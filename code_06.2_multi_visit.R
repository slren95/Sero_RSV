# 状态：未使用，已放弃；当前分析无必要采用本方案。
# 与 code_06.1_multi_visit.R 的区别：
# 1. 06.1 使用 B_epi 的后验均值；本脚本最多抽取 1000 组 B_epi 后验样本，逐组计算衰减率。
# 2. 本脚本在筛选 include == 1 后重新排序并计算随访间隔，使间隔与纳入记录的前后配对一致。
# 3. 06.1 要求五个表位同时高于背景，再计算 lag；本脚本按表位独立处理，先配对，再筛选前后均高于背景的记录。
# 4. 06.1 拟合衰减率分布并导出 epi_prior_stan.rds；本脚本仅汇总、绘图，不生成该先验文件。
# 当前分析沿用 06.1 生成的先验文件；本脚本保留作为历史探索记录。

rm(list = ls())

library(tidyverse)
library(rio)
library(psych)
library(tidyxl)
library(janitor)
library(readxl)
library(patchwork)
library(openxlsx)
library(furrr)

# ============================================================
# 0. 全局设置
# ============================================================
site_cols <- c("site_iii", "site_ii", "site_v", "site_ph", "site_iv")

N_DRAWS_USE <- 1000
set.seed(20260619)

future::plan(multisession, workers = max(1, parallel::detectCores() - 1))


# ============================================================
# 1. 读取并清洗主数据
# ============================================================
df_new <- read_excel("data/KKH-V2-RSV Clinical Dataset for Model Analysis-20260619_new.xlsx") %>%
  clean_names() %>%
  mutate(across(6:9, ~ as.Date(as.numeric(.x), origin = "1899-12-30"))) %>%
  mutate(visit_1_date = ymd(visit_1_date)) %>%
  mutate(subject_ID = str_sub(sample_visit, 1, 8), .before = 1) %>%
  separate(age, into = c("age_value", "age_unit"), remove = FALSE) %>%
  mutate(
    age_unit = tolower(age_unit),
    age_value = as.integer(age_value),
    DOB = case_when(
      age_unit == "years" ~ visit_1_date %m-% years(age_value),
      age_unit == "months" ~ visit_1_date %m-% months(age_value),
      age_unit == "weeks" ~ visit_1_date %m-% weeks(age_value),
      age_unit == "days" ~ visit_1_date %m-% days(age_value),
      TRUE ~ NA_Date_
    )
  ) %>%
  mutate(
    age_month_1 = interval(DOB, visit_1_date) / months(1),
    age_month_2 = interval(DOB, visit_2_date) / months(1),
    age_month_3 = interval(DOB, visit_3_date) / months(1),
    age_month_4 = interval(DOB, visit_4_date) / months(1),
    age_month_5 = interval(DOB, visit_5_date) / months(1)
  ) %>%
  mutate(
    age_visit = case_when(
      visit == 1 ~ age_month_1,
      visit == 2 ~ age_month_2,
      visit == 3 ~ age_month_3,
      visit == 4 ~ age_month_4,
      visit == 5 ~ age_month_5,
      TRUE ~ NA_real_
    ),
    .after = age_unit
  ) %>%
  group_by(subject_ID) %>%
  mutate(
    visit_interval = age_visit - lag(age_visit),
    subject_include_site = as.integer(sum(include, na.rm = TRUE) > 0)
  ) %>%
  ungroup() %>%
  rename_with(~ str_extract(.x, "^[^_]+_[^_]+"), .cols = starts_with("site")) %>%
  relocate(DOB, starts_with("age_month_"), .after = visit_5_date) %>%
  relocate(subject_include_site, .after = include)


# ============================================================
# 2. 分析数据
#
# 筛选后重新计算 visit_interval，
# 保证 interval 与实际进入分析的相邻 visit 一致
# ============================================================
df_new.2 <- df_new %>%
  filter(age_visit <= 30, include == 1) %>%
  arrange(subject_ID, age_visit) %>%
  group_by(subject_ID) %>%
  mutate(visit_interval = age_visit - lag(age_visit)) %>%
  ungroup() %>%
  add_count(subject_ID)


# ============================================================
# 3. 提取 B_epi posterior draws
# ============================================================
fit_simple <- readRDS("rda/fit_simple.rds")

B_epi_draws_full <- rstan::extract(fit_simple, pars = "B_epi")$B_epi
colnames(B_epi_draws_full) <- site_cols

S_all <- nrow(B_epi_draws_full)
idx <- sample.int(S_all, size = min(N_DRAWS_USE, S_all))

B_epi_draws <- B_epi_draws_full[idx, , drop = FALSE]

# posterior mean，仅用于参考
B_epi_mean <- colMeans(B_epi_draws_full)

B_epi_mean
qlogis(B_epi_mean)


# ============================================================
# 4. 转成长格式
#
# 每个 site 独立处理
# ============================================================
df_long <- df_new.2 %>%
  select(subject_ID, visit, age_visit, visit_interval, all_of(site_cols)) %>%
  pivot_longer(cols = all_of(site_cols), names_to = "site", values_to = "value") %>%
  arrange(subject_ID, site, age_visit)


# ============================================================
# 5. 单个 posterior BG draw 下计算所有 decay rates
#
# 关键：
# 1. 每个 site 使用自己的 background
# 2. 每个 site 独立判断是否 > background
# 3. lag 在 background 筛选之前计算
# 4. visit_interval 对应同一对实际相邻 visit
# ============================================================
compute_decay_one_draw <- function(B_epi_vec, data_long = df_long) {

  # B_epi_vec：当前一个 posterior draw 中 5 个 site 的 background values
  # data_long：长格式数据，每行对应一个 subject × visit × site

  data_long %>%

    # 1. 为每个 site 匹配当前 posterior draw 对应的 background，
    #    并计算 background-adjusted logit response
    mutate(
      bg_value = B_epi_vec[site],
      adjusted_logit = qlogis(value / 100) - qlogis(bg_value)
    ) %>%

    # 2. 按 subject、site 和年龄排序，确保后续 lag() 对应真实的前一次 visit
    arrange(subject_ID, site, age_visit) %>%

    # 3. 每个 subject × site 独立计算
    group_by(subject_ID, site) %>%

    # 4. 获取当前 visit 对应的前一次 visit 信息
    #    previous_adjusted_logit 就是前一次 visit 的 adjusted_logit
    mutate(
      previous_visit = lag(visit),
      previous_age_visit = lag(age_visit),
      previous_adjusted_logit = lag(adjusted_logit)
    ) %>%

    # 5. 只有相邻两次 visit 的 adjusted logit 均 > 0 时才计算
    #    adjusted_logit > 0 等价于原始 response 高于当前 draw 的 background
    filter(
      !is.na(previous_adjusted_logit),
      !is.na(adjusted_logit),
      previous_adjusted_logit > 0,
      adjusted_logit > 0,
      !is.na(visit_interval),
      visit_interval > 0
    ) %>%

    # 6. 根据相邻两次 visit 的 response 变化计算衰减率
    #
    #    previous_adjusted_logit = 前一次 visit 的 response
    #    adjusted_logit          = 当前 visit 的 response
    #
    #    k = -[log(R_current) - log(R_previous)] / Δt
    #
    #    k > 0：response 下降
    #    k < 0：response 上升
    mutate(
      k = -(log(adjusted_logit) - log(previous_adjusted_logit)) / visit_interval,
      t12 = log(2) / k
    ) %>%

    ungroup() %>%

    # 7. 保留每一个 subject × site × 相邻 visit pair 的计算结果
    select(
      subject_ID, site,
      previous_visit, visit,
      previous_age_visit, age_visit, visit_interval,
      bg_value,
      previous_adjusted_logit, adjusted_logit,
      k, t12
    )
}


# ============================================================
# 6. 检查第 1 个 posterior draw
# ============================================================
B_epi_draw_1 <- B_epi_draws[1, ]

df_decay_draw1 <- compute_decay_one_draw(B_epi_draw_1)

df_decay_draw1


# ============================================================
# 7. 第 1 个 draw：decay rate 分布
# ============================================================
ggplot(df_decay_draw1, aes(x = k)) +
  geom_histogram(bins = 30, color = "black") +
  facet_wrap(~ site, scales = "free") +
  theme_bw() +
  labs(
    title = "Distribution of individual decay rates (Draw 1)",
    x = "Decay rate k (per month)",
    y = "Count"
  )


# ============================================================
# 8. 第 1 个 draw：half-life 分布
# ============================================================
ggplot(df_decay_draw1, aes(x = t12)) +
  geom_histogram(bins = 30, color = "black") +
  facet_wrap(~ site, scales = "free") +
  theme_bw() +
  labs(
    title = "Distribution of individual half-lives (Draw 1)",
    x = "Half-life (months)",
    y = "Count"
  )


# ============================================================
# 9. 对 1000 个 posterior BG draws 全部计算
#
# index = 1:1000
# posterior_draw = 对应原始 Stan posterior draw 编号
#
# 不做任何 individual-level 汇总
# ============================================================
decay_draw_list <- furrr::future_map(
  seq_len(nrow(B_epi_draws)),
  function(s) {
    compute_decay_one_draw(B_epi_draws[s, ]) %>%
      mutate(index = s, posterior_draw = idx[s], .before = 1)
  },
  .options = furrr::furrr_options(seed = TRUE)
)

DF_decay_post <- bind_rows(decay_draw_list)

# ============================================================
# 1. 每个 site 中非正 decay rate 的数量和比例
# ============================================================
DF_negative_summary <- DF_decay_post %>%
  group_by(site) %>%
  summarise(
    n_total = n(),
    n_k_nonpositive = sum(k <= 0, na.rm = TRUE),
    prop_k_nonpositive = mean(k <= 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(site = factor(site, levels = site_cols)) %>%
  arrange(site)

DF_negative_summary


# ============================================================
# 2. 仅保留真正的 waning observations：k > 0
# ============================================================
DF_decay_positive <- DF_decay_post %>%
  filter(k > 0)


# ============================================================
# 3. 对 k > 0 的数据计算 k 和 t12 的分布
# ============================================================
DF_decay_summary <- DF_decay_positive %>%
  group_by(site) %>%
  summarise(
    n = n(),
    k_median = median(k, na.rm = TRUE),
    k_lower = quantile(k, 0.025, na.rm = TRUE),
    k_upper = quantile(k, 0.975, na.rm = TRUE),
    t12_median = median(t12, na.rm = TRUE),
    t12_lower = quantile(t12, 0.025, na.rm = TRUE),
    t12_upper = quantile(t12, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = factor(site, levels = site_cols)) %>%
  arrange(site)

DF_decay_summary


# ============================================================
# 4. 合并：同时展示排除比例和最终结果
# ============================================================
DF_decay_summary_all <- DF_negative_summary %>%
  left_join(DF_decay_summary, by = "site")

DF_decay_summary_all


# 绘图 ----
# k
p_k <- DF_decay_positive %>%
  mutate(site = factor(site, levels = site_cols)) %>%
  ggplot(aes(x = k)) +
  geom_histogram(bins = 50, color = "black") +
  facet_wrap(~ site, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Decay rate (per month)", y = "Count")

p_k


# t12
p_t12 <- DF_decay_positive %>%
  mutate(site = factor(site, levels = site_cols)) %>%
  ggplot(aes(x = t12)) +
  geom_histogram(bins = 50, color = "black") +
  facet_wrap(~ site, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Half-life (months)", y = "Count")

p_t12
