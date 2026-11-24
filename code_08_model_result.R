# 单独查看新版主分析结果；从项目根目录 D:/NJMU/Sero_RSV 运行。
# 模型定义：stan_new/fit_main_new.stan
# 拟合脚本：code_02_new/code_02_fit_main_new.R
# 拟合结果：rda_new/fit_main_new.rds（由上述拟合脚本保存）
# 感染概率输入：code_03_new/fit_main_new/df_individual_wide_all.rds，
# 由 code_03_new/code_03.1_stan_plot_new.R 生成；需与上述主模型拟合结果对应。
# 本脚本只读取已有结果，不重新拟合；不依赖 code_03/code_04 的运行环境。
library(rstan)
library(dplyr)
library(tidyr)

MODEL_NAME <- "fit_main_new"
fit_path <- file.path("rda_new", paste0(MODEL_NAME, ".rds"))
if (!file.exists(fit_path)) {
  stop("找不到主分析结果：", fit_path, "；请检查项目工作目录和拟合结果文件。")
}
fit_truncat <- readRDS(fit_path)

# 提取 pi 的后验抽样
pi_draws <- rstan::extract(fit_truncat, pars = "pi")$pi

# 计算后验中位数和 95% CrI
pi_summary <- t(apply(
  pi_draws, 2, quantile,
  probs = c(0.5, 0.025, 0.975)
))

# 转换为百分比
pi_summary <- round(pi_summary * 100, 1)

# 设置名称
rownames(pi_summary) <- c(
  "IgA", "IgG", "Site III", "Site II",
  "Site V", "Site Ø", "Site IV"
)

colnames(pi_summary) <- c("Median", "Lower", "Upper")

# 查看结果
print(pi_summary)

# Age-specific relative susceptibility ----
# w[1] 对应月龄 0，w[25] 对应月龄 24；报告后验中位数及等尾 95% CrI。
# Figure 2 的易感性曲线原代码使用后验均值；此处按正文要求使用中位数。
w_draws <- rstan::extract(fit_truncat, pars = "w")$w
stopifnot(ncol(w_draws) >= 25)
w_quantiles <- t(apply(w_draws, 2, quantile, probs = c(0.5, 0.025, 0.975)))
w_summary <- data.frame(
  age_month = seq_len(ncol(w_draws)) - 1L,
  median = w_quantiles[, 1], lower = w_quantiles[, 2], upper = w_quantiles[, 3]
)

# 最高/最低点按 0–23 月龄的后验中位数曲线确定；并列时保留全部月龄。
# 这些点的区间是对应月龄的点位 CrI，不是极值位置本身的可信区间。
w_first24 <- filter(w_summary, age_month <= 23)
w_peak <- filter(w_first24, median == max(median))
w_nadir <- filter(w_first24, median == min(median))
w_key <- bind_rows(
  mutate(filter(w_summary, age_month == 0), item = "Month 0"),
  mutate(w_peak, item = "Highest median within months 0–23"),
  mutate(w_nadir, item = "Lowest median within months 0–23"),
  mutate(filter(w_summary, age_month == 24), item = "Month 24")
) %>%
  mutate(report = sprintf("%.2f (95%% CrI, %.2f–%.2f)", median, lower, upper)) %>%
  select(item, age_month, report)

cat("\nAge-specific relative susceptibility to RSV infection\n")
print(w_key, row.names = FALSE)

# Age-dependent probability of RSV infection ----
# 沿用 code_04_new/code_04.2_plot_figure2_new.R 的口径：
# 每个后验抽样内，对该月的个体概率求均值（按原代码将 0 视为未纳入）；
# 随后按年龄排序，逐抽样计算 1 - cumprod(1 - mean_prob)，最后求分位数。
# 先对宽表分组汇总再转长表，避免把全部个体抽样展开成很大的长表。
risk_path <- file.path("code_03_new", MODEL_NAME, "df_individual_wide_all.rds")
if (!file.exists(risk_path)) {
  stop("找不到感染概率后处理结果：", risk_path,
       "；请先为该模型运行 code_03_new/code_03.1_stan_plot_new.R。")
}
risk_data <- readRDS(risk_path)
prob_columns <- paste0("prob_infect_age_", 0:23)
stopifnot(all(c("post_id", prob_columns) %in% names(risk_data)))

risk_draws <- risk_data %>%
  group_by(post_id) %>%
  summarise(across(all_of(prob_columns), ~ mean(na_if(.x, 0), na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_longer(cols = all_of(prob_columns), names_to = "age_month", values_to = "mean_prob") %>%
  mutate(age_month = as.integer(sub("prob_infect_age_", "", age_month))) %>%
  arrange(post_id, age_month)
if (any(!is.finite(risk_draws$mean_prob))) {
  stop("0–23 月龄内存在无有效个体概率的后验抽样，无法完整计算累计概率。")
}
risk_draws <- risk_draws %>%
  group_by(post_id) %>%
  mutate(cum_prob = 1 - cumprod(1 - mean_prob)) %>%
  ungroup()

risk_summary <- risk_draws %>%
  group_by(age_month) %>%
  summarise(across(c(mean_prob, cum_prob), list(
    median = ~ quantile(.x, 0.5),
    lower = ~ quantile(.x, 0.025),
    upper = ~ quantile(.x, 0.975)
  )), .groups = "drop") %>%
  arrange(age_month)

monthly_peak <- filter(risk_summary, mean_prob_median == max(mean_prob_median))
monthly_nadir <- filter(risk_summary, mean_prob_median == min(mean_prob_median))
monthly_key <- bind_rows(
  mutate(filter(risk_summary, age_month == 0), item = "First month"),
  mutate(monthly_peak, item = "Highest median within months 0–23"),
  mutate(monthly_nadir, item = "Lowest median within months 0–23"),
  mutate(filter(risk_summary, age_month == 23), item = "23–<24 months")
) %>%
  mutate(age_interval = sprintf("%d–<%d months", age_month, age_month + 1L),
         report = sprintf("%.1f%% (95%% CrI, %.1f–%.1f%%)",
                          100 * mean_prob_median, 100 * mean_prob_lower, 100 * mean_prob_upper)) %>%
  select(item, age_interval, report)

# 截至 6、12、24 月：分别累积 0–5、0–11、0–23 的月份区间。
cumulative_key <- risk_summary %>%
  filter(age_month %in% c(5, 11, 23)) %>%
  mutate(by_age_month = age_month + 1L,
         report = sprintf("%.1f%% (95%% CrI, %.1f–%.1f%%)",
                          100 * cum_prob_median, 100 * cum_prob_lower, 100 * cum_prob_upper)) %>%
  select(by_age_month, report)

cat("\nMonthly probability of RSV infection\n")
print(as.data.frame(monthly_key), row.names = FALSE)
cat("\nCumulative probability of RSV infection\n")
print(as.data.frame(cumulative_key), row.names = FALSE)

# by_age_month                      report
# 6 66.4% (95% CrI, 59.1–73.0%)
# 12 81.2% (95% CrI, 76.8–85.0%)
# 24 95.7% (95% CrI, 94.6–96.7%)

# 工作区保留完整数字表：w_summary、risk_summary（未四舍五入、概率为 0–1）。
# 正文所需摘要：w_key、monthly_key、cumulative_key；pi_summary 为指标敏感度。

