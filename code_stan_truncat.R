rm(list=ls())
library(ggsci)
library(rstan)
library(splines)
library(tidyverse)
library(bayesplot)

# --- 全局配置 ---
sample_size <- 530  # 修改此处即可自由切换：50, 100, 530
run_type <- ifelse(sample_size < 500, "test", "full")

n_iter <- ifelse(run_type == "test", 500, 2000)
n_warmup <- n_iter / 2
target_delta <- ifelse(run_type == "test", 0.8, 0.95)

# --- 数据处理 ---
source('code_data_explore.R')

df_sero_sub <- df_sero.3 %>%
  mutate(DOB = as.Date(DOB),
         DOC = as.Date(DOC),
         DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1) %>%
  slice_head(n = sample_size)

df_pos.2 <- df_pos.2 %>% arrange(MONTH)
lambda_vec <- df_pos.2$pos_scale
M_total <- length(lambda_vec)

# --- B-Spline 与 索引矩阵 ---
knots_internal <- c(4, 8, 12, 16)
B_basis <- bs(0:23, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, 23), intercept = TRUE)

N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = 24)

for (i in 1:N_current) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:23)
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1

# --Bg value ---
df_sero_sub[, c("IgA-MFI","IgG-MFI","SiteIII-MPE8", "Site II-Motavizumab",
                "Site V-Suptavumab", "Site Φ-Nirsevimab",
                "Site IV-Clesrovimab")] %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  map_dfr(~quantile(.x,c(1:10/10),na.rm = T))

df_sero_sub %>%
  select(`IgA-MFI`, `IgG-MFI`, `SiteIII-MPE8`, `Site II-Motavizumab`,
         `Site V-Suptavumab`, `Site Φ-Nirsevimab`, `Site IV-Clesrovimab`) %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  facet_wrap(~Variable, scales = "free", ncol = 3) +
  scale_x_log10() +  # Log transform due to wide range
  labs(title = "Distribution of Antibody/Antigen Measurements",
       x = "Value (log10 scale)", y = "Count") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8, face = "bold"))

BG_VALUE<-df_sero_sub[, c("IgA-MFI","IgG-MFI","SiteIII-MPE8", "Site II-Motavizumab",
                "Site V-Suptavumab", "Site Φ-Nirsevimab",
                "Site IV-Clesrovimab")] %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  map_dbl(~quantile(.x,c(0.1),na.rm = T))


# --- 构造 Stan Data ---
stan_data <- list(
  N = N_current,
  M = M_total,
  A = as.integer(df_sero_sub$age_value),
  cal_idx = cal_idx_matrix,
  Y_IgA = df_sero_sub$`IgA-MFI`,
  Y_IgG = df_sero_sub$`IgG-MFI`,
  Y_epi = as.matrix(df_sero_sub[, c("SiteIII-MPE8", "Site II-Motavizumab",
                                    "Site V-Suptavumab", "Site Φ-Nirsevimab",
                                    "Site IV-Clesrovimab")]),
  lambda = lambda_vec,
  r_mat = 20.8/36,
  r1=20.8/(36/2),
  r2=20.8/34,
  r_epi=rep(20.8/36,5),

  J = ncol(B_basis),
  B_basis = B_basis,

  B_IgA_base=BG_VALUE[1],
  B_IgG_base=BG_VALUE[2],
  B_epi_base=BG_VALUE[3:7]
)

# --- 运行模型 ---
options(mc.cores = parallel::detectCores())
cat("🚀 启动样本量:", N_current, " | 模式:", run_type, "\n")

# Truncat model ----

my_inits <- function() {
  list(
    beta = rep(0.5, 8),
    tau = 0.5,
    pi = rep(0.85, 7),                             # 敏感度初始设为 85%，很合理
    peak = c(500, 5000, 1.5, 1.5, 1.5, 1.5, 1.5),  # 区分 MFI 尺度和 Logit 尺度
    Mat0 = 150,
    Mat0_k = rep(10, 5),
    sigma1 = 1.5,
    sigma2 = 1.5,
    phi = rep(5, 5),
    B_IgA_raw = 0,                                 # 扰动参数初始不扰动
    B_IgG_raw = 0,
    B_epi_raw = rep(0, 5)
  )
}


# fit_truncat <- stan(
#   file = "rsv_model_truncat.stan",
#   data = stan_data,
#   chains = 4,
#   iter = n_iter,
#   warmup = n_warmup,
#   init = my_inits,
#   control = list(adapt_delta = target_delta, max_treedepth = 12)
# )
#
# saveRDS(fit_truncat,'rda/fit_truncat.rds')

fit_truncat<-readRDS('rda/fit_truncat.rds')
# --- 结果分析 ---
fit_truncat

cat(capture.output(print(fit_truncat)), file = "fit_truncat_output.txt", sep = "\n")

# 易感性曲线 w(a) 提取
posterior_w_truncat <- rstan::extract(fit_truncat, "w")$w
w_res_truncat <- data.frame(
  age = 0:23,
  mean = colMeans(posterior_w_truncat),
  lower = apply(posterior_w_truncat, 2, quantile, 0.025),
  upper = apply(posterior_w_truncat, 2, quantile, 0.975)
)

# 可视化
ggplot(w_res_truncat, aes(x = age, y = mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  labs(#title = paste("RSV Susceptibility (N =", N_current, ")"),
       x = "Age (months)", y = "Relative Susceptibility") +
  scale_x_continuous(breaks=seq(0,23,by=3),expand = expansion())+
  scale_y_continuous(expand = expansion())+
  theme_bw()+
  theme()

ggsave('plot/truncat/w_by_age.jpg',width=8,height=6)

fit_truncat@par_dims %>% map_dbl(~length(.x))

color_scheme_set("viridisA")
walk2(fit_truncat@model_pars, fit_truncat@par_dims, ~{

  # print(.x)
  # print(.y)
  # print(length(.y))

  # 标量参数
  if (length(.y) == 0) {
    p <- mcmc_trace(
      fit_truncat,
      pars = .x
    )
    ggsave(
      filename = paste0("plot/truncat/", .x, ".jpg"),
      plot = p,
      width = 8,
      height = 5
    )

  } else {
    # 非标量参数，生成 Stan 风格参数名
    idx <- expand.grid(lapply(.y, seq_len))
    pars_use <- apply(idx, 1, function(i) {
      paste0(.x, "[", paste(i, collapse = ","), "]")
    })
    p <- mcmc_trace(
      fit_truncat,
      pars = pars_use
    )
    ggsave(
      filename = paste0("plot/truncat/", .x, ".jpg"),
      plot = p,
      width = 10,
      height = 6
    )
  }
})

df_means<-enframe(get_posterior_mean(fit_truncat)[, 5])

df_real_peaks <- df_means %>%
  pivot_wider(names_from = name, values_from = value) %>%
  pivot_longer(cols = starts_with("peak["), names_to = "peak_name", values_to = "stan_peak_delta") %>%
  mutate(
    idx = as.integer(str_extract(peak_name, "\\d+")),
    marker = c("IgA-MFI", "IgG-MFI", "SiteIII-MPE8", "Site II-Motavizumab",
               "Site V-Suptavumab", "Site Φ-Nirsevimab", "Site IV-Clesrovimab")[idx]
  ) %>%
  rowwise() %>%
  mutate(
    real_peak_value = case_when(
      idx == 1 ~ B_IgA + stan_peak_delta,
      idx == 2 ~ B_IgG + Mat0 + stan_peak_delta,
      idx == 3 ~ boot::inv.logit(boot::logit(`B_epi[1]`) + stan_peak_delta),
      idx == 4 ~ boot::inv.logit(boot::logit(`B_epi[2]`) + stan_peak_delta),
      idx == 5 ~ boot::inv.logit(boot::logit(`B_epi[3]`) + stan_peak_delta),
      idx == 6 ~ boot::inv.logit(boot::logit(`B_epi[4]`) + stan_peak_delta),
      idx == 7 ~ boot::inv.logit(boot::logit(`B_epi[5]`) + stan_peak_delta)
    )
  ) %>%
  ungroup() %>%
  select(marker,B_IgA,B_IgG,starts_with('B_epi['),Mat0,stan_peak_delta, real_peak_value) %>%
  bind_cols(BG_VALUE=BG_VALUE)

print(df_real_peaks)


# Risk of infection ----
df_individual_raw <- df_sero_sub %>%
  mutate(ID = row_number()) %>%
  bind_cols(
    df_means %>%
      filter(str_detect(name, "^w\\[")) %>%
      mutate(m = as.integer(str_extract(name, "\\d+"))) %>%
      select(m, value) %>%
      pivot_wider(names_from = m, names_prefix = "w_", values_from = value) %>%
      slice(rep(1, nrow(df_sero_sub)))
  ) %>%
  rowwise() %>%
  mutate(
    age_i = as.integer(age_value),
    start_m = max(1, age_i - 6),
    eff_len = age_i - start_m
  ) %>%
  filter(eff_len > 0) %>%
  mutate(
    timeline_data = list({
      w_vec <- c(w_1, w_2, w_3, w_4, w_5, w_6, w_7, w_8, w_9, w_10, w_11, w_12,
                 w_13, w_14, w_15, w_16, w_17, w_18, w_19, w_20, w_21, w_22, w_23, w_24)
      h <- numeric(eff_len)
      S_surv_win <- numeric(eff_len + 1)
      S_surv_win[1] <- 1.0
      prob_t <- numeric(eff_len)

      for (idx in 1:eff_len) {
        m <- start_m + idx - 1
        c_idx <- cal_idx_matrix[cur_group_rows(), m]
        h[idx] <- lambda_vec[c_idx] * w_vec[m]
        S_surv_win[idx + 1] <- S_surv_win[idx] * exp(-h[idx])
        prob_t[idx] <- S_surv_win[idx] * -expm1(-h[idx])
      }

      tibble(
        absolute_age_month = start_m : (age_i - 1),
        prob_infection_at_this_month = prob_t
      )
    }),
    prob_uninfected = {
      w_vec <- c(w_1, w_2, w_3, w_4, w_5, w_6, w_7, w_8, w_9, w_10, w_11, w_12,
                 w_13, w_14, w_15, w_16, w_17, w_18, w_19, w_20, w_21, w_22, w_23, w_24)
      h_sum <- 0
      for (idx in 1:eff_len) {
        m <- start_m + idx - 1
        c_idx <- cal_idx_matrix[cur_group_rows(), m]
        h_sum <- h_sum + lambda_vec[c_idx] * w_vec[m]
      }
      exp(-h_sum)
    }
  ) %>%
  ungroup()

df_individual_absolute_timeline<-df_individual_raw %>%
  select(ID, age_value, prob_uninfected, timeline_data) %>%
  unnest(timeline_data)

print(df_individual_absolute_timeline)

df_individual_wide <- df_individual_absolute_timeline %>%
  mutate(absolute_age_month = factor(absolute_age_month, levels = 0:23)) %>%
  pivot_wider(
    names_from = absolute_age_month,         # 把绝对月份变成新列名
    names_prefix = "prob_infect_age_",       # 给列名加个清爽的前缀
    values_from = prob_infection_at_this_month, # 填充的数值
    values_fill = NA_integer_                       # 不在窗口内的月份自动填 0
  ) %>%
  arrange(age_value)

# 打印查看一人一行的完美大宽表
print(df_individual_wide)

export(df_individual_wide,'docs/df_individual_wide.xlsx')

df_individual_wide %>%
  summarise(expected_uninfected_prob = mean(prob_uninfected, na.rm = TRUE),.by=age_value) %>%
  ggplot(aes(age_value,1-expected_uninfected_prob))+
  geom_point()+
  geom_line(group=1)

df_individual_wide %>%
  summarise(
    across(starts_with("prob_infect_age_"), ~mean(.x, na.rm = TRUE), .names = "{.col}")
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "expected_prob"
  ) %>%
  mutate(age = as.integer(str_remove(age, "prob_infect_age_"))) %>%
  arrange(age) %>%
  ggplot(aes(x = age, y = expected_prob)) +
  geom_line(color = "#2471A3", size = 1.2) +
  geom_point(color = "#E74C3C", size = 3) +
  scale_x_continuous(breaks = 1:24) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Expected Primary RSV Infection Probability by Absolute Age",
    subtitle = paste("Based on Latent Survival Analysis (N =", nrow(df_individual_wide), ")"),
    x = "Absolute Age (months)",
    y = "Expected Probability of Primary Infection"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5),
    panel.grid.minor = element_blank()
  )
