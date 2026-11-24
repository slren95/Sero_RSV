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
B_basis <- bs(1:24, knots = knots_internal, degree = 3,
              Boundary.knots = c(1, 24), intercept = TRUE)

N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = 23)

for (i in 1:N_current) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:22)
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1

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
  r_mat = 0.58,
  J = ncol(B_basis),
  B_basis = B_basis
)

# --- 运行模型 ---
options(mc.cores = parallel::detectCores())
cat("🚀 启动样本量:", N_current, " | 模式:", run_type, "\n")

# Stan model ----

# fit <- stan(
#   file = "rsv_model.stan",
#   data = stan_data,
#   chains = 4,
#   iter = n_iter,
#   warmup = n_warmup,
#   init = 0,
#   control = list(adapt_delta = target_delta, max_treedepth = 12)
# )

fit<-readRDS('rda/fit.rds')
#saveRDS(fit,'rda/fit.RDS')

# --- 结果分析 ---
print(fit, pars = c("w","pi","B_IgG","B_IgA","peak", "r1", "r2","Mat0","Mat0_k","pi"))

# 易感性曲线 w(a) 提取
posterior_w <- rstan::extract(fit, "w")$w
w_res <- data.frame(
  age = 1:24,
  mean = colMeans(posterior_w),
  lower = apply(posterior_w, 2, quantile, 0.025),
  upper = apply(posterior_w, 2, quantile, 0.975)
)

# 可视化
ggplot(w_res, aes(x = age, y = mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  labs(title = paste("RSV Susceptibility (N =", N_current, ")"),
       x = "Age (months)", y = "Relative Susceptibility") +
  theme_minimal()

#ggsave('plot/w_by_age.jpg',width=10,height=6)




color_scheme_set("viridis")
walk2(fit@model_pars, fit@par_dims, ~{

  # print(.x)
  # print(.y)
  # print(length(.y))

  # 标量参数
  if (length(.y) == 0 || all(.y == 1)) {
    p <- mcmc_trace(
      fit,
      pars = .x
    )
    ggsave(
      filename = paste0("plot/", .x, ".jpg"),
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
      fit,
      pars = pars_use
    )
    ggsave(
      filename = paste0("plot/", .x, ".jpg"),
      plot = p,
      width = 10,
      height = 6
    )
  }
})

# Truncat model ----

# fit_truncat <- stan(
#   file = "rsv_model_truncat.stan",
#   data = stan_data,
#   chains = 4,
#   iter = n_iter,
#   warmup = n_warmup,
#   init = 0,
#   control = list(adapt_delta = target_delta, max_treedepth = 12)
# )
#
# saveRDS(fit_truncat,'rda/fit_truncat.rds')

fit_truncat<-readRDS('rda/fit_truncat.rds')

writeLines(capture.output(print(fit_truncat)), 'fit_truncat_output.txt')
# --- 结果分析 ---
print(fit_truncat, pars = c("w","pi","B_IgG","B_IgA","peak", "r1", "r2","Mat0","Mat0_k","pi"))

# 易感性曲线 w(a) 提取
posterior_w_truncat <- rstan::extract(fit_truncat, "w")$w
w_res_truncat <- data.frame(
  age = 1:24,
  mean = colMeans(posterior_w_truncat),
  lower = apply(posterior_w_truncat, 2, quantile, 0.025),
  upper = apply(posterior_w_truncat, 2, quantile, 0.975)
)

# 可视化
ggplot(w_res_truncat, aes(x = age, y = mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  labs(title = paste("RSV Susceptibility (N =", N_current, ")"),
       x = "Age (months)", y = "Relative Susceptibility") +
  theme_minimal()

ggsave('plot/truncat_w_by_age.jpg',width=10,height=6)

walk2(fit_truncat@model_pars, fit_truncat@par_dims, ~{

  # print(.x)
  # print(.y)
  # print(length(.y))

  # 标量参数
  if (length(.y) == 0 || all(.y == 1)) {
    p <- mcmc_trace(
      fit_truncat,
      pars = .x
    )
    ggsave(
      filename = paste0("plot/truncat_", .x, ".jpg"),
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
      filename = paste0("plot/truncat_", .x, ".jpg"),
      plot = p,
      width = 10,
      height = 6
    )
  }
})
