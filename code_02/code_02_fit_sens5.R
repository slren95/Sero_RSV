#rm(list=ls())
library(ggsci)
library(rstan)
library(splines)
library(tidyverse)
library(bayesplot)
library(posterior)
library(rio)

# Model Detail ---
# Model Name: fit_sens5 (v6.0 极简版)
# Description: Sensitivity analysis - IgA + IgG ONLY (No Epitopes, No Toggles)
# Age range: 0-29

# --- 1. 数据处理 ---
load('rda/code_01_data_explore.RData')

df_sero_sub <- df_sero.4 %>%
  mutate(DOB = as.Date(DOB),
         DOC = as.Date(DOC),
         DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1) %>%
  filter(age_value < 30)

df_pos.2 <- df_pos.2 %>% arrange(MONTH)
lambda_vec <- df_pos.2$pos_scale
M_total <- length(lambda_vec)

# --- 2. B-Spline 与 索引矩阵 ---
knots_internal <- c(4, 8, 12, 16, 20, 24)
B_basis <- bs(0:29, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, 29), intercept = TRUE)

N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = 30)

for (i in 1:N_current) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:29)
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1


# --- 3. 构造 适配 v6.0 Stan 代码的 stan_data ----
# 剔除了: Y_epi, use_IgA, use_IgG, use_epi, r_epi 以及 DF_t12
stan_data <- list(
  N = N_current,
  M = M_total,
  A = as.integer(df_sero_sub$age_value),   cal_idx = cal_idx_matrix,      # 观测数据 (仅双抗体)
  Y_IgA = df_sero_sub$`IgA-MFI`,
  Y_IgG = df_sero_sub$`IgG-MFI`,

  # 外部已知参数与衰减率
  lambda = lambda_vec,
  r_mat = 20.8 / 36,
  r1    = 20.8 / (34 / 2),
  r2    = 20.8 / 34,

  J = ncol(B_basis),
  B_basis = B_basis
)


# --- 4. 适配 v6.0 Stan 代码的 Initial Value Generator ----
# 剔除了: Mat0_k, phi, B_epi；将 pi 和 peak 数组长度缩减为 2
my_inits <- function() {
  list(
    beta   = rep(0.5, ncol(B_basis)),
    tau    = 0.5,
    pi     = c(0.85, 0.85),        # 仅保留 IgA, IgG (长度为 2)
    peak   = c(500, 5000),         # 仅保留 IgA, IgG (长度为 2)
    Mat0   = 150,
    sigma1 = 1.5,
    sigma2 = 1.5,
    B_IgA  = runif(1, 4, 6),
    B_IgG  = runif(1, 40, 60)
  )
}

# --- 5. Run Bayes Model ---
options(mc.cores = parallel::detectCores())

# fit_sens5 <- stan(
#   file = "stan/fit_sens5.stan",  # 指向极简版的 v6.0 stan 文件
#   data = stan_data,
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   init = my_inits,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )

# 方案：单核串行跑 4 条链（日志实时输出，适合观察每条链的运行状态）
# fit_sens5 <- stan(
#   file = "stan/fit_sens5.stan",
#   data = stan_data,
#   chains = 4,        # 跑满 4 条链（确保能计算 Rhat）
#   cores = 1,         # 强制单核运行，日志实时打印在控制台
#   iter = 2000,
#   warmup = 1000,
#   refresh = 50,      # 每 50 步刷新一次进度
#   init = my_inits,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
#
# saveRDS(fit_sens5, 'rda/fit_sens5.rds')
