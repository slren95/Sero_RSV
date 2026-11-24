#rm(list=ls())
library(ggsci)
library(rstan)
library(splines)
library(tidyverse)
library(bayesplot)
library(posterior)
library(rio)

# Model Detail ---
# Data: df_sero.4
# Age range: 0-29
# Waning of epitope: from longitudinal data DF_t12

# Sens1: 衰减速率允许变异

# --- 数据处理 ---
load('rda/code_01_data_explore.RData')
DF_t12<-import('rda/DF_t12.rds',trust=T)

df_sero_sub <- df_sero.4 %>%
  mutate(DOB = as.Date(DOB),
         DOC = as.Date(DOC),
         DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1) %>%
  filter(age_value < 30)


df_pos.2 <- df_pos.2 %>% arrange(MONTH)
lambda_vec <- df_pos.2$pos_scale
M_total <- length(lambda_vec)

# --- B-Spline 与 索引矩阵 ---
knots_internal <- c(4, 8, 12, 16, 20, 24)
B_basis <- bs(0:29, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, 29), intercept = TRUE)

N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = 30)

for (i in 1:N_current) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:29)
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1


# --- 1. 基础 Stan Data 模板 ----
stan_data_base <- list(
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

  # 先验均值配置
  prior_r_mat_mean = 20.8 / 36,
  prior_r1_mean     = 20.8 / (34 / 2),
  prior_r2_mean     = 20.8 / 34,
  prior_r_epi_mean  = log(2) / DF_t12$median,

  J = ncol(B_basis),
  B_basis = B_basis
)

# --- 2. 构造 3 种敏感性分析 Data 列表 ----
sens_scenarios <- list(
  "sens1_0.05_0.05" = c(main = 0.05, epi = 0.05), # 基线场景
  "sens1_0.05_0.15" = c(main = 0.05, epi = 0.15), # 仅放宽表位先验 SD
  "sens1_0.15_0.15" = c(main = 0.15, epi = 0.15),  # 主抗体与表位均放宽
  "sens1_0.10_0.10" = c(main = 0.10, epi = 0.10),  # 主抗体与表位均放宽
  "sens1_0.10_0.15" = c(main = 0.10, epi = 0.15),  # 主抗体与表位均放宽
  "sens1_0.15_0.10" = c(main = 0.15, epi = 0.10)  # 主抗体与表位均放宽
)

# 动态生成 stan_data 集合 (保存在列表中)
stan_data_list <- lapply(sens_scenarios, function(pars) {
  d <- stan_data_base
  d$prior_sd_main <- pars["main"]
  d$prior_sd_epi  <- pars["epi"]
  return(d)
})

# 此时直接通过命名提取：
# stan_data_list[["sens1_0.05_0.05"]]

# --- 3. Initial Values 函数 ----
my_inits <- function() {
  list(
    beta = rep(0.5, ncol(B_basis)),
    tau = 0.5,
    pi = rep(0.85, 7),
    peak = c(500, 5000, 1.5, 1.5, 1.5, 1.5, 1.5),
    Mat0 = 150,
    Mat0_k = rep(10, 5),
    sigma1 = 1.5,
    sigma2 = 1.5,
    phi = rep(5, 5),
    B_IgA = runif(1, 4, 6),
    B_IgG = runif(1, 40, 60),
    B_epi = runif(5, 0.03, 0.07)
  )
}

# --- 4. 运行单个模型或批量循环采样 ----
options(mc.cores = parallel::detectCores())

# 示例：运行基线场景 sens1_0.05_0.05
# fit_truncat_extend_sens1_0.05_0.05 <- stan(
#   file = "stan/fit_truncat_extend_sens1.stan",
#   data = stan_data_list[["sens1_0.05_0.05"]],
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   init = my_inits,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
#
# saveRDS(fit_truncat_extend_sens1_0.05_0.05,'rda/fit_truncat_extend_sens1_0.05_0.05.rds')



