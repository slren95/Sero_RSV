# sens3 公共数据准备与初值脚本；由 3、9、12 月窗口入口调用，不执行拟合或保存。
stopifnot(exists("window_size_val"), exists("MAX_AGE_val"))
#rm(list=ls())
library(ggsci)
library(rstan)
library(splines)
library(tidyverse)
library(bayesplot)
library(posterior)
library(rio)

# Main
# Model Detail ---
# Data: df_sero.4
# Age range: 0-29
# 新版主模型：r_mat、r1、r2、r_epi 作为带对数正态先验的待估参数。
# 从项目根目录运行；输入读取现有 rda/，拟合结果保存到 rda_new/。

# --- 数据处理 ---
load('rda/code_01_data_explore.RData')
epi_prior_stan<-import('rda/epi_prior_stan.rds',trust=T)

# 与 Y_epi 的列顺序一致：III、II、V、Phi、IV。
epi_order <- c("t12_k_site_iii_logit", "t12_k_site_ii_logit",
               "t12_k_site_v_logit", "t12_k_site_ph_logit",
               "t12_k_site_iv_logit")
stopifnot(all(c("site", "meanlog", "sdlog") %in% names(epi_prior_stan)),
          !anyDuplicated(epi_prior_stan$site),
          all(epi_order %in% epi_prior_stan$site))
epi_prior_stan <- epi_prior_stan[match(epi_order, epi_prior_stan$site), ]
stopifnot(all(is.finite(epi_prior_stan$meanlog)),
          all(is.finite(epi_prior_stan$sdlog)),
          all(epi_prior_stan$sdlog > 0))

df_sero_sub <- df_sero.4 %>%
  mutate(DOB = as.Date(DOB),
         DOC = as.Date(DOC),
         DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1) %>%
  filter(age_value < MAX_AGE_val)


df_pos.2 <- df_pos.2 %>% arrange(MONTH)
lambda_vec <- df_pos.2$pos_scale
M_total <- length(lambda_vec)

# --- B-Spline 与 索引矩阵 ---
knots_internal <- c(4, 8, 12, 16, 20, 24)
B_basis <- bs(0:(MAX_AGE_val - 1), knots = knots_internal, degree = 3,
              Boundary.knots = c(0, MAX_AGE_val - 1), intercept = TRUE)

N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = MAX_AGE_val)

for (i in 1:N_current) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:(MAX_AGE_val - 1))
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1


# --- 构造 Stan Data ----
stan_data <- list(
  MAX_AGE = MAX_AGE_val,
  window_size = window_size_val,
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

  prior_r_mat = 20.8/36,
  prior_r1 = 20.8/(34/2),
  prior_r2 = 20.8/34,
  prior_meanlog_epi = epi_prior_stan$meanlog,
  prior_sdlog_epi = epi_prior_stan$sdlog,

  J = ncol(B_basis),
  B_basis = B_basis
)

# --- Run Bayes Model ----
options(mc.cores = 4)

my_inits <- function() {
  list(
    beta = rep(0.5, ncol(B_basis)),
    tau = 0.5,
    r_mat = stan_data$prior_r_mat,
    r1 = stan_data$prior_r1,
    r2 = stan_data$prior_r2,
    r_epi = exp(stan_data$prior_meanlog_epi),
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
