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
# Waning of epitope: r_mat、r1、r2 和 r_epi 并不是死死固定死的常数，而是被当作待估计的未知参数（Parameters），并在模型中为其赋予了带有均值和标准差的先验分布（即你所说的层次先验/随机先验形式）。
# epi_prior_stan

# --- 数据处理 ---
load('rda/code_01_data_explore.RData')
epi_prior_stan<-import('rda/epi_prior_stan.rds',trust=T)

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


# --- 构造 Stan Data ----
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

  prior_r_mat = 20.8/36,
  prior_r1 = 20.8/(34/2),
  prior_r2 = 20.8/34,
  prior_meanlog_epi = epi_prior_stan$meanlog,
  prior_sdlog_epi = epi_prior_stan$sdlog,

  J = ncol(B_basis),
  B_basis = B_basis
)

# --- Run Bayes Model ----
options(mc.cores = parallel::detectCores())

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

# fit_sens1_2 <- stan(
#   file = "stan/fit_sens1_2.stan",
#   data = stan_data,
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   init = my_inits,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
#
# saveRDS(fit_sens1_2,'rda/fit_sens1_2.rds')



