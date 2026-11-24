# --- 💡 动态参数设置（在此处自由调整） ---
window_size_val <- 9  # 💡 可以自由改为 3, 6, 9, 12 等
MAX_AGE_val     <- 30 # 💡 对应 0-29 月龄 (共 30 个月龄单元)

source('code_02/code_02_fit_sens3_base.R')

# fit_truncat_extend_sens3_9 <- stan(
#   file = "stan/fit_truncat_extend_sens3_base.stan",
#   data = stan_data,
#   chains = 4,
#   iter = 2000,
#   warmup = 1000,
#   init = my_inits,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
#
# saveRDS(fit_truncat_extend_sens3_9, 'rda/fit_truncat_extend_sens3_9.rds')
