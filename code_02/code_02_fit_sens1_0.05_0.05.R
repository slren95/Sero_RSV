source('code_02/code_02_fit_sens1_base.R')

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
