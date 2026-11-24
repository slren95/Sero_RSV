# sens3：3 个月追溯窗口；从项目根目录运行。
window_size_val <- 3
MAX_AGE_val <- 30

source("code_02_new/code_02_fit_sens3_base_new.R", local = TRUE, encoding = "UTF-8")

if (!dir.exists("rda_new")) dir.create("rda_new")

fit_sens3_3_new <- stan(
  file = "stan_new/fit_sens3_base_new.stan",
  data = stan_data,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  init = my_inits,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

saveRDS(fit_sens3_3_new, "rda_new/fit_sens3_3_new.rds")
