# Prepare age, calendar-time and observed-serology data for code_03_new.
# Run from the project root.
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)
library(rio)

load("rda/code_01_data_explore.RData")

df_sero_sub <- df_sero.4 %>%
  mutate(
    DOB = as.Date(DOB),
    DOC = as.Date(DOC),
    DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1
  ) %>%
  filter(age_value < 30)

df_pos.2 <- df_pos.2 %>% arrange(MONTH)
lambda_vec <- df_pos.2$pos_scale
M_total <- length(lambda_vec)
N_current <- nrow(df_sero_sub)
cal_idx_matrix <- matrix(1, nrow = N_current, ncol = 30)
for (i in seq_len(N_current)) {
  cal_idx_matrix[i, ] <- df_sero_sub$DOM_index[i] + (0:29)
}
cal_idx_matrix[cal_idx_matrix > M_total] <- 1
