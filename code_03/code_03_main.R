rm(list=ls())
library(tidyverse)

model_names <- c(
  "fit_main",
  "fit_simple",
  "fit_sens3_3",
  "fit_sens3_9",
  "fit_sens3_12",
  "fit_sens2",
  "fit_sens1_0.05_0.05",
  "fit_sens1_0.05_0.15",
  "fit_sens1_0.15_0.15",
  "fit_sens1_0.10_0.10",
  "fit_sens1_0.10_0.15",
  "fit_sens1_0.15_0.10",
  "fit_sens4_5_5",
  "fit_sens4_1_1",
  "fit_sens5",
  "fit_main_new"
)

model_names[16] %>%
  walk(~{
    file_name1 <- paste0('stan/', .x, '.stan')
    file_name2 <-  file.path('code_02', paste0('code_02_', .x, '.R'))
    if(!file.exists(file_name1)) file.create(file_name1) else cat(file_name1,'已存在\n')
    if(!file.exists(file_name2)) file.create(file_name2) else cat(file_name2,'已存在\n')
    output_dir <- file.path('code_03', .x)
    if(!dir.exists(output_dir)) dir.create(output_dir) else cat(output_dir,'文件夹已存在\n')
  })

model_names[16] %>%
  purrr::walk(~{
    cat('precessing',.x,'\n\n')
    assign("MODEL_NAME", .x, envir = .GlobalEnv)
    source('code_03/code_03.1_stan_plot.R')
    source('code_03/code_03.2_infection_risk.R')
    source('code_03/code_03.3_dynamic.R')
  })
