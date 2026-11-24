# Run the existing code_03 post-processing workflow for saved new-model fits.
# Run this script from the project root; it does not fit or recompile Stan models.

model_names <- c(
  "fit_main_new",
  "fit_sen1.1_0.05_new",
  "fit_sen1.1_0.15_new",
  "fit_sen1.1_0.2_new",
  "fit_sens1.2_0.5_new",
  "fit_sens1.2_0.2_new",
  "fit_sens2_new",
  "fit_sens3_3_new",
  "fit_sens3_9_new",
  "fit_sens3_12_new",
  "fit_sens4_1_1_new",
  "fit_sens4_5_5_new",
  "fit_sens5_new"
)

for (i in seq_along(model_names)[3]) {
  model_name <- model_names[[i]]
  fit_path <- file.path("rda_new", paste0(model_name, ".rds"))

  if (!file.exists(fit_path)) {
    stop("找不到已保存的拟合结果：", fit_path, call. = FALSE)
  }

  message(sprintf("[%d/%d] 后处理 %s", i, length(model_names), model_name))
  model_env <- new.env(parent = globalenv())
  assign("MODEL_NAME", model_name, envir = model_env)

  source("code_03_new/code_03.1_stan_plot_new.R",
         local = model_env, encoding = "UTF-8")
  source("code_03_new/code_03.2_infection_risk_new.R",
         local = model_env, encoding = "UTF-8")

  dynamic_script <- if (identical(model_name, "fit_sens5_new")) {
    "code_03_new/code_03.3_dynamic_sens5_new.R"
  } else {
    "code_03_new/code_03.3_dynamic_new.R"
  }
  source(dynamic_script, local = model_env, encoding = "UTF-8")
}

message("全部新模型后处理完成。")
