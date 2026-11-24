# 从项目根目录运行：
# source("code_02_new/run_all_model.R", encoding = "UTF-8")
# 依次拟合 model_names 中启用的模型；每个模型 4 条链、4 核，结果保存到 rda_new/。
# 任一模型报错时停止；重新运行本脚本会重新拟合并覆盖对应 RDS。

local({
  model_names <- c(
    #"fit_main_new",            # 主分析
    #"fit_sen1.1_0.05_new",     # 前三个衰减率：sdlog = 0.05
    #"fit_sen1.1_0.15_new",     # 前三个衰减率：sdlog = 0.15；新增场景，取消注释可启用
    #"fit_sen1.1_0.2_new",      # 前三个衰减率：sdlog = 0.20
    #"fit_sens1.2_0.5_new",     # 五个表位衰减率：sdlog = 0.50
    #"fit_sens1.2_0.2_new",     # 五个表位衰减率：sdlog = 0.20
    #"fit_sens2_new",           # 样条内部节点
    #"fit_sens3_3_new",         # 3 个月追溯窗口
    #"fit_sens3_9_new",         # 9 个月追溯窗口
    "fit_sens3_12_new",        # 12 个月追溯窗口
    #"fit_sens4_1_1_new",       # pi ~ beta(1, 1)
    "fit_sens4_5_5_new"       # pi ~ beta(5, 5)
    #"fit_sens5_new"            # 仅 IgA 和 IgG
  )

  scripts <- file.path("code_02_new", paste0("code_02_", model_names, ".R"))
  if (!all(file.exists(scripts))) {
    stop("请从项目根目录运行；缺少脚本：",
         paste(scripts[!file.exists(scripts)], collapse = ", "))
  }

  for (i in seq_along(scripts)) {
    started_at <- Sys.time()
    message(sprintf("[%d/%d] 开始 %s：%s",
                    i, length(scripts), model_names[i], format(started_at)))

    # 每个模型独立准备数据；保存完成后释放对象，避免累计占用内存。
    model_env <- new.env(parent = globalenv())
    source(scripts[i], local = model_env, encoding = "UTF-8")
    rm(model_env)
    invisible(gc())

    message(sprintf("[%d/%d] 完成 %s，耗时 %.1f 分钟；结果：rda_new/%s.rds",
                    i, length(scripts), model_names[i],
                    as.numeric(difftime(Sys.time(), started_at, units = "mins")),
                    model_names[i]))
  }

  message(sprintf("所选 %d 个模型拟合及保存完成。", length(model_names)))
})
