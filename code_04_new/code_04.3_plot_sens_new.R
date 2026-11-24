library(tidyverse)
library(patchwork)
library(rio)

color_palette <- c(
  # 01–10: Tableau 10 (original palette)
  "#1F77B4",  # Blue
  "#FF7F0E",  # Orange
  "#2CA02C",  # Green
  "#D62728",  # Red
  "#9467BD",  # Purple
  "#8C564B",  # Brown
  "#E377C2",  # Pink
  "#7F7F7F",  # Gray
  "#BCBD22",  # Olive
  "#17BECF",  # Cyan

  # 11–20: Extended scientific palette
  "#005F73",  # Deep teal
  "#E76F51",  # Coral
  "#6A4C93",  # Deep violet
  "#219EBC",  # Ocean blue
  "#CA6702",  # Burnt orange
  "#40916C",  # Forest green
  "#B56576",  # Dusty rose
  "#3A86FF",  # Bright blue
  "#9B5DE5",  # Orchid purple
  "#F4A261",  # Apricot

  # 21–30: Additional contrasting colors
  "#264653",  # Dark slate
  "#F72585",  # Magenta
  "#4D908E",  # Muted teal
  "#F8961E",  # Amber
  "#577590",  # Steel blue
  "#90BE6D",  # Light green
  "#B5179E",  # Plum
  "#277DA1",  # Azure
  "#A44A3F",  # Brick red
  "#43AA8B",  # Sea green

  # 31–40: Further extensions
  "#7209B7",  # Royal purple
  "#F3722C",  # Tangerine
  "#184E77",  # Navy blue
  "#D45087",  # Raspberry
  "#588157",  # Sage green
  "#C77DFF",  # Lavender
  "#9A6324",  # Ochre brown
  "#00A6A6",  # Turquoise
  "#B23A48",  # Crimson
  "#4361EE"   # Royal blue
)

dir.create('code_04_new/plot', recursive = TRUE, showWarnings = FALSE)

model_names <- c(
  "fit_main_new",
  "fit_sen1.1_0.05_new",
  "fit_sen1.1_0.15_new",
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


# 四类图共用分组；仅保留 model_names 中启用的模型，自动跳过空组。
# 后续恢复模型时，只需取消上面 model_names 中相应行的注释。
sens_groups <- list(
  #c("fit_sen1.1_0.05_new"),
  #c("fit_sens1.2_0.5_new", "fit_sens1.2_0.2_new"),
  "fit_sens2_new",
  c("fit_sens3_3_new", "fit_sens3_9_new"),
  c("fit_sens4_5_5_new"),
  "fit_sens5_new"
) %>%
  map(~ intersect(.x, model_names)) %>%
  keep(~ length(.x) > 0)

color_map <- setNames(
  color_palette[seq_along(model_names)],
  model_names
)

label_map <- c(
  "fit_main_new" = "Main analysis",
  "fit_sen1.1_0.05_new" = "Sensitivity 1.1 (0.05)",
  "fit_sen1.1_0.15_new" = "Sensitivity 1.1 (0.15)",
  "fit_sens1.2_0.5_new" = "Sensitivity 1.2 (0.5)",
  "fit_sens1.2_0.2_new" = "Sensitivity 1.2 (0.2)",
  "fit_sens2_new" = "Sensitivity analysis 1",
  "fit_sens3_3_new" = "Sensitivity analysis 2 (3 months)",
  "fit_sens3_9_new" = "Sensitivity analysis 2 (9 months)",
  "fit_sens3_12_new" = "Sensitivity analysis 2 (12 months)",
  "fit_sens4_1_1_new" = "Sensitivity analysis 3 (1, 1)",
  "fit_sens4_5_5_new" = "Sensitivity analysis 3",
  "fit_sens5_new" = "Sensitivity analysis 4"
)


# 易感性曲线 w(a) ----

df_w<-map_dfr(model_names,~{
  fit_truncat <- readRDS(file.path('rda_new', paste0(.x, '.rds')))
  posterior_w_truncat <- rstan::extract(fit_truncat, "w")$w
  w_res_truncat <- data.frame(
    age = 0:(ncol(posterior_w_truncat)-1),
    mean = colMeans(posterior_w_truncat),
    lower = apply(posterior_w_truncat, 2, quantile, 0.025),
    upper = apply(posterior_w_truncat, 2, quantile, 0.975),
    model = .x
  )
}
)

df_w %>%
  filter(str_detect(model,'(main|sen1)')) %>%
  ggplot(aes(x = age, y = mean, color = model, fill = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x = "Age (months)",
    y = "Relative susceptibility for RSV infection",
    color = "Model", fill = "Model"
  ) +
  scale_x_continuous(
    breaks = seq(0, 23, by = 1),
    expand = expansion(add = c(0, .5)),
    labels = function(x) sprintf('%s-<%s', x, x + 1)
  ) +
  scale_y_continuous(expand = expansion()) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(color = "gray90", fill = NA, linewidth = .5),
    axis.line = element_line(linewidth = .2, arrow = arrow(length = unit(0, 'mm'))),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = .5, size = 6)
  ) +
  coord_cartesian(xlim = c(0, 23))



# 绘图函数：ggplot 只写一次，支持多个 sens 模型
plot_w <- function(dat, main_model = "fit_main_new", sens_models) {
  dat_sub <- dat %>% filter(model %in% c(main_model, sens_models)) %>%
    filter(model != "fit_sens3_3_new" | age >= 2)

  ggplot(dat_sub, aes(x = age, y = mean, color = model, fill = model)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = color_map, labels = label_map, name = NULL) +
    scale_fill_manual(values = color_map, labels = label_map, name = NULL) +
    labs(
      x = "Age (months)",
      y = "Relative susceptibility for RSV infection",
      color = "Model", fill = "Model"
    ) +
    scale_x_continuous(
      limits = c(0,23),
      breaks = seq(0, 23, by = 1),
      expand = expansion(add = c(0.5, .5)),
      labels = function(x) sprintf('%s-<%s', x, x + 1)
    ) +
    scale_y_continuous(limits = c(0,3)) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linewidth = .2),
      axis.line = element_line(linewidth = .2, arrow = arrow(length = unit(0, "mm"))),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = .5, size = 6),
      legend.position = 'top'
    ) +
    coord_cartesian(xlim = c(0, 23))
}

plots_w <- sens_groups %>%
  map(~ plot_w(df_w, main_model = "fit_main_new", sens_models = .x))

wrap_plots(plots_w, nrow = 2, axis_titles = 'collect')

ggsave("code_04_new/plot/Sens_w.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")


# 灵敏度图 ----
df_s<-map_dfr(model_names,~{
  fit_truncat <- readRDS(file.path('rda_new', paste0(.x, '.rds')))
  pi_draws <- rstan::extract(fit_truncat, pars = "pi")$pi
  if(ncol(pi_draws)==2){
    pi_draws <- cbind(pi_draws, matrix(NA, nrow = nrow(pi_draws), ncol = 5))
  }
  df_markers <- tibble(
    marker_idx = 1:7,
    # 对应标书文本：注入具体的抗体与标记物生物学标签
    biological_meaning = c("IgA","IgG","SiteIII-MPE8", "Site II-Motavizumab",
                           "Site V-Suptavumab", "Site Φ-Nirsevimab",
                           "Site IV-Clesrovimab"),
    # 💡 核心对账步：直接从 Stan 后验样本矩阵中切片计算分位数
    median = apply(pi_draws, 2, quantile, 0.50, na.rm = TRUE),  # 后验中位数
    lower  = apply(pi_draws, 2, quantile, 0.025, na.rm = TRUE), # 95% 信用区间下界
    upper  = apply(pi_draws, 2, quantile, 0.975, na.rm = TRUE)  # 95% 信用区间上界
  ) %>%
    mutate(biological_meaning = factor(biological_meaning, levels = rev(biological_meaning)),
           model=.x)
})

shape_map <- setNames(ifelse(model_names == "fit_main_new", 10, 16), model_names)

plot_s <- function(dat, main_model = "fit_main_new", sens_models,
                   dodge_width = 0.6, show_label = TRUE) {
  dat_sub <- dat %>%
    filter(model %in% c(main_model, sens_models)) %>%
    mutate(model = factor(model, levels = c(main_model, sens_models)))

  pos <- position_dodge(width = dodge_width)
  lv  <- levels(dat_sub$model)              # 当前图里实际出现的模型

  p <- ggplot(dat_sub, aes(x = median, y = biological_meaning,
                           color = model, fill = model)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0.15, linewidth = 0.9,
                   position = pos) +
    geom_point(aes(shape = model), size = 2, position = pos) +
    scale_color_manual(values = color_map,labels = label_map, name = NULL) +
    scale_fill_manual(values  = color_map,labels = label_map, name = NULL) +
    scale_shape_manual(values = shape_map,labels = label_map, name = NULL) +
    guides(
      color = guide_legend(override.aes = list(
        fill  = unname(color_map[lv]),      # ← 只取当前 levels
        shape = unname(shape_map[lv])
      )),
      fill  = "none",
      shape = "none"
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0.7, 1.02)) +
    labs(x = 'Sensitivity to RSV infection', y = NULL) +
    theme_bw()+
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
          panel.grid.minor = element_blank(),
          panel.grid = element_line(linewidth = .2),
          axis.line = element_line(linewidth = .2, arrow = arrow(length = unit(0, "mm"))),
          axis.title = element_text(face = "bold", size = 10),
          axis.text = element_text(color = "black"),
          legend.position = 'top')

  if (show_label) {
    p <- p + geom_text(
      aes(label = sprintf("%.1f%% (%.1f%%-%.1f%%)",
                          median * 100, lower * 100, upper * 100)),
      position = pos, hjust = 0.5, vjust = -1.0,
      size = 2.6, color = "black", show.legend = FALSE
    )
  }
  p
}


# 4. 配对并 map
plots_s <- sens_groups %>%
  map(~ plot_s(df_s, main_model = "fit_main_new", sens_models = .x, show_label = FALSE))

wrap_plots(plots_s, nrow = 2, axis_titles = 'collect')

ggsave("code_04_new/plot/Sens_sens.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")

# 感染概率 ----
# 计算并保存每个模型的感染概率汇总结果 --

if (F) {
  walk(model_names, ~{
    MODEL_NAME <- .x
    message("Processing: ", MODEL_NAME)

    df_individual_wide_all <- import(file.path('code_03_new',MODEL_NAME,'df_individual_wide_all.rds'), trust = TRUE)

    df_c_model <- df_individual_wide_all %>%
      dplyr::select(post_id, ID, starts_with("prob_infect_age_")) %>%
      mutate(across(starts_with("prob_infect_age_"), ~ na_if(.x, 0))) %>%
      pivot_longer(
        cols = starts_with("prob_infect_age_"),
        names_to = "age",
        values_to = "prob"
      ) %>%
      mutate(age = as.integer(str_remove(age, "prob_infect_age_"))) %>%
      summarise(mean_prob = mean(prob, na.rm = TRUE), .by = c(post_id, age)) %>%
      arrange(post_id, age) %>%
      mutate(cum_prob = 1 - cumprod(1 - mean_prob), .by = post_id) %>%
      summarise(
        across(
          c(mean_prob, cum_prob),
          list(
            est = ~ quantile(.x, 0.5, na.rm = TRUE),
            lci = ~ quantile(.x, 0.025, na.rm = TRUE),
            uci = ~ quantile(.x, 0.975, na.rm = TRUE)
          )
        ),
        .by = age
      ) %>%
      mutate(model = MODEL_NAME)

    saveRDS(df_c_model, file.path('rda_new', sprintf("df_c_%s.rds", MODEL_NAME)))
  })
}


df_c <- map_dfr(model_names, ~ import(file.path("rda_new", sprintf("df_c_%s.rds", .x)), trust = TRUE))

# 累计感染概率绘图函数 --
plot_c <- function(dat, main_model = "fit_main_new", sens_models, MAX_AGE = 23) {

  dat_sub <- dat %>%
    filter(model %in% c(main_model, sens_models)) %>%
    mutate(model = factor(model, levels = c(main_model, sens_models))) %>%
    filter(model != "fit_sens3_3_new" | age >= 2)

  ggplot(dat_sub, aes(x = age, y = cum_prob_est, color = model, fill = model)) +
    geom_ribbon(aes(ymin = cum_prob_lci, ymax = cum_prob_uci), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = color_map, labels = label_map, name = NULL) +
    scale_fill_manual(values = color_map, labels = label_map, name = NULL) +
    scale_x_continuous(
      limits = c(0,MAX_AGE),
      breaks = seq(0, MAX_AGE, by = 1),
      expand = expansion(add = c(0.5, .7)),
      labels = function(x) sprintf("%s-<%s", x, x + 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linewidth = .2),
      axis.line = element_line(linewidth = .2, arrow = arrow(length = unit(0, "mm"))),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = .5, size = 6),
      legend.position = 'top'
    ) +
    labs(
      x = "Age (months)",
      y = "Cumulative probability of RSV infection",
      color = "Model",
      fill = "Model"
    ) +
    coord_cartesian(xlim = c(0, MAX_AGE))
}


plots_c <- sens_groups %>%
  map(~ plot_c(df_c, main_model = "fit_main_new", sens_models = .x, MAX_AGE = 23))

wrap_plots(plots_c, nrow = 2, axis_titles = 'collect')

ggsave("code_04_new/plot/Sens_cum.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")

plot_i <- function(dat, main_model = "fit_main_new", sens_models, MAX_AGE = 23) {

  dat_sub <- dat %>%
    filter(model %in% c(main_model, sens_models)) %>%
    mutate(model = factor(model, levels = c(main_model, sens_models))) %>%
    filter(model != "fit_sens3_3_new" | age >= 2)

  ggplot(dat_sub, aes(x = age, y = mean_prob_est, color = model, fill = model)) +
    geom_ribbon(aes(ymin = mean_prob_lci, ymax = mean_prob_uci), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = color_map, labels = label_map, name = NULL) +
    scale_fill_manual(values = color_map, labels = label_map, name = NULL) +
    scale_x_continuous(
      limits = c(0,MAX_AGE),
      breaks = seq(0, MAX_AGE, by = 1),
      expand = expansion(add = c(0.5, .7)),
      labels = function(x) sprintf("%s-<%s", x, x + 1)
    ) +
    scale_y_continuous(
      limits = c(0, 0.6),
      breaks = seq(0, 0.6, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linewidth = .2),
      axis.line = element_line(linewidth = .2, arrow = arrow(length = unit(0, "mm"))),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = .5, size = 6),
      legend.position = 'top'
    ) +
    labs(
      x = "Age (months)",
      y = "Monthly probability of RSV infection",
      color = "Model",
      fill = "Model"
    ) +
    coord_cartesian(xlim = c(0, MAX_AGE))
}

plots_i <- sens_groups %>%
  map(~ plot_i(df_c, main_model = "fit_main_new", sens_models = .x, MAX_AGE = 23))

wrap_plots(plots_i, nrow = 2, axis_titles = 'collect')

ggsave("code_04_new/plot/Sens_ins.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")

# Sens 5 igA igA only ----
labels_alt <- label_map
labels_alt["fit_sens5_new"] <- "IgA and IgG only"

# 易感性
plots_w[[7]] <- plots_w[[4]] +
  scale_color_manual(values = color_map, labels = labels_alt, name = NULL) +
  scale_fill_manual(values = color_map, labels = labels_alt, name = NULL)

# 指标敏感度：还需同步修改 shape 图例
plots_s[[7]] <- plots_s[[4]] +
  scale_color_manual(values = color_map, labels = labels_alt, name = NULL) +
  scale_fill_manual(values = color_map, labels = labels_alt, name = NULL) +
  scale_shape_manual(values = shape_map, labels = labels_alt, name = NULL)

# 累计感染概率
plots_c[[7]] <- plots_c[[4]] +
  scale_color_manual(values = color_map, labels = labels_alt, name = NULL) +
  scale_fill_manual(values = color_map, labels = labels_alt, name = NULL)

# 月感染概率
plots_i[[7]] <- plots_i[[4]] +
  scale_color_manual(values = color_map, labels = labels_alt, name = NULL) +
  scale_fill_manual(values = color_map, labels = labels_alt, name = NULL)

