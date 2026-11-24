library(tidyverse)
library(lubridate)
library(GGally)
library(patchwork)
library(scales)
library(grid)

font_theme <- theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(size = 9, face = "plain"),
  axis.text = element_text(size = 7, color = "black"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  strip.text = element_text(size = 8)
)

# 1. 数据准备 ----
load("rda/code_01_data_explore.RData")

df_sero_sub <- df_sero.4 %>%
  mutate(
    DOB = as.Date(DOB), DOC = as.Date(DOC),
    DOM_index = interval(ymd("2020-06-01"), DOB) %/% months(1) + 1,
    iga_log = log(`IgA-MFI`), igg_log = log(`IgG-MFI`)
  ) %>%
  filter(age_value < 30)

site_order <- c(
  "iga_log", "igg_log",
  "Site II-Motavizumab", "SiteIII-MPE8",
  "Site IV-Clesrovimab", "Site V-Suptavumab",
  "Site Φ-Nirsevimab"
)

df_plot <- df_sero_sub %>% select(all_of(site_order))
variable_labels <- c("ln(IgA MFI)", "ln(IgG MFI)", "Site II", "Site III",
                     "Site IV", "Site V", "Site Φ")

# 2. 坐标轴配置 ----
# Site inhibition 为 0–1，显示为 0%–100%
axis_config <- list(
  iga_log = list(limits = c(0, 12), breaks = seq(0, 12, 4), labels = label_number()),
  igg_log = list(limits = c(0, 12), breaks = seq(0, 12, 4), labels = label_number()),
  `Site II-Motavizumab` = list(limits = c(0, 1), breaks = seq(0, 1, 0.5), labels = label_percent()),
  `SiteIII-MPE8` = list(limits = c(0, 1), breaks = seq(0, 1, 0.5), labels = label_percent()),
  `Site IV-Clesrovimab` = list(limits = c(0, 1), breaks = seq(0, 1, 0.5), labels = label_percent()),
  `Site V-Suptavumab` = list(limits = c(0, 1), breaks = seq(0, 1, 0.5), labels = label_percent()),
  `Site Φ-Nirsevimab` = list(limits = c(0, 1), breaks = seq(0, 1, 0.5), labels = label_percent())
)

add_axis_config <- function(p, vx, vy = NULL) {
  cx <- axis_config[[vx]]
  cy <- if (!is.null(vy)) axis_config[[vy]] else NULL

  if (!is.null(cx)) p <- p + scale_x_continuous(breaks = cx$breaks, labels = cx$labels)
  if (!is.null(cy)) p <- p + scale_y_continuous(breaks = cy$breaks, labels = cy$labels)

  if (!is.null(cx)) {
    p <- p + scale_x_continuous(
      breaks = cx$breaks,
      labels = cx$labels,
      expand = expansion(mult = 0.06)
    )
  }

  if (!is.null(cy)) {
    p <- p + scale_y_continuous(
      breaks = cy$breaks,
      labels = cy$labels,
      expand = expansion(mult = 0.06)
    )
  }
  p + coord_cartesian(
    xlim = if (!is.null(cx)) cx$limits else NULL,
    ylim = if (!is.null(cy)) cy$limits else NULL,
    expand = TRUE
  )
}

# 3. Upper：Spearman 相关系数 + 显著性 ----
panel_cor_tile <- function(data, mapping, ...) {
  x <- rlang::eval_tidy(mapping$x, data)
  y <- rlang::eval_tidy(mapping$y, data)
  valid <- is.finite(x) & is.finite(y)

  if (sum(valid) >= 3 && length(unique(x[valid])) > 1 && length(unique(y[valid])) > 1) {
    test <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman", exact = FALSE))
    r <- unname(test$estimate)
    pval <- test$p.value
  } else {
    r <- NA_real_
    pval <- NA_real_
  }

  stars <- case_when(
    is.na(pval) ~ "",
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
  )

  lab <- if (is.na(r)) "NA" else paste0(sprintf("%.2f", r), stars)

  ggplot(tibble(x = 1, y = 1, r = r, label = lab), aes(x, y)) +
    geom_tile(aes(fill = r), width = 1, height = 1) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_gradient2(low = "#6D9EC1", mid = "white", high = "salmon",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90", guide = "none") +
    # scale_fill_gradient2(
    #   low = "#547DA5", mid = "white", high = "#C76570",
    #   midpoint = 0, limits = c(-1, 1),
    #   na.value = "grey90", guide = "none"
    # )+
    coord_cartesian(xlim = c(0.5, 1.5), ylim = c(0.5, 1.5), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3))
}

# 4. Diagonal：核密度 + Rug ----
panel_density <- function(data, mapping, ...) {
  vx <- rlang::as_name(mapping$x)

  p <- ggplot(data, mapping) +
    geom_density(fill = "#9BB8D0", alpha = 0.7, color = "#365F80",
                 linewidth = 0.4, na.rm = TRUE) +
    geom_rug(sides = "b", alpha = 0.18, linewidth = 0.2, na.rm = TRUE) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  add_axis_config(p, vx)
}

# 5. Lower：散点图 + LOESS ----
panel_scatter <- function(data, mapping, ...) {
  vx <- rlang::as_name(mapping$x)
  vy <- rlang::as_name(mapping$y)

  p <- ggplot(data, mapping) +
    geom_point(alpha = 0.18, size = 0.55, color = "grey35", na.rm = TRUE) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                color = "#C44E52", linewidth = 0.6, na.rm = TRUE) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank())

  add_axis_config(p, vx, vy)
}

# 6. ggpairs 主图 ----
p <- ggpairs(
  df_plot,
  upper = list(continuous = panel_cor_tile),
  diag = list(continuous = panel_density),
  lower = list(continuous = panel_scatter),
  columnLabels = variable_labels,
  axisLabels = "show",
  progress = FALSE
) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_blank(),
    #strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7, color = "grey30"),
    axis.text.x = element_text(margin = margin(t = 3)),
    axis.text.y = element_text(margin = margin(r = 3)),
    axis.ticks.length = unit(1.5, "pt")
  ) + font_theme

p[1, 1] <- p[1, 1] +
  scale_y_continuous(
    breaks = c(0, 0.1, 0.2),
    labels = c(0,0.1,0.2)
  )

print(p)



# 7. Spearman 相关系数颜色图例（-1 至 1）----
legend_data <- tibble(x = seq(-1, 1, length.out = 201), y = 1)

p_legend <- ggplot(legend_data, aes(x, y, fill = x)) +
  geom_tile(width = 0.011, height = 0.4) +
  scale_fill_gradient2(low = "#6D9EC1", mid = "white", high = "#E46726",
                       midpoint = 0, limits = c(-1, 1), guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0, 0)) +
  labs(x = "Spearman's rho", y = NULL) +
  theme_void(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9, margin = margin(t = 2)),
    axis.ticks.x = element_line(color = "grey40"),
    axis.ticks.length.x = unit(2, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  )

# 8. 拼接主图与顶部 Legend ----
p_grob <- grid.grabExpr(print(p), wrap = TRUE)

p_legend_top <- plot_spacer() + p_legend + plot_spacer() +
  plot_layout(widths = c(1.8, 1, 1.8))

pairs_final <- p_legend_top / wrap_elements(full = p_grob) +
  plot_layout(heights = c(0.7, 18))

print(pairs_final)

# 9. 导出 ----

ggsave("plot/serology_correlation_matrix.tiff", p, width = 12, height = 10, dpi = 300, bg = "white")


