dir.create('code_04/plot', recursive = TRUE, showWarnings = FALSE)

library(posterior)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(grid)

model_names <- c(
  "fit_main",
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
  "fit_sens1_2"
)

model_names %>%
  walk(~{
    MODEL_NAME <- "fit_main"
    MODEL_NAME <- .x
    fit_main <- readRDS(sprintf("rda/%s.rds", MODEL_NAME))

    # ============================================================
    # 1. MCMC diagnostics
    # ============================================================

    draws <- as_draws_array(fit_main)
    diag <- summarise_draws(draws, rhat, ess_bulk, ess_tail) %>% filter(variable != "lp__")

    # ============================================================
    # 2. Parameter groups
    # ============================================================

    diag <- diag %>% mutate(group = case_when(
      variable %in% c("B_IgA", "B_IgG") | str_detect(variable, "^B_epi\\[") ~ "Background antibody level",
      str_detect(variable, "^beta\\[") ~ "B-spline coefficients",
      variable == "Mat0" | str_detect(variable, "^Mat0_k\\[") ~ "Maternally transferred antibody at birth",
      str_detect(variable, "^peak\\[") ~ "Post-infection antibody peak increment",
      str_detect(variable, "^sigma") | str_detect(variable, "^phi\\[") ~ "Observation model variability",
      str_detect(variable, "^pi\\[") ~ "Marker sensitivity",
      variable == "tau" ~ "Spline smoothing parameter",
      str_detect(variable, "^w\\[") ~ "Age-specific relative susceptibility"
    ))

    # ============================================================
    # 3. Colours
    # ============================================================

    cols <- c(
      "Background antibody level" = "#2878B5",
      "B-spline coefficients" = "#2CA02C",
      "Maternally transferred antibody at birth" = "#E64B35",
      "Post-infection antibody peak increment" = "#7E57C2",
      "Observation model variability" = "#D97706",
      "Marker sensitivity" = "#159D93",
      "Spline smoothing parameter" = "#666666",
      "Age-specific relative susceptibility" = "#174A9C"
    )

    fills <- c(
      "Background antibody level" = "#F1F7FB",
      "B-spline coefficients" = "#F2F8F2",
      "Maternally transferred antibody at birth" = "#FDF3F1",
      "Post-infection antibody peak increment" = "#F6F2FB",
      "Observation model variability" = "#FDF6ED",
      "Marker sensitivity" = "#EFF8F7",
      "Spline smoothing parameter" = "#F4F4F4",
      "Age-specific relative susceptibility" = "#F0F4FA"
    )

    # ============================================================
    # 4. Left panel: 38 parameters
    # ============================================================

    order_left <- c("B_IgA", "B_IgG", paste0("B_epi[", 1:5, "]"), paste0("beta[", 1:10, "]"), "Mat0", paste0("Mat0_k[", 1:5, "]"), paste0("peak[", 1:7, "]"), "sigma1", "sigma2", paste0("phi[", 1:5, "]"), "tau")

    left <- diag %>% filter(variable %in% order_left) %>% mutate(variable = factor(variable, levels = order_left)) %>% arrange(variable) %>% mutate(y = rev(seq_len(n())))

    left_groups <- left %>% group_by(group) %>% summarise(ymin = min(y), ymax = max(y), n_par = n(), .groups = "drop") %>% mutate(label_y = ifelse(n_par == 1, ymax + 0.38, ymax - 0.08))

    left_sep <- left_groups$ymin[left_groups$ymin != min(left_groups$ymin)] - 0.5

    # ============================================================
    # 5. Right panel: 37 parameters
    # ============================================================

    order_right <- c(paste0("pi[", 1:7, "]"), paste0("w[", 1:30, "]"))

    right <- diag %>% filter(variable %in% order_right) %>% mutate(variable = factor(variable, levels = order_right)) %>% arrange(variable) %>% mutate(y = rev(seq_len(n())))

    right_groups <- right %>% group_by(group) %>% summarise(ymin = min(y), ymax = max(y), n_par = n(), .groups = "drop") %>% mutate(label_y = ymax - 0.08)

    right_sep <- right_groups$ymin[right_groups$ymin != min(right_groups$ymin)] - 0.5

    # ============================================================
    # 6. Common axis
    # ============================================================

    xmin <- min(diag$rhat, na.rm = TRUE) - 0.0003
    xmax <- 1.0106
    xbreaks <- c(0.999, 1.000, 1.002, 1.004, 1.006, 1.008, 1.010)
    label_x <- 1.00985

    # ============================================================
    # 7. Left plot
    # ============================================================

    p_left <- ggplot(left, aes(x = rhat, y = y, colour = group)) +
      geom_hline(yintercept = left$y, colour = "grey94", linewidth = 0.22) +
      geom_hline(yintercept = left_sep, colour = "grey72", linewidth = 0.40) +
      geom_vline(xintercept = 1.01, linetype = "dashed", colour = "#D73027", linewidth = 0.50) +
      geom_point(size = 1.8) +
      geom_label(data = left_groups, aes(x = label_x, y = label_y, label = group, colour = group, fill = group), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.05, fontface = "bold", label.padding = unit(0.09, "lines"), label.r = unit(0.05, "lines"), linewidth = 0) +
      scale_colour_manual(values = cols, guide = "none") +
      scale_fill_manual(values = fills, guide = "none") +
      scale_y_continuous(breaks = left$y, labels = as.character(left$variable), expand = expansion(add = c(0.5, 0.5))) +
      scale_x_continuous(breaks = xbreaks, labels = function(x) sprintf("%.3f", x), limits = c(xmin, xmax), expand = c(0, 0)) +
      labs(x = expression(hat(R)), y = NULL) +
      theme_classic(base_size = 9) +
      theme(axis.text.y = element_text(size = 7, colour = "black", margin = margin(r = 4)), axis.text.x = element_text(size = 7.5, colour = "black"), axis.title.x = element_text(size = 10, margin = margin(t = 2)), axis.ticks.y = element_blank(), panel.grid = element_blank(), plot.margin = margin(2, 8, 0, 5))

    # ============================================================
    # 8. Right plot
    # ============================================================

    p_right <- ggplot(right, aes(x = rhat, y = y, colour = group)) +
      geom_hline(yintercept = right$y, colour = "grey94", linewidth = 0.22) +
      geom_hline(yintercept = right_sep, colour = "grey72", linewidth = 0.40) +
      geom_vline(xintercept = 1.01, linetype = "dashed", colour = "#D73027", linewidth = 0.50) +
      geom_point(size = 1.8) +
      geom_label(data = right_groups, aes(x = label_x, y = label_y, label = group, colour = group, fill = group), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.05, fontface = "bold", label.padding = unit(0.09, "lines"), label.r = unit(0.05, "lines"), linewidth = 0) +
      scale_colour_manual(values = cols, guide = "none") +
      scale_fill_manual(values = fills, guide = "none") +
      scale_y_continuous(breaks = right$y, labels = as.character(right$variable), expand = expansion(add = c(0.5, 0.5))) +
      scale_x_continuous(breaks = xbreaks, labels = function(x) sprintf("%.3f", x), limits = c(xmin, xmax), expand = c(0, 0)) +
      labs(x = expression(hat(R)), y = NULL) +
      theme_classic(base_size = 9) +
      theme(axis.text.y = element_text(size = 7, colour = "black", margin = margin(r = 4)), axis.text.x = element_text(size = 7.5, colour = "black"), axis.title.x = element_text(size = 10, margin = margin(t = 2)), axis.ticks.y = element_blank(), panel.grid = element_blank(), plot.margin = margin(2, 5, 0, 8))

    # ============================================================
    # 9. Combine
    # ============================================================

    p_final <- p_left + p_right +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 0))))

    p_final

    # ============================================================
    # 10. Export
    # ============================================================

    ggsave(sprintf("code_04/plot/Rhat_convergence_%s.tiff",.x), p_final, device = "tiff", width = 12, height = 6, units = "in", dpi = 600, compression = "lzw")

  })

