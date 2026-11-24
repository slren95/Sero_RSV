#rm(list=ls())
library(tidyverse)
library(boot)
library(posterior)

MODEL_NAME <- 'fit_sens5_new'
OUTPUT_DIR <- file.path('code_03_new',MODEL_NAME)
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

fit_truncat <- readRDS(file.path('rda_new', paste0(MODEL_NAME, '.rds')))

df_means <- enframe(get_posterior_mean(fit_truncat)[, 5])

## peak value ----
df_real_peaks <- df_means %>%
  pivot_wider(names_from = name, values_from = value) %>%
  pivot_longer(cols = starts_with("peak["), names_to = "peak_name", values_to = "inf_peak_delta") %>%
  mutate(
    idx = as.integer(str_extract(peak_name, "\\d+")),
    marker = c("IgA-MFI", "IgG-MFI")[idx]
  ) %>%
  rowwise() %>%
  mutate(
    inf_peak_value = case_when(
      idx == 1 ~ B_IgA + inf_peak_delta,
      idx == 2 ~ B_IgG + inf_peak_delta
    ),
    mat_peak_value = case_when(
      idx == 1 ~ B_IgA,
      idx == 2 ~ B_IgG + Mat0
    )
  ) %>%
  ungroup() %>%
  select(marker, B_IgA, B_IgG, Mat0, inf_peak_delta, inf_peak_value, mat_peak_value)

print(df_real_peaks)


# --- 1. 从数据框中动态提取所需的标量参数 ---
p_base <- df_real_peaks %>% slice(1)

b_iga   <- p_base$B_IgA
b_igg   <- p_base$B_IgG
mat0    <- p_base$Mat0

# 提取 2 个抗体的峰值增量 (1:IgA, 2:IgG)
peaks_vec <- df_real_peaks$inf_peak_delta

# 使用新版 fit_sens5 的后验衰减率均值
rate_draws <- rstan::extract(fit_truncat, pars=c("r_mat","r1","r2"))
r_mat <- mean(rate_draws$r_mat)
r1 <- mean(rate_draws$r1)
r2 <- mean(rate_draws$r2)


# --- 2. 构建模拟函数 (设定具体的感染满月龄 m_inf，例如 6 月龄) ---
get_kinetics_trajectory <- function(m_inf = 6) {
  # 生成从 0 到 23 月龄的高密度平滑时间点 (步长 0.05 个月)
  time_seq <- seq(0, 23, by = 0.05)

  map_df(time_seq, function(t) {
    # 生物学机制：m_inf 感染，m_inf + 1 达峰
    has_peaked <- t >= (m_inf + 1)
    dt <- if (has_peaked) t - m_inf - 1 else 0

    # 1. 计算 IgA 轨迹
    mu_IgA <- if (has_peaked) b_iga + peaks_vec[1] * exp(-r1 * dt) else b_iga

    # 2. 计算 IgG 轨迹 (母传 + 感染后突跃)
    mu_IgG_base <- b_igg + mat0 * exp(-r_mat * t)
    mu_IgG <- if (has_peaked) mu_IgG_base + peaks_vec[2] * exp(-r2 * dt) else mu_IgG_base

    # 返回当前时间点的结果 (仅包含 IgA 和 IgG)
    tibble(
      Age = t,
      `IgA-MFI` = mu_IgA,
      `IgG-MFI` = mu_IgG
    )
  })
}

# --- 3. 运行模拟并转换数据格式 ---
df_plot_data <- get_kinetics_trajectory(m_inf = 6) %>%
  pivot_longer(-Age, names_to = "Marker", values_to = "Value")

# --- 4. 优雅的 Apple-style 极简可视化 (双图分面) ---
ggplot(df_plot_data, aes(x = Age, y = Value, color = Marker)) +
  geom_vline(xintercept = 6, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 7, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  facet_wrap(~Marker, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = seq(0, 23, by = 2)) +
  labs(
    title = "Model-Inferred RSV Antibody Kinetics (IgA & IgG Only)",
    subtitle = "Simulated primary infection at 6 months (peaking at 7 months)",
    x = "Age (Months)",
    y = "Antibody Concentration (MFI)",
    caption = "Dashed line = Infection onset (6m) | Dotted line = Peak concentration (7m)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 11, colour = "grey20"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, margin = margin(b=4)),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b=15))
  )

ggsave(file.path(OUTPUT_DIR,'antibody_dynamic.tiff'), width = 10, height = 5, dpi = 300)
