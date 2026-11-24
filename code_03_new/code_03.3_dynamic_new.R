#rm(list=ls())
library(tidyverse)
library(boot)
library(posterior)


# MODEL_NAME="fit_main"

source('code_03_new/prepare_model_data_new.R', local=TRUE, encoding='UTF-8')

fit_truncat <- readRDS(file.path('rda_new', paste0(MODEL_NAME, '.rds')))
OUTPUT_DIR <- file.path('code_03_new',MODEL_NAME)
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

DF_t12<-import('rda/DF_t12.rds',trust=T)

df_means<-enframe(get_posterior_mean(fit_truncat)[, 5])

## peak value ----
df_real_peaks <- df_means %>%
  pivot_wider(names_from = name, values_from = value) %>%
  pivot_longer(cols = starts_with("peak["), names_to = "peak_name", values_to = "inf_peak_delta") %>%
  mutate(
    idx = as.integer(str_extract(peak_name, "\\d+")),
    marker = c("IgA-MFI", "IgG-MFI", "SiteIII-MPE8", "Site II-Motavizumab",
               "Site V-Suptavumab", "Site Φ-Nirsevimab", "Site IV-Clesrovimab")[idx]
  ) %>%
  rowwise() %>%
  mutate(
    inf_peak_value = case_when(
      idx == 1 ~ B_IgA + inf_peak_delta,
      idx == 2 ~ B_IgG + inf_peak_delta,
      idx == 3 ~ boot::inv.logit(boot::logit(`B_epi[1]`) + inf_peak_delta),
      idx == 4 ~ boot::inv.logit(boot::logit(`B_epi[2]`) + inf_peak_delta),
      idx == 5 ~ boot::inv.logit(boot::logit(`B_epi[3]`) + inf_peak_delta),
      idx == 6 ~ boot::inv.logit(boot::logit(`B_epi[4]`) + inf_peak_delta),
      idx == 7 ~ boot::inv.logit(boot::logit(`B_epi[5]`) + inf_peak_delta)
    ),
    mat_peak_value = case_when(
      idx == 1 ~ B_IgA,
      idx == 2 ~ B_IgG + Mat0,
      idx == 3 ~ boot::inv.logit(boot::logit(`B_epi[1]`) + `Mat0_k[1]`),
      idx == 4 ~ boot::inv.logit(boot::logit(`B_epi[2]`) + `Mat0_k[2]`),
      idx == 5 ~ boot::inv.logit(boot::logit(`B_epi[3]`) + `Mat0_k[3]`),
      idx == 6 ~ boot::inv.logit(boot::logit(`B_epi[4]`) + `Mat0_k[4]`),
      idx == 7 ~ boot::inv.logit(boot::logit(`B_epi[5]`) + `Mat0_k[5]`)
    )
  ) %>%
  ungroup() %>%
  select(marker,B_IgA,B_IgG,starts_with('B_epi['),Mat0,starts_with('Mat0_k['),inf_peak_delta, inf_peak_value,mat_peak_value)

print(df_real_peaks)


# --- 1. 从数据框中动态提取所需的标量参数 ---
# 由于每一行这些全局参数都一样，我们直接取第一行
p_base <- df_real_peaks %>% slice(1)

b_iga   <- p_base$B_IgA
b_igg   <- p_base$B_IgG
mat0    <- p_base$Mat0

# 提取 5 个表位的背景 B_epi
b_epi_vec <- c(p_base$`B_epi[1]`, p_base$`B_epi[2]`, p_base$`B_epi[3]`, p_base$`B_epi[4]`, p_base$`B_epi[5]`)
# 提取 5 个表位的母传贡献 Mat0_k
mat0_k_vec <- c(p_base$`Mat0_k[1]`, p_base$`Mat0_k[2]`, p_base$`Mat0_k[3]`, p_base$`Mat0_k[4]`, p_base$`Mat0_k[5]`)

# 从 7 行数据中依次提取 peak 的增量 (inf_peak_delta)
# 对应顺序：1:IgA, 2:IgG, 3-7:表位
peaks_vec <- df_real_peaks$inf_peak_delta

# 衰减速率常量
rate_draws <- rstan::extract(fit_truncat, pars=c("r_mat","r1","r2"))
r_mat <- mean(rate_draws$r_mat)
r1 <- mean(rate_draws$r1)
r2 <- mean(rate_draws$r2)

r_epi <- apply(rstan::extract(fit_truncat, pars="r_epi")$r_epi, 2, mean)



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

    # 3. 计算 5 个表位的抑制率轨迹 (Logit 尺度下做加和，再 inv.logit 转换)
    mu_epi <- numeric(5)
    for (k in 1:5) {
      epi_S0 <- inv.logit(logit(b_epi_vec[k]) + mat0_k_vec[k] * exp(-r_mat * t))
      mu_epi[k] <- if (has_peaked) {
        inv.logit(logit(epi_S0) + peaks_vec[k+2] * exp(-r_epi[k] * dt))
      } else {
        epi_S0
      }
    }

    # 返回当前时间点的结果
    tibble(
      Age = t,
      `IgA-MFI` = mu_IgA,
      `IgG-MFI` = mu_IgG,
      `SiteIII-MPE8` = mu_epi[1],
      `Site II-Motavizumab` = mu_epi[2],
      `Site V-Suptavumab` = mu_epi[3],
      `Site Φ-Nirsevimab` = mu_epi[4],
      `Site IV-Clesrovimab` = mu_epi[5]
    )
  })
}

# --- 3. 运行模拟并转换数据格式 ---
# 假设出生后第 6 个月感染
df_plot_data <- get_kinetics_trajectory(m_inf = 6) %>%
  pivot_longer(-Age, names_to = "Marker", values_to = "Value") %>%
  # 为了分面美观，区分 MFI 尺度和抑制率百分比尺度
  mutate(Scale_Group = ifelse(str_detect(Marker, "MFI"), "Total Antibody (MFI)", "Epitope Inhibition Rate (0-1)"))

# --- 4. 优雅的 Apple-style 极简可视化 ---
ggplot(df_plot_data, aes(x = Age, y = Value, color = Marker)) +
  # 辅助线表明生物学逻辑进程
  geom_vline(xintercept = 6, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 7, linetype = "dotted", color = "grey50") +
  geom_line(size = 1.2, alpha = 0.85) +
  facet_wrap(~Marker, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = seq(0, 23, by = 2)) +
  labs(
    title = "Model-Inferred RSV Antibody & Epitope Kinetics",
    subtitle = "Simulated primary infection at 6 months (peaking at 7 months)",
    x = "Age (Months)",
    y = "Antibody (MFI) / Epitope Inhibition Rate",
    caption = "Dashed line = Infection onset (6m) | Dotted line = Peak concentration (7m)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10, colour = "grey20"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=4)),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b=15))
  )

ggsave(file.path(OUTPUT_DIR,'antibody_dynamic.tiff'), width = 12, height = 6, dpi = 300)
