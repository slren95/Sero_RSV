rm(list=ls())
library(tidyverse)
library(rio)
library(psych)
library(tidyxl)
library(janitor)
library(readxl)
library(patchwork)
library(openxlsx)
library(fitdistrplus)


## Manual select Include data----
df_new<-read_excel('data/KKH-V2-RSV Clinical Dataset for Model Analysis-20260619_new.xlsx') %>%
  clean_names() %>%
  mutate(across(6:9,~as.Date(as.numeric(.x), origin = "1899-12-30"))) %>%
  mutate(visit_1_date=ymd(visit_1_date)) %>%
  mutate(subject_ID=str_sub(sample_visit,1,8),.before = 1) %>%
  separate(age,into=c('age_value','age_unit'),remove = F) %>%
  mutate(
    age_unit = tolower(age_unit),
    age_value = as.integer(age_value),
    DOB = case_when(
      age_unit == "years"  ~ visit_1_date %m-% years(age_value),
      age_unit == "months" ~ visit_1_date %m-% months(age_value),
      age_unit == "weeks"  ~ visit_1_date %m-% weeks(age_value),
      age_unit == "days"   ~ visit_1_date %m-% days(age_value),
      TRUE ~ NA_Date_
    )
  ) %>%
  mutate(age_month_1=interval(DOB,visit_1_date)/months(1),
         age_month_2=interval(DOB,visit_2_date)/months(1),
         age_month_3=interval(DOB,visit_3_date)/months(1),
         age_month_4=interval(DOB,visit_4_date)/months(1),
         age_month_5=interval(DOB,visit_5_date)/months(1)) %>%
  mutate(
    age_visit = case_when(
      visit == 1 ~ age_month_1,
      visit == 2 ~ age_month_2,
      visit == 3 ~ age_month_3,
      visit == 4 ~ age_month_4,
      visit == 5 ~ age_month_5,
      TRUE ~ NA_real_
    ),.after=age_unit
  ) %>%
  group_by(subject_ID) %>%
  mutate(visit_interval=age_visit-lag(age_visit),.after=age_visit,
         subject_include_site=as.integer(sum(include,na.rm=T)>0)) %>%
  ungroup() %>%
  rename_with(~str_extract(.x, "^[^_]+_[^_]+"),.cols = starts_with("site")) %>%
  relocate(DOB,starts_with('age_month_'),.after = visit_5_date) %>%
  relocate(subject_include_site,.after = include)

df_new %>%
  filter(age_visit<=30)  %>%
  distinct(subject_ID)

df_new.2<-df_new %>%
  filter(age_visit<=30) %>% filter(include==1) %>%
  add_count(subject_ID)

rio::export(df_new.2,'docs/df_new.2.xlsx')

df_new.2 %>% distinct(subject_ID)



# Fit Wanning ----
fit_simple<-readRDS('rda/fit_simple.rds')

site_cols <- c('site_iii', 'site_ii','site_v', 'site_ph','site_iv')


B_epi<-rstan::summary(fit_simple, pars = "B_epi")$summary[, "mean"] %>%
  setNames(site_cols)

qlogis(0.5)
logit(0.5)
qlogis(B_epi)

df_logit<-df_new.2 %>%
  arrange(subject_ID, age_visit) %>%
  group_by(subject_ID) %>%
  dplyr::select(subject_ID, visit,age_visit, visit_interval,all_of(site_cols)) %>%
  mutate(across(all_of(site_cols),~{
    current_col <- cur_column()
    #print(current_col)
    qlogis(.x/100) - qlogis(B_epi[current_col])
  },.names = '{.col}_logit')) %>%
  filter(if_all(ends_with('_logit'), ~ .x > 0)) %>%
  mutate(across(ends_with('_logit'),~{
    k=-(log(.x)-log(lag(.x)))/visit_interval
  },.names = 'k_{.col}')) %>%
  ungroup() %>%
  mutate(across(
    starts_with("k_"),
    ~ log(2) / .x,
    .names = "t12_{.col}"
  ))

boot::inv.logit(df_logit[[1,'site_iii_logit']]+logit(B_epi[1]))

summary(df_logit %>% dplyr::select(starts_with('t12_'),starts_with('k_')))


df_t12 <- df_logit %>%
  dplyr::select(subject_ID, age_visit, visit,starts_with("t12_k_")) %>%
  pivot_longer(
    cols = starts_with("t12_k_"),
    names_to = "site",
    values_to = "t12_months"
  ) %>%
  mutate(
    age_group = cut(
      age_visit,
      breaks = c(0, 6, 12, 24, 30),
      labels = c("0-6m", "6-12m", "12-24m", "24-30m")
    )
  ) %>%
  filter(!is.na(t12_months))


ggplot(df_t12, aes(x = t12_months, fill = factor(visit))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ site, scales = "free") +
  theme_bw() +
  labs(x = "Half-life (days)", y = "Density")

ggplot(df_t12, aes(x = t12_months)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~ site, scales = "free") +
  theme_bw() +
  labs(x = "Half-life (days)")

## DF_t12 ----
DF_t12<-df_logit %>%
  dplyr::select(starts_with("t12_")) %>%
  pivot_longer(cols = everything(),
               names_to = "site",
               values_to = "t12") %>%
  group_by(site) %>%
  summarise(
    mean = mean(t12, na.rm = TRUE),
    median = median(t12, na.rm = TRUE),
    q1= quantile(t12,0.25,na.rm=T),
    q3= quantile(t12,0.75,na.rm=T),
    .groups = "drop"
  ) %>%
  mutate(Site=str_extract(site, "(?<=t12_k_).*(?=_logit)"),
         Site=factor(Site,levels=c("site_iii","site_ii","site_v","site_ph","site_iv"))) %>%
  arrange(Site)

DF_t12

# # A tibble: 5 × 4
# site                  mean median Site
# <chr>                <dbl>  <dbl> <fct>
#   1 t12_k_site_iii_logit 14.5   11.0  site_iii
# 2 t12_k_site_ii_logit  14.9   13.8  site_ii
# 3 t12_k_site_v_logit   21.4   15.1  site_v
# 4 t12_k_site_ph_logit   5.09   4.28 site_ph
# 5 t12_k_site_iv_logit  22.8    8.74 site_iv

export(DF_t12,'rda/DF_t12.rds')

# Excel with Plot----


# 需要格式化为保留 1 位小数的列名
target_cols <- c(
  "age_visit", "visit_interval", "ig_g_net_mfi_new_beads", "ig_a_net_mfi_new_beads",
  "ig_a_normalized_net_mfi", "ig_g_normalized_net_mfi", "site_iii", "site_ii",
  "site_v", "site_ph", "site_iv","age_month_2","age_month_3","age_month_4","age_month_5"
)

# 2. 数据预处理
df_clean <- df_new %>%
  mutate(across(all_of(c("age_visit", target_cols)), as.numeric)) %>%
  # 对指定列保留 1 位小数
  mutate(across(all_of(target_cols), ~ round(.x, 1))) %>%
  # 按照 subject_ID 和 age_visit 排序以确保折线图顺序正确
  arrange(subject_ID, age_visit)

# 3. 计算每个 subject_ID 对应的行合并区间（Start Row 到 End Row）
subject_layout <- df_clean %>%
  mutate(row_idx = row_number() + 1) %>%
  group_by(subject_ID) %>%
  summarise(
    start_row = min(row_idx),
    end_row = max(row_idx),
    .groups = 'drop'
  )

# 4. 创建工作簿并写入数据
wb <- createWorkbook()
sheet_name <- "with_plot"
addWorksheet(wb, sheet_name)
writeData(wb, sheet_name, df_clean)

# 追加三个绘图列
total_cols <- ncol(df_clean)
col_igg <- total_cols + 2
col_iga <- total_cols + 3
col_sites <- total_cols + 1

writeData(wb, sheet_name, "Plot: IgG Net MFI", startCol = col_igg, startRow = 1)
writeData(wb, sheet_name, "Plot: IgA Net MFI", startCol = col_iga, startRow = 1)
writeData(wb, sheet_name, "Plot: 5 Sites Trend", startCol = col_sites, startRow = 1)

# 设置绘图列宽
setColWidths(wb, sheet_name, cols = c(col_igg, col_iga, col_sites), widths = 40)

# 5. 合并单元格、调整高度，并为每个患者的最后一行添加下划线 (Bottom Border)
# 定义下划线样式（中等灰色，优雅大方）
bottom_line_style <- createStyle(
  border = "bottom",
  borderColour = "#808080",
  borderStyle = "medium"
)

for (i in 1:nrow(subject_layout)) {
  s_row <- subject_layout$start_row[i]
  e_row <- subject_layout$end_row[i]

  mergeCells(wb, sheet_name, cols = col_igg, rows = s_row:e_row)
  mergeCells(wb, sheet_name, cols = col_iga, rows = s_row:e_row)
  mergeCells(wb, sheet_name, cols = col_sites, rows = s_row:e_row)

  num_rows <- e_row - s_row + 1
  row_height <- max(120 / num_rows, 20)
  setRowHeights(wb, sheet_name, rows = s_row:e_row, heights = rep(row_height, num_rows))

  # 💡 核心改动：为该患者的末尾行 (e_row) 的所有列（数据列 + 绘图列）添加底部边框
  # 这样可以在不同 subject_ID 之间拉起一条清晰的“下划隔离线”
  addStyle(
    wb, sheet_name,
    style = bottom_line_style,
    rows = e_row,
    cols = 1:col_sites,
    gridExpand = TRUE,
    stack = TRUE # stack = TRUE 允许保留原单元格的其他样式（例如色阶填充等）
  )
}

# 6. 定义一个生成“X 轴 Label 动态着色”的辅助函数
# 传入该 Subject 的所有 age_visit 刻度，返回对应的颜色向量
get_label_colors <- function(breaks) {
  ifelse(breaks < 30, "#FF0000", "#333333")
}

# 7. 循环绘制 ggplot 并插入
for (i in 1:nrow(subject_layout)) {
  s_id <- subject_layout$subject_ID[i]
  s_row <- subject_layout$start_row[i]

  # 提取画图数据并确保 age_visit 是保留 1 位小数的唯一值序列（作为 X 轴刻度）
  sub_data <- df_clean %>% filter(subject_ID == s_id)
  x_breaks <- unique(sub_data$age_visit)
  x_colors <- get_label_colors(x_breaks)

  # 定义带动态刻度着色的基础主题
  dynamic_theme <- theme_minimal(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8, color = "#333333"),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5),
      # 核心：将算出来的颜色向量应用到 X 轴刻度标签上
      axis.text.x = element_text(color = x_colors, face = "bold")
    )

  ### 图 1: IgG 趋势图 ------
  p_igg <- ggplot(sub_data, aes(x = age_visit, y = ig_g_net_mfi_new_beads)) +
    geom_line(color = "#005F73", size = 1) +
    geom_point(color = "#005F73", size = 2) +
    scale_x_continuous(breaks = x_breaks) + # 强制显示所有随访点的 age_visit
    labs(title = paste("IgG:", s_id), x = "Age at Visit", y = "MFI") +
    dynamic_theme

 ## 图 2: IgA 趋势图 ------
  p_iga <- ggplot(sub_data, aes(x = age_visit, y = ig_a_net_mfi_new_beads)) +
    geom_line(color = "#AE2012", size = 1) +
    geom_point(color = "#AE2012", size = 2) +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = paste("IgA:", s_id), x = "Age at Visit", y = "MFI") +
    dynamic_theme

  ### ------ 图 3: 5 Sites 趋势图 ------
  sub_sites_long <- sub_data %>%
    dplyr::select(age_visit, site_iii, site_ii, site_v, site_ph, site_iv) %>%
    pivot_longer(cols = starts_with("site_"), names_to = "Site", values_to = "Value")

  p_sites <- ggplot(sub_sites_long, aes(x = age_visit, y = Value, color = Site)) +
    geom_line(size = 0.8) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "5 Sites Co-plot", x = "Age at Visit", y = "Value") +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8, color = "#333333"),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      plot.margin = margin(5, 5, 5, 5),
      # 5 sites 图同样应用动态刻度着色
      axis.text.x = element_text(color = x_colors, face = "bold")
    )

  # 打印并插入
  print(p_igg)
  insertPlot(wb, sheet_name, startCol = col_igg, startRow = s_row,
             width = 2.8, height = 1.6, fileType = "png", units = "in")

  print(p_iga)
  insertPlot(wb, sheet_name, startCol = col_iga, startRow = s_row,
             width = 2.8, height = 1.6, fileType = "png", units = "in")

  print(p_sites)
  insertPlot(wb, sheet_name, startCol = col_sites, startRow = s_row,
             width = 2.8, height = 1.6, fileType = "png", units = "in")
}

# 8. 为目标数据列应用二元色（Dual-Color）条件格式化
col_indices <- which(names(df_clean) %in% c(target_cols,'include'))
dual_palette <- c("#FFFFFF", "#005F73") # 优雅的珍珠白到淡天蓝色

for (col_idx in col_indices) {
  conditionalFormatting(
    wb, sheet_name,
    cols = col_idx,
    rows = 2:(nrow(df_clean) + 1),
    style = dual_palette,
    rule = NULL,
    type = "colourScale"
  )
}

# 9. 保存文件
saveWorkbook(wb, "data/kkh_with_plot.xlsx", overwrite = TRUE)

# Plot Waning ----
names(df_new.2)[c(19:23)] %>%
  map(~{
    df_new.2 %>%
      ggplot(aes(age_visit,!!sym(.x),color=subject_ID))+
      geom_line(aes(group=subject_ID))+
      geom_point()+
      labs(title=.x)+
      theme_bw()+
      theme(legend.position = 'none')
  }) %>%
  wrap_plots(ncol=3)

ggsave("plot/Waning_site.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")

df_waning<-names(df_new.2)[c(19:23)] %>%
  map_dfr(~{
    df_new.2 %>%
      group_by(subject_ID) %>%
      arrange(visit) %>%
      slice_head(n=1) %>%
      ungroup()
  }) %>%
  summarise(across(all_of(site_cols),~mean(.x))) %>%
  pivot_longer(cols=everything(),names_to = 'Site',values_to = 'start') %>%
  mutate(base=B_epi[Site]) %>%
  left_join(DF_t12 %>% transmute(Site,t12_month_median=median)) %>%
  pmap_dfr(
    function(Site, start, base, t12_month_median, ...) {
      data.frame(
        Site = Site,
        t = seq(0, 24, by = .5)
      ) %>%
        mutate(
          t12=t12_month_median,
          value = boot::inv.logit(
            (qlogis(start / 100) - qlogis(base)) * exp(-log(2) / t12_month_median * t) + qlogis(base)
          )
        )
    }
  )

df_waning %>%
  ggplot(aes(t,value,color=Site,group=Site))+
  geom_point()+
  geom_line()

names(df_new.2)[c(25:29)] %>%
  map(~{
    df_new.2 %>%
      group_by(subject_ID) %>%
      arrange(visit) %>%
      mutate(time=age_visit-first(age_visit)) %>%
      ggplot(aes(time,!!sym(.x)/100,color=subject_ID))+
      geom_line(aes(group=subject_ID),alpha = .5)+
      geom_point(alpha = .5)+
      theme_bw(base_size=14)+
      theme(legend.position = 'right')+
      geom_line(
        data = df_waning %>% filter(Site==.x),
        aes(x = t, y = value),
        color = "red",
        linewidth = 1,
        inherit.aes = FALSE
      ) +
      scale_y_continuous(labels = scales::percent_format(),limits = c(0,NA),expand = expansion(add = c(0,0),mult=c(0,.1)))+
      labs(
        title = .x,color='Subject ID',
        caption = sprintf('half-life period: %s months',round(DF_t12[DF_t12$Site==.x,'median'],1)),
        x = "Time since baseline (months)",
        y = "Inhibition (%)"
      )+
      theme(legend.title = element_text(hjust = .5))
  }) %>%
  wrap_plots(ncol=3) +
  plot_layout(guides = 'collect',axis_titles = 'collect')

ggsave("plot/Waning_site_zero.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")

# Fitting ----
head(df_t12)

# Waning rate (/month)
df_decay <- df_t12 %>% filter(is.finite(t12_months), t12_months > 0) %>% mutate(r = log(2) / t12_months)

# Sample size
df_decay %>% group_by(site) %>% summarise(n = n(), n_subject = n_distinct(subject_ID), median = median(r), Q1 = quantile(r, .25), Q3 = quantile(r, .75))

# Fit distributions
fits <- df_decay %>% split(.$site) %>% map(~list(
  lnorm = fitdist(.x$r, "lnorm"),
  gamma = fitdist(.x$r, "gamma"),
  weibull = fitdist(.x$r, "weibull")
))

# Compare AIC
fit_summary <- imap_dfr(fits, ~tibble(
  site = .y,
  distribution = c("Lognormal", "Gamma", "Weibull"),
  AIC = c(.x$lnorm$aic, .x$gamma$aic, .x$weibull$aic)
)) %>% group_by(site) %>% mutate(delta_AIC = AIC - min(AIC)) %>% arrange(site, AIC)

fit_summary

# Fitted parameters
imap(fits, ~list(site = .y, lognormal = .x$lnorm$estimate, gamma = .x$gamma$estimate, weibull = .x$weibull$estimate))



epi_prior <- imap_dfr(fits, ~tibble(site=.y, meanlog=.x$lnorm$estimate["meanlog"], sdlog=.x$lnorm$estimate["sdlog"])) %>%
  mutate(
    # 采用对数正态分布的中位数: exp(meanlog)
    rate_median = exp(meanlog),

    # 采用对数正态分布的均值 (期望值): exp(meanlog + (sdlog^2) / 2)
    rate_mean = exp(meanlog + (sdlog^2) / 2),

    # 2. 对应的原始半衰期 (Half-life = ln(2) / rate)
    # 根据中位数速率计算的半衰期
    half_life_median = log(2) / rate_median,

    # 根据均值速率计算的半衰期
    half_life_mean = log(2) / rate_mean
  )

epi_prior

epi_order <- c("t12_k_site_iii_logit", "t12_k_site_ii_logit", "t12_k_site_v_logit",
               "t12_k_site_ph_logit", "t12_k_site_iv_logit")

epi_prior_stan <- epi_prior %>% slice(match(epi_order, site))
epi_prior_stan

DF_t12

export(epi_prior_stan,'rda/epi_prior_stan.rds')


# Plot 5 visit ----
df_new_lt_30<-df_new %>%
  filter(age_visit<=30)

head(df_new_lt_30)

names(df_new_lt_30)

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. 整理数据
# ============================================================

markers <- c(
  "ig_a_normalized_net_mfi",
  "ig_g_normalized_net_mfi",
  "site_iii",
  "site_ii",
  "site_v",
  "site_ph",
  "site_iv"
)

marker_labels <- c(
  ig_a_normalized_net_mfi = "IgA",
  ig_g_normalized_net_mfi = "IgG",
  site_iii = "Site III–MPE8",
  site_ii  = "Site II–motavizumab",
  site_v   = "Site V–suptavumab",
  site_ph  = "Site Φ–nirsevimab",
  site_iv  = "Site IV–clesrovimab"
)

df_long <- df_new_lt_30 %>%
  dplyr::select(
    subject_ID, sample_visit, visit, age_visit,
    subject_include_site, include,
    all_of(markers)
  ) %>%
  pivot_longer(
    cols = all_of(markers),
    names_to = "marker",
    values_to = "value"
  ) %>%
  mutate(
    marker = factor(
      marker,
      levels = markers,
      labels = marker_labels
    )
  ) %>%
  arrange(marker, subject_ID, age_visit)

# ============================================================
# 2. 提取 include == 1 对应的相邻两次检测
# ============================================================

# 假设 include == 1 标记的是下降线段的终点
df_segments <- df_long %>%
  dplyr::arrange(marker, subject_ID, visit) %>%
  dplyr::group_by(marker, subject_ID) %>%
  dplyr::mutate(
    x_start = age_visit,
    y_start = value,
    x_end = dplyr::lead(age_visit),
    y_end = dplyr::lead(value),
    next_include = dplyr::lead(include)
  ) %>%
  dplyr::filter(
    subject_include_site == 1,
    include == 1,
    next_include == 1,
    !is.na(y_start),
    !is.na(y_end)
  ) %>%
  dplyr::ungroup()

# ============================================================
# 3. 绘图函数：一次绘制一个 marker
# ============================================================
# 为所有高亮儿童分配固定颜色，确保不同 marker 颜色一致
ids <- sort(unique(as.character(df_segments$subject_ID)))
id_cols <- setNames(
  grDevices::hcl.colors(length(ids), "Dark 3"),
  ids
)

plot_marker <- function(marker_name) {
  dat <- dplyr::filter(df_long, marker == marker_name, !is.na(value))
  seg <- dplyr::filter(df_segments, marker == marker_name)

  ggplot() +
    geom_line(data = dat, aes(age_visit, value, group = subject_ID),
              colour = "grey78", linewidth = 0.35, alpha = 0.55) +
    geom_point(data = dat, aes(age_visit, value),
               colour = "grey70", size = 1.3, alpha = 0.65) +
    geom_segment(data = seg,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                     colour = as.character(subject_ID)), linewidth = 1) +
    geom_point(data = seg, aes(x_start, y_start, colour = as.character(subject_ID)),
               size = 2) +
    geom_point(data = seg, aes(x_end, y_end, colour = as.character(subject_ID)),
               size = 2) +
    scale_colour_manual(values = id_cols) +
    scale_x_continuous(breaks = seq(0, 30, 6), limits = c(0, 30)) +
    labs(title = marker_name, x = "Age at visit (months)",
         y = "Epitope-specific inhibition (%)",
         colour = "Subject ID") +
    theme_bw() +
    theme(legend.position = 'right', panel.grid.minor = element_blank(),
          plot.title = element_text())
}

# ============================================================
# 4. 单独绘图
# ============================================================

plots<-c("IgA","IgG","Site II–motavizumab","Site III–MPE8",
         "Site IV–clesrovimab","Site V–suptavumab","Site Φ–nirsevimab") %>%
  map(~plot_marker(.x))

wrap_plots(plots[3:7])+
  plot_layout(axis_titles = 'collect',guides = 'collect')

ggsave("plot/Figure_muti_subject.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")


