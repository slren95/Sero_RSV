library(patchwork)
library(cowplot)
library(GGally)
library(rio)

MODEL_NAME<-'fit_main_new'

source('code_03_new/prepare_model_data_new.R', encoding = 'UTF-8')

fit_truncat <- readRDS(file.path('rda_new', paste0(MODEL_NAME, '.rds')))
fit_truncat

dir.create('code_04_new/plot', recursive = TRUE, showWarnings = FALSE)
cat(capture.output(print(fit_truncat)), file = file.path('code_04_new/plot','fit_main_new_output.txt'), sep = "\n")

df_means<-enframe(get_posterior_mean(fit_truncat)[, 5])

# 易感性曲线 w(a) ----
posterior_w_truncat <- rstan::extract(fit_truncat, "w")$w
w_res_truncat <- data.frame(
  age = 0:(ncol(posterior_w_truncat)-1),
  mean = colMeans(posterior_w_truncat),
  lower = apply(posterior_w_truncat, 2, quantile, 0.025),
  upper = apply(posterior_w_truncat, 2, quantile, 0.975)
)

plots_w<-c(23,ncol(cal_idx_matrix)-1) %>% map(~{
  MAX_AGE<-.x
  ggplot(w_res_truncat, aes(x = age, y = mean)) +
    geom_line(color = "#1A5276", linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#1A5276") +
    labs(#title = paste("RSV Susceptibility (N =", N_current, ")"),
      x = "Age (months)", y = "Relative susceptibility to RSV infection") +
    scale_x_continuous(breaks=seq(0,MAX_AGE,by=1),expand = expansion(add=c(0,.5)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)})+
    scale_y_continuous(expand = expansion())+
    theme_bw()+
    theme(      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
                panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
                panel.grid.minor = element_blank(),
                panel.grid = element_line(linewidth = .2),
                axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
                axis.title = element_text(face = "bold",size=10),
                axis.text = element_text(color = "black"),
          axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6))+
    coord_cartesian(xlim=c(0,MAX_AGE))

  #ggsave(sprintf('%s/w_by_age_%s.jpg',MODEL_NAME,MAX_AGE),width=8,height=4)
})

plots_w[[1]]




# 2. 累计感染概率 ----
df_individual_wide_all <- import(file.path('code_03_new',MODEL_NAME,'df_individual_wide_all.rds'), trust = TRUE)

df_plot_data <- df_individual_wide_all %>%
  select(post_id, ID, starts_with("prob_infect_age_")) %>%
  mutate(across(starts_with("prob_infect_age_"),~na_if(.x, 0))) %>%
  pivot_longer(
    cols = starts_with("prob_infect_age_"),
    names_to = "age",
    values_to = "prob"
  ) %>%
  mutate(age = as.integer(str_remove(age, "prob_infect_age_"))) %>%
  summarise(
    mean_prob = mean(prob, na.rm = TRUE),
    .by = c(post_id, age)
  ) %>%
  mutate(cum_prob = 1-cumprod(1-mean_prob),.by = post_id) %>%
  summarise(across(c(mean_prob,cum_prob),list(est=~quantile(.x,0.5,na.rm=T),
                                              lci=~quantile(.x,0.025,na.rm=T),
                                              uci=~quantile(.x,0.975,na.rm=T))),.by=age)


plot_cums<-c(23,nrow(df_plot_data)-1) %>% map(~{
  MAX_AGE<-.x
  ggplot(df_plot_data, aes(x = age, y = cum_prob_est)) +
    geom_ribbon(aes(ymin = cum_prob_lci, ymax = cum_prob_uci), fill = "#2471A3", alpha = 0.15) +
    geom_line(color = "#2471A3", linewidth = 1.2) +
    geom_point(color = "#E74C3C", size = 2.5, shape = 16) +
    scale_x_continuous(limits = c(NA,MAX_AGE),breaks = seq(0, MAX_AGE, by = 1),
                       expand = expansion(add=c(0.1,.7)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)}) +
    scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linewidth = .2),
      axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
      axis.title = element_text(face = "bold",size=10),
      axis.text = element_text(color = "black"),
      axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
    ) +
    labs(
      x = "Age (months)",
      y = "Cumulative probability of RSV infection"
    )

  #ggsave(sprintf('%s/cumulative_infection_with_CI_%s.jpg',MODEL_NAME,MAX_AGE), width = 8, height = 4, dpi = 300)

})


plot_insts<-c(23,nrow(df_plot_data)-1) %>% map(~{
  MAX_AGE<-.x
  MAX_Y<-ceiling(max(df_plot_data$mean_prob_uci,na.rm=T) / 0.05) * 0.05
  ggplot(df_plot_data, aes(x = age, y = mean_prob_est)) +
    geom_ribbon(aes(ymin = mean_prob_lci, ymax = mean_prob_uci), fill = "#2471A3", alpha = 0.15) +
    geom_line(color = "#2471A3", linewidth = 1.2) +
    geom_point(color = "#E74C3C", size = 2.5, shape = 16) +
    scale_x_continuous(limits = c(0,MAX_AGE),breaks = seq(0, MAX_AGE, by = 1),
                       expand = expansion(add=c(0.1,.7)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)}) +
    scale_y_continuous(limits = c(0, MAX_Y),breaks = seq(0, MAX_Y, by = 0.05), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linewidth = .2),
      axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
      axis.title = element_text(face = "bold",size=10),
      axis.text = element_text(color = "black"),
      axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
    ) +
    labs(
      x = "Age (months)",
      y = "Monthly probability of RSV infection"
    )

  #ggsave(sprintf('%s/instantaneous_infection_with_CI_%s.jpg',MODEL_NAME,MAX_AGE), width = 8, height = 4, dpi = 300)

})

# 灵敏度曲线 ----

pi_draws <- rstan::extract(fit_truncat, pars = "pi")$pi

# 2. 纵向收拢 4000 个世界，计算贝叶斯分位数
df_markers <- tibble(
  marker_idx = 1:7,
  # 对应标书文本：注入具体的抗体与标记物生物学标签
  biological_meaning = c("IgA","IgG","SiteIII-MPE8", "Site II-Motavizumab",
                         "Site V-Suptavumab", "Site Φ-Nirsevimab",
                         "Site IV-Clesrovimab"),
  # 💡 核心对账步：直接从 Stan 后验样本矩阵中切片计算分位数
  median = apply(pi_draws, 2, quantile, 0.50),  # 后验中位数
  lower  = apply(pi_draws, 2, quantile, 0.025), # 95% 信用区间下界
  upper  = apply(pi_draws, 2, quantile, 0.975)  # 95% 信用区间上界
) %>%
  # 严格锁定 Y 轴排列顺序，使 pi[1] 盘踞在图形最上方
  mutate(biological_meaning = factor(biological_meaning, levels = rev(biological_meaning)))

# 3. 顶刊学术风森林图可视化
plot_sens<-ggplot(df_markers, aes(x = median, y = biological_meaning)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.18, color = "#1A5276", linewidth = 0.9) +
  geom_point(color = "#E74C3C", size = 3, shape = 16) +
  geom_text(aes(label = sprintf("%.1f%% (%.1f%%-%.1f%%)", median * 100, lower * 100, upper * 100)),
            hjust = 0.5, vjust = -0.9, size = 3, color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),breaks = c(0.85,0.9,0.95,1),
                     limits = c(NA, 1.02)) +
  theme_bw()+
  labs(x='Sensitivity to RSV infection',y=NULL)+
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(linewidth = .2),
    axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
    axis.title = element_text(face = "bold",size=10),
    axis.text = element_text(color = "black")
  )


plot_grid(
  plots_w[[1]], plot_insts[[1]], plot_sens, plot_cums[[1]],
  ncol = 2,
  labels = c('a', 'c', 'b', 'd')   # 按 左到右、上到下 的顺序对应
)

ggsave("code_04_new/plot/Figure_2.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")


# Parirwise Plot ----
get_p_pairwise<-function(){
  source('code_07_plot_pairwise.R',local = T)
  p
}
p_pairwise<-get_p_pairwise()


plot_sens /
  wrap_elements(full = GGally::ggmatrix_gtable(p_pairwise)) +
  plot_layout(heights = c(2, 5))

# Sens Plot ----
get_p_sens<-function(){
  source('code_04_new/code_04.3_plot_sens_new.R',local = T)
  list(plots_w,plots_s,plots_c,plots_i)
}

p_sens<-get_p_sens()


# Figure patch ----
p_sens[[1]][[7]]

p_sens[[2]][[7]]

p_sens[[3]][[7]]

plot_sens

p_pairwise

font_theme <- theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(size = 9, face = "plain"),
  axis.text = element_text(size = 7, color = "black"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  strip.text = element_text(size = 8)
)

p1 <- p_sens[[1]][[7]] +
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, 0.5)
  )

p_left <- wrap_plots(
  p1,
  p_sens[[4]][[6]] + theme(legend.position = "none") +
    scale_y_continuous(limits = c(0,.4),breaks = seq(0,0.4,by=.1),labels = scales::percent_format()),
  p_sens[[3]][[6]] + theme(legend.position = "none") +
    scale_y_continuous(breaks = seq(0,1,by=.2),labels = scales::percent_format()),
  ncol = 1
) + plot_layout(axis_titles = 'collect') +
  plot_annotation(tag_levels = list(letters[1:3])) & font_theme &
  theme(axis.title.y = element_text(margin = margin(r = 2, l = 12)),
        plot.margin = margin(r=10))

p2 <- p_sens[[2]][[7]]

# 修改图中的点
for (i in seq_along(p2$layers)) {
  if (inherits(p2$layers[[i]]$geom, "GeomPoint")) {
    p2$layers[[i]]$aes_params$shape <- 16
  }
}

# 修改图例中的点
p2 <- p2 + guides(
  colour = guide_legend(override.aes = list(shape = 16))
)


p2

# 右侧：2张
p_right <- p2 /
  wrap_elements(full = GGally::ggmatrix_gtable(p_pairwise)) +
  plot_annotation(tag_levels = list(letters[4:5])) +
  plot_layout(heights = c(2,3)) & font_theme


p_final <- (p_left | p_right) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = "Arial",
      size = 10,
      face = "bold"
    ),
    plot.tag.position = c(0, 1)
  )


p_final

ggsave("plot/Figure_2_new.tiff", width = 10, height = 9, dpi = 300,compression = "lzw")

save.image('code_04_new/code_04.2_plot_figure2_new.RData')

# Sens Plot New ----
p2 +
  wrap_elements(full = GGally::ggmatrix_gtable(p_pairwise)) +
  plot_annotation(tag_levels = list(letters[1:2])) +
  plot_layout(widths = c(2,3.5)) & font_theme

ggsave("code_04_new/plot/Figure_3_pair.tiff", width = 10, height =5, dpi = 300,compression = "lzw")


wrap_plots(
  p1+theme(legend.position = 'null'),
  p_sens[[4]][[7]] +theme(legend.position = 'top')+
    scale_y_continuous(limits = c(0,.4),breaks = seq(0,0.4,by=.1),labels = scales::percent_format()),
  p_sens[[3]][[7]] +theme(legend.position = 'null')+
    scale_y_continuous(breaks = seq(0,1,by=.2),labels = scales::percent_format()),
  nrow = 1
) + plot_layout(axis_titles = 'collect') +
  plot_annotation(tag_levels = list(letters[1:3])) & font_theme

ggsave("code_04_new/plot/Figure_4.tiff", width = 12, height =5, dpi = 300,compression = "lzw")
