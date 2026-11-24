dir.create('code_04/plot', recursive = TRUE, showWarnings = FALSE)

library(patchwork)
library(cowplot)

MODEL_NAME<-'fit_main'

source(sprintf('code_02/code_02_%s.R',MODEL_NAME))

fit_truncat <- readRDS(sprintf('rda/%s.rds',MODEL_NAME))
fit_truncat

cat(capture.output(print(fit_truncat)), file = paste0('code_03/',MODEL_NAME,'/fit_truncat_output.txt'), sep = "\n")

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

  #ggsave(sprintf('code_03/%s/w_by_age_%s.jpg',MODEL_NAME,MAX_AGE),width=8,height=4)
})

plots_w[[1]]




# 2. 累计感染概率 ----
df_individual_wide_all <- import(paste0('code_03/',MODEL_NAME,'/df_individual_wide_all.rds'), trust = TRUE)

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

  #ggsave(sprintf('code_03/%s/cumulative_infection_with_CI_%s.jpg',MODEL_NAME,MAX_AGE), width = 8, height = 4, dpi = 300)

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

  #ggsave(sprintf('code_03/%s/instantaneous_infection_with_CI_%s.jpg',MODEL_NAME,MAX_AGE), width = 8, height = 4, dpi = 300)

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

ggsave("code_04/plot/Figure_2.tiff", width = 12, height = 8, dpi = 300,compression = "lzw")

