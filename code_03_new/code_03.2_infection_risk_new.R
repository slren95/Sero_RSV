library(tidyverse)
library(scales)
library(rio)


# MODEL_NAME="fit_sens3_3"



df_individual_wide_all <- import(file.path('code_03_new',MODEL_NAME,'df_individual_wide_all.rds'), trust=TRUE)



# 2. 累计感染概率 ----

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

if(MODEL_NAME=="fit_sens3_3_new") {
  df_plot_data<-bind_rows(
    tibble(age=0:1),
    df_plot_data
  )
}

df_plot_data %>%
  mutate(across(2:7,~round(.*100,4))) %>%
  filter(age %in% c(0,5,11,23)) %>%
  export(file.path('code_03_new',MODEL_NAME,'cum.txt'))
export(df_plot_data,file.path('code_03_new',MODEL_NAME,'infection_risk_data.txt'))

# 3. 图表渲染 ----

c(23,nrow(df_plot_data)-1) %>% walk(~{
  MAX_AGE<-.x
  ggplot(df_plot_data, aes(x = age, y = cum_prob_est)) +
    geom_ribbon(aes(ymin = cum_prob_lci, ymax = cum_prob_uci), fill = "#2471A3", alpha = 0.15) +
    geom_line(color = "#2471A3", linewidth = 1.2) +
    geom_point(color = "#E74C3C", size = 2.5, shape = 16) +
    scale_x_continuous(limits = c(NA,MAX_AGE),breaks = seq(0, MAX_AGE, by = 1),
                       expand = expansion(add=c(0.1,.5)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)}) +
    scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
    ) +
    labs(
      x = "Age (months)",
      y = "Cumulative probability of RSV infection"
    )

  ggsave(file.path('code_03_new',MODEL_NAME,sprintf('cumulative_infection_with_CI_%s.jpg',MAX_AGE)), width = 8, height = 4, dpi = 300)

})


c(23,nrow(df_plot_data)-1) %>% walk(~{
  MAX_AGE<-.x
  MAX_Y<-ceiling(max(df_plot_data$mean_prob_uci,na.rm=T) / 0.05) * 0.05
  ggplot(df_plot_data, aes(x = age, y = mean_prob_est)) +
    geom_ribbon(aes(ymin = mean_prob_lci, ymax = mean_prob_uci), fill = "#2471A3", alpha = 0.15) +
    geom_line(color = "#2471A3", linewidth = 1.2) +
    geom_point(color = "#E74C3C", size = 2.5, shape = 16) +
    scale_x_continuous(limits = c(0,MAX_AGE),breaks = seq(0, MAX_AGE, by = 1),
                       expand = expansion(add=c(0.1,.5)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)}) +
    scale_y_continuous(limits = c(0, MAX_Y),breaks = seq(0, MAX_Y, by = 0.05), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
    ) +
    labs(
      x = "Age (months)",
      y = "Instantaneous probability of RSV infection"
    )

  ggsave(file.path('code_03_new',MODEL_NAME,sprintf('instantaneous_infection_with_CI_%s.jpg',MAX_AGE)), width = 8, height = 4, dpi = 300)

})
