dir.create('code_04/plot', recursive = TRUE, showWarnings = FALSE)

library(rio)
library(janitor)
library(scales)
library(ggpubr)

df_sero.4_29<-import('rda/df_sero.4_29months_705.xlsx') %>%
  clean_names() %>%
  mutate(age_f=factor(age_value,levels=5:29))

names(df_sero.4_29)

plot_iga<-ggviolin(df_sero.4_29, x = "age_f", y = "ig_a_mfi", palette = "npg", alpha = 0.35) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 0.5, alpha = 0.4, color = "#E74C3C") +
  scale_y_log10(labels = scales::label_number(accuracy = 1),
                limits = c(0.2, NA),
                expand = expansion(mult = c(0.001, 0.05))) +
  scale_x_discrete(limits = as.character(5:29), breaks = as.character(5:29),
                   expand = expansion(add = c(0.1, 0.7)),
                   labels = function(x) sprintf("%s-<%s", x, as.numeric(x) + 1)) +
  labs(x = "Age (months)", y = "IgA MFI") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
        axis.title = element_text(face = "bold",size=10),
        axis.text = element_text(color = "black"),
        axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6))

plot_igg<-ggviolin(df_sero.4_29, x = "age_f", y = "ig_g_mfi", palette = "npg", alpha = 0.35) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 0.5, alpha = 0.4, color =  "#3498DB") +
  scale_y_log10(labels = scales::label_number(accuracy = 1),
                limits = c(0.2, NA),
                expand = expansion(mult = c(0.001, 0.05))) +
  scale_x_discrete(limits = as.character(5:29), breaks = as.character(5:29),
                   expand = expansion(add = c(0.1, 0.7)),
                   labels = function(x) sprintf("%s-<%s", x, as.numeric(x) + 1)) +
  labs(x = "Age (months)", y = "IgA MFI") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
        axis.title = element_text(face = "bold",size=10),
        axis.text = element_text(color = "black"),
        axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6))

plot_iga+plot_igg+
  plot_layout(nrow=2,axes = 'collect_x')+
  plot_annotation(tag_levels = 'A')

ggsave("code_04/plot/Figure_1_ig.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")



plot_info <- tribble(
  ~var,                    ~label,                    ~color,
  "site_iii_mpe8",         "Site III–MPE8",           "#2C7FB8",
  "site_ii_motavizumab",   "Site II–Motavizumab",     "#E67E22",
  "site_v_suptavumab",     "Site V–Suptavumab",       "#27AE60",
  "site_ph_nirsevimab",    "Site Φ–Nirsevimab",       "#8E44AD",
  "site_iv_clesrovimab",   "Site IV–Clesrovimab",     "#D35400"
)


p_epi <- pmap(
  plot_info,
  function(var, label, color) {
    ggviolin(df_sero.4_29, x = "age_f", y = var, palette = "npg", alpha = 0.35) +
      geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 0.5, alpha = 0.4, color = color) +
      scale_x_discrete(limits = as.character(5:29), breaks = as.character(5:29),
                       expand = expansion(add = c(0.1, 0.7)),
                       labels = function(x) sprintf("%s-<%s", x, as.numeric(x) + 1)) +
      scale_y_continuous(limits = c(0,1),labels=percent_format(),expand = expansion(add=c(0.1,0.1)))+
      labs(x = "Age (months)", y = label) +
      theme_bw() +
      theme(legend.position = "none")
  }
)

wrap_plots(p_epi)+
  plot_layout(nrow=5,axes = 'collect_x')+
  plot_annotation(tag_levels = 'A')

# NEW Cut ----

df_sero.4_cut <- df_sero.4_29 %>%
  mutate(age_f = cut(age_value, breaks = c(5, 9, 12, 15, 18, 21, 24, 27, 30), right = FALSE,
                     labels = c("5-<9","9-<12", "12-<15", "15-<18", "18-<21", "21-<24", "24-<27", "27-<30")))

plot_iga<-ggviolin(df_sero.4_cut, x = "age_f", y = "ig_a_mfi", fill = "#F5B7B1", color = "black",
         alpha = 0.35, width = 0.9, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4, color = "#E74C3C") +
  scale_y_log10(labels = scales::label_number(accuracy = 1),
                limits = c(0.2, NA),
                expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "Age (months)", y = "IgA MFI") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
        axis.title = element_text(size=10),
        axis.text = element_text(color = "black"),
        #axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
        )

plot_igg<-ggviolin(df_sero.4_cut, x = "age_f", y = "ig_g_mfi", fill = "#9ECAE1", color = "black",
                   alpha = 0.35, width = 0.9, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4, color = "#2C7FB8") +
  scale_y_log10(labels = scales::label_number(accuracy = 1),
                limits = c(0.2, NA),
                expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "Age (months)", y = "IgG MFI") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
        axis.title = element_text(size=10),
        axis.text = element_text(color = "black"),
        #axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6)
  )

plot_iga+plot_igg+
  plot_layout(nrow=2,axes = 'collect_x')+
  plot_annotation(tag_levels = 'A')

plot_info <- tribble(
  ~var,                   ~label,      ~color,    ~fill,
  "site_ii_motavizumab",  "Site II",   "#E69F00", "#F8DFA6",  # 橙黄
  "site_iii_mpe8",        "Site III",  "#009E73", "#A8DCCB",  # 青绿
  "site_iv_clesrovimab",  "Site IV",   "#CC79A7", "#EBC6DC",  # 粉紫
  "site_v_suptavumab",    "Site V",    "#7A9A01", "#D5DFA5",  # 橄榄绿
  "site_ph_nirsevimab",   "Site Φ",    "#7B5AA6", "#D4C8E2"   # 紫
)


p_epi <- pmap(
  plot_info,
  function(var, label, color, fill) {
    ggplot(df_sero.4_cut, aes(x = age_f, y = .data[[var]])) +
      geom_jitter(width = 0.18, size = 0.8, alpha = 0.5, color = color) +
      geom_boxplot(width = 0.3, fill = fill, color = "black", outlier.shape = NA, alpha = 0.6) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format(),
                         expand = expansion(add = c(0.1, 0.1))) +
      labs(x = "Age (months)", y = label) +
      theme_bw()+
      theme(        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                    plot.subtitle = element_text(color = "gray40", hjust = 0.5, size = 11),
                    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.6),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
                    axis.title = element_text(size=10),
                    axis.text = element_text(color = "black"),)
  }
)


wrap_plots(p_epi) +
  plot_layout(nrow = 5, axes = "collect_x")

p_ab <- wrap_elements(
  plot_iga + plot_igg +
    plot_layout(nrow = 2, axes = "collect_x")
)

p_epi_all <- wrap_elements(
  wrap_plots(p_epi) +
    plot_layout(nrow = 5, axes = "collect_x")
)

p_ab + p_epi_all +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("code_04/plot/Figure_1.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")

