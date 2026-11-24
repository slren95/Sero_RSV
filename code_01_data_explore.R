library(tidyverse)
library(rio)
library(ggsci)
library(DataExplorer)
library(patchwork)
library(psych)
library(scales)


df_sero<-import('data/RSV Clinical Dataset for Model Analysis (n = 2320)-17042026.xlsx')
df_pos<-import('data/RSV_June 2020 to May 2025_Leo.xlsx')

# Batch 1 Sero data ----
df_sero %>% count(Age)


df_sero.1<-df_sero %>%
  separate(Age,into=c('age_value','age_unit'),remove = F) %>%
  mutate(
    age_unit = tolower(age_unit),
    age_value = as.integer(age_value),
    DOC = dmy(`Date of Collection`),
    DOB = case_when(
      age_unit == "years"  ~ DOC %m-% years(age_value),
      age_unit == "months" ~ DOC %m-% months(age_value),
      age_unit == "weeks"  ~ DOC %m-% weeks(age_value),
      age_unit == "days"   ~ DOC %m-% days(age_value),
      TRUE ~ NA_Date_
    )
  )

df_sero.1 %>% count(age_unit,Age)

#create_report(df_sero.1)

names(df_sero.1)[8:14] %>%
  map(~{
    ggplot(df_sero.1, aes(x = .data[[.x]])) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      labs(title = .x) +
      theme_minimal()
  }) %>%
  wrap_plots(ncol=3)

summary(df_sero.1)

df_sero.2<-df_sero.1 %>%
  mutate(across(contains('Site'),~if_else(.x<=0,1e-6,.x/100)))

summary(df_sero.2)

df_sero.3<-df_sero.2 %>%
  filter(str_detect(Age,'Months')) %>%
  select(`Subject No.`,DOB,DOC,Gender,age_value,age_unit,contains('Site'),contains('Ig'))

lapply(df_sero.3[,c('DOB','DOC')],range)
# $DOB
# [1] "2021-04-21" "2024-05-15"
#
# $DOC
# [1] "2023-03-13" "2024-11-15"



# Positive rate ----

df_pos.1<-df_pos %>%
  fill(Year, .direction='down') %>%
  mutate(MONTH = sprintf('%s-%02d', Year, match(Month, month.abb)),
         pos = rowSums(select(., 3:12), na.rm = TRUE),
         pos_scale=(pos-min(pos))/(max(pos)-min(pos))
  )

df_pos.2<-df_pos.1 %>%
  select(MONTH,pos_scale)


head(df_pos.2)



## RSV positity Curve ----
Sys.setlocale("LC_TIME", "English")
df_pos.1 %>%
  mutate(Date = as.Date(paste0(MONTH, "-01"))) %>%
  ggplot(., aes(x = Date, y = pos)) +
  # 绘制折线和点
  geom_line(color = "salmon", linewidth = 1) +
  geom_point(color = "salmon", fill = "white", size = 2, shape = 21, stroke = 1) +
  # 设置横轴标签格式为“年-月”
  scale_x_date(
    date_breaks = "3 months",     # 每3个月显示一个刻度
    date_labels = "%Y-%b",         # 格式示例：2013-Jan
    expand = expansion(mult = c(.01,.01))
  ) +
  # 细节修饰
  labs(
    x = NULL,
    y = "Number of RSV-positive patients"
  ) +
  theme_bw(base_size = 12)+
  theme(panel.border = element_rect(color = "gray90", fill = NA, linewidth = .5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"))

ggsave('plot/RSV_posivity.jpg', width = 8, height = 4, dpi = 300)

# Susceptibility Demo Plot ----

# 1. 模拟 0-23 月龄数据
# 使用 Gamma 样分布函数模拟易感性：初期保护 -> 峰值 -> 缓慢下降
df_omega <- data.frame(month = 0:23) %>%
  mutate(
    # 模拟函数：(m+2) * exp(-0.25 * m) 产生一个在 2-3 月龄左右的峰值
    w_raw = (month + 2) * exp(-0.25 * month),
    # 归一化，使最大易感性为 1
    w_estimated = w_raw / max(w_raw)
  )

# 2. 绘图
ggplot(df_omega, aes(x = month, y = w_estimated)) +
  # 绘制推断曲线 (Estimated Curve)
  geom_line(color = "#C00000", size = 1.2, linetype = "solid") +
  # 绘制空心点并遮挡背景线
  #geom_point(color = "#C00000", fill = "white", size = 3, shape = 21, stroke = 1.2) +
  # 坐标轴设置
  scale_x_continuous(breaks = seq(0, 23, 2)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  # 标签：呼应你 PPT 中 "Estimated/Unknown" 的标记需求
  labs(
    x = "Age in Months",
    y = "RSV Susceptibility"
  ) +
  # 主题美化
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(color = "#C00000", face = "bold", size = 14)
  )


# Distribution of IgA/IgG epi ----

names(df_sero.3)[7:13] %>%
  walk(~{
    print(ggplot(df_sero.3, aes(!!sym(.x))) +
            geom_histogram(aes(y = ..density..),
                           bins = 30,
                           fill = "steelblue",
                           color = "white",
                           alpha = 0.7) +
            geom_density(color = "red",
                         size = 1.2,
                         kernel = "gaussian") +
            labs(title = paste("Distribution of", .x),
                 x = .x,
                 y = "Density") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  panel.grid.minor = element_blank()))
  })


df_sero.3 %>% dplyr::select(7:13) %>%
  mutate(across(c(everything()),~.x)) %>%
  pairs.panels(method = "pearson", # correlation method
               hist.col = "#00AFBB",
               density = TRUE,  # show density plots
               ellipses = TRUE # show correlation ellipses
  )

df_sero.3 %>% dplyr::select(7:13) %>% skimr::skim() %>%
  capture.output(file = "docs/data_explore_skim.txt")


df_sero.3 %>%
  dplyr::select(7:13) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = "", y = value)) +
  geom_violin(fill = "#00AFBB", color = "darkblue", alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.7, outlier.size = 0.8) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(x = "", y = "Value",
       title = "Violin + Boxplot for Variables 4-8",
       caption = "Violin shows density, box shows summary statistics") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "#f5f5f5", color = NA),
        strip.text = element_text(face = "bold", size = 10))

df_sero.3 %>% count(month(DOB)) %>% plot()

# New identified 24-30 months subjects (2 batch)----

df_extend1<-import('data/2) 24 to 30 Months for subjects (Residual).xlsx') %>% mutate(`Date of Collection`=dmy(`Date of Collection`))
df_extend2<-import('data/KKH-V3-RSV Clinical Dataset for Model Analysis (n=202) 20260722.xlsx') %>% rename(`Subject No.`=`Sample ID`)
df_extend<-bind_rows(list(batch1=df_extend1,batch2=df_extend2),.id = 'batch') %>% add_count(`Subject No.`)

df_extend1$`Subject No.` %in% df_extend2$`Subject No.` # 1 equal 2

df_extend2$`Subject No.` %in% df_sero.3$`Subject No.`
df_extend2$`Subject No.` %in% df_sero$`Subject No.`

df_tem<-df_extend2 %>% left_join(df_sero,by='Subject No.') %>%
  relocate(starts_with('IgA'),starts_with('IgG'))

df_tem %>%
  filter(!is.na(`IgA-MFI`)) %>%
  select(1:6) %>%
  pairs.panels(method = "pearson", # correlation method
              hist.col = "#00AFBB",
              density = TRUE,  # show density plots
              ellipses = TRUE # show correlation ellipses
  )

df_tem %>%
  filter(!is.na(`IgA-MFI`)) %>%
  select(1:6) %>%
  summary()

# use df_extend2 Normalized value as df_sero.IgX-MFI
#myvenn(names(df_extend2),names(df_sero.3))


df_extend2_new<-df_extend2 %>%
  rename(age_value=`Age (months)`,`IgA-MFI`=`IgA Normalized Net MFI`,`IgG-MFI`=`IgG Normalized Net MFI`) %>%
  mutate(age_unit='months',
         DOC = ymd(`Date of Collection`),
         DOB = case_when(
           age_unit == "years"  ~ DOC %m-% years(age_value),
           age_unit == "months" ~ DOC %m-% months(age_value),
           age_unit == "weeks"  ~ DOC %m-% weeks(age_value),
           age_unit == "days"   ~ DOC %m-% days(age_value),
           TRUE ~ NA_Date_
         )) %>%
  select(-matches('(Normalized|Raw)'),-`Date of Collection`) %>%
  mutate(across(contains('Site'),~if_else(.x<=0,1e-6,.x/100))) %>%
  mutate(across(contains('MFI'),~if_else(.x<=0,1e-6,.x)))

#myvenn(names(df_extend2_new),names(df_sero.3))

lapply(df_extend2_new[,c('DOB','DOC')],range)
# $DOB
# [1] "2020-09-13" "2023-04-22"
#
# $DOC
# [1] "2023-03-13" "2025-04-22"

## df_sero.4 ----

df_sero.4<-bind_rows(list(df_sero.3=df_sero.3,df_extend2_new=df_extend2_new),.id = 'from')

export(df_sero.4,'rda/df_sero.4.rds')

# Subject Age ----
df_sero.4 %>%
  filter(age_value<30) %>%
  ggplot(aes(age_value)) +
  geom_bar(fill = "#3182ce", color = "#2c5282", alpha = 0.85, width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  labs(x = "Age at collection", y = "Count") +
  scale_x_continuous(breaks=5:29,labels = function(x) sprintf('%s-<%sm',x-1,x))+
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 30,size=8))

ggsave("plot/Age_at_collection.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")

# Result Population 0-29----
df_sero.4_29<-df_sero.4 %>% filter(age_value<30) %>%
  rename(`Site III-MPE8`=`SiteIII-MPE8`)

df_sero.4_29 %>% count(Gender)
df_sero.4_29 %>% count(age_value)

df_sero.4_29 %>%
  dplyr::select(age_value,starts_with('Site'),ends_with('MFI')) %>%
  summary()

df_sero.4_29 %>%
  count(month(DOC))

lapply(df_sero.4_29[,c('DOB','DOC')],range)


p_iga<-df_sero.4_29 %>%
  ggplot(aes(age_value, `IgA-MFI`, color = Gender)) +
  geom_point() +
  geom_smooth(se = TRUE, method = "loess") +  # 或 method = "gam"
  scale_y_log10(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0.001, 0.05)),  # 底部少扩展，上部稍多
    labels = label_number(accuracy = 1)       # 显示为普通数字，不科学计数法
  ) +
  scale_x_continuous(breaks = seq(5,29,by=2))+
  labs(
    x = "Age (Months)",
    y = "IgA-MFI",
    color = "Gender"
  ) +
  theme_bw()

p_igg<-df_sero.4_29 %>%
  ggplot(aes(age_value, `IgG-MFI`, color = Gender)) +
  geom_point() +
  geom_smooth(se = TRUE, method = "loess") +  # 或 method = "gam"
  scale_y_log10(
    limits = c(0.2, NA),
    expand = expansion(mult = c(0.001, 0.05)),  # 底部少扩展，上部稍多
    labels = label_number(accuracy = 1)       # 显示为普通数字，不科学计数法
  ) +
  scale_x_continuous(breaks = seq(5,29,by=2))+
  labs(
    x = "Age (Months)",
    y = "IgG-MFI",
    color = "Gender"
  ) +
  theme_bw()

p_iga + p_igg + plot_layout(guides = 'collect',axis='collect') &
  theme(legend.position = 'top')

# Figure 1 ----
ggsave("plot/Ig_GA.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")

names(df_sero.4_29)[c(9,8,12,10,11)]

df_sero.4_29 %>%
  pivot_longer(starts_with('Site')) %>%
  #mutate(age_value=factor(age_value,levels=5:29)) %>%
  ggplot(aes(age_value, value,color=Gender)) +
  geom_point(alpha=.5)+
  geom_smooth()+
  facet_wrap(vars(name))+
  scale_y_continuous(labels=percent_format())+
  scale_x_continuous(breaks = seq(5,29,by=2),labels = function(x){sprintf('%s-<%s',x,x+1)})+
  labs(
    x = "Age (Months)",
    y = "Inhibition(%)",
    color = "Gender"
  ) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=30,hjust = 1),
        legend.position = 'top')

ggsave("plot/Site.tiff", width = 12, height = 6, dpi = 300,compression = "lzw")


library(tidyverse)
library(ggridges)

# 定义位点顺序（Φ 放最后）
site_order <- c("Site II-Motavizumab",
                "SiteIII-MPE8",
                "Site IV-Clesrovimab",
                "Site V-Suptavumab",
                "Site Φ-Nirsevimab")

df_sero.4_29 %>%
  pivot_longer(starts_with('Site'),
               names_to = "Site",
               values_to = "MFI") %>%
  mutate(
    age_value = factor(age_value, levels = 5:29),
    Site = factor(Site, levels = site_order)
  ) %>%
  ggplot(aes(x = MFI, y = age_value, fill = Site)) +
  geom_density_ridges(
    alpha = 0.6,
    scale = 1.2,
    rel_min_height = 0.01
  ) +
  scale_x_log10(
    limits = c(0.1, NA),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  facet_wrap(~ Site, ncol = 3, scales = "free_x") +  # 分面！每个 Site 单独一个图
  labs(
    title = "Antibody MFI Distribution by Age and Epitope Site",
    x = "MFI (log10 scale)",
    y = "Age (years)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",           # 分面后不需要图例
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

export(df_sero.4_29,'rda/df_sero.4_29months_705.xlsx')

df_sero.4_29 %>% count(from)

df_sero.4_29 %>%
  select(8:14) %>%
  pairs.panels(method = "spearman", # correlation method
              hist.col = "#00AFBB",
              density = TRUE,  # show density plots
              ellipses = TRUE # show correlation ellipses
  )

rm(list=setdiff(ls(),c('df_sero.4','df_sero.3','df_pos.2')))

save.image('rda/code_01_data_explore.RData')


