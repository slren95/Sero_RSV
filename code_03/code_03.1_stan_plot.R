# Post Bayes Model ----
#rm(list=ls())

# MODEL_NAME="fit_sens5"


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

c(23,ncol(cal_idx_matrix)-1) %>% walk(~{
  MAX_AGE<-.x
  ggplot(w_res_truncat, aes(x = age, y = mean)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
    labs(#title = paste("RSV Susceptibility (N =", N_current, ")"),
      x = "Age (months)", y = "Relative susceptibility for RSV infection") +
    scale_x_continuous(breaks=seq(0,MAX_AGE,by=1),expand = expansion(add=c(0,.5)),
                       labels=function(x){sprintf('%s-<%s',x,x+1)})+
    scale_y_continuous(expand = expansion())+
    theme_bw(base_size = 12)+
    theme(panel.border = element_rect(color = "gray90", fill = NA, linewidth = .5),
          #axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          axis.text.x= element_text(angle=-45,hjust=0,vjust=.5,size=6))+
    coord_cartesian(xlim=c(0,MAX_AGE))

  ggsave(sprintf('code_03/%s/w_by_age_%s.jpg',MODEL_NAME,MAX_AGE),width=8,height=4)
})


# MCMC trace Plot ----

fit_truncat@par_dims %>% map_dbl(~length(.x))

color_scheme_set("viridisA")
walk2(fit_truncat@model_pars, fit_truncat@par_dims, ~{
  # 标量参数
  if (length(.y) == 0) {
    p <- mcmc_trace(fit_truncat,pars = .x)
    ggsave(filename = paste0('code_03/',MODEL_NAME,'/par_', .x, ".jpg"),plot = p,
           width = 8,height = 5)

  } else {
    # 非标量参数，生成 Stan 风格参数名
    idx <- expand.grid(lapply(.y, seq_len))
    pars_use <- apply(idx, 1, function(i) {
      paste0(.x, "[", paste(i, collapse = ","), "]")
    })
    p <- mcmc_trace(fit_truncat,pars = pars_use)
    ggsave(filename = paste0('code_03/',MODEL_NAME,'/par_',.x, ".jpg"),plot = p,
           width = 10, height = 6)
  }
})

# Sensitivity Plot ----

pi_draws <- rstan::extract(fit_truncat, pars = "pi")$pi

num_markers <- ncol(pi_draws)

# 1. 动态定义生物学标签
if (num_markers == 2) {
  marker_meanings <- c("IgA", "IgG")
} else {
  marker_meanings <- c("IgA", "IgG", "SiteIII-MPE8", "Site II-Motavizumab",
                       "Site V-Suptavumab", "Site Φ-Nirsevimab",
                       "Site IV-Clesrovimab")
}

# 2. 计算分位数
df_markers <- tibble(
  marker_idx = 1:num_markers,
  biological_meaning = marker_meanings,
  median = apply(pi_draws, 2, quantile, 0.50),
  lower  = apply(pi_draws, 2, quantile, 0.025),
  upper  = apply(pi_draws, 2, quantile, 0.975)
) %>%
  mutate(biological_meaning = factor(biological_meaning, levels = rev(biological_meaning)))

df_markers

df_markers

# 3. 顶刊学术风森林图可视化
ggplot(df_markers, aes(x = median, y = biological_meaning)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.18, color = "#1A5276", linewidth = 0.9) +
  geom_point(color = "#E74C3C", size = 3, shape = 16) +
  geom_text(aes(label = sprintf("%.1f%% (%.1f%%-%.1f%%)", median * 100, lower * 100, upper * 100)),
            hjust = 0.5, vjust = -0.9, size = 3, color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(NA, 1.05)) +
  theme_bw(base_size = 12)+
  labs(x='Sensitivity to RSV infection',y=NULL)

ggsave(paste0('code_03/',MODEL_NAME,'/pi_sensitivity.jpg'), width = 7.5, height = 4.5, dpi = 300)


mcmc_trace(fit_truncat,pars = c('pi[1]','pi[2]'))

marker_labels <- c(
  "pi[1]" = "IgA ~ sensitivity ~ (pi[1])",
  "pi[2]" = "IgG ~ sensitivity ~ (pi[2])"
)

#color_scheme_set("brewer-Spectral")
mcmc_trace(fit_truncat, pars = c('pi[1]', 'pi[2]')) +
  # 💡 核心魔法：强行覆盖原有的分面标签，并开启数学表达式解析 (label_parsed)
  facet_wrap(~ parameter,
             labeller = labeller(parameter = as_labeller(marker_labels, default = label_parsed))) +
  # 使用你最习惯的顶刊低调学术风主题
  theme_bw(base_size = 12)+
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  )+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.minor  = element_blank(),
    axis.line = element_line(linewidth = .2,arrow = arrow(length = unit(0,'mm'))),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )+
  labs(x = "Iteration step",y='Sensitivity for RSV infection')

ggsave(paste0('code_03/',MODEL_NAME,'/marker_mcmc_trace.jpg'), width = 8, height = 4, dpi = 300)


# Risk of infection ----
## Use Mean ----
(MAX_AGE<-ncol(cal_idx_matrix)-1)

WINDOW_MONTH<-case_when(
  MODEL_NAME=='fit_sens3_3'~3,
  MODEL_NAME=='fit_sens3_9'~9,
  MODEL_NAME=='fit_sens3_12'~12,
  T~6
)

get_df_individual_raw <- function(df_sero_sub, df_means, WINDOW_MONTH) {

  # 1. 动态提取 w_* 列的名称前缀，比如 "w_1", "w_2", ... "w_30"
  w_names <- df_means %>%
    filter(str_detect(name, "^w\\[")) %>%
    mutate(m = as.integer(str_extract(name, "\\d+"))) %>%
    pull(m) %>%
    unique() %>%
    sort() %>%
    paste0("w_", .) # 生成形如 c("w_1", "w_2", ..., "w_30") 的字符向量

  df_sero_sub %>%
    mutate(ID = row_number()) %>%
    # 💡 联动修改 1：动态拼入 MCMC 迭代生成的 w 列
    bind_cols(
      df_means %>%
        filter(str_detect(name, "^w\\[")) %>%
        mutate(m = as.integer(str_extract(name, "\\d+"))) %>%
        select(m, value) %>%
        pivot_wider(names_from = m, names_prefix = "w_", values_from = value) %>%
        slice(rep(1, nrow(df_sero_sub)))
    ) %>%
    rowwise() %>%
    mutate(
      age_i = as.integer(age_value),
      #start_m = max(0, age_i - 6), # 满月龄尺度，最小触底到 0
      start_m = max(0, age_i - WINDOW_MONTH), # 满月龄尺度，最小触底到 0
      eff_len = age_i - start_m     # 严格的随访格子数，不加 1
    ) %>%
    filter(eff_len > 0) %>%
    mutate(
      # 🛠️ 1. 上面这个块彻底改写：利用 c_across 动态捕获 w_1 到 w_N
      timeline_data = list({
        # 💡 使用 c_across 配合提前生成的列名向量，动态转成 w_vec 向量（兼容 24 或 30）
        w_vec <- c_across(all_of(w_names))

        h <- numeric(eff_len)
        S_surv_win <- numeric(eff_len + 1)
        S_surv_win[1] <- 1.0
        prob_t <- numeric(eff_len)

        for (idx in 1:eff_len) {
          m <- start_m + idx - 1       # 纯粹的物理满月龄 (0, 1, 2...)
          stan_idx <- m + 1            # 核心安全保护：R 下标从 1 开始 (1, 2, 3...)

          c_idx <- cal_idx_matrix[ID, stan_idx] # 精准锁定当前婴儿 (ID) 和当前月份格子

          # 💡 动态 w_vec 已经自适应长度，直接通过 stan_idx 索引
          h[idx] <- lambda_vec[c_idx] * w_vec[stan_idx]

          S_surv_win[idx + 1] <- S_surv_win[idx] * exp(-h[idx])
          prob_t[idx] <- S_surv_win[idx] * -expm1(-h[idx])
        }

        tibble(
          absolute_age_month = start_m : (age_i - 1),
          prob_infection_at_this_month = prob_t
        )
      }),

      # 🛠️ 2. 下面这个块同步对齐：同样使用动态捕获的 w_vec
      prob_uninfected = {
        w_vec <- c_across(all_of(w_names))

        h_sum <- 0
        for (idx in 1:eff_len) {
          m <- start_m + idx - 1
          stan_idx <- m + 1
          c_idx <- cal_idx_matrix[ID, stan_idx]

          h_sum <- h_sum + lambda_vec[c_idx] * w_vec[stan_idx]
        }
        exp(h_sum * -1) # 相当于 exp(-h_sum)，但写法对一些解析器更友好
      }
    ) %>%
    ungroup()
}

df_individual_raw <- get_df_individual_raw(df_sero_sub,df_means, WINDOW_MONTH)

df_individual_wide<-df_individual_raw %>%
  select(ID, age_value, prob_uninfected, timeline_data) %>%
  unnest(timeline_data) %>%
  mutate(absolute_age_month = factor(absolute_age_month, levels = 0:MAX_AGE)) %>%
  pivot_wider(
    names_from = absolute_age_month,         # 把绝对月份变成新列名
    names_prefix = "prob_infect_age_",       # 给列名加个清爽的前缀
    values_from = prob_infection_at_this_month, # 填充的数值
    values_fill = NA_integer_                       # 不在窗口内的月份自动填 0
  ) %>%
  arrange(age_value) %>%
  relocate(all_of(paste0('prob_infect_age_',max(5-WINDOW_MONTH,0):(MAX_AGE-1))),.after = prob_uninfected) %>%
  rowwise() %>%
  mutate(PROB=sum(c_across(starts_with('prob_')),na.rm=T)) %>%
  ungroup()


export(df_individual_wide,paste0('code_03/',MODEL_NAME,'/df_individual_wide.xlsx'))

df_individual_wide %>%
  summarise(expected_uninfected_prob = mean(prob_uninfected, na.rm = TRUE),.by=age_value) %>%
  ggplot(aes(age_value,1-expected_uninfected_prob))+
  geom_point()+
  geom_line(group=1)+
  scale_x_continuous(breaks = 5:MAX_AGE)

df_individual_wide %>%
  summarise(
    across(starts_with("prob_infect_age_"), ~mean(.x, na.rm = TRUE), .names = "{.col}")
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "expected_prob"
  ) %>%
  mutate(age = as.integer(str_remove(age, "prob_infect_age_"))) %>%
  arrange(age) %>%
  mutate(cum_prob=cumsum(expected_prob)) %>%
  ggplot(aes(x = age, y = expected_prob)) +
  #ggplot(aes(x = age, y = cum_prob)) +
  geom_line(color = "#2471A3", size = 1.2) +
  geom_point(color = "#E74C3C", size = 3) +
  scale_x_continuous(breaks = 0:MAX_AGE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_minimal() +
  labs(
    # title = "Expected Primary RSV Infection Probability by Absolute Age",
    # subtitle = paste("Based on Latent Survival Analysis (N =", nrow(df_individual_wide), ")"),
    x = "Absolute Age (months)",
    y = "Expected probability of RSV infection"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "gray40", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

#ggsave('plot/truncat/expected_prob_infection_mean.jpg', width = 8, height = 4, dpi = 300)
#ggsave('plot/truncat/cum_prob_infection_mean.jpg', width = 8, height = 4, dpi = 300)

## Use Posterior sample----
df_all_draws <- as_draws_df(fit_truncat)

# df_individual_wide_all<-1:4 %>%
#   map_dfr(~{
#     cat('--',.x,'\n')
#     post_id=.x
#     df_pseudo_means<-as.data.frame(df_all_draws[post_id,])  %>% t() %>% enframe() %>%
#       mutate(value = as.numeric(value))
#     df_individual_raw <- get_df_individual_raw(df_sero_sub,df_pseudo_means)
#     df_individual_wide<-df_individual_raw %>%
#       select(ID, age_value, prob_uninfected, timeline_data) %>%
#       unnest(timeline_data) %>%
#       mutate(absolute_age_month = factor(absolute_age_month, levels = 0:23)) %>%
#       pivot_wider(
#         names_from = absolute_age_month,         # 把绝对月份变成新列名
#         names_prefix = "prob_infect_age_",       # 给列名加个清爽的前缀
#         values_from = prob_infection_at_this_month, # 填充的数值
#         values_fill = NA_integer_                       # 不在窗口内的月份自动填 0
#       ) %>%
#       mutate(post_id=post_id,.before = 1)
#   }) %>%
#   relocate(all_of(paste0('prob_infect_age_',0:22)),.after = prob_uninfected) %>%
#   arrange(ID,post_id)

### Fast Cauculate----

get_df_individual_wide_fast <- function(df_sero_sub, df_all_draws, cal_idx_matrix, lambda_vec,WINDOW_MONTH) {

  # 1. 精准、干净地切出 4000 行 w[1] 到 w[24] 的纯数值矩阵
  w_col_indices <- grep("^w\\[", colnames(df_all_draws))
  mat_w <- as.matrix(df_all_draws[, w_col_indices]) # [4000 × 24]

  n_draws <- nrow(mat_w)
  n_infants <- nrow(df_sero_sub)
  age_vec <- as.integer(df_sero_sub$age_value)

  # 2. 预先在内存中一次性开辟 4000 × 530 的大矩阵空间，拒绝碎片化
  mat_uninfected <- matrix(1, nrow = n_draws, ncol = n_infants)

  # 建立 24 个月份的超大空列表，用来存分月概率
  mat_monthly_list <- lapply(0:MAX_AGE, function(x) matrix(0, nrow = n_draws, ncol = n_infants))

  cat("⚡ 矩阵并行计算核心启动...\n")

  # 3. 💡 绝妙的逆向思维：外层循环只跑 530 个婴儿，而不是 4000 次！
  # 因为 530 次循环极其轻量，而内部可以直接利用 R 的矩阵乘法让 4000 行全样本并行爆发！
  for (i in 1:n_infants) {
    age_i <- age_vec[i]
    start_m <- max(0, age_i - WINDOW_MONTH)
    eff_len <- age_i - start_m

    if (eff_len > 0) {
      m_range <- start_m : (age_i - 1)
      stan_idx_vec <- m_range + 1

      # 提取当前宝宝回溯窗口内各月的真实社区流行度
      c_indices <- cal_idx_matrix[i, stan_idx_vec]
      lambdas <- lambda_vec[c_indices] # 长度为 eff_len

      # 💡 核心魔法：让当前宝宝的 lambdas 与 4000 行 w 对应的列直接进行向量化矩阵乘法！
      # mat_w[, stan_idx_vec] 是一张 [4000 × eff_len] 的矩阵
      # sweep 函数让这 4000 行在一瞬间全部乘以 lambdas 向量，速度快到飞起！
      h_matrix <- sweep(mat_w[, stan_idx_vec, drop = FALSE], 2, lambdas, `*`) # [4000 × eff_len]

      # 计算未感染概率（行求和后指数化，一次性得到 4000 个世界的结果）
      h_sum_vec <- rowSums(h_matrix)
      mat_uninfected[, i] <- exp(-h_sum_vec)

      # 递推计算每个历史月份的初次感染概率
      # 利用累计求和矩阵并行递推生存状态
      if (eff_len == 1) {
        mat_monthly_list[[stan_idx_vec[1]]][, i] <- -expm1(-h_matrix[, 1])
      } else {
        # 算出每一格左侧的累积生存概率矩阵
        h_cum <- t(apply(h_matrix, 1, cumsum))
        S_surv <- cbind(1, exp(-h_cum[, -eff_len, drop = FALSE])) # [4000 × eff_len]

        prob_t_matrix <- S_surv * (-expm1(-h_matrix))

        # 完美的把这 4000 行结果批量倒进对应的月份大账本里
        for (idx in 1:eff_len) {
          curr_m_idx <- stan_idx_vec[idx]
          mat_monthly_list[[curr_m_idx]][, i] <- prob_t_matrix[, idx]
        }
      }
    }
  }

  cat("📊 核心计算完成！正在进行高速扁平化网格组装...\n")

  # 4. 利用 expand.grid 瞬间平铺出 4000 × 530 = 2,120,000 行的大骨架
  df_base_grid <- expand.grid(
    post_id = 1:n_draws,
    ID = 1:n_infants
  ) %>%
    mutate(age_value = age_vec[ID])

  # 5. 矩阵直接降维拍扁（as.vector），一秒填充满所有数据列
  df_base_grid$prob_uninfected <- as.vector(mat_uninfected)

  for(m in max(5-WINDOW_MONTH,0):MAX_AGE) {
    col_name <- paste0("prob_infect_age_", m)
    df_base_grid[[col_name]] <- as.vector(mat_monthly_list[[m + 1]])
  }

  # 6. 重定位并完美排序
  df_final_wide <- df_base_grid %>%
    relocate(all_of(paste0('prob_infect_age_', max(5-WINDOW_MONTH,0):(MAX_AGE-1))), .after = prob_uninfected) %>%
    arrange(ID, post_id)

  return(df_final_wide)
}
# 💡 以前需要跑几个小时，现在一敲回车，十几秒直接收网！
df_individual_wide_all <- get_df_individual_wide_fast(
  df_sero_sub    = df_sero_sub,
  df_all_draws   = df_all_draws,
  cal_idx_matrix = cal_idx_matrix,
  lambda_vec     = lambda_vec,
  WINDOW_MONTH   = WINDOW_MONTH
)

df_individual_wide_all %>%
  summary()


rio::export(df_individual_wide_all,paste0('code_03/',MODEL_NAME,'/df_individual_wide_all.rds'))


