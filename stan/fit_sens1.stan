// RSV 初次感染年龄推断模型 - 基于横断面血清学数据的潜变量生存分析
// 版本: v5.3 (2026年7月 - 满月龄0起点 - 拓宽至 0-29 月龄安全版)

data {
  int<lower=1> N;                // 样本量
  int<lower=1> M;                // 流行强度序列的总月数
  int<lower=5, upper=29> A[N];   // 💡 修改：上限扩展至 29 (对应小于30个月)

  // 💡 修改：列数拓宽为 30
  int<lower=1, upper=M> cal_idx[N, 30];

  // 观测指标数据
  vector[N] Y_IgA;
  vector[N] Y_IgG;
  matrix[N, 5] Y_epi;            // 5个表位特异性抑制率

  // 外部已知量
  vector[M] lambda;              // 全局按月相对流行强度 (0-1)

  // 💡 修改：B-Spline 基函数预计算矩阵 (30行对应满月龄 0:29)
  int<lower=1> J;
  matrix[30, J] B_basis;

  // 💡 新增：外部输入的先验中心值
  real<lower=0> prior_r_mat_mean;
  real<lower=0> prior_r1_mean;
  real<lower=0> prior_r2_mean;
  vector<lower=0>[5] prior_r_epi_mean;

  // 💡 参数重构：按抗体类型划分先验标准差 (用于敏感性分析)
  real<lower=0> prior_sd_main;   // 作用于主抗体: r_mat, r1, r2
  real<lower=0> prior_sd_epi;    // 作用于 5 个表位: r_epi[1:5]
}

parameters {
  // 1. 年龄别易感性权重参数 (30个格子由 J 个基底控制)
  vector<lower=0>[J] beta;
  real<lower=0> tau;             // 平滑惩罚超参数

  // 2. 诊断敏感度
  vector<lower=0, upper=1>[7] pi;

  // 3. 感染后抗体峰值增量
  vector<lower=0>[7] peak;

  real<lower=0> Mat0;            // 出生时母源 IgG
  vector<lower=0>[5] Mat0_k;     // 出生时母源对各表位的贡献

  // 5. 背景与基线
  real<lower=0> B_IgA;
  real<lower=0> B_IgG;
  vector<lower=0, upper=1>[5] B_epi;

  // 6. 观测误差与精度
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector<lower=0>[5] phi;        // Beta 分布精度参数

  // 💡 将衰减速率改为待估参数（加 lower=0 约束）
  real<lower=0> r_mat;  // 母源抗体衰减速率
  real<lower=0> r1;     // IgA 衰减速率
  real<lower=0> r2;     // IgG 衰减速率
  vector<lower=0>[5] r_epi; // 表位抗体衰减速率
}

transformed parameters {
  vector[30] w;                  // 💡 修改：长度拓宽为 30
  // 生成年龄权重曲线 (对应满月龄 0, 1, 2 ... 29)
  w = B_basis * beta;
}

model {
  // --- 使用重构后的 SD 参数配置先验 ---
  // 主抗体 (r_mat, r1, r2) 共享 prior_sd_main
  r_mat ~ lognormal(log(prior_r_mat_mean), prior_sd_main);
  r1    ~ lognormal(log(prior_r1_mean),    prior_sd_main);
  r2    ~ lognormal(log(prior_r2_mean),    prior_sd_main);

  // 5 个表位抗体共享 prior_sd_epi
  for (k in 1:5) {
    r_epi[k] ~ lognormal(log(prior_r_epi_mean[k]), prior_sd_epi);
  }

  // --- 先验部分 ---
  beta ~ normal(0, 10);
  for (j in 2:(J-1)) {
    target += normal_lpdf(beta[j-1] - 2*beta[j] + beta[j+1] | 0, tau);
  }
  tau ~ exponential(1);

  pi ~ beta(8, 2);
  peak ~ lognormal(2, 1);
  Mat0 ~ lognormal(5, 1);
  sigma1 ~ exponential(1);
  sigma2 ~ exponential(1);
  phi ~ exponential(0.1);

  B_IgA ~ lognormal(log(5), 0.5);
  B_IgG ~ lognormal(log(50), 0.5);
  B_epi ~ beta(0.05, 0.95);

  // --- 似然函数 ---
  // 💡 注意：由于使用了动态月龄 A[i] 边界，以下部分无需修改，自动兼容 0-29 尺度
  for (i in 1:N) {
    int window_size = 6;

    int start_m = max(0, A[i] - window_size);
    int eff_len = A[i] - start_m;

    if (eff_len > 0) {
      vector[eff_len] h;
      vector[eff_len + 1] S_surv_win;
      vector[eff_len] prob_t;
      vector[eff_len] log_lik_t;

      real prob_S0;
      real log_lik_S0;

      S_surv_win[1] = 1.0;

      for (idx in 1:eff_len) {
        int m = start_m + idx - 1;
        int stan_idx = m + 1; // 自动映射到 1:30 的安全数组下标

        h[idx] = lambda[cal_idx[i, stan_idx]] * w[stan_idx];
        S_surv_win[idx + 1] = S_surv_win[idx] * exp(-h[idx]);
        prob_t[idx] = S_surv_win[idx] * -expm1(-h[idx]);
      }

      prob_S0 = S_surv_win[eff_len + 1];

      // 情况 A: 窗口内未感染 (S=0)
      real mu_IgA_S0 = B_IgA;
      real mu_IgG_S0 = B_IgG + Mat0 * exp(-r_mat * A[i]);

      log_lik_S0 = log(prob_S0) +
                   lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S0), sigma1) +
                   lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2);

      for(k in 1:5) {
        real mu_epi_S0 = inv_logit(logit(B_epi[k]) + Mat0_k[k] * exp(-r_mat * A[i]));
        log_lik_S0 += beta_lpdf(Y_epi[i, k] | mu_epi_S0 * phi[k], (1-mu_epi_S0) * phi[k]);
      }

      // 情况 B: 窗口内已感染 (S=1)
      for (idx in 1:eff_len) {
        int m = start_m + idx - 1;
        real dt = A[i] - m - 1;

        // --- IgA 计算 ---
        real mu_IgA_S1 = mu_IgA_S0 + peak[1] * exp(-r1 * dt);
        log_lik_t[idx] = log_sum_exp(
          log(pi[1]) + lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S1), sigma1),
          log(1 - pi[1]) + lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S0), sigma1)
        );

        // --- IgG 计算 ---
        real mu_IgG_S1 = mu_IgG_S0 + peak[2] * exp(-r2 * dt);
        log_lik_t[idx] += log_sum_exp(
          log(pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S1), sigma2),
          log(1 - pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2)
        );

        // --- 表位计算 ---
        for (k in 1:5) {
          real mu_epi_S0_t = inv_logit(logit(B_epi[k]) + Mat0_k[k] * exp(-r_mat * A[i]));
          real mu_epi_S1_t = inv_logit(logit(mu_epi_S0_t) + peak[k+2] * exp(-r_epi[k] * dt));
          log_lik_t[idx] += log_sum_exp(
            log(pi[k+2]) + beta_lpdf(Y_epi[i, k] | mu_epi_S1_t * phi[k], (1-mu_epi_S1_t) * phi[k]),
            log(1 - pi[k+2]) + beta_lpdf(Y_epi[i, k] | mu_epi_S0_t * phi[k], (1-mu_epi_S0_t) * phi[k])
          );
        }

        log_lik_t[idx] += log(prob_t[idx]);
      }

      target += log_sum_exp(log_lik_S0, log_sum_exp(log_lik_t));
    }
  }
}
