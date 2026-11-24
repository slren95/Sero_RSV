// fit_sens5_new: IgA and IgG only; all three waning rates estimated with lognormal priors.
// RSV 初次感染年龄推断模型 - 基于横断面血清学数据的潜变量生存分析
// 版本: v6.0 (双抗体 IgA + IgG 极简优化版 - 完整代码)

// 与 stan_new/fit_main_new.stan 对照：
// 1. 仅保留 IgA/IgG 观测似然，移除 Y_epi 和五个表位的观测信息。
// 2. pi、peak 的长度由 7 改为 2；移除 Mat0_k、B_epi、phi、r_epi 及其相关先验。
// 3. 移除表位先验输入 prior_meanlog_epi/prior_sdlog_epi；对应 R 脚本无需读取 epi_prior_stan。
// 4. r_mat、r1、r2 仍按主模型使用 sdlog = 0.10 的对数正态先验，其余保留部分的统计设定相同。

data {
  int<lower=1> N;                // 样本量
  int<lower=1> M;                // 流行强度序列的总月数
  int<lower=5, upper=29> A[N];   // 满月龄 (对应 0-29，上限 29)

  int<lower=1, upper=M> cal_idx[N, 30]; // 对应 30 个月龄的日历月份索引

  // 观测指标数据 (仅保留 IgA 与 IgG)
  vector[N] Y_IgA;
  vector[N] Y_IgG;

  // 外部已知量
  vector[M] lambda;              // 全局按月相对流行强度 (0-1)
  real<lower=0> prior_r_mat;
  real<lower=0> prior_r1;
  real<lower=0> prior_r2;

  // B-Spline 基函数预计算矩阵 (30行对应满月龄 0:29)
  int<lower=1> J;
  matrix[30, J] B_basis;
}

parameters {
  real<lower=0> r_mat;
  real<lower=0> r1;
  real<lower=0> r2;
  // 1. 年龄别易感性权重参数 (30个格子由 J 个基底控制)
  vector<lower=0>[J] beta;
  real<lower=0> tau;              // 平滑惩罚超参数

  // 2. 诊断敏感度 (仅保留 IgA 和 IgG 两个指标)
  vector<lower=0, upper=1>[2] pi; // pi[1]: IgA, pi[2]: IgG

  // 3. 感染后抗体峰值增量
  vector<lower=0>[2] peak;       // peak[1]: IgA, peak[2]: IgG

  // 4. 母源抗体
  real<lower=0> Mat0;            // 出生时母源 IgG

  // 5. 背景与基线
  real<lower=0> B_IgA;
  real<lower=0> B_IgG;

  // 6. 观测误差
  real<lower=0> sigma1;          // IgA 测量残差
  real<lower=0> sigma2;          // IgG 测量残差
}

transformed parameters {
  vector[30] w;                  // 生成年龄权重曲线 (对应满月龄 0, 1, 2 ... 29)
  w = B_basis * beta;
}

model {
  r_mat ~ lognormal(log(prior_r_mat), 0.10);
  r1 ~ lognormal(log(prior_r1), 0.10);
  r2 ~ lognormal(log(prior_r2), 0.10);
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

  B_IgA ~ lognormal(log(5), 0.5);
  B_IgG ~ lognormal(log(50), 0.5);

  // --- 似然函数 ---
  for (i in 1:N) {
    int window_size = 6;

    int start_m = max(0, A[i] - window_size);
    int eff_len = A[i] - start_m;

    if (eff_len > 0) {
      vector[eff_len] h;
      vector[eff_len + 1] S_surv_win;
      vector[eff_len] prob_t;
      vector[eff_len] log_lik_t;

      real log_lik_S0;

      S_surv_win[1] = 1.0;

      for (idx in 1:eff_len) {
        int m = start_m + idx - 1;
        int stan_idx = m + 1; // 自动映射到 1:30 的安全数组下标

        h[idx] = lambda[cal_idx[i, stan_idx]] * w[stan_idx];
        S_surv_win[idx + 1] = S_surv_win[idx] * exp(-h[idx]);
        prob_t[idx] = S_surv_win[idx] * -expm1(-h[idx]);
      }

      real prob_S0 = S_surv_win[eff_len + 1];

      // ==================== 情况 A: 窗口内未感染 (S=0) ====================
      real mu_IgA_S0 = B_IgA;
      real mu_IgG_S0 = B_IgG + Mat0 * exp(-r_mat * A[i]);

      log_lik_S0 = log(prob_S0)
                 + lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S0), sigma1)
                 + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2);

      // ==================== 情况 B: 窗口内已感染 (S=1) ====================
      for (idx in 1:eff_len) {
        int m = start_m + idx - 1;
        real dt = A[i] - m - 1;

        real mu_IgA_S1 = mu_IgA_S0 + peak[1] * exp(-r1 * dt);
        real mu_IgG_S1 = mu_IgG_S0 + peak[2] * exp(-r2 * dt);

        // IgA 混合项
        real log_lik_IgA = log_sum_exp(
          log(pi[1]) + lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S1), sigma1),
          log(1 - pi[1]) + lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S0), sigma1)
        );

        // IgG 混合项
        real log_lik_IgG = log_sum_exp(
          log(pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S1), sigma2),
          log(1 - pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2)
        );

        log_lik_t[idx] = log_lik_IgA + log_lik_IgG + log(prob_t[idx]);
      }

      target += log_sum_exp(log_lik_S0, log_sum_exp(log_lik_t));
    }
  }
}


