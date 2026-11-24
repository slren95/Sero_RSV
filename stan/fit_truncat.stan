// RSV 初次感染年龄推断模型 - 基于横断面血清学数据的潜变量生存分析
// 版本: v5.2 (2026年6月 - 满月龄0起点 + 抗体延迟显现安全版)

data {
  int<lower=1> N;                // 样本量
  int<lower=1> M;                // 流行强度序列的总月数
  int<lower=5, upper=23> A[N];   // 每个婴儿采样时的满月龄 (0-23尺度)

  // 索引矩阵：列数拓宽为 24，代表满 0 月龄到满 23 月龄的历史足迹
  int<lower=1, upper=M> cal_idx[N, 24];

  // 观测指标数据
  vector[N] Y_IgA;
  vector[N] Y_IgG;
  matrix[N, 5] Y_epi;            // 5个表位特异性抑制率

  // 外部已知量
  vector[M] lambda;              // 全局按月相对流行强度 (0-1)
  real r_mat;                    // 母源抗体月衰减速率
  real r1;                       // IgA 衰减
  real r2;                       // IgG 衰减
  vector[5] r_epi;               // 表位衰减

  // B-Spline 基函数预计算矩阵 (24行对应满月龄 0:23)
  int<lower=1> J;
  matrix[24, J] B_basis;
}

parameters {
  // 1. 年龄别易感性权重参数 (24个格子由 J=8 个基底控制)
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
}

transformed parameters {
  vector[24] w;
  // 生成年龄权重曲线 (对应满月龄 0, 1, 2 ... 23)
  w = B_basis * beta;
}

model {
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

  B_IgA ~ lognormal(log(5), 0.5); // 背景 MFI 约 5 [cite: 165]
  B_IgG ~ lognormal(log(50), 0.5); // 基线 MFI 约 50 [cite: 168]
  B_epi ~ beta(0.05, 0.95);

  // --- 似然函数 ---
  for (i in 1:N) {
    int window_size = 6;

    // start_m 代表回溯的起始满月龄，最小取到 0
    int start_m = max(0, A[i] - window_size);

    // 💡 正确：去掉 +1！eff_len 严格代表区间内的格子数
    int eff_len = A[i] - start_m;

    if (eff_len > 0) {
      vector[eff_len] h;
      vector[eff_len + 1] S_surv_win;
      vector[eff_len] prob_t;
      vector[eff_len] log_lik_t;

      real prob_S0;
      real log_lik_S0;

      S_surv_win[1] = 1.0;

      // 递推条件生存概率
      for (idx in 1:eff_len) {
        // m 此时为纯粹的生物学【满月龄】(最小为 0，代表出生那月)
        int m = start_m + idx - 1;

        // 💡 核心保护机制：转换为从 1 开始的 Stan 数组安全下标
        int stan_idx = m + 1;

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
        int m = start_m + idx - 1; // 当前模拟感染发生时的【满月龄】(0, 1, ... 23)

        // 💡 核心生物学假设：在 m 月龄感染，抗体在 m+1 月龄到顶并衰减
        // 间隔时间 dt = A[i] - (m + 1)。
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

        // 加上当前历史月份的初次感染概率贡献
        log_lik_t[idx] += log(prob_t[idx]);
      }

      // 【安全混合同步】
      target += log_sum_exp(log_lik_S0, log_sum_exp(log_lik_t));
    }
  }
}
