// RSV 初次感染年龄推断模型 - 基于横断面血清学数据的潜变量生存分析
// 版本: v4.1 (2026年4月)

data {
  int<lower=1> N;                // 样本量 [cite: 19]
  int<lower=1> M;                // 流行强度序列的总月数
  int<lower=5, upper=23> A[N];   // 每个婴儿采样时的月龄 [cite: 20]

  // 索引矩阵：婴儿 i 在月龄 m 时，对应 lambda 序列中的绝对月份索引
  int<lower=1, upper=M> cal_idx[N, 23];

  // 观测指标数据 [cite: 20]
  vector[N] Y_IgA;
  vector[N] Y_IgG;
  matrix[N, 5] Y_epi;            // 5个表位特异性抑制率 (Site III, II, V, I, IV) [cite: 21-26]

  // 外部已知量 [cite: 30]
  vector[M] lambda;              // 全局按月相对流行强度 (0-1)
  real r_mat;                    // 母源抗体月衰减速率 (固定值 0.58) [cite: 32]

  // B-Spline 基函数预计算矩阵 (由 R 预计算，J=8 对应 4个内部节点) [cite: 36, 48]
  int<lower=1> J;
  matrix[24, J] B_basis;
}

parameters {
  // 1. 年龄别易感性权重参数 [cite: 36, 42]
  vector<lower=0>[J] beta;
  real<lower=0> tau;             // 平滑惩罚超参数 [cite: 36, 52]

  // 2. 诊断敏感度 (pi[1] 固定为 1) [cite: 36, 37, 149]
  vector<lower=0, upper=1>[6] pi_rest;

  // 3. 感染后抗体峰值增量 [cite: 36, 154]
  vector<lower=0>[7] peak;

  // 4. 抗体衰减速率 [cite: 32, 177-182]
  real<lower=0> r1;              // IgA 衰减
  real<lower=0> r2;              // IgG 衰减
  vector<lower=0>[5] r_epi;      // 表位衰减 (层次先验)

  // 5. 背景与基线 [cite: 36, 163-171]
  real<lower=0> B_IgA;
  real<lower=0> B_IgG;
  vector<lower=0, upper=1>[5] B_epi;
  real<lower=0> Mat0;            // 出生时母源 IgG [cite: 172]
  vector<lower=0>[5] Mat0_k;     // 出生时母源对各表位的贡献 [cite: 175]

  // 6. 观测误差与精度 [cite: 36, 157-162]
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector<lower=0>[5] phi;        // Beta 分布精度参数
}

transformed parameters {
  vector[7] pi;
  vector[24] w;

  // 锚定约束：IgA 敏感度设为 1 [cite: 37, 149]
  pi[1] = 1.0;
  for(k in 2:7) pi[k] = pi_rest[k-1];

  // 生成年龄权重曲线 [cite: 42, 189]
  w = B_basis * beta;
}

model {
  // --- 先验部分 ---
  beta ~ normal(0, 10);
  for (j in 2:(J-1)) {
    // 二阶差分平滑惩罚 [cite: 55, 143]
    target += normal_lpdf(beta[j-1] - 2*beta[j] + beta[j+1] | 0, tau);
  }
  tau ~ exponential(1);          // [cite: 52, 146]

  pi_rest ~ beta(8, 2);          // 敏感度均值 0.8 [cite: 152]
  peak ~ lognormal(2, 1);        // 峰值弱信息先验 [cite: 155]

  r1 ~ lognormal(log(0.15), 0.5); // [cite: 179]
  r2 ~ lognormal(log(0.07), 0.5); // [cite: 179]
  r_epi ~ lognormal(-2.0, 0.5);   // 层次先验

  B_IgA ~ lognormal(log(5), 0.5); // 背景 MFI 约 5 [cite: 165]
  B_IgG ~ lognormal(log(50), 0.5); // 基线 MFI 约 50 [cite: 168]
  Mat0 ~ lognormal(5, 1);

  sigma1 ~ exponential(1);
  sigma2 ~ exponential(1);
  phi ~ exponential(0.1);

  // --- 似然函数 ---
  for (i in 1:N) {
    vector[A[i]] h;
    vector[A[i] + 1] S_surv;
    vector[A[i]] prob_t;
    real prob_S0;

    // 1. 计算每个月的感染风险与生存概率 [cite: 56, 64, 70]
    S_surv[1] = 1.0;
    for (m in 1:A[i]) {
      // 瞬时风险 = 流行强度 * 年龄权重 [cite: 56, 189]
      h[m] = lambda[cal_idx[i, m]] * w[m];
      S_surv[m+1] = S_surv[m] * exp(-h[m]);
      // 使用 expm1 提高数值精度
      prob_t[m] = S_surv[m] * -expm1(-h[m]);
    }
    prob_S0 = S_surv[A[i]+1];   // 采样时未感染的概率 [cite: 67, 189]

    // 2. 边际似然汇总 [cite: 125-138]

    // 情况 A: 采样时未感染 (S=0) [cite: 93-106]
    real log_lik_S0;
    real mu_IgG_S0 = B_IgG + Mat0 * exp(-r_mat * A[i]); // [cite: 100, 189]

    log_lik_S0 = log(prob_S0) +
                 lognormal_lpdf(Y_IgA[i] | log(B_IgA), sigma1) +
                 lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2);

    for(k in 1:5) {
      // 表位未感染期望值计算 [cite: 106, 189]
      real mu_epi_S0 = inv_logit(logit(B_epi[k]) + Mat0_k[k] * exp(-r_mat * A[i]));
      log_lik_S0 += beta_lpdf(Y_epi[i, k] | mu_epi_S0 * phi[k], (1-mu_epi_S0) * phi[k]);
    }

    // 情况 B: 采样前已感染 (S=1)，对所有可能的首次感染时间 t 进行边际化 [cite: 107-124]
    vector[A[i]] log_lik_t;
    for (m in 1:A[i]) {
      real dt = A[i] - m; // 感染至今的时长 [cite: 112]

      // 指标 1: IgA (pi_1=1) [cite: 114, 124]
      real mu_IgA_S1 = B_IgA + peak[1] * exp(-r1 * dt);
      log_lik_t[m] = lognormal_lpdf(Y_IgA[i] | log(mu_IgA_S1), sigma1);

      // 指标 2: IgG (混合模型) [cite: 117, 126]
      real mu_IgG_S1 = mu_IgG_S0 + peak[2] * exp(-r2 * dt);
      log_lik_t[m] += log_sum_exp(
        log(pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S1), sigma2),
        log(1 - pi[2]) + lognormal_lpdf(Y_IgG[i] | log(mu_IgG_S0), sigma2)
      );

      // 指标 3-7: 各表位 (混合模型) [cite: 120, 126, 189]
      for (k in 1:5) {
        real mu_epi_S0_t = inv_logit(logit(B_epi[k]) + Mat0_k[k] * exp(-r_mat * A[i]));
        real mu_epi_S1_t = inv_logit(logit(mu_epi_S0_t) + peak[k+2] * exp(-r_epi[k] * dt));
        log_lik_t[m] += log_sum_exp(
          log(pi[k+2]) + beta_lpdf(Y_epi[i, k] | mu_epi_S1_t * phi[k], (1-mu_epi_S1_t) * phi[k]),
          log(1 - pi[k+2]) + beta_lpdf(Y_epi[i, k] | mu_epi_S0_t * phi[k], (1-mu_epi_S0_t) * phi[k])
        );
      }
      log_lik_t[m] += log(prob_t[m]); // 乘以在该时间点感染的概率分布 [cite: 133]
    }

    // 最终汇总个体对数似然：Log-Sum-Exp 技巧处理混合分布 [cite: 138]
    target += log_sum_exp(log_lik_S0, log_sum_exp(log_lik_t));
  }
}
