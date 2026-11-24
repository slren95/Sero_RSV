 # stan/ 模型文件说明

本项目所有模型都是同一套核心框架的不同版本/变体：
**基于横断面血清学的 RSV 初次感染年龄潜变量生存模型** —— 离散时间风险 `h(m)=lambda[cal_idx]*w[m]`，
`w` 由 B 样条基底 `B_basis * beta` 生成；对"窗口内未感染 (S0)"与"某月龄 m 感染 (S1)"两种潜在状态做混合边际化
（IgA/IgG 用对数正态，5 个表位抑制率用 Beta，全部经 log_sum_exp 汇总）。

## 文件总览

| 文件 | 版本/用途 | 状态 |
|---|---|---|
| `fit_truncat_extend.stan` | v5.3（2026-07），0–29 月龄主模型 | ✅ 现行基线 |
| `fit_truncat_extend_sens3.stan` | v5.4（2026-09），window_size 与 MAX_AGE 动态化 | ✅ 现行（敏感性用） |
| `fit_truncat.stan` | v5.2（2026-06），0–23 月龄 | 🔶 历史版（24 列网格版） |
| `rsv_model_truncat_old.stan` | v5.2，0–23 月龄 | 🔶 历史版（含废弃的背景扰动重参数化） |
| `fit_truncat_extend_sens1.stan` | v5.3 + 衰减率改为待估参数 | 🔶 敏感性分析 1 |
| `rsv_model.stan` | v4.1（2026-04），原始全史模型 | ⛔ 废弃 |
| `fit_truncat_extend_sens3_3/9/12.stan`、`fit_truncat_extend_sens4.stan` | 空占位文件（utils.R 脚手架生成） | ⛔ 空文件，无代码 |

## 主要异同

### 共同点（所有非空文件）
- 参数集一致：`beta`（易感性权重，B 样条）+ `tau`（二阶差分平滑惩罚）、`pi[7]`（7 个指标诊断敏感度）、
  `peak[7]`（感染后抗体峰值增量）、`Mat0`/`Mat0_k`（母源抗体）、`B_IgA/B_IgG/B_epi`（背景基线）、
  `sigma1/sigma2/phi`（观测误差）；
- 先验一致：`beta~N(0,10)`、`pi~beta(8,2)`、`peak~lognormal(2,1)`、`Mat0~lognormal(5,1)`、
  `tau/sigma~exp(1)`、`phi~exp(0.1)`；
- 似然结构一致：逐个体在采样月龄前的小窗口内递推条件生存概率，S0/S1 两分支做 log_sum_exp 边际化。

### 差异点
| 维度 | rsv_model (v4.1) | rsv_model_truncat_old / fit_truncat (v5.2) | fit_truncat_extend (v5.3) | sens1 | sens3 (v5.4) |
|---|---|---|---|---|---|
| 月龄尺度/列数 | 采样月龄 1–23，`cal_idx[N,23]` | 满月龄 0 起点，`cal_idx[N,24]` | 满月龄 0 起点，**拓宽至 0–29**，`cal_idx[N,30]` | 同 extend | `MAX_AGE` 由 data 传入，`cal_idx[N,MAX_AGE]` |
| 追溯窗口 | 全史（无窗口截断） | 硬编码 `window_size=6` | 硬编码 `window_size=6` | 同 extend | `window_size` 由 data 传入（3/6/9/12） |
| IgA 敏感度 `pi[1]` | 锚定 =1（只估 6 个） | 7 个全部估计 | 7 个全部估计 | 同 extend | 同 extend |
| 衰减率 r1/r2/r_epi/r_mat | r1/r2/r_epi 估计，r_mat 固定 | 全部固定为外部 data | 全部固定为外部 data | **改为待估参数**，用 data 传入的均值做对数正态信息先验 | 全部固定 |
| 背景 B 参数化 | 直接参数 | `_old` 版用 base+raw 扰动重参数化；`fit_truncat` 为直接参数 | 直接参数 | 直接参数 | 直接参数 |
| dt 定义（抗体延迟假设） | `dt=A-m` | `dt=A-m-1`（m 月感染、m+1 月抗体达峰） | 同 v5.2 | 同 extend | 同 v5.2 |

### 关键区别解读
- **v4.1 → v5.2**：改为"满月龄 0 起点"并引入 **追溯窗口截断 + 抗体延迟显现假设**（`dt=A-m-1`），
  是 v5.2 之后的系列与旧全史模型最本质的区别；
- **v5.2 → v5.3**：仅把年龄上限从 23 拓宽到 29（24 列→30 列），似然主体无需改动；
- **sens1**：唯一把抗体衰减率 `r_mat/r1/r2/r_epi` 从"外部固定值"改为"待估参数"的模型
  （先验中心值来自纵向数据 `DF_t12`，经 data 传入）；
- **sens3**：唯一不硬编码窗口与年龄上界的模型，专为 `window_size ∈ {3,6,9,12}` 敏感性扫描设计；
- **sens2 / simple 变体没有独立 .stan 文件**：它们复用 `fit_truncat_extend.stan`，
  仅在 R 端改变样条内节点（sens2：`c(6,12,18,24)`）或表位衰减取值（simple：与 IgG 相同，不读 `DF_t12`）。

## 用法建议
- 新分析以 `fit_truncat_extend.stan`（基线）和 `fit_truncat_extend_sens3.stan`（扫 window）为准；
- `fit_truncat_extend_sens3_*.stan` 与 `sens4.stan` 为空占位，勿直接编译；
- 配对 R 脚本见 `code_02/code_02_*.R`，仍从仓库根目录运行（编译调用多被注释，rda/ 中已有拟合结果可直接复用）。
