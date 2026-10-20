# threeSGE

`threeSGE` 将二元网络与加权/概率网络下的 ISS 子图数量估计、VISS
方差估计和 TSGE 模体检验整理为统一的 R 包接口。

## 安装与帮助

从 GitHub 安装：

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("Sergio-Jensen/threeSGE")
library(threeSGE)
```

也可以从本地构建的源码包安装并查看帮助：

```r
install.packages("threeSGE_0.1.0.tar.gz", repos = NULL, type = "source")
library(threeSGE)

help(package = "threeSGE")
?est.subgraph
?est.variance
?tsge
```

## 输入数据约定

- 主网络和模体均可使用方阵；安装 `igraph` 或 `tidygraph` 后也可以传入
  `igraph`/`tbl_graph` 对象。
- `H` 是模体邻接矩阵：非零值表示边，对角线会被忽略。
- 无向分析建议始终传入对称、对角线为 0 的矩阵。
- 二元网络的非对角元素应为 0/1；加权网络必须非负。
- 当 `normalize_weight = TRUE` 时，使用
  `P[i,j] = W[i,j] / max(W)`；为 `FALSE` 时，输入将直接视为 `[0,1]`
  概率矩阵。
- `network_type = "auto"` 根据非对角元素是否全部为 0/1 自动选择二元
  bitset 或加权 ISS。必要时可用 `"binary"` 或 `"weighted"` 强制指定。

## 1. `est.subgraph()`：子图数量或期望数量

```r
est.subgraph(G, H,
             network_type = c("auto", "binary", "weighted"),
             normalize_weight = TRUE, ...)
```

主要参数：

| 参数 | 含义与要求 |
|---|---|
| `G` | 主网络方阵或图对象；二元网络表示观测图，加权网络表示权重或边概率。 |
| `H` | 目标模体方阵或图对象，节点数不得超过 `G`。 |
| `network_type` | `"auto"`、`"binary"` 或 `"weighted"`。 |
| `normalize_weight` | 是否把加权网络除以最大权重；仅影响 weighted 分支。 |
| `B` | ISS 抽样次数，正整数，默认 10000。 |
| `directed` | 是否为有向网络；binary 默认 `FALSE`，weighted 默认 `TRUE`。 |
| `seed` | 可复现随机种子。 |
| `order` | 模体节点编号的一个排列，用于指定贯序。 |
| `order_method` | binary 专用：`"mot.order"`、`"reverse"` 或 `"random"`。 |
| `keep_mappings` | 是否保留成功抽样得到的节点映射。 |
| `verbose` | 是否打印抽样进度。 |
| `aut_count` | 已知的模体自同构数量；`NULL` 时自动计算。 |
| `bootstrap_ci` | binary 专用，是否计算非参数 bootstrap 区间。 |
| `bootstrap_B` | bootstrap 重抽样次数，默认 1000。 |
| `bootstrap_seed` | bootstrap 随机种子。 |
| `verbose.boot` | 是否打印 bootstrap 进度。 |

主要输出：

- 通用：`B`、`directed`、`graph_nodes`、`network_type`、抽样贯序、
  `aut.morph`、`successes` 和 `success_rate`。
- binary：`labeled_estimate`、`unlabeled_estimate` 及相应 `*_se`；Wald、
  log-Wald 和可选 bootstrap 区间；`weight_diagnostics`；原始 `weights`；
  可选 `mappings`。
- weighted：`labeled_expectation`、`unlabeled_expectation` 及相应 SE；
  `probability_matrix`、`probability_density`、`probability_rule`、
  `motif_order` 和 `log_weights`。

```r
H <- matrix(1, 3, 3); diag(H) <- 0
G <- matrix(1, 5, 5); diag(G) <- 0
est.subgraph(G, H, B = 500, directed = FALSE,
             seed = 1, verbose = FALSE)

P <- matrix(0.5, 5, 5); diag(P) <- 0
est.subgraph(P, H, B = 500, directed = FALSE,
             normalize_weight = FALSE, seed = 1, verbose = FALSE)
```

## 2. `iss.sensitivity()`：贯序敏感性分析

参数：`G`、`H`、`B`、`directed`、`seed` 和 `aut_count` 含义同上；
`order` 指定基准贯序；`orders` 可包含 `"mot.order"`、`"reverse"`、
`"given"`；`random_orders` 指定附加随机贯序数量；`verbose` 控制每次 ISS
进度；`bootstrap_ci`、`bootstrap_B` 和 `verbose.boot` 控制各次估计的
bootstrap；`keep_results` 决定是否保留每个完整结果；`verbose.sens` 控制
贯序比较进度。

输出 `sensitivity_table`（各贯序的估计、SE、区间、成功率、ESS、权重 CV
和最大权重份额）、可选的 `results`，以及抽样规模和网络密度等元数据。

```r
iss.sensitivity(G, H, B = 200, random_orders = 2,
                seed = 2, verbose.sens = FALSE)
```

## 3. `est.pairunion()`：FPUG/二阶重叠结构

```r
est.pairunion(X, H, B = 10000, directed = FALSE, induced = TRUE,
              normalize_weight = FALSE, overlap_sizes = NULL,
              mode = c("count", "expectation"),
              allocation = c("equal", "proportional"), ...)
```

| 参数 | 含义与要求 |
|---|---|
| `X` | `mode="count"` 时为二元观测图；`mode="expectation"` 时为概率/权重矩阵。 |
| `H` | 目标模体。 |
| `B` | pairwise-union ISS 总抽样次数。 |
| `directed` | 是否有向。 |
| `induced` | 是否同时施加模体非边约束。 |
| `normalize_weight` | expectation 模式下是否将权重归一化。 |
| `overlap_sizes` | 两个模体共享节点数；默认 `2:(k-1)`。 |
| `mode` | `"count"` 或 `"expectation"`。 |
| `allocation` | 各重叠层采用 `"equal"` 或 `"proportional"` 抽样分配。 |
| `mu_labeled`, `mu_unlabeled` | 可选的已知均值，用于标准化 clumping 指标。 |
| `aut_count` | 可选的已知自同构数量。 |
| `include_self_overlap` | 是否纳入完全重叠 `s=k`；二阶 union 默认不纳入。 |
| `seed` | 随机种子。 |
| `verbose`, `verbose.B` | 分别控制重叠层和层内抽样进度。 |
| `preprocessed` | 是否将 `X` 视为已核验的概率矩阵。 |
| `ci_level` | 置信水平。 |
| `bootstrap_ci`, `bootstrap_B`, `bootstrap_seed`, `verbose.boot` | 分层 bootstrap 控制项。 |

输出包括 `labeled_total`、`unlabeled_union_total`、
`unlabeled_pair_total` 及各自 SE；`layer_table` 给出各重叠层的模式数量、
抽样分配和估计；`clumping` 给出按重叠节点数和共享边数划分的结构指标；
`inference` 给出区间与抽样权重诊断。

```r
est.pairunion(P, H, B = 300, mode = "expectation",
              directed = FALSE, seed = 3,
              verbose = FALSE, verbose.B = FALSE)
```

## 4. `est.variance()`：VISS 子图数量方差

| 参数 | 含义与要求 |
|---|---|
| `P` | 独立边概率矩阵；若 `normalize_weight=TRUE`，也可为非负权重矩阵。 |
| `H` | 目标模体。 |
| `B` | FPUG/IPUM pairwise-union 总抽样次数。 |
| `B_mean` | 未提供均值时用于均值估计的 ISS 次数，默认等于 `B`。 |
| `directed` | 是否有向。 |
| `induced` | 是否进行诱导子图估计。 |
| `normalize_weight` | 是否将原始权重除以最大值。 |
| `overlap_sizes` | 重叠层；默认 `2:k`。 |
| `allocation` | `"equal"` 或 `"proportional"`。 |
| `mu_labeled`, `mu_unlabeled` | 可选已知均值；提供后可避免相应均值估计。 |
| `aut_count` | 可选已知自同构数量。 |
| `seed` | 随机种子。 |
| `verbose`, `verbose.B` | 重叠层与层内进度开关。 |
| `include_incompatible` | 是否纳入不相容模式的 IPUM；完整方差估计通常应保持 `TRUE`。 |
| `preprocessed` | 是否把 `P` 视为已核验的概率矩阵。 |
| `ci_level` | 置信水平。 |
| `bootstrap_ci`, `bootstrap_B`, `bootstrap_seed`, `verbose.boot` | 分层 bootstrap 控制项。 |

输出结构：

- `mean`：labeled/unlabeled 均值以及是否在本次调用中估计；
- `variance`：labeled/unlabeled 原始方差、截断至非负后的方差及 Monte
  Carlo SE；
- `layer_table`：每个重叠层的 FPUG joint mass、IPUM marginal-product
  mass、二者之差及 SE；
- `overdispersion`：原始和非负过度离散系数；
- `clumping`：整体及分层重叠结构；
- `inference`：Wald/log-Wald/bootstrap 区间及 FPUG、IPUM 权重诊断。

```r
v <- est.variance(P, H, B = 500, B_mean = 500,
                  directed = FALSE, seed = 4,
                  verbose = FALSE, verbose.B = FALSE)
v$variance$unlabeled
```

## 5. `sec.moment()`：二阶矩

`sec.moment(...)` 接受 `est.variance()` 的全部参数。输出保留完整方差结果，
并增加 `second_moment`，其中含 labeled/unlabeled 的 `E[N^2]` 和相应 SE。

```r
sec.moment(P = P, H = H, B = 300, B_mean = 300,
           directed = FALSE, seed = 5,
           verbose = FALSE, verbose.B = FALSE)$second_moment
```

## 6. `tsge()`：模体检验

通用参数为 `G`、`H`、`network_type` 和 `normalize_weight`。通过 `...`
传递：

| 参数 | 含义与选项 |
|---|---|
| `B_observed` | 观测网络 ISS 次数。 |
| `R_null` | 零模型期望 ISS 次数。 |
| `directed` | 是否有向。 |
| `null_model` | binary：`chung_lu`、`erdos_renyi`、`given`；weighted：`weighted_chung_lu`、`weighted_erdos_renyi`、`given`。 |
| `P_null` | `null_model="given"` 时必须提供的概率矩阵。 |
| `null_variance` | `poisson`、`sample`、`direct`（`sample` 别名）、`none`、`provided`；weighted 还可用 `calibrated`。 |
| `tau2_null` | `null_variance="provided"` 时的非负方差。 |
| `K_null_variance` | weighted calibrated 分支模拟网络数量，至少 2。 |
| `B_null_variance` | sample/calibrated 分支中每次使用的抽样次数。 |
| `null_variance_overlap_sizes` | 传给 `est.variance()` 的重叠层。 |
| `null_variance_allocation` | `equal` 或 `proportional`。 |
| `alternative` | `greater`、`less` 或 `two.sided`。 |
| `seed`, `order`, `automorphism_count_value`, `verbose` | 随机种子、贯序、自同构数量和进度控制。 |
| `null_expectation`, `observed` | binary 和 weighted 分支均可复用的既有 `est.subgraph` 结果。 |

输出 `coefficient_table`（log enrichment）、`auxiliary_table`（数量差异）、
`estimates`、底层 `observed`/`null_expectation` 估计，以及选择 sample 或
calibrated 方差时的完整方差结果或随机网络计数。

```r
tsge(P, H, network_type = "weighted", normalize_weight = FALSE,
     directed = FALSE, null_model = "given", P_null = P,
     null_variance = "none", B_observed = 300, R_null = 300,
     seed = 6, verbose = FALSE)
```

## 7. `tsge.cal()`：随机网络校准

| 参数 | 含义与选项 |
|---|---|
| `enrichment_result` | 可复用的 `tsge()` 结果。 |
| `observed` | 可复用的观测网络 `est.subgraph()` 结果；binary 与 weighted 均支持。 |
| `null_expectation` | 可复用的零模型 `est.subgraph()` 结果；binary 与 weighted 均支持。 |
| `K` | 随机网络数量；最小经验 p 值为 `1/(K+1)`。 |
| `R_null` | 每个随机网络的 ISS 次数。 |
| `B_observed` | 未提供既有结果时的观测网络 ISS 次数。 |
| `directed` | 是否有向。 |
| `null_model` | binary 还支持 `degree_swap`；weighted 使用加权 Chung–Lu、加权 ER 或 `given`。 |
| `P_null` | `given` 模型的概率矩阵。 |
| `statistic` | `log_enrichment`、`count` 或 `excess`。 |
| `alternative` | `greater`、`less` 或 `two.sided`。 |
| `seed`, `order`, `automorphism_count_value` | 随机种子、贯序、自同构数量。 |
| `null_expectation_value` | 兼容原接口的可选数值型零假设期望；若同时提供 `null_expectation`，以后者为准。 |
| `degree_swap_steps` | binary degree-swap 分支的成功换边次数。 |
| `log_offset` | log enrichment 的非负偏移，可能出现零计数时建议设为正数。 |
| `verbose`, `verbose.K` | 单次 ISS 和随机网络循环的进度开关。 |

显式传入的 `observed` 和 `null_expectation` 优先于 `enrichment_result`
中的对应对象或数值。输出 `calibration_table`、有限模拟经验 p 值 `empirical.p`、
`random_counts`、`random_count_ses`、`random_statistics`、观测数量及 SE、
零假设期望、复用后的 `observed` / `null_expectation_result` 对象、模体贯序
与自同构数量等。

```r
cal <- tsge.cal(P, H, network_type = "weighted",
                normalize_weight = FALSE, directed = FALSE,
                null_model = "given", P_null = P,
                K = 9, R_null = 200, B_observed = 200,
                statistic = "count", seed = 7,
                verbose = FALSE, verbose.K = FALSE)
cal$calibration_table
cal$empirical.p
```

已有估计可在两类网络中直接复用，例如：

```r
obs <- est.subgraph(P, H, network_type = "weighted",
                    normalize_weight = FALSE, directed = FALSE,
                    B = 200, seed = 8, verbose = FALSE)
tsge.cal(P, H, network_type = "weighted", normalize_weight = FALSE,
         directed = FALSE, null_model = "given", P_null = P,
         observed = obs, null_expectation = obs,
         K = 9, R_null = 200, statistic = "count",
         seed = 9, verbose = FALSE, verbose.K = FALSE)
```

## 8. `PairUnion()` 与 `PairVar()`：精确枚举

这两个函数不进行 ISS 抽样，结果没有 Monte Carlo 误差，适合小型网络、
理论值计算以及对 `est.pairunion()`/`est.variance()` 的核验。一般枚举的
计算量会随网络、模体及 union 规模迅速增长；齐次概率网络会自动采用快速
闭式计算。

`PairUnion()` 参数：

| 参数 | 含义与要求 |
|---|---|
| `X` | `mode="count"` 时为 0/1 方阵；`mode="expectation"` 时为 `[0,1]` 概率方阵。 |
| `H` | 目标模体方阵。 |
| `directed` | 是否为有向网络。 |
| `induced` | 是否同时施加非边约束。 |
| `overlap_sizes` | 两份模体共享的节点数，默认 `2:nrow(H)`。 |
| `mode` | `count` 或 `expectation`。 |
| `include_incompatible` | 是否在 `by_pattern` 中显示不相容模式；这些模式对 joint total 的贡献始终为 0。 |
| `method` | `auto`、`enumerate` 或 `homogeneous`。 |
| `homogeneous_tol` | `auto` 判断齐次概率时使用的数值容差。 |

输出 `labeled_total`、`unlabeled_union_total`、
`unlabeled_pair_total`、模体自同构数量、逐模式结果 `by_pattern` 以及实际
采用的计算方法。

`PairVar()` 的 `P` 必须是 `[0,1]` 概率方阵；其余 `H`、`directed`、
`induced`、`overlap_sizes`、`method` 和 `homogeneous_tol` 与上述含义相同。
`include_incompatible` 仅决定是否在 `by_pattern` 中显示不相容模式；其
IPUM marginal-product mass 始终纳入总方差。

`PairVar()` 输出 `labeled_variance`、`unlabeled_variance`、
`labeled_joint_total`、`labeled_marginal_product_total`、逐模式 joint/IPUM/
covariance 分解和自同构数量。

```r
exact_union <- PairUnion(P, H, mode = "expectation", directed = FALSE)
exact_var <- PairVar(P, H, directed = FALSE)

exact_union$by_pattern
exact_var$unlabeled_variance
```
