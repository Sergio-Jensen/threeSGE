library(threeSGE)

H3 <- matrix(1, 3, 3)
diag(H3) <- 0

G5 <- matrix(1, 5, 5)
diag(G5) <- 0

P6 <- matrix(0.5, 6, 6)
diag(P6) <- 0

binary_fit <- est.subgraph(
  G5, H3, B = 100, directed = FALSE, seed = 1, verbose = FALSE
)
stopifnot(binary_fit$unlabeled_estimate == choose(5, 3))

weighted_fit <- est.subgraph(
  P6, H3, B = 100, directed = FALSE, normalize_weight = FALSE,
  seed = 1, verbose = FALSE
)
stopifnot(abs(weighted_fit$unlabeled_expectation - 2.5) < 1e-10)

variance_fit <- est.variance(
  P6, H3, B = 100, B_mean = 100, directed = FALSE,
  seed = 2, verbose = FALSE, verbose.B = FALSE
)
stopifnot(abs(variance_fit$variance$unlabeled - 5) < 1e-10)

union_fit <- est.pairunion(
  P6, H3, B = 100, mode = "expectation", directed = FALSE,
  seed = 3, verbose = FALSE, verbose.B = FALSE
)
stopifnot(inherits(union_fit, "est.pairunion"))

tsge_fit <- tsge(
  P6, H3, network_type = "weighted", normalize_weight = FALSE,
  directed = FALSE, null_model = "given", P_null = P6,
  null_variance = "none", B_observed = 100, R_null = 100,
  seed = 4, verbose = FALSE
)
stopifnot(inherits(tsge_fit, "tsge"))

P5 <- P6[seq_len(5), seq_len(5)]
binary_null <- est.subgraph(
  P5, H3, network_type = "weighted", normalize_weight = FALSE,
  B = 100, directed = FALSE, seed = 40, verbose = FALSE
)
reused_binary_fit <- tsge(
  G5, H3, observed = binary_fit, null_expectation = binary_null,
  null_model = "given", P_null = P5,
  null_variance = "none", directed = FALSE, verbose = FALSE
)
stopifnot(
  identical(reused_binary_fit$observed, binary_fit),
  identical(reused_binary_fit$null_expectation, binary_null)
)

cal_fit <- tsge.cal(
  P6, H3, network_type = "weighted", normalize_weight = FALSE,
  directed = FALSE, null_model = "given", P_null = P6,
  K = 3, B_observed = 50, R_null = 50, statistic = "count",
  seed = 5, verbose = FALSE, verbose.K = FALSE
)
stopifnot(length(cal_fit$random_counts) == 3)

cal_binary_reuse <- tsge.cal(
  G5, H3, observed = binary_fit, null_expectation = binary_null,
  directed = FALSE, null_model = "given", P_null = P5,
  K = 3, R_null = 50, statistic = "count",
  seed = 6, verbose = FALSE, verbose.K = FALSE
)
stopifnot(
  identical(cal_binary_reuse$observed, binary_fit),
  identical(cal_binary_reuse$null_expectation_result, binary_null),
  cal_binary_reuse$observed_count == binary_fit$unlabeled_estimate,
  cal_binary_reuse$null_expectation == binary_null$unlabeled_expectation
)

cal_weighted_reuse <- tsge.cal(
  P6, H3, network_type = "weighted", normalize_weight = FALSE,
  observed = weighted_fit, null_expectation = weighted_fit,
  directed = FALSE, null_model = "given", P_null = P6,
  K = 3, R_null = 50, statistic = "count",
  seed = 7, verbose = FALSE, verbose.K = FALSE
)
stopifnot(
  identical(cal_weighted_reuse$observed, weighted_fit),
  identical(cal_weighted_reuse$null_expectation_result, weighted_fit),
  cal_weighted_reuse$observed_count == weighted_fit$unlabeled_expectation,
  cal_weighted_reuse$null_expectation == weighted_fit$unlabeled_expectation
)

exact_union <- PairUnion(
  P6, H3, mode = "expectation", directed = FALSE
)
exact_variance <- PairVar(P6, H3, directed = FALSE)
stopifnot(
  abs(exact_union$unlabeled_pair_total - 8.125) < 1e-10,
  abs(exact_variance$unlabeled_variance - 5) < 1e-10
)
