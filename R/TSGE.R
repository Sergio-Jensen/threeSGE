#==========================================================================
## Three-Stage Graph Estimation (TSGE) motif verification
#==========================================================================
#--------------------------------------------------------------------------
## Unweighted Implementation
#--------------------------------------------------------------------------

is.compact.probmat <- function(P){
  is.list(P) && !is.null(P$compact_probmat) && isTRUE(P$compact_probmat)
}
compact.to.matrix <- function(P, max_entries = 5e7){
  if(!is.compact.probmat(P)){return(P)}
  if(P$n * P$n > max_entries){
    stop("Sampled null variance currently needs an explicit probability matrix; compact null matrix is too large to expand safely.")
  }
  idx <- expand.grid(from = seq_len(P$n), to = seq_len(P$n))
  M <- matrix(prob.values(P, idx$from, idx$to, eps = 0), P$n, P$n)
  diag(M) <- 0
  if(!P$directed){
    M <- (M + t(M))/2
    diag(M) <- 0
  }
  M
}
make.compact.probmat <- function(type, n, directed = FALSE, density = NA_real_,
                                 p = NULL, degree = NULL, out_degree = NULL,
                                 in_degree = NULL, m = NULL, cap = 1 - 1e-12){
  structure(
    list(
      compact_probmat = TRUE,
      type = type,
      n = as.integer(n),
      directed = directed,
      density = density,
      p = p,
      degree = degree,
      out_degree = out_degree,
      in_degree = in_degree,
      m = m,
      cap = cap
    ),
    class = c("compact.probmat", "list")
  )
}
prob.values <- function(P, from, to, eps = 1e-300){
  if(!is.compact.probmat(P)){
    return(pmin(pmax(P[cbind(from, to)], eps), 1 - 1e-12))
  }
  if(P$type == "erdos_renyi"){
    p <- rep(P$p, length(from))
  }else if(P$directed){
    p <- P$out_degree[from] * P$in_degree[to]/P$m
  }else{
    p <- P$degree[from] * P$degree[to]/(2 * P$m)
  }
  pmin(pmax(p, eps), P$cap)
}
norm.probmat <- function(P, directed = FALSE, cap = 1 - 1e-12){
  if(is.compact.probmat(P)){return(P)}
  motif.norm.probmat(P, directed = directed, cap = cap)
}
make.null <- function(G, method = c("chung_lu", "erdos_renyi"), directed = FALSE,
                      cap = 1 - 1e-12, compact_threshold = 5e7){
  method <- match.arg(method)
  G <- motif.graph.to.matrix(G, weighted = FALSE)
  G <- motif.norm.graph(G, directed = directed)
  n <- nrow(G)
  use_compact <- (n * n) > compact_threshold
  if(method == "erdos_renyi"){
    p <- motif.graph.density(G, directed = directed)
    if(use_compact){
      return(make.compact.probmat(
        type = "erdos_renyi",
        n = n,
        directed = directed,
        density = p,
        p = min(max(p, 0), cap),
        cap = cap
      ))
    }
    P <- matrix(p, n, n)
    diag(P) <- 0
    return(norm.probmat(P, directed = directed, cap = cap))
  }
  if(directed){
    m <- sum(G)
    if(m == 0){
      if(use_compact){
        return(make.compact.probmat(
          type = "erdos_renyi",
          n = n,
          directed = TRUE,
          density = 0,
          p = 0,
          cap = cap
        ))
      }
      return(matrix(0, n, n))
    }
    out_degree <- rowSums(G)
    in_degree <- colSums(G)
    if(use_compact){
      return(make.compact.probmat(
        type = "chung_lu",
        n = n,
        directed = TRUE,
        density = m/(n * (n - 1)),
        out_degree = out_degree,
        in_degree = in_degree,
        m = m,
        cap = cap
      ))
    }
    P <- outer(out_degree, in_degree, "*")/m
    diag(P) <- 0
    return(norm.probmat(P, directed = TRUE, cap = cap))
  }
  m <- sum(G)/2
  if(m == 0){
    if(use_compact){
      return(make.compact.probmat(
        type = "erdos_renyi",
        n = n,
        directed = FALSE,
        density = 0,
        p = 0,
        cap = cap
      ))
    }
    return(matrix(0, n, n))
  }
  degree <- rowSums(G)
  if(use_compact){
    return(make.compact.probmat(
      type = "chung_lu",
      n = n,
      directed = FALSE,
      density = m/choose(n, 2),
      degree = degree,
      m = m,
      cap = cap
    ))
  }
  P <- outer(degree, degree, "*")/(2 * m)
  diag(P) <- 0
  norm.probmat(P, directed = FALSE, cap = cap)
}
.tsge.result <- function(call, theta, theta_se, mu0, mu0_se,
                         alternative, null_model, null_variance,
                         motif_order, automorphism_count,
                         observed, null_expectation,
                         sample_null_variance = NULL,
                         calibrated_null_counts = NULL,
                         calibrated_null_ses = NULL,
                         calibrated_tau2 = NULL,
                         tau2_null = NULL,
                         probability_rule = NULL,
                         result_class = "tsge"){
  tau2 <- if(null_variance == "calibrated"){
    calibrated_tau2
  }else if(null_variance == "sample"){
    sample_null_variance$variance$unlabeled
  }else if(null_variance == "poisson"){
    mu0
  }else if(null_variance == "provided"){
    if(is.null(tau2_null) || !is.numeric(tau2_null) || length(tau2_null) != 1 ||
       is.na(tau2_null) || tau2_null < 0){
      stop("tau2_null must be a non-negative number when null_variance = 'provided'.")
    }
    tau2_null
  }else{
    0
  }
  if(theta <= 0 || mu0 <= 0){
    warning("theta or mu0 is non-positive; log enrichment is undefined. Returning NA coefficient.")
    beta <- beta_se_mc <- beta_se_total <- z <- p <- NA_real_
  }else{
    beta <- log(theta/mu0)
    beta_se_mc <- sqrt((theta_se/theta)^2 + (mu0_se/mu0)^2)
    beta_se_total <- sqrt(beta_se_mc^2 + tau2/(mu0^2))
    z <- beta/beta_se_total
    p <- motif.wald.p(z, alternative = alternative)
  }
  excess <- theta - mu0
  excess_se_mc <- sqrt(theta_se^2 + mu0_se^2)
  excess_se_total <- sqrt(excess_se_mc^2 + tau2)
  excess_z <- excess/excess_se_total
  excess_p <- motif.wald.p(excess_z, alternative = alternative)
  result <- list(
    call = call,
    coefficient_table = data.frame(
      term = "log_enrichment", estimate = beta, std.error = beta_se_total,
      statistic = z, p.value = p, row.names = NULL
    ),
    auxiliary_table = data.frame(
      term = "excess_count", estimate = excess, std.error = excess_se_total,
      statistic = excess_z, p.value = excess_p, row.names = NULL
    ),
    alternative = alternative,
    null_model = null_model,
    null_variance = null_variance,
    probability_rule = probability_rule,
    motif_order = motif_order,
    automorphism_count = automorphism_count,
    observed = observed,
    null_expectation = null_expectation,
    sample_null_variance = sample_null_variance,
    direct_null_variance = sample_null_variance,
    calibrated_null_counts = calibrated_null_counts,
    calibrated_null_ses = calibrated_null_ses,
    calibrated_tau2 = calibrated_tau2,
    estimates = list(
      observed_count = theta,
      observed_se = theta_se,
      null_expectation = mu0,
      null_expectation_se = mu0_se,
      enrichment_ratio = if(mu0 > 0) theta/mu0 else NA_real_,
      log_enrichment = beta,
      log_enrichment_se_mc = beta_se_mc,
      log_enrichment_se_total = beta_se_total,
      excess = excess,
      excess_se_mc = excess_se_mc,
      excess_se_total = excess_se_total,
      excess_z = excess_z,
      excess_p_value = excess_p
    )
  )
  class(result) <- result_class
  result
}

.tsge.binary <- function(G, H, B_observed = 10000, R_null = 10000, directed = FALSE,
                       null_model = c("chung_lu", "erdos_renyi", "given"),
                       P_null = NULL, null_expectation = NULL,
                       null_variance = c("poisson", "sample", "direct", "none", "provided"), tau2_null = NULL,
                       B_null_variance = R_null,
                       null_variance_overlap_sizes = NULL,
                       null_variance_allocation = c("equal", "proportional"),
                       alternative = c("greater", "less", "two.sided"), seed = NULL,
                       order = NULL, automorphism_count_value = NULL,
                       observed = NULL, verbose = TRUE){
  null_model <- match.arg(null_model)
  null_variance <- match.arg(null_variance)
  if(null_variance == "direct"){null_variance <- "sample"}
  null_variance_allocation <- match.arg(null_variance_allocation)
  alternative <- match.arg(alternative)
  if(!is.null(seed)){set.seed(seed)}
  G <- motif.graph.to.matrix(G, weighted = FALSE)
  H <- motif.graph.to.matrix(H, weighted = FALSE)
  G <- motif.norm.graph(G, directed = directed)
  H <- motif.norm.graph(H, directed = directed)
  if(is.null(order)){
    density <- motif.graph.density(G, directed = directed)
    order <- motif.order(H, density, directed = directed)
  }
  if(is.null(automorphism_count_value)){
    automorphism_count_value <- motif.automorphism.count(H, directed = directed)
  }
  if(is.null(observed)){
    observed_args <- list(
      G = G, H = H, B = B_observed, directed = directed,
      seed = if(is.null(seed)) NULL else seed + 1L,
      order = order, keep_mappings = FALSE, verbose = verbose
    )
    observed <- do.call(est.subgraph, observed_args)
  }
  P <- if(null_model == "given"){
    if(is.null(P_null)){stop("P_null must be supplied when null_model = 'given'.")}
    norm.probmat(P_null, directed = directed)
  }else{
    make.null(G, method = null_model, directed = directed)
  }
  P_explicit <- if(is.compact.probmat(P)){compact.to.matrix(P)}else{P}
  if(is.null(null_expectation)){
    null_expectation <- est.subgraph(
      G = P_explicit,
      H = H,
      B = R_null,
      directed = directed,
      aut_count = automorphism_count_value,
      normalize_weight = FALSE,
      seed = if(is.null(seed)) NULL else seed + 2L,
      verbose = verbose,
      network_type = "weighted"
    )
  }
  sample_null_variance <- NULL
  theta <- motif.est.value(observed)
  theta_se <- motif.est.se(observed)
  mu0 <- motif.est.value(null_expectation)
  mu0_se <- motif.est.se(null_expectation)
  if(null_variance == "sample"){
    if(!exists("est.variance", mode = "function")){
      stop("Source VISS.R before requesting sampled null variance.")
    }
    sample_null_variance <- est.variance(
      P = if(is.compact.probmat(P)){compact.to.matrix(P)}else{P},
      H = H,
      B = B_null_variance,
      directed = directed,
      induced = TRUE,
      normalize_weight = FALSE,
      overlap_sizes = null_variance_overlap_sizes,
      allocation = null_variance_allocation,
      include_incompatible = TRUE,
      mu_unlabeled = mu0,
      aut_count = automorphism_count_value,
      seed = if(is.null(seed)) NULL else seed + 3L,
      verbose = verbose
    )
  }
  .tsge.result(
    call = match.call(),
    theta = theta,
    theta_se = theta_se,
    mu0 = mu0,
    mu0_se = mu0_se,
    alternative = alternative,
    null_model = null_model,
    null_variance = null_variance,
    tau2_null = tau2_null,
    motif_order = order,
    automorphism_count = automorphism_count_value,
    observed = observed,
    null_expectation = null_expectation,
    sample_null_variance = sample_null_variance,
    result_class = "tsge"
  )
}
print.tsge <- function(x, ...){
  weighted <- identical(x$network_type, "weighted")
  cat(if(weighted){"Weighted TSGE motif test\n"}else{"TSGE motif test\n"})
  if(weighted){cat("  probability rule: ", x$probability_rule, "\n", sep = "")}
  cat("  null model:       ", x$null_model, "\n", sep = "")
  cat("  null variance:    ", x$null_variance, "\n", sep = "")
  cat("  alternative:      ", x$alternative, "\n", sep = "")
  cat("  motif order:      ", paste(x$motif_order, collapse = " -> "), "\n", sep = "")
  cat("  automorphisms:    ", x$automorphism_count, "\n", sep = "")
  cat("\nCounts\n")
  cat(if(weighted){"  observed expected count: "}else{"  observed count:   "},
      signif(x$estimates$observed_count, 6),
      " (SE ", signif(x$estimates$observed_se, 4), ")\n", sep = "")
  cat("  null expectation: ", signif(x$estimates$null_expectation, 6),
      " (SE ", signif(x$estimates$null_expectation_se, 4), ")\n", sep = "")
  if(!is.null(x$calibrated_tau2)){
    cat("  calibrated null SD: ", signif(sqrt(x$calibrated_tau2), 6), "\n", sep = "")
  }
  cat("  beta = log(theta/mu0): ", signif(x$estimates$log_enrichment, 6),
      " (total SE ", signif(x$estimates$log_enrichment_se_total, 4), ")\n", sep = "")
  cat("  exp(beta):        ", signif(x$estimates$enrichment_ratio, 6), "\n", sep = "")
  cat("\nCoefficient table\n")
  print(x$coefficient_table, row.names = FALSE)
  ## cat("\nAuxiliary count-difference table\n")
  ## print(x$auxiliary_table, row.names = FALSE)
  invisible(x)
}
summary.tsge <- function(object, ...){
  list(coefficients = object$coefficient_table,
       auxiliary = object$auxiliary_table)
}
randomize <- function(G, n_swap = NULL, max_tries = NULL,
                      directed = FALSE){
  G <- motif.graph.to.matrix(G, weighted = FALSE)
  G <- motif.norm.graph(G, directed = directed)
  n <- nrow(G)
  if(directed){
    edges <- which(G, arr.ind = TRUE)
    m <- nrow(edges)
    if(m < 2){return(G)}
    if(is.null(n_swap)){n_swap <- 10L * m}
    if(is.null(max_tries)){max_tries <- 50L * n_swap}
    swaps <- 0L
    tries <- 0L
    while (swaps < n_swap && tries < max_tries){
      tries <- tries + 1L
      idx <- sample.int(m, 2)
      a <- edges[idx[1], 1]
      b <- edges[idx[1], 2]
      c <- edges[idx[2], 1]
      d <- edges[idx[2], 2]
      if(a == c || b == d || a == d || c == b){next}
      if(G[a, d] || G[c, b]){next}
      G[a, b] <- FALSE
      G[c, d] <- FALSE
      G[a, d] <- TRUE
      G[c, b] <- TRUE
      edges[idx[1], ] <- c(a, d)
      edges[idx[2], ] <- c(c, b)
      swaps <- swaps + 1L
    }
    return(G)
  }
  edges <- which(upper.tri(G) & G, arr.ind = TRUE)
  m <- nrow(edges)
  if(m < 2){return(G)}
  if(is.null(n_swap)){n_swap <- 10L * m}
  if(is.null(max_tries)){max_tries <- 50L * n_swap}
  sorted_edge <- function(x, y){
    if(x < y) c(x, y) else c(y, x)
  }
  swaps <- 0L
  tries <- 0L
  while (swaps < n_swap && tries < max_tries){
    tries <- tries + 1L
    idx <- sample.int(m, 2)
    a <- edges[idx[1], 1]
    b <- edges[idx[1], 2]
    c <- edges[idx[2], 1]
    d <- edges[idx[2], 2]
    if(length(unique(c(a, b, c, d))) < 4){next}
    if(stats::runif(1) < 0.5){
      e1 <- sorted_edge(a, d)
      e2 <- sorted_edge(c, b)
    }else{
      e1 <- sorted_edge(a, c)
      e2 <- sorted_edge(b, d)
    }
    if(identical(e1, e2)){next}
    if(G[e1[1], e1[2]] || G[e2[1], e2[2]]){next}
    G[a, b] <- G[b, a] <- FALSE
    G[c, d] <- G[d, c] <- FALSE
    G[e1[1], e1[2]] <- G[e1[2], e1[1]] <- TRUE
    G[e2[1], e2[2]] <- G[e2[2], e2[1]] <- TRUE
    edges[idx[1], ] <- e1
    edges[idx[2], ] <- e2
    swaps <- swaps + 1L
  }
  return(G)
}
cal.stat <- function(count, mu0, statistic = c("log_enrichment", "count", "excess"), log_offset = 0){
  statistic <- match.arg(statistic)
  if(statistic == "count"){return(count)}
  if(statistic == "excess"){return(count - mu0)}
  log((count + log_offset)/(mu0 + log_offset))
}
emp.p <- function(observed, null_values, center, alternative = c("greater", "less", "two.sided")){
  alternative <- match.arg(alternative)
  if(alternative == "greater"){
    return((1 + sum(null_values >= observed, na.rm = TRUE))/(length(null_values) + 1))
  }
  if(alternative == "less"){
    return((1 + sum(null_values <= observed, na.rm = TRUE))/(length(null_values) + 1))
  }
  (1 + sum(abs(null_values - center) >= abs(observed - center), na.rm = TRUE)) /
    (length(null_values) + 1)
}

.tsge.cal.result <- function(call, observed_count, observed_count_se,
                             null_expectation, random_counts, random_ses,
                             statistic, alternative, null_model,
                             K, R_null, motif_order, automorphism_count,
                             log_offset = 0, probability_rule = NULL,
                             observed_result = NULL,
                             null_expectation_result = NULL,
                             result_class = "tsge.cal"){
  observed_stat <- cal.stat(
    observed_count, null_expectation,
    statistic = statistic, log_offset = log_offset
  )
  random_stats <- vapply(
    random_counts, cal.stat, numeric(1), mu0 = null_expectation,
    statistic = statistic, log_offset = log_offset
  )
  center <- switch(statistic,
                   count = null_expectation,
                   log_enrichment = 0,
                   excess = 0)
  finite_stats <- random_stats[is.finite(random_stats)]
  if(length(finite_stats) < length(random_stats)){
    warning("Some calibrated statistics are non-finite. Consider setting log_offset > 0.")
  }
  null_mean <- if(length(finite_stats) > 0) mean(finite_stats) else NA_real_
  null_sd <- if(length(finite_stats) > 1) stats::sd(finite_stats) else NA_real_
  calibrated_z <- if(is.finite(observed_stat) && is.finite(null_sd) && null_sd > 0){
    (observed_stat - null_mean)/null_sd
  }else{
    NA_real_
  }
  result <- list(
    call = call,
    calibration_table = data.frame(
      term = paste0("calibrated_", statistic),
      estimate = observed_stat,
      null.mean = null_mean,
      null.sd = null_sd,
      statistic = calibrated_z,
      p.value = motif.wald.p(calibrated_z, alternative = alternative),
      row.names = NULL
    ),
    statistic = statistic,
    alternative = alternative,
    null_model = null_model,
    probability_rule = probability_rule,
    K = K,
    R_null = R_null,
    observed_count = observed_count,
    observed_count_se = observed_count_se,
    null_expectation = null_expectation,
    observed = observed_result,
    null_expectation_result = null_expectation_result,
    random_counts = random_counts,
    random_count_ses = random_ses,
    observed_statistic = observed_stat,
    random_statistics = random_stats,
    motif_order = motif_order,
    automorphism_count = automorphism_count,
    empirical.p = emp.p(
      observed = observed_stat,
      null_values = random_stats,
      center = center,
      alternative = alternative
    ),
    min_p_value = 1/(K + 1)
  )
  class(result) <- result_class
  result
}

.tsge.binary.cal <- function(G, H, enrichment_result = NULL,
                           observed = NULL, null_expectation = NULL,
                           K = 99, R_null = 10000, B_observed = 10000,
                           directed = FALSE, null_model = NULL, P_null = NULL, statistic = c("log_enrichment", "count", "excess"),
                           alternative = c("greater", "less", "two.sided"), seed = NULL, order = NULL,
                           automorphism_count_value = NULL, null_expectation_value = NULL,
                           degree_swap_steps = NULL, log_offset = 0, verbose = TRUE, verbose.K = TRUE){
  statistic <- match.arg(statistic)
  alternative <- match.arg(alternative)
  if(!is.numeric(K) || length(K) != 1 || is.na(K) || K < 1 || K != as.integer(K)){
    stop("K must be a positive integer.")
  }
  if(!is.numeric(R_null) || length(R_null) != 1 || is.na(R_null) ||
     R_null < 1 || R_null != as.integer(R_null)){
    stop("R_null must be a positive integer.")
  }
  K <- as.integer(K)
  R_null <- as.integer(R_null)
  if(!is.null(seed)){set.seed(seed)}
  G <- motif.graph.to.matrix(G, weighted = FALSE)
  H <- motif.graph.to.matrix(H, weighted = FALSE)
  G <- motif.norm.graph(G, directed = directed)
  H <- motif.norm.graph(H, directed = directed)
  observed_result <- observed
  null_expectation_result <- null_expectation
  observed_count <- NULL
  observed_count_se <- NULL
  if(!is.null(enrichment_result)){
    if(is.null(null_model)) null_model <- enrichment_result$null_model
    if(is.null(order)) order <- enrichment_result$motif_order
    if(is.null(automorphism_count_value)){
      automorphism_count_value <- enrichment_result$automorphism_count
    }
    if(is.null(observed_result)){
      observed_result <- enrichment_result$observed
    }
    if(is.null(null_expectation_result)){
      null_expectation_result <- enrichment_result$null_expectation
    }
    if(is.null(null_expectation_value)){
      null_expectation_value <- enrichment_result$estimates$null_expectation
    }
    if(is.null(observed_result)){
      observed_count <- enrichment_result$estimates$observed_count
      observed_count_se <- enrichment_result$estimates$observed_se
    }
  }
  if(!is.null(observed_result)){
    observed_count <- motif.est.value(observed_result)
    observed_count_se <- motif.est.se(observed_result)
    if(is.null(order)){
      order <- motif.est.order(observed_result)
    }
  }
  if(!is.null(null_expectation_result)){
    null_expectation_value <- motif.est.value(null_expectation_result)
  }
  if(is.null(null_model)) null_model <- "chung_lu"
  null_model <- match.arg(null_model, c("chung_lu", "erdos_renyi", "given", "degree_swap"))
  if(is.null(order)){
    order <- motif.order(
      H, motif.graph.density(G, directed = directed), directed = directed
    )
  }
  if(is.null(automorphism_count_value)){
    automorphism_count_value <- motif.automorphism.count(H, directed = directed)
  }
  if(is.null(observed_count)){
    observed_args <- list(G = G, H = H, B = B_observed,
                          directed = directed, seed = if(is.null(seed)) NULL else seed + 1L,
                          order = order, keep_mappings = FALSE, verbose = verbose)
    observed_result <- do.call(est.subgraph, observed_args)
    observed_count <- motif.est.value(observed_result)
    observed_count_se <- motif.est.se(observed_result)
  }
  P <- NULL
  if(null_model == "given"){
    if(is.null(P_null)) stop("P_null must be supplied when null_model = 'given'.")
    P <- norm.probmat(P_null, directed = directed)
  }else if(null_model %in% c("chung_lu", "erdos_renyi")){
    P <- make.null(G, method = null_model, directed = directed)
  }
  random_counts <- numeric(K)
  random_ses <- numeric(K)
  for(i in seq_len(K)){
    G_star <- if(null_model == "degree_swap"){
      randomize(G, n_swap = degree_swap_steps, directed = directed)
    }else{
      motif.sample.graph(P, directed = directed)
    }
    random_args <- list(G = G_star,
                        H = H,
                        B = R_null,
                        directed = directed,
                        seed = NULL,
                        order = order,
                        keep_mappings = FALSE,
                        verbose = verbose)
    random_est <- do.call(est.subgraph, random_args)
    random_counts[i] <- random_est$unlabeled_estimate
    random_ses[i] <- random_est$unlabeled_se
    if(verbose.K && K >= 10 && i %% max(1, floor(K/10)) == 0){
      cat("Finished", i, "of", K, "random-network calibrations\n")
    }
  }
  if(is.null(null_expectation_value)){null_expectation_value <- mean(random_counts)}
  .tsge.cal.result(
    call = match.call(),
    observed_count = observed_count,
    observed_count_se = observed_count_se,
    null_expectation = null_expectation_value,
    random_counts = random_counts,
    random_ses = random_ses,
    statistic = statistic,
    alternative = alternative,
    null_model = null_model,
    K = K,
    R_null = R_null,
    motif_order = order,
    automorphism_count = automorphism_count_value,
    log_offset = log_offset,
    observed_result = observed_result,
    null_expectation_result = null_expectation_result,
    result_class = "tsge.cal"
  )
}
print.tsge.cal <- function(x, ...){
  weighted <- identical(x$network_type, "weighted")
  cat(if(weighted){"Weighted random-network motif calibration\n"}else{"Random-network motif calibration\n"})
  if(weighted){cat("  probability rule: ", x$probability_rule, "\n", sep = "")}
  cat("  null model:       ", x$null_model, "\n", sep = "")
  cat("  statistic:        ", x$statistic, "\n", sep = "")
  cat("  alternative:      ", x$alternative, "\n", sep = "")
  cat("  random networks:  ", x$K, "\n", sep = "")
  cat("  min p-value:      ", signif(x$min_p_value, 4), "\n", sep = "")
  cat("  observed count:   ", signif(x$observed_count, 6),
      " (SE ", signif(x$observed_count_se, 4), ")\n", sep = "")
  cat("  null expectation: ", signif(x$null_expectation, 6), "\n", sep = "")
  cat("\nCalibration table\n")
  print(x$calibration_table, row.names = FALSE)
  invisible(x)
}
summary.tsge.cal <- function(object, ...){
  object$calibration_table
}

#--------------------------------------------------------------------------
## Weighted Implementation
#--------------------------------------------------------------------------
make.null.weighted <- function(P,
                               method = c("weighted_chung_lu",
                                          "weighted_erdos_renyi",
                                          "chung_lu",
                                          "erdos_renyi"),
                               directed = TRUE){
  method <- match.arg(method)
  if(method == "chung_lu"){method <- "weighted_chung_lu"}
  if(method == "erdos_renyi"){method <- "weighted_erdos_renyi"}
  P <- motif.norm.probmat(P, directed = directed, strict = FALSE)
  n <- nrow(P)
  if(method == "weighted_erdos_renyi"){
    p <- motif.prob.density(P, directed = directed)
    P0 <- matrix(p, n, n)
    diag(P0) <- 0
    return(motif.norm.probmat(P0, directed = directed, strict = FALSE))
  }
  if(directed){
    m <- sum(P)
    if(m <= 0){return(matrix(0, n, n))}
    out_strength <- rowSums(P)
    in_strength <- colSums(P)
    P0 <- outer(out_strength, in_strength, "*")/m
    diag(P0) <- 0
    return(motif.norm.probmat(P0, directed = TRUE, strict = FALSE))
  }
  m <- sum(P)/2
  if(m <= 0){return(matrix(0, n, n))}
  strength <- rowSums(P)
  P0 <- outer(strength, strength, "*")/(2 * m)
  diag(P0) <- 0
  motif.norm.probmat(P0, directed = FALSE, strict = FALSE)
}

.tsge.weighted <- function(W, H,
                                B_observed = 10000,
                                R_null = 10000,
                                directed = TRUE,
                                normalize_weight = TRUE,
                                null_model = c("weighted_chung_lu",
                                               "weighted_erdos_renyi",
                                               "given",
                                               "chung_lu",
                                               "erdos_renyi"),
                                P_null = NULL,
                                null_expectation = NULL,
                                null_variance = c("calibrated", "sample", "direct", "poisson", "none", "provided"),
                                tau2_null = NULL,
                                K_null_variance = 99,
                                B_null_variance = R_null,
                                null_variance_overlap_sizes = NULL,
                                null_variance_allocation = c("equal", "proportional"),
                                alternative = c("greater", "less", "two.sided"),
                                seed = NULL,
                                order = NULL,
                                automorphism_count_value = NULL,
                                observed = NULL,
                                verbose = TRUE){
  null_model <- match.arg(null_model)
  null_variance <- match.arg(null_variance)
  if(null_variance == "direct"){null_variance <- "sample"}
  null_variance_allocation <- match.arg(null_variance_allocation)
  alternative <- match.arg(alternative)
  if(!is.null(seed)){set.seed(seed)}
  P_obs <- motif.matrix.to.probmat.max(
    W,
    directed = directed,
    normalize_weight = normalize_weight,
    preprocessed = FALSE
  )
  H <- motif.graph.to.matrix(H, weighted = FALSE)
  H <- motif.norm.motif(H, directed = directed)
  if(is.null(order)){
    order <- motif.order(
      H, motif.prob.density(P_obs, directed = directed), directed = directed
    )
  }
  if(is.null(automorphism_count_value)){
    automorphism_count_value <- motif.automorphism.count(H, directed = directed)
  }
  if(is.null(observed)){
    observed <- est.subgraph(
      G = P_obs,
      H = H,
      network_type = "weighted",
      B = B_observed,
      directed = directed,
      seed = if(is.null(seed)) NULL else seed + 1L,
      order = order,
      aut_count = automorphism_count_value,
      keep_mappings = FALSE,
      verbose = verbose,
      normalize_weight = FALSE
    )
  }
  P0 <- if(null_model == "given"){
    if(is.null(P_null)){stop("P_null must be supplied when null_model = 'given'.")}
    motif.norm.probmat(P_null, directed = directed, strict = FALSE)
  }else{
    make.null.weighted(P_obs, method = null_model, directed = directed)
  }
  if(is.null(null_expectation)){
    null_expectation <- est.subgraph(
      G = P0,
      H = H,
      network_type = "weighted",
      B = R_null,
      directed = directed,
      seed = if(is.null(seed)) NULL else seed + 2L,
      order = order,
      aut_count = automorphism_count_value,
      keep_mappings = FALSE,
      verbose = verbose,
      normalize_weight = FALSE
    )
  }
  calibrated_null_counts <- NULL
  calibrated_null_ses <- NULL
  calibrated_tau2 <- NULL
  sample_null_variance <- NULL
  if(null_variance == "calibrated"){
    if(!is.numeric(K_null_variance) || length(K_null_variance) != 1 ||
        is.na(K_null_variance) || K_null_variance < 2 ||
        K_null_variance != as.integer(K_null_variance)){
      stop("K_null_variance must be an integer greater than or equal to 2 when null_variance = 'calibrated'.")
    }
    if(!is.numeric(B_null_variance) || length(B_null_variance) != 1 ||
        is.na(B_null_variance) || B_null_variance < 1 ||
        B_null_variance != as.integer(B_null_variance)){
      stop("B_null_variance must be a positive integer when null_variance = 'calibrated'.")
    }
    K_null_variance <- as.integer(K_null_variance)
    B_null_variance <- as.integer(B_null_variance)
    calibrated_null_counts <- numeric(K_null_variance)
    calibrated_null_ses <- numeric(K_null_variance)
    for(i in seq_len(K_null_variance)){
      G_star <- motif.sample.graph(P0, directed = directed)
      random_est <- est.subgraph(
        G = G_star,
        H = H,
        B = B_null_variance,
        directed = directed,
        seed = NULL,
        order = order,
        aut_count = automorphism_count_value,
        keep_mappings = FALSE,
        verbose = FALSE,
        normalize_weight = FALSE
      )
      calibrated_null_counts[i] <- motif.est.value(random_est)
      calibrated_null_ses[i] <- motif.est.se(random_est)
      if(verbose && K_null_variance >= 10 &&
          i %% max(1, floor(K_null_variance/10)) == 0){
        cat("Finished", i, "of", K_null_variance, "null-variance calibrations\n")
      }
    }
    calibrated_tau2 <- stats::var(calibrated_null_counts)
  }else if(null_variance == "sample"){
    if(!exists("est.variance", mode = "function")){
      stop("Source VISS.R before requesting sampled null variance.")
    }
    sample_null_variance <- est.variance(
      P = P0,
      H = H,
      B = B_null_variance,
      directed = directed,
      induced = TRUE,
      normalize_weight = FALSE,
      overlap_sizes = null_variance_overlap_sizes,
      allocation = null_variance_allocation,
      include_incompatible = TRUE,
      mu_unlabeled = motif.est.value(null_expectation),
      aut_count = automorphism_count_value,
      seed = if(is.null(seed)) NULL else seed + 3L,
      verbose = verbose
    )
  }
  theta <- motif.est.value(observed)
  theta_se <- motif.est.se(observed)
  mu0 <- motif.est.value(null_expectation)
  mu0_se <- motif.est.se(null_expectation)
  .tsge.result(
    call = match.call(),
    theta = theta,
    theta_se = theta_se,
    mu0 = mu0,
    mu0_se = mu0_se,
    alternative = alternative,
    null_model = null_model,
    null_variance = null_variance,
    tau2_null = tau2_null,
    probability_rule = if(normalize_weight){"P_obs,ij = W_ij/max(W)"}else{"P_obs,ij = W_ij"},
    motif_order = motif.est.order(observed),
    automorphism_count = automorphism_count_value,
    observed = observed,
    null_expectation = null_expectation,
    sample_null_variance = sample_null_variance,
    calibrated_null_counts = calibrated_null_counts,
    calibrated_null_ses = calibrated_null_ses,
    calibrated_tau2 = calibrated_tau2,
    result_class = c("tsge.weighted", "tsge")
  )
}

.tsge.weighted.cal <- function(W, H,
                                    enrichment_result = NULL,
                                    observed = NULL,
                                    null_expectation = NULL,
                                    K = 99,
                                    R_null = 10000,
                                    B_observed = 10000,
                                    directed = TRUE,
                                    normalize_weight = TRUE,
                                    null_model = NULL,
                                    P_null = NULL,
                                    statistic = c("log_enrichment", "count", "excess"),
                                    alternative = c("greater", "less", "two.sided"),
                                    seed = NULL,
                                    order = NULL,
                                    automorphism_count_value = NULL,
                                    null_expectation_value = NULL,
                                    log_offset = 0,
                                    verbose = TRUE,
                                    verbose.K = TRUE){
  statistic <- match.arg(statistic)
  alternative <- match.arg(alternative)
  if(!is.numeric(K) || length(K) != 1 || is.na(K) || K < 1 || K != as.integer(K)){
    stop("K must be a positive integer.")
  }
  if(!is.numeric(R_null) || length(R_null) != 1 || is.na(R_null) ||
      R_null < 1 || R_null != as.integer(R_null)){
    stop("R_null must be a positive integer.")
  }
  K <- as.integer(K)
  R_null <- as.integer(R_null)
  if(!is.null(seed)){set.seed(seed)}
  P_obs <- motif.matrix.to.probmat.max(
    W,
    directed = directed,
    normalize_weight = normalize_weight,
    preprocessed = FALSE
  )
  H <- motif.graph.to.matrix(H, weighted = FALSE)
  H <- motif.norm.motif(H, directed = directed)
  observed_result <- observed
  null_expectation_result <- null_expectation
  observed_count <- NULL
  observed_count_se <- NULL
  if(!is.null(enrichment_result)){
    if(is.null(null_model)){null_model <- enrichment_result$null_model}
    if(is.null(order)){order <- enrichment_result$motif_order}
    if(is.null(automorphism_count_value)){
      automorphism_count_value <- enrichment_result$automorphism_count
    }
    if(is.null(observed_result)){
      observed_result <- enrichment_result$observed
    }
    if(is.null(null_expectation_result)){
      null_expectation_result <- enrichment_result$null_expectation
    }
    if(is.null(null_expectation_value)){
      null_expectation_value <- enrichment_result$estimates$null_expectation
    }
    if(is.null(observed_result)){
      observed_count <- enrichment_result$estimates$observed_count
      observed_count_se <- enrichment_result$estimates$observed_se
    }
  }
  if(!is.null(observed_result)){
    observed_count <- motif.est.value(observed_result)
    observed_count_se <- motif.est.se(observed_result)
    if(is.null(order)){
      order <- motif.est.order(observed_result)
    }
  }
  if(!is.null(null_expectation_result)){
    null_expectation_value <- motif.est.value(null_expectation_result)
  }
  if(is.null(null_model)){null_model <- "weighted_chung_lu"}
  null_model <- match.arg(null_model, c("weighted_chung_lu",
                                        "weighted_erdos_renyi",
                                        "given",
                                        "chung_lu",
                                        "erdos_renyi"))
  if(is.null(order)){
    order <- motif.order(
      H, motif.prob.density(P_obs, directed = directed), directed = directed
    )
  }
  if(is.null(automorphism_count_value)){
    automorphism_count_value <- motif.automorphism.count(H, directed = directed)
  }
  if(is.null(observed_count)){
    observed_result <- est.subgraph(
      G = P_obs,
      H = H,
      network_type = "weighted",
      B = B_observed,
      directed = directed,
      seed = if(is.null(seed)) NULL else seed + 1L,
      order = order,
      aut_count = automorphism_count_value,
      keep_mappings = FALSE,
      verbose = verbose,
      normalize_weight = FALSE
    )
    observed_count <- motif.est.value(observed_result)
    observed_count_se <- motif.est.se(observed_result)
    order <- motif.est.order(observed_result)
  }
  P0 <- if(null_model == "given"){
    if(is.null(P_null)){stop("P_null must be supplied when null_model = 'given'.")}
    motif.norm.probmat(P_null, directed = directed, strict = FALSE)
  }else{
    make.null.weighted(P_obs, method = null_model, directed = directed)
  }
  random_counts <- numeric(K)
  random_ses <- numeric(K)
  for(i in seq_len(K)){
    G_star <- motif.sample.graph(P0, directed = directed)
    random_est <- est.subgraph(
      G = G_star,
      H = H,
      B = R_null,
      directed = directed,
      seed = NULL,
      order = order,
      aut_count = automorphism_count_value,
      keep_mappings = FALSE,
      verbose = verbose,
      normalize_weight = FALSE
    )
    random_counts[i] <- motif.est.value(random_est)
    random_ses[i] <- motif.est.se(random_est)
    if(verbose.K && K >= 10 && i %% max(1, floor(K/10)) == 0){
      cat("Finished", i, "of", K, "weighted random-network calibrations\n")
    }
  }
  if(is.null(null_expectation_value)){null_expectation_value <- mean(random_counts)}
  .tsge.cal.result(
    call = match.call(),
    observed_count = observed_count,
    observed_count_se = observed_count_se,
    null_expectation = null_expectation_value,
    random_counts = random_counts,
    random_ses = random_ses,
    statistic = statistic,
    alternative = alternative,
    null_model = null_model,
    probability_rule = if(normalize_weight){"P_obs,ij = W_ij/max(W)"}else{"P_obs,ij = W_ij"},
    K = K,
    R_null = R_null,
    motif_order = order,
    automorphism_count = automorphism_count_value,
    log_offset = log_offset,
    observed_result = observed_result,
    null_expectation_result = null_expectation_result,
    result_class = c("tsge.weighted.cal", "tsge.cal")
  )
}

tsge <- function(G, H,
                 network_type = c("auto", "binary", "weighted"),
                 normalize_weight = TRUE, ...){
  network_type <- motif.network.type(G, network_type = network_type)
  result <- if(network_type == "binary"){
    .tsge.binary(G = G, H = H, ...)
  }else{
    .tsge.weighted(
      W = G,
      H = H,
      normalize_weight = normalize_weight,
      ...
    )
  }
  result$network_type <- network_type
  class(result) <- "tsge"
  result
}

tsge.cal <- function(G, H,
                     network_type = c("auto", "binary", "weighted"),
                     normalize_weight = TRUE, ...){
  network_type <- motif.network.type(G, network_type = network_type)
  result <- if(network_type == "binary"){
    .tsge.binary.cal(G = G, H = H, ...)
  }else{
    .tsge.weighted.cal(
      W = G,
      H = H,
      normalize_weight = normalize_weight,
      ...
    )
  }
  result$network_type <- network_type
  class(result) <- "tsge.cal"
  result
}
