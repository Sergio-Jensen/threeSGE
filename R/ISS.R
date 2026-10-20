#==========================================================================
## Induced-Sequential-Sampling (ISS) subgraph-count estimation
#==========================================================================
bootCI <- function(weights = NULL, log_weights = NULL, scale = 1,
                                  B_boot = 1000, level = 0.95, seed = NULL,
                                  verbose = TRUE){
  if(!is.null(weights)){
    n <- length(weights)
  }else if(!is.null(log_weights)){
    n <- length(log_weights)
  }else{
    stop("Either weights or log_weights must be supplied.")
  }
  if(n < 1 || B_boot < 1){return(c(NA_real_, NA_real_))}
  if(!is.null(seed)){set.seed(seed)}
  estimates <- numeric(B_boot)
  for(b in seq_len(B_boot)){
    idx <- sample.int(n, n, replace = TRUE)
    estimates[b] <- if(!is.null(weights)){
      mean(weights[idx])/scale
    }else{
      motif.mean.logs(log_weights[idx])/scale
    }
    if(verbose && B_boot >= 10 && b %% max(1, floor(B_boot/10)) == 0){
      cat("Finished", b, "of", B_boot, "Bootstrap samples\n")
    }
  }
  stats::quantile(estimates, probs = c((1 - level)/2, 1 - (1 - level)/2),
                  na.rm = TRUE, names = FALSE, type = 7)
}
valid.o <- function(order, H, name = "order"){
  k <- nrow(H)
  if(!is.numeric(order) && !is.integer(order)){
    stop(name, " must be a numeric permutation of motif node indices.")
  }
  if(length(order) != k || !setequal(order, seq_len(k))){
    stop(name, " must be a permutation of motif node indices.")
  }
  as.integer(order)
}
orders <- function(H, density, directed = FALSE, order = NULL,
                              order_method = c("mot.order", "reverse", "random")){
  if(!is.null(order)){
    return(valid.o(order, H, name = "order"))
  }
  order_method <- match.arg(order_method)
  base_order <- motif.order(H, density, directed = directed)
  if(order_method == "reverse"){
    return(rev(base_order))
  }
  if(order_method == "random"){
    return(sample(seq_len(nrow(H))))
  }
  base_order
}
bitsets <- function(G){
  n <- nrow(G)
  nonG <- !G
  diag(nonG) <- FALSE
  return(list(adj = G, non_adj = nonG, all = rep(TRUE, n)))
}
get.cand <- function(H, order, step, mapping, used, masks, directed = FALSE){
  subgraph_node <- order[step]
  candidates <- masks$all & !used
  if(step == 1){return(which(candidates))}
  for(prev_step in seq_len(step - 1)){
    prev_subgraph <- order[prev_step]
    prev_graph <- mapping[prev_step]
    if(directed){
      if(H[subgraph_node, prev_subgraph]){
        candidates <- candidates & masks$adj[, prev_graph]
      }else{
        candidates <- candidates & masks$non_adj[, prev_graph]
      }
      if(H[prev_subgraph, subgraph_node]){
        candidates <- candidates & masks$adj[prev_graph, ]
      }else{
        candidates <- candidates & masks$non_adj[prev_graph, ]
      }
    }else{
      if(H[subgraph_node, prev_subgraph]){
        candidates <- candidates & masks$adj[, prev_graph]
      }else{
        candidates <- candidates & masks$non_adj[, prev_graph]
      }
    }
    if(!any(candidates)){
      return(integer(0))
    }
  }
  which(candidates)
}
sample.iss <- function(H, order, masks, directed = FALSE){
  n <- nrow(masks$adj)
  k <- nrow(H)
  mapping <- integer(k)
  used <- rep(FALSE, n)
  weight <- 1
  for(step in seq_len(k)){
    candidates <- get.cand(
      H = H,
      order = order,
      step = step,
      mapping = mapping,
      used = used,
      masks = masks,
      directed = directed
    )
    candidate_count <- length(candidates)
    if(candidate_count == 0){
      return(list(weight = 0, mapping = integer(0), success = FALSE))
    }
    chosen <- if(candidate_count == 1){candidates}else{sample(candidates, 1)}
    mapping[step] <- chosen
    used[chosen] <- TRUE
    weight <- weight * candidate_count
  }
  full_mapping <- integer(k)
  full_mapping[order] <- mapping
  return(list(weight = weight, mapping = full_mapping, success = TRUE))
}
.est.subgraph.binary <- function(G, H, B = 10000, directed = FALSE, seed = NULL, order = NULL,
                         order_method = c("mot.order", "reverse", "random"),
                         keep_mappings = FALSE, verbose = TRUE, aut_count = NULL,
                         bootstrap_ci = FALSE, bootstrap_B = 1000, verbose.boot = FALSE,
                         bootstrap_seed = NULL){
  if(!is.null(seed)){set.seed(seed)}
  G <- motif.norm.graph(G, directed = directed)
  H <- motif.norm.graph(H, directed = directed)
  if(nrow(G) != ncol(G)){stop("G must be a square adjacency matrix.")}
  if(nrow(H) != ncol(H)){stop("H must be a square subgraph adjacency matrix.")}
  if(nrow(H) > nrow(G)){stop("The subgraph cannot have more nodes than the graph.")}
  density <- motif.graph.density(G, directed = directed)
  resolved_order_method <- if(is.null(order)){
    match.arg(order_method)
  }else{
    "given"
  }
  order <- orders(
    H = H,
    density = density,
    directed = directed,
    order = order,
    order_method = if(resolved_order_method == "given") "mot.order" else resolved_order_method
  )
  masks <- bitsets(G)
  weights <- numeric(B)
  successes <- logical(B)
  mappings <- if(keep_mappings){vector("list", B)}else{NULL}
  for(b in seq_len(B)){
    draw <- sample.iss(H, order, masks, directed = directed)
    weights[b] <- draw$weight
    successes[b] <- draw$success
    if(keep_mappings && draw$success){
      mappings[[b]] <- draw$mapping
    }
    if(verbose && B >= 10 && b %% max(1, floor(B/10)) == 0){
      cat("Finished", b, "of", B, "samples\n")
    }
  }
  labeled_estimate <- mean(weights)
  labeled_se <- if(B > 1){stats::sd(weights)/sqrt(B)}else{NA_real_}
  if(is.null(aut_count)){
    aut_count <- motif.automorphism.count(H, directed = directed)
  }
  unlabeled_estimate <- labeled_estimate/aut_count
  unlabeled_se <- labeled_se/aut_count
  labeled_wald_ci_95 <- motif.wald.ci(labeled_estimate, labeled_se)
  unlabeled_wald_ci_95 <- motif.wald.ci(unlabeled_estimate, unlabeled_se)
  labeled_log_wald_ci_95 <- motif.log.wald.ci(labeled_estimate, labeled_se)
  unlabeled_log_wald_ci_95 <- motif.log.wald.ci(unlabeled_estimate, unlabeled_se)
  diagnostics <- motif.weight.diagnostics(weights)
  bootstrap <- NULL
  if(isTRUE(bootstrap_ci)){
    bootstrap <- list(
      labeled_ci_95 = bootCI(
        weights = weights,
        scale = 1,
        B_boot = bootstrap_B,
        seed = bootstrap_seed,
        verbose = verbose.boot
      ),
      unlabeled_ci_95 = bootCI(
        weights = weights,
        scale = aut_count,
        B_boot = bootstrap_B,
        seed = if(is.null(bootstrap_seed)) NULL else bootstrap_seed + 1L,
        verbose = verbose.boot
      ),
      B = bootstrap_B,
      method = "nonparametric bootstrap"
    )
  }
  result <- list(B = B,
    directed = directed,
    graph_nodes = nrow(G),
    subgraph_nodes = nrow(H),
    graph_density = density,
    subgraph_order = order,
    order_method = resolved_order_method,
    aut.morph = aut_count,
    successes = sum(successes),
    success_rate = mean(successes),
    labeled_estimate = labeled_estimate,
    labeled_se = labeled_se,
    labeled_ci_95 = labeled_wald_ci_95,
    labeled_wald_ci_95 = labeled_wald_ci_95,
    labeled_log_wald_ci_95 = labeled_log_wald_ci_95,
    unlabeled_estimate = unlabeled_estimate,
    unlabeled_se = unlabeled_se,
    unlabeled_ci_95 = unlabeled_wald_ci_95,
    unlabeled_wald_ci_95 = unlabeled_wald_ci_95,
    unlabeled_log_wald_ci_95 = unlabeled_log_wald_ci_95,
    ci_method = "Wald",
    supplemental_ci_method = "log-scale Wald",
    weight_diagnostics = diagnostics,
    bootstrap_ci = bootstrap,
    weights = weights)
  if(keep_mappings){result$mappings <- mappings}
  class(result) <- "est.subgraph"
  return(result)
}

print.est.subgraph <- function(x, ...){
  if(identical(x$network_type, "weighted")){
    cat("Weighted induced-SIS expectation estimate\n")
    cat("  graph nodes:      ", x$graph_nodes, "\n", sep = "")
    cat("  motif nodes:      ", x$motif_nodes, "\n", sep = "")
    cat("  samples:          ", x$B, "\n", sep = "")
    cat("  directed:         ", x$directed, "\n", sep = "")
    cat("  normalize weight: ", x$normalize_weight, "\n", sep = "")
    cat("  probability rule: ", x$probability_rule, "\n", sep = "")
    cat("  probability dens: ", signif(x$probability_density, 6), "\n", sep = "")
    cat("  motif order:      ", paste(x$motif_order, collapse = " -> "), "\n", sep = "")
    cat("  automorphisms:    ", x$aut.morph, "\n", sep = "")
    cat("  labeled expected: ", signif(x$labeled_expectation, 6),
        " (SE ", signif(x$labeled_se, 4), ")\n", sep = "")
    cat("  unlabeled expected: ", signif(x$unlabeled_expectation, 6),
        " (SE ", signif(x$unlabeled_se, 4), ")\n", sep = "")
    return(invisible(x))
  }
  cat("Induced-Sequential-Sampling estimate\n")
  cat("  graph nodes:      ", x$graph_nodes, "\n", sep = "")
  cat("  subgraph nodes:      ", x$subgraph_nodes, "\n", sep = "")
  cat("  samples:          ", x$B, "\n", sep = "")
  cat("  success rate:     ", signif(x$success_rate, 4), "\n", sep = "")
  cat("  graph density:    ", signif(x$graph_density, 4), "\n", sep = "")
  cat("  subgraph order:      ", paste(x$subgraph_order, collapse = " -> "), "\n", sep = "")
  cat("  order method:     ", x$order_method, "\n", sep = "")
  cat("  automorphisms:    ", x$aut.morph, "\n", sep = "")
  cat("  labeled estimate: ", signif(x$labeled_estimate, 6),
      " (SE ", signif(x$labeled_se, 4), ")\n", sep = "")
  cat("  labeled Wald 95% CI: [", signif(x$labeled_wald_ci_95[1], 6),
      ", ", signif(x$labeled_wald_ci_95[2], 6), "]\n", sep = "")
  cat("  unlabeled count:  ", signif(x$unlabeled_estimate, 6),
      " (SE ", signif(x$unlabeled_se, 4), ")\n", sep = "")
  cat("  unlabeled Wald 95% CI: [", signif(x$unlabeled_wald_ci_95[1], 6),
      ", ", signif(x$unlabeled_wald_ci_95[2], 6), "]\n", sep = "")
  if(!is.null(x$weight_diagnostics)){
    cat("  ESS ratio:        ", signif(x$weight_diagnostics$ess_ratio, 4), "\n", sep = "")
    cat("  max weight share: ", signif(x$weight_diagnostics$max_weight_share, 4), "\n", sep = "")
  }
  if(motif.rare.case(x$weight_diagnostics, success_rate = x$success_rate) &&
      all(is.finite(x$unlabeled_log_wald_ci_95))){
    cat("  supplemental log-Wald unlabeled 95% CI: [",
        signif(x$unlabeled_log_wald_ci_95[1], 6), ", ",
        signif(x$unlabeled_log_wald_ci_95[2], 6), "]\n", sep = "")
  }
  if(!is.null(x$bootstrap_ci)){
    cat("  bootstrap unlabeled 95% CI: [", signif(x$bootstrap_ci$unlabeled_ci_95[1], 6),
        ", ", signif(x$bootstrap_ci$unlabeled_ci_95[2], 6), "]\n", sep = "")
  }
  invisible(x)
}
iss.sensitivity <- function(G, H, B = 10000, directed = FALSE, seed = NULL,
                                  order = NULL,
                                  orders = c("mot.order", "reverse"),
                                  random_orders = 0,
                                  aut_count = NULL,
                                  verbose = FALSE,
                                  bootstrap_ci = FALSE,
                                  bootstrap_B = 1000,
                                  verbose.boot = FALSE,
                                  keep_results = TRUE,
                                  verbose.sens = TRUE){
  if(!is.null(seed)){set.seed(seed)}
  G <- motif.norm.graph(G, directed = directed)
  H <- motif.norm.graph(H, directed = directed)
  density <- motif.graph.density(G, directed = directed)
  if(is.null(aut_count)){
    aut_count <- motif.automorphism.count(H, directed = directed)
  }
  order_list <- list()
  add_order <- function(name, value){
    order_list[[name]] <<- valid.o(value, H, name = name)
    invisible(NULL)
  }
  base_order <- if(is.null(order)){
    motif.order(H, density, directed = directed)
  }else{
    valid.o(order, H, name = "order")
  }
  if(length(orders) > 0){
    for(spec in orders){
      if(spec %in% c("mot.order", "default", "rarity")){
        add_order("mot.order", base_order)
      }else if(spec == "reverse"){
        add_order("reverse", rev(base_order))
      }else if(spec == "given"){
        add_order("given", base_order)
      }else{
        stop("Unknown ISS order sensitivity option: ", spec)
      }
    }
  }
  if(random_orders > 0){
    for(i in seq_len(random_orders)){
      add_order(paste0("random_", i), sample(seq_len(nrow(H))))
    }
  }
  if(length(order_list) == 0){stop("No valid ISS orders were supplied.")}
  results <- if(keep_results){vector("list", length(order_list))}else{NULL}
  if(keep_results){names(results) <- names(order_list)}
  table <- vector("list", length(order_list))
  i <- 0L
  for(name in names(order_list)){
    i <- i + 1L
    res <- est.subgraph(
      G = G,
      H = H,
      B = B,
      directed = directed,
      seed = if(is.null(seed)) NULL else seed + 100L * i,
      order = order_list[[name]],
      keep_mappings = FALSE,
      verbose = verbose,
      aut_count = aut_count,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      verbose.boot = verbose.boot,
      bootstrap_seed = if(is.null(seed)) NULL else seed + 10000L + 100L * i
    )
    if(keep_results){results[[name]] <- res}
    table[[i]] <- data.frame(
      order_name = name,
      order = paste(res$subgraph_order, collapse = " -> "),
      unlabeled_estimate = res$unlabeled_estimate,
      unlabeled_se = res$unlabeled_se,
      wald_lower = res$unlabeled_wald_ci_95[1],
      wald_upper = res$unlabeled_wald_ci_95[2],
      log_wald_lower = res$unlabeled_log_wald_ci_95[1],
      log_wald_upper = res$unlabeled_log_wald_ci_95[2],
      success_rate = res$success_rate,
      ess_ratio = res$weight_diagnostics$ess_ratio,
      cv_weight = res$weight_diagnostics$cv_weight,
      max_weight_share = res$weight_diagnostics$max_weight_share,
      row.names = NULL
    )
    if(verbose.sens){
      cat("Finished", name, "sensitivity test\n")
    }
  }
  out <- list(
    call = match.call(),
    sensitivity_table = do.call(rbind, table),
    results = results,
    B = B,
    directed = directed,
    graph_density = density,
    random_orders = random_orders
  )
  class(out) <- "iss.sensitivity"
  out
}
print.iss.sensitivity <- function(x, ...){
  cat("ISS order sensitivity analysis\n")
  cat("  samples:       ", x$B, "\n", sep = "")
  cat("  graph density: ", signif(x$graph_density, 4), "\n", sep = "")
  cat("\nSensitivity table\n")
  print(x$sensitivity_table, row.names = FALSE)
  invisible(x)
}

.weighted.candidate.logp <- function(P, H, order, step, mapping, candidates,
                                     directed = TRUE){
  if(step == 1){return(rep(0, length(candidates)))}
  motif_node <- order[step]
  log_prob <- numeric(length(candidates))
  for(prev_step in seq_len(step - 1)){
    prev_motif <- order[prev_step]
    prev_graph <- mapping[prev_step]
    if(directed){
      log_prob <- log_prob + motif.edge.logprob(
        P[candidates, prev_graph], H[motif_node, prev_motif]
      )
      log_prob <- log_prob + motif.edge.logprob(
        P[prev_graph, candidates], H[prev_motif, motif_node]
      )
    }else{
      log_prob <- log_prob + motif.edge.logprob(
        P[candidates, prev_graph], H[motif_node, prev_motif]
      )
    }
  }
  log_prob
}

.sample.iss.weighted <- function(P, H, order, directed = TRUE){
  n <- nrow(P)
  k <- nrow(H)
  mapping <- integer(k)
  used <- rep(FALSE, n)
  log_weight <- 0
  for(step in seq_len(k)){
    candidates <- which(!used)
    log_g <- .weighted.candidate.logp(
      P = P,
      H = H,
      order = order,
      step = step,
      mapping = mapping,
      candidates = candidates,
      directed = directed
    )
    log_z <- motif.log.sum(log_g)
    if(!is.finite(log_z)){
      return(list(log_weight = -Inf, mapping = integer(0), success = FALSE))
    }
    probs <- exp(log_g - log_z)
    chosen <- if(length(candidates) == 1){
      candidates
    }else{
      sample(candidates, 1, prob = probs)
    }
    mapping[step] <- chosen
    used[chosen] <- TRUE
    log_weight <- log_weight + log_z
  }
  full_mapping <- integer(k)
  full_mapping[order] <- mapping
  list(log_weight = log_weight, mapping = full_mapping, success = TRUE)
}

.est.subgraph.weighted <- function(W, H, B = 10000, directed = TRUE,
                                   seed = NULL, order = NULL,
                                   keep_mappings = FALSE, verbose = TRUE,
                                   aut_count = NULL,
                                   normalize_weight = TRUE){
  if(!is.numeric(B) || length(B) != 1 || is.na(B) ||
     B < 1 || B != as.integer(B)){
    stop("B must be a positive integer.")
  }
  B <- as.integer(B)
  if(!is.null(seed)){set.seed(seed)}
  P <- motif.matrix.to.probmat(
    W,
    directed = directed,
    normalize_weight = normalize_weight
  )
  H <- motif.norm.motif(H, directed = directed)
  if(nrow(H) > nrow(P)){stop("H cannot have more nodes than W.")}
  density <- motif.prob.density(P, directed = directed)
  if(is.null(order)){
    order <- motif.order(H, density, directed = directed)
  }
  order <- valid.o(order, H, name = "order")
  if(is.null(aut_count)){
    aut_count <- motif.automorphism.count(H, directed = directed)
  }
  log_weights <- numeric(B)
  successes <- logical(B)
  mappings <- if(keep_mappings){vector("list", B)}else{NULL}
  for(b in seq_len(B)){
    draw <- .sample.iss.weighted(P, H, order, directed = directed)
    log_weights[b] <- draw$log_weight
    successes[b] <- draw$success
    if(keep_mappings && draw$success){mappings[[b]] <- draw$mapping}
    if(verbose && B >= 10 && b %% max(1, floor(B/10)) == 0){
      cat("Finished", b, "of", B, "weighted samples\n")
    }
  }
  labeled_estimate <- motif.mean.logs(log_weights)
  labeled_var <- motif.var.logs(log_weights)
  labeled_se <- sqrt(labeled_var/B)
  unlabeled_estimate <- labeled_estimate/aut_count
  unlabeled_se <- labeled_se/aut_count
  result <- list(
    B = B,
    directed = directed,
    graph_nodes = nrow(P),
    motif_nodes = nrow(H),
    probability_density = density,
    motif_order = order,
    aut.morph = aut_count,
    normalize_weight = normalize_weight,
    probability_rule = if(normalize_weight){
      "P_ij = W_ij / max(W)"
    }else{
      "P_ij = W_ij"
    },
    successes = sum(successes),
    success_rate = mean(successes),
    labeled_expectation = labeled_estimate,
    labeled_se = labeled_se,
    unlabeled_expectation = unlabeled_estimate,
    unlabeled_se = unlabeled_se,
    log_weights = log_weights,
    probability_matrix = P
  )
  if(keep_mappings){result$mappings <- mappings}
  result
}

est.subgraph <- function(G, H,
                         network_type = c("auto", "binary", "weighted"),
                         normalize_weight = TRUE, ...){
  network_type <- motif.network.type(G, network_type = network_type)
  result <- if(network_type == "binary"){
    .est.subgraph.binary(G = G, H = H, ...)
  }else{
    .est.subgraph.weighted(
      W = G,
      H = H,
      normalize_weight = normalize_weight,
      ...
    )
  }
  result$network_type <- network_type
  class(result) <- "est.subgraph"
  result
}
