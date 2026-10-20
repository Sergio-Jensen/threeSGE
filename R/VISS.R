#==========================================================================
## Variance ISS (VISS) subgraph-count variance estimation
#==========================================================================
sm.falling <- function(n, k){
  if(k < 0 || k > n){return(0)}
  if(k == 0){return(1)}
  prod(n - seq.int(0, k - 1))
}

sm.layer <- function(n, k, s){
  sm.falling(n, k) * choose(k, s) * sm.falling(k, s) *
    sm.falling(n - k, k - s)
}

sm.log.var <- function(log_w){
  finite <- is.finite(log_w)
  if(sum(finite) <= 1){return(NA_real_)}
  m <- max(log_w[finite])
  w <- exp(log_w[finite] - m)
  stats::var(w) * exp(2 * m)
}

sm.strata.diagns <- function(weight_layers, use_abs = FALSE){
  if(length(weight_layers) == 0L){
    return(list(pooled = motif.weight.diagnostics(numeric(0)), by_layer = data.frame()))
  }
  clean_layers <- lapply(weight_layers, function(w){
    w <- as.numeric(w)
    if(use_abs){w <- abs(w)}
    w
  })
  by_layer <- do.call(rbind, lapply(seq_along(clean_layers), function(i){
    d <- motif.weight.diagnostics(clean_layers[[i]])
    data.frame(
      layer = i,
      n = d$n,
      finite = d$finite,
      zero_rate = d$zero_rate,
      positive_rate = d$positive_rate,
      ess = d$ess,
      ess_ratio = d$ess_ratio,
      cv_weight = d$cv_weight,
      max_weight_share = d$max_weight_share,
      row.names = NULL
    )
  }))
  pooled <- motif.weight.diagnostics(unlist(clean_layers, use.names = FALSE))
  list(
    pooled = pooled,
    by_layer = by_layer,
    use_abs = use_abs,
    min_ess_ratio = if(nrow(by_layer) > 0){min(by_layer$ess_ratio, na.rm = TRUE)}else{NA_real_},
    min_positive_rate = if(nrow(by_layer) > 0){min(by_layer$positive_rate, na.rm = TRUE)}else{NA_real_},
    max_weight_share = if(nrow(by_layer) > 0){max(by_layer$max_weight_share, na.rm = TRUE)}else{NA_real_}
  )
}

sm.bootCI.strata <- function(weight_layers, B_boot = 1000, level = 0.95,
                             seed = NULL, verbose = FALSE){
  if(length(weight_layers) == 0L || B_boot < 1){
    return(c(NA_real_, NA_real_))
  }
  if(!is.null(seed)){set.seed(seed)}
  estimates <- numeric(B_boot)
  for(b in seq_len(B_boot)){
    estimates[b] <- sum(vapply(weight_layers, function(w){
      n <- length(w)
      if(n < 1){return(0)}
      mean(w[sample.int(n, n, replace = TRUE)])
    }, numeric(1)))
    if(verbose && B_boot >= 10 && b %% max(1, floor(B_boot/10)) == 0){
      cat("Finished", b, "of", B_boot, "Bootstrap samples\n")
    }
  }
  stats::quantile(estimates, probs = c((1 - level)/2, 1 - (1 - level)/2),
                  na.rm = TRUE, names = FALSE, type = 7)
}

sm.inference <- function(weight_layers, estimate, se, level = 0.95,
                         bootstrap_ci = FALSE, bootstrap_B = 1000,
                         bootstrap_seed = NULL, verbose.boot = FALSE,
                         use_abs_diagnostics = FALSE){
  diagnostics <- sm.strata.diagns(weight_layers, use_abs = use_abs_diagnostics)
  boot <- NULL
  if(isTRUE(bootstrap_ci)){
    boot <- sm.bootCI.strata(
      weight_layers = weight_layers,
      B_boot = bootstrap_B,
      level = level,
      seed = bootstrap_seed,
      verbose = verbose.boot
    )
  }
  list(
    estimate = estimate,
    se = se,
    level = level,
    wald_ci = motif.wald.ci(estimate, se, level = level),
    log_wald_ci = motif.log.wald.ci(estimate, se, level = level),
    ci_method = "Wald",
    supplemental_ci_method = "log-scale Wald",
    bootstrap_ci = boot,
    bootstrap_method = if(isTRUE(bootstrap_ci)){"stratified nonparametric bootstrap"}else{NULL},
    weight_diagnostics = diagnostics,
    rare_case = motif.rare.case(diagnostics$pooled)
  )
}

sm.add.count <- function(edge_counts, nonedge_counts, i, j, require_edge,
                         directed = FALSE){
  if(i == j){return(list(edge = edge_counts, nonedge = nonedge_counts))}
  if(!directed && i > j){tmp <- i; i <- j; j <- tmp}
  if(require_edge){
    edge_counts[i, j] <- edge_counts[i, j] + 1L
    if(!directed){edge_counts[j, i] <- edge_counts[j, i] + 1L}
  }else{
    nonedge_counts[i, j] <- nonedge_counts[i, j] + 1L
    if(!directed){nonedge_counts[j, i] <- nonedge_counts[j, i] + 1L}
  }
  list(edge = edge_counts, nonedge = nonedge_counts)
}

sm.add.joint <- function(C, i, j, require_edge, directed = FALSE){
  if(i == j){return(C)}
  if(!directed && i > j){tmp <- i; i <- j; j <- tmp}
  new_value <- if(require_edge){1}else{0}
  old <- C[i, j]
  if(!is.na(old) && old != new_value){
    attr(C, "conflict") <- TRUE
    return(C)
  }
  C[i, j] <- new_value
  if(!directed){C[j, i] <- new_value}
  C
}

sm.union.pattern <- function(H, phi_for_psi, directed = FALSE, induced = TRUE){
  k <- nrow(H)
  phi_nodes <- seq_len(k)
  psi_nodes <- integer(k)
  next_node <- k
  for(b in seq_len(k)){
    if(!is.na(phi_for_psi[b])){
      psi_nodes[b] <- phi_for_psi[b]
    }else{
      next_node <- next_node + 1L
      psi_nodes[b] <- next_node
    }
  }
  m <- next_node
  joint <- matrix(NA_real_, m, m)
  edge_counts <- matrix(0L, m, m)
  nonedge_counts <- matrix(0L, m, m)
  add.copy <- function(nodes){
    if(directed){
      for(a in seq_len(k)){
        for(b in seq_len(k)){
          if(a == b){next}
          if(H[a, b] || induced){
            cnt <- sm.add.count(edge_counts, nonedge_counts, nodes[a], nodes[b],
                                H[a, b], directed = TRUE)
            edge_counts <<- cnt$edge
            nonedge_counts <<- cnt$nonedge
            joint <<- sm.add.joint(joint, nodes[a], nodes[b], H[a, b],
                                   directed = TRUE)
          }
        }
      }
    }else if(k > 1){
      for(a in seq_len(k - 1)){
        for(b in seq.int(a + 1, k)){
          if(H[a, b] || induced){
            cnt <- sm.add.count(edge_counts, nonedge_counts, nodes[a], nodes[b],
                                H[a, b], directed = FALSE)
            edge_counts <<- cnt$edge
            nonedge_counts <<- cnt$nonedge
            joint <<- sm.add.joint(joint, nodes[a], nodes[b], H[a, b],
                                   directed = FALSE)
          }
        }
      }
    }
    invisible(NULL)
  }
  add.copy(phi_nodes)
  add.copy(psi_nodes)
  diag(joint) <- NA_real_
  diag(edge_counts) <- 0L
  diag(nonedge_counts) <- 0L
  list(
    conflict = isTRUE(attr(joint, "conflict")),
    constraints = joint,
    product_edge_counts = edge_counts,
    product_nonedge_counts = nonedge_counts,
    phi_nodes = phi_nodes,
    psi_nodes = psi_nodes,
    union_nodes = m
  )
}

sm.n4mat <- function(C, directed = FALSE){
  m <- nrow(C)
  edge_counts <- matrix(0L, m, m)
  nonedge_counts <- matrix(0L, m, m)
  if(directed){
    idx <- which(!is.na(C) & row(C) != col(C), arr.ind = TRUE)
  }else{
    idx <- which(!is.na(C) & upper.tri(C), arr.ind = TRUE)
  }
  if(nrow(idx) > 0){
    for(i in seq_len(nrow(idx))){
      a <- idx[i, 1]
      b <- idx[i, 2]
      cnt <- sm.add.count(edge_counts, nonedge_counts, a, b, C[a, b] == 1,
                          directed = directed)
      edge_counts <- cnt$edge
      nonedge_counts <- cnt$nonedge
    }
  }
  list(edge = edge_counts, nonedge = nonedge_counts)
}

sm.constns <- function(edge_counts, nonedge_counts, directed = FALSE){
  if(directed){
    idx <- which((edge_counts > 0L | nonedge_counts > 0L) &
                   row(edge_counts) != col(edge_counts), arr.ind = TRUE)
  }else{
    idx <- which((edge_counts > 0L | nonedge_counts > 0L) & upper.tri(edge_counts),
                 arr.ind = TRUE)
  }
  if(nrow(idx) == 0){
    return(list(from = integer(0), to = integer(0), edge_count = integer(0),
                nonedge_count = integer(0), n = 0L))
  }
  list(
    from = idx[, 1],
    to = idx[, 2],
    edge_count = as.integer(edge_counts[idx]),
    nonedge_count = as.integer(nonedge_counts[idx]),
    n = nrow(idx)
  )
}

sm.constnp <- function(P, mapping, constraints){
  if(constraints$n == 0L){return(1)}
  p <- P[cbind(mapping[constraints$from], mapping[constraints$to])]
  prod((p ^ constraints$edge_count) * ((1 - p) ^ constraints$nonedge_count))
}

sm.edge.count <- function(H, directed = FALSE){
  if(directed){
    sum(H[row(H) != col(H)])
  }else{
    sum(H[upper.tri(H)])
  }
}

sm.shared.edge.count <- function(compiled){
  constraints <- compiled$product_constraint_list
  if(is.null(constraints) || constraints$n == 0L){return(0L)}
  sum(constraints$edge_count >= 2L)
}

sm.consts <- function(H, directed = FALSE, induced = TRUE){
  k <- nrow(H)
  edge_counts <- matrix(0L, k, k)
  nonedge_counts <- matrix(0L, k, k)
  if(directed){
    for(a in seq_len(k)){
      for(b in seq_len(k)){
        if(a == b){next}
        if(H[a, b] || induced){
          cnt <- sm.add.count(edge_counts, nonedge_counts, a, b, H[a, b],
                              directed = TRUE)
          edge_counts <- cnt$edge
          nonedge_counts <- cnt$nonedge
        }
      }
    }
  }else if(k > 1){
    for(a in seq_len(k - 1)){
      for(b in seq.int(a + 1, k)){
        if(H[a, b] || induced){
          cnt <- sm.add.count(edge_counts, nonedge_counts, a, b, H[a, b],
                              directed = FALSE)
          edge_counts <- cnt$edge
          nonedge_counts <- cnt$nonedge
        }
      }
    }
  }
  sm.constns(edge_counts, nonedge_counts, directed = directed)
}

sm.partial.o <- function(C = NULL, directed = FALSE, graph_density = 0.5,
                         edge_counts = NULL, nonedge_counts = NULL,
                         eps = 1e-12){
  if(is.null(edge_counts) || is.null(nonedge_counts)){
    if(is.null(C)){stop("Either C or count matrices must be supplied.")}
    counts <- sm.n4mat(C, directed = directed)
    edge_counts <- counts$edge
    nonedge_counts <- counts$nonedge
  }
  m <- nrow(edge_counts)
  p <- min(max(graph_density, eps), 1 - eps)
  edge_w <- -log(p)
  nonedge_w <- -log1p(-p)
  score.mat <- edge_counts * edge_w + nonedge_counts * nonedge_w
  total_score <- vapply(seq_len(m), function(u){
    if(directed){
      sum(score.mat[u, -u, drop = TRUE]) + sum(score.mat[-u, u, drop = TRUE])
    }else{
      sum(score.mat[u, -u, drop = TRUE])
    }
  }, numeric(1))
  if(all(total_score == 0)){return(seq_len(m))}
  first <- which.max(total_score)
  order <- first
  remaining <- setdiff(seq_len(m), first)
  while(length(remaining) > 0){
    scores <- vapply(remaining, function(u){
      imposed <- sum(vapply(order, function(v){
        if(directed){score.mat[u, v] + score.mat[v, u]}else{score.mat[u, v]}
      }, numeric(1)))
      imposed + 1e-9 * total_score[u]
    }, numeric(1))
    next_node <- remaining[which.max(scores)]
    order <- c(order, next_node)
    remaining <- setdiff(remaining, next_node)
  }
  order
}

sm.permutations <- function(x){
  n <- length(x)
  if(n == 0){return(matrix(integer(0), nrow = 1))}
  if(n == 1){return(matrix(x, nrow = 1))}
  out <- NULL
  for(i in seq_along(x)){
    sub <- sm.permutations(x[-i])
    out <- rbind(out, cbind(x[i], sub))
  }
  out
}

sm.constraint.aut <- function(edge_counts, nonedge_counts, directed = FALSE){
  m <- nrow(edge_counts)
  if(m <= 1L){return(1L)}
  state <- matrix(NA_integer_, m, m)
  state[edge_counts > 0L] <- 1L
  state[nonedge_counts > 0L] <- 0L
  diag(state) <- NA_integer_
  if(!directed){
    state[lower.tri(state)] <- t(state)[lower.tri(state)]
  }
  perms <- sm.permutations(seq_len(m))
  sum(vapply(seq_len(nrow(perms)), function(i){
    p <- perms[i, ]
    identical(state, state[p, p, drop = FALSE])
  }, logical(1)))
}

sm.overlap <- function(k, s){
  if(s < 0 || s > k){stop("s must be between 0 and k.")}
  if(s == 0){return(list(rep(NA_integer_, k)))}
  out <- list()
  psi_subsets <- utils::combn(seq_len(k), s, simplify = FALSE)
  phi_subsets <- utils::combn(seq_len(k), s, simplify = FALSE)
  for(psi_pos in psi_subsets){
    for(phi_pos in phi_subsets){
      perm_mat <- if(s == 1){matrix(phi_pos, nrow = 1)}else{sm.permutations(phi_pos)}
      for(i in seq_len(nrow(perm_mat))){
        pattern <- rep(NA_integer_, k)
        pattern[psi_pos] <- perm_mat[i, ]
        out[[length(out) + 1L]] <- pattern
      }
    }
  }
  out
}

sm.compilen <- function(edge_counts, nonedge_counts, directed = FALSE,
                        graph_density = 0.5, phi_nodes = integer(0),
                        psi_nodes = integer(0), union_nodes = NULL,
                        compatible = TRUE, target = c("joint", "product")){
  target <- match.arg(target)
  m <- union_nodes %||% nrow(edge_counts)
  order <- sm.partial.o(
    directed = directed,
    graph_density = graph_density,
    edge_counts = edge_counts,
    nonedge_counts = nonedge_counts
  )
  step_constraints <- vector("list", m)
  for(step in seq_len(m)){
    if(step == 1){
      step_constraints[[step]] <- list(prev_step = integer(0), direction = integer(0),
                                       edge_count = integer(0), nonedge_count = integer(0),
                                       n = 0L)
      next
    }
    subgraph_node <- order[step]
    prev_steps <- integer(0)
    directions <- integer(0)
    edge_cs <- integer(0)
    nonedge_cs <- integer(0)
    for(prev_step in seq_len(step - 1)){
      prev_subgraph <- order[prev_step]
      ec <- edge_counts[subgraph_node, prev_subgraph]
      nc <- nonedge_counts[subgraph_node, prev_subgraph]
      if(ec > 0L || nc > 0L){
        prev_steps <- c(prev_steps, prev_step)
        directions <- c(directions, 1L)
        edge_cs <- c(edge_cs, ec)
        nonedge_cs <- c(nonedge_cs, nc)
      }
      if(directed){
        ec <- edge_counts[prev_subgraph, subgraph_node]
        nc <- nonedge_counts[prev_subgraph, subgraph_node]
        if(ec > 0L || nc > 0L){
          prev_steps <- c(prev_steps, prev_step)
          directions <- c(directions, -1L)
          edge_cs <- c(edge_cs, ec)
          nonedge_cs <- c(nonedge_cs, nc)
        }
      }
    }
    step_constraints[[step]] <- list(prev_step = prev_steps, direction = directions,
                                     edge_count = edge_cs, nonedge_count = nonedge_cs,
                                     n = length(prev_steps))
  }
  list(
    constraint_list = sm.constns(edge_counts, nonedge_counts,
                                 directed = directed),
    edge_counts = edge_counts,
    nonedge_counts = nonedge_counts,
    order = order,
    step_constraints = step_constraints,
    phi_nodes = phi_nodes,
    psi_nodes = psi_nodes,
    union_nodes = m,
    compatible = compatible,
    target = target
  )
}

sm.compileu <- function(union, directed = FALSE, graph_density = 0.5,
                        target = NULL){
  if(is.null(target)){target <- if(isTRUE(union$conflict)){"product"}else{"joint"}}
  if(target == "joint"){
    counts <- sm.n4mat(union$constraints, directed = directed)
    out <- sm.compilen(
      counts$edge, counts$nonedge, directed = directed,
      graph_density = graph_density, phi_nodes = union$phi_nodes,
      psi_nodes = union$psi_nodes, union_nodes = union$union_nodes,
      compatible = !isTRUE(union$conflict), target = "joint"
    )
  }else{
    out <- sm.compilen(
      union$product_edge_counts, union$product_nonedge_counts,
      directed = directed, graph_density = graph_density,
      phi_nodes = union$phi_nodes, psi_nodes = union$psi_nodes,
      union_nodes = union$union_nodes, compatible = !isTRUE(union$conflict),
      target = "product"
    )
  }
  out$product_constraint_list <- sm.constns(
    union$product_edge_counts, union$product_nonedge_counts, directed = directed
  )
  out
}

sm.compile.overlap <- function(H, s, directed = FALSE, induced = TRUE,
                               n = Inf, graph_density = 0.5,
                               include_incompatible = TRUE){
  raw_patterns <- sm.overlap(nrow(H), s)
  compiled <- list()
  for(pattern in raw_patterns){
    union <- sm.union.pattern(H = H, phi_for_psi = pattern,
                              directed = directed, induced = induced)
    if(union$union_nodes > n){next}
    if(isTRUE(union$conflict) && !include_incompatible){next}
    compiled[[length(compiled) + 1L]] <- sm.compileu(
      union, directed = directed, graph_density = graph_density,
      target = if(isTRUE(union$conflict)){"product"}else{"joint"}
    )
  }
  compiled
}

sm.pcand.logp <- function(P, compiled, step, mapping,
                          candidates){
  constraints <- compiled$step_constraints[[step]]
  if(constraints$n == 0L){return(rep(0, length(candidates)))}
  log_prob <- numeric(length(candidates))
  for(i in seq_len(constraints$n)){
    prev_graph <- mapping[constraints$prev_step[i]]
    p <- if(constraints$direction[i] == 1L){
      P[candidates, prev_graph]
    }else{
      P[prev_graph, candidates]
    }
    if(constraints$edge_count[i] > 0L){
      log_prob <- log_prob +
        constraints$edge_count[i] * motif.edge.logprob(p, TRUE)
    }
    if(constraints$nonedge_count[i] > 0L){
      log_prob <- log_prob +
        constraints$nonedge_count[i] * motif.edge.logprob(p, FALSE)
    }
  }
  log_prob
}

sm.sample.compile <- function(P, compiled, fixed_mapping = NULL){
  n <- nrow(P)
  m <- compiled$union_nodes
  order <- compiled$order
  mapping <- integer(m)
  used <- rep(FALSE, n)
  log_weight <- 0
  if(!is.null(fixed_mapping)){
    if(length(fixed_mapping) != m){stop("fixed_mapping must have length union_nodes.")}
    fixed_mapping <- as.integer(fixed_mapping)
  }
  for(step in seq_len(m)){
    union_node <- order[step]
    fixed_node <- if(is.null(fixed_mapping)){NA_integer_}else{fixed_mapping[union_node]}
    if(!is.na(fixed_node)){
      if(fixed_node < 1 || fixed_node > n || used[fixed_node]){
        return(list(success = FALSE, mapping = integer(0), log_weight = -Inf))
      }
      candidates <- fixed_node
      log_g <- sm.pcand.logp(P, compiled, step, mapping, candidates)
      if(!is.finite(log_g[1])){
        return(list(success = FALSE, mapping = integer(0), log_weight = -Inf))
      }
      chosen <- fixed_node
      log_weight <- log_weight + log_g[1]
    }else{
      candidates <- which(!used)
      log_g <- sm.pcand.logp(P, compiled, step, mapping, candidates)
      log_z <- motif.log.sum(log_g)
      if(!is.finite(log_z)){
        return(list(success = FALSE, mapping = integer(0), log_weight = -Inf))
      }
      probs <- exp(log_g - log_z)
      chosen <- if(length(candidates) == 1){candidates}else{sample(candidates, 1, prob = probs)}
      log_weight <- log_weight + log_z
    }
    mapping[step] <- chosen
    used[chosen] <- TRUE
  }
  full_mapping <- integer(m)
  full_mapping[order] <- mapping
  list(success = TRUE, mapping = full_mapping, log_weight = log_weight)
}

sm.alloc <- function(B, weights, allocation = c("equal", "proportional")){
  allocation <- match.arg(allocation)
  m <- length(weights)
  if(m == 0){return(integer(0))}
  if(B < m){stop("B must be at least the number of valid overlap layers.")}
  if(allocation == "equal"){
    out <- rep(floor(B/m), m)
    remainder <- B - sum(out)
    if(remainder > 0){out[seq_len(remainder)] <- out[seq_len(remainder)] + 1L}
    return(out)
  }
  if(sum(weights) <= 0){return(sm.alloc(B, rep(1, m), allocation = "equal"))}
  raw <- B * weights/sum(weights)
  out <- pmax(1L, floor(raw))
  while(sum(out) < B){
    i <- which.max(raw - out)
    out[i] <- out[i] + 1L
  }
  while(sum(out) > B){
    candidates <- which(out > 1L)
    i <- candidates[which.min(raw[candidates] - out[candidates])]
    out[i] <- out[i] - 1L
  }
  out
}

sm.est.mean <- function(P, H, B = 10000, directed = FALSE,
                        induced = TRUE, aut_count = NULL,
                        normalize_weight = FALSE, seed = NULL,
                        verbose = TRUE, preprocessed = FALSE){
  if(!is.numeric(B) || length(B) != 1 || is.na(B) || B < 1 || B != as.integer(B)){
    stop("B must be a positive integer.")
  }
  B <- as.integer(B)
  if(!is.null(seed)){set.seed(seed)}
  P <- motif.matrix.to.probmat.max(P, directed = directed,
                normalize_weight = normalize_weight,
                preprocessed = preprocessed)
  H <- motif.norm.motif(H, directed = directed)
  n <- nrow(P)
  k <- nrow(H)
  if(k > n){stop("H cannot have more nodes than P.")}
  if(is.null(aut_count)){aut_count <- motif.automorphism.count(H, directed = directed)}
  density <- motif.prob.density(P, directed = directed)
  subgraph_counts <- sm.consts(H, directed = directed, induced = induced)
  edge_counts <- matrix(0L, k, k)
  nonedge_counts <- matrix(0L, k, k)
  for(i in seq_len(subgraph_counts$n)){
    cnt <- sm.add.count(edge_counts, nonedge_counts, subgraph_counts$from[i],
                        subgraph_counts$to[i], subgraph_counts$edge_count[i] > 0L,
                        directed = directed)
    edge_counts <- cnt$edge
    nonedge_counts <- cnt$nonedge
  }
  compiled <- sm.compilen(edge_counts, nonedge_counts,
                          directed = directed, graph_density = density,
                          phi_nodes = seq_len(k), union_nodes = k,
                          compatible = TRUE, target = "joint")
  log_weights <- numeric(B)
  successes <- logical(B)
  for(b in seq_len(B)){
    draw <- sm.sample.compile(P, compiled)
    log_weights[b] <- draw$log_weight
    successes[b] <- draw$success
    if(verbose && B >= 10 && b %% max(1, floor(B/10)) == 0){
      cat("Finished", b, "of", B, "mean ISS samples\n")
    }
  }
  labeled_mean <- motif.mean.logs(log_weights)
  labeled_se <- sqrt(sm.log.var(log_weights)/B)
  list(
    B = B,
    labeled_mean = labeled_mean,
    labeled_se = labeled_se,
    unlabeled_mean = labeled_mean/aut_count,
    unlabeled_se = labeled_se/aut_count,
    aut_count = aut_count,
    success_rate = mean(successes),
    subgraph_order = compiled$order
  )
}

sm.layer.weight <- function(P, compiled, subgraph_constraints, pattern_count){
  draw <- sm.sample.compile(P = P, compiled = compiled)
  if(!draw$success){
    return(list(weight = 0, joint_weight = 0, product_weight = 0,
                success = FALSE))
  }
  if(isTRUE(compiled$compatible)){
    q_joint <- sm.constnp(P, draw$mapping, compiled$constraint_list)
    if(q_joint <= 0 || !is.finite(q_joint)){
      return(list(weight = 0, joint_weight = 0, product_weight = 0,
                  success = FALSE))
    }
    q_product <- sm.constnp(P, draw$mapping,
                            compiled$product_constraint_list)
    joint_weight <- pattern_count * exp(draw$log_weight)
    product_weight <- joint_weight * q_product/q_joint
    weight <- joint_weight - product_weight
  }else{
    joint_weight <- 0
    product_weight <- pattern_count * exp(draw$log_weight)
    weight <- -product_weight
  }
  list(weight = weight, joint_weight = joint_weight,
       product_weight = product_weight, success = TRUE)
}

sm.pairwise.union.weight <- function(P, compiled, pattern_count){
  draw <- sm.sample.compile(P = P, compiled = compiled)
  if(!draw$success){return(list(weight = 0, success = FALSE))}
  list(weight = pattern_count * exp(draw$log_weight), success = TRUE)
}

est.pairunion <- function(X, H,
                          B = 10000,
                          directed = FALSE,
                          induced = TRUE,
                          normalize_weight = FALSE,
                          overlap_sizes = NULL,
                          mode = c("count", "expectation"),
                          allocation = c("equal", "proportional"),
                          mu_labeled = NULL,
                          mu_unlabeled = NULL,
                          aut_count = NULL,
                          include_self_overlap = FALSE,
                          seed = NULL,
                          verbose = TRUE,
                          verbose.B = TRUE,
                          preprocessed = FALSE,
                          ci_level = 0.95,
                          bootstrap_ci = FALSE,
                          bootstrap_B = 1000,
                          bootstrap_seed = NULL,
                          verbose.boot = FALSE){
  mode <- match.arg(mode)
  allocation <- match.arg(allocation)
  if(!is.numeric(B) || length(B) != 1 || is.na(B) || B < 1 || B != as.integer(B)){
    stop("B must be a positive integer.")
  }
  B <- as.integer(B)
  if(!is.null(seed)){set.seed(seed)}
  
  if(mode == "count"){
    P <- motif.graph.to.matrix(X, weighted = FALSE)
    motif.assert.square(P, "X")
    P <- 1L * (P != 0)
    if(!directed){P <- 1L * ((P + t(P)) > 0)}
    diag(P) <- 0L
  }else{
    P <- motif.matrix.to.probmat.max(X, directed = directed,
                  normalize_weight = normalize_weight,
                  preprocessed = preprocessed)
  }
  
  H <- motif.norm.motif(H, directed = directed)
  n <- nrow(P)
  k <- nrow(H)
  H_edge_count <- sm.edge.count(H, directed = directed)
  if(k > n){stop("H cannot have more nodes than X.")}
  if(is.null(aut_count)){aut_count <- motif.automorphism.count(H, directed = directed)}
  if(is.null(overlap_sizes)){
    overlap_sizes <- if(k <= 2L){integer(0)}else{seq.int(2L, k - 1L)}
  }
  overlap_sizes <- sort(unique(as.integer(overlap_sizes)))
  if(any(overlap_sizes < 0 | overlap_sizes > k)){
    stop("overlap_sizes must be between 0 and the number of subgraph nodes.")
  }
  if(!isTRUE(include_self_overlap)){
    overlap_sizes <- overlap_sizes[overlap_sizes < k]
  }
  if(length(overlap_sizes) == 0L){
    stop("No proper pairwise-union overlap size is available. ",
         "For second-order structure, overlap_sizes must exclude the full-overlap case s = k.")
  }
  
  density <- motif.prob.density(P, directed = directed)
  compiled_layers <- vector("list", length(overlap_sizes))
  compatible_counts <- integer(length(overlap_sizes))
  layer_embedding_counts <- numeric(length(overlap_sizes))
  for(idx in seq_along(overlap_sizes)){
    s <- overlap_sizes[idx]
    compiled_layers[[idx]] <- sm.compile.overlap(
      H = H, s = s, directed = directed, induced = induced,
      n = n, graph_density = density,
      include_incompatible = FALSE
    )
    compatible_counts[idx] <- length(compiled_layers[[idx]])
    layer_embedding_counts[idx] <- compatible_counts[idx] * sm.falling(n, 2L * k - s)
  }
  
  valid <- compatible_counts > 0 & layer_embedding_counts > 0
  overlap_sizes <- overlap_sizes[valid]
  compiled_layers <- compiled_layers[valid]
  compatible_counts <- compatible_counts[valid]
  layer_embedding_counts <- layer_embedding_counts[valid]
  if(length(overlap_sizes) == 0){
    stop("No compatible pairwise union pattern is available.")
  }
  
  B_layers <- sm.alloc(B, layer_embedding_counts, allocation = allocation)
  layer_labeled <- numeric(length(overlap_sizes))
  layer_labeled_se <- numeric(length(overlap_sizes))
  layer_unlabeled_union <- numeric(length(overlap_sizes))
  layer_unlabeled_union_se <- numeric(length(overlap_sizes))
  layer_unlabeled_pair <- numeric(length(overlap_sizes))
  layer_unlabeled_pair_se <- numeric(length(overlap_sizes))
  layer_success_rate <- numeric(length(overlap_sizes))
  shared_edge_tables <- vector("list", length(overlap_sizes))
  layer_labeled_weights <- vector("list", length(overlap_sizes))
  layer_unlabeled_union_weights <- vector("list", length(overlap_sizes))
  layer_unlabeled_pair_weights <- vector("list", length(overlap_sizes))
  
  for(idx in seq_along(overlap_sizes)){
    b_s <- B_layers[idx]
    patterns <- compiled_layers[[idx]]
    pattern_count <- compatible_counts[idx]
    union_aut <- vapply(patterns, function(compiled){
      sm.constraint.aut(compiled$edge_counts, compiled$nonedge_counts,
                        directed = directed)
    }, numeric(1))
    shared_edges <- vapply(patterns, sm.shared.edge.count, integer(1))
    shared_edge_levels <- sort(unique(shared_edges))
    
    labeled_w <- numeric(b_s)
    unlabeled_union_w <- numeric(b_s)
    unlabeled_pair_w <- numeric(b_s)
    successes <- logical(b_s)
    shared_edge_labeled <- matrix(0, nrow = b_s, ncol = length(shared_edge_levels))
    shared_edge_unlabeled_pair <- matrix(0, nrow = b_s,
                                         ncol = length(shared_edge_levels))
    
    for(b in seq_len(b_s)){
      j <- sample.int(pattern_count, 1L)
      draw <- sm.pairwise.union.weight(P, patterns[[j]], pattern_count)
      labeled_w[b] <- draw$weight
      unlabeled_union_w[b] <- draw$weight / union_aut[j]
      unlabeled_pair_w[b] <- draw$weight / (aut_count^2)
      shared_idx <- match(shared_edges[j], shared_edge_levels)
      shared_edge_labeled[b, shared_idx] <- draw$weight
      shared_edge_unlabeled_pair[b, shared_idx] <- draw$weight/(aut_count^2)
      successes[b] <- draw$success
      if(verbose && verbose.B && b_s >= 10 && b %% max(1, floor(b_s/10)) == 0){
        cat("Finished", b, "of", b_s, "pairwise-union ISS samples\n")
      }
    }
    
    layer_labeled[idx] <- mean(labeled_w)
    layer_labeled_se[idx] <- if(b_s > 1){stats::sd(labeled_w)/sqrt(b_s)}else{NA_real_}
    layer_unlabeled_union[idx] <- mean(unlabeled_union_w)
    layer_unlabeled_union_se[idx] <- if(b_s > 1){stats::sd(unlabeled_union_w)/sqrt(b_s)}else{NA_real_}
    layer_unlabeled_pair[idx] <- mean(unlabeled_pair_w)
    layer_unlabeled_pair_se[idx] <- if(b_s > 1){stats::sd(unlabeled_pair_w)/sqrt(b_s)}else{NA_real_}
    layer_success_rate[idx] <- mean(successes)
    layer_labeled_weights[[idx]] <- labeled_w
    layer_unlabeled_union_weights[[idx]] <- unlabeled_union_w
    layer_unlabeled_pair_weights[[idx]] <- unlabeled_pair_w
    shared_edge_tables[[idx]] <- data.frame(
      overlap_size = overlap_sizes[idx],
      shared_edge_count = shared_edge_levels,
      labeled_joint_overlap_mass = colMeans(shared_edge_labeled),
      unlabeled_pair_overlap_mass = colMeans(shared_edge_unlabeled_pair),
      row.names = NULL
    )
    
    if(verbose){
      cat("Finished overlap size", overlap_sizes[idx], "with", b_s,
          "compatible pairwise-union ISS samples\n")
    }
  }
  
  labeled_total <- sum(layer_labeled)
  labeled_se <- sqrt(sum(layer_labeled_se^2, na.rm = TRUE))
  unlabeled_union_total <- sum(layer_unlabeled_union)
  unlabeled_union_se <- sqrt(sum(layer_unlabeled_union_se^2, na.rm = TRUE))
  unlabeled_pair_total <- sum(layer_unlabeled_pair)
  unlabeled_pair_se <- sqrt(sum(layer_unlabeled_pair_se^2, na.rm = TRUE))
  
  suffix <- if(mode == "count"){"count"}else{"expectation"}
  layer_table <- data.frame(
    overlap_size = overlap_sizes,
    compatible_pattern_count = compatible_counts,
    compatible_embedding_count = layer_embedding_counts,
    samples = B_layers,
    success_rate = layer_success_rate,
    labeled_value = layer_labeled,
    labeled_se = layer_labeled_se,
    unlabeled_union_value = layer_unlabeled_union,
    unlabeled_union_se = layer_unlabeled_union_se,
    unlabeled_pair_value = layer_unlabeled_pair,
    unlabeled_pair_se = layer_unlabeled_pair_se,
    row.names = NULL
  )
  names(layer_table)[names(layer_table) == "labeled_value"] <- paste0("labeled_", suffix)
  names(layer_table)[names(layer_table) == "unlabeled_union_value"] <- paste0("unlabeled_union_", suffix)
  names(layer_table)[names(layer_table) == "unlabeled_pair_value"] <- paste0("unlabeled_pair_", suffix)
  if(is.null(mu_labeled) && !is.null(mu_unlabeled)){
    mu_labeled <- mu_unlabeled * aut_count
  }else if(is.null(mu_unlabeled) && !is.null(mu_labeled)){
    mu_unlabeled <- mu_labeled/aut_count
  }
  clumping <- sm.pu.clumping(
    layer_table = layer_table,
    suffix = suffix,
    aut_count = aut_count,
    subgraph_nodes = k,
    subgraph_edges = H_edge_count,
    shared_edge_table = do.call(rbind, shared_edge_tables),
    mu_labeled = mu_labeled,
    mu_unlabeled = mu_unlabeled
  )
  inference <- list(
    level = ci_level,
    labeled = sm.inference(
      layer_labeled_weights,
      estimate = labeled_total,
      se = labeled_se,
      level = ci_level,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      verbose.boot = verbose.boot
    ),
    unlabeled_union = sm.inference(
      layer_unlabeled_union_weights,
      estimate = unlabeled_union_total,
      se = unlabeled_union_se,
      level = ci_level,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = if(is.null(bootstrap_seed)) NULL else bootstrap_seed + 1L,
      verbose.boot = verbose.boot
    ),
    unlabeled_pair = sm.inference(
      layer_unlabeled_pair_weights,
      estimate = unlabeled_pair_total,
      se = unlabeled_pair_se,
      level = ci_level,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = if(is.null(bootstrap_seed)) NULL else bootstrap_seed + 2L,
      verbose.boot = verbose.boot
    )
  )
  
  result <- list(
    call = match.call(),
    B = B,
    proposal = "iss_compatible_pairwise_union",
    ordering = "pairwise-union mot.order",
    mode = mode,
    directed = directed,
    induced = induced,
    normalize_weight = normalize_weight,
    include_self_overlap = include_self_overlap,
    preprocessed = preprocessed,
    graph_nodes = n,
    subgraph_nodes = k,
    probability_density = density,
    aut_count = aut_count,
    overlap_sizes = overlap_sizes,
    allocation = allocation,
    total = labeled_total,
    labeled_total = labeled_total,
    labeled_se = labeled_se,
    unlabeled_union_total = unlabeled_union_total,
    unlabeled_union_se = unlabeled_union_se,
    unlabeled_pair_total = unlabeled_pair_total,
    unlabeled_pair_se = unlabeled_pair_se,
    layer_table = layer_table,
    mean = list(labeled = mu_labeled, unlabeled = mu_unlabeled),
    clumping = clumping,
    inference = inference
  )
  class(result) <- c("est.pairunion", "pairwise.union")
  result
}

sm.safe.div <- function(num, den){
  if(length(den) != 1L || !is.finite(den) || den == 0){return(NA_real_)}
  num/den
}

sm.pu.clumping <- function(layer_table, suffix, aut_count,
                           subgraph_nodes = NULL,
                           subgraph_edges = NULL,
                           shared_edge_table = NULL,
                           mu_labeled = NULL,
                           mu_unlabeled = NULL){
  value_col <- paste0("labeled_", suffix)
  union_col <- paste0("unlabeled_union_", suffix)
  pair_col <- paste0("unlabeled_pair_", suffix)
  if(is.null(layer_table[[value_col]])){return(NULL)}
  
  joint <- layer_table[[value_col]]
  s <- layer_table$overlap_size
  joint_total <- sum(joint, na.rm = TRUE)
  joint_share <- if(joint_total > 0){
    joint/joint_total
  }else{
    rep(NA_real_, length(joint))
  }
  
  unlabeled_union_total <- if(!is.null(layer_table[[union_col]])){
    sum(layer_table[[union_col]], na.rm = TRUE)
  }else{
    NA_real_
  }
  unlabeled_pair_total <- if(!is.null(layer_table[[pair_col]])){
    sum(layer_table[[pair_col]], na.rm = TRUE)
  }else{
    joint_total/(aut_count^2)
  }
  node_overlap_mass <- sum(s * joint, na.rm = TRUE)
  unlabeled_node_overlap_mass <- node_overlap_mass/(aut_count^2)
  shared_edge_total <- if(!is.null(shared_edge_table) &&
                          nrow(shared_edge_table) > 0L){
    sum(shared_edge_table$labeled_joint_overlap_mass, na.rm = TRUE)
  }else{
    NA_real_
  }
  edge_overlap_mass <- if(!is.null(shared_edge_table) &&
                          nrow(shared_edge_table) > 0L){
    sum(shared_edge_table$shared_edge_count *
          shared_edge_table$labeled_joint_overlap_mass, na.rm = TRUE)
  }else{
    NA_real_
  }
  unlabeled_edge_overlap_mass <- edge_overlap_mass/(aut_count^2)
  
  by_overlap <- data.frame(
    overlap_size = s,
    labeled_joint_overlap_mass = joint,
    unlabeled_union_overlap_mass = layer_table[[union_col]],
    unlabeled_pair_overlap_mass = layer_table[[pair_col]],
    joint_overlap_share = joint_share,
    row.names = NULL
  )
  by_shared_edge <- NULL
  if(!is.null(shared_edge_table) && nrow(shared_edge_table) > 0L){
    by_shared_edge <- stats::aggregate(
      cbind(labeled_joint_overlap_mass, unlabeled_pair_overlap_mass) ~
        shared_edge_count,
      data = shared_edge_table,
      FUN = sum
    )
    by_shared_edge$edge_overlap_share <- if(shared_edge_total > 0){
      by_shared_edge$labeled_joint_overlap_mass/shared_edge_total
    }else{
      NA_real_
    }
  }
  
  summary <- list(
    labeled_joint_overlap_mass = joint_total,
    unlabeled_union_overlap_mass = unlabeled_union_total,
    unlabeled_pair_overlap_mass = unlabeled_pair_total,
    unlabeled_joint_overlap_mass = unlabeled_pair_total,
    node_overlap_mass = node_overlap_mass,
    unlabeled_node_overlap_mass = unlabeled_node_overlap_mass,
    edge_overlap_mass = edge_overlap_mass,
    unlabeled_edge_overlap_mass = unlabeled_edge_overlap_mass,
    weighted_mean_overlap_size = sm.safe.div(node_overlap_mass, joint_total),
    overlap_tau = if(joint_total > 0){sum(joint_share^2, na.rm = TRUE)}else{NA_real_},
    max_overlap_share = if(joint_total > 0){max(joint_share, na.rm = TRUE)}else{NA_real_},
    node_overlap_intensity = if(!is.null(mu_labeled) && is.finite(mu_labeled) &&
                                !is.null(subgraph_nodes) && subgraph_nodes > 1){
      sm.safe.div(node_overlap_mass, (subgraph_nodes - 1) * (mu_labeled^2))
    }else{NA_real_},
    edge_overlap_intensity = if(!is.null(mu_labeled) && is.finite(mu_labeled) &&
                                !is.null(subgraph_edges) && subgraph_edges > 0){
      sm.safe.div(edge_overlap_mass, subgraph_edges * (mu_labeled^2))
    }else{NA_real_},
    unlabeled_node_overlap_intensity = if(!is.null(mu_unlabeled) &&
                                          is.finite(mu_unlabeled) &&
                                          !is.null(subgraph_nodes) &&
                                          subgraph_nodes > 1){
      sm.safe.div(unlabeled_node_overlap_mass,
                  (subgraph_nodes - 1) * (mu_unlabeled^2))
    }else{NA_real_},
    unlabeled_edge_overlap_intensity = if(!is.null(mu_unlabeled) &&
                                          is.finite(mu_unlabeled) &&
                                          !is.null(subgraph_edges) &&
                                          subgraph_edges > 0){
      sm.safe.div(unlabeled_edge_overlap_mass,
                  subgraph_edges * (mu_unlabeled^2))
    }else{NA_real_},
    effective_subgraph_count = if(!is.null(mu_labeled) && is.finite(mu_labeled)){
      sm.safe.div(mu_labeled^2, joint_total)
    }else{NA_real_},
    effective_unlabeled_subgraph_count = if(!is.null(mu_unlabeled) && is.finite(mu_unlabeled)){
      sm.safe.div(mu_unlabeled^2, unlabeled_pair_total)
    }else{NA_real_}
  )
  
  list(summary = summary, by_overlap = by_overlap, by_shared_edge = by_shared_edge)
}

sm.clumping <- function(layer_table, mu_labeled, mu_unlabeled,
                        labeled_variance_raw, unlabeled_variance_raw,
                        aut_count){
  if(is.null(layer_table$labeled_joint_overlap_mass)){
    return(NULL)
  }
  joint <- layer_table$labeled_joint_overlap_mass
  product <- layer_table$labeled_marginal_product_mass
  cov <- layer_table$labeled_variance_contribution
  s <- layer_table$overlap_size
  
  joint_total <- sum(joint, na.rm = TRUE)
  product_total <- sum(product, na.rm = TRUE)
  abs_cov_total <- sum(abs(cov), na.rm = TRUE)
  second_labeled <- mu_labeled^2 + labeled_variance_raw
  second_unlabeled <- mu_unlabeled^2 + unlabeled_variance_raw
  unlabeled_joint_total <- joint_total/(aut_count^2)
  unlabeled_product_total <- product_total/(aut_count^2)
  
  joint_share <- if(joint_total > 0){joint/joint_total}else{rep(NA_real_, length(joint))}
  abs_cov_share <- if(abs_cov_total > 0){abs(cov)/abs_cov_total}else{rep(NA_real_, length(cov))}
  
  by_overlap <- data.frame(
    overlap_size = s,
    labeled_joint_overlap_mass = joint,
    labeled_marginal_product_mass = product,
    labeled_excess_overlap_mass = cov,
    joint_overlap_share = joint_share,
    abs_covariance_share = abs_cov_share,
    row.names = NULL
  )
  
  summary <- list(
    labeled_joint_overlap_mass = joint_total,
    unlabeled_joint_overlap_mass = unlabeled_joint_total,
    labeled_marginal_product_mass = product_total,
    unlabeled_marginal_product_mass = unlabeled_product_total,
    labeled_excess_overlap_mass = labeled_variance_raw,
    unlabeled_excess_overlap_mass = unlabeled_variance_raw,
    joint_to_product_overlap_ratio = sm.safe.div(joint_total, product_total),
    overlap_fraction_of_second_moment = sm.safe.div(joint_total, second_labeled),
    effective_subgraph_count = sm.safe.div(mu_labeled^2, joint_total),
    weighted_mean_overlap_size = sm.safe.div(sum(s * joint, na.rm = TRUE), joint_total),
    overlap_tau = if(joint_total > 0){sum(joint_share^2, na.rm = TRUE)}else{NA_real_},
    max_overlap_share = if(joint_total > 0){max(joint_share, na.rm = TRUE)}else{NA_real_},
    abs_covariance_tau = if(abs_cov_total > 0){sum(abs_cov_share^2, na.rm = TRUE)}else{NA_real_},
    variance_to_mean_ratio = sm.safe.div(unlabeled_variance_raw, mu_unlabeled),
    overdispersion_c = if(is.finite(mu_unlabeled) && mu_unlabeled > 0){
      (unlabeled_variance_raw - mu_unlabeled)/(mu_unlabeled^2)
    }else{NA_real_},
    labeled_second_moment = second_labeled,
    unlabeled_second_moment = second_unlabeled
  )
  
  list(summary = summary, by_overlap = by_overlap)
}

est.variance <- function(P, H,
                         B = 10000,
                         B_mean = B,
                         directed = FALSE,
                         induced = TRUE,
                         normalize_weight = FALSE,
                         overlap_sizes = NULL,
                         allocation = c("equal", "proportional"),
                         mu_labeled = NULL,
                         mu_unlabeled = NULL,
                         aut_count = NULL,
                         seed = NULL,
                         verbose = TRUE,
                         verbose.B = TRUE,
                         include_incompatible = TRUE,
                         preprocessed = FALSE,
                         ci_level = 0.95,
                         bootstrap_ci = FALSE,
                         bootstrap_B = 1000,
                         bootstrap_seed = NULL,
                         verbose.boot = FALSE){
  allocation <- match.arg(allocation)
  if(!is.numeric(B) || length(B) != 1 || is.na(B) || B < 1 || B != as.integer(B)){
    stop("B must be a positive integer.")
  }
  B <- as.integer(B)
  if(!is.null(seed)){set.seed(seed)}
  P <- motif.matrix.to.probmat.max(P, directed = directed,
                normalize_weight = normalize_weight,
                preprocessed = preprocessed)
  H <- motif.norm.motif(H, directed = directed)
  n <- nrow(P)
  k <- nrow(H)
  if(k > n){stop("H cannot have more nodes than P.")}
  if(is.null(aut_count)){aut_count <- motif.automorphism.count(H, directed = directed)}
  if(is.null(overlap_sizes)){overlap_sizes <- if(k <= 1){k}else{seq.int(2L, k)}}
  overlap_sizes <- sort(unique(as.integer(overlap_sizes)))
  if(any(overlap_sizes < 0 | overlap_sizes > k)){
    stop("overlap_sizes must be between 0 and the number of subgraph nodes.")
  }
  density <- motif.prob.density(P, directed = directed)
  layer_pair_counts <- vapply(overlap_sizes, function(s){sm.layer(n, k, s)}, numeric(1))
  valid <- layer_pair_counts > 0
  overlap_sizes <- overlap_sizes[valid]
  layer_pair_counts <- layer_pair_counts[valid]
  if(length(overlap_sizes) == 0){stop("No valid overlap layer is available.")}
  compiled_layers <- vector("list", length(overlap_sizes))
  compatible_counts <- integer(length(overlap_sizes))
  incompatible_counts <- integer(length(overlap_sizes))
  pattern_counts <- integer(length(overlap_sizes))
  for(idx in seq_along(overlap_sizes)){
    compiled_layers[[idx]] <- sm.compile.overlap(
      H = H, s = overlap_sizes[idx], directed = directed, induced = induced,
      n = n, graph_density = density,
      include_incompatible = include_incompatible
    )
    pattern_counts[idx] <- length(compiled_layers[[idx]])
    compatible_counts[idx] <- sum(vapply(compiled_layers[[idx]], function(x){isTRUE(x$compatible)}, logical(1)))
    incompatible_counts[idx] <- pattern_counts[idx] - compatible_counts[idx]
    if(pattern_counts[idx] == 0){
      stop("No valid overlap pattern is available for overlap size ", overlap_sizes[idx], ".")
    }
  }
  B_layers <- sm.alloc(B, layer_pair_counts, allocation = allocation)
  subgraph_constraints <- sm.consts(H, directed = directed, induced = induced)
  layer_covariances <- numeric(length(overlap_sizes))
  layer_ses <- numeric(length(overlap_sizes))
  layer_mean_weight <- numeric(length(overlap_sizes))
  layer_success_rate <- numeric(length(overlap_sizes))
  layer_joint_overlap_mass <- numeric(length(overlap_sizes))
  layer_joint_overlap_ses <- numeric(length(overlap_sizes))
  layer_product_mass <- numeric(length(overlap_sizes))
  layer_product_ses <- numeric(length(overlap_sizes))
  layer_labeled_variance_weights <- vector("list", length(overlap_sizes))
  layer_unlabeled_variance_weights <- vector("list", length(overlap_sizes))
  layer_joint_weights <- vector("list", length(overlap_sizes))
  layer_product_weights <- vector("list", length(overlap_sizes))
  for(idx in seq_along(overlap_sizes)){
    b_s <- B_layers[idx]
    patterns <- compiled_layers[[idx]]
    pattern_count <- pattern_counts[idx]
    cov_weights <- numeric(b_s)
    joint_weights <- numeric(b_s)
    product_weights <- numeric(b_s)
    successes <- logical(b_s)
    for(b in seq_len(b_s)){
      compiled <- patterns[[sample.int(pattern_count, 1)]]
      draw <- sm.layer.weight(P, compiled, subgraph_constraints, pattern_count)
      cov_weights[b] <- draw$weight
      joint_weights[b] <- draw$joint_weight
      product_weights[b] <- draw$product_weight
      successes[b] <- draw$success
      if(verbose && verbose.B && b_s >= 10 && b %% max(1, floor(b_s/10)) == 0){
        cat("Finished", b, "of", b_s, "union-ISS covariance samples\n")
      }
    }
    layer_mean_weight[idx] <- mean(cov_weights)
    layer_covariances[idx] <- layer_mean_weight[idx]
    layer_ses[idx] <- if(b_s > 1){stats::sd(cov_weights)/sqrt(b_s)}else{NA_real_}
    layer_joint_overlap_mass[idx] <- mean(joint_weights)
    layer_joint_overlap_ses[idx] <- if(b_s > 1){stats::sd(joint_weights)/sqrt(b_s)}else{NA_real_}
    layer_product_mass[idx] <- mean(product_weights)
    layer_product_ses[idx] <- if(b_s > 1){stats::sd(product_weights)/sqrt(b_s)}else{NA_real_}
    layer_success_rate[idx] <- mean(successes)
    layer_labeled_variance_weights[[idx]] <- cov_weights
    layer_unlabeled_variance_weights[[idx]] <- cov_weights/(aut_count^2)
    layer_joint_weights[[idx]] <- joint_weights
    layer_product_weights[[idx]] <- product_weights
    if(verbose){
      cat("Finished overlap size", overlap_sizes[idx], "with", b_s,
          "union-ISS covariance samples\n")
    }
  }
  labeled_variance_raw <- sum(layer_covariances)
  labeled_variance_se <- sqrt(sum(layer_ses^2, na.rm = TRUE))
  unlabeled_variance_raw <- labeled_variance_raw/(aut_count^2)
  unlabeled_variance_se <- labeled_variance_se/(aut_count^2)
  labeled_variance <- max(labeled_variance_raw, 0)
  unlabeled_variance <- max(unlabeled_variance_raw, 0)
  mean_result <- NULL
  if(is.null(mu_labeled) && is.null(mu_unlabeled)){
    mean_result <- sm.est.mean(
      P = P, H = H, B = B_mean, directed = directed, induced = induced,
      aut_count = aut_count, normalize_weight = FALSE, seed = NULL,
      verbose = verbose, preprocessed = TRUE
    )
    mu_labeled <- mean_result$labeled_mean
    mu_unlabeled <- mean_result$unlabeled_mean
  }else if(is.null(mu_labeled)){
    mu_labeled <- mu_unlabeled * aut_count
  }else if(is.null(mu_unlabeled)){
    mu_unlabeled <- mu_labeled/aut_count
  }
  c_raw <- if(mu_unlabeled > 0){(unlabeled_variance_raw - mu_unlabeled)/(mu_unlabeled^2)}else{NA_real_}
  c_overdispersed <- if(is.finite(c_raw)){max(c_raw, 0)}else{NA_real_}
  layer_table <- data.frame(
    overlap_size = overlap_sizes,
    pair_count = layer_pair_counts,
    pattern_count = pattern_counts,
    compatible_pattern_count = compatible_counts,
    incompatible_pattern_count = incompatible_counts,
    samples = B_layers,
    success_rate = layer_success_rate,
    mean_importance_weight = layer_mean_weight,
    labeled_joint_overlap_mass = layer_joint_overlap_mass,
    labeled_joint_overlap_mass_se = layer_joint_overlap_ses,
    labeled_marginal_product_mass = layer_product_mass,
    labeled_marginal_product_mass_se = layer_product_ses,
    labeled_variance_contribution = layer_covariances,
    labeled_variance_contribution_se = layer_ses,
    row.names = NULL
  )
  clumping <- sm.clumping(
    layer_table = layer_table,
    mu_labeled = mu_labeled,
    mu_unlabeled = mu_unlabeled,
    labeled_variance_raw = labeled_variance_raw,
    unlabeled_variance_raw = unlabeled_variance_raw,
    aut_count = aut_count
  )
  inference <- list(
    level = ci_level,
    labeled_raw = sm.inference(
      layer_labeled_variance_weights,
      estimate = labeled_variance_raw,
      se = labeled_variance_se,
      level = ci_level,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      verbose.boot = verbose.boot,
      use_abs_diagnostics = TRUE
    ),
    unlabeled_raw = sm.inference(
      layer_unlabeled_variance_weights,
      estimate = unlabeled_variance_raw,
      se = unlabeled_variance_se,
      level = ci_level,
      bootstrap_ci = bootstrap_ci,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = if(is.null(bootstrap_seed)) NULL else bootstrap_seed + 1L,
      verbose.boot = verbose.boot,
      use_abs_diagnostics = TRUE
    ),
    labeled_wald = list(
      estimate = labeled_variance,
      se = labeled_variance_se,
      level = ci_level,
      wald_ci = motif.wald.ci(labeled_variance, labeled_variance_se, level = ci_level),
      log_wald_ci = motif.log.wald.ci(labeled_variance, labeled_variance_se, level = ci_level)
    ),
    unlabeled_wald = list(
      estimate = unlabeled_variance,
      se = unlabeled_variance_se,
      level = ci_level,
      wald_ci = motif.wald.ci(unlabeled_variance, unlabeled_variance_se, level = ci_level),
      log_wald_ci = motif.log.wald.ci(unlabeled_variance, unlabeled_variance_se, level = ci_level)
    ),
    diagnostics_note = paste(
      "Variance-contribution weights may be negative;",
      "diagnostics for raw variance estimators are computed from absolute weights."
    ),
    joint_overlap_diagnostics = sm.strata.diagns(layer_joint_weights, use_abs = FALSE),
    marginal_product_diagnostics = sm.strata.diagns(layer_product_weights, use_abs = FALSE)
  )
  result <- list(
    call = match.call(),
    B = B,
    B_mean = B_mean,
    proposal = "iss_pairwise_union",
    ordering = "pairwise-union mot.order",
    directed = directed,
    induced = induced,
    normalize_weight = normalize_weight,
    preprocessed = preprocessed,
    include_incompatible = include_incompatible,
    probability_rule = if(normalize_weight){"P_ij = W_ij / max(W)"}else{"P_ij = W_ij"},
    graph_nodes = n,
    subgraph_nodes = k,
    probability_density = density,
    aut_count = aut_count,
    overlap_sizes = overlap_sizes,
    allocation = allocation,
    layer_table = layer_table,
    mean = list(labeled = mu_labeled, unlabeled = mu_unlabeled,
                estimated = !is.null(mean_result), mean_result = mean_result),
    variance = list(
      labeled_raw = labeled_variance_raw,
      labeled = labeled_variance,
      labeled_se = labeled_variance_se,
      unlabeled_raw = unlabeled_variance_raw,
      unlabeled = unlabeled_variance,
      unlabeled_se = unlabeled_variance_se
    ),
    overdispersion = list(
      c_raw = c_raw,
      c_overdispersed = c_overdispersed,
      variance_poisson_overdispersed = if(mu_unlabeled > 0){
        mu_unlabeled + c_overdispersed * mu_unlabeled^2
      }else{NA_real_}
    ),
    clumping = clumping,
    inference = inference
  )
  class(result) <- c("variance", "subgraph.var")
  result
}

sec.moment <- function(...){
  res <- est.variance(...)
  res$second_moment <- list(
    labeled = res$variance$labeled_raw + res$mean$labeled^2,
    labeled_se = res$variance$labeled_se,
    unlabeled = res$variance$unlabeled_raw + res$mean$unlabeled^2,
    unlabeled_se = res$variance$unlabeled_se
  )
  class(res) <- c("sec.moment", class(res))
  res
}

print.sec.moment <- function(x, ...){
  cat("Sampled subgraph second-moment estimate\n")
  cat("  graph nodes:        ", x$graph_nodes, "\n", sep = "")
  cat("  subgraph nodes:        ", x$subgraph_nodes, "\n", sep = "")
  cat("  directed:           ", x$directed, "\n", sep = "")
  cat("  induced:            ", x$induced, "\n", sep = "")
  cat("  pair samples:       ", x$B, "\n", sep = "")
  cat("  proposal:           ", x$proposal, "\n", sep = "")
  cat("\nMean\n")
  cat("  unlabeled mu:       ", signif(x$mean$unlabeled, 6), "\n", sep = "")
  cat("  labeled mu:         ", signif(x$mean$labeled, 6), "\n", sep = "")
  cat("\nSecond moment\n")
  cat("  unlabeled E[N^2]:   ", signif(x$second_moment$unlabeled, 6),
      " (SE ", signif(x$second_moment$unlabeled_se, 4), ")\n", sep = "")
  cat("  labeled E[N^2]:     ", signif(x$second_moment$labeled, 6),
      " (SE ", signif(x$second_moment$labeled_se, 4), ")\n", sep = "")
  invisible(x)
}

print.est.pairunion <- function(x, ...){
  label <- if(x$mode == "count"){"count"}else{"expectation"}
  cat("ISS compatible pairwise-union estimate\n")
  cat("  graph nodes:         ", x$graph_nodes, "\n", sep = "")
  cat("  subgraph nodes:         ", x$subgraph_nodes, "\n", sep = "")
  cat("  directed:            ", x$directed, "\n", sep = "")
  cat("  induced:             ", x$induced, "\n", sep = "")
  cat("  mode:                ", x$mode, "\n", sep = "")
  cat("  pair samples:        ", x$B, "\n", sep = "")
  cat("  proposal:            ", x$proposal, "\n", sep = "")
  cat("  ordering:            ", x$ordering, "\n", sep = "")
  cat("  overlap sizes:       ", paste(x$overlap_sizes, collapse = ", "), "\n", sep = "")
  cat("  overlap allocation:  ", x$allocation, "\n", sep = "")
  cat("  automorphisms:       ", x$aut_count, "\n", sep = "")
  cat("\nEstimate\n")
  cat("  labeled ", label, ":          ", signif(x$labeled_total, 6),
      " (SE ", signif(x$labeled_se, 4), ")\n", sep = "")
  cat("  unlabeled union ", label, ":  ", signif(x$unlabeled_union_total, 6),
      " (SE ", signif(x$unlabeled_union_se, 4), ")\n", sep = "")
  cat("  unlabeled pair ", label, ":   ", signif(x$unlabeled_pair_total, 6),
      " (SE ", signif(x$unlabeled_pair_se, 4), ")\n", sep = "")
  if(!is.null(x$clumping) && !is.null(x$clumping$summary)){
    cs <- x$clumping$summary
    cat("\nCompatible clumping / overlap structure\n")
    cat("  pair overlap mass:  ", signif(cs$unlabeled_pair_overlap_mass, 6), "\n", sep = "")
    cat("  union overlap mass: ", signif(cs$unlabeled_union_overlap_mass, 6), "\n", sep = "")
    cat("  mean overlap size:  ", signif(cs$weighted_mean_overlap_size, 6), "\n", sep = "")
    cat("  overlap concentration:        ", signif(cs$overlap_tau, 6), "\n", sep = "")
    cat("  max overlap share:  ", signif(cs$max_overlap_share, 6), "\n", sep = "")
    cat("  node intensity:     ", signif(cs$unlabeled_node_overlap_intensity, 6), "\n", sep = "")
    cat("  edge intensity:     ", signif(cs$unlabeled_edge_overlap_intensity, 6), "\n", sep = "")
    if(is.finite(cs$effective_unlabeled_subgraph_count)){
      cat("  effective subgraphs:   ", signif(cs$effective_unlabeled_subgraph_count, 6), "\n", sep = "")
    }
  }
  invisible(x)
}

print.subgraph.var <- function(x, ...){
  cat("Sampled subgraph variance estimate\n")
  cat("  graph nodes:        ", x$graph_nodes, "\n", sep = "")
  cat("  subgraph nodes:        ", x$subgraph_nodes, "\n", sep = "")
  cat("  directed:           ", x$directed, "\n", sep = "")
  cat("  induced:            ", x$induced, "\n", sep = "")
  cat("  normalize weight:   ", x$normalize_weight %||% FALSE, "\n", sep = "")
  cat("  pair samples:       ", x$B, "\n", sep = "")
  cat("  proposal:           ", x$proposal %||% "iss_pairwise_union", "\n", sep = "")
  cat("  ordering:           ", x$ordering %||% "pairwise-union mot.order", "\n", sep = "")
  cat("  incompatible used:  ", x$include_incompatible %||% TRUE, "\n", sep = "")
  cat("  overlap sizes:      ", paste(x$overlap_sizes, collapse = ", "), "\n", sep = "")
  cat("  overlap allocation: ", x$allocation, "\n", sep = "")
  cat("  automorphisms:      ", x$aut_count, "\n", sep = "")
  cat("\nVariance\n")
  cat("  unlabeled raw Var:  ", signif(x$variance$unlabeled_raw, 6),
      " (SE ", signif(x$variance$unlabeled_se, 4), ")\n", sep = "")
  cat("  unlabeled Var:      ", signif(x$variance$unlabeled, 6), "\n", sep = "")
  cat("  labeled raw Var:    ", signif(x$variance$labeled_raw, 6),
      " (SE ", signif(x$variance$labeled_se, 4), ")\n", sep = "")
  cat("\nOverdispersion\n")
  cat("  c raw:              ", signif(x$overdispersion$c_raw, 6), "\n", sep = "")
  cat("  c overdispersed:    ", signif(x$overdispersion$c_overdispersed, 6), "\n", sep = "")
  if(!is.null(x$clumping) && !is.null(x$clumping$summary)){
    cs <- x$clumping$summary
    cat("\nClumping / overlap structure\n")
    cat("  joint overlap mass: ", signif(cs$unlabeled_joint_overlap_mass, 6), "\n", sep = "")
    cat("  joint/product ratio:", signif(cs$joint_to_product_overlap_ratio, 6), "\n", sep = "")
    cat("  overlap fraction:   ", signif(cs$overlap_fraction_of_second_moment, 6), "\n", sep = "")
    cat("  effective subgraphs:   ", signif(cs$effective_subgraph_count, 6), "\n", sep = "")
    cat("  mean overlap size:  ", signif(cs$weighted_mean_overlap_size, 6), "\n", sep = "")
    cat("  overlap concentration:        ", signif(cs$overlap_tau, 6), "\n", sep = "")
  }
  invisible(x)
}
