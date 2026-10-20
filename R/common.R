#==========================================================================
## Basic helpers shared by ISS, VISS, and TSGE
#==========================================================================

`%||%` <- function(x, y){if(is.null(x)) y else x}

motif.graph.to.matrix <- function(x, weighted = FALSE){
  if(inherits(x, "tbl_graph")){
    if(!requireNamespace("tidygraph", quietly = TRUE)){
      stop("Package 'tidygraph' is required to use tbl_graph inputs.")
    }
    x <- tidygraph::as.igraph(x)
  }
  if(inherits(x, "igraph")){
    if(!requireNamespace("igraph", quietly = TRUE)){
      stop("Package 'igraph' is required to use igraph or tidygraph inputs.")
    }
    attr_name <- if(weighted && "weight" %in% igraph::edge_attr_names(x)){
      "weight"
    }else{
      NULL
    }
    x <- igraph::as_adjacency_matrix(x, attr = attr_name, sparse = FALSE)
  }
  x <- as.matrix(x)
  if(!weighted){x <- x != 0}
  x
}

motif.assert.square <- function(A, name = "matrix"){
  if(is.null(dim(A)) || length(dim(A)) != 2 || nrow(A) != ncol(A)){
    stop(name, " must be a square matrix.")
  }
  invisible(TRUE)
}

motif.network.type <- function(x,
                               network_type = c("auto", "binary", "weighted"),
                               tolerance = sqrt(.Machine$double.eps)){
  network_type <- match.arg(network_type)
  if(network_type != "auto"){return(network_type)}
  M <- motif.graph.to.matrix(x, weighted = TRUE)
  if(is.null(dim(M)) || length(dim(M)) != 2 || nrow(M) != ncol(M)){
    stop("The observed network must be a square matrix or graph object.")
  }
  if(!is.numeric(M) && !is.logical(M)){
    stop("The observed network must contain numeric or logical edge values.")
  }
  storage.mode(M) <- "numeric"
  M[is.na(M)] <- 0
  if(any(!is.finite(M))){stop("The observed network must contain only finite edge values.")}
  values <- M[row(M) != col(M)]
  is_binary <- all(abs(values) <= tolerance | abs(values - 1) <= tolerance)
  if(is_binary){"binary"}else{"weighted"}
}

motif.est.value <- function(x){
  value <- x$unlabeled_estimate %||% x$unlabeled_expectation %||% x$unlabeled_mean
  if(is.null(value)){stop("The subgraph estimator returned no unlabeled estimate.")}
  value
}

motif.est.se <- function(x){
  value <- x$unlabeled_se
  if(is.null(value)){stop("The subgraph estimator returned no unlabeled standard error.")}
  value
}

motif.est.order <- function(x){
  x$subgraph_order %||% x$motif_order %||% x$order
}

motif.norm.graph <- function(A, directed = FALSE){
  if(isTRUE(attr(A, "norm.graph.directed") == directed)){
    return(A)
  }
  A <- motif.graph.to.matrix(A, weighted = FALSE)
  if(!is.logical(A)){A <- A != 0}
  diag(A) <- FALSE
  if(!directed){
    A <- A | t(A)
    diag(A) <- FALSE
  }
  attr(A, "norm.graph.directed") <- directed
  A
}

motif.norm.motif <- function(H, directed = FALSE){
  H <- motif.graph.to.matrix(H, weighted = FALSE)
  H <- H != 0
  diag(H) <- FALSE
  if(!directed){
    H <- H | t(H)
    diag(H) <- FALSE
  }
  H
}

motif.w.to.probmat <- function(W, directed = FALSE){
  W <- motif.graph.to.matrix(W, weighted = TRUE)
  if(nrow(W) != ncol(W)){stop("W must be a square weighted matrix.")}
  if(any(W < 0, na.rm = TRUE)){stop("W must be non-negative.")}
  W[is.na(W)] <- 0
  diag(W) <- 0
  if(!directed){
    W <- (W + t(W))/2
    diag(W) <- 0
  }
  max_w <- max(W)
  if(!is.finite(max_w) || max_w <= 0){
    stop("W must contain at least one positive off-diagonal weight.")
  }
  P <- W/max_w
  diag(P) <- 0
  P
}

motif.norm.probmat <- function(P, directed = FALSE, cap = 1,
                               strict = FALSE, label = "P"){
  P <- motif.graph.to.matrix(P, weighted = TRUE)
  if(nrow(P) != ncol(P)){stop(label, " must be a square probability matrix.")}
  if(strict && any(P < 0 | P > 1, na.rm = TRUE)){
    stop(label, " must have entries in [0, 1].")
  }
  P[is.na(P)] <- 0
  P <- pmin(pmax(P, 0), cap)
  diag(P) <- 0
  if(!directed){
    P <- (P + t(P))/2
    P <- pmin(pmax(P, 0), cap)
    diag(P) <- 0
  }
  P
}

motif.matrix.to.probmat <- function(W, directed = FALSE, normalize_weight = FALSE,
                                    strict_probability = TRUE){
  if(normalize_weight){
    motif.w.to.probmat(W, directed = directed)
  }else{
    motif.norm.probmat(
      W,
      directed = directed,
      strict = strict_probability,
      label = "W"
    )
  }
}

motif.norm.probmat.max <- function(P, directed = FALSE, strict = TRUE,
                                   preprocessed = FALSE){
  if(isTRUE(preprocessed)){
    motif.assert.square(P, "P")
    return(P)
  }
  rn <- rownames(P)
  cn <- colnames(P)
  P <- motif.graph.to.matrix(P, weighted = TRUE)
  motif.assert.square(P, "P")
  storage.mode(P) <- "numeric"
  if(any(!is.finite(P))){stop("P must contain only finite values.")}
  if(strict && (min(P) < 0 || max(P) > 1)){
    stop("P must have entries in [0, 1]. Use normalize_weight = TRUE for raw weights.")
  }
  P[P < 0] <- 0
  P[P > 1] <- 1
  if(!directed && !isTRUE(all(P == t(P)))){P <- pmax(P, t(P))}
  if(any(diag(P) != 0)){diag(P) <- 0}
  rownames(P) <- rn %||% rownames(P)
  colnames(P) <- cn %||% colnames(P)
  P
}

motif.w.to.probmat.max <- function(W, directed = FALSE,
                                   preprocessed = FALSE){
  if(isTRUE(preprocessed)){
    return(motif.norm.probmat.max(
      W, directed = directed, strict = TRUE, preprocessed = TRUE
    ))
  }
  rn <- rownames(W)
  cn <- colnames(W)
  W <- motif.graph.to.matrix(W, weighted = TRUE)
  motif.assert.square(W, "W")
  storage.mode(W) <- "numeric"
  if(any(!is.finite(W))){stop("W must contain only finite values.")}
  if(min(W) < 0){stop("W must be nonnegative.")}
  if(!directed && !isTRUE(all(W == t(W)))){W <- pmax(W, t(W))}
  if(any(diag(W) != 0)){diag(W) <- 0}
  scale <- max(W)
  P <- if(is.finite(scale) && scale > 0){W/scale}else{W}
  rownames(P) <- rn %||% rownames(P)
  colnames(P) <- cn %||% colnames(P)
  motif.norm.probmat.max(P, directed = directed, strict = TRUE)
}

motif.matrix.to.probmat.max <- function(W, directed = FALSE,
                                        normalize_weight = FALSE,
                                        preprocessed = FALSE){
  if(isTRUE(normalize_weight)){
    motif.w.to.probmat.max(W, directed = directed, preprocessed = FALSE)
  }else{
    motif.norm.probmat.max(
      W, directed = directed, strict = TRUE, preprocessed = preprocessed
    )
  }
}

motif.graph.density <- function(G, directed = FALSE){
  G <- motif.norm.graph(G, directed = directed)
  n <- nrow(G)
  if(n <= 1){return(0)}
  if(directed){sum(G)/(n * (n - 1))}else{sum(G)/(n * (n - 1))}
}

motif.prob.density <- function(P, directed = FALSE){
  n <- nrow(P)
  if(n <= 1){return(0)}
  if(directed){sum(P)/(n * (n - 1))}else{sum(P[upper.tri(P)])/choose(n, 2)}
}

motif.log.sum <- function(x){
  finite <- is.finite(x)
  if(!any(finite)){return(-Inf)}
  m <- max(x[finite])
  m + log(sum(exp(x[finite] - m)))
}

motif.mean.logs <- function(log_w){
  exp(motif.log.sum(log_w) - log(length(log_w)))
}

motif.var.logs <- function(log_w){
  n <- length(log_w)
  if(n <= 1){return(NA_real_)}
  log_sum_w <- motif.log.sum(log_w)
  log_sum_w2 <- motif.log.sum(2 * log_w)
  mean_w <- exp(log_sum_w - log(n))
  sum_w2 <- exp(log_sum_w2)
  max((sum_w2 - n * mean_w^2)/(n - 1), 0)
}

motif.order <- function(H, graph_density, directed = FALSE, eps = 1e-12){
  k <- nrow(H)
  p <- min(max(graph_density, eps), 1 - eps)
  edge_w <- -log(p)
  nonedge_w <- -log(1 - p)
  pair_weight <- function(has_edge){if(has_edge) edge_w else nonedge_w}
  global_score <- numeric(k)
  for(u in seq_len(k)){
    others <- setdiff(seq_len(k), u)
    global_score[u] <- sum(vapply(others, function(v){
      if(directed){
        pair_weight(H[u, v]) + pair_weight(H[v, u])
      }else{
        pair_weight(H[u, v])
      }
    }, numeric(1)))
  }
  first <- which.max(global_score)
  order <- first
  remaining <- setdiff(seq_len(k), first)
  while(length(remaining) > 0){
    scores <- vapply(remaining, function(u){
      imposed <- sum(vapply(order, function(s){
        if(directed){
          pair_weight(H[u, s]) + pair_weight(H[s, u])
        }else{
          pair_weight(H[u, s])
        }
      }, numeric(1)))
      imposed + 1e-9 * global_score[u]
    }, numeric(1))
    next_node <- remaining[which.max(scores)]
    order <- c(order, next_node)
    remaining <- setdiff(remaining, next_node)
  }
  order
}

motif.automorphism.count <- function(H, directed = FALSE){
  H <- motif.norm.motif(H, directed = directed)
  count <- 0L
  visit <- function(prefix, remaining){
    if(length(remaining) == 0){
      if(identical(H[prefix, prefix, drop = FALSE], H)){
        count <<- count + 1L
      }
      return(invisible(NULL))
    }
    for(v in remaining){visit(c(prefix, v), setdiff(remaining, v))}
    invisible(NULL)
  }
  visit(integer(0), seq_len(nrow(H)))
  count
}

motif.weight.diagnostics <- function(weights){
  finite <- is.finite(weights)
  w <- weights[finite]
  n <- length(weights)
  finite_n <- length(w)
  zero_n <- sum(w == 0)
  total <- sum(w)
  sum_w2 <- sum(w^2)
  mean_w <- if(finite_n > 0){mean(w)}else{NA_real_}
  sd_w <- if(finite_n > 1){stats::sd(w)}else{NA_real_}
  ess <- if(is.finite(total) && total > 0 && is.finite(sum_w2) && sum_w2 > 0){
    total^2/sum_w2
  }else{0}
  list(
    n = n,
    finite = finite_n,
    zero_count = zero_n,
    zero_rate = if(finite_n > 0){zero_n/finite_n}else{NA_real_},
    positive_rate = if(finite_n > 0){mean(w > 0)}else{NA_real_},
    ess = ess,
    ess_ratio = if(finite_n > 0){ess/finite_n}else{NA_real_},
    cv_weight = if(is.finite(mean_w) && mean_w > 0 && is.finite(sd_w)){
      sd_w/mean_w
    }else{NA_real_},
    max_weight = if(finite_n > 0){max(w)}else{NA_real_},
    max_weight_share = if(is.finite(total) && total > 0){max(w)/total}else{NA_real_}
  )
}

motif.wald.ci <- function(estimate, se, level = 0.95){
  if(!is.finite(estimate) || !is.finite(se) || se < 0){
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level)/2)
  estimate + c(-1, 1) * z * se
}

motif.log.wald.ci <- function(estimate, se, level = 0.95){
  if(!is.finite(estimate) || estimate <= 0 || !is.finite(se) || se < 0){
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level)/2)
  exp(log(estimate) + c(-1, 1) * z * se/estimate)
}

motif.rare.case <- function(diagnostics, success_rate = NULL,
                            ess_ratio_threshold = 0.1,
                            positive_rate_threshold = 0.1,
                            max_weight_share_threshold = 0.1){
  rare <- FALSE
  if(!is.null(success_rate) && is.finite(success_rate)){
    rare <- success_rate < positive_rate_threshold
  }
  if(!is.null(diagnostics)){
    rare <- rare ||
      (is.finite(diagnostics$ess_ratio) &&
         diagnostics$ess_ratio < ess_ratio_threshold) ||
      (is.finite(diagnostics$positive_rate) &&
         diagnostics$positive_rate < positive_rate_threshold) ||
      (is.finite(diagnostics$max_weight_share) &&
         diagnostics$max_weight_share > max_weight_share_threshold)
  }
  rare
}

motif.edge.logprob <- function(p, has_edge){
  if(has_edge){
    out <- log(p)
    out[p <= 0] <- -Inf
  }else{
    out <- log1p(-p)
    out[p >= 1] <- -Inf
  }
  out
}

motif.sample.graph <- function(P, directed = FALSE){
  P <- motif.norm.probmat(P, directed = directed)
  n <- nrow(P)
  G <- matrix(FALSE, n, n)
  if(directed){
    mask <- row(P) != col(P)
    G[mask] <- as.logical(stats::rbinom(sum(mask), 1, P[mask]))
    diag(G) <- FALSE
    return(G)
  }
  mask <- upper.tri(P)
  G[mask] <- as.logical(stats::rbinom(sum(mask), 1, P[mask]))
  G <- G | t(G)
  diag(G) <- FALSE
  G
}

motif.wald.p <- function(z, alternative = c("greater", "less", "two.sided")){
  alternative <- match.arg(alternative)
  if(alternative == "greater"){return(stats::pnorm(z, lower.tail = FALSE))}
  if(alternative == "less"){return(stats::pnorm(z, lower.tail = TRUE))}
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}
