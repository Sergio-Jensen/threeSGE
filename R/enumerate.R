vecs <- function(x) {
  x <- as.integer(x)
  if(length(x) == 0L) return(list(integer(0)))
  if(length(x) == 1L) return(list(x))
  out <- list()
  for(i in seq_along(x)) {
    for(p in vecs(x[-i])) out[[length(out) + 1L]] <- c(x[i], p)
  }
  out
}

combs <- function(x, m) {
  if(m == 0L) return(list(integer(0)))
  utils::combn(x, m, simplify = FALSE)
}

consts <- function(H, directed = FALSE, induced = TRUE) {
  H <- as.matrix(H)
  k <- nrow(H)
  A <- 1L * (H != 0)
  diag(A) <- 0L
  
  if(!directed) {
    A <- 1L * ((A + t(A)) > 0)
    ij <- which(upper.tri(A), arr.ind = TRUE)
  } else {
    ij <- which(row(A) != col(A), arr.ind = TRUE)
  }
  
  cval <- A[ij]
  if(!induced) {
    keep <- cval == 1L
    ij <- ij[keep, , drop = FALSE]
    cval <- cval[keep]
  }
  
  data.frame(a = ij[, 1], b = ij[, 2], c = as.integer(cval))
}

const.aut <- function(cons, u, directed = FALSE) {
  if(u <= 1L) return(1L)
  perms <- vecs(seq_len(u))

  if(nrow(cons) == 0L) return(length(perms))

  if(!directed) {
    key0 <- paste(pmin(cons$a, cons$b), pmax(cons$a, cons$b), cons$c, sep = "_")
  } else {
    key0 <- paste(cons$a, cons$b, cons$c, sep = "_")
  }
  key0 <- sort(key0)

  sum(vapply(perms, function(p) {
    aa <- p[cons$a]
    bb <- p[cons$b]
    if(!directed) {
      keyp <- paste(pmin(aa, bb), pmax(aa, bb), cons$c, sep = "_")
    } else {
      keyp <- paste(aa, bb, cons$c, sep = "_")
    }
    identical(sort(keyp), key0)
  }, logical(1)))
}

aut <- function(H, directed = FALSE) {
  H <- as.matrix(H)
  k <- nrow(H)
  A <- 1L * (H != 0)
  diag(A) <- 0L
  if(!directed) A <- 1L * ((A + t(A)) > 0)

  perms <- vecs(seq_len(k))
  sum(vapply(perms, function(p) {
    identical(A, A[p, p, drop = FALSE])
  }, logical(1)))
}

overlap <- function(k, overlap_sizes = 2:k) {
  patterns <- list()
  for(s in overlap_sizes) {
    for(S1 in combs(seq_len(k), s)) {
      for(S2 in combs(seq_len(k), s)) {
        for(S2p in vecs(S2)) {
          rho <- stats::setNames(as.integer(S2p), as.character(S1))
          patterns[[length(patterns) + 1L]] <- list(s = s, rho = rho)
        }
      }
    }
  }
  patterns
}

u.consts <- function(H, rho, directed = FALSE, induced = TRUE) {
  k <- nrow(H)
  parent <- seq_len(2L * k)
  
  find <- function(x) {
    while(parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  
  unite <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if(ra != rb) parent[rb] <<- ra
  }
  
  if(length(rho) > 0L) {
    for(i in seq_along(rho)) {
      a <- as.integer(names(rho)[i])
      b <- as.integer(rho[i])
      unite(a, k + b)
    }
  }
  
  reps <- vapply(seq_len(2L * k), find, integer(1))
  ids <- stats::setNames(seq_along(sort(unique(reps))), sort(unique(reps)))
  node_id <- as.integer(ids[as.character(reps)])
  
  copy1 <- node_id[seq_len(k)]
  copy2 <- node_id[k + seq_len(k)]
  u <- length(unique(node_id))
  
  base <- consts(H, directed = directed, induced = induced)
  
  raw <- rbind(
    data.frame(a = copy1[base$a], b = copy1[base$b], c = base$c),
    data.frame(a = copy2[base$a], b = copy2[base$b], c = base$c)
  )
  
  if(!directed && nrow(raw) > 0L) {
    aa <- pmin(raw$a, raw$b)
    bb <- pmax(raw$a, raw$b)
    raw$a <- aa
    raw$b <- bb
  }
  
  if(nrow(raw) == 0L) {
    return(list(
      u = u,
      constraints = raw,
      multiset_constraints = raw,
      incompatible = FALSE
    ))
  }
  
  key <- paste(raw$a, raw$b, sep = "_")
  parts <- split(seq_len(nrow(raw)), key)
  
  incompatible <- any(vapply(parts, function(id) {
    length(unique(raw$c[id])) > 1L
  }, logical(1)))
  
  if(incompatible) {
    return(list(
      u = u,
      constraints = raw[0, ],
      multiset_constraints = raw,
      incompatible = TRUE
    ))
  }
  
  keep <- vapply(parts, function(id) id[1L], integer(1))
  cons <- raw[keep, , drop = FALSE]
  rownames(cons) <- NULL
  
  list(
    u = u,
    constraints = cons,
    multiset_constraints = raw,
    incompatible = FALSE
  )
}

consts.adj <- function(cons, u, directed = FALSE, value_col = "c") {
  A <- matrix(0L, u, u)
  if(nrow(cons) == 0L) return(A)
  
  edge_rows <- which(cons[[value_col]] == 1L)
  if(length(edge_rows) == 0L) return(A)
  
  for(ii in edge_rows) {
    a <- cons$a[ii]
    b <- cons$b[ii]
    A[a, b] <- 1L
    if(!directed) A[b, a] <- 1L
  }
  
  diag(A) <- 0L
  A
}

overlapG <- function(H,
                            directed = FALSE,
                            induced = TRUE,
                            overlap_sizes = 2:nrow(H),
                            include_incompatible = FALSE) {
  H <- as.matrix(H)
  k <- nrow(H)
  patterns <- overlap(k, overlap_sizes)
  
  out <- vector("list", length(patterns))
  keep <- rep(TRUE, length(patterns))
  
  for(id in seq_along(patterns)) {
    pat <- patterns[[id]]
    U <- u.consts(H, pat$rho, directed = directed, induced = induced)
    
    if(U$incompatible && !include_incompatible) {
      keep[id] <- FALSE
      next
    }
    
    adjacency <- if(U$incompatible) {
      NULL
    } else {
      consts.adj(U$constraints, U$u, directed = directed)
    }
    
    edge_constraints <- U$multiset_constraints[
      U$multiset_constraints$c == 1L, , drop = FALSE
    ]
    nonedge_constraints <- U$multiset_constraints[
      U$multiset_constraints$c == 0L, , drop = FALSE
    ]
    
    out[[id]] <- list(
      pattern_id = id,
      overlap_size = pat$s,
      rho = pat$rho,
      union_nodes = U$u,
      incompatible = U$incompatible,
      adjacency = adjacency,
      constraints = U$constraints,
      multiset_constraints = U$multiset_constraints,
      edge_constraints = edge_constraints,
      nonedge_constraints = nonedge_constraints
    )
  }
  
  out <- out[keep]
  names(out) <- paste0("pattern_", vapply(out, function(x) x$pattern_id, integer(1)))
  out
}

enum.u <- function(X, cons, u, directed = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  if(u > n) return(0)
  if(nrow(cons) == 0L) return(prod(n - seq_len(u) + 1))
  
  deg <- tabulate(c(cons$a, cons$b), nbins = u)
  ord <- order(-deg, seq_len(u))
  
  involved <- lapply(seq_len(u), function(r) which(cons$a == r | cons$b == r))
  assigned <- integer(u)
  used <- rep(FALSE, n)
  total <- 0
  
  get.val <- function(i, j) {
    if(!directed && i > j) X[j, i] else X[i, j]
  }
  
  rec <- function(t, weight) {
    if(t > u) {
      total <<- total + weight
      return(NULL)
    }
    
    r <- ord[t]
    for(x in which(!used)) {
      inc <- 1
      
      for(ii in involved[[r]]) {
        a <- cons$a[ii]; b <- cons$b[ii]
        other <- if(a == r) b else a
        if(assigned[other] == 0L) next
        
        i <- if(a == r) x else assigned[a]
        j <- if(b == r) x else assigned[b]
        
        pij <- get.val(i, j)
        inc <- inc * if(cons$c[ii] == 1L) pij else (1 - pij)
        if(inc == 0) break
      }
      
      if(inc != 0) {
        assigned[r] <<- x
        used[x] <<- TRUE
        rec(t + 1L, weight * inc)
        used[x] <<- FALSE
        assigned[r] <<- 0L
      }
    }
  }
  
  rec(1L, 1)
  total
}

PairUnion <- function(X, H,
                                directed = FALSE,
                                induced = TRUE,
                                overlap_sizes = 2:nrow(H),
                                mode = c("count", "expectation"),
                                include_incompatible = FALSE,
                                method = c("auto", "enumerate", "homogeneous"),
                                homogeneous_tol = 1e-12) {
  mode <- match.arg(mode)
  method <- match.arg(method)
  X <- as.matrix(X)
  
  if(mode == "count" && any(!(X %in% c(0, 1)))) {
    stop("For mode = 'count', X must be a 0/1 adjacency matrix.")
  }
  if(mode == "expectation" && any(X < 0 | X > 1)) {
    stop("For mode = 'expectation', X must be a probability matrix with entries in [0,1].")
  }
  
  k <- nrow(H)
  n <- nrow(X)
  offdiag <- if(directed) X[row(X) != col(X)] else X[upper.tri(X)]
  is_homogeneous <- length(offdiag) > 0L &&
    all(is.finite(offdiag)) &&
    max(abs(offdiag - offdiag[1L])) <= homogeneous_tol
  if(method == "homogeneous" && !is_homogeneous) {
    stop("method = 'homogeneous' requires all relevant off-diagonal entries of X to be equal.")
  }
  use_homogeneous <- method == "homogeneous" || (method == "auto" && is_homogeneous)
  homogeneous_p <- if(use_homogeneous) offdiag[1L] else NA_real_

  patterns <- overlap(k, overlap_sizes)
  motif_aut <- aut(H, directed = directed)

  out <- list()
  labeled_total <- 0
  unlabeled_union_total <- 0
  unlabeled_pair_total <- 0
  
  for(id in seq_along(patterns)) {
    pat <- patterns[[id]]
    U <- u.consts(H, pat$rho, directed, induced)
    
    value <- 0
    union_aut <- NA_integer_
    if(!U$incompatible) {
      if(use_homogeneous) {
        edge_n <- sum(U$constraints$c == 1L)
        nonedge_n <- sum(U$constraints$c == 0L)
        embedding_n <- if(U$u > n) {
          0
        } else if(U$u == 0L) {
          1
        } else {
          prod(n - seq.int(0L, U$u - 1L))
        }
        value <- embedding_n * (homogeneous_p ^ edge_n) *
          ((1 - homogeneous_p) ^ nonedge_n)
      } else {
        value <- enum.u(X, U$constraints, U$u, directed)
      }
      union_aut <- const.aut(U$constraints, U$u, directed)
    }

    unlabeled_union_value <- if(is.na(union_aut)) 0 else value / union_aut
    unlabeled_pair_value <- value / (motif_aut^2)

    labeled_total <- labeled_total + value
    unlabeled_union_total <- unlabeled_union_total + unlabeled_union_value
    unlabeled_pair_total <- unlabeled_pair_total + unlabeled_pair_value
    
    if(include_incompatible || !U$incompatible) {
      out[[length(out) + 1L]] <- data.frame(
        pattern_id = id,
        overlap_size = pat$s,
        union_nodes = U$u,
        incompatible = U$incompatible,
        union_automorphism_count = union_aut,
        motif_automorphism_count = motif_aut,
        labeled_value = value,
        unlabeled_union_value = unlabeled_union_value,
        unlabeled_pair_value = unlabeled_pair_value
      )
    }
  }
  
  by_pattern <- do.call(rbind, out)
  suffix <- if(mode == "count") "count" else "expectation"
  names(by_pattern)[names(by_pattern) == "labeled_value"] <- paste0("labeled_", suffix)
  names(by_pattern)[names(by_pattern) == "unlabeled_union_value"] <- paste0("unlabeled_union_", suffix)
  names(by_pattern)[names(by_pattern) == "unlabeled_pair_value"] <- paste0("unlabeled_pair_", suffix)
  
  list(
    total = labeled_total,
    labeled_total = labeled_total,
    unlabeled_union_total = unlabeled_union_total,
    unlabeled_pair_total = unlabeled_pair_total,
    motif_automorphism_count = motif_aut,
    by_pattern = by_pattern,
    mode = mode,
    method = if(use_homogeneous) "homogeneous" else "enumerate",
    homogeneous_probability = homogeneous_p,
    directed = directed,
    induced = induced,
    overlap_sizes = overlap_sizes
  )
}

PairVar <- function(P, H,
                                       directed = FALSE,
                                       induced = TRUE,
                                       overlap_sizes = 2:nrow(H),
                                       include_incompatible = TRUE,
                                       method = c("auto", "enumerate", "homogeneous"),
                                       homogeneous_tol = 1e-12) {
  method <- match.arg(method)
  P <- as.matrix(P)
  if(any(P < 0 | P > 1)) {
    stop("P must be a probability matrix with entries in [0,1].")
  }

  k <- nrow(H)
  n <- nrow(P)
  offdiag <- if(directed) P[row(P) != col(P)] else P[upper.tri(P)]
  is_homogeneous <- length(offdiag) > 0L &&
    all(is.finite(offdiag)) &&
    max(abs(offdiag - offdiag[1L])) <= homogeneous_tol
  if(method == "homogeneous" && !is_homogeneous) {
    stop("method = 'homogeneous' requires all relevant off-diagonal entries of P to be equal.")
  }
  use_homogeneous <- method == "homogeneous" || (method == "auto" && is_homogeneous)
  homogeneous_p <- if(use_homogeneous) offdiag[1L] else NA_real_

  patterns <- overlap(k, overlap_sizes)
  motif_aut <- aut(H, directed = directed)

  out <- list()
  labeled_variance <- 0
  joint_total <- 0
  marginal_product_total <- 0

  for(id in seq_along(patterns)) {
    pat <- patterns[[id]]
    U <- u.consts(H, pat$rho, directed, induced)

    joint_sum <- 0
    if(!U$incompatible) {
      if(use_homogeneous) {
        joint_edges <- sum(U$constraints$c == 1L)
        joint_nonedges <- sum(U$constraints$c == 0L)
        embedding_n <- if(U$u > n) {
          0
        } else if(U$u == 0L) {
          1
        } else {
          prod(n - seq.int(0L, U$u - 1L))
        }
        joint_sum <- embedding_n * (homogeneous_p ^ joint_edges) *
          ((1 - homogeneous_p) ^ joint_nonedges)
      } else {
        joint_sum <- enum.u(P, U$constraints, U$u, directed)
      }
    }

    if(use_homogeneous) {
      product_edges <- sum(U$multiset_constraints$c == 1L)
      product_nonedges <- sum(U$multiset_constraints$c == 0L)
      embedding_n <- if(U$u > n) {
        0
      } else if(U$u == 0L) {
        1
      } else {
        prod(n - seq.int(0L, U$u - 1L))
      }
      marginal_product_sum <- embedding_n * (homogeneous_p ^ product_edges) *
        ((1 - homogeneous_p) ^ product_nonedges)
    } else {
      marginal_product_sum <- enum.u(
        P,
        U$multiset_constraints,
        U$u,
        directed
      )
    }

    covariance_sum <- joint_sum - marginal_product_sum
    unlabeled_covariance_sum <- covariance_sum / (motif_aut^2)

    labeled_variance <- labeled_variance + covariance_sum
    joint_total <- joint_total + joint_sum
    marginal_product_total <- marginal_product_total + marginal_product_sum

    if(include_incompatible || !U$incompatible) {
      out[[length(out) + 1L]] <- data.frame(
        pattern_id = id,
        overlap_size = pat$s,
        union_nodes = U$u,
        incompatible = U$incompatible,
        motif_automorphism_count = motif_aut,
        joint_sum = joint_sum,
        marginal_product_sum = marginal_product_sum,
        labeled_covariance_sum = covariance_sum,
        unlabeled_covariance_sum = unlabeled_covariance_sum
      )
    }
  }

  by_pattern <- do.call(rbind, out)
  unlabeled_variance <- labeled_variance / (motif_aut^2)

  list(
    labeled_variance = labeled_variance,
    unlabeled_variance = unlabeled_variance,
    labeled_joint_total = joint_total,
    labeled_marginal_product_total = marginal_product_total,
    motif_automorphism_count = motif_aut,
    method = if(use_homogeneous) "homogeneous" else "enumerate",
    homogeneous_probability = homogeneous_p,
    by_pattern = by_pattern,
    directed = directed,
    induced = induced,
    overlap_sizes = overlap_sizes
  )
}
