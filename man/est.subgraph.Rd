\name{est.subgraph}
\alias{est.subgraph}
\alias{print.est.subgraph}
\title{Estimate an Induced Subgraph Count or Expectation}
\description{
Uses induced sequential importance sampling (ISS) to estimate the number of
copies of a motif. Binary observations use the bitset implementation;
non-binary observations are interpreted as edge weights or probabilities.
}
\usage{
est.subgraph(G, H, network_type = c("auto", "binary", "weighted"),
             normalize_weight = TRUE, ...)
}
\arguments{
  \item{G}{Observed host network. Supply a square numeric or logical matrix,
  an \code{igraph} object, or a \code{tbl_graph}. For a binary network,
  off-diagonal entries must be 0/1. For a weighted network, entries must be
  nonnegative when \code{normalize_weight = TRUE}, or probabilities in
  \eqn{[0,1]} when it is \code{FALSE}.}
  \item{H}{Square motif adjacency matrix or graph object. Nonzero entries
  denote motif edges; the diagonal is ignored. The motif cannot have more
  vertices than \code{G}.}
  \item{network_type}{\code{"auto"} detects binary versus weighted input
  from off-diagonal values. \code{"binary"} and \code{"weighted"} force a
  backend.}
  \item{normalize_weight}{Weighted input only. If \code{TRUE}, construct
  \eqn{P_{ij}=W_{ij}/\max(W)}. If \code{FALSE}, treat \code{G} directly as
  a probability matrix.}
  \item{...}{Backend controls. \code{B} is a positive number of ISS draws
  (default 10000); \code{directed} selects directed motifs (binary default
  \code{FALSE}, weighted default \code{TRUE}); \code{seed} controls random
  sampling; \code{order} is a permutation of motif vertex indices;
  \code{keep_mappings} retains successful sampled embeddings;
  \code{verbose} prints progress; and \code{aut_count} supplies a known
  automorphism count. Binary-only controls are \code{order_method}, one of
  \code{"mot.order"}, \code{"reverse"}, or \code{"random"};
  \code{bootstrap_ci}; \code{bootstrap_B}; \code{bootstrap_seed}; and
  \code{verbose.boot}.}
}
\value{
An object of class \code{"est.subgraph"}. Common fields are \code{B},
\code{directed}, \code{graph_nodes}, \code{network_type}, the selected motif
order, \code{aut.morph}, \code{successes}, and \code{success_rate}.

For binary input, \code{labeled_estimate} and \code{unlabeled_estimate} are
the estimated counts, with \code{labeled_se} and \code{unlabeled_se}; Wald
and log-Wald intervals are returned in the corresponding \code{*_ci_95}
fields. \code{weight_diagnostics} contains ESS, coefficient-of-variation,
zero-rate, and maximum-weight-share diagnostics. \code{weights} stores the
ISS weights, \code{bootstrap_ci} stores optional bootstrap results, and
\code{mappings} is present when requested.

For weighted input, \code{labeled_expectation} and
\code{unlabeled_expectation} are expected counts, with standard errors in
\code{labeled_se} and \code{unlabeled_se}. \code{probability_matrix},
\code{probability_density}, \code{probability_rule}, \code{motif_order}, and
\code{log_weights} describe the probability-network calculation.
}
\details{
With \code{directed = FALSE}, symmetric matrices with zero diagonal are
strongly recommended. The unlabeled estimate equals the labeled estimate
divided by the motif automorphism count.
}
\examples{
H <- matrix(1, 3, 3)
diag(H) <- 0

G <- matrix(1, 5, 5)
diag(G) <- 0
fit_binary <- est.subgraph(G, H, B = 200, directed = FALSE,
                           seed = 1, verbose = FALSE)
fit_binary$unlabeled_estimate

P <- matrix(0.5, 5, 5)
diag(P) <- 0
fit_weighted <- est.subgraph(P, H, B = 200, directed = FALSE,
                             normalize_weight = FALSE,
                             seed = 1, verbose = FALSE)
fit_weighted$unlabeled_expectation
}
\seealso{\code{\link{iss.sensitivity}}, \code{\link{est.variance}}}
