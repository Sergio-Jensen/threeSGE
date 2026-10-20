\name{tsge.cal}
\alias{tsge.cal}
\alias{print.tsge.cal}
\alias{summary.tsge.cal}
\title{Calibrate a TSGE Motif Statistic with Random Networks}
\description{
Simulates random networks under a fitted, supplied, or degree-preserving null
model, re-estimates the motif count in every network, and reports empirical
and normal-reference calibration results.
}
\usage{
tsge.cal(G, H, network_type = c("auto", "binary", "weighted"),
         normalize_weight = TRUE, ...)
}
\arguments{
  \item{G}{Observed square network matrix, \code{igraph}, or
  \code{tbl_graph}.}
  \item{H}{Square motif adjacency matrix or graph object.}
  \item{network_type}{\code{"auto"}, \code{"binary"}, or
  \code{"weighted"}.}
  \item{normalize_weight}{Weighted input only; divide raw weights by their
  maximum when \code{TRUE}, otherwise use them as probabilities.}
  \item{...}{Calibration controls described below.}
}
\details{
Arguments passed through \code{...} are:
\describe{
  \item{\code{enrichment_result}}{Optional existing result from
  \code{tsge()}. Its observed estimate, null model, motif order,
  automorphism count, and null expectation are reused unless overridden.}
  \item{\code{observed}}{Optional existing \code{est.subgraph()} result for
  the observed network. This is supported for both binary and weighted
  networks and avoids recomputing the observed count.}
  \item{\code{null_expectation}}{Optional existing \code{est.subgraph()}
  result for the null-network expectation. This is supported for both binary
  and weighted networks and avoids recomputing the scalar null expectation.}
  \item{\code{K}}{Positive number of simulated random networks; default 99.
  The smallest attainable empirical p-value is \eqn{1/(K+1)}.}
  \item{\code{R_null}}{Positive ISS sample size for each random network;
  default 10000.}
  \item{\code{B_observed}}{Positive ISS sample size for the observed network
  when no reusable \code{enrichment_result} is supplied; default 10000.}
  \item{\code{directed}}{Logical; binary default \code{FALSE}, weighted
  default \code{TRUE}.}
  \item{\code{null_model}}{Binary choices are \code{"chung_lu"},
  \code{"erdos_renyi"}, \code{"given"}, and \code{"degree_swap"}; default
  \code{"chung_lu"}. Weighted choices are \code{"weighted_chung_lu"},
  \code{"weighted_erdos_renyi"}, and \code{"given"}, with unprefixed model
  names accepted as aliases; default \code{"weighted_chung_lu"}.}
  \item{\code{P_null}}{Square probability matrix required by
  \code{null_model = "given"}.}
  \item{\code{statistic}}{\code{"log_enrichment"}, \code{"count"}, or
  \code{"excess"}.}
  \item{\code{alternative}}{\code{"greater"}, \code{"less"}, or
  \code{"two.sided"}.}
  \item{\code{seed}}{Optional integer random seed.}
  \item{\code{order}}{Optional motif-vertex permutation.}
  \item{\code{automorphism_count_value}}{Optional known automorphism count.}
  \item{\code{null_expectation_value}}{Optional numeric null expectation.
  If absent, the average random-network count is used. A supplied
  \code{null_expectation} result takes precedence over this legacy numeric
  argument.}
  \item{\code{degree_swap_steps}}{Binary degree-swap mode only; requested
  number of successful edge swaps per random network. By default it is ten
  times the observed edge count.}
  \item{\code{log_offset}}{Nonnegative offset added before taking logarithms.
  Use a positive value when zero motif counts are possible.}
  \item{\code{verbose}}{Logical; print ISS progress for observed and random
  networks.}
  \item{\code{verbose.K}}{Logical; print progress across the \code{K}
  random networks.}
}
Explicit \code{observed} and \code{null_expectation} objects take precedence
over the corresponding objects or scalar values in \code{enrichment_result}.
}
\value{
An object of class \code{"tsge.cal"}. \code{calibration_table} contains the
observed statistic, random-network mean and SD, calibrated z statistic, and
normal-reference p-value. \code{empirical.p} is the finite-simulation
empirical p-value and \code{min_p_value} is \eqn{1/(K+1)}.
\code{random_counts}, \code{random_count_ses}, and
\code{random_statistics} retain all calibration draws. Other fields include
the observed count and SE, scalar null expectation, retained \code{observed}
estimate, retained \code{null_expectation_result} estimate, statistic,
alternative, null model, \code{K}, \code{R_null}, motif order, automorphism
count, probability rule, and detected network type.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 5, 5); diag(P) <- 0
weighted_est <- est.subgraph(P, H, network_type = "weighted",
                             normalize_weight = FALSE, directed = FALSE,
                             B = 100, seed = 6, verbose = FALSE)
cal <- tsge.cal(P, H, network_type = "weighted",
                normalize_weight = FALSE, directed = FALSE,
                null_model = "given", P_null = P,
                observed = weighted_est,
                null_expectation = weighted_est,
                K = 3, R_null = 100,
                statistic = "count", seed = 7,
                verbose = FALSE, verbose.K = FALSE)
cal$calibration_table
cal$empirical.p

G <- matrix(1, 5, 5); diag(G) <- 0
binary_est <- est.subgraph(G, H, directed = FALSE, B = 100,
                           seed = 8, verbose = FALSE)
binary_cal <- tsge.cal(G, H, directed = FALSE,
                       null_model = "given", P_null = P,
                       observed = binary_est,
                       null_expectation = weighted_est,
                       K = 3, R_null = 100, statistic = "count",
                       seed = 9, verbose = FALSE, verbose.K = FALSE)
binary_cal$observed_count
}
\seealso{\code{\link{tsge}}}
