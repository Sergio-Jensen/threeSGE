\name{tsge}
\alias{tsge}
\alias{print.tsge}
\alias{summary.tsge}
\title{Three-Stage Graph Estimation Motif Test}
\description{
Compares an observed motif count with its expectation under a fitted or
supplied null probability network. Binary and weighted observations use one
interface and are dispatched automatically.
}
\usage{
tsge(G, H, network_type = c("auto", "binary", "weighted"),
     normalize_weight = TRUE, ...)
}
\arguments{
  \item{G}{Observed square network matrix, \code{igraph}, or
  \code{tbl_graph}. Binary off-diagonal values select the bitset counting
  path; other numeric values select the weighted expectation path.}
  \item{H}{Square motif adjacency matrix or graph object. Nonzero entries
  denote required motif edges.}
  \item{network_type}{\code{"auto"}, \code{"binary"}, or
  \code{"weighted"}. Use an explicit value to override detection.}
  \item{normalize_weight}{Weighted observations only. If \code{TRUE}, divide
  raw nonnegative weights by their maximum. If \code{FALSE}, use values as
  probabilities in \eqn{[0,1]}.}
  \item{...}{Test controls described below.}
}
\details{
The following arguments may be passed through \code{...}:
\describe{
  \item{\code{B_observed}}{Positive ISS sample size for the observed count or
  expectation; default 10000.}
  \item{\code{R_null}}{Positive ISS sample size for the null expectation;
  default 10000.}
  \item{\code{directed}}{Logical; binary default \code{FALSE}, weighted
  default \code{TRUE}.}
  \item{\code{null_model}}{Binary choices are \code{"chung_lu"},
  \code{"erdos_renyi"}, and \code{"given"}. Weighted choices are
  \code{"weighted_chung_lu"}, \code{"weighted_erdos_renyi"},
  \code{"given"}, with \code{"chung_lu"} and \code{"erdos_renyi"}
  accepted as aliases.}
  \item{\code{P_null}}{Square probability matrix required when
  \code{null_model = "given"}.}
  \item{\code{null_variance}}{How to represent null count variation.
  \code{"poisson"} uses the null mean; \code{"sample"} calls
  \code{est.variance()}; \code{"direct"} is an alias for \code{"sample"};
  \code{"none"} uses zero; and \code{"provided"} uses \code{tau2_null}.
  Weighted tests additionally allow \code{"calibrated"}, which estimates
  variance from simulated null networks. Binary default is \code{"poisson"};
  weighted default is \code{"calibrated"}.}
  \item{\code{tau2_null}}{Nonnegative supplied null variance required for
  \code{null_variance = "provided"}.}
  \item{\code{K_null_variance}}{Weighted calibrated mode only; integer at
  least 2 giving the number of simulated null networks; default 99.}
  \item{\code{B_null_variance}}{Positive ISS sample size for sampled or
  calibrated null variance; defaults to \code{R_null}.}
  \item{\code{null_variance_overlap_sizes}}{Optional overlap-size vector
  passed to \code{est.variance()}.}
  \item{\code{null_variance_allocation}}{\code{"equal"} or
  \code{"proportional"}; passed to \code{est.variance()}.}
  \item{\code{alternative}}{\code{"greater"}, \code{"less"}, or
  \code{"two.sided"}.}
  \item{\code{seed}}{Optional integer random seed.}
  \item{\code{order}}{Optional permutation of motif vertices for sequential
  sampling.}
  \item{\code{automorphism_count_value}}{Optional known motif automorphism
  count.}
  \item{\code{verbose}}{Logical; print ISS and calibration progress.}
  \item{\code{null_expectation}}{Optional previously computed
  \code{est.subgraph()} result for the null network. It can be reused in both
  binary and weighted TSGE.}
  \item{\code{observed}}{Optional previously computed
  \code{est.subgraph()} result for the observed network. It can be reused in
  both binary and weighted TSGE.}
}
The primary coefficient is
\eqn{\log(\widehat\theta/\widehat\mu_0)}. Its total standard error combines
Monte Carlo error from both ISS estimates and the selected null-count
variance. The auxiliary statistic is the count difference
\eqn{\widehat\theta-\widehat\mu_0}.
}
\value{
An object of class \code{"tsge"}. \code{coefficient_table} reports log
enrichment, total SE, z statistic, and p-value. \code{auxiliary_table}
reports the count-difference analogue. \code{estimates} contains observed
and null counts, their Monte Carlo SEs, enrichment ratio, log enrichment,
Monte Carlo and total SEs, and excess-count statistics. \code{observed} and
\code{null_expectation} retain the underlying ISS results.
\code{sample_null_variance}, \code{calibrated_null_counts},
\code{calibrated_null_ses}, and \code{calibrated_tau2} are populated when
their corresponding variance method is used. Other fields record the null
model, variance method, alternative, probability rule, motif order,
automorphism count, and detected \code{network_type}.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
G <- matrix(c(
  0,1,1,0,0,
  1,0,1,1,0,
  1,1,0,1,1,
  0,1,1,0,1,
  0,0,1,1,0
), 5, 5, byrow = TRUE)

fit <- tsge(G, H, B_observed = 200, R_null = 200,
            directed = FALSE, null_model = "erdos_renyi",
            null_variance = "poisson", seed = 6, verbose = FALSE)
fit$coefficient_table

# Reuse both binary-branch estimates without repeating ISS.
fit_reused <- tsge(G, H, observed = fit$observed,
                   null_expectation = fit$null_expectation,
                   directed = FALSE, null_model = "erdos_renyi",
                   null_variance = "poisson", verbose = FALSE)

W <- matrix(0.5, 5, 5); diag(W) <- 0
fit_w <- tsge(W, H, network_type = "weighted",
              normalize_weight = FALSE, directed = FALSE,
              null_model = "given", P_null = W,
              null_variance = "none", B_observed = 200, R_null = 200,
              seed = 6, verbose = FALSE)
fit_w$estimates
}
\seealso{\code{\link{tsge.cal}}, \code{\link{est.subgraph}},
\code{\link{est.variance}}}
