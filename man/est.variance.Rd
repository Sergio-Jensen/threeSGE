\name{est.variance}
\alias{est.variance}
\alias{print.subgraph.var}
\title{Estimate the Variance of a Subgraph Count}
\description{
Uses pairwise-union ISS to estimate the variance of the labeled and
unlabeled count of a motif in an independent-edge probability network.
Both compatible joint mass (FPUG) and marginal-product mass (IPUM) are
sampled within each overlap layer.
}
\usage{
est.variance(P, H, B = 10000, B_mean = B, directed = FALSE,
             induced = TRUE, normalize_weight = FALSE,
             overlap_sizes = NULL,
             allocation = c("equal", "proportional"),
             mu_labeled = NULL, mu_unlabeled = NULL, aut_count = NULL,
             seed = NULL, verbose = TRUE, verbose.B = TRUE,
             include_incompatible = TRUE, preprocessed = FALSE,
             ci_level = 0.95, bootstrap_ci = FALSE, bootstrap_B = 1000,
             bootstrap_seed = NULL, verbose.boot = FALSE)
}
\arguments{
  \item{P}{Square independent-edge probability matrix, \code{igraph}, or
  \code{tbl_graph}. Entries must be in \eqn{[0,1]} unless
  \code{normalize_weight = TRUE}.}
  \item{H}{Square motif adjacency matrix or graph object; nonzero entries
  indicate motif edges.}
  \item{B}{Positive total number of pairwise-union ISS draws.}
  \item{B_mean}{Positive number of ISS draws used to estimate the motif mean
  when neither \code{mu_labeled} nor \code{mu_unlabeled} is supplied.}
  \item{directed}{Logical; whether edges are directed.}
  \item{induced}{Logical; impose both motif-edge and motif-nonedge
  constraints. If \code{FALSE}, only motif edges are required.}
  \item{normalize_weight}{If \code{TRUE}, divide nonnegative raw weights by
  their maximum. If \code{FALSE}, use \code{P} directly as probabilities.}
  \item{overlap_sizes}{Integer vector of shared-vertex counts. The default is
  \code{2:k} for a \eqn{k}-vertex motif (or \code{k} when \eqn{k<=1}).}
  \item{allocation}{\code{"equal"} or \code{"proportional"} allocation of
  \code{B} across overlap layers.}
  \item{mu_labeled, mu_unlabeled}{Optional known labeled or unlabeled motif
  mean. Supply either or both to avoid or reduce separate mean estimation.}
  \item{aut_count}{Optional known motif automorphism count.}
  \item{seed}{Optional integer random seed.}
  \item{verbose}{Logical; print overlap-layer completion messages.}
  \item{verbose.B}{Logical; print within-layer sampling progress.}
  \item{include_incompatible}{Logical; include incompatible overlap patterns
  when estimating marginal-product mass. This should normally remain
  \code{TRUE} for complete variance estimation.}
  \item{preprocessed}{Logical; treat \code{P} as already validated and
  normalized. Intended for internal repeated calculations.}
  \item{ci_level}{Confidence level in \eqn{(0,1)}.}
  \item{bootstrap_ci}{Logical; add stratified nonparametric bootstrap
  intervals.}
  \item{bootstrap_B}{Positive number of bootstrap resamples.}
  \item{bootstrap_seed}{Optional bootstrap seed.}
  \item{verbose.boot}{Logical; print bootstrap progress.}
}
\value{
An object with classes \code{"variance"} and \code{"subgraph.var"}.
\code{mean} contains labeled/unlabeled means and the optional internal mean
result. \code{variance} contains raw and nonnegative-truncated labeled and
unlabeled variance estimates plus Monte Carlo SEs. \code{overdispersion}
contains raw and nonnegative overdispersion coefficients and the resulting
Poisson-overdispersed variance. \code{layer_table} reports pair counts,
compatible/incompatible pattern counts, allocations, success rates, joint
overlap mass, marginal-product mass, and variance contributions by overlap.
\code{clumping} contains overall and by-overlap structural summaries.
\code{inference} contains confidence intervals and separate diagnostics for
variance, joint-overlap, and marginal-product weights. Metadata fields record
\code{B}, \code{B_mean}, graph sizes, probability rule, density,
automorphism count, overlap sizes, and allocation.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 6, 6); diag(P) <- 0
ans <- est.variance(P, H, B = 300, B_mean = 300,
                    directed = FALSE, seed = 4,
                    verbose = FALSE, verbose.B = FALSE)
ans$variance$unlabeled  # theoretical value is 5 for this example
ans$layer_table
}
\seealso{\code{\link{est.pairunion}}, \code{\link{sec.moment}}}
