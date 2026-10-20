\name{est.pairunion}
\alias{est.pairunion}
\alias{print.est.pairunion}
\title{Estimate Compatible Pairwise-Union Structure}
\description{
Estimates first- and second-order motif overlap mass by sampling compatible
pairwise unions at selected overlap sizes.
}
\usage{
est.pairunion(X, H, B = 10000, directed = FALSE, induced = TRUE,
              normalize_weight = FALSE, overlap_sizes = NULL,
              mode = c("count", "expectation"),
              allocation = c("equal", "proportional"),
              mu_labeled = NULL, mu_unlabeled = NULL, aut_count = NULL,
              include_self_overlap = FALSE, seed = NULL, verbose = TRUE,
              verbose.B = TRUE, preprocessed = FALSE, ci_level = 0.95,
              bootstrap_ci = FALSE, bootstrap_B = 1000,
              bootstrap_seed = NULL, verbose.boot = FALSE)
}
\arguments{
  \item{X}{Square host matrix or graph object. In \code{mode = "count"},
  nonzero entries are binary edges. In \code{mode = "expectation"}, entries
  are edge probabilities unless \code{normalize_weight = TRUE}.}
  \item{H}{Square motif adjacency matrix or graph object; nonzero entries
  are edges.}
  \item{B}{Positive total number of pairwise-union ISS draws.}
  \item{directed}{Logical; whether host and motif edges are directed.}
  \item{induced}{Logical; include motif non-edge constraints when compiling
  each union pattern.}
  \item{normalize_weight}{Expectation mode only. Divide nonnegative raw
  weights by their maximum before estimation. Otherwise entries must be in
  \eqn{[0,1]}.}
  \item{overlap_sizes}{Integer vector of numbers of shared motif vertices.
  The default is \code{2:(k - 1)} for a \eqn{k}-vertex motif.}
  \item{mode}{\code{"count"} for a realized binary graph or
  \code{"expectation"} for a probability network.}
  \item{allocation}{\code{"equal"} divides \code{B} approximately equally
  across overlap layers; \code{"proportional"} allocates according to the
  layer embedding counts.}
  \item{mu_labeled, mu_unlabeled}{Optional known motif means used to compute
  normalized clumping summaries. If one is supplied, the other is derived
  using \code{aut_count}.}
  \item{aut_count}{Optional known motif automorphism count.}
  \item{include_self_overlap}{Logical; allow full overlap \eqn{s=k}. The
  default excludes self-overlap for genuine pairwise-union structure.}
  \item{seed}{Optional integer random seed.}
  \item{verbose}{Logical; print layer completion messages.}
  \item{verbose.B}{Logical; print within-layer sampling progress.}
  \item{preprocessed}{Expectation mode only. If \code{TRUE}, regard \code{X}
  as an already validated probability matrix.}
  \item{ci_level}{Confidence level in \eqn{(0,1)} for reported intervals.}
  \item{bootstrap_ci}{Logical; compute stratified nonparametric bootstrap
  intervals.}
  \item{bootstrap_B}{Positive number of bootstrap resamples.}
  \item{bootstrap_seed}{Optional bootstrap seed.}
  \item{verbose.boot}{Logical; print bootstrap progress.}
}
\value{
An object of class \code{"est.pairunion"}. Scalar totals are returned in
\code{labeled_total}, \code{unlabeled_union_total}, and
\code{unlabeled_pair_total}, with corresponding \code{*_se} fields.
\code{layer_table} reports pattern counts, sample allocation, success rates,
and estimates by overlap size. \code{mean} records supplied means.
\code{clumping$summary} gives overlap mass, mean overlap size, concentration,
node/edge intensities, and effective subgraph counts;
\code{clumping$by_overlap} and \code{clumping$by_shared_edge} give detailed
breakdowns. \code{inference} contains Wald, log-Wald, optional bootstrap, and
importance-weight diagnostics for each reported total. Remaining fields
record the call, proposal, data mode, graph sizes, density, overlap sizes,
allocation, and automorphism count.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 6, 6); diag(P) <- 0
ans <- est.pairunion(P, H, B = 200, mode = "expectation",
                     directed = FALSE, seed = 3,
                     verbose = FALSE, verbose.B = FALSE)
ans$layer_table
}
\seealso{\code{\link{est.variance}}, \code{\link{sec.moment}}}
