\name{iss.sensitivity}
\alias{iss.sensitivity}
\alias{print.iss.sensitivity}
\title{Assess ISS Sensitivity to the Sequential Vertex Order}
\description{
Repeats binary ISS estimation under several motif-vertex orders and compares
the estimates, standard errors, success rates, and weight diagnostics.
}
\usage{
iss.sensitivity(G, H, B = 10000, directed = FALSE, seed = NULL,
                order = NULL, orders = c("mot.order", "reverse"),
                random_orders = 0, aut_count = NULL, verbose = FALSE,
                bootstrap_ci = FALSE, bootstrap_B = 1000,
                verbose.boot = FALSE, keep_results = TRUE,
                verbose.sens = TRUE)
}
\arguments{
  \item{G}{Square binary host adjacency matrix or graph object.}
  \item{H}{Square binary motif adjacency matrix or graph object.}
  \item{B}{Positive number of ISS draws for each tested order.}
  \item{directed}{Logical; whether edges are directed.}
  \item{seed}{Optional integer random seed. Each order receives a distinct,
  reproducible derived seed.}
  \item{order}{Optional baseline permutation of \code{seq_len(nrow(H))}.}
  \item{orders}{Character vector containing \code{"mot.order"},
  \code{"reverse"}, or \code{"given"}. \code{"default"} and
  \code{"rarity"} are accepted synonyms for \code{"mot.order"}.}
  \item{random_orders}{Nonnegative integer number of additional random
  permutations.}
  \item{aut_count}{Optional known motif automorphism count.}
  \item{verbose}{Logical; print sampling progress inside each ISS run.}
  \item{bootstrap_ci}{Logical; compute bootstrap intervals in each run.}
  \item{bootstrap_B}{Positive number of bootstrap resamples.}
  \item{verbose.boot}{Logical; print bootstrap progress.}
  \item{keep_results}{Logical; retain every full \code{est.subgraph} result.}
  \item{verbose.sens}{Logical; print completion messages for tested orders.}
}
\value{
An object of class \code{"iss.sensitivity"} with
\code{sensitivity_table}, one row per tested order; \code{results}, a named
list of full estimates when \code{keep_results = TRUE}; and metadata
\code{B}, \code{directed}, \code{graph_density}, and \code{random_orders}.
The table reports the order, unlabeled estimate and SE, Wald and log-Wald
limits, success rate, ESS ratio, weight CV, and maximum-weight share.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
G <- matrix(1, 6, 6); diag(G) <- 0
ans <- iss.sensitivity(G, H, B = 100, random_orders = 1,
                       seed = 2, verbose.sens = FALSE)
ans$sensitivity_table
}
\seealso{\code{\link{est.subgraph}}}
