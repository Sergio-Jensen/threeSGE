\name{PairUnion}
\alias{PairUnion}
\title{Exactly Enumerate Pairwise-Union Motif Structure}
\description{
Exactly computes labeled, unlabeled-union, and unlabeled-pair totals for two
copies of a motif sharing selected numbers of vertices. Unlike
\code{est.pairunion()}, this function performs no Monte Carlo sampling.
}
\usage{
PairUnion(X, H, directed = FALSE, induced = TRUE,
          overlap_sizes = 2:nrow(H),
          mode = c("count", "expectation"),
          include_incompatible = FALSE,
          method = c("auto", "enumerate", "homogeneous"),
          homogeneous_tol = 1e-12)
}
\arguments{
  \item{X}{Square numeric host matrix. With \code{mode = "count"}, every
  entry must be 0 or 1. With \code{mode = "expectation"}, every entry must
  be an edge probability in \eqn{[0,1]}. For undirected analyses, supply a
  symmetric matrix with zero diagonal.}
  \item{H}{Square motif adjacency matrix. Nonzero entries denote motif
  edges; the diagonal is ignored.}
  \item{directed}{Logical; whether ordered directed edges are used.}
  \item{induced}{Logical; if \code{TRUE}, both edges and nonedges of each
  motif copy are constrained. If \code{FALSE}, only motif edges are
  constrained.}
  \item{overlap_sizes}{Integer vector of shared-vertex counts. Values should
  be between 0 and \code{nrow(H)}; the default is \code{2:nrow(H)}.}
  \item{mode}{\code{"count"} enumerates copies in a realized binary graph;
  \code{"expectation"} sums probabilities in an independent-edge network.}
  \item{include_incompatible}{Logical; include incompatible overlap patterns
  as zero-valued rows in \code{by_pattern}. Such patterns never contribute
  to the totals.}
  \item{method}{\code{"enumerate"} always uses recursive embedding
  enumeration. \code{"homogeneous"} uses an exact closed form and requires
  all relevant off-diagonal entries of \code{X} to be equal.
  \code{"auto"} selects the homogeneous formula when possible and otherwise
  enumerates embeddings.}
  \item{homogeneous_tol}{Nonnegative numerical tolerance used by
  \code{method = "auto"} to identify a homogeneous matrix.}
}
\value{
A list with \code{total} and \code{labeled_total};
\code{unlabeled_union_total}; \code{unlabeled_pair_total};
\code{motif_automorphism_count}; and \code{by_pattern}, a data frame giving
the overlap size, union size, compatibility, union and motif automorphism
counts, and labeled/unlabeled value for every reported pattern. The fields
\code{mode}, \code{method}, \code{homogeneous_probability},
\code{directed}, \code{induced}, and \code{overlap_sizes} record the
calculation settings.
}
\details{
The result is exact up to floating-point arithmetic. General enumeration
grows rapidly with network size, motif size, and union size; it is intended
for small networks or for validating ISS estimates. The homogeneous path is
usually much faster.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 6, 6); diag(P) <- 0

exact_union <- PairUnion(P, H, mode = "expectation",
                         directed = FALSE, method = "auto")
exact_union$unlabeled_pair_total
exact_union$by_pattern
}
\seealso{\code{\link{PairVar}}, \code{\link{est.pairunion}}}

