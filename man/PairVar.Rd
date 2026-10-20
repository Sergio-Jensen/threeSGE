\name{PairVar}
\alias{PairVar}
\title{Exactly Enumerate the Variance of a Motif Count}
\description{
Exactly computes motif-count variance by enumerating every requested overlap
pattern and subtracting marginal-product mass (IPUM) from compatible joint
mass (FPUG). Unlike \code{est.variance()}, no importance sampling is used.
}
\usage{
PairVar(P, H, directed = FALSE, induced = TRUE,
        overlap_sizes = 2:nrow(H), include_incompatible = TRUE,
        method = c("auto", "enumerate", "homogeneous"),
        homogeneous_tol = 1e-12)
}
\arguments{
  \item{P}{Square independent-edge probability matrix with every entry in
  \eqn{[0,1]}. For undirected analyses, supply a symmetric matrix with zero
  diagonal. Raw weights are not normalized by this function.}
  \item{H}{Square motif adjacency matrix; nonzero entries denote edges.}
  \item{directed}{Logical; whether ordered directed edges are used.}
  \item{induced}{Logical; if \code{TRUE}, impose motif-edge and motif-nonedge
  constraints. If \code{FALSE}, impose edge constraints only.}
  \item{overlap_sizes}{Integer vector of shared-vertex counts between 0 and
  \code{nrow(H)}. The default is \code{2:nrow(H)}.}
  \item{include_incompatible}{Logical; include incompatible patterns in
  \code{by_pattern}. Regardless of this display choice, their zero joint
  mass and nonzero marginal-product mass remain part of the variance total.}
  \item{method}{\code{"enumerate"} recursively enumerates embeddings;
  \code{"homogeneous"} uses the exact homogeneous closed form; and
  \code{"auto"} chooses the closed form when all relevant off-diagonal
  probabilities agree within \code{homogeneous_tol}.}
  \item{homogeneous_tol}{Nonnegative tolerance for homogeneous-probability
  detection.}
}
\value{
A list containing \code{labeled_variance} and
\code{unlabeled_variance}; \code{labeled_joint_total};
\code{labeled_marginal_product_total}; \code{motif_automorphism_count}; and
\code{by_pattern}. The pattern table gives overlap size, union size,
compatibility, joint sum, marginal-product sum, and labeled/unlabeled
covariance contributions. \code{method}, \code{homogeneous_probability},
\code{directed}, \code{induced}, and \code{overlap_sizes} record the
settings.
}
\details{
The result is exact up to floating-point arithmetic. Recursive enumeration
can be expensive and is primarily useful for small networks, theoretical
benchmarks, and validation of VISS. Homogeneous networks use a much faster
closed-form calculation.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 6, 6); diag(P) <- 0

exact_var <- PairVar(P, H, directed = FALSE, method = "auto")
exact_var$unlabeled_variance  # 5
exact_var$by_pattern
}
\seealso{\code{\link{PairUnion}}, \code{\link{est.variance}}}
