\name{threeSGE-package}
\alias{threeSGE}
\alias{threeSGE-package}
\docType{package}
\title{Induced Sequential Sampling and Three-Stage Graph Estimation}
\description{
The \pkg{threeSGE} package provides unified interfaces for binary and
weighted induced-subgraph estimation, pairwise-union and variance
estimation, and motif verification under parametric or calibrated null
networks. Exact enumeration is available for small networks when Monte Carlo
estimation is not desired.
}
\details{
The main workflow is:
\enumerate{
  \item Use \code{est.subgraph()} to estimate an observed count or a
  probability-network expectation.
  \item Use \code{est.variance()} when the sampling variance of a subgraph
  count is required.
  \item Use \code{tsge()} for an analytical motif test, or \code{tsge.cal()}
  for random-network calibration.
  \item Use \code{PairUnion()} or \code{PairVar()} for exact enumeration on
  computationally manageable networks.
}
For undirected analyses, symmetric input matrices with a zero diagonal are
recommended. Nonzero entries of the motif matrix denote required edges.
}
\seealso{
\code{\link{est.subgraph}}, \code{\link{est.variance}},
\code{\link{PairUnion}}, \code{\link{PairVar}}, \code{\link{tsge}},
\code{\link{tsge.cal}}
}
\author{Sergio Developer}
\keyword{package}
