\name{sec.moment}
\alias{sec.moment}
\alias{print.sec.moment}
\title{Estimate the Second Moment of a Subgraph Count}
\description{
Runs \code{est.variance()} and adds labeled and unlabeled second moments using
\eqn{E(N^2)=Var(N)+E(N)^2}.
}
\usage{sec.moment(...)}
\arguments{
  \item{...}{All arguments accepted by \code{est.variance()}, including
  \code{P}, \code{H}, \code{B}, \code{B_mean}, direction and inducedness,
  overlap selection, mean inputs, allocation, seeds, and inference controls.}
}
\value{
The full \code{est.variance()} result with class \code{"sec.moment"} added
and a \code{second_moment} list. That list contains \code{labeled},
\code{labeled_se}, \code{unlabeled}, and \code{unlabeled_se}.
}
\examples{
H <- matrix(1, 3, 3); diag(H) <- 0
P <- matrix(0.5, 5, 5); diag(P) <- 0
ans <- sec.moment(P = P, H = H, B = 200, B_mean = 200,
                  directed = FALSE, seed = 5,
                  verbose = FALSE, verbose.B = FALSE)
ans$second_moment
}
\seealso{\code{\link{est.variance}}}
