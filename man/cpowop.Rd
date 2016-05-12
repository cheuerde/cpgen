\name{\%**\%}
\alias{\%**\%}
\title{
Square matrix power operator
}
\description{
This operator computes an arbitrary power of a positive definite square matrix
using an Eigen-decomposition: \eqn{\mathbf{X}^p = \mathbf{UD}^{p}\mathbf{U}'}
}
\usage{
X \%**\% power
}

\arguments{
\item{X}{Positive definite square matrix}
\item{power}{numeric scalar - desired power of X}
}

\value{
Matrix \code{X} to the power \code{p}
}

\examples{
# Inverse Square Root of a positive definite square matrix
X <- matrix(rnorm(100*500),100,500)

XX <- ccross(X)

XX_InvSqrt <- XX \%**\% -0.5

# check result: ((XX')^-0.5 (XX')^-0.5)^-1 = XX'
table(round(csolve(XX_InvSqrt \%c\% XX_InvSqrt),digits=2) == round(XX,digits=2) )
}

\keyword{Tools}


