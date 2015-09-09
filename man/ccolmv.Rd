\name{ccolmv}
\alias{ccolmv}
\title{
Colwise means or variances
}
\description{
Computes the colwise means or variances of a matrix - internal use
}
\usage{
ccolmv(X,compute_var=FALSE)
}

\arguments{
  \item{X}{matrix of type: \code{matrix} or \code{dgCMatrix}}
  \item{compute_var}{boolean, defines whether the colwise variances rather than the means will be returned}
}

\value{
Numeric Vector of colwise means or variances of X
}

\examples{

X <- matrix(rnorm(1000*500),1000,500)
means <- ccolmv(X)
vars <- ccolmv(X,compute_var=TRUE)

}

\keyword{Tools}
