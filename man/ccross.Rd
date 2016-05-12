\name{ccross}
\alias{ccross}
\title{
ccross
}
\description{
Computation of the following matrix-product: \eqn{\mathbf{XDX}'}
Where \eqn{\mathbf{D}} is a diagonal matrix, which is being passed to the function as a vector.
}
\usage{
ccross(X,D=NULL)	
}

\arguments{
  \item{X}{matrix}
  \item{D}{numeric vector, will be used as a weighting diagonal matrix of dimension ncol(X). 
If omitted an identity matrix will be assigned.}
}

\value{
Square matrix of dimension nrow(X)
}

\examples{

# Computing the matrix-square-root of a positive definite square matrix:
# generate random data
rand_data(100,500)

W <- ccross(M)

# this is the implementation of the matrix power-operator '\%**\%'
W_sqrt <- with(eigen(W), ccross(vectors,values**0.5))
}

\keyword{Tools}


