\name{csolve}
\alias{csolve}
\title{
csolve
}
\description{
This is a wrapper for the Cholesky-solver 'llt' from Eigen.
The function computes the solution: \deqn{\mathbf{b} = \mathbf{X}^{-1} \mathbf{y}}
If no vector \code{y} is passed, an identity matrix will be assigned
and the functions returns the inverse of \eqn{\mathbf{X}}.
}
\usage{
csolve(X,y=NULL)	
}

\arguments{
  \item{X}{positive definite square matrix}
  \item{y}{vector of length equal to columns/rows of \code{X}}
}

\value{
Solution vector/matrix
}

\examples{

# Least Squares Solving

# Generate random data 

n = 1000
p = 500

M <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)

# least squares solution:

b <- csolve(t(M) \%c\% M, t(M) \%c\% y)


}

\keyword{Tools}


