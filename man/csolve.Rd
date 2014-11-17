\name{csolve}
\alias{csolve}
\title{
csolve
}
\description{
This is a wrapper for the Cholesky-solvers 'LLT' (dense case) or 'Simplicial-LLT' (sparse case) from Eigen.
The function computes the solution: \deqn{\mathbf{b} = \mathbf{X}^{-1} \mathbf{y}}
If no vector \code{y} is passed, an identity matrix will be assigned
and the functions returns the inverse of \eqn{\mathbf{X}}.
In the case of multiple right hand sides (as is the case when computing an inverse matrix)
multiple threads will solve equal parts of it.
}
\usage{
csolve(X,y=NULL)	
}

\arguments{
  \item{X}{positive definite square matrix of type \code{matrix} or \code{dgCMatrix}}
  \item{y}{numeric vector of length equal to columns/rows of \code{X}}
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


