\name{\%c\%}
\alias{\%c\%}
\title{
(Parallel) Matrix product operator
}
\description{
This operator computes the crossproduct between two matrices.
It can be used as a replacement for \code{\%*\%} in many cases. The operator
only accepts matrices of types: \code{matrix} or \code{dgCMatrix}. 
In the case of two dense matrices the operator will compute the crossproduct 
in parallel (Eigen + OpenMP)
}
\usage{
X\%c\%Y	
}

\arguments{
\item{X}{Matrix or vector (treated as column-vector) of type: \code{matrix} or \code{dgCMatrix}}
\item{Y}{as \code{X}}
}

\value{
Matrix of type: \code{matrix} or \code{dgCMatrix}
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


