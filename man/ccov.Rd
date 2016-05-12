\name{ccov}
\alias{ccov}
\title{
ccov
}
\description{
Computation of covariance- or correlation-matrix. 
Shrinkage estimate through the use of 'lambda'. Weights for observations
can be passed.

}
\usage{
ccov(X,lambda=0, w=NULL, compute_cor=FALSE)	
}

\arguments{
  \item{X}{matrix}
  \item{lambda}{numeric scalar, shrinkage parameter}
  \item{w}{numeric vector of weights with same lengths as rows in X}
  \item{compute_cor}{boolean - defines whether the functions returns a correlation- rather than a covariance matrix}
}

\value{
Covariance matrix with dimension ncol(X)
}


\examples{
# generate random data
rand_data(100,500)

# compute correlation matrix of t(M)
corM <- ccov(t(M),compute_cor=TRUE)
}

\keyword{Tools}
