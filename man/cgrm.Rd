\name{cgrm}
\alias{cgrm}

\title{Genomic Relationship Matrices}

\description{
  Based on a coefficient-matrix (i.e. marker matrix) \eqn{\mathbf{X}} that will be scaled column-wise, a weight-vector \eqn{\mathbf{w}} and a shrinkage parameter \eqn{\lambda}, \code{cgrm} returns
  the following similarity matrix: 
                             \deqn{\mathbf{G} = (1-\lambda) \frac{\mathbf{X D X}^{'}}{\sum\mathbf{w}} + \mathbf{I}\lambda}

  where \eqn{\mathbf{D}} = \eqn{diag(\mathbf{w})}. 
  A weighted genomic relationship matrix allows running TA-BLUP as described in Zhang et al. (2010).

}



\usage{
cgrm(X, w = NULL, lambda=0)
}
\arguments{
  \item{X}{coefficient matrix}
  \item{w}{numeric vector of weights for every column in X}
  \item{lambda}{numeric scalar, shrinkage parameter}
}

\details{
...
}

\value{
Similarity matrix with dimension nrow(X)
}


\author{
Claas Heuer
}

\references{
de los Campos, G., Vazquez, A.I., Fernando, R., Klimentidis, Y.C., Sorensen, D., 2013. "Prediction of Complex Human Traits Using the Genomic Best Linear Unbiased Predictor". PLoS Genetics 9, e1003608. doi:10.1371/journal.pgen.1003608

Zhang Z, Liu J, Ding X, Bijma P, de Koning D-J, et al. (2010) "Best Linear Unbiased Prediction of Genomic Breeding Values Using a Trait-Specific Marker-Derived Relationship Matrix". PLoS ONE 5(9): e12648. doi:10.1371/journal.pone.0012648
}


\seealso{\code{\link{cgrm.A}}, \code{\link{cgrm.D}}.}

\examples{
# generate random data
rand_data(100,500)

weights <- (cor(M,y)**2)[,1]

G <- cgrm(M,weights,lambda=0.01)
}
\keyword{Genomic Relationship}
