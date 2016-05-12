\name{cgrm.D}
\alias{cgrm.D}

\title{Dominance Genomic Relationship Matrix}

\description{
  Based on a marker matrix \eqn{\mathbf{X}} with \{-1,0,1\} - out of which a column-wise centered dominance coefficient matrix will be constructed and a shrinkage parameter \eqn{\lambda}, \code{cgrm.D} returns
  the following dominance genomic relationship matrix according to Su et al. (2012):
                             \deqn{\mathbf{G} = (1-\lambda) \frac{\mathbf{X X}^{'}}{\sum\limits_{i=1}^n 2 p_i q_i(1-2 p_i q_i) }  + \mathbf{I}\lambda}

  The additive marker coefficients will be used to compute dominance coefficients as: \code{1-abs(X)}

}



\usage{
cgrm.D(X, lambda=0)
}
\arguments{
  \item{X}{marker matrix}
  \item{lambda}{numeric scalar, shrinkage parameter}
}

\details{
...
}

\value{
Dominance relationship matrix with dimension nrow(X)
}


\author{
Claas Heuer
}

\references{
Su G, Christensen OF, Ostersen T, Henryon M, Lund MS (2012) "Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers". PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293
}



\seealso{\code{\link{cgrm}}, \code{\link{cgrm.A}}.}

\examples{
# generate random data
rand_data(100,500)

D <- cgrm.D(M,lambda=0.01)
}
\keyword{Genomic Relationship}
