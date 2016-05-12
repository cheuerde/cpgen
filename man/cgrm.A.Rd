\name{cgrm.A}
\alias{cgrm.A}

\title{Additive Genomic Relationship Matrix}

\description{
  Based on a marker matrix \eqn{\mathbf{X}} with \{-1,0,1\} - coding that will be centered column-wise and a shrinkage parameter \eqn{\lambda}, \code{cgrm.A} returns
  the following additive genomic relationship matrix according to VanRaden (2008):
                             \deqn{\mathbf{G} = (1-\lambda) \frac{\mathbf{X X}^{'}}{\sum\limits_{i=1}^n 2 p_i q_i}  + \mathbf{I}\lambda}

}



\usage{
cgrm.A(X, lambda=0, yang=FALSE)
}
\arguments{
  \item{X}{marker matrix}
  \item{lambda}{numeric scalar, shrinkage parameter}
  \item{yang}{boolean, diagonal elements of A according to Yang et al. (2010)}
}

\details{
...
}

\value{
Additive genomic relationship matrix with dimension nrow(X)
}


\author{
Claas Heuer
}

\references{
VanRaden, P.M. "Efficient Methods to Compute Genomic Predictions". Journal of Dairy Science 91, no. 11 (November 2008): 4414-23. doi:10.3168/jds.2007-0980.

Yang, Jian, Beben Benyamin, Brian P McEvoy, Scott Gordon, Anjali K Henders, Dale R Nyholt, Pamela A Madden, et al. "Common SNPs Explain a Large Proportion of the Heritability for Human Height". Nature Genetics 42, no. 7 (July 2010): 565-69. doi:10.1038/ng.608.
}

\seealso{\code{\link{cgrm}}, \code{\link{cgrm.D}}}

\examples{
# generate random data
rand_data(100,500)

### compute the additive genomic relationship matrix
A <- cgrm.A(M,lambda=0.01)
}
\keyword{Genomic Relationship}
