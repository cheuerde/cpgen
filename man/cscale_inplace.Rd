\name{cscale_inplace}
\alias{cscale_inplace}
\title{
cscale_inplace
}
\description{
Center (and scale) a matrix 'inplace'. The function is meant for
big matrices that shall be scaled inplace, hence without creating a copy

}
\usage{
cscale_inplace(X, means = NULL, vars = NULL, scale=FALSE)	
}

\arguments{
  \item{X}{numeric matrix}
  \item{means}{numeric vector, if ommitted will be computed using code{ccolmv}}
  \item{vars}{numeric vector, if ommitted will be computed using \code{ccolmv}}
  \item{scale}{boolean - scale the matrix}
}

\value{
nothing, function works 'inplace'
}


\examples{
\dontrun{
# generate random data
rand_data(500,5000)

# scale matrix
cscale_inplace(M,scale=TRUE)
}
}

\keyword{Tools}
