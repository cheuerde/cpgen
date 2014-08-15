\name{get_pred}
\alias{get_pred}
\title{
Extract predictions vectors of an object returned by \code{\link{clmm.CV}}
}
\description{
Takes an object returned by \code{\link{clmm.CV}} and returns a matrix
of predicted values from every model. Every columns represents the prediction vector of one model
}
\usage{
get_pred(mod)
}

\arguments{
  \item{mod}{List returned by \code{\link{clmm.CV}}}
}

\value{
Matrix of prediction vectors in columns
}

\examples{

# see example of clmm.CV

}

\keyword{Tools}
