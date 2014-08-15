\name{get_cor}
\alias{get_cor}
\title{
Compute the prediction accuracy from Cross Validition
}
\description{
Takes a matrix of predictions returned by \code{\link{get_pred}}, a
list of masked phenotypes returned by \code{\link{cCV}} and the original
phenotype vector and returns the correlation between predicted and observed
values
}
\usage{
get_cor(predictions,cv_pheno,y)
}

\arguments{
  \item{predictions}{Prediction matrix returned by \code{\link{get_pred}}}
  \item{cv_pheno}{List of masked phenotypes returned by \code{\link{cCV}}}
  \item{y}{Original unmasked phenotype vector that has been used in \code{\link{cCV}}}
}

\value{
Numeric scalar - Mean prediction accuracy measured as correlation between predicted and observed phenotypes
}

\examples{

# see example of clmm.CV

}

\keyword{Tools}
