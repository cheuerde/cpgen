\name{cCV}
\alias{cCV}
\title{
Generate phenotype vectors for cross validation
}
\description{
This function takes a phenotype vector and generates \code{folds * reps} masked 
vectors for cross validation. Every vector has as many additional missing values as \code{length(y) / folds}.

}
\usage{
cCV(y,folds=5,reps=1,matrix=FALSE,seed=NULL)
}

\arguments{
  \item{y}{vector of phenotypes - may already contain missing values}
  \item{folds}{integer, number of folds}
  \item{reps}{integer, number of replications}
  \item{matrix}{boolean, if TRUE function returns a matrix rather than a list}
  \item{seed}{numeric scalar, seed for \code{sample}}
}

\value{
List (matrix) with as many items (columns) as \code{folds * reps}
}

\seealso{\code{\link{clmm}, \link{get_pred}, \link{get_cor}}}

\examples{
# generate random data
rand_data(500,1000)

y_CV <- cCV(y,folds=5,reps=20)
}

\keyword{Genomic Prediction}
