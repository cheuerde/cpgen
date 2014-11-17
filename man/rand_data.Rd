\name{rand_data}
\alias{rand_data}
\title{
Generate random data for test purposes
}
\description{
Generates a random marker-matrix in \{-1,0,1\} coding and a phenotype vector.
Phenotypic variance times h2 (variance explained by markers) is equally spreaded among
all markers (sampled from uniform distribution).
}
\usage{
rand_data(n=500,p_marker=10000,h2=0.3,prop_qtl=0.01,seed=NULL)
}

\arguments{
  \item{n}{Number of oberservations}
  \item{p_marker}{Number of markers}
  \item{h2}{Heritability of the trait}
  \item{prop_qtl}{Proportion of QTL of total number of markers}
  \item{seed}{Seed for RNG}
}

\value{
No return value. Generates two objects globally (\code{M} and \code{y}) that can be used after the execution of the function.
\code{M} is the marker matrix and \code{y} the phenotype vector
}

\examples{

# Generate random data with 100 observations and 500 markers
rand_data(100,500)

# check that objects have been created
str(M)
str(y)

}

\keyword{Genomic Prediction}
