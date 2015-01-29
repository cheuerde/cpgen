\name{cpgen-package}
\alias{cpgen-package}
\docType{package}
\title{
cpgen - Parallel genomic evaluations
}
\description{
The package offers a variety of functions that are frequently being used in genomic prediction
and genomewide association studies. The package is based on Rcpp and RcppEigen, hence all routines
are implemented using the matrix algebra library Eigen.
The main emphasis of the package lies in parallel computing which is realized by C++ functions making
use of OpenMP. 
}
\details{
\tabular{ll}{
Package: \tab cpgen\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-07-10\cr
License: \tab License: GPL (>= 2)\cr
}
}

\author{
Claas Heuer

Maintainer: Claas Heuer <cheuer@tierzucht.uni-kiel.de>
}

\references{

Guennebaud, G., Jacob, B., et al.: "Eigen v3". http://eigen.tuxfamily.org (2010)

Dirk Eddelbuettel and Romain Francois (2011). "Rcpp: Seamless R and
C++ Integration". Journal of Statistical Software, 40(8), 1-18. URL
http://www.jstatsoft.org/v40/i08/.

Douglas Bates, Dirk Eddelbuettel (2013). "Fast and Elegant Numerical
Linear Algebra Using the RcppEigen Package". Journal of Statistical
Software, 52(5), 1-24. URL http://www.jstatsoft.org/v52/i05/.




}

\keyword{ package }

