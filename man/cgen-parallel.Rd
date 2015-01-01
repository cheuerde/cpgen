\name{Parallelization}
\alias{cpgen-parallel}
\docType{package}
\title{
Multithreading using \code{cpgen}
}
\description{
The package \code{cpgen} makes use of shared memory multi-threading using
OpenMP. R is of single-threaded nature, hence almost the entire package is written
in C++. The package offers a variety of functions that lets you control and check
the number of threads that are being used by the functions of the package.
Internally every function uses the global variable \code{cpgen.threads} which is stored in
\code{options()$cpgen.threads}.
The value can be changed using the function \code{set_num_threads()}. When the package is loaded
in an R-session \code{cpgen.threads} will be set to the value returned by \code{get_max_threads()} which
is a wrapper for the OpenMP-header function \code{omp_get_max_threads()}}

\details{
The following functions are multithreaded and access the variable \code{cpgen.threads}:
\itemize{
\item{cGWAS}
\item{cGWAS.emmax}
\item{clmm}
\item{cGBLUP}
\item{ccross}
\item{\%c\%}
\item{cgrm}
\item{cgrm.A}
\item{cgrm.D}
\item{ccov}
\item{csolve}
\item{cSSBR.setup}
\item{cSSBR}
}}

\seealso{\code{\link{set_num_threads}, \link{get_num_threads}, \link{get_max_threads}, \link{check_openmp}}}

\keyword{Parallelization}

