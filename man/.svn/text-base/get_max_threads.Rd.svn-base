\name{get_max_threads}
\alias{get_max_threads}
\title{
Get the maximum number of threads available
}
\description{
This is a wrapper for the OpenMP-function \code{omp_get_max_threads()}, hence the function will
report the result of the according omp-function.
Note: The returned value does not necessarily reflect the number of physical cores present but in most
cases it will.
}
\usage{
get_max_threads()
}


\value{
Returns the value reported by \code{omp_get_max_threads()}
}

\seealso{\code{\link{set_num_threads}, \link{get_num_threads}, \link{check_openmp}}}

\examples{

# set number of threads to the value reported by get_max_threads()
n_threads <- get_max_threads()
set_num_threads(n_threads)

# check
get_num_threads()


}

\keyword{Parallelization}


