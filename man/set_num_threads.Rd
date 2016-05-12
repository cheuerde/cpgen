\name{set_num_threads}
\alias{set_num_threads}
\title{
Set the number of OpenMP threads used by the functions of package \code{cpgen}
}
\description{
This function sets the value of the global variable stored in \code{options()$cpgen.threads} to the assigned integer.
Note_1: The assigned value may exceed the number of physical cores present but that might lead to dramatical
decrease in performance.
Note_2: The function can override the global variable 'OMP_NUM_THREADS' (if \code{global=TRUE} and hence also other non-cpgen functions
are affected by a call to \code{set_num_threads()}.
}
\usage{
set_num_threads(x,silent=FALSE, global=FALSE)
}

\arguments{
\item{x}{Integer scalar that specifies the number of threads to be used by \code{cpgen}-functions}
\item{silent}{boolean, controls whether to print a message}
\item{global}{boolean, change openmp threads globally (might effect other libraries)}
}

\value{
Changes the global variable \code{cpgen.threads} to the value in \code{x}
}

\seealso{\code{\link{get_num_threads}, \link{get_max_threads}, \link{check_openmp}}}

\examples{

# Control the number of threads being used in an R-session:

# set the number of threads to 1
set_num_threads(1)

#### Use a parallelized cpgen-function

# generate random data
rand_data(100,500)

# check single-threaded performance
system.time(W <- M\%c\%t(M))

# set number of threads to 2

set_num_threads(2)

# check multi-threaded performance
system.time(W <- M\%c\%t(M))
}


\keyword{Parallelization}


