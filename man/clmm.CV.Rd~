\name{clmm.CV}
\alias{clmm.CV}

\title{Cross Validation with Linear Mixed Models using Gibbs Sampling}

\description{
This function allows running Cross Validation using \code{\link{clmm}}.
All the formulation is equal to \code{\link{clmm}} except for the phenotypes, which
are being passed as a list of equally sized vectors.
The main advantage of the function is that several threads can access the very same 
data once assigned, which means that the design matrices only have to be allocated once.
The parallel scaling of this function is almost linear.

}


\usage{
clmm.CV(y, X = NULL , random = NULL, par_random = NULL, niter=10000,
burnin=5000,scale_e=0,df_e=-2, verbose = TRUE, seed = NULL)
}
\arguments{
\item{y}{list of phenotype vectors}
\item{X}{Fixed effects design matrix of type: \code{matrix} or \code{dgCMatrix}. If omitted a column-vector of ones will be assigned }
\item{random}{list of design matrices for random effects - every element of the list represents one random effect and may be of type: \code{matrix} or \code{dgCMatrix}}
\item{par_random}{list of options for random effects. If passed, the list must have as many elements as \code{random}. Every element must be a list of 3:
\itemize{
\item{scale - scale parameter for the inverse chi-square prior}
\item{df - degrees of freedom for the inverse chi-square prior}
\item{method - method to be used for the random effects, may be: \code{random} or \code{BayesA}}
}
}
\item{niter}{number of iterations}
\item{burnin}{number of iterations to be discarded as burnin}
\item{verbose}{prints progress to the screen}
\item{scale_e}{scale parameter for the inverse chi-square prior for the residuals}
\item{df_e}{degrees of freedom for the inverse chi-square prior for the residuals}
\item{seed}{seed for the random number generator. If omitted, a seed will be generated based on machine and time}
}



\details{
In C++:
For every element of the phenotype list a new instance of an MCMC-object will be created.
All the memory allocation needed for running the model is done by the major thread. 
The function then iterates over all objects and runs the gibbs sampler. This step is parallelized,
which means that as many models are being run at the same time as threads available.
All MCMC-objects are totally independent from each other, they only share the same design-matrices.
Every object has its own random-number generator with its own seed which allows perfectly reproducible
results. 

}

\value{
List of \code{length(Y)} with elements equal to the output of \code{\link{clmm}}
}


\author{
Claas Heuer
}


\seealso{\code{\link{clmm}}}

\examples{

### Running a 4-fold cross-validation with one repetition:
\dontrun{

# generate random data
rand_data(500,5000)

### compute the list of masked phenotype-vectors for CV
y_CV <- cCV(y,fold=4,reps=1)


### Cross Validation using GBLUP
G.A <- cgrm.A(M,lambda=0.01)


### generate the list of design matrices for clmm
random = list(t(chol(G.A)))

### specify options
par_random = list(list(method="random",scale=var(y)/2,df=5))

### run 
fit <- clmm.CV(y_CV,random=random,par_random=par_random,niter=5000,burnin=2500)

### inspect results
str(fit)

### obtain predictions
pred <- get_pred(fit)

### prediction accuracy
get_cor(pred,y_CV,y)
}

}
\keyword{Genomic Prediction}
