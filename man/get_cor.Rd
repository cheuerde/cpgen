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

\seealso{\code{\link{clmm}, \link{get_pred}, \link{cCV}}}

\examples{
### Running a 4-fold cross-validation with one repetition:

# generate random data
rand_data(100,500)

### compute the list of masked phenotype-vectors for CV
y_CV <- cCV(y,fold=4,reps=1)


### Cross Validation using GBLUP
G.A <- cgrm.A(M,lambda=0.01)


### generate the list of design matrices for clmm
Z_list = list(t(chol(G.A)))

### specify options
h2 = 0.3
scale = unlist(lapply(y_CV,function(x)var(x,na.rm=TRUE))) * h2
df = rep(5,length(y_CV))
par_random = list(list(method="ridge",scale=scale,df=df))

### run 
fit <- clmm(y_CV, Z=Z_list, par_random=par_random, niter=5000, burnin=2500)

### inspect results
str(fit)

### obtain predictions
pred <- get_pred(fit)

### prediction accuracy
get_cor(pred,y_CV,y)
}

\keyword{Genomic Prediction}
\keyword{Tools}
