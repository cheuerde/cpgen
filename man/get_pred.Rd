\name{get_pred}
\alias{get_pred}
\title{
Extract predictions vectors of an object returned by \code{\link{clmm}} using multpile phenotypes
}
\description{
Takes an object returned by \code{\link{clmm}} using multpile phenotypes and returns a matrix
of predicted values from every model. Every columns represents the prediction vector of one model
}
\usage{
get_pred(mod)
}

\arguments{
  \item{mod}{List returned by \code{\link{clmm}} using multpile phenotypes}
}

\value{
Matrix of prediction vectors in columns
}

\seealso{\code{\link{clmm}, \link{get_cor}, \link{cCV}}}

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
h2 = 0.3
scale = unlist(lapply(y_CV,function(x)var(x,na.rm=T))) * h2
df = rep(5,length(y_CV))
par_random = list(list(method="ridge",scale=scale,df=df))

### run 
fit <- clmm(y_CV,random=random,par_random=par_random,niter=5000,burnin=2500)

### inspect results
str(fit)

### obtain predictions
pred <- get_pred(fit)

### prediction accuracy
get_cor(pred,y_CV,y)
}

}
\keyword{Genomic Prediction}
\keyword{Tools}

