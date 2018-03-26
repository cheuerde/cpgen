\name{cSSBR.setup}
\alias{cSSBR.setup}

\title{Preparing Model terms for Single Step Bayesian Regression}

\description{
This function prepares all model terms for SSBR using pedigree and marker information.
The function is particularly useful for using the reported model terms on multiple 
phenotypes, for cross validation (\code{\link{clmm}}), for genomewide association studies
or to pass them to alternative software.
}


\usage{
cSSBR.setup(data, M, M.id, verbose=TRUE, returnAll = FALSE)
}
\arguments{
\item{data}{\code{data.frame} with four columns: \code{id, sire, dam, y}}
\item{M}{Marker Matrix for genotyped individuals}
\item{M.id}{Vector of length \code{nrow(M)} representing rownames for M}
\item{verbose}{Prints progress to the screen}
\item{returnAll}{Return all individuals present in \code{M} and pedigree. If false, function only return individuals with phenotypes}
}



\details{
...

}

\value{
List of 5:
  \item{ids}{ids for the model (ordered as in other model terms)}
  \item{y}{phenotype vector}
  \item{Marker_Matrix}{Combined Marker Matrix including imputed and genotyped individuals}
  \item{Z_residual}{Design Matrix used to model the residual error for the imputed individuals}
  \item{\code{ginverse_residual}}{Submatrix of the inverse of the numerator relationship matrix. 
                                  Used to model the residual error for the imputed individuals}
  \item{J}{Design Matrix for modelling differences in means for genotyped and non-genotyped individuals}
  \item{imputed}{List of imputed individuals}
}


\author{
Claas Heuer
}

\references{

Fernando, R.L., Dekkers, J.C., Garrick, D.J.: A class of bayesian methods to combine large numbers of
genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution 46(1), 50 (2014)

}

\seealso{\code{\link{cSSBR.setup}, \link{clmm}}}

\examples{

# example dataset

id <- 1:6
sire <- c(rep(NA,3),rep(1,3))
dam <- c(rep(NA,3),2,2,3)

# phenotypes
y <- c(NA, 0.45, 0.87, 1.26, 1.03, 0.67)

dat <- data.frame(id=id,sire=sire,dam=dam,y=y)


# Marker genotypes
M <- rbind(c(1,2,1,1,0,0,1,2,1,0),
           c(2,1,1,1,2,0,1,1,1,1),
           c(0,1,0,0,2,1,2,1,1,1))

M.id <- 1:3

model_terms <- cSSBR.setup(dat,M, M.id)

var_y <- var(y,na.rm=TRUE)
var_e <- (10*var_y / 21)
var_a <- var_e 
var_m <- var_e / 10

# put emphasis on the prior
df = 500

par_random=list(list(method="ridge",scale=var_m,df = df),list(method="ridge",scale=var_a,df=df))

set_num_threads(1)

# passing model terms to 'clmm'
mod<-clmm(y=model_terms$y, 
          Z=list(model_terms$Marker_Matrix,model_terms$Z_residual),
          ginverse = list(NULL, model_terms$ginverse_residual),
          par_random=par_random,
          scale_e = var_e,
          df_e=df,
          niter=50000,
          burnin=30000)
          
# check marker effects
print(round(mod[[4]]$posterior$estimates_mean,digits=2))


}
\keyword{Genomic Prediction}
