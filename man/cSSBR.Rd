\name{cSSBR}
\alias{cSSBR}

\title{Single Step Bayesian Regression}

\description{
This function runs Single Step Bayesian Regression (SSBR) for the prediction of breeding values
in a unified model that incorporates genotyped and ungenotyped individuals (Fernando et al., 2014).
}


\usage{
cSSBR(data, M, M.id, X=NULL, par_random=NULL, scale_e=0, df_e=0, 
      niter=5000, burnin=2500, seed=NULL, verbose=TRUE)
}
\arguments{
\item{data}{\code{data.frame} with four columns: id, sire, dam, y}
\item{M}{Marker Matrix for genotyped individuals. MUST contain rownames}
\item{M.id}{Vector of length \code{nrow(M)} representing rownames for M}
\item{X}{Fixed effects design matrix of type: \code{matrix} or \code{dgCMatrix}.
         If omitted a column-vector of ones will be assigned. Must have as many rows as data}
\item{par_random}{as in \code{\link{clmm}}}
\item{niter}{as in \code{\link{clmm}}}
\item{burnin}{as in \code{\link{clmm}}}
\item{verbose}{as in \code{\link{clmm}}}
\item{scale_e}{as in \code{\link{clmm}}}
\item{df_e}{as in \code{\link{clmm}}}
\item{seed}{as in \code{\link{clmm}}}
}



\details{
yet to be written ...

}

\value{
List of 3 + number of random effects as in \code{\link{clmm}} +
  \item{SSBR}{List of 4:
     \itemize{
       \item{ids - ids used in the model (ordered as in other model terms)}
       \item{y - phenotype vector}
       \item{Marker_Matrix - Combined Marker Matrix including imputed and genotyped individuals}
       \item{Z_residual - Design Matrix used to model the residual error for the imputed individuals}
       \item{Breeding_Values - Predicted Breeding Values for all animals in \code{data} that have genotypes and/or phenotypes}
      }
   }
}


\author{
Claas Heuer
}

\references{

Rohan L. Fernando, Jack C. M. Dekkers and Dorian J. Garrick. "A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole genome
analyses". 2014

}

\seealso{\code{\link{cSSBR.setup}, \link{clmm}}}

\examples{

# this is the example dataset given in Fernando et al., 2014

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

var_y <- var(y,na.rm=TRUE)
var_e <- (10*var_y / 21)
var_a <- var_e 
var_m <- var_e / 10

# put emphasis on the prior (reproducing results given in SSBR-paper)
df = 500

par_random=list(list(method="random",scale=var_m,df = df),list(method="random",scale=var_a,df=df))

set_num_threads(1)
mod<-cSSBR(data = dat,M=M, M.id=M.id, par_random=par_random, scale_e = var_e, df_e=df,niter=50000,burnin=30000)

# check marker effects
print(round(mod[[4]]$posterior$estimates_mean,digits=2))

# check breeding value prediction:
print(round(mod$SSBR$Breeding_Values,digits=2))

}
\keyword{Genomic Prediction}
