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
\item{M}{Marker Matrix for genotyped individuals}
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
The function sets up the following model using \code{\link{cSSBR.setup}}:
\deqn{
\mathbf{y} = \mathbf{Xb} + \mathbf{M\alpha} + \mathbf{Z\epsilon} + \mathbf{e}
}
The matrix \eqn{\mathbf{M}} denotes a combined marker matrix consisting of actual and imputed marker covariates.
Best linear predictions of gene content (Gengler et al., 2007) for the non-genotyped individuals are obtained using: \eqn{\mathbf{A}^{11}\hat{\mathbf{M}_1} = -\mathbf{A}^{12}\mathbf{M}_2} (Fernando et al., 2014).
\eqn{\mathbf{A}^{11}} and \eqn{\mathbf{A}^{12}} are submatrices of the inverse of the numerator relationship matrix, which is easily obtained (Henderson, 1976). The subscripts 1 and 2 denote non genotyped and genotyped individuals respectively. The very sparse equation system is being solved using a sparse cholesky solver provided by the Eigen library.
The residual imputation error has variance: \eqn{(\mathbf{A}^{11})^{-1}\sigma_{\epsilon}^2} and is modelled by constructing the design matrix as \eqn{\mathbf{Z} = \mathbf{L}_{11}}, where \eqn{\mathbf{A}=\mathbf{LL}'} and \eqn{(\mathbf{A}^{11})^{-1} = \mathbf{L}_{11}\mathbf{L'}_{11}}.
}

\value{
List of 3 + number of random effects as in \code{\link{clmm}} +
  \item{SSBR}{List of 4:
     \itemize{
       \item{ids - ids used in the model (ordered as in other model terms)}
       \item{y - phenotype vector}
       \item{X - Design matrix for fixed effects}
       \item{Marker_Matrix - Combined Marker Matrix including imputed and genotyped individuals}
       \item{Z_residual - Design Matrix used to model the residual error for the imputed individuals}
       \item{Breeding_Values - Predicted Breeding Values for all animals in \code{data} that have
             genotypes and/or phenotypes}
      }
   }
}


\author{
Claas Heuer
}

\references{

Fernando, R.L., Dekkers, J.C., Garrick, D.J.: A class of bayesian methods to combine large numbers of
genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution 46(1), 50 (2014)

Gengler, N., Mayeres, P., Szydlowski, M.: A simple method to approximate gene content in large pedigree
populations: application to the myostatin gene in dual-purpose belgian blue cattle. animal 1(01), 21 (2007)

Henderson, C.R.: A simple method for computing the inverse of a numerator relationship matrix used in
prediction of breeding values. Biometrics 32(1), 69–83 (1976)

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
mod<-cSSBR(data = dat,
           M=M,
           M.id=M.id,
           par_random=par_random,
           scale_e = var_e,
           df_e=df,
           niter=50000,
           burnin=30000)

# check marker effects
print(round(mod[[4]]$posterior$estimates_mean,digits=2))

# check breeding value prediction:
print(round(mod$SSBR$Breeding_Values,digits=2))

}
\keyword{Genomic Prediction}
