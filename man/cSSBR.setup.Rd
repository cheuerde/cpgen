\name{cSSBR.setup}
\alias{cSSBR.setup}

\title{Preparing Model terms for Single Step Bayesian Regression}

\description{
This function prepares all model terms for SSBR using pedigree and marker information.
The function is particularly useful for using the reported model terms on multiple 
phenotypes or for cross validation (\code{\link{clmm.CV}}).
}


\usage{
cSSBR.setup(data, M, verbose=TRUE)
}
\arguments{
\item{data}{\code{data.frame} with four columns: id, sire, dam, y}
\item{M}{Marker Matrix for genotyped individuals. MUST contain rownames}
\item{verbose}{Prints progress to the screen}
}



\details{
yet to be written ...

}

\value{
List of 4:
  \item{ids}{ids for the model (ordered as in other model terms)}
  \item{y}{phenotype vector}
  \item{Marker_Matrix}{Combined Marker Matrix including imputed and genotyped individuals}
  \item{Z_residual}{Design Matrix used to model the residual error for the imputed individuals}
}


\author{
Claas Heuer
}

\references{

Rohan L Fernando, Jack C M Dekkers and Dorian J Garrick. "A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole genome
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

rownames(M) <- 1:3

model_terms <- cSSBR.setup(dat,M)

var_y <- var(y,na.rm=TRUE)
var_e <- (10*var_y / 21)
var_a <- var_e 
var_m <- var_e / 10

# put emphasis on the prior (reproducing results given in SSBR-paper)
df = 500

par_random=list(list(method="random",scale=var_m,df = df),list(method="random",scale=var_a,df=df))

set_num_threads(1)

# passing model terms to 'clmm'
mod<-clmm(y=model_terms$y, 
          random=list(model_terms$Marker_Matrix,model_terms$Z_residual),
          par_random=par_random,
          scale_e = var_e,
          df_e=df,
          niter=50000,
          burnin=30000)
          
# check marker effects
print(round(mod[[4]]$posterior$estimates_mean,digits=2))


}
\keyword{Genomic Prediction}
