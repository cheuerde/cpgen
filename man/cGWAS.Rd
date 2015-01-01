\name{cGWAS}
\alias{cGWAS}

\title{Genomewide Association Study}

\description{
This function runs GWAS for continuous traits. Population structure
that can lead to false positive association signals can be accounted for
by passing a Variance-covariance matrix of the phenotype vector (Kang et al., 2010).
The GLS-solution for fixed effects is computed as:

\deqn{\hat{\boldymbol{\beta}} = (\mathbf{X'V}^{-1}\mathbf{X})^{-1}\mathbf{X'V}^{-1}\mathbf{y}}

Equivalent solutions are obtained by premultiplying the design matrix \eqn{\mathbf{X}} for fixed effects and the phenotype vector \eqn{\mathbf{y}}  by \eqn{\mathbf{V}^{-1/2}} :

\deqn{\hat{\boldymbol{\beta}} = (\mathbf{X}^{\ast\prime}\mathbf{X}^{\ast})^{-1}\mathbf{X}^{\ast\prime}\mathbf{y}^{\ast}} with \deqn{\mathbf{X}^{\ast} =\mathbf{V}^{-1/2}\mathbf{X}}    \deqn{\mathbf{y}^{\ast} =\mathbf{V}^{-1/2}\mathbf{y}}


}



\usage{
cGWAS(y,M,X=NULL,V=NULL,dom=FALSE, verbose=TRUE)
}
\arguments{
\item{y}{vector of phenotypes}
\item{M}{Marker matrix}
\item{X}{Optional Design Matrix for additional fixed effects. If omitted a column-vector of ones will be assigned}
\item{V}{Inverse square root of the Variance-covariance matrix for the phenotype vector of type: \code{matrix} or \code{dgCMatrix}. Used for computing the GLS-solution of fixed effects. If omitted an identity-matrix will be assigned}
\item{dom}{Defines whether to include an additional dominance coefficient for every marker. Note: only useful if the genotype-coding in \code{M} follows \{-1,0,1\}
The dominance coefficient is computed as: \code{1-abs(M)}}
\item{verbose}{prints progress to the screen}
}

\details{
...
}

\value{
List of 3 vectors or matrices. If \code{dom=TRUE} every element of the list will be a matrix with two columns. First column additive, second dominance:
  \item{p-value}{Vector of p-values for every marker}
  \item{beta}{GLS solution for fixed marker effects}
  \item{se}{Standard Errors for values in \code{beta}}
}


\author{
Claas Heuer
}


\references{

Kang, Hyun Min, Jae Hoon Sul, Susan K Service, Noah A Zaitlen, Sit-yee Kong, Nelson B Freimer, Chiara Sabatti, and Eleazar Eskin. "Variance Component Model to Account for Sample Structure in Genome-Wide Association Studies." Nature Genetics 42, no. 4 (April 2010): 348-54. doi:10.1038/ng.548.
}


\seealso{\code{\link{cGWAS.emmax}}}

\examples{
\dontrun{
# generate random data
rand_data(500,5000)


### GWAS without accounting for population structure
mod <- cGWAS(y,M)

### GWAS - accounting for population structure
## Estimate variance covariance matrix of y

G <- cgrm.A(M,lambda=0.01)

fit <- cGBLUP(y,G,verbose=FALSE)

### construct V
V <- G*fit$var_a + diag(length(y))*fit$var_e

### get the inverse square root of V
V2inv <- V \%**\% -0.5

### run GWAS again
mod2 <- cGWAS(y,M,V=V2inv,verbose=TRUE)
}


}
\keyword{GWAS}
