\name{cGWAS.emmax}
\alias{cGWAS.emmax}

\title{Genomewide Association Study - EMMAX}

\description{
This is a convenience function that uses the function \code{\link{cGWAS}} but 
estimates the variance-covariance matrix of the phenotype vector in advance 
using \code{\link{clmm}}. This method was termed EMMAX (Kang et al., 2010).

}



\usage{
cGWAS.emmax(y,M,A=NULL,X=NULL,dom=FALSE,verbose=TRUE,scale_a = 0, df_a = -2,
		   scale_e = 0, df_e = -2,niter=15000,burnin=7500,seed=NULL)
}
\arguments{
\item{y}{vector of phenotypes}
\item{M}{Marker matrix}
\item{A}{Relationship matrix that is being used to estimate \eqn{V} - if omitted, \code{A} will be constructed using \code{M} and \code{\link{cgrm}}}
\item{X}{Optional Design Matrix for additional fixed effects. If omitted a column-vector of ones will be assigned}
\item{dom}{Defines whether to include an additional dominance coefficient for every marker. Note: only useful if the genotype-coding in \code{M} follows \{-1,0,1\}
The dominance coefficient is computed as: \code{1-abs(M)}}
\item{verbose}{Prints progress to the screen}
\item{scale_a}{prior scale parameter for \eqn{a}}
\item{df_a}{prior degrees of freedom for \eqn{a}}
\item{scale_e}{prior scale parameter for \eqn{e}}
\item{df_e}{prior degrees of freedom for \eqn{e}}
\item{niter}{Number of iterations}
\item{burnin}{Burnin}
\item{seed}{Seed}
}


\details{
...
}

\value{
List of 3 vectors or matrices. If \code{dom=TRUE} every element of the list will be a matrix with two columns. First column additive, second dominance:
  \item{p-value}{Vector of p-values for every marker}
  \item{beta}{GLS solution for fixed marker effects}
  \item{se}{Standard Errors for values in \code{beta}}
  \item{marker_variance}{Estimate of the marker variance reported by \code{\link{clmm}}}
  \item{residual_variance}{Estimate of the residual variance reported by \code{\link{clmm}}}
}


\author{
Claas Heuer
}

\references{

Kang, H. M., N. A. Zaitlen, C. M. Wade, A. Kirby, D. Heckerman, M. J. Daly, and E. Eskin. "Efficient Control of Population Structure in Model Organism Association Mapping." Genetics 178, no. 3 (February 1, 2008): 1709-23. doi:10.1534/genetics.107.080101.

Kang, Hyun Min, Jae Hoon Sul, Susan K Service, Noah A Zaitlen, Sit-yee Kong, Nelson B Freimer, Chiara Sabatti, and Eleazar Eskin. "Variance Component Model to Account for Sample Structure in Genome-Wide Association Studies." Nature Genetics 42, no. 4 (April 2010): 348-54. doi:10.1038/ng.548.

}

\seealso{\code{\link{cGWAS}}}

\examples{
\dontrun{
# generate random data
rand_data(500,5000)

# run EMMAX
res <- cGWAS.emmax(y,M,verbose=TRUE)
}
}
\keyword{GWAS}






