\name{cGBLUP}
\alias{cGBLUP}

\title{Genomic BLUP}

\description{
This function allows fitting a mixed model with one random effect besides the residual 
using \code{\link{clmm}}. The random effect \eqn{\mathbf{a}} follows some
covariance-structure \eqn{\mathbf{G}}.

}



\usage{
cGBLUP(y,G,X=NULL, scale_a = 0, df_a = -2, scale_e = 0, df_e = -2,
          niter = 10000, burnin = 5000, seed = NULL, verbose=TRUE)
}
\arguments{
\item{y}{vector of phenotypes}
\item{G}{Relationship matrix / covariance structure for random effects}
\item{X}{Optional Design Matrix for fixed effects. If omitted a column-vector of ones will be assigned}
\item{scale_a}{prior scale parameter for \eqn{a}}
\item{df_a}{prior degrees of freedom for \eqn{a}}
\item{scale_e}{prior scale parameter for \eqn{e}}
\item{df_e}{prior degrees of freedom for \eqn{e}}
\item{niter}{Number of iterations}
\item{burnin}{Burnin}
\item{seed}{Seed}
\item{verbose}{Prints progress to the screen}
}


\details{
Kang et al. (2008):
\deqn{
\mathbf{y} = \mathbf{Xb} + \mathbf{a} + \mathbf{e}  \textrm{   with:   } \mathbf{a} \sim MVN(\mathbf{0},\mathbf{G}\sigma^2_a)}

By finding the decomposition: \eqn{\mathbf{G = UDU'}} and premultiplying the model equation by \eqn{\mathbf{U'}} we get:
\deqn{\mathbf{U'y = U'Xb + U'a + U'e}}
with:
\deqn{Var(\mathbf{U'y}) = \mathbf{U'G'U} \sigma^2_a + \mathbf{U'U} \sigma^2_e }
\deqn{\mathbf{U'UDU'U}\sigma^2_a + \mathbf{I}\sigma^2_e}
\deqn{\mathbf{D}\sigma^2_a + \mathbf{I}\sigma^2_e}

After diagonalization of the variance-covariance structure the transformed model is being fitted by passing \eqn{\mathbf{D}^{1/2}}
as the design matrix for the random effects to \code{\link{clmm}}.
The results are subsequently backtransformed and returned by the function.

}

\value{
List of 6:
  \item{var_e}{Posterior mean of the residual variance}
  \item{var_a}{Posterior mean of the random-effect variance}
  \item{b}{Posterior means of the fixed effects}
  \item{a}{Posterior means of the random effects}
  \item{posterior_var_e}{Posterior of the residual variance}
  \item{posterior_var_u}{Posterior of the random variance}
}


\author{
Claas Heuer
}

\references{

Kang, H. M., N. A. Zaitlen, C. M. Wade, A. Kirby, D. Heckerman, M. J. Daly, and E. Eskin. "Efficient Control of Population Structure in Model Organism Association Mapping." Genetics 178, no. 3 (February 1, 2008): 1709-23. doi:10.1534/genetics.107.080101.

}


\seealso{\code{\link{clmm}}, \code{\link{cgrm}}, \code{\link{cGWAS.emmax}}}

\examples{
\dontrun{
# generate random data
rand_data(500,5000)

# compute a genomic relationship-matrix
G <- cgrm(M,lambda=0.01)

# run model
mod <- cGBLUP(y,G)
}
}
\keyword{Genomic Prediction}
