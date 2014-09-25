\name{clmm}
\alias{clmm}

\title{Linear Mixed Models using Gibbs Sampling}

\description{
This function runs linear mixed models of the following form:
\deqn{
\mathbf{y} = \mathbf{Xb} + \mathbf{Z}_{1}\mathbf{u}_1 + \mathbf{Z}_{2}\mathbf{u}_2 + \mathbf{Z}_{3}\mathbf{u}_3 + ... + \mathbf{Z}_{k}\mathbf{u}_k + \mathbf{e}
}
The function allows to include an arbitrary number of independent random effects with each of them being assumed to follow: \eqn{MVN(\mathbf{0},\mathbf{I}\sigma^2_{u_k})}. If the covariance structure of one random effect is assumed to follow some \eqn{\mathbf{G}} then it it necessary to construct the design matrix for
that random effect as described in Waldmann et al. (2008): \eqn{\mathbf{F} = \mathbf{ZG}^{1/2}}.


}


\usage{
clmm(y, X = NULL , random = NULL, par_random = NULL, niter=10000, 
burnin=5000,scale_e=0,df_e=-2, verbose = TRUE, seed = NULL)
}
\arguments{
\item{y}{vector of phenotypes}
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
At this point the function allows to specify the method for any random term as: 'random' or 'BayesA'. 'random' assumes a common variance for all levels of the 
random effect, 'BayesA' assumes every level of the random effect to have its own distribution and variance as described in Meuwissen et al. (2001).
A wider range of methods is available in the excellent BGLR-package, which also allows phenotypes to be discrete (de los Campos et al. 2013). 

The focus of this function is to allow solving high-dimensional problems that are mixtures of sparse and dense features in the design matrices.
The computational expensive parts of the Gibbs Sampler are parallelized as described in Fernando et al. (2014).
Note that the parallel performance highly depends on the number of observations and features present in the design matrices. 
It is highly recommended to set the number of threads for less than 10000 observations (length of phenotype vector) to 1 using: \code{set_num_threads(1)} before
running a model. Even for larger sample sizes the parallel performance still depends on the dimension of the feature matrices.
Good results in terms of parallel scaling were observed starting from 50000 observations and more than 10000 features (i.e. number of markers).
Single threaded performance is very good thanks to smart computations during gibbs sampling (Fernando, 2013 (personal communication), de los Campos et al., 2009) 
and the use of efficient Eigen-methods for dense and sparse algebra.


}

\value{
List of 3 + number of random effects:
  \item{Residual_Variance}{List of 4:
     \itemize{
       \item{Posterior_Mean - Mean estimate of the residual variance}
       \item{Posterior - Distribution of residual variance}
       \item{scale_prior - scale parameter that has been assigned}
       \item{df_prior - degrees of freedom that have been assigned}
      }
   }

  \item{Predicted}{numeric vector of predicted values}
  \item{Effect_1}{List of 4:
     \itemize{
       \item{type - dense or sparse design matrix}
       \item{method - method that has been used = "fixed"}
       \item{scale_prior - scale parameter that has been assigned}
       \item{df_prior - degrees of freedom that have been assigned}
       \item{posterior - list of 1 = mean of the solution for fixed effects}
      }
   }

Susequently as many additional items as random effects of the following form
  \item{Effect_k}{List of 4:
     \itemize{
       \item{type - dense or sparse design matrix}
       \item{method - method that has been used}
       \item{scale_prior - scale parameter that has been assigned}
       \item{df_prior - degrees of freedom that have been assigned}
       \item{posterior - list of 3 
          \itemize{
            \item{estimates_mean - mean solutions for random effects}
            \item{variance_mean - mean variance}
            \item{variance - distribution of variance}
          }
        }
      }
    }


}


\author{
Claas Heuer

Credits:
Xiaochen Sun (Iowa State University, Ames) gave strong assistance in the theoretical parts and contributed in the very first implementation
of the Gibbs Sampler. Essential parts were adopted from the BayesC-implementation of Rohan Fernando and the BLR-package of Gustavo de los Campos.
The idea of how to parallelize the single site Gibbs Sampler came from Rohan Fernando (2013). 

}

\references{
de los Campos, G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel, and J. M. Cotes. "Predicting Quantitative Traits With Regression Models for Dense Molecular Markers and Pedigree." Genetics 182, no. 1 (May 1, 2009): 375-85. doi:10.1534/genetics.109.101501.

Waldmann, Patrik, Jon Hallander, Fabian Hoti, and Mikko J. Sillanpaa. "Efficient Markov Chain Monte Carlo Implementation of Bayesian Analysis of Additive and Dominance Genetic Variances in Noninbred Pedigrees." Genetics 179, no. 2 (June 1, 2008): 1101-12. doi:10.1534/genetics.107.084160.

Meuwissen, T., B. J. Hayes, and M. E. Goddard. "Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps." Genetics 157, no. 4 (2001): 1819-29.

de los Campos, Gustavo, Paulino Perez Rodriguez, and Maintainer Paulino Perez Rodriguez. "Package 'BGLR,'" 2013. ftp://128.31.0.28/pub/CRAN/web/packages/BGLR/BGLR.pdf.


}

\seealso{\code{\link{clmm.CV}, \link{cGBLUP}, \link{cGWAS.emmax}}}

\examples{

### Running a model with an additive and dominance effect
\dontrun{
# generate random data
rand_data(500,5000)

### compute the relationship matrices
G.A <- cgrm.A(M,lambda=0.01)
G.D <- cgrm.D(M,lambda=0.01)

### generate the list of design matrices for clmm
random = list(t(chol(G.A)),t(chol(G.D)))

### specify options
par_random = list(list(method="random",scale=var(y)/2 ,df=5),
		  list(method="random",scale=var(y)/10,df=5))

### run 

set_num_threads(1)
fit <- clmm(y,random=random,par_random=par_random,niter=5000,burnin=2500)

### inspect results
str(fit)
}
}
\keyword{Genomic Prediction}
