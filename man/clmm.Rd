\name{clmm}
\alias{clmm}

\title{Linear Mixed Models using Gibbs Sampling}

\description{
This function runs linear mixed models of the following form:
\deqn{
\mathbf{y} = \mathbf{Xb} + \mathbf{Z}_{1}\mathbf{u}_1 + \mathbf{Z}_{2}\mathbf{u}_2 + \mathbf{Z}_{3}\mathbf{u}_3 + ... + \mathbf{Z}_{k}\mathbf{u}_k + \mathbf{e}
}
The function allows to include an arbitrary number of independent random effects with each of them being assumed to follow: \eqn{MVN(\mathbf{0},\mathbf{I}\sigma^2_{u_k})}. If the covariance structure of one random effect is assumed to follow some \eqn{\mathbf{G}} then the design matrix for 
that random effect can be constructed as described in Waldmann et al. (2008): \eqn{\mathbf{F} = \mathbf{ZG}^{1/2}}.
Alternatively, the argument \code{ginverse} can be passed.
}


\usage{
clmm(y, X = NULL, Z = NULL, ginverse = NULL, par_random = NULL,
niter=10000, burnin=5000, scale_e=0, df_e=-2, beta_posterior = FALSE,
verbose = TRUE, timings = FALSE, seed = NULL, use_BLAS=FALSE)
}
\arguments{
\item{y}{vector or list of phenotypes}
\item{X}{fixed effects design matrix of type: \code{matrix} or \code{dgCMatrix}. If omitted a column-vector of ones will be assigned }
\item{Z}{list of design matrices for random effects - every element of the list represents one random effect and may be of type: \code{matrix} or \code{dgCMatrix}}
\item{ginverse}{list of inverse covariance matrices for random effects in \code{Z} (e.g. Inverse of numerator relationship matrix).
                Every element of the list represents one random effect and may be of type: \code{matrix} or \code{dgCMatrix}. 
                Note: If passed, \code{ginverse} must have as many items as \code{Z}. For no ginverse assign \code{NULL} for a particular random effect.}
\item{par_random}{list of options for random effects. If passed, the list must have as many elements as \code{random}. Every element may be a list of 4:
  \itemize{
    \item{\code{scale} - (vector of) scale parameters for the inverse chi-square prior}
    \item{\code{df} - (vector of) degrees of freedom for the inverse chi-square prior}
    \item{\code{method} - method to be used for the random effects, may be: \code{ridge} or \code{BayesA}}
    \item{\code{name} - name for that effect}
    \item{\code{GWAS} - list of options for conducting GWAS using window variance proportions (Fernando et al, 2013):
    \itemize{
      \item{\code{window_size} - number of markers used to form a single window}
      \item{\code{threshold} - window porportion of total variance, used to determine presents of association}
    } }
  }
} 
\item{niter}{number of iterations}
\item{burnin}{number of iterations to be discarded as burnin}
\item{verbose}{prints progress to the screen}
\item{beta_posterior}{save all posterior samples of regression coefficients}
\item{timings}{prints time per iteration to the screen - sets \code{verbose = FALSE}}
\item{scale_e}{scale parameter for the inverse chi-square prior for the residuals}
\item{df_e}{degrees of freedom for the inverse chi-square prior for the residuals}
\item{seed}{seed for the random number generator. If omitted, a seed will be generated based on machine and time}
\item{use_BLAS}{use BLAS library instead of Eigen}
}



\details{
\strong{Single Model run}

At this point the function allows to specify the method for any random term as: 'ridge' or 'BayesA'.
In Ridge Regression the prior distribution of the random
effects is assumed to be normal with a constant variance component, while in
BayesA the marginal prior is a t-distribution, resulting from locus specific
variances with inverse chi-square priors (Gianola et al., 2009).
A wider range of methods is available in the excellent BGLR-package, which also allows 
phenotypes to be discrete (de los Campos et al. 2013). 

The focus of this function is to allow solving high-dimensional problems that are mixtures of sparse and dense features in the design matrices.
The computational expensive parts of the Gibbs Sampler are parallelized as described in Fernando et al. (2014).
Note that the parallel performance highly depends on the number of observations and features present in the design matrices. 
It is highly recommended to set the number of threads for less than 10000 observations (length of phenotype vector) to 1 using: \code{set_num_threads(1)} before
running a model. Even for larger sample sizes the parallel performance still depends on the dimension of the feature matrices.
Good results in terms of parallel scaling were observed starting from 50000 observations and more than 10000 features (i.e. number of markers).
Single threaded performance is very good thanks to smart computations during gibbs sampling (Fernando, 2013 (personal communication), de los Campos et al., 2009) 
and the use of efficient Eigen-methods for dense and sparse algebra.

\strong{Parallel Model runs}

In the case of multiple phenotypes passed to the function as a list,
the main advantage of the function is that several threads can access the very same 
data once assigned, which means that the design matrices only have to be allocated once.
The parallel scaling of this function using multiple phenotypes is almost linear.

In C++:

For every element of the phenotype list a new instance of an MCMC-object will be created.
All the memory allocation needed for running the model is done by the major thread. 
The function then iterates over all objects and runs the gibbs sampler. This step is parallelized,
which means that as many models are being run at the same time as threads available.
All MCMC-objects are totally independent from each other, they only share the same design-matrices.
Every object has its own random-number generator with its own seed which allows perfectly reproducible
results.

\strong{GWAS using genomic windows}

The function allows to specify options to any random effect for conducting genomewide association
studies using prediction vector variances of marker windows as described in Fernando et al. (2013).
In every effective sample of the Gibbs Sampler the sampling variance of the
genotypic value vector \eqn{\mathbf{g} = \mathbf{Zu}} of the particular random 
effect is computed as:
\eqn{\tilde{\sigma}_{g}^2 = \left(\sum_{j=1}^{n} (g_j - \mu_g)^2\right) (n -1)^{-1}}
, with \eqn{\mu_g} being the mean of \eqn{\mathbf{g}} and \eqn{n} the number of observations. 
Then for any window \eqn{w} the sampling variance of \eqn{\mathbf{g_w} =
\mathbf{Z}_w\mathbf{u}_w} is obtained as:
\eqn{\tilde{\sigma}_{g_{w}}^2 = \left(\sum_{j=1}^{n} (g_{w_j} - \mu_{g_w})^2\right) (n -1)^{-1}}, where \eqn{w} indicates the range over the
columns of \eqn{\mathbf{Z}} that forms the window \eqn{w}. The posterior
probability that a window exceeds a specified proportion \eqn{\delta} of
the total variance is estimated by the number of samples in which
\eqn{\frac{\tilde{\sigma}_{g_{w}}^2}{\tilde{\sigma}_{g}^2} > \delta}
divided by the total number of samples. It can be shown that among
marker windows that have a posterior probability \eqn{p} or greater for
having a variance greater than \eqn{\delta} of the total variance, the
proportion of false positive signals (PFP) are expected to be lower
than \eqn{1-p} (Fernando et al. 2004, Fernando et al., 2013). 

}

\value{
List of 4 + number of random effects:
  \item{Residual_Variance}{List of 4:
     \itemize{
       \item{\code{Posterior_Mean} - Mean estimate of the residual variance}
       \item{\code{Posterior} - posterior samples of residual variance}
       \item{\code{scale_prior} - scale parameter that has been assigned}
       \item{\code{df_prior} - degrees of freedom that have been assigned}
      }
   }

  \item{Predicted}{numeric vector of predicted values}
  \item{fixed_effects}{List of 4:
     \itemize{
       \item{\code{type} - dense or sparse design matrix}
       \item{\code{method} - method that has been used = "fixed"}
       \item{\code{posterior} - list of 1 + 1 (if \code{beta_posterior=TRUE}) 
          \itemize{
            \item{\code{estimates_mean} - mean solutions for random effects}
            \item{\code{estimates} - posterior samples of random effects}
          }
        }
      }
   }

Susequently as many additional items as random effects of the following form
  \item{Effect_k}{List of 4  + 1 (if \code{GWAS} options were specified):
     \itemize{
       \item{\code{type} - dense or sparse design matrix}
       \item{\code{method} - method that has been used}
       \item{\code{scale_prior} - scale parameter that has been assigned}
       \item{\code{df_prior} - degrees of freedom that have been assigned}
       \item{\code{posterior} - list of 3 + 1 (if \code{beta_posterior=TRUE}) 
          \itemize{
            \item{\code{estimates_mean} - mean solutions for random effects}
            \item{\code{variance_mean} - mean variance}
            \item{\code{variance} - posterior samples of variance}
            \item{\code{estimates} - posterior samples of random effects}
          }
        }

       \item{\code{GWAS} - list of 9 (if specified)
          \itemize{
            \item{\code{window_size} - number of features (markers) used to form a single window}
            \item{\code{threshold} - window porportion of total variance, used to determine presents of association}
            \item{\code{mean_variance} - mean variance of prediction vector using all windows}
            \item{\code{windows} - identifier}
            \item{\code{start} - starting column for window}
            \item{\code{end} - ending column for window}
            \item{\code{window_variance} - mean variance of prediction vector using this window}
            \item{\code{window_variance_proportion} - mean window proportion of total variance}
            \item{\code{prob_window_var_bigger_threshold} - mean probability that window variance exceeds threshold}
          }
        }
     }
  }


  \item{mcmc}{List of 4 + 1 (if \code{timings=TRUE}):
     \itemize{
       \item{\code{niter} - number of iterations}
       \item{\code{burnin} - number of samples discarded as burnin}
       \item{\code{number_of_samples} - number of samples used to estimate posterior means}
       \item{\code{seed} - seed used for the random number generator}
       \item{\code{time_per_iter} - average seconds per iteration}
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

Gianola, D., de Los Campos, G., Hill, W.G., Manfredi, E., Fernando, R.: "Additive genetic variability and the
bayesian alphabet." Genetics 183(1), 347-363 (2009)

de los Campos, G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel, and J. M. Cotes. "Predicting Quantitative Traits With Regression Models for Dense Molecular Markers and Pedigree." Genetics 182, no. 1 (May 1, 2009): 375-85. doi:10.1534/genetics.109.101501.

Waldmann, Patrik, Jon Hallander, Fabian Hoti, and Mikko J. Sillanpaa. "Efficient Markov Chain Monte Carlo Implementation of Bayesian Analysis of Additive and Dominance Genetic Variances in Noninbred Pedigrees." Genetics 179, no. 2 (June 1, 2008): 1101-12. doi:10.1534/genetics.107.084160.

Meuwissen, T., B. J. Hayes, and M. E. Goddard. "Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps." Genetics 157, no. 4 (2001): 1819-29.

de los Campos, Gustavo, Paulino Perez Rodriguez, and Maintainer Paulino Perez Rodriguez. "Package 'BGLR,'" 2013. ftp://128.31.0.28/pub/CRAN/web/packages/BGLR/BGLR.pdf.

Fernando, R.L., Dekkers, J.C., Garrick, D.J. "A class of bayesian methods to combine large numbers of
genotyped and non-genotyped animals for whole-genome analyses." Genetics Selection Evolution 46(1), 50 (2014)

Fernando, R., Nettleton, D., Southey, B., Dekkers, J., Rothschild, M., Soller, M. "Controlling the proportion of
false positives in multiple dependent tests." Genetics 166(1), 611-619 (2004)

Fernando, Rohan L., and Dorian Garrick. "Bayesian methods applied to GWAS." Genome-Wide Association Studies and Genomic Prediction. Humana Press, 2013. 237-274.


}

\seealso{\code{\link{cGBLUP}, \link{cSSBR}, \link{cGWAS.emmax}}}

\examples{

\dontrun{

############################################################# 
### Running a model with an additive and dominance effect ###
############################################################# 

# generate random data
rand_data(500,5000)

### compute the relationship matrices
G.A <- cgrm.A(M,lambda=0.01)
G.D <- cgrm.D(M,lambda=0.01)

### generate the list of design matrices for clmm
Z_list = list(t(chol(G.A)),t(chol(G.D)))

### specify options
par_random = list(list(method="ridge",scale=var(y)/2 ,df=5, name="add"),
		  list(method="ridge",scale=var(y)/10,df=5, name="dom"))

### run 

set_num_threads(1)
fit <- clmm(y = y, Z = Z_list, par_random=par_random, niter=5000, burnin=2500)

### inspect results
str(fit)


########################
### Cross Validation ###
########################
 
### 4-fold cross-validation with one repetition:
# generate random data
rand_data(500,5000)

### compute the list of masked phenotype-vectors for CV
y_CV <- cCV(y, fold=4, reps=1)


### Cross Validation using GBLUP
G.A <- cgrm.A(M,lambda=0.01)


### generate the list of design matrices for clmm
Z_list = list(t(chol(G.A)))

### specify options
h2 = 0.3
scale = unlist(lapply(y_CV,function(x)var(x,na.rm=T))) * h2
df = rep(5,length(y_CV))
par_random = list(list(method="ridge",scale=scale,df=df, name="animal"))

### run model with 4 threads
set_num_threads(4)
fit <- clmm(y = y_CV, Z = Z_list, par_random=par_random, niter=5000, burnin=2500)

### inspect results
str(fit)

### obtain predictions
pred <- get_pred(fit)

### prediction accuracy
get_cor(pred,y_CV,y)

########################################################
### GWAS using Bayesian Regression on marker windows ###
########################################################
 
# generate random data
rand_data(500,5000)

### generate the list of design matrices for clmm
Z_list = list(M)

### specify options
h2 = 0.3
scale = var(y) * h2
df = 5
# specifying the model
# Here we use ridge regression on the marker covariates
# and define a window size of 100 and a threshold of 0.01
# which defines the proportion of genetic variance accounted
# for by a single window
par_random = list(list(method="ridge",scale=scale,df=df,
             GWAS=list(window_size=100, threshold=0.01), name="markers"))

### run 
set_num_threads(1)
fit <- clmm(y =y, Z = Z_list, par_random=par_random, niter=2000, burnin=1000)

### inspect results
str(fit)

### extract GWAS part
gwas <- fit$markers$GWAS

# plot window variance proportions
plot(gwas$window_variance_proportion)

##########################################################
### Sparse Animal Model using the pedigreemm milk data ###
##########################################################

# cpgen offers two ways of running models with random effects that
# are assumed to follow some covariance structure.
# 1) Construct the Covariance matrix and pass the cholesky of that
#    as design matrix for that random effect
# 2) Construct the inverse of the covariance matrix (ginverse) and pass the design
#    matrix 'Z' that relates observations to factors in ginverse in conjunction
#    with the symmetric ginverse.

# This is approach 2) which is more convenient for pedigree based
# animal models, as ginverse (Inverse of numerator relationship matrix) is 
# very sparse and can be obtained very efficiently

set_num_threads(1)

# load the data
data(milk)

# get Ainverse
# Ainv <- as(getAInv(pedCows),"dgCMatrix")

T_Inv <- as(pedCows, "sparseMatrix")
D_Inv <- Diagonal(x=1/Dmat(pedCows))
Ainv<-t(T_Inv) %*% D_Inv %*% T_Inv
dimnames(Ainv)[[1]]<-dimnames(Ainv)[[2]] <-pedCows@label
Ainv <- as(Ainv, "dgCMatrix")

# We need to construct the design matrix.
# Therefore we create a second id column with factor levels
# equal to the animals in the pedigree
milk$id2 <- factor(as.character(milk$id), levels = pedCows@label)

# set up the design matrix
Z <- sparse.model.matrix(~ -1 + id2, data = milk, drop.unused.levels=FALSE)

# run the model
niter = 5000
burnin = 2500

modAinv <- clmm(y = as.numeric(milk$milk), Z = list(Z), ginverse = list(Ainv), 
                niter = niter, burnin = burnin)

# This is approach 1) run an equivalent model using the cholesky of A

# get L from A = LL'
L <- as(t(relfactor(pedCows)),"dgCMatrix")

# match with ids
ZL <- L[match(milk$id, pedCows@label),]

# run the model
modL <- clmm(as.numeric(milk$milk), Z= list(ZL), 
             niter = niter, burnin = burnin)


### a more advanced model

# y = Xb + Zu + a + e
#
# u = permanent environment of animal
# a = additive genetic effect of animal

Zpe <- sparse.model.matrix(~ -1 + id2, drop.unused.levels = TRUE, data = milk)

# make X and account for lactation and herd
X <- sparse.model.matrix(~ 1 + lact + herd, data = milk)

niter = 10000
burnin = 2500

mod2 <- clmm(as.numeric(milk$milk), X = X, Z = list(Zpe,Z), ginverse = list(NULL, Ainv), 
                        niter = niter, burnin = burnin))


# run all phenotypes in the milk dataset at once in parallel

Y <- list(as.numeric(milk$milk),as.numeric(milk$fat),as.numeric(milk$prot),as.numeric(milk$scs))

set_num_threads(4)

# ginverse version
model <- clmm(Y, X = X, Z = list(Zpe,Z), ginverse = list(NULL, Ainv), 
              niter = niter, burnin = burnin)

# get heritabilities and repeatabilities with their standard deviations

heritabilities <- array(0, dim=c(length(Y),2))
colnames(heritabilities) <- c("h2","sd")

# only use post-burnin samples
range <- (burnin+1):niter

# h2
heritabilities[,1] <- unlist(lapply(model, function(x)
                                    mean(
                                         x$Effect_2$posterior$variance[range] /
                                         (x$Effect_1$posterior$variance[range] + 
                                         x$Effect_2$posterior$variance[range] + 
                                         x$Residual_Variance$Posterior[range]))
                                        )
                             )

# standard deviation of h2
heritabilities[,2] <- unlist(lapply(model, function(x)
                                    sd(
                                       x$Effect_2$posterior$variance[range] /
                                       (x$Effect_1$posterior$variance[range] + 
                                       x$Effect_2$posterior$variance[range] + 
                                       x$Residual_Variance$Posterior[range]))
                                      )
                             )

}
}
\keyword{Genomic Prediction}
\keyword{GWAS}
