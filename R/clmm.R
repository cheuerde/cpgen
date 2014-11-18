#
# clmm.R
# Claas Heuer, June 2014
#
# Copyright (C)  2014 Claas Heuer
#
# This file is part of cpgen.
#
# cpgen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cpgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# in file.path(R.home("share"), "licenses").  If not, see
# <http://www.gnu.org/licenses/>.
#

# clmm

clmm <- function(y, X = NULL , random = NULL, par_random = NULL, niter=10000, burnin=5000,scale_e=0,df_e=-2, verbose = TRUE, timings = FALSE, seed = NULL){

default_scale = 0
default_df = -2

allowed=c("numeric")
a = class(y)

if(sum(a%in%allowed)!=1) { stop("phenotypes must match one of the following types: 'numeric'") }

if(class(y) == "list") {
  p = length(y)
  if (sum(unlist(lapply(y,is.vector))) != p) stop("phenotypes must be supplied as vectors")
  n = unlist(lapply(y,length))
  if (sum(n[1]!=n) > 0) stop("phenoytpe vectors must have same length")
  n = n[1]
  if(is.null(names(y))) names(y) <- paste("phenotype: ",1:p,sep="") 
  } else {
    if(!is.vector(y)) stop("phenotype must be supplied as vector")
    n <- length(y) 
    y <- list(y)
    names(y) <- "phenotype: 1"
  }
#isy <- (1:n)[!is.na(y)]

if(is.null(X)) {
  X = array(1,dim=c(n,1))
  par_fixed <- list(scale=default_scale,df=default_df,sparse_or_dense="dense",method="fixed")
  } else {                   
    if(X_is_ok(X,n,"fixed")) {
      if(class(X) == "matrix") { type = "dense" } else { type = "sparse" }
      par_fixed <- list(scale=default_scale,df=default_df,sparse_or_dense=type,method="fixed")    
    }
  }

if(is.null(random)) {
  random = list()
  par_random = list()
  } else {

  if(is.null(par_random)) {
    par_random<-list(length(random))
    if(is.null(names(random))) names(random) = 1:length(random)
    for(i in 1:length(random)){

      if(X_is_ok(random[[i]],n,names(random)[i])) {
      method = "random"
      if(class(random[[i]]) == "matrix") type = "dense"
      if(class(random[[i]]) == "dgCMatrix") type = "sparse"
      par_random[[i]] = list(scale=default_scale,df=default_df,sparse_or_dense=type,method=method) }
    }

    } else {

      if(length(par_random) != length(random)) stop(" 'par_effects' must have as many items as 'random' ")

      for(i in 1:length(par_random)) {
        X_is_ok(random[[i]],n,names(random)[i])
        allowed_methods = c("fixed","random","BayesA")
        if(is.null(par_random[[i]]$method)) stop(paste("Define a method for random-effect: ",i,sep=""))
        if(par_random[[i]]$method %in% allowed_methods == FALSE) stop(paste("Method must be one of: ",allowed_methods,sep=""))

        if(is.null(par_random[[i]]$df) | !is.numeric(par_random[[i]]$df) | length(par_random[[i]]$df) > 1)  {

          if(par_random[[i]]$method == "BayesA") { 

            par_random[[i]]$df = 4.0 } else { 

              par_random[[i]]$df = default_df }
        
          }

        if(is.null(par_random[[i]]$scale) | !is.numeric(par_random[[i]]$scale) | length(par_random[[i]]$scale) > 1) {

          if(par_random[[i]]$method == "BayesA") { 

## FIXME the scale for BayesA doesnt make sense for more than one phenotype
            dfA = par_random[[i]]$df
            par_random[[i]]$scale = (((dfA-2)/dfA) * (var(y[[1]],na.rm=T) / 5)) / (ncol(random[[i]]) * sum(ccolmv(random[[i]])**2,na.rm=T)) } else {

            par_random[[i]]$scale = default_scale }
     
        }



          if(class(random[[i]]) == "matrix") { type = "dense" } else { type = "sparse" }
          par_random[[i]]$sparse_or_dense = type 

        }   
      }

    }
  



# RNG Seed based on system time and process id
# Taken from: http://stackoverflow.com/questions/8810338/same-random-numbers-every-time
if(is.null(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

## FIXME check arguments

if(timings) verbose = FALSE
par_mcmc = list(niter=niter, burnin=burnin, full_output=TRUE, verbose=verbose, timings = timings, scale_e = scale_e, df_e = df_e, seed = as.character(seed))

#set.seed(seed)
#cat(paste("\n seed: ",seed,"\n",sep="")) 
 .Call("clmm",y, X , par_fixed ,random, par_random ,par_mcmc, verbose=verbose, options()$cpgen.threads, PACKAGE = "cpgen" )[[1]]

#return(list(X=X,par_fixed=par_fixed,random=random,par_random=par_random,par_mcmc=par_mcmc))

}





# clmm.CV

clmm.CV <- function(y, X = NULL , random = NULL, par_random = NULL, niter=10000, burnin=5000,scale_e=0,df_e=-2, verbose = TRUE, timings = FALSE, seed = NULL){

default_scale = 0
default_df = -2

allowed=c("numeric","list")
a = class(y)

if(sum(a%in%allowed)!=1)  stop(paste("phenotypes must match one of the following types: ",allowed,sep="")) 

if(class(y) == "list") {
  p = length(y)
  if (sum(unlist(lapply(y,is.vector))) != p) stop("phenotypes must be supplied as vectors")
  n = unlist(lapply(y,length))
  if (sum(n[1]!=n) > 0) stop("phenoytpe vectors must have same length")
  n = n[1]
  if(is.null(names(y))) names(y) <- paste("phenotype: ",1:p,sep="") 
  } else {
    if(!is.vector(y)) stop("phenotype must be supplied as vector")
    n <- length(y) 
    y <- list(y)
    names(y) <- "phenotype: 1"
  }
#isy <- (1:n)[!is.na(y)]

if(is.null(X)) {
  X = array(1,dim=c(n,1))
  par_fixed <- list(scale=0,df=0,sparse_or_dense="dense",method="fixed")
  } else {                   
    if(X_is_ok(X,n,"fixed")) {
      if(class(X) == "matrix") { type = "dense" } else { type = "sparse" }
      par_fixed <- list(scale=default_scale,df=default_df,sparse_or_dense=type,method="fixed")    
    }
  }

if(is.null(random)) {
  random = list()
  par_random = list()
  } else {

  if(is.null(par_random)) {
    par_random<-list(length(random))
    if(is.null(names(random))) names(random) = 1:length(random)
    for(i in 1:length(random)){

      if(X_is_ok(random[[i]],n,names(random)[i])) {
      method = "random"
      if(class(random[[i]]) == "matrix") type = "dense"
      if(class(random[[i]]) == "dgCMatrix") type = "sparse"
      par_random[[i]] = list(scale=default_scale,df=default_df,sparse_or_dense=type,method=method) }
    }

    } else {

      if(length(par_random) != length(random)) stop(" 'par_effects' must have as many items as 'random' ")

      for(i in 1:length(par_random)) {
        X_is_ok(random[[i]],n,names(random)[i])
        allowed_methods = c("fixed","random","BayesA")
        if(is.null(par_random[[i]]$method)) stop(paste("Define a method for random-effect: ",i,sep=""))
        if(par_random[[i]]$method %in% allowed_methods == FALSE) stop(paste("Method must be one of: ",allowed_methods,sep=""))

        if(is.null(par_random[[i]]$df) | !is.numeric(par_random[[i]]$df) | length(par_random[[i]]$df) > 1)  {

          if(par_random[[i]]$method == "BayesA") { 

            par_random[[i]]$df = 4.0 } else { 

              par_random[[i]]$df = default_df }
        
          }

        if(is.null(par_random[[i]]$scale) | !is.numeric(par_random[[i]]$scale) | length(par_random[[i]]$scale) > 1) {

          if(par_random[[i]]$method == "BayesA") { 

## FIXME the scale for BayesA doesnt make sense for more than one phenotype
            dfA = par_random[[i]]$df
            par_random[[i]]$scale = (((dfA-2)/dfA) * (var(y[[1]],na.rm=T) / 5)) / (ncol(random[[i]]) * sum(ccolmv(random[[i]])**2,na.rm=T)) } else {

            par_random[[i]]$scale = default_scale }
     
        }

          if(class(random[[i]]) == "matrix") { type = "dense" } else { type = "sparse" }
          par_random[[i]]$sparse_or_dense = type 

        }   

      }

    }
  



# RNG Seed based on system time and process id
# Taken from: http://stackoverflow.com/questions/8810338/same-random-numbers-every-time
if(is.null(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

## Timings are not allowed here
timings=FALSE
par_mcmc = list(niter=niter, burnin=burnin, full_output=TRUE, verbose=FALSE, timings = timings, scale_e = scale_e, df_e = df_e, seed = as.character(seed))


#set.seed(seed)
#cat(paste("\n seed: ",seed,"\n",sep="")) 
.Call("clmm",y, X , par_fixed ,random, par_random ,par_mcmc, verbose=verbose, options()$cpgen.threads, PACKAGE = "cpgen" )

#return(list(X=X,par_fixed=par_fixed,random=random,par_random=par_random,par_mcmc=par_mcmc))

}


get_pred <- function(mod) {

return(matrix(unlist(lapply(mod,function(x)x$Predicted)),ncol=length(mod),nrow=length(mod[[1]]$Predicted)))

}


get_cor <- function(predictions,cv_pheno,y) {

cv_vec <- matrix(unlist(cv_pheno),nrow=length(y),ncol=length(cv_pheno),byrow=FALSE)

mean_pred <- rep(NA,nrow(predictions))

for(i in 1:nrow(predictions)) {

mean_pred[i] <- mean(predictions[i,which(is.na(cv_vec[i,]))])

}

return(cor(mean_pred,y,use="pairwise.complete.obs"))

}

# cGBLUP


cGBLUP <- function(y,G,X=NULL, scale_a = 0, df_a = -2, scale_e = 0, df_e = -2,niter = 10000, burnin = 5000, seed = NULL, verbose=TRUE){

isy <- (1:length(y))[!is.na(y)]

if(verbose) cat("\nComputing Eigen Decomposition\n")

if(length(isy) < length(y)) {
  UD <- eigen(G[isy,isy]) } else {
    UD <- eigen(G) }

n <- length(isy)

if(is.null(X)) X = rep(1,length(y[isy]))

Uy <- (t(UD$vectors)%c%y[isy])[,1]
UX <- t(UD$vectors)%c%X

D_sqrt <- sqrt(UD$values)
Z<-sparseMatrix(i=1:n,j=1:n,x=D_sqrt)

par_random <- list(list(scale=scale_a,df=df_a,sparse_or_dense="sparse",method="random"))

if(verbose) cat("Running Model\n")
if(is.null(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

# set the number of threads to 1 for clmm
old_threads <- get_num_threads()
set_num_threads(1,silent=TRUE)

mod <- clmm(Uy, UX , list(Z), par_random ,scale_e=scale_e,df_e=df_e, verbose=verbose, niter=niter,burnin=burnin,seed=seed)

# set number of threads to old value
set_num_threads(old_threads,silent=TRUE)

u <- rep(NA,length(y))


u[isy] <- UD$vectors %c% (D_sqrt * mod[[4]]$posterior$estimates_mean)
if(length(isy) < length(y)) { u[-isy] <- G[-isy,isy] %c% csolve(G[isy,isy],u[isy]) }

e<-mod$Residual_Variance$Posterior

return(list(var_e = mod$Residual_Variance$Posterior_Mean,
	    var_a = mod[[4]]$posterior$variance_mean, 
            b = mod[[3]]$posterior$estimates_mean,
	    a = u,
	    posterior_var_e = mod$Residual_Variance$Posterior,
	    posterior_var_a = mod[[4]]$posterior$variance))


}






X_is_ok <- function(X,n,name) {

allowed=c("matrix","dgCMatrix")
a = class(X)
#if(sum(a%in%allowed)!=1) stop(paste(c("lol","rofl"))) 
if(sum(a%in%allowed)!=1) stop(paste("design matrix '",name,"' must match one of the following types: ",paste(allowed,collapse=" or "),sep="")) 

if(anyNA(X)) { stop(paste("No NAs allowed in design matrix '", name,"'", sep="")) } 
      
if(a=="matrix" | a=="dgCMatrix") { if(nrow(X) != n) stop(paste("Number of rows in design matrix '",name,"' doesnt match number of observations in y",sep="")) }


return(1) 

}
















