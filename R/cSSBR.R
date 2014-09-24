#
# cSSBR.R
# Claas Heuer, September 2014
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

# cSSBR

cSSBR <- function(data, M, X=NULL , par_random=NULL, scale_e=0, df_e=0, niter=5000, burnin=2500, seed=NULL, verbose=TRUE) {


# double check some stuff to prevent clmm from failing 
# after having already done some heavy computations
if(!is.null(X)) {
  if(class(X)!="matrix") stop("'X' has to be of type 'matrix'") 
  if(anyNA(X)) stop("No NAs allowed in X")
} 

if(!is.null(par_random)) {
  if(length(par_random)!=2) stop("par_random has to be of length 2")
  for(i in 1:length(par_random)) {
    allowed_methods = c("fixed","random","BayesA")
    if(is.null(par_random[[i]]$method)) stop(paste("Define a method for random-effect: ",i,sep=""))
    if(par_random[[i]]$method %in% allowed_methods == FALSE) stop(paste("Method must be one of: ",allowed_methods,sep=""))
  }
}

# obtain the model terms
model_terms <- cSSBR.setup(data,M,verbose)

#################
### Run Model ###
#################

if(verbose) cat(" Running Model\n")
mod <- clmm(y=model_terms$y,
	    X=X,
	    random=list(model_terms$Marker_Matrix,model_terms$Z_residual),
	    par_random=par_random,
	    scale_e=scale_e,
	    df_e=df_e,
	    niter=niter,
	    burnin=burnin,
	    seed=seed,
	    verbose=verbose) 

# obtain breeding values
# FIXME there might be more than just the intercept
bv <- mod$Predicted - mod$Effect_1$posterior$estimates_mean

# add animals that were not in model to the bv-vector
gt_bv <- (M %c% mod$Effect_2$posterior$estimates_mean)[,1]
gt_not_in_model <- rownames(M)[!rownames(M)%in%model_terms$ids]
bv <- c(bv,gt_bv[!rownames(M)%in%model_terms$ids])
names(bv) <- c(model_terms$ids,gt_not_in_model)

# add some important information concerning SSBR
mod$SSBR <- list(ids = model_terms$ids, y=model_terms$y, Marker_Matrix=model_terms$Marker_Matrix, Z_residual = model_terms$Z_residual, Breeding_Values = bv)

return(mod)

}


