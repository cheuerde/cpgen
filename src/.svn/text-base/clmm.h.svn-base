#ifndef _cpgen_clmm_H
#define _cpgen_clmm_H

#ifdef _OPENMP
  #define has_openmp 1
  #include <omp.h>
#else 
  #define has_openmp 0
  #define omp_get_num_threads() 1
  #define omp_set_num_threads(x) 1
  #define omp_get_max_threads() 1
  #define omp_get_thread_limit() 1
  #define omp_set_dynamic(x) 1
  #define omp_get_thread_num() 0
#endif

#include <RcppEigen.h>
#include <vector>
#include "clmm/mcmc.h"


RcppExport SEXP clmm(SEXP yR, SEXP XR, SEXP par_XR, SEXP list_of_design_matricesR, SEXP par_design_matricesR, SEXP par_mcmcR, SEXP verboseR, SEXP threadsR);

#endif


