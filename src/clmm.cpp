/*
// clmm.cpp
// Claas Heuer, June 2014
//
// Copyright (C)  2014 Claas Heuer
//
// This file is part of cpgen.
//
// cpgen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// cpgen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/

//
// Here we use two different ways of parallelization.
//
// 1) Using parallel computations in every iteration of 
//    a gibbs-sampler as in Fernando et al., 2014. This is already
//    implemented in the class defined in 'mcmc.h' through the classes
//    defined in 'base_methods.h'
// 2) Running several independent models at the same time, with all of them using
//    the very same design-matrices. This is for cross-validation or simply running
//    one model on a large number of phenotypes. This usually scales perfectly.
//    The most important thing here is, that every running instance of an mcmc-object
//    has its very own random-number-generator. This is achieved by encapsulating
//    the RNG-engine in the class defined in 'mt_sampler.h'. C++11 already comes with
//    a Mersenne-Twister Engine and distribution functions (adopted from the Boost-libraries). 
// 

#include "clmm.h"
#include "printer.h"

typedef vector<MCMC<base_methods_st> > mcmc_st;
typedef vector<MCMC<base_methods_mp> > mcmc_mp;

SEXP clmm(SEXP yR, SEXP XR, SEXP par_XR, SEXP list_of_design_matricesR, SEXP par_design_matricesR, SEXP par_mcmcR, SEXP verboseR, SEXP threadsR){

int threads = as<int>(threadsR); 
int verbose = as<int>(verboseR); 

Rcpp::List list_of_phenotypes(yR);
int p = list_of_phenotypes.size();

// if number of threads is larger than p
// we set the number of threads to p
if(threads > p) threads = p;

//omp_set_dynamic(0);
omp_set_num_threads(threads);

Eigen::setNbThreads(1);
Eigen::initParallel();

printer prog(p / threads);

Rcpp::List summary_out;

mcmc_st vec_mcmc_st;
mcmc_mp vec_mcmc_mp;


if((p>1) | (threads==1)) {

// fill the container with mcmc_objects
  for(int i=0;i<p;i++) {
  

    vec_mcmc_st.push_back(MCMC<base_methods_st>(list_of_phenotypes[i], XR, par_XR, list_of_design_matricesR ,par_design_matricesR, par_mcmcR,i));

  }

  for(mcmc_st::iterator it = vec_mcmc_st.begin(); it != vec_mcmc_st.end(); it++) {

    it->initialize();

  }

  if ((p > 1) & verbose) { prog.initialize(); }

// this looks easy - the work was to allow this step to be parallelized
#pragma omp parallel for 
  for(unsigned int i=0;i<vec_mcmc_st.size();i++){

    vec_mcmc_st.at(i).gibbs();

// verbose

  if((p>1) & verbose) { 
   
    if(omp_get_thread_num()==0) {

      prog.DoProgress(); 

    }

  }



  }


  Rcpp::List Summary;
  for(mcmc_st::iterator it = vec_mcmc_st.begin(); it != vec_mcmc_st.end(); it++) {

    it->summary();
    Summary[it->get_name()] = it->get_summary();     

  }

  summary_out = Summary;

} else {

// if the number of threads is larger than and the number of phenotypes is equal to 1 
// the function runs one parallelized Gibbs Sampler 

// fill the container with mcmc_objects
    for(int i=0;i<1;i++) {
  
      vec_mcmc_mp.push_back(MCMC<base_methods_mp>(list_of_phenotypes[i], XR, par_XR, list_of_design_matricesR ,par_design_matricesR, par_mcmcR,i));

    }

    for(mcmc_mp::iterator it = vec_mcmc_mp.begin(); it != vec_mcmc_mp.end(); it++) {

      it->initialize();
      it->gibbs();

    }


    Rcpp::List Summary;
    for(mcmc_mp::iterator it = vec_mcmc_mp.begin(); it != vec_mcmc_mp.end(); it++) {

      it->summary();
      Summary[it->get_name()] = it->get_summary();     

    }


    summary_out = Summary;

  }


return summary_out;

}



