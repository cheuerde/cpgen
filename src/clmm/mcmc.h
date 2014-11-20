/*
// mcmc.h
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
// This is the main class for 'clmm'. It did not necessarily have to be 
// a template class, but I wanted compiletime polymorphism in order to
// be able to surpass OpenMP entirely. 
// The reason for two different classes is simply the enormous overhead
// induced by OpenMP. In small datasets it does not make any sense to run
// things in parallel, actually it would decrease the performance dramatically.
//
// There are three layers of polymorphism which are achieved
// by a strategy pattern using virtual classes or by templates:
//
// 1) Single or Multithreaded:
//    The behaviour is defined by the base-class: base_methods_abstract ('base_methods_abstract.h')
//    from which the two derived classes in 'base_methods.h' inherit.
// 2) Method to apply for a fixed or random effect.
//    This is implemented using the base class in 'function_base.h' from which 
//    , as of now, three classes inherit ('functions.h').
// 3) Sparse or dense design matrices.
//    This is achieved by explicitely specializing the functions in 'base_methods.h' for the two classes:
//    Eigen::Map<Eigen::MatrixXd> and Eigen::MappedSparseMatrix<double> .
//    The derived template-classes in 'functions.h' hold those objects, once assigned from a SEXP pointer.
//    

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatrixXd;
typedef Eigen::MappedSparseMatrix<double> MapSparseMatrixXd;
typedef Eigen::Map<Eigen::VectorXd> MapVectorXd;
typedef Eigen::Map<Eigen::ArrayXd> MapArrayXd;

// Taken from: http://gallery.rcpp.org/articles/sparse-iterators/
typedef MapSparseMatrixXd::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

typedef std::vector<std::map<std::string, int> > mp_container;

//#ifndef _CRAN
// #include "mt_sampler.h"
//#else
// #include "R_sampler.h"
//#endif

#include "mt_sampler.h"
#include "effects.h"
#include "../printer.h"
// only available from gcc 4.7 onwards.
// Rtools uses 4.6.3, so we can't include chrono
// #include <chrono>  //FIXME TIMING
// tradidtional C-Way:
#include <sys/time.h>

////////////////
// MCMC class //
////////////////

template<class F>
class MCMC {

private:

  int niter;
  int  burnin;
  bool full_output;
  bool verbose;
  bool timings;
  bool initialized;

  double scale_e;
  double df_e;

  sampler mcmc_sampler;

public:

  Rcpp::List list_of_design_matrices;
  SEXP X_design_matrix;
  Rcpp::List summary_list;
  Rcpp::List par_fixed;
  Rcpp::List par_random;

  double var_y;
  double post_var_e;
  double mu;
 
  int effiter;
  int n_random;

  int n;

  std::string seed;
 
  std::string name;

  VectorXd y;
  VectorXd ycorr;
  VectorXd mean_var_e;

  mp_container thread_vec;

// this obejct determines whethter we run in parallel or not
  base_methods_abstract * my_base_functions;

  vector<effects> model_effects;

  double var_e;

  bool has_na;
  vector<int> isna;

//  for timings
//  std::chrono::high_resolution_clock::time_point t0; //FIXME TIMING
//  std::chrono::high_resolution_clock::time_point t1; //FIXME TIMING
  struct timeval t0;
  struct timeval t1;
  double time_temp;
  double mean_time_per_iter = 0;

//  void populate(SEXP y_from_R, SEXP X_from_R, SEXP par_fixed_from_R, SEXP list_of_design_matrices_from_R, SEXP par_design_matrices_from_R, SEXP par_from_R, int phenotype_number);
  inline void initialize();
  inline void gibbs();
  inline void summary();
  std::string get_name();
  Rcpp::List get_summary();
// taken from http://stackoverflow.com/questions/307082/cleaning-up-an-stl-list-vector-of-pointers
//  ~MCMC() { while(!model_effects.empty()) delete model_effects.back(), model_effects.pop_back() ;};
//  MCMC() : my_base_functions(new F) {};
   MCMC(SEXP y_from_R, SEXP X_from_R, SEXP par_fixed_from_R, SEXP list_of_design_matrices_from_R, 
   SEXP par_design_matrices_from_R, SEXP par_from_R, int phenotype_number);
//   MCMC(const MyString& mcmc_source);
  ~MCMC(){ delete my_base_functions; };

};

template<class F>
MCMC<F>::MCMC(SEXP y_from_R, SEXP X_from_R, SEXP par_fixed_from_R, SEXP list_of_design_matrices_from_R, SEXP par_design_matrices_from_R, SEXP par_from_R, int phenotype_number) : my_base_functions(0)  {


// Initializing all members
  Rcpp::List par(par_from_R);
  list_of_design_matrices = Rcpp::List(list_of_design_matrices_from_R);
  X_design_matrix = X_from_R;
  MapVectorXd y_temp = MapVectorXd(as<MapVectorXd> (y_from_R));

  niter = Rcpp::as<int>(par["niter"]);
  burnin = Rcpp::as<int>(par["burnin"]);
  full_output = Rcpp::as<bool>(par["full_output"]);
  verbose = Rcpp::as<bool>(par["verbose"]);
  timings = Rcpp::as<bool>(par["timings"]);
  scale_e = Rcpp::as<double>(par["scale_e"]);
  df_e = Rcpp::as<double>(par["df_e"]);
  seed = Rcpp::as<std::string>(par["seed"]);

  par_fixed = Rcpp::List(par_fixed_from_R);
  par_random = Rcpp::List(par_design_matrices_from_R);

  std::ostringstream oss; 
  oss << phenotype_number + 1; 
  name = "PHENOTYPE_" + oss.str();

// this is to ensure that if a couple of models are run at the same time
// that every instance has got a unique seed - openmp
  seed.append(oss.str());
  mcmc_sampler = sampler(seed);
//  my_base_functions = new F;

  y = y_temp;

  mean_var_e = VectorXd::Zero(niter);

  effiter=0;
  post_var_e=0;

  n = y.size();
  n_random = list_of_design_matrices.size() -1;
  if(n_random<1) n_random=1;

// check for NAs in pehnotype vector
// comparison x!=x yields TRUE if x=NA

  has_na = true;
  mu = 0;
  var_y=0;
  isna = vector<int>();
  model_effects = vector<effects>();
  initialized = false;

}



/*

///////////////////////
// Copy Constructor ///
//////////////////////

template<class F>
MCMC<F>::MCMC(const MyString& mcmc_source) {

  niter = mcmc_source.niter;
  burnin = mcmc_source.burnin;
  full_output = mcmc_source.full_output;
  verbose = mcmc_source.verbose;
  initialized = mcmc_source.initialized;
  scale_e = mcmc_source.scale_e;
  df_e = mcmc_source.df_e;

// default copy constructor of class 'sampler' is sufficent (no dynamic objects)
  mcmc_sampler = mcmc_source.mcmc_sampler;


  list_of_design_matrices = mcmc_source.list_of_design_matrices;
  summary_list = mcmc_source.summary_list;
  var_y = mcmc_source.var_y;
  post_var_e = mcmc_source.post_var_e;
  mu = mcmc_source.mu;
 
  effiter = mcmc_source.effiter;
  n_random = mcmc_source.n_random;

  n = mcmc_source.n;

  seed = mcmc_source.seed;
 
  name = mcmc_source.name;

  y = mcmc_source.y;
  ycorr = mcmc_source.ycorr;
  mean_var_e = mcmc_source.mean_var_e;

  thread_vec = mcmc_source.thread_vec;

// First reason for explicit copy constructor
  my_base_functions = new base_methods_abstract(*mcmc_source.my_base_functions);

// FIXME check this
  vector<effects> model_effects;

  var_e = mcmc_source.var_e;

  has_na = mcmc_source.has_na;
  isna = mcmc_source.isna;

}

*/


template<class F>
void MCMC<F>::initialize() {

// here we assign the base_function class outside
// the constructor to prevent the nightmare of having to 
// code copy constructors for every single virtual and 
// derived class

  delete my_base_functions;

  my_base_functions = new F; 

//  isna.clear();

  for(int i=0;i<y.size();i++) { if( y(i)!=y(i) ){ isna.push_back(i);} else { mu+=y(i);} }
//Rcout << endl << "na size: " << isna.size() << endl;
  has_na = isna.size() > 0 ? true : false;

  mu = mu / (y.size() - isna.size());  

// compute variance; missing at random for NAs -- FIXME crucial part
// FIXME assign: if( y(i)!=y(i)  or if( y(i)!=0 ??
  for(int i=0;i<y.size();i++) { if( y(i)!=y(i) ) {y(i) = mu;} else {var_y += (y(i) - mu) * (y(i) - mu);} }
  var_y = var_y / (y.size() - isna.size() - 1);
  var_e = var_y - (var_y / (2 * n_random));

  ycorr = y;

/////////////////////
// Multithreading //
////////////////////

  int n_threads;

// container to store start and length for OpenMP threads - Credit: Hao Cheng (Iowa State University)
// get the number of threads
#pragma omp parallel
{
if(omp_get_thread_num()==0) { n_threads = omp_get_num_threads(); }
}


thread_vec = mp_container(n_threads);


for(int i=0;i<n_threads;i++) {

  thread_vec.at(i)["start"] = i * n / n_threads;
  if(i==(n_threads-1)){thread_vec.at(i)["end"] = n;} else {
  thread_vec.at(i)["end"] = (i+1) * n / n_threads;}
  thread_vec.at(i)["length"] = thread_vec.at(i)["end"] - thread_vec.at(i)["start"];

}


// populate with effects
// include fixed effect

  model_effects.clear();
  model_effects.push_back(effects(X_design_matrix, ycorr.data(),par_fixed, &var_e, var_y,n_random,niter, burnin, full_output)); 

// random effects

  for(int i=0;i<list_of_design_matrices.size();i++){

    SEXP temp_list_sexp = par_random[i];

    Rcpp::List temp_list(temp_list_sexp);

    model_effects.push_back(effects(list_of_design_matrices[i],ycorr.data(),temp_list, &var_e, var_y,n_random,niter, burnin, full_output)); 

  }




//  Rcout << endl << "MCMC I am single_thread" << endl;




// initializing effects

    for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {

    it->initialize(my_base_functions);

    }

    initialized = true;

}


template<class F>
void MCMC<F>::gibbs() {

if(!initialized) initialize();

vector<effects>::iterator it;  

  printer prog(niter);
  if(verbose) { prog.initialize(); }

  for(int gibbs_iter=0; gibbs_iter<niter; gibbs_iter++){

// timings
//    if(timings) t0 = std::chrono::high_resolution_clock::now(); //FIXME TIMING
    if(timings) gettimeofday(&t0, NULL);
    for(it = model_effects.begin(); it != model_effects.end(); it++) {

    it->sample_effects(mcmc_sampler,my_base_functions,thread_vec);

    }

    for(it = model_effects.begin(); it != model_effects.end(); it++) {

    it->sample_variance(mcmc_sampler, gibbs_iter);

    }



// sample residual variance

   var_e = (ycorr.matrix().squaredNorm() + scale_e * df_e) / mcmc_sampler.rchisq(n + df_e);
   mean_var_e(gibbs_iter)=var_e;

// residual noise to NAs -- Reference: de los Campos (2009) - BLR
    if(has_na) 
     {
       for (unsigned int i=0;i<isna.size();i++) {ycorr(isna[i]) = mcmc_sampler.rnorm(0,sqrt(var_e));}
     }

// posterior means

    if(gibbs_iter>burnin) {

      for(it = model_effects.begin(); it != model_effects.end(); it++) {

      it->update_means();

      }

    post_var_e += var_e;
    effiter++;

    }

    if (verbose) {

      prog.DoProgress();

    }

// timings
    if(timings) {

//      t1 = std::chrono::high_resolution_clock::now(); //FIXME TIMING
//      time_temp = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() / 1000.0; //FIXME TIMING
      if(timings) gettimeofday(&t1, NULL);
      time_temp = (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1000;
      mean_time_per_iter += time_temp;
      Rcpp::Rcout << std::endl 
                  << " Iteration: |"   << gibbs_iter + 1 
//                << "|  var_e: |"     << var_e  
                  << "|  secs/iter: |" << time_temp
                  << "|"               << std::endl;

    }

//    if(verbose) { Rcout << endl << "Iteration: " << gibbs_iter << "   s2e: " << var_e; }

  }

//  if(verbose) { Rcout << endl; }


}


template<class F>
void MCMC<F>::summary() {


  summary_list["Residual_Variance"] = Rcpp::List::create(Rcpp::Named("Posterior_Mean") = post_var_e / effiter,
			      Rcpp::Named("Posterior") = mean_var_e,
			      Rcpp::Named("scale_prior") = scale_e,
			      Rcpp::Named("df_prior") = df_e);

  VectorXd yhat = VectorXd::Zero(n);
  for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {

    yhat += it->predict(effiter);

  }
  

  summary_list["Predicted"] = yhat;

  int count = 1;

  std::string list_name;

  for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {
  
    std::ostringstream oss; 
    oss << count;
    list_name = "Effect_" + oss.str();
    summary_list[list_name] = it->get_summary(effiter);
    count++;
  }

  if(timings) summary_list["time_per_iter"] = mean_time_per_iter / niter;
  
}


template<class F>
std::string MCMC<F>::get_name() {


  return name;
  

}


template<class F>
Rcpp::List MCMC<F>::get_summary() {

 return summary_list;

}

