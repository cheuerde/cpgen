/*
// functions.h
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


// default copy constructor is enough. it will call base and the only data
// member is 'design_matrix' which copy-constructor is being called (Eigen).

#include "function_base.h"

/////////////////////
// function_fixed //
////////////////////


template<class T>
class function_fixed: public function_base {
public:


T design_matrix;

inline void initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects);
inline void sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec);
inline void sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter){};
inline void update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var); 
inline Eigen::VectorXd predict(VectorXd& mean_estimates, int start, int length);
inline Rcpp::List summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior);
function_fixed(SEXP design_matrix_pointer) : design_matrix(as<T>(design_matrix_pointer)){};
~function_fixed(){};
//~function_fixed(){Rcout << endl << "function_fixed destructor" << endl;};

};

template<class T>
void function_fixed<T>::initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects) {

  base_fun->initialize(design_matrix,xtx, columns);
  var = VectorXd::Constant(columns,var_y * 1000);
  mean_var = VectorXd::Constant(1,0);
  initialized = 1;

}

template<class T>
void function_fixed<T>::sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec) {

  base_fun->sample_effects(design_matrix,xtx,estimates,ycorr,var,var_e, mcmc_sampler, thread_vec);

}

template<class T>
void function_fixed<T>::update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var) {

  mean_estimates += estimates;
   
}

template<class T>
Eigen::VectorXd function_fixed<T>::predict(VectorXd& mean_estimates, int start, int length) {

  return design_matrix.block(0,start,design_matrix.rows(),length) * mean_estimates.segment(start,length);
   
}

template<class T>
Rcpp::List function_fixed<T>::summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior) {

  return Rcpp::List::create(Rcpp::Named("estimates_mean") = mean_estimates/effiter);
 
}


/////////////////////
// function_random //
/////////////////////

template<class T>
class function_random: public function_base {

public:

T design_matrix;
inline void initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects);
inline void sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec);
inline void sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter);
inline void update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var); 
inline Eigen::VectorXd predict(VectorXd& mean_estimates, int start, int length);
inline Rcpp::List summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior);
function_random(SEXP design_matrix_pointer) : design_matrix(as<T>(design_matrix_pointer)){};
~function_random(){};

};


template<class T>
void function_random<T>::initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects) {

  base_fun->initialize(design_matrix,xtx, columns);
// FIXME : find appropiate starting value
  var = VectorXd::Constant(columns,var_y / 2 * total_number_effects);
  mean_var = VectorXd::Constant(1,0);
  initialized = 1;

}

template<class T>
void function_random<T>::sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec) {

  base_fun->sample_effects(design_matrix,xtx,estimates,ycorr,var,var_e, mcmc_sampler, thread_vec);

}

template<class T>
void function_random<T>::sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter) {

   double var_scalar = (scale * df + estimates.squaredNorm())/mcmc_sampler.rchisq(columns + df);
   var.noalias() = VectorXd::Constant(columns, var_scalar);  
   var_posterior.at(niter) = var_scalar;

}


template<class T>
void function_random<T>::update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var) {

  mean_estimates += estimates;
  mean_var(0) += var(0);
   
}


template<class T>
Eigen::VectorXd function_random<T>::predict(VectorXd& mean_estimates, int start, int length) {

  return design_matrix.block(0,start,design_matrix.rows(),length) * mean_estimates.segment(start,length);
   
}

template<class T>
Rcpp::List function_random<T>::summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior) {

  return Rcpp::List::create(Rcpp::Named("estimates_mean") = mean_estimates/effiter,
		            Rcpp::Named("variance_mean") = mean_var(0)/effiter,
                            Rcpp::Named("variance") = var_posterior);
   
}



/////////////////////
// function_bayesA //
/////////////////////

template<class T>
class function_bayesA: public function_base {
public:

T design_matrix;

inline void initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects);
inline void sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec);
inline void sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter);
inline void update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var); 
inline Eigen::VectorXd predict(VectorXd& mean_estimates, int start, int length);
inline Rcpp::List summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior);
function_bayesA(SEXP design_matrix_pointer) : design_matrix(as<T>(design_matrix_pointer)){};
~function_bayesA(){};

};


template<class T>
void function_bayesA<T>::initialize(base_methods_abstract*& base_fun,VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects) {

  base_fun->initialize(design_matrix,xtx, columns);
// FIXME : find appropiate starting value
  var = VectorXd::Constant(columns,var_y / 2 * total_number_effects / columns);
  mean_var = VectorXd::Constant(columns,0);
  initialized = 1;

}

template<class T>
void function_bayesA<T>::sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec) {

  base_fun->sample_effects(design_matrix,xtx,estimates,ycorr,var,var_e, mcmc_sampler, thread_vec);

}

template<class T>
void function_bayesA<T>::sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter) {

  for (int i=0;i<columns;i++) { var(i) = (scale * df + estimates(i) * estimates(i)) / mcmc_sampler.rchisq(1 + df); }

}

template<class T>
void function_bayesA<T>::update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var) {

  mean_estimates += estimates;
  mean_var += var;
   
}


template<class T>
Eigen::VectorXd function_bayesA<T>::predict(VectorXd& mean_estimates, int start, int length) {

  return design_matrix.block(0,start,design_matrix.rows(),length) * mean_estimates.segment(start,length);
   
}

template<class T>
Rcpp::List function_bayesA<T>::summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior) {

  return Rcpp::List::create(Rcpp::Named("estimates_mean") = mean_estimates/effiter,
		            Rcpp::Named("variances_mean") = mean_var/effiter);
   
}



