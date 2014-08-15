/*
// effects.h
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

#include "functions.h"


////////////////////
// Effects class //
///////////////////


class effects {
public:

  function_base * my_functions;

  SEXP design_matrix;
  int total_number_effects;
  int columns;
  int niter;
  int burnin;
  double scale;
  double df;
  double * var_e;
  VectorXd mean_var;
//  double var_scalar;
  VectorXd var;
  vector<double> var_posterior;
  VectorXd xtx;
  double * ycorr;
//  MatrixXd mean_var;

  std::string sparse_or_dense;
  std::string method;
  bool initialized;
  bool full_output;
  double var_y;
  VectorXd estimates;
  VectorXd mean_estimates;

  effects(SEXP design_matrix_from_MCMC, double * ycorr_from_MCMC, Rcpp::List list_from_MCMC, double * var_e_from_MCMC, double var_y_from_MCMC,int total_number_effects_from_MCMC,int niter_from_MCMC, int burnin_from_MCMC, bool full_output_from_MCMC);
  inline void initialize(base_methods_abstract*& base_fun);  
  inline void sample_effects(sampler& mcmc_sampler,base_methods_abstract*& base_fun, mp_container& thread_vec); 
  inline void sample_variance(sampler& mcmc_sampler, int& niter_from_MCMC); 
  inline void update_means();
  inline Eigen::VectorXd predict(int effiter);
  Rcpp::List get_summary(int effiter);

  effects() : my_functions(0){};
// deleting the null-pointer does nothing, so we dont have to check whether a new function class has been assigned
// but if: it will be deleted
  ~effects() {delete my_functions;};
//  ~effects() {delete my_functions;Rcout << endl << "effects_destructor" << endl;};


};


effects::effects(SEXP design_matrix_from_MCMC, double * ycorr_from_MCMC, Rcpp::List list_from_MCMC, double * var_e_from_MCMC, double var_y_from_MCMC,int total_number_effects_from_MCMC,int niter_from_MCMC, int burnin_from_MCMC, bool full_output_from_MCMC) : my_functions(0) {
 
  design_matrix = design_matrix_from_MCMC;
  scale = as<double>(list_from_MCMC["scale"]);
  df = as<double>(list_from_MCMC["df"]);
  sparse_or_dense = as<std::string>(list_from_MCMC["sparse_or_dense"]);
  method = as<std::string>(list_from_MCMC["method"]);
  total_number_effects = total_number_effects_from_MCMC;
  ycorr = ycorr_from_MCMC;
  var_e = var_e_from_MCMC;
  var_y = var_y_from_MCMC;
  full_output = full_output_from_MCMC;
  niter = niter_from_MCMC;
  burnin = burnin_from_MCMC;

  initialized = false;


//  Rcout << endl << "effects I am single_thread" << endl;

}



void effects::initialize(base_methods_abstract*& base_fun){

  delete my_functions;

  if(sparse_or_dense == "sparse") {
    if(method == "fixed") my_functions = new function_fixed<MapSparseMatrixXd>(design_matrix); 
    if(method == "random") my_functions = new function_random<MapSparseMatrixXd>(design_matrix); 
    if(method == "BayesA") my_functions = new function_bayesA<MapSparseMatrixXd>(design_matrix); 
  } else {
      if(method == "fixed") my_functions = new function_fixed<MapMatrixXd>(design_matrix); 
      if(method == "random") my_functions = new function_random<MapMatrixXd>(design_matrix); 
      if(method == "BayesA") my_functions = new function_bayesA<MapMatrixXd>(design_matrix); 
    }


  my_functions->initialize(base_fun,xtx,columns,var, mean_var,var_y,initialized,total_number_effects);
  estimates = VectorXd::Zero(columns);
  mean_estimates = VectorXd::Zero(columns);
  var_posterior.resize(niter);
  initialized = true;
  

}

void effects::sample_effects(sampler& mcmc_sampler,base_methods_abstract*& base_fun, mp_container& thread_vec){

  my_functions->sample_effects(base_fun,xtx,estimates, ycorr, var,var_e, mcmc_sampler,thread_vec);

}


void effects::sample_variance(sampler& mcmc_sampler, int& niter_from_MCMC){

  my_functions->sample_variance(var, scale, df, estimates, columns, mcmc_sampler, var_posterior, niter_from_MCMC);

}



void effects::update_means(){


  my_functions->update_means(mean_estimates,estimates,mean_var,var);

//  mean_estimates += estimates;
//  mean_var += var(0);
//  effiter++;

}


Eigen::VectorXd effects::predict(int effiter){


  VectorXd predicted;
  if(sparse_or_dense == "sparse") {
  MapSparseMatrixXd Z = MapSparseMatrixXd(as<MapSparseMatrixXd> (design_matrix));
  predicted = Z * (mean_estimates / effiter); 
  } else {
    MapMatrixXd Z = MapMatrixXd(as<MapMatrixXd> (design_matrix));
    predicted = Z * (mean_estimates / effiter); 
    }
 
  return predicted;

}


Rcpp::List effects::get_summary(int effiter){

  Rcpp::List out =  Rcpp::List::create(Rcpp::Named("type") = sparse_or_dense,
			  	       Rcpp::Named("method") = method,
			               Rcpp::Named("scale_prior") = scale,
			               Rcpp::Named("df_prior") = df);
  out["posterior"] = my_functions->summary(mean_estimates,mean_var,effiter, var_posterior);

  return out;

}


