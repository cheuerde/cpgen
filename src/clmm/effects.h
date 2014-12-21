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
typedef std::vector<std::map<std::string, int> > mp_container;


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
  int n_samples;
  int burnin;
  double scale;
  double df;
  double * var_e;
  VectorXd mean_var;
//  double var_scalar;
  VectorXd var;
  vector<double> var_posterior;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> effects_posterior;
  VectorXd xtx;
  double * ycorr;
//  MatrixXd mean_var;

  std::string sparse_or_dense;
  std::string method;
  std::string name;
  bool initialized;
  bool full_output;
  bool GWAS;
  int GWAS_window_size;
  int GWAS_n_windows;
  double GWAS_threshold;
// we use the mp_container for multithreading also for the windows in GWAS
  mp_container GWAS_window_container;
  Eigen::VectorXd GWAS_window_genetic_variance;
  Eigen::VectorXd GWAS_window_genetic_proportion;
  Eigen::VectorXd GWAS_window_bigger_threshold;
  Eigen::VectorXd genetic_values;
  double mean_gen_var;

  double var_y;
  VectorXd estimates;
  VectorXd mean_estimates;

  effects(SEXP design_matrix_from_MCMC, double * ycorr_from_MCMC, Rcpp::List list_from_MCMC, double * var_e_from_MCMC, double var_y_from_MCMC,int total_number_effects_from_MCMC,int niter_from_MCMC, int burnin_from_MCMC, bool full_output_from_MCMC);
  inline void initialize(base_methods_abstract*& base_fun);  
  inline void sample_effects(sampler& mcmc_sampler,base_methods_abstract*& base_fun, mp_container& thread_vec, int& niter_from_MCMC); 
  inline void sample_variance(sampler& mcmc_sampler, int& niter_from_MCMC); 
  inline void update_means();
  inline Eigen::VectorXd predict(int effiter);
  inline void update_GWAS();
  Rcpp::List get_summary(int effiter);
  std::string get_name() { return name; };

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
  name = as<std::string>(list_from_MCMC["name"]);
  GWAS = as<bool>(list_from_MCMC["GWAS"]);
  GWAS_threshold = as<double>(list_from_MCMC["GWAS_threshold"]);
  GWAS_window_size = as<int>(list_from_MCMC["GWAS_window_size"]);
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
    if(method == "ridge") my_functions = new function_random<MapSparseMatrixXd>(design_matrix); 
    if(method == "BayesA") my_functions = new function_bayesA<MapSparseMatrixXd>(design_matrix); 
  } else {
      if(method == "fixed") my_functions = new function_fixed<MapMatrixXd>(design_matrix); 
      if(method == "ridge") my_functions = new function_random<MapMatrixXd>(design_matrix); 
      if(method == "BayesA") my_functions = new function_bayesA<MapMatrixXd>(design_matrix); 
    }


  my_functions->initialize(base_fun,xtx,columns,var, mean_var,var_y,initialized,total_number_effects);
  estimates = VectorXd::Zero(columns);
  mean_estimates = VectorXd::Zero(columns);
  var_posterior.resize(niter);
// store posterior distribution of estimates
  if(full_output) effects_posterior = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(niter,columns);

// prepare GWAS stuff
  if(GWAS) {

    GWAS_n_windows = columns / GWAS_window_size;
    if(GWAS_n_windows < 1) GWAS_n_windows = 1;

// we use the mp_container for multithreading also for the windows in GWAS
    GWAS_window_container = mp_container(GWAS_n_windows);

    for(int i=0;i<GWAS_n_windows;i++) {

      GWAS_window_container.at(i)["start"] = i * columns / GWAS_n_windows;
      if(i==(GWAS_n_windows-1)){GWAS_window_container.at(i)["end"] = columns;} else {
      GWAS_window_container.at(i)["end"] = (i+1) * columns / GWAS_n_windows;}
      GWAS_window_container.at(i)["length"] = GWAS_window_container.at(i)["end"] - GWAS_window_container.at(i)["start"];

    }

    GWAS_window_genetic_variance = Eigen::VectorXd::Zero(GWAS_window_container.size());
    GWAS_window_genetic_proportion = Eigen::VectorXd::Zero(GWAS_window_container.size());
    GWAS_window_bigger_threshold = Eigen::VectorXd::Zero(GWAS_window_container.size());
    genetic_values = Eigen::VectorXd::Zero(columns);
    mean_gen_var = 0;

  }

  n_samples = 0;
  initialized = true;
  

}

void effects::sample_effects(sampler& mcmc_sampler,base_methods_abstract*& base_fun, mp_container& thread_vec, int& niter_from_MCMC){

  my_functions->sample_effects(base_fun,xtx,estimates, ycorr, var,var_e, mcmc_sampler,thread_vec);
  if(full_output) effects_posterior.row(niter_from_MCMC).noalias() = estimates;

}


void effects::sample_variance(sampler& mcmc_sampler, int& niter_from_MCMC){

  my_functions->sample_variance(var, scale, df, estimates, columns, mcmc_sampler, var_posterior, niter_from_MCMC);

}



void effects::update_means(){


  my_functions->update_means(mean_estimates,estimates,mean_var,var);

  if(GWAS) update_GWAS();

//  mean_estimates += estimates;
//  mean_var += var(0);
//  effiter++;

}


Eigen::VectorXd effects::predict(int effiter){

  Eigen::VectorXd temp = mean_estimates / effiter;
  return my_functions->predict(temp, 0, columns);

}


void effects::update_GWAS(){

//  auto get_var = [] (Eigen::Map<Eigen::VectorXd> vec) -> double { 

//    return (vec.array() - (vec.sum() / vec.size())).matrix.squaredNorm() / (vec.size() - 1);

//  };

// total genetic variance
  genetic_values.noalias() = my_functions->predict(estimates, 0, columns);
  double mu = genetic_values.sum() / genetic_values.size();
  double genetic_variance = (genetic_values.array() - mu).matrix().squaredNorm() / (genetic_values.size() - 1);
  mean_gen_var += genetic_variance;

  for(int i=0; i < GWAS_window_container.size(); i++) {

// variance for windows
    genetic_values.noalias() = my_functions->predict(estimates, GWAS_window_container.at(i)["start"], GWAS_window_container.at(i)["length"]);
    mu = genetic_values.sum() / genetic_values.size();
    double window_genetic_variance = (genetic_values.array() - mu).matrix().squaredNorm() / (genetic_values.size() - 1);
    GWAS_window_genetic_variance(i) += window_genetic_variance;
    double window_genetic_proportion = window_genetic_variance / genetic_variance;
    GWAS_window_genetic_proportion(i) += window_genetic_proportion;
    GWAS_window_bigger_threshold(i) += window_genetic_proportion > GWAS_threshold ? 1.0 : 0.0;

  }

}




Rcpp::List effects::get_summary(int effiter){

  Rcpp::List out =  Rcpp::List::create(Rcpp::Named("type") = sparse_or_dense,
			  	       Rcpp::Named("method") = method,
			               Rcpp::Named("scale_prior") = scale,
			               Rcpp::Named("df_prior") = df);
  Rcpp::List out_posterior = my_functions->summary(mean_estimates,mean_var,effiter, var_posterior);
  if(full_output) out_posterior["estimates"] = effects_posterior;
  out["posterior"] = out_posterior;

  if(GWAS) {

    Rcpp::List gwas_out;
    Eigen::VectorXi window(GWAS_window_container.size());
    Eigen::VectorXi start(GWAS_window_container.size());
    Eigen::VectorXi end(GWAS_window_container.size());

    for(int i=0; i < GWAS_window_container.size(); i++) {

      window(i) = i + 1;
      start(i) = GWAS_window_container.at(i)["start"] + 1;
      end(i) = GWAS_window_container.at(i)["end"];

    }

    gwas_out["window_size"] = GWAS_window_size;
    gwas_out["threshold"] = GWAS_threshold;
    gwas_out["mean_variance"] = mean_gen_var / effiter;
    gwas_out["window"] = window;
    gwas_out["start"] = start;
    gwas_out["end"] = end;
    gwas_out["window_variance"] = GWAS_window_genetic_variance / effiter;
    gwas_out["window_variance_proportion"] = GWAS_window_genetic_proportion / effiter;
    gwas_out["prob_window_var_bigger_threshold"] = GWAS_window_bigger_threshold / effiter;

    out["GWAS"] = gwas_out;

  }

  return out;

}


