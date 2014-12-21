/*
// function_base.h
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


#include "base_methods.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

/////////////////////////////
// Abstract function class //
/////////////////////////////

class function_base {
public:
virtual void initialize(base_methods_abstract*& base_fun, VectorXd& xtx, int& columns, VectorXd& var, VectorXd& mean_var, double& var_y, bool& initialized, int total_number_effects) = 0;
virtual void sample_effects(base_methods_abstract*& base_fun,VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec) = 0;
virtual void sample_variance(VectorXd& var, double& scale, double& df, VectorXd& estimates, int& columns, sampler& mcmc_sampler, vector<double>& var_posterior,int niter) = 0;
virtual void update_means(VectorXd& mean_estimates, VectorXd& estimates, VectorXd& mean_var, VectorXd& var) = 0; 
virtual Eigen::VectorXd predict(VectorXd& mean_estimates, int start, int length) = 0; 
virtual Rcpp::List summary(VectorXd& mean_estimates, VectorXd& mean_var, int effiter, vector<double>& var_posterior) = 0;
virtual ~function_base(){};

};



