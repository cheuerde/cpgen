/*
// base_methods_abstract.h
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

#include <R_ext/Lapack.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;


class base_methods_abstract{

	public:
		virtual void initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns) = 0;
		virtual void initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns) = 0;
		virtual void sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, 
					    double * ycorr, VectorXd& var, double * var_e, 
					    sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1) = 0;
		virtual void sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, 
					    double * ycorr, VectorXd& var, double * var_e, 
					    sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1) = 0;
		virtual ~base_methods_abstract(){};

};
