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

#include <progress.hpp>

using namespace Eigen;
using namespace Rcpp;
using namespace std;


class MCMC_abstract{

	public:
		virtual int get_niter()  = 0;
		virtual bool get_verbose() = 0;
		virtual void initialize() = 0;
		virtual void gibbs(Progress * prog) = 0;
		virtual void gibbs() = 0;
		virtual void summary() = 0;
		virtual std::string get_name() = 0;
		virtual Rcpp::List get_summary() = 0;
		virtual ~MCMC_abstract(){};

};
