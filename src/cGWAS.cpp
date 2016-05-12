/*
// cGWAS.cpp
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

#include "cGWAS.h"

SEXP cGWAS(SEXP yR,SEXP MR, SEXP VR, SEXP V2R, SEXP XR, SEXP domR, SEXP second_transformR, SEXP sparseR, SEXP verboseR, SEXP threadsR) {

	bool dom = as<bool> (domR); 
	bool sparse = as<bool> (sparseR);
	bool second_transform = as<bool> (second_transformR);
	bool verbose = as<bool> (verboseR);
	int threads = as<int>(threadsR); 

	omp_set_num_threads(threads);

	Eigen::setNbThreads(1);
	Eigen::initParallel();

	GWA<MapSparseMatrixXd> W_sparse;
	GWA<MapMatrixXd> W_dense;

	if(sparse) {
		 
		W_sparse.initialize(yR,MR,VR,V2R,XR,dom, second_transform,verbose);
		W_sparse.run_GWAS();
		return W_sparse.summary(); 

	} else {

		W_dense.initialize(yR,MR,VR,V2R,XR,dom, second_transform,verbose);
		W_dense.run_GWAS();
		return W_dense.summary(); 

	}

}

