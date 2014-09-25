/*
// ctools.h
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

#ifndef _cpgen_ctools_H
#define _cpgen_ctools_H

#include <RcppEigen.h>

#ifdef _OPENMP
  #define has_openmp 1
  #include <omp.h>
  #define OMP_VERSION _OPENMP
#else 
  #define has_openmp 0
  #define omp_get_num_threads() 1
  #define omp_set_num_threads(x) 1
  #define omp_get_max_threads() 1
  #define omp_get_thread_limit() 1
  #define omp_set_dynamic(x) 1
  #define omp_get_thread_num() 0
  #define OMP_VERSION 0
#endif

using namespace Rcpp;
using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatrixXd;
typedef Eigen::MappedSparseMatrix<double> MapSparseMatrixXd;
typedef Eigen::Map<Eigen::VectorXd> MapVectorXd;
typedef Eigen::Map<Eigen::ArrayXd> MapArrayXd;

template<class T1,class T2>
SEXP ccrossproduct(SEXP XR, SEXP ZR)

{

  T1 X(as<T1> (XR));
  T2 Z(as<T2> (ZR));

// allocate R-matrix
  int n = X.rows();
  int p = Z.cols();
  Rcpp::NumericMatrix sol(n,p);

// map that R-matrix to an Eigen-class
  MapMatrixXd sol_map(sol.begin(),n,p);

  sol_map.noalias() = X*Z;

  return sol;
  
};


// template for eigensolvers using dense
// and/or sparse matrices.
// return type is always a dense matrix here.
template<class T1, class T2, class T3>
SEXP eigensolver(SEXP XR, SEXP yR, SEXP threadsR) {

  int threads = as<int>(threadsR);
  omp_set_num_threads(threads);
  Eigen::setNbThreads(1);
  Eigen::initParallel();

  T1 X(as<T1> (XR));
  T2 y(as<T2> (yR));

  T3 W; 
  
  W.compute(X);

// allocate R-matrix
  int n = X.rows();
  int p = y.cols();
  Rcpp::NumericMatrix sol(n,p);

// map that R-matrix to an Eigen-class
  MapMatrixXd sol_map(sol.begin(),n,p);

// only use omp if threads > 1
if(threads >1) {

# pragma omp parallel for
    for(size_t i=0;i<p;++i) { 
      sol_map.col(i).noalias() = W.solve(y.col(i));
    } 
    
} else {
  
    sol.map.noalias() = W.solve(y);
    
  }
  
  return sol;

};



RcppExport SEXP check_openmp();
RcppExport SEXP check_openmp_version();
RcppExport SEXP get_max_threads();
RcppExport SEXP set_num_threads(SEXP n);
//RcppExport SEXP get_limit_threads();

RcppExport SEXP ccp_dense_dense(SEXP X, SEXP Z, SEXP threadsR);
RcppExport SEXP ccp_dense_sparse(SEXP X, SEXP Z);
RcppExport SEXP ccp_sparse_dense(SEXP X, SEXP Z);
RcppExport SEXP ccp_sparse_sparse(SEXP X, SEXP Z);
RcppExport SEXP ctrace(SEXP XR);

RcppExport SEXP ccov(SEXP Xa,SEXP lambdaR, SEXP wR, SEXP corR, SEXP threadsR);
RcppExport SEXP cscanx(SEXP path);
RcppExport SEXP camat(SEXP Xa,SEXP lambdaR, SEXP yangR, SEXP threadsR);
RcppExport SEXP cdmat(SEXP Xa, SEXP lambdaR, SEXP threadsR);
RcppExport SEXP cgrm(SEXP XR,SEXP wR, SEXP iswR, SEXP lambdaR, SEXP threadsR);
RcppExport SEXP ccross(SEXP Xa, SEXP Da, SEXP threadsR);
RcppExport SEXP csolve(SEXP XR, SEXP yR, SEXP threadsR);
RcppExport SEXP csolve_sparse(SEXP XR, SEXP yR, SEXP threadsR);
RcppExport SEXP cmaf(SEXP Xa);
RcppExport SEXP ccolmv(SEXP XR,SEXP varR);

#endif
