/*
// cGWAS.h
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

#ifndef _cpgen_cGWAS_H
#define _cpgen_cGWAS_H

#include <RcppEigen.h>
//#include <R.h>
//#include <Rmath.h>

#ifdef _OPENMP
  #define has_openmp 1
  #include <omp.h>
#else 
  #define has_openmp 0
  #define omp_get_num_threads() 1
  #define omp_set_num_threads(x) 1
  #define omp_get_max_threads() 1
  #define omp_get_thread_limit() 1
  #define omp_set_dynamic(x) 1
  #define omp_get_thread_num() 0
#endif

#include "printer.h"


using namespace Rcpp;
using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatrixXd;
typedef Eigen::MappedSparseMatrix<double> MapSparseMatrixXd;
typedef Eigen::Map<Eigen::VectorXd> MapVectorXd;
typedef Eigen::Map<Eigen::ArrayXd> MapArrayXd;

RcppExport SEXP cGWAS(SEXP yR,SEXP MR, SEXP VR, SEXP V2R, SEXP XR, SEXP domR, SEXP second_transformR, SEXP sparseR, SEXP verboseR, SEXP threadsR);


template<class T>
class GWA {
private:

SEXP V_R;
SEXP y_R;
SEXP X_R;
SEXP M_R;
SEXP V2_R;

public:

MatrixXd beta;
MatrixXd p_value;
MatrixXd se;

VectorXd ystar;
MatrixXd Xstar;


int p;
bool dom;
bool second_transform;
bool initialized;
bool verbose;
int number_effects;

void initialize(SEXP y_R_in, SEXP M_R_in, SEXP V_R_in, SEXP V2_R_in, SEXP X_R_in, bool has_dom, bool second_transform_in, bool verbose_R_in);
void run_GWAS();
Rcpp::List summary();

};

template<class T> void GWA<T>::initialize(SEXP y_R_in, SEXP M_R_in, SEXP V_R_in, SEXP V2_R_in, SEXP X_R_in, bool has_dom, bool second_transform_in, bool verbose_R_in){

V_R = V_R_in;
y_R = y_R_in;
M_R = M_R_in;
X_R = X_R_in;
V2_R = V2_R_in;

dom = has_dom;

T V = T(as<T> (V_R));
MapVectorXd y = MapVectorXd(as<MapVectorXd> (y_R));
MapMatrixXd X = MapMatrixXd(as<MapMatrixXd> (X_R));
MapMatrixXd M = MapMatrixXd(as<MapMatrixXd> (M_R));

p = M.cols();

ystar = V*y;
Xstar = V*X;

number_effects = (dom) ? 2 : 1;

beta = MatrixXd(p,number_effects);
p_value = MatrixXd(p,number_effects);
se = MatrixXd(p,number_effects);

second_transform=second_transform_in;
verbose=verbose_R_in;

initialized=1;

};


template<class T> void GWA<T>::run_GWAS(){

MapMatrixXd M = MapMatrixXd(as<MapMatrixXd> (M_R));
T V = T(as<T> (V_R));
MapMatrixXd V2 = MapMatrixXd(as<MapMatrixXd> (V2_R));


#pragma omp parallel
{
  MatrixXd xtx = MatrixXd(Xstar.cols()+number_effects,Xstar.cols()+number_effects);
  MatrixXd xtx_inv = MatrixXd(Xstar.cols()+number_effects,Xstar.cols()+number_effects);
  VectorXd xty = VectorXd(Xstar.cols()+number_effects);
  VectorXd e(ystar.size());
  VectorXd b(Xstar.cols()+number_effects);
  MatrixXd X_marker(ystar.size(),number_effects);
  MatrixXd xtx_temp = MatrixXd(number_effects,Xstar.cols());


  double df;
  double n = ystar.size();
  double res_se;

  df = n - Xstar.cols() - number_effects;

  xtx.block(0,0,Xstar.cols(),Xstar.cols()) = Xstar.transpose()*Xstar;
  xty.segment(0,Xstar.cols()) = Xstar.transpose()*ystar;

// for verbose
  int n_threads = omp_get_num_threads();
  int max = M.cols() / n_threads;
  printer prog(max);

  


#pragma omp for
  for(int i=0;i<M.cols();i++){

// this is for emmax - we need to perform the transformation: U'X on the markers before V2inv*Marker
    if(second_transform){
      X_marker.col(0) = V * (V2 * M.col(i));
      if(dom) {X_marker.col(1) = V * (V2 * (1.0 - abs(M.col(i).array())).matrix());}
    } else {
        X_marker.col(0) = V * M.col(i);
        if(dom) {X_marker.col(1) = V * (1.0 - abs(M.col(i).array())).matrix();}
      }

    xtx_temp = X_marker.transpose()*Xstar;
    xtx.block(Xstar.cols(),Xstar.cols(),number_effects,number_effects) = X_marker.transpose()*X_marker;
    xtx.block(Xstar.cols(),0,number_effects,Xstar.cols()) = xtx_temp;
    xtx.block(0,Xstar.cols(),Xstar.cols(),number_effects) = xtx_temp.transpose();
    xty.segment(Xstar.cols(),number_effects) = X_marker.transpose()*ystar;

    xtx_inv = xtx.inverse();
    b = xtx_inv*xty;
    e = ystar - Xstar*b.segment(0,Xstar.cols()) - X_marker*b.segment(Xstar.cols(),number_effects);
    res_se = e.squaredNorm()/df;
    se.row(i) = (xtx_inv.diagonal().segment(Xstar.cols(),number_effects).array() * res_se).sqrt();

    for(int j = 0;j<number_effects;j++){

      p_value(i,j) = R::pt(abs(b(Xstar.cols()+j))/se(i,j),df,0,0) * 2.0;
      beta(i,j) = b(Xstar.cols()+j);

      }  


    if(omp_get_thread_num()==0) {

      if (verbose) {

        prog.DoProgress();

      }

    }


//    if(omp_get_thread_num()==0) {
//        if(verbose) { 
//          progress++;
//          DoProgress(progress*n_threads, M.cols() );
//        }
//      }  

    }

  }

}

template<class T> Rcpp::List GWA<T>::summary(){

return Rcpp::List::create(Rcpp::Named("p_value") = p_value,
                          Rcpp::Named("beta") = beta,
			  Rcpp::Named("se") = se);
}


#endif











