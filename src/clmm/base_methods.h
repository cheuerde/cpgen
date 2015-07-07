/*
// base_methods.h
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

#include "base_methods_abstract.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

class base_methods_st: public base_methods_abstract {
public:
inline void initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns);
inline void initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns);
inline void sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);
inline void sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);

};


void base_methods_st::initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
  xtx = Z.colwise().squaredNorm();
  columns = Z.cols();

}


void base_methods_st::initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
// the loop is necessary as Eigen::Sparse doesnt support colwise() 
  for(int i=0;i<Z.cols();i++) { xtx(i) = Z.col(i).squaredNorm(); }
  columns = Z.cols();

};


void base_methods_st::sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){


  double b_temp,rhs,lhs,inv_lhs,mean;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    b_temp=estimates(i);
// this saves a lot of computation time, as we don't have to adjust ycorr
// Taken from Rohan Fernando's BayesC implementation
    rhs = Z.col(i).dot(ycorr_map);
    rhs += xtx(i) * b_temp;
// here we can have fixed or random effects behaviour by doing shrinkage or not (lambda)
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);
    ycorr_map += Z.col(i) * b_adj;

  }

}





// sparse specialization
void base_methods_st::sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){

  double b_temp,rhs,lhs,inv_lhs,mean;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    b_temp=estimates(i);
    rhs=0;

// Iterate over the non-zero values in the i-th column of Z
    for (InIterMat it_(Z, i); it_; ++it_){

      rhs += it_.value() * ycorr_map(it_.index()); 

    }

    rhs += xtx(i) * b_temp;
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);

    for (InIterMat it_(Z, i); it_; ++it_){

      ycorr_map(it_.index()) += it_.value() * b_adj; 

    }

  }

}



///////////////////////
/// Parallel /////////
//////////////////////


class base_methods_mp: public base_methods_abstract {
public:
inline void initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns);
inline void initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns);
inline void sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);
inline void sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);

};


void base_methods_mp::initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
#pragma omp parallel for
  for(int i=0;i<Z.cols();i++) { xtx(i) = Z.col(i).squaredNorm(); }
  columns = Z.cols();

}


void base_methods_mp::initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
// the loop is necessary as Eigen::Sparse doesnt support colwise() 
#pragma omp parallel for
  for(int i=0;i<Z.cols();i++) { xtx(i) = Z.col(i).squaredNorm(); }
  columns = Z.cols();

};


void base_methods_mp::sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){


  double b_temp,rhs,lhs,inv_lhs,mean,sum_mp;
  int who;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    sum_mp=0;
    b_temp=estimates(i);

#pragma omp parallel reduction(+:sum_mp) private(who)
{
      who=omp_get_thread_num();
      MapVectorXd x_mp(MapVectorXd(&Z(thread_vec.at(who)["start"],i),thread_vec.at(who)["length"],1));
      MapVectorXd y_mp(MapVectorXd(&ycorr_map(thread_vec.at(who)["start"]),thread_vec.at(who)["length"])); 
      sum_mp+=x_mp.dot(y_mp);
}

    rhs = sum_mp + xtx(i) * b_temp;
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);

#pragma omp parallel  private(who)
{
      who=omp_get_thread_num();
      MapVectorXd x_mp(MapVectorXd(&Z(thread_vec.at(who)["start"],i),thread_vec.at(who)["length"],1));
      MapVectorXd y_mp(MapVectorXd(&ycorr_map(thread_vec.at(who)["start"]),thread_vec.at(who)["length"])); 
      y_mp += x_mp * b_adj;
}


  }

}



// sparse specialization
void base_methods_mp::sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){

  double b_temp,rhs,lhs,inv_lhs,mean;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    b_temp=estimates(i);
    rhs=0;

// Iterate over the non-zero values in the i-th column of Z
    for (InIterMat it_(Z, i); it_; ++it_){

      rhs += it_.value() * ycorr_map(it_.index()); 

    }

    rhs += xtx(i) * b_temp;
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);

    for (InIterMat it_(Z, i); it_; ++it_){

      ycorr_map(it_.index()) += it_.value() * b_adj; 

    }

  }

}


/////////////////
/// For BLAS ////
/////////////////


class base_methods_BLAS: public base_methods_abstract {
public:
inline void initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns);
inline void initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns);
inline void sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);
inline void sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler,mp_container& thread_vec, double lambda = 1);

};


void base_methods_BLAS::initialize(MapMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
  xtx = Z.colwise().squaredNorm();
  columns = Z.cols();

}


void base_methods_BLAS::initialize(MapSparseMatrixXd& Z, VectorXd& xtx, int& columns){

  xtx = VectorXd(Z.cols());
// the loop is necessary as Eigen::Sparse doesnt support colwise() 
  for(int i=0;i<Z.cols();i++) { xtx(i) = Z.col(i).squaredNorm(); }
  columns = Z.cols();

};


void base_methods_BLAS::sample_effects(MapMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){


  double b_temp,rhs,lhs,inv_lhs,mean;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  int obs = ycorr_map.size();
  int incr = 1;
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    b_temp=estimates(i);
    rhs=F77_NAME(ddot)(&obs,&Z(0,i),&incr, ycorr_map.data(),&incr);
    rhs += xtx(i) * b_temp;
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);
    F77_NAME(daxpy)(&obs, &b_adj,&Z(0,i),&incr, ycorr_map.data(),&incr);

  }

}




// sparse specialization
void base_methods_BLAS::sample_effects(MapSparseMatrixXd& Z, VectorXd& xtx, VectorXd& estimates, double * ycorr, VectorXd& var, double * var_e, sampler& mcmc_sampler, mp_container& thread_vec, double lambda){

  double b_temp,rhs,lhs,inv_lhs,mean;
  MapVectorXd ycorr_map(ycorr,Z.rows());
  double b_adj;

  for(int i=0;i<Z.cols();i++) { 

    b_temp=estimates(i);
    rhs=0;

// Iterate over the non-zero values in the i-th column of Z
    for (InIterMat it_(Z, i); it_; ++it_){

      rhs += it_.value() * ycorr_map(it_.index()); 

    }

    rhs += xtx(i) * b_temp;
    lhs = xtx(i) + (*var_e / var(i) * lambda);
    inv_lhs = 1.0 / lhs;
    mean = inv_lhs*rhs;
    estimates(i) = mcmc_sampler.rnorm(mean,sqrt(inv_lhs * *var_e));
    b_adj = b_temp - estimates(i);

    for (InIterMat it_(Z, i); it_; ++it_){

      ycorr_map(it_.index()) += it_.value() * b_adj; 

    }

  }

}


