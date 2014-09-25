/*
// ctools.cpp
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


#include "ctools.h"

// check if openmp is available

SEXP check_openmp() { return wrap(has_openmp); }


// check openmp_version

SEXP check_openmp_version() { return wrap(OMP_VERSION); }


// check_max_threads

SEXP get_max_threads(){ return wrap(omp_get_max_threads()); } 

// set_num_threads

SEXP set_num_threads(SEXP n){ omp_set_num_threads(as<int>(n)); } 

// check_limit_threads

// SEXP get_limit_threads(){ return wrap(omp_get_thread_limit()); } 



// ccrossproduct


SEXP ccp_dense_dense(SEXP XR, SEXP ZR, SEXP threadsR){

  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  return(ccrossproduct<MapMatrixXd,MapMatrixXd>(XR,ZR));

}

SEXP ccp_dense_sparse(SEXP XR, SEXP ZR){


  return(ccrossproduct<MapMatrixXd,MapSparseMatrixXd>(XR,ZR));

}

SEXP ccp_sparse_dense(SEXP XR, SEXP ZR){


  return(ccrossproduct<MapSparseMatrixXd,MapMatrixXd>(XR,ZR));

}

SEXP ccp_sparse_sparse(SEXP XR, SEXP ZR){


  MapSparseMatrixXd X(as<MapSparseMatrixXd> (XR));
  MapSparseMatrixXd Z(as<MapSparseMatrixXd> (ZR));

// On this occasion I use wrap, as it is convenient
// and memory might not matter that much here
  return(wrap(X*Z));

}


// ctrace

SEXP ctrace(SEXP XR){
 
  MapMatrixXd X(as<MapMatrixXd> (XR));
  return wrap(X.trace());

}


// ccov

SEXP ccov(SEXP Xa,SEXP lambdaR, SEXP wR, SEXP corR, SEXP threadsR)

{

  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  const Map<MatrixXd> X(as<Map<MatrixXd> >(Xa));
  double lambda = as<double>(lambdaR);
  int cor = as<int>(corR);
  const Map<VectorXd> w(as<Map<VectorXd> >(wR));


  RowVectorXd mu = (X*w.asDiagonal()).colwise().sum();
  VectorXd sd;
  int p = X.cols();
  Rcpp::NumericMatrix covX_out(p,p);
  Eigen::Map<Eigen::MatrixXd> covX(covX_out.begin(),p,p);

  covX = (1-lambda) * ((X.rowwise()-mu).transpose()*w.asDiagonal()*(X.rowwise()-mu)).array() / (1-(w.dot(w))); 

  if(cor!=0) 
  
  { 

    sd = 1 / covX.diagonal().array().sqrt();
    covX = sd.asDiagonal()*covX*sd.asDiagonal();

  }
    
  covX.diagonal().array() += lambda;

  return covX_out;  
    

}


// cscanx


SEXP cscanx(SEXP path){

  int row=0,col=0;


  std::string line;

  std::string c_path = as<std::string>(path);
  std::ifstream f(c_path.c_str());
  if (!f.is_open()) { perror(("error while opening file " + c_path).c_str()); }



  std::vector<double> X;
  std::string n;


  while(getline(f, line)) 

	{
	std::stringstream t(line);	
	if(row==0) 
		{ while (t >> n) 
			{ X.push_back(atof(n.c_str())); col++;};
	} else {
		while (t >> n)
			{

			X.push_back(atof(n.c_str())); 

			}
		}
	row++; 
	}	

  Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> >Z(X.data(),row,col);

  return wrap(Z);


}

// camat


SEXP camat(SEXP Xa,SEXP lambdaR, SEXP yangR, SEXP threadsR)

{

  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  const Map<MatrixXd> X(as<Map<MatrixXd> >(Xa));
  int yang = as<int>(yangR);
  double lambda = as<double>(lambdaR);
  double p = X.cols();
  int n = X.rows();

  RowVectorXd maf = (X.array()+1).colwise().sum() / (n*2);
  Map<ArrayXd> mafa(maf.data(),X.cols());

  Rcpp::NumericMatrix A_out(n,n);
  Eigen::Map<Eigen::MatrixXd> A(A_out.begin(),n,n);

 

  if(yang == 0) {

    double c = (maf.array()*(1-maf.array()).array()).sum() * 2;

    A.noalias() = (1-lambda) * ((((X.array()+1).matrix().rowwise()-2*maf)*((X.array()+1).matrix().rowwise()-2*maf).transpose()) / c).array();
    A.diagonal().array() += lambda;

  } else {

      RowVectorXd var = (2*maf).array()*(1-maf.array());
      A = (((X.array()+1).matrix().rowwise()-2*maf).matrix()*(1/var.array()).matrix().asDiagonal()
                   *((X.array()+1).matrix().rowwise()-2*maf).matrix().transpose()).array() / p;

// Diagonals


      ArrayXd x(X.cols());

      for (int i=0;i<n;i++) {

        x = X.row(i).array()+1;

        A(i,i) =  ((x*x - (1+2*mafa)*x + 2*(mafa*mafa)) / (2*mafa*(1-mafa))).sum() / p + 1;

      }

    }
    
  return A_out;  

}


// cdmat 


SEXP cdmat(SEXP Xa, SEXP lambdaR, SEXP threadsR){

 
  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  const Map<MatrixXd> X(as<Map<MatrixXd> >(Xa));

  double lambda = as<double>(lambdaR);

  int n = X.rows();
  Rcpp::NumericMatrix D_out(n,n);
  Eigen::Map<Eigen::MatrixXd> D(D_out.begin(),n,n);

// compute allele-frequencies
  RowVectorXd p = (X.array()+1).matrix().colwise().sum() / (n*2);

  RowVectorXd q = 1 - p.array();

  double c = (2 * p.array() * q.array() * (1 - 2 * p.array() * q.array())).sum();

  D.noalias() = (1-lambda) * ((((1 - X.array().abs()).matrix().rowwise() - (2*p.array()*q.array()).matrix())
                * ((1 - X.array().abs()).matrix().rowwise() - (2*p.array()*q.array()).matrix()).transpose()) / c);

  D.diagonal().array() += lambda;

  return D_out;


}


// cgrm



SEXP cgrm(SEXP XR,SEXP wR, SEXP iswR, SEXP lambdaR, SEXP threadsR){


  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  Map<MatrixXd> X(as<Map<MatrixXd> >(XR));
  double lambda = as<double>(lambdaR);
  bool isw = as<bool>(iswR);
  
  int n = X.rows();
  Rcpp::NumericMatrix G_out(n,n);
  Eigen::Map<Eigen::MatrixXd> G(G_out.begin(),n,n);

// colwise means and variances
  RowVectorXd mu = X.colwise().sum() / n;
  ArrayXd var = (X.rowwise() - mu).colwise().squaredNorm() / (n-1);
 
  if(isw) {
  	
    Map<ArrayXd> w(as<Map<ArrayXd> >(wR));
    G.noalias() = (1-lambda) * (((X.rowwise() - mu) * (w / var).matrix().asDiagonal() * (X.rowwise() - mu).transpose()) / w.sum() );
    
  } else {
  	
      G.noalias() = (1-lambda) * (((X.rowwise() - mu) * (X.rowwise() - mu).transpose()) / var.sum());
      
    }

  G.diagonal().array() += lambda;

  return G_out;

}


// ccross


SEXP ccross(SEXP Xa, SEXP Da, SEXP threadsR){


  int nt = as<int>(threadsR);
  omp_set_num_threads(nt);
  Eigen::setNbThreads(nt);

  const Map<MatrixXd> X(as<Map<MatrixXd> >(Xa));
  const Map<VectorXd> D(as<Map<VectorXd> >(Da));
  
  int n = X.rows();
  Rcpp::NumericMatrix A_out(n,n);
  Eigen::Map<Eigen::MatrixXd> A(A_out.begin(),n,n);

  A.noalias() = X*D.asDiagonal()*X.transpose();

  return A_out;

}



// csolve


SEXP csolve(SEXP XR, SEXP yR,SEXP threadsR) {

  return(eigensolver<MapMatrixXd, MapMatrixXd, Eigen::LLT<Eigen::MatrixXd> >(XR,yR,threadsR));

}

// csolve_sparse

SEXP csolve_sparse(SEXP XR, SEXP yR, SEXP threadsR) {

  return(eigensolver<MapSparseMatrixXd,MapMatrixXd, Eigen::SimplicialLLT<Eigen::SparseMatrix<double, Eigen::ColMajor> > >(XR,yR, threadsR));

}



// cmaf

SEXP cmaf(SEXP Xa){


  const Map<MatrixXd> X(as<Map<MatrixXd> >(Xa));

  double n = X.rows();

  const RowVectorXd maf = (X.array()+1).matrix().colwise().sum() / (n*2);

  return wrap(maf);

}


// ccolmv


SEXP ccolmv(SEXP XR,SEXP varR){


  Map<MatrixXd> X(as<Map<MatrixXd> >(XR));
  bool var = as<bool>(varR);

  double n = X.rows();

// colwise means and variances
  RowVectorXd mu = X.colwise().sum() / n;
  if(var) {
    RowVectorXd var = (X.rowwise() - mu).colwise().squaredNorm() / (n-1);
    return wrap(var); 
  } else {
    
      return wrap(mu);
	
    }
}



			

