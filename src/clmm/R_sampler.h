/*
// R_sampler.h
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

#include <R.h>
#include <Rmath.h>

class sampler {

public:
double  rnorm(double mean, double sd){

  return R::rnorm(mean,sd);

  };

double  rchisq(int df){

  return R::rchisq(df);

  };

  void set_seed(uint32_t seed){};

  void set_seed(){};

void check_sampler(){ 

  std::cout << std::endl << " I am the R-RNG " << std::endl;

  };

};
