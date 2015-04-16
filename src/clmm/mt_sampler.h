/*
// mt_sampler.h
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

#include <iostream>
#include <random>

// source: http://www.cplusplus.com/reference/random/mersenne_twister_engine/seed/

class sampler {

private:

  std::mt19937 gen;

public:
  double rnorm(double mean, double sd){

    return std::normal_distribution<double>(mean,sd)(gen);

  };

  double rchisq(double df){

    return std::chi_squared_distribution<double>(df)(gen);

  };

  int rbinom(size_t trials, double prob){

    return std::binomial_distribution<int>(trials, prob)(gen);

  };



  void set_seed(std::string string_seed){

  std::seed_seq my_seed(string_seed.begin(),string_seed.end());
  gen.seed(my_seed);

  };


// not nice
  void set_seed(){

  std::random_device rd;
  gen.seed(static_cast<uint32_t> (rd())); 

  };

  void check_sampler(){ std::cout << std::endl << " I am mt19937_C++ " << std::endl;};

  sampler() {

  std::random_device rd;
  gen = std::mt19937(static_cast<uint32_t> (rd())); 

  }

  sampler(std::string string_seed) {

  std::seed_seq my_seed(string_seed.begin(),string_seed.end());
 // cout << endl << my_seed << endl;
  gen = std::mt19937(my_seed);

  }

//  ~sampler(){ delete gen; }


};

