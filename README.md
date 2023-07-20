[![R build status](https://github.com/r-lib/actions/workflows/R-CMD-check/badge.svg)](https://github.com/cheuerde/cpgen/actions/actions?workflow=R-CMD-check)

cpgen
=====

## Parallel Genomic Evaluations in R

The package offers a variety of functions that are frequently being used in genomic prediction
and genomewide association studies. The package is based on`Rcpp` and`RcppEigen`, hence all routines
are implemented using the matrix algebra library `Eigen`.

## Parallelization

The main emphasis of the package lies in parallel computing which is realized by C++ functions making
use of shared memory multi-threading implemented through `OpenMP`. 
`R` is of single-threaded nature, hence almost the entire package is written
in C++. The package offers a variety of functions that lets you control and check
the number of threads that are being used by the functions of the package.
Internally every function uses the global variable `cpgen.threads` which is stored in
`options()$cpgen.threads`.
The value can be changed using the function `set_num_threads()`. When the package is loaded
in an R-session `cpgen.threads` will be set to the value returned by `get_max_threads()` which
is a wrapper for the OpenMP-header function `omp_get_max_threads()`

## Installation


```r
devtools::install_github("cheuerde/cpgen", ref = "master", build_vignettes=FALSE)
```
