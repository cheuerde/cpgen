## Package update submission
* besides functionality and bug fixes I made changes according
to notification emails from CRAN maintainers.
Notably, "package=" was changed to "PACKAGE=" at one instance
as well as a CXX1X reference in a commented make file was removed.
 
## Test environments
* local ubuntu 16.04, R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 4 NOTEs:

checking top-level files ... NOTE
File
  LICENSE
is not mentioned in the DESCRIPTION file.

checking Rd line widths ... NOTE
Rd file 'clmm.Rd':
  \usage lines wider than 90 characters:
     verbose = TRUE, timings = FALSE, seed = NULL, use_BLAS=FALSE, beta_posterior_fixed = FALSE)

These lines will be truncated in the PDF manual.

checking compiled code ... NOTE
File ‘cpgen/libs/cpgen.so’:
  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

It is good practice to register native routines and to disable symbol
search.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.


* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.
