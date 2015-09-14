
// Note: This is copied and modified form the sources of the 
// cran package: 'pedigreemm'

#include "pedigree_stuff.h"

/**
 * This is a creplacement for the former R-function "getGenAncestors"
 * It takes integer pointers to the sire, dam, id and generation vectors
 * and it changes the generations directly on the passed R-object, so there is no
 * return value = void.
 */

void calc_generation(int* sire, int* dam, int* id, int* gene, int this_id) {

/* j is the individual (label id) we are currently working with */
    int j = this_id;
/* integers to store the generations of dam and sire */    
    int gene_sire, gene_dam;
    
/* If there are no parents assign 0 */
    gene[j] = 0;

/* check first parent */ 
    if(sire[j] != NA_INTEGER) {

/* if the parent is not a founder continue */
      if (gene[sire[j]-1] == NA_INTEGER) {

/* we call this function again to get the generation of this parent - recursion */
/* The function doesnt return anything, it changes 'gene' directly */
/* we need the '-1' for the indices because R starts at 1, C at 0 */
        calc_generation(sire, dam, id, gene, sire[j]-1);
        gene_sire = 1 + gene[sire[j]-1];

      } else {
/* if the parent is a founder my generation is 1 */
          gene_sire = 1 + gene[sire[j]-1];      

        } 
/* we assign the generation derived from father for this individual */
/* but this might get overriden by the dam's derived generation if it is larger */
      gene[j] = gene_sire;

    }

/* check second parent */ 
/* all the same as above */
    if(dam[j] != NA_INTEGER) {

      if (gene[dam[j]-1] == NA_INTEGER) {

        calc_generation(sire, dam, id, gene, dam[j]-1);
        gene_dam = 1 + gene[dam[j]-1];

      } else {

          gene_dam = 1 + gene[dam[j]-1];      

        }
/* if the dam's derived generation is higher then override gene[j] with it */        
      if (gene_dam > gene[j])  gene[j] = gene_dam;  

    }

}

/**
 * This function iterates over every element of the pedigree and let
 * "calc_generation" calculate the generations. Again without a return value.
 */


void get_generation(SEXP sire_in, SEXP dam_in, SEXP id_in, SEXP gene_in, SEXP verbose_in) {

// get R-vectors for sire, dam, id and generation
    int *sire = INTEGER(sire_in),
	*dam = INTEGER(dam_in),
	*id = INTEGER(id_in),
	*gene = INTEGER(gene_in);
 
    int verbose = *INTEGER(verbose_in);
    int n = LENGTH(sire_in);

    for(int i=0; i<n; i++) {
      
      if(verbose) Rprintf("%i \n", i); 

// we only need to calculate the generation for non-founders, which
// is the case when we assigned 'NA' to generation (founders have 0)
      if(gene[i] == NA_INTEGER) {

// call the function above to get the generation for individual 'i'.
// we pass the vectors as pointers to calc_generation()
        calc_generation(sire, dam, id, gene, i);

      }

    }

}

