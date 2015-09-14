/*
// pedigree_stuff.h
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

// Note: This is copied and modified form the sources of the cran
// package: 'pedigreemm'

#ifndef pedigree_stuff_H
#define pedigree_stuff_H

#include <R.h>
#include <Rdefines.h>

void get_generation(SEXP sire_in, SEXP dam_in, SEXP id_in, SEXP gene_in, SEXP verbose_in);
void calc_generation(int* sire, int* dam, int* id, int* gene, int this_id);

#endif 
