/* ************************************************************************** */
/*  Copyright (C) 2018 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <stdio.h>
#include "compRat_poly.h"

int compRat_poly_fprint(FILE * file, const compRat_poly_t z){
    int r;
    r = fprintf(file, "real part: ");
    r = realRat_poly_fprint(file, compRat_poly_realref(z));
    r = fprintf(file, "\nimag part: ");
    r = realRat_poly_fprint(file, compRat_poly_imagref(z));
    return r;
}

int compRat_poly_fprint_pretty(FILE * file, const compRat_poly_t z, const char * var){
    int r;
    r = fprintf(file, "real part: ");
    r = realRat_poly_fprint_pretty(file, compRat_poly_realref(z), var);
    r = fprintf(file, "\nimag part: ");
    r = realRat_poly_fprint_pretty(file, compRat_poly_imagref(z), var);
    return r;
}