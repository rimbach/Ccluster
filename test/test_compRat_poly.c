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
#include "polynomials/compRat_poly.h"

int main() {
    
    compRat_poly_t p;
    compRat_poly_init(p);
    
    compRat_poly_zero(p);
    printf("--p: \n");
    compRat_poly_print_pretty(p, "x");
    printf("\n");
    
    realRat_poly_t pbern, pmign;
    realRat_poly_init(pbern);
    realRat_poly_init(pmign);
    bernoulli_polynomial(pbern, 5);
    mignotte_polynomial(pmign, 10, (int)64/2-1);
    compRat_poly_set_realRat_poly(p,pbern);
    printf("--p: \n");
    compRat_poly_print_pretty(p, "x");
    printf("\n");
    
    compRat_poly_set2_realRat_poly(p,pmign,pbern);
    printf("--p: \n");
    compRat_poly_print_pretty(p, "x");
    printf("\n");
    
    realRat_poly_clear(pbern);
    realRat_poly_clear(pmign);
    compRat_poly_clear(p);
    return 0;
}