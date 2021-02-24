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
#include "polynomials/realRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include "flint/flint.h"

int main() {
        
    compApp_poly_t test;
    compApp_poly_init(test);
    compApp_poly_clear(test);
    
    realRat_poly_t poly;
    realRat_poly_t poly2;
    realRat_poly_init(poly);
    realRat_poly_init(poly2);
    realRat_poly_fit_length(poly,10);
    
    printf("length: %d\n", realRat_poly_length(poly));
    printf("degree: %d\n", realRat_poly_degree(poly));
    
    realRat_t r;
    realRat_init(r);
    realRat_set_si(r, 1, 3);
    realRat_poly_set_coeff_realRat(poly, 9, r);
    
//     realRat_poly_print(poly);
//     printf("\n");
    
    
    realRat_poly_set(poly2, poly);
    
    realRat_set_si(r, 1, 2);
    realRat_poly_set_coeff_realRat(poly, 8, r);
    
    realRat_poly_print_pretty(poly, "x");
    printf("\n");
    realRat_poly_print_pretty(poly2, "x");
    printf("\n");
    
    realRat_poly_get_coeff_realRat(r, poly, 8);
    realRat_print(r);
    printf("\n");
    
//     arith_bernoulli_polynomial(poly , 10);
    bernoulli_polynomial(poly , 10);
    printf("10-th Bernoulli polynomial: ");
    realRat_poly_print_pretty(poly, "x");
    printf("\n");
    
    int deg = 64;
    int bitsize = 64;
    realRat_poly_init(poly2);
    mignotte_polynomial(poly2 , deg, (int)bitsize/2-1);
    printf("Mignotte polynomial: ");
    realRat_poly_print_pretty(poly2, "x");
    printf("\n");
    
//     realRat_clear(coeff);
    realRat_clear(r);
    realRat_poly_clear(poly);
    realRat_poly_clear(poly2);
    return 0;   
}
