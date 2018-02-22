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
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"

int main() {
    
    compApp_poly_t p;
    compApp_poly_init(p);
    compApp_poly_clear(p);
    
    compApp_poly_init2(p, 10);
    printf("p: "); compApp_poly_printd(p, 10); printf("\n");
    compApp_poly_clear(p);
    
    realRat_poly_t pbern, pmign;
    realRat_poly_init(pbern);
    realRat_poly_init(pmign);
    bernoulli_polynomial(pbern, 5);
    mignotte_polynomial(pmign, 5, (int)64/2-1);
    
    printf("--pbern: \n"); realRat_poly_print_pretty(pbern, "x"); printf("\n\n");
    printf("--pmign: \n"); realRat_poly_print_pretty(pmign, "x"); printf("\n\n");
    
    compApp_poly_init(p);
    compApp_poly_set_fmpq_poly( p, pbern, 5);
    printf("p bern: "); compApp_poly_printd(p, 10); printf("\n\n");
    compApp_poly_clear(p);
    compApp_poly_init(p);
    compApp_poly_set_fmpq_poly( p, pmign, 5);
    printf("p mign: "); compApp_poly_printd(p, 10); printf("\n\n");
    compApp_poly_clear(p);
    
    compApp_poly_init(p);
    compApp_poly_set2_fmpq_poly( p, pmign, pbern, 5);
    printf("p : "); compApp_poly_printd(p, 10); printf("\n\n");
    compApp_poly_clear(p);
        
    realRat_poly_clear(pbern);
    realRat_poly_clear(pmign);
    return 0;
}