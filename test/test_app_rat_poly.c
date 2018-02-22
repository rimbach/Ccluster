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
#include "polynomials/app_rat_poly.h"

int main() {
    compApp_poly_t p;
    
    realRat_poly_t pbern, pmign;
    realRat_poly_init(pbern);
    realRat_poly_init(pmign);
    bernoulli_polynomial(pbern, 5);
    mignotte_polynomial(pmign, 5, (int)64/2-1);
    
    printf("--pbern: \n"); realRat_poly_print_pretty(pbern, "x"); printf("\n\n");
    printf("--pmign: \n"); realRat_poly_print_pretty(pmign, "x"); printf("\n\n");
    
    compApp_poly_init(p);
    compApp_poly_set_realRat_poly( p, pbern, 5);
    printf("p bern: "); compApp_poly_printd(p, 10); printf("\n\n");
    compApp_poly_clear(p);
    compApp_poly_init(p);
    compApp_poly_set_realRat_poly( p, pmign, 5);
    printf("p mign: "); compApp_poly_printd(p, 10); printf("\n\n");
    compApp_poly_clear(p);
    
    compRat_poly_t pRat;
    compRat_poly_init(pRat);
    compRat_poly_set2_realRat_poly(pRat,pmign,pbern);
    compApp_poly_init(p);
    compApp_poly_set_compRat_poly(p, pRat, 10);
    printf("p : "); compApp_poly_printd(p, 10); printf("\n\n");
    
    realRat_t q;
    realRat_init(q);
    realRat_set_si(q,1,2);
    compApp_poly_scale_realRat_in_place( p->coeffs, q, p->length, 10 );
    printf("p scaled by 1/2: "); compApp_poly_printd(p, 10); printf("\n\n");
    realRat_clear(q);
    
    compRat_t center;
    realRat_t radius;
    compRat_init(center);
    realRat_init(radius);
    compRat_set_sisi(center, 1,1,1,1);
    realRat_set_si(radius, 1,1);
    compApp_poly_taylorShift_in_place( p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, 10);
    printf("p shifted in 1+i with radius 1: "); compApp_poly_printd(p, 10); printf("\n\n");
    compRat_clear(center);
    realRat_clear(radius);
    
    compApp_poly_clear(p);
    realRat_poly_clear(pbern);
    realRat_poly_clear(pmign);
    compRat_poly_clear(pRat);
    return 0;
    
}