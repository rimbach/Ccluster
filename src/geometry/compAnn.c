/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <stdio.h>
#include "geometry/compAnn.h"

void compAnn_fprintd( FILE * file, const compAnn_t x, slong digits ){
    fprintf(file, "#indMax: %ld, indMin: %ld, rrInPo: %d, rrInNe: %d \n", 
            compAnn_indMaxref(x), compAnn_indMinref(x), compAnn_rrInPoref(x), compAnn_rrInNeref(x) );
    if (compAnn_centerReref(x)!=0)
        fprintf(file, "#center: %ld, ", compAnn_centerReref(x) );
    if (compAnn_centerImref(x)!=0)
        fprintf(file, "#center: i%ld, ", compAnn_centerImref(x) );
    fprintf(file, "radInf: ");
    realApp_fprintd(file, compAnn_radInfref(x), digits);
    fprintf(file, "  radSup: ");
    realApp_fprintd(file, compAnn_radSupref(x), digits);
//     fprintf(file, "\n");
}

void compAnn_fprint( FILE * file, const compAnn_t x){
    fprintf(file, "#indMax: %ld, indMin: %ld, rrInPo: %d, rrInNe: %d \n", 
            compAnn_indMaxref(x), compAnn_indMinref(x), compAnn_rrInPoref(x), compAnn_rrInNeref(x) );
    if (compAnn_centerReref(x)!=0)
        fprintf(file, "#center: %ld, ", compAnn_centerReref(x) );
    if (compAnn_centerImref(x)!=0)
        fprintf(file, "#center: i%ld, ", compAnn_centerImref(x) );
    fprintf(file, "radInf: ");
    realApp_fprint(file, compAnn_radInfref(x));
    fprintf(file, "  radSup: ");
    realApp_fprint(file, compAnn_radSupref(x));
//     fprintf(file, "\n");
}

/* assume a1 and a2 are centered on the real line; a1 = c1, r1; a2 = c2, r2 */
/* assume c1 = 0 */
/* returns 0 only if a1 and a2 have no intersection */
/* otherwise returns the intersection that is imaginary positive */
int compAnn_intersect_realCenter( compApp_t intersection, const compAnn_t a1, const compAnn_t a2, slong prec){
    
    int res = 1;
    realApp_t r1, r2, c2, zero;
    
    realApp_init(r1);
    realApp_init(r2);
    realApp_init(c2);
    realApp_init(zero);
    realApp_zero(zero);
#define x compApp_realref(intersection)
#define y compApp_imagref(intersection)    
    
    realApp_set_si(c2,   compAnn_centerReref( a2 ) );
    realApp_set   (r1,  compAnn_radSupref( a1 ) );
    realApp_set   (r2,  compAnn_radSupref( a2 ) );
    realApp_union (r1, r1, compAnn_radInfref( a1 ),  prec  );
    realApp_union (r2, r2, compAnn_radInfref( a2 ),  prec  );
    
    /* compute x = ( c2^2 + r1^2 - r2^2 )/(2c2) */
    realApp_sqr   ( r1, r1, prec );
    realApp_sqr   ( r2, r2, prec );
    realApp_sqr   ( x,  c2, prec );
    realApp_add   ( x,  x,  r1  , prec );
    realApp_sub   ( x,  x,  r2  , prec );
    realApp_mul_si( c2, c2, 2,    prec );
    realApp_div   ( x,  x,  c2  , prec );
    /* compute y^2 = r1^2 - x^2 */
    /* put x^2 in c2 */
    realApp_sqr   ( c2, x, prec );
    realApp_sub   ( y,  r1, c2  , prec );
    /* check if y < 0; in that case, no intersection */
    if (realApp_lt( y, zero )==1) {
        res = 0;
    }
//     printf("# intersection before square root: ");
//     compApp_printd( intersection, 10 );
//     printf("\n");
    /* compute y = sqrt(r1^2 - x^2) */
    if (realApp_contains_zero (y) ==1 ) { /* center y in 0 */
//         realApp_union (y, zero, y,  prec  );
        realApp_add_error(zero, y);
        realApp_set(y, zero);
        mag_sqrt( arb_radref(y), arb_radref(y) );
    }
    else {
        realApp_sqrt  ( y,  y, prec );
    }
//     printf("# intersection after square root: ");
//     compApp_printd( intersection, 10 );
//     printf("\n");
    
    realApp_clear(r1);
    realApp_clear(r2);
    realApp_clear(c2);
    realApp_clear(zero);
    
#undef x
#undef y
    
    return res;
}
