/* ************************************************************************** */
/*  Copyright (C) 2019 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <stdio.h>
#include "numbers/realApp.h"
#include "polynomials/realRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "caches/cacheApp.h"
#include "rootRadii/realIntRootRadii.h"

int main() {
    
        realApp_t absPi;
        realApp_t absPj;
        realApp_t absPk;
        realApp_t rad;
        
        realApp_init(absPi);
        realApp_init(absPj);
        realApp_init(absPk);
        realApp_init(rad);
        
        slong i = 1;
        slong j = 2;
        slong k = 3;
        realApp_set_d(absPi, 1.);
        realApp_set_d(absPj, 2.);
        realApp_set_d(absPk, 3.);
        realApp_set_d(rad, 0.1);
        realApp_add_error(absPi, rad);
        realApp_add_error(absPj, rad);
        realApp_add_error(absPk, rad);
        
        int res = realIntRootRadii_liesBelow( i, absPi, j, absPj, k, absPk, 53 );
        printf( "res: %d \n", res);
        
        realApp_set_d(absPi, 2.);
        realApp_set_d(absPj, 4.);
        realApp_set_d(absPk, 8.);
        realApp_set_d(rad, 0.001);
        realApp_add_error(absPi, rad);
        realApp_add_error(absPj, rad);
        realApp_add_error(absPk, rad);
        
        res = realIntRootRadii_liesBelow( i, absPi, j, absPj, k, absPk, 53 );
        printf( "res: %d \n", res);
        
        realApp_set_d(absPi, 2.);
        realApp_set_d(absPj, 3.5);
        realApp_set_d(absPk, 8.);
        realApp_set_d(rad, 0.001);
        realApp_add_error(absPi, rad);
        realApp_add_error(absPj, rad);
        realApp_add_error(absPk, rad);
        
        res = realIntRootRadii_liesBelow( i, absPi, j, absPj, k, absPk, 53 );
        printf( "res: %d \n", res);
        
        realApp_set_d(absPi, 2.);
        realApp_set_d(absPj, 4.5);
        realApp_set_d(absPk, 8.);
        realApp_set_d(rad, 0.001);
        realApp_add_error(absPi, rad);
        realApp_add_error(absPj, rad);
        realApp_add_error(absPk, rad);
        
        res = realIntRootRadii_liesBelow( i, absPi, j, absPj, k, absPk, 53 );
        printf( "res: %d \n", res);
        
        i=1;
        j=3;
        k=5;
        realApp_set_d(absPi, 2.);
        realApp_set_d(absPj, 8.);
        realApp_set_d(absPk, 32.);
        realApp_set_d(rad, 1e-10);
        realApp_add_error(absPi, rad);
        realApp_add_error(absPj, rad);
        realApp_add_error(absPk, rad);
        
        res = realIntRootRadii_liesBelow( i, absPi, j, absPj, k, absPk, 53 );
        printf( "res: %d \n", res);
        
        
        
        realApp_set_d(rad, 1e-10);
        
        slong len = 4;
        realApp_ptr absCoeffs = (realApp_ptr) ccluster_malloc ( (len)*sizeof(realApp) );
        for (slong ind = 0; ind<len; ind ++)
            realApp_init(absCoeffs + ind);
        realApp_set_d(absCoeffs + 0, 2.);   realApp_add_error(absCoeffs + 0, rad);
        realApp_set_d(absCoeffs + 1, 5.);   realApp_add_error(absCoeffs + 1, rad);
        realApp_set_d(absCoeffs + 2, 8.);   realApp_add_error(absCoeffs + 2, rad);
        realApp_set_d(absCoeffs + 3, 15.);  realApp_add_error(absCoeffs + 3, rad);
        
        printf(" absCoeffs[0]: "); realApp_printd( absCoeffs + 0, 10); printf("\n");
        printf(" absCoeffs[1]: "); realApp_printd( absCoeffs + 1, 10); printf("\n");
        printf(" absCoeffs[2]: "); realApp_printd( absCoeffs + 2, 10); printf("\n");
        printf(" absCoeffs[3]: "); realApp_printd( absCoeffs + 3, 10); printf("\n");
        
//         res = realRootRadii_liesBelow( 0, absCoeffs + 0, 1, absCoeffs + 1, 2, absCoeffs + 2, 53 );
//         printf( "res: %d \n", res);
        
        slong lenCh = 0;
        slong * convexHull = (slong *) ccluster_malloc ( (len)*sizeof(slong) );
        
        lenCh = realIntRootRadii_convexHull( convexHull, absCoeffs, len, 53 );
        
        printf(" Convex hull: %ld vertices: ", lenCh );
        for (slong ind = 0; ind < lenCh; ind++)
            printf("%ld, ", convexHull[ind]);
        printf("\n");
    
        realApp_clear(absPi);
        realApp_clear(absPj);
        realApp_clear(absPk);
        realApp_clear(rad);
        
        for (slong ind = 0; ind<len; ind ++)
            realApp_clear(absCoeffs + ind);
        ccluster_free(absCoeffs);
        ccluster_free(convexHull);
        
        slong degree = 512;
        realRat_poly_t pmig;
        realRat_poly_init(pmig);
        mignotte_polynomial(pmig , degree, 14);
//         wilkinson_polynomial(pmig , degree);
//         bernoulliInt_polynomial(pmig , degree);
        cacheApp_t cache;
        cacheApp_init_realRat_poly(cache, pmig);
        
        realRat_t delta;
        realRat_init(delta);
        realRat_set_si(delta, 1, degree*degree);
//         realRat_set_si(delta, 1, degree);
        
        compAnn_list_t annulii;
        compAnn_list_init(annulii);
        
        slong prec = realIntRootRadii_rootRadii( annulii, cache, delta );
        
        compAnn_list_printd(annulii, 10);
        printf("\n\n");
        
        realIntRootRadii_connectedComponents( annulii, prec );
        
        compAnn_list_printd(annulii, 10);
        printf("\n\n");
        
        realIntRootRadii_containsRealRoot( annulii, cache, prec );
        
        compAnn_list_printd(annulii, 10);
        printf("\n\n");
        
        compAnn_list_clear(annulii);
        
        realRat_clear(delta);
        realRat_poly_clear(pmig);
        cacheApp_clear(cache);
        
        return 0;
    
}
