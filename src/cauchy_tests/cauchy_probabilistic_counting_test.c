/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "cauchy_tests/cauchy_tests.h"

/* version for separated input disk:
 * D(center, radius) and D(center, 4*radius) 
 * are assumed to contain the same number of roots */
cauchyTest_res cauchyTest_probabilistic_counting_test( const compRat_t center,
                                                       const realRat_t radius,           
                                                       cacheApp_t cache,
                                                       cacheCauchy_t cacheCau,
                                                       slong prec,
                                                       metadatas_t meta, int depth){
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
//     realRat_t RadInf, RadSup, Rad, ExRad, ratio;
//     realApp_t nbCentersApp;
//     realRat_init(RadInf);
//     realRat_init(RadSup);
//     realRat_init(Rad);
//     realRat_init(ExRad);
//     realRat_init(ratio);
//     realApp_init(nbCentersApp);
    
//     realRat_set(RadInf, radius);
//     realRat_mul_si(RadSup, RadInf, 4);
//     realRat_add(Rad, RadInf, RadSup);
//     realRat_div_ui(Rad, Rad, 2);
//     realRat_sub(ExRad, RadSup, RadInf);
//     realRat_div_ui(ExRad, ExRad, 2);
//     realRat_div(ExRad, ExRad, cacheCauchy_isoRatioref(cacheCau));
//     realRat_div(ratio, Rad, ExRad);
//     realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
//     slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    
    cacheCauchy_set_bounds( cacheCau, radius, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#---cauchy probabilistic counting test: \n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
    compApp_t s0;
    compApp_init(s0);
    int alreadyEvaluated = 0;
    
    /* apply counting test to D(center, 2*radius) */
    res = cauchyTest_computeS0Approx(s0, center, radius,
                                     NULL, 0, 0, 
                                    0, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_UNCERTIFI, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
    
    res.nbOfSol = ( res.nbOfSol==-2? -1:0 );
    
    if (res.nbOfSol==0) {
        realApp_add_error( compApp_realref(s0), wP + 0 );
        realApp_add_error( compApp_imagref(s0), wP + 0 );
//         res.nbOfSol = ( compApp_contains_zero(s0)==1? 0:1 );
        slong nbOfSol = -1;
        int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
        int containsZero = realApp_contains_zero( compApp_imagref(s0) );
        if (! ( unique && containsZero) ){
            res.nbOfSol = -1;
//         printf("--- ps counting test returns -1; prec: %d \n", (int) res.appPrec);
//         printf("------ s0 with 1/2 error: "); compApp_printd(s0, 20); printf("\n");
//         printf("------ unique: %d\n", unique);
//         printf("------ contains 0: %d\n", containsZero);
        }
        else {
            res.nbOfSol = (int) nbOfSol;
        }
    }
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#------ res: %i\n", res.nbOfSol);
        printf("#------ s0: "); compApp_printd(s0,10); 
        slong nbOfSol = -1;
        int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
        int containsZero = realApp_contains_zero( compApp_imagref(s0) );
        if (unique)
            printf(" real part contains a unique integer: YES, %d", (int) nbOfSol);
        else
            printf(" real part contains a unique integer: NO");
        if (containsZero)
            printf(" imaginary part contains zero: YES");
        else
            printf(" imaginary part contains zero: NO");
        printf("\n");
    }
    
    compApp_clear(s0);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
        metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
    }
    
    return res;
    
}

