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

cauchyTest_res cauchyTest_probabilistic_counting( const compDsk_t Delta,         
                                                  cacheApp_t cache,
                                                  cacheCauchy_t cacheCau,
                                                  slong prec,
                                                  metadatas_t meta, int depth){
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    cacheCauchy_set_bounds( cacheCau, compDsk_radiusref(Delta), CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#---cauchy probabilistic counting: \n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
    compApp_t s0;
    compApp_init(s0);
    int alreadyEvaluated = 0;
    
    res = cauchyTest_computeS0Approx(s0, compDsk_centerref(Delta), compDsk_radiusref(Delta),
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

