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

cauchyTest_res cauchyTest_deterministic_exclusion_test( const compRat_t center,
                                          const realRat_t radius,
                                         const realRat_t radius2,
                                         slong vangle,           
                                         slong vindex,           
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth){
    
    int level = 4;
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    cacheCauchy_set_bounds( cacheCau, radius, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cauchy deterministic exclusion test: \n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
    slong q1 = cacheCauchy_nbEvalExref(cacheCau);
//     slong q2 = cacheCauchy_nbEvalCeref(cacheCau);
    slong quo = cacheCauchy_quotientref(cacheCau);
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    realApp_ptr wP2 = cacheCauchy_wanErrCeref(cacheCau);
    
    slong nbPowerSums = cacheCauchy_nbPwSuExref(cacheCau);
    
//     if (metadatas_getVerbo(meta)>=3) {
// //         printf("------ number of eval points q1: %ld\n", q1);
//         printf("------ number of eval points q1: %ld\n", cacheCauchy_nbEvalExref(cacheCau));
// //         printf("------ number of eval points q2: %ld\n", q2);
//         printf("------ number of eval points q2: %ld\n", cacheCauchy_nbEvalCeref(cacheCau));
// //         printf("------ error wP1: "); realApp_printd(wP, 10); printf("\n");
//         printf("------ error wP1: "); realApp_printd(cacheCauchy_wanErrExref(cacheCau), 10); printf("\n");
// //         printf("------ error wP2: "); realApp_printd(wP2, 10); printf("\n");
//         printf("------ error wP2: "); realApp_printd(cacheCauchy_wanErrCeref(cacheCau), 10); printf("\n");
// //         printf("------ lowerBound: "); realApp_printd(lbApp, 10); printf("\n");
//         printf("------ lowerBound: "); realApp_printd(cacheCauchy_lBoundApref(cacheCau), 10); printf("\n");
// //         printf("------ upperBound: "); realApp_printd(ubApp, 10); printf("\n");
//         printf("------ upperBound: "); realApp_printd(cacheCauchy_uBoundApref(cacheCau), 10); printf("\n");
//     }
    
    compApp_ptr ps = (compApp_ptr) ccluster_malloc( nbPowerSums*sizeof(compApp) );
    for (int j=0; j<nbPowerSums; j++)
        compApp_init( ps +j );
    res = cauchyTest_computeSsApprox(ps, center, radius, 
                                     radius2, vangle, vindex, 
                                     cache, cacheCau, CAUCHYTEST_UNCERTIFI, prec, inCounting, meta, depth );
    
    res.nbOfSol = ( res.nbOfSol==-2? 1:0 );
    int j=0;
    while ( (j<nbPowerSums) && (res.nbOfSol==0) ) {
        realApp_add_error( compApp_realref(ps + j), wP + j );
        realApp_add_error( compApp_imagref(ps + j), wP + j );
        res.nbOfSol = ( compApp_contains_zero(ps+j)==1? 0:1 );
        j++;
    }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------ res: %i\n", res.nbOfSol);
        for (int j=0; j<nbPowerSums; j++) {
            printf("#------ s%d: ",j); compApp_printd(ps+j,10); 
            if (compApp_contains_zero(ps+j))
                printf(" contains zero: YES");
            else
                printf(" contains zero: NO");
            printf("\n");
        }
    }
    
    int alreadyEvaluated = 0;
    compApp_t s0;
    compApp_init(s0);
    
    if (res.nbOfSol==0) {
        
        alreadyEvaluated = 0;
        slong g = 0;
        while ((g<q1)&&(res.nbOfSol==0)) {
            
            res = cauchyTest_computeS0Approx(s0, center, radius,
                                             radius2, vangle, vindex,
                                             quo*g, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_CERTIFIED, res.appPrec, inCounting, meta, depth );
            
            
            /* test */
            realRat_t arg;
            compApp_t rot;
            realRat_init(arg);
            compApp_init(rot);
            slong q2 = cacheCauchy_nbEvalCeref(cacheCau);
            realRat_set_si (arg, -2*quo*g, q2);
            compApp_set_realRat(rot, arg, 2*prec);
            acb_exp_pi_i (rot,    rot,  2*prec);
            acb_mul(rot, rot, s0, 2*prec);
            printf("---g: %ld, s0: ", g); compApp_printd(s0, 10); printf("\n");
            printf("      rot s0: "); compApp_printd(rot, 10); printf("\n");
            realRat_clear(arg);
            compApp_clear(rot);
            
            res.nbOfSol = ( res.nbOfSol==-2? 1:0 );
            
            if (res.nbOfSol==0) {
                realApp_add_error( compApp_realref(s0), wP2 );
                realApp_add_error( compApp_imagref(s0), wP2 );
                res.nbOfSol = ( compApp_contains_zero(s0)==1? 0:1 );
            }
            
            g = g+1;
        }
        
        if (metadatas_getVerbo(meta)>=2) {
            if (res.nbOfSol!=0) {
                printf("#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
                printf("#------ res: %i, g:%ld\n", res.nbOfSol, g);
                printf("#------ s0: "); compApp_printd(s0,10); printf("\n");
            }
        }
        
    }
    
    compApp_clear(s0);
    
    for (int j=0; j<nbPowerSums; j++)
        compApp_clear( ps +j );
    ccluster_free(ps);
    
    if ( metadatas_haveToCount(meta) && (inCounting==0) ) {
        metadatas_add_time_CauExTo(meta, (double) (clock() - start));
        metadatas_add_CauchyExTest(meta, depth, res.appPrec);
    }
    
    return res;
}
