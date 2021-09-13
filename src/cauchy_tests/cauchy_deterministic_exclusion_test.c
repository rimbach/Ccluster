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
    
//     slong q1 = cacheCauchy_nbEvalExref(cacheCau);
    slong q2 = cacheCauchy_nbEvalCeref(cacheCau);
//     slong quo = cacheCauchy_quotientref(cacheCau);
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
        while ((g<q2)&&(res.nbOfSol==0)) {
            
            res = cauchyTest_computeS0Approx(s0, center, radius,
                                             radius2, vangle, vindex,
                                             g, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_CERTIFIED, res.appPrec, inCounting, meta, depth );
            
            res.nbOfSol = ( res.nbOfSol==-2? 1:0 );
            
            if (res.nbOfSol==0) {
                realApp_add_error( compApp_realref(s0), wP2 );
                realApp_add_error( compApp_imagref(s0), wP2 );
                res.nbOfSol = ( compApp_contains_zero(s0)==1? 0:1 );
            }
            
//             if (g==0) {
//                 
//                 realApp_t div, max, mu;
//                 realApp_init(div);
//                 realApp_init(max);
//                 realApp_init(mu);
//                 // compute the max of fdivs 
//                 compApp_abs( max, cacheCauchy_fdivsCeref(cacheCau) + 0, res.appPrec);
//                 for (slong gg = 1; gg < q2; gg++) {
//                     compApp_abs( div, cacheCauchy_fdivsCeref(cacheCau) + gg, res.appPrec);
//                     realApp_max( max, max, div, res.appPrec);
//                 }
//                 // if the max of fdivs is strictly less than 1/(2q2 +1)
//                 realApp_one(mu);
//                 realApp_div_ui(mu, mu, 2*q2 + 1, res.appPrec);
//                 
//                 if (metadatas_getVerbo(meta)>=2) {
//                     printf("#--- max of fdivs:"); realApp_printd( max, 10); printf("\n");
//                     printf("#--- mu          :"); realApp_printd( mu,  10); printf("\n");
//                     printf("#--- max < mu?  :%d\n", realApp_lt(max, mu));
//                 }
//             
//                 if ( realApp_lt(max, mu) ) {
//                     res.nbOfSol = 0;
//                     g = q2;
//                 }
//                 realApp_clear(div);
//                 realApp_clear(max);
//                 realApp_clear(mu);
//             }
            
            g = g+1;
        }
        
        if (metadatas_getVerbo(meta)>=level) {
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

cauchyTest_res cauchyTest_deterministic_exclusion_testNEW( const compRat_t center,
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
    
    int level = 3;
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    cacheCauchy_set_bounds( cacheCau, radius, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cauchy deterministic exclusion test: BEGIN\n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
//     slong q1 = cacheCauchy_nbEvalExref(cacheCau);
    slong q2 = cacheCauchy_nbEvalCeref(cacheCau);
//     slong quo = cacheCauchy_quotientref(cacheCau);
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
//     realApp_ptr wP2 = cacheCauchy_wanErrCeref(cacheCau);
    
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
    
//     int alreadyEvaluated = 0;

    if (res.nbOfSol==0) {
        
        compApp_t sstar, sh;
        realApp_t abssstar, errabssstar, temp, max, e, mu;
        compApp_init(sstar);
        compApp_init(sh);
        realApp_init(abssstar);
        realApp_init(errabssstar);
        realApp_init(temp);
        realApp_init(max);
        realApp_init(e);
        realApp_init(mu);
        
//         realApp_one(e);
//         realApp_div_ui(e, e, 4*q2, res.appPrec);
//         realApp_div_ui(e, e, q2, res.appPrec);
//         realApp_div_ui(e, e, q2, res.appPrec);
//     
//         compApp_set_sisi(sstar, 1, 1);
//         realApp_one(abssstar);
//         realApp_add_error( compApp_realref(sstar), abssstar );
//         realApp_add_error( compApp_imagref(sstar), abssstar );
//         compApp_abs( abssstar, sstar, res.appPrec );
//         realApp_get_rad_realApp( abssstar, abssstar );
//         realApp_mul_2exp_si( abssstar, abssstar, 1 );
        
        res.nbOfSol = -1;
        
//         while ( compApp_contains_zero( sstar ) && realApp_ge( abssstar, e ) ) {
        while (res.nbOfSol==-1){
            // get evaluation points
            cauchyTest_getEvaluationPoints( center, radius,
                                        radius2, vangle, vindex,
                                        cacheCau, CAUCHYTEST_CERTIFIED, res.appPrec);
            // evaluate at points
            cauchyTest_evaluateAtPoints( cache, cacheCau, CAUCHYTEST_CERTIFIED, res.appPrec, inCounting, meta, depth);
            // compute fdivs
            res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, CAUCHYTEST_CERTIFIED, res.appPrec, inCounting, meta);
            
            if (metadatas_getVerbo(meta)>=level) {
                    printf("#--- computeFdivs_fromVals res:  %d\n", res.nbOfSol);
            }
                
            if (res.nbOfSol == 1) {
        
                // compute the max of fdivs 
                compApp_abs( max, cacheCauchy_fdivsCeref(cacheCau) + 0, res.appPrec);
                for (slong g = 1; g < q2; g++) {
                    compApp_abs( temp, cacheCauchy_fdivsCeref(cacheCau) + g, res.appPrec);
                    realApp_max( max, max, temp, res.appPrec);
                }
                // if the max of fdivs is strictly less than 1/(2q2 +1)
                realApp_one(mu);
                realApp_div_ui(mu, mu, 2*q2 + 1, res.appPrec);
                
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#--- max of fdivs:"); realApp_printd( max, 10); printf("\n");
                    printf("#--- mu          :"); realApp_printd( mu,  10); printf("\n");
                    printf("#--- max < mu?  :%d\n", realApp_lt(max, mu));
                }
            
                if ( realApp_lt(max, mu) )
                    res.nbOfSol = 0;
                else {
                    realApp_one(e);
                    realApp_div_ui(e, e, 4*q2, res.appPrec);
                    realApp_div_ui(e, e, q2, res.appPrec);
                    realApp_div_ui(e, e, q2, res.appPrec);
                    // compute all the power sums
                    cauchyTest_computeShApprox_fromVals( sstar, 0, cacheCau, res.appPrec, inCounting, meta );
                    compApp_pow_si( sstar, sstar, 2, res.appPrec );
                    compApp_abs( abssstar, sstar, res.appPrec );
                    realApp_get_rad_realApp( errabssstar, abssstar );
                    realApp_mul_2exp_si( errabssstar, errabssstar, 1 );
                    slong h = 1;
                    while ( (h<=q2-1) 
//                           && compApp_contains_zero( sstar ) 
                          && realApp_lt( abssstar, e ) ) {
                        cauchyTest_computeShApprox_fromVals( sh, h, cacheCau, res.appPrec, inCounting, meta );
                        compApp_pow_si( sh, sh, 2, res.appPrec );
                        compApp_add( sstar, sstar, sh, res.appPrec );
                        compApp_abs( abssstar, sstar, res.appPrec );
                        realApp_get_rad_realApp( errabssstar, abssstar );
                        realApp_mul_2exp_si( errabssstar, errabssstar, 1 );
                        h = h+1;
                    }
                    
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#--- h               :%ld\n", h);
                        printf("#--- sstar           :"); compApp_printd( sstar, 10); printf("\n");
                        printf("#--- abssstar        :"); realApp_printd( abssstar,  10); printf("\n");
                        printf("#--- errabssstar     :"); realApp_printd( errabssstar,  10); printf("\n");
                        printf("#--- e               :"); realApp_printd( e,  10); printf("\n");
                        printf("#--- abssstar < e?   :%d\n", realApp_lt( abssstar, e ));
                        printf("#--- abssstar > e?   :%d\n", realApp_gt( abssstar, e ));
                        printf("#--- errabssstar < e?:%d\n", realApp_lt( errabssstar, e ));
                    }
                
                    if ( realApp_lt( abssstar, e ) ){
                        // then contains no root
                        res.nbOfSol = 0;
                    } else if ( realApp_gt( abssstar, e ) ) {
                        res.nbOfSol = 1;
                    } else if ( realApp_lt( errabssstar, e ) ) {
                        res.nbOfSol = 1;
                    } else {
                        res.nbOfSol = -1;
                    }
                    
// //                     if ( compApp_contains_zero( sstar ) && realApp_lt( abssstar, e ) ){
//                     if ( realApp_lt( abssstar, e ) ){
//                         // then contains no root
//                         res.nbOfSol = 0;
// //                     } else if (!(compApp_contains_zero( sstar ))) {
//                     } else if (realApp_gt( abssstar, e ) ) {
//                             res.nbOfSol = 1;
//                     } else
//                         res.nbOfSol = -1;
                }
                
                if (res.nbOfSol == -2)
                    res.nbOfSol = 1;
                
                if (res.nbOfSol == -1)
                    res.appPrec = 2*res.appPrec;
            }
        }
        
    
        compApp_clear(sstar);
        compApp_clear(sh);
        realApp_clear(abssstar);
        realApp_clear(errabssstar);
        realApp_clear(temp);
        realApp_clear(max);
        realApp_clear(e);
        realApp_clear(mu);
    
    }
    
    for (int j=0; j<nbPowerSums; j++)
        compApp_clear( ps +j );
    ccluster_free(ps);
    
    if ( metadatas_haveToCount(meta) && (inCounting==0) ) {
        metadatas_add_time_CauExTo(meta, (double) (clock() - start));
        metadatas_add_CauchyExTest(meta, depth, res.appPrec);
    }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cauchy deterministic exclusion test: END, res.nbOfSol: %d \n", res.nbOfSol);
    }
    return res;
}
