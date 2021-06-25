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

cauchyTest_res cauchyTest_probabilistic_verification( const compDsk_t Delta,
                                                      slong nbOfRoots,
                                                      const realRat_t a,
                                                      cacheApp_t cache,
                                                      cacheCauchy_t cacheCau,
                                                      slong prec,
                                                      metadatas_t meta, int depth){
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    /* call probabilistic counting test */
    if (realRat_cmp(a, cacheCauchy_isoRatioref(cacheCau))==0) {
        res = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, res.appPrec, meta, depth);
    }
    else {
        res = cauchyTest_probabilistic_counting_withIsoRatio( a, Delta, cache, cacheCau, res.appPrec, meta, depth);
    }
    
    if (metadatas_getVerbo(meta)>=3) {
            printf("#------------------cauchyTest_probabilistic_verification: res of proba counting is %d\n", res.nbOfSol);
    }

    if ( res.nbOfSol != nbOfRoots )
        res.nbOfSol = -1;
    else
        res.nbOfSol = nbOfRoots;
    
    return res;
}

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

/* this version works for any iso ratio, not necessarily the one
 * defined in CacheCau 
 * the disk is not necessarily isoRatio-isolated
 * can fail, otherwise
 * returns the number of roots in Delta */
cauchyTest_res cauchyTest_probabilistic_counting_withIsoRatio( const realRat_t isoRatio,
                                                               const compDsk_t Delta,
                                                               cacheApp_t cache,
                                                               cacheCauchy_t cacheCau,
                                                               slong prec,
                                                               metadatas_t meta, int depth){
    cauchyTest_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    
    clock_t start = clock();
    double evalTime=0;
    
    /* want |s0*-s0| less than 1/4 */
    /* => number of evaluation points = ceil ( log_isoRatio (4*degree +1) ) + 1*/
    slong q = cacheCauchy_get_NbOfEvalPoints( cacheCauchy_degreeref(cacheCau), isoRatio, 1, CCLUSTER_DEFAULT_PREC );
    if (metadatas_getVerbo(meta)>=3)
        printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: nb of eval points: %ld\n", q);
    
    /* want w(s0*)<1/4 */
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_d(errAp, 0.25);
    
    compApp_t point       ;
    compApp_t pointShifted;
    compApp_t fval        ;
    compApp_t fderval     ;
    compApp_t fdiv        ;
    compApp_t s0star      ;
    
    compApp_init( point        );
    compApp_init( pointShifted );
    compApp_init( fval         );
    compApp_init( fderval      );
    compApp_init( fdiv         );
    compApp_init( s0star       );
    
    realApp_t radRe;
    realApp_t radIm;
    realApp_init(radRe);
    realApp_init(radIm);
    
    realApp_t lb, ub, modulus;
    realApp_init(lb);
    realApp_init(ub);
    realApp_init(modulus);
    
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
    compApp_poly_ptr app = NULL;
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=3)
            printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: res.appPrec: %ld\n", res.appPrec);
        
        cacheCauchy_lBoundApp( lb, cacheCauchy_degreeref(cacheCau), isoRatio, compDsk_radiusref(Delta), res.appPrec);
        cacheCauchy_uBoundApp( ub, cacheCauchy_degreeref(cacheCau), isoRatio, compDsk_radiusref(Delta), res.appPrec);
    
        enoughPrec = 1;
        compApp_zero(s0star);
        compApp_set_compRat(c, compDsk_centerref(Delta), res.appPrec);
        
        if (cacheCauchy_evalFastref(cacheCau) == NULL)
            app = cacheApp_getApproximation ( cache, res.appPrec );
        
        for(slong i=0; i<q; i++) {
            realRat_set_si(argu, 2*i, q);
            compApp_set_realRat(a, argu, res.appPrec);
            acb_exp_pi_i( point, a, res.appPrec);
            compApp_mul_realRat(pointShifted, point, compDsk_radiusref(Delta), res.appPrec);
            compApp_add( pointShifted, c, pointShifted, res.appPrec);
            
            start = clock();
            if (cacheCauchy_evalFastref(cacheCau) == NULL)
                compApp_poly_evaluate2_rectangular(fval, fderval, app, pointShifted, res.appPrec);
            else
                (cacheCauchy_evalFastref(cacheCau))( fval, fderval, pointShifted, res.appPrec);
            evalTime += (double) (clock() - start);
            
            compApp_abs(modulus, fval, res.appPrec);
            
            if (compApp_contains_zero( fval )){
//                 if (metadatas_getVerbo(meta)>=3)
//                     printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: fval contains 0, i: %ld\n", i);
                if (realApp_lt( modulus, lb )){
                    enoughPrec=-2;
                    break;
                }
                else {
                    res.appPrec = 2*res.appPrec;
                    enoughPrec = 0;
                    break;
                }
            } else if (realApp_lt( modulus, lb )) {
                enoughPrec=-2;
                break;
            } else if (!realApp_ge( modulus, lb )){
                res.appPrec = 2*res.appPrec;
                enoughPrec = 0;
                break;
            }
            
            compApp_div(fdiv, fderval, fval, res.appPrec);
            /* check if the ratio is less than the upper bound */
            compApp_abs(modulus, fdiv, res.appPrec);
            if (realApp_gt( modulus, ub )) {
                enoughPrec=-2;
                break;
            }
            
            compApp_mul(fdiv, fdiv, point, res.appPrec);
            compApp_mul_realRat(fdiv, fdiv, compDsk_radiusref(Delta), res.appPrec);
            compApp_div_si( fdiv, fdiv, q, res.appPrec);
            compApp_add(s0star, s0star, fdiv, res.appPrec);
            
            /* check error of s0star */
            realApp_get_rad_realApp( radRe, compApp_realref(s0star) );
            realApp_get_rad_realApp( radIm, compApp_imagref(s0star) );
//             realApp_mul_realRat( radRe, radRe, compDsk_radiusref(Delta), appPrec );
//             realApp_mul_realRat( radIm, radIm, compDsk_radiusref(Delta), appPrec );
            if ( realApp_ge( radRe, errAp ) || realApp_ge( radIm, errAp ) ){
//                 if (metadatas_getVerbo(meta)>=3)
//                     printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: s0star not precise enough, i: %ld\n", i);
                res.appPrec = 2*res.appPrec;
                enoughPrec = 0;
                break;
            }
                
        }
//         if (metadatas_getVerbo(meta)>=3)
//             if (enoughPrec == 1)
//                 printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: s0star precise enough\n");
        
    }
    
    if (enoughPrec!=-2) {
    
        /* add error */
        realApp_add_error( compApp_realref(s0star), errAp );
        realApp_add_error( compApp_imagref(s0star), errAp );
    
        slong nbOfSol = -1;
        int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0star) );
        int containsZero = realApp_contains_zero( compApp_imagref(s0star) );
        if (! ( unique && containsZero) ){
            res.nbOfSol = -1;
        }
        else {
            res.nbOfSol = (int) nbOfSol;
        }
    } else {
        res.nbOfSol = -1;
    }
        
//     realApp_clear(liR);
//     realApp_clear(qApp);
    realApp_clear(lb);
    realApp_clear(ub);
    realApp_clear(modulus);
    realApp_clear(errAp);
    compApp_clear( point        );
    compApp_clear( pointShifted );
    compApp_clear( fval         );
    compApp_clear( fderval      );
    compApp_clear( fdiv         );
    compApp_clear( s0star       );
    realApp_clear(radRe);
    realApp_clear(radIm);
    compApp_clear(c);
    compApp_clear(a);
    realRat_clear(argu);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoED(meta, evalTime);
        metadatas_add_CauchyCoEvalsD(meta, depth, q);
    }
    
    return res;
}

