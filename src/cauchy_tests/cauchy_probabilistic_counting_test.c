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
                                                  cacheCauchy_t cacheCau,
                                                  slong prec,
                                                  metadatas_t meta, int depth){
    
    int level = 4;
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    cacheCauchy_set_bounds( cacheCau, compDsk_radiusref(Delta), CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cauchy probabilistic counting: \n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
    compApp_t s0;
    compApp_init(s0);
    int alreadyEvaluated = 0;
    
    res = cauchyTest_computeS0Approx(s0, compDsk_centerref(Delta), compDsk_radiusref(Delta),
                                     NULL, 0, 0, &alreadyEvaluated, cacheCau, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
    
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
    
    if (metadatas_getVerbo(meta)>=level) {
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
                                                               cacheCauchy_t cacheCau,
                                                               slong prec,
                                                               metadatas_t meta, int depth){
    int level = 4;
    cauchyTest_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    
    clock_t start = clock();
    double evalTime=0;
    
    /* want |s0*-s0| less than 1/4 */
    /* => number of evaluation points = ceil ( log_isoRatio (4*degree +1) ) + 1*/
    slong q = cacheCauchy_get_NbOfEvalPoints( cacheCauchy_degreeref(cacheCau), isoRatio, 1, CCLUSTER_DEFAULT_PREC );
    if (metadatas_getVerbo(meta)>=level)
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
    
    realApp_t lb, ub;
    realApp_init(lb);
    realApp_init(ub);
    
    compApp_t c;
    
    compApp_init(c);
    
//     compApp_poly_ptr app = NULL;
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=level)
            printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: res.appPrec: %ld\n", res.appPrec);
        
        cacheCauchy_lBoundApp( lb, cacheCauchy_degreeref(cacheCau), isoRatio, compDsk_radiusref(Delta), res.appPrec);
        cacheCauchy_uBoundApp( ub, cacheCauchy_degreeref(cacheCau), isoRatio, compDsk_radiusref(Delta), res.appPrec);
    
        enoughPrec = 1;
        compApp_zero(s0star);
        compApp_set_compRat(c, compDsk_centerref(Delta), res.appPrec);
        
        for(slong i=0; (i<q)&&(enoughPrec==1); i++) {
            
            cauchyTest_computePointPointShifted( point, pointShifted, c, q, i, compDsk_radiusref(Delta), res.appPrec);
            
            start = clock();
            cacheCauchy_eval( fval, fderval, pointShifted, 1, cacheCau, res.appPrec, meta);
            evalTime += (double) (clock() - start);
            
            enoughPrec = cauchyTest_compute_fdiv_checkPrecAndBounds( fdiv, fval, fderval, lb, ub, res.appPrec );
            if (enoughPrec==-1) {
                res.appPrec = 2*res.appPrec;
                enoughPrec = 0;
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
//                 if (metadatas_getVerbo(meta)>=level)
//                     printf("#------------------cauchyTest_probabilistic_counting_withIsoRatio: s0star not precise enough, i: %ld\n", i);
                res.appPrec = 2*res.appPrec;
                enoughPrec = 0;
                break;
            }
                
        }
        
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
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoEv(meta, evalTime);
        metadatas_add_CauchyCoEvals(meta, depth, q);
    }
    
    return res;
}

