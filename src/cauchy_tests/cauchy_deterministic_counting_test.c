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
cauchyTest_res cauchyTest_deterministic_counting_test( const compRat_t center,
                                                       const realRat_t radius,           
                                                       cacheApp_t cache,
                                                       cacheCauchy_t cacheCau,
                                                       slong prec,
                                                       metadatas_t meta, int depth){
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    realRat_t RadInf, RadSup, Rad, ExRad, ratio;
    realApp_t nbCentersApp;
    realRat_init(RadInf);
    realRat_init(RadSup);
    realRat_init(Rad);
    realRat_init(ExRad);
    realRat_init(ratio);
    realApp_init(nbCentersApp);
    
    realRat_set(RadInf, radius);
    realRat_mul_si(RadSup, RadInf, 4);
    realRat_add(Rad, RadInf, RadSup);
    realRat_div_ui(Rad, Rad, 2);
    realRat_sub(ExRad, RadSup, RadInf);
    realRat_div_ui(ExRad, ExRad, 2);
    realRat_div(ExRad, ExRad, cacheCauchy_isoRatioref(cacheCau));
    realRat_div(ratio, Rad, ExRad);
    realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
    realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
    slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    
    cacheCauchy_set_bounds( cacheCau, Rad, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#---cauchy deterministic counting test: \n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
        printf("#------ nb of discs on contour for certification: %d \n", (int) nbCenters);
        printf("#------ radius of discs on contour for certification: "); realRat_print( ExRad ); printf("\n");
    }
    
    compApp_t s0;
    compApp_init(s0);
    int alreadyEvaluated = 0;
    
    /* apply counting test to D(center, 2*radius) */
    res = cauchyTest_computeS0Approx(s0, center, Rad,
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
    
    if ( res.nbOfSol > 0 ) { 
        /* try certification */
        cauchyTest_res resEx;
        resEx.nbOfSol = 0;
        resEx.appPrec = res.appPrec;
        for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
            resEx = cauchyTest_deterministic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cache, cacheCau, resEx.appPrec, CAUCHYTEST_INCOUNTIN, meta, depth );
        }
        
        if (resEx.nbOfSol != 0) {
            res.nbOfSol = -1;
            if (metadatas_getVerbo(meta)>3) {
                printf("#------ certification failed!\n");
            }
        }
    }
    
    
    compApp_clear(s0);
    realRat_clear(RadInf);
    realRat_clear(RadSup);
    realRat_clear(Rad);
    realRat_clear(ExRad);
    realRat_clear(ratio);
    realApp_clear(nbCentersApp);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
        metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
    }
    
    return res;
    
}

/* version for newton:
 * D(center, radSup) contains nbOfRoots roots */
cauchyTest_res cauchyTest_deterministic_counting_test_for_newton( const compRat_t center,
                                                                  realRat_t radInf,  
                                                                  const realRat_t radSup,
                                                                  slong nbOfRoots,
                                                                  cacheApp_t cache,
                                                                  cacheCauchy_t cacheCau,
                                                                  slong prec,
                                                                  metadatas_t meta, int depth){
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    realRat_t Rad, ExRad, ratio, nRad;
    realApp_t nbCentersApp;
    realRat_init(Rad);
    realRat_init(ExRad);
    realRat_init(ratio);
    realRat_init(nRad);
    realApp_init(nbCentersApp);
    
    realRat_add(Rad, radInf, radSup);
    realRat_div_ui(Rad, Rad, 2);
    realRat_sub(ExRad, radSup, radInf);
    realRat_div_ui(ExRad, ExRad, 2);
//     realRat_div(ExRad, ExRad, cacheCauchy_isoRatioref(cacheCau));
    realRat_div(ratio, Rad, ExRad);
    realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
    realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
    slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    
    realRat_set(nRad, radInf);
    realRat_mul_si(nRad, nRad, 2);
    realRat_div(nRad, ExRad, nRad);
    realRat_add_si(nRad, nRad, 1);
    realRat_div(nRad, Rad, nRad);
    
    cacheCauchy_set_bounds( cacheCau, Rad, CCLUSTER_DEFAULT_PREC);
//     cacheCauchy_set_bounds( cacheCau, radInf, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#---cauchy deterministic counting test in Newton iteration: \n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
        printf("#------ ratio radSup/radInf: ");
        realApp_t rat;
        realApp_init(rat);
        realApp_set_realRat(rat, radSup, CCLUSTER_DEFAULT_PREC);
        realApp_div_realRat(rat, rat, radInf, CCLUSTER_DEFAULT_PREC);
        realApp_printd(rat, 10); printf("\n");
        realApp_clear(rat);
        printf("#------ ratio mu/rho': "); realRat_print(ratio); printf("\n");
        printf("#------ ratio nRad/radInf: ");
        realApp_init(rat);
        realApp_set_realRat(rat, nRad, CCLUSTER_DEFAULT_PREC);
        realApp_div_realRat(rat, rat, radInf, CCLUSTER_DEFAULT_PREC);
        realApp_div_ui(rat, rat, 3, CCLUSTER_DEFAULT_PREC);
        realApp_mul_si(rat, rat, 4, CCLUSTER_DEFAULT_PREC);
        realApp_printd(rat, 10); printf("\n");
        realApp_clear(rat);
        printf("#------ nb of discs on contour for certification: %d \n", (int) nbCenters);
        printf("#------ radius of discs on contour for certification: "); realRat_print( ExRad ); printf("\n");
        printf("#------ number of roots in disc: %ld\n", nbOfRoots);
    }
    
    compApp_t s0;
    compApp_init(s0);
    int alreadyEvaluated = 0;
    
    /* apply counting test to D(center, radInf) */
    res = cauchyTest_computeS0Approx(s0, center, Rad,
                                     NULL, 0, 0, 
                                    0, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_UNCERTIFI, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
//     res = cauchyTest_computeS0Approx(s0, center, radInf,
//                                      NULL, 0, 0, 
//                                     0, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_UNCERTIFI, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
    
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
    if (res.nbOfSol != nbOfRoots) {
        res.nbOfSol = -1;
        if (metadatas_getVerbo(meta)>3) {
                printf("#------ not all solutions!\n");
            }
    } else if ( res.nbOfSol > 0 ) { 
        /* try certification */
        cauchyTest_res resEx;
        resEx.nbOfSol = 0;
        resEx.appPrec = res.appPrec;
        for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
            resEx = cauchyTest_deterministic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cache, cacheCau, resEx.appPrec, CAUCHYTEST_INCOUNTIN, meta, depth );
        }
        
        if (resEx.nbOfSol != 0) {
            res.nbOfSol = -1;
            if (metadatas_getVerbo(meta)>3) {
                printf("#------ certification failed!\n");
            }
        }
    }
    
    if (res.nbOfSol>0)
        realRat_set(radInf, nRad);
    
    compApp_clear(s0);
    realRat_clear(Rad);
    realRat_clear(ExRad);
    realRat_clear(ratio);
    realRat_clear(nRad);
    realApp_clear(nbCentersApp);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
        metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
    }
    
    return res;
}

cauchyTest_res cauchyTest_deterministic_counting_test_combinatorial( const compRat_t center,
                                                                     realRat_t radius,
                                                                     slong nbOfRoots,
                                                                     cacheApp_t cache,
                                                                     cacheCauchy_t cacheCau,
                                                                     slong prec,
                                                                     metadatas_t meta, int depth){
    
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    realRat_t gamma, rad;
    realRat_init(gamma);
    realRat_init(rad);
    realRat_mul(gamma, cacheCauchy_isoRatioref(cacheCau), cacheCauchy_isoRatioref(cacheCau));
    fmpz_add_si( realRat_numref(gamma), realRat_numref(gamma), 1);
    realRat_set(rad, radius);
    
    if (metadatas_getVerbo(meta)>3) {
        printf("#---cauchy deterministic counting test in Newton iteration, combinatorial version: \n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
        printf("#------ max number of roots in disc: %ld\n", nbOfRoots);
        printf("#------ gamma: "); realRat_print(gamma); printf("\n");
    }
    
    slong i=0;
    res.nbOfSol = nbOfRoots;
    compApp_t s0;
    compApp_init(s0);
    
    while ((i < nbOfRoots+1) && (res.nbOfSol == nbOfRoots) ) {
        
        cacheCauchy_set_bounds( cacheCau, rad, CCLUSTER_DEFAULT_PREC);
        int alreadyEvaluated = 0;
        /* apply counting test to D(center, rad) */
        res = cauchyTest_computeS0Approx(s0, center, rad,
                                         NULL, 0, 0, 
                                         0, &alreadyEvaluated, cache, cacheCau, CAUCHYTEST_UNCERTIFI, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
    
        res.nbOfSol = ( res.nbOfSol==-2? -1:0 );
        if (res.nbOfSol==0) {
            realApp_add_error( compApp_realref(s0), wP + 0 );
            realApp_add_error( compApp_imagref(s0), wP + 0 );
            slong nbOfSol = -1;
            int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
            int containsZero = realApp_contains_zero( compApp_imagref(s0) );
            if (! ( unique && containsZero) ){
                res.nbOfSol = -1;
            }
            else {
                res.nbOfSol = (int) nbOfSol;
            }
        }
        
        if (metadatas_getVerbo(meta)>3) {
            if ((res.nbOfSol == -1)||(res.nbOfSol != nbOfRoots)){
                printf("# ------%ld-th iteration, certification failed; res.nbOfSol: %d\n", i, res.nbOfSol);
            }
        }
        
        i++;
        realRat_div( rad, rad, gamma );
    }
    
    realRat_clear(gamma);
    realRat_clear(rad);
    compApp_clear(s0);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
        metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
    }
    
    return res;
    
}

