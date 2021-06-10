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

/* version for any disk */
cauchyTest_res cauchyTest_deterministic_counting( const compDsk_t Delta,
                                                  const realRat_t a,
                                                  cacheApp_t cache,
                                                  cacheCauchy_t cacheCau,
                                                  slong prec,
                                                  metadatas_t meta, int depth){
    
    realRat_t radInf, radSup;
    realRat_init(radInf);
    realRat_init(radSup);
    
    realRat_set(radInf, compDsk_radiusref(Delta));
    realRat_div(radInf, radInf, a);
    realRat_set(radSup, compDsk_radiusref(Delta));
    realRat_mul(radSup, radSup, a);
    
    cauchyTest_res res = cauchyTest_rootFreeAnnulus( compDsk_centerref(Delta), radInf, radSup,
                                                     cache, cacheCau, prec, meta, depth);
    
    if (res.nbOfSol==0) {
        /* Delta is a-isolated */
        res.nbOfSol = cauchyTest_computeS0compDsk( a, Delta, cache, cacheCau, meta, depth);
    }
    
    realRat_clear(radInf);
    realRat_clear(radSup);
    
    return res;
                
}

cauchyTest_res cauchyTest_deterministic_verification( const compDsk_t Delta,
                                                      slong nbOfRoots,
                                                      cacheApp_t cache,
                                                      cacheCauchy_t cacheCau,
                                                      slong prec,
                                                      metadatas_t meta, int depth){
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    realRat_t radInf, radSup, a;
    realRat_init(radInf);
    realRat_init(radSup);
    realRat_init(a);
    
    realRat_set(a, cacheCauchy_isoRatioref(cacheCau));
    
    realRat_set(radInf, compDsk_radiusref(Delta));
    realRat_div(radInf, radInf, a);
    realRat_set(radSup, compDsk_radiusref(Delta));
    realRat_mul(radSup, radSup, a);
    
    /* first call deterministic test */
    res = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, res.appPrec, meta, depth);
    
    if ( res.nbOfSol != nbOfRoots )
        res.nbOfSol = -1;
    else {
        res = cauchyTest_rootFreeAnnulus( compDsk_centerref(Delta), radInf, radSup,
                                          cache, cacheCau, prec, meta, depth);
        if (res.nbOfSol == 0)
           res.nbOfSol = nbOfRoots;
    }
    
    realRat_clear(radInf);
    realRat_clear(radSup);
    realRat_clear(a);
    
    return res;
}

/* certification that A(center, radInf, radSup) contains no root: */
/* if  0 then A(center, radInf, radSup) contains no root */
/* if -1 then A(center, (radSup+radInf)/2 - (5/4)*isoRatio*(radSup-radInf)/2,  */
/*                      (radSup+radInf)/2 + (5/4)*isoRatio*(radSup-radInf)/2 ) */
/*             contains a root */
cauchyTest_res cauchyTest_rootFreeAnnulus( const compRat_t center,
                                           const realRat_t radInf,  
                                           const realRat_t radSup,
                                           cacheApp_t cache,
                                           cacheCauchy_t cacheCau,
                                           slong prec,
                                           metadatas_t meta, int depth){
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    realRat_t Rad, ExRad, ratio;
    realApp_t nbCentersApp;
    realRat_init(Rad);
    realRat_init(ExRad);
    realRat_init(ratio);
    realApp_init(nbCentersApp);
    
    /* Rad = (radSup + radInf)/2 */
    realRat_add(Rad, radInf, radSup);
    realRat_div_ui(Rad, Rad, 2);
    /* ExRad = (radSup - radInf)/2 */
    realRat_sub(ExRad, radSup, radInf);
    realRat_div_ui(ExRad, ExRad, 2);
    /* ratio = Rad/ExRad */
    realRat_div(ratio, Rad, ExRad);
    /* nbCenters = ceil ( 2*pi*ratio ) */
    realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
    realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
    slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("#---cauchy RootFreeAnnulus: \n");
        printf("#------ nb of discs on contour for certification: %d \n", (int) nbCenters);
        printf("#------ radius of discs on contour for certification: 5/4*"); realRat_print( ExRad ); printf("\n");
        printf("#------ isoRatio for exclusion: "); realRat_print( cacheCauchy_isoRatioref(cacheCau) ); printf("\n");
    }
    
    /* try to prove that A(c, Rad-ExRad, Rad+ExRad) contains no root */
    /*ExRad = (5/4)*ExRad */
    realRat_mul_si( ExRad, ExRad, 5);
    realRat_div_ui( ExRad, ExRad, 4);
    
    cauchyTest_res resEx;
    resEx.nbOfSol = 0;
    resEx.appPrec = res.appPrec;
    for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
        resEx = cauchyTest_deterministic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cache, cacheCau, resEx.appPrec, CAUCHYTEST_INCOUNTIN, meta, depth );
    }
    
    if (resEx.nbOfSol != 0) {
        res.nbOfSol = -1;
        if (metadatas_getVerbo(meta)>=3) {
            printf("#------ certification failed!\n");
        }
    } else {
    /* D(c,Rad) is isoRatio-isolated */
        res.nbOfSol = 0;
    }
    
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

cauchyTest_res cauchyTest_deterministic_counting_combinatorial( const compRat_t center,
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

/* OLD version: */
// cauchyTest_res cauchyTest_deterministic_counting( const compDsk_t Delta,
//                                                   cacheApp_t cache,
//                                                   cacheCauchy_t cacheCau,
//                                                   slong prec,
//                                                   metadatas_t meta, int depth){
//     
//     realRat_t radInf, radSup;
//     realRat_init(radInf);
//     realRat_init(radSup);
//     
//     realRat_set(radInf, compDsk_radiusref(Delta));
//     realRat_mul_si(radInf, radInf, 3);
//     realRat_div_ui(radInf, radInf, 4);
//     realRat_set(radSup, compDsk_radiusref(Delta));
//     realRat_mul_si(radSup, radSup, 2);
//     
//     cauchyTest_res res = cauchyTest_deterministic_counting_with_radInf_radSup( compDsk_centerref(Delta), radInf, radSup, 
//                                                                                cache, cacheCau, prec,
//                                                                                meta, depth);
//     
//     realRat_clear(radInf);
//     realRat_clear(radSup);
//     
//     return res;
//                 
// }
// 
// cauchyTest_res cauchyTest_deterministic_counting_with_radInf_radSup( const compRat_t center,
//                                                                      const realRat_t radInf,  
//                                                                      const realRat_t radSup,
//                                                                      cacheApp_t cache,
//                                                                      cacheCauchy_t cacheCau,
//                                                                      slong prec,
//                                                                      metadatas_t meta, int depth){
//     clock_t start = clock();
//     
//     cauchyTest_res res;
//     res.appPrec = prec;
//     
//     realRat_t Rad, ExRad, ratio, isoRatio;
//     realApp_t nbCentersApp;
//     realRat_init(Rad);
//     realRat_init(ExRad);
//     realRat_init(ratio);
//     realApp_init(nbCentersApp);
//     realRat_init(isoRatio);
//     
//     realRat_add(Rad, radInf, radSup);
//     realRat_div_ui(Rad, Rad, 2);
//     realRat_sub(ExRad, radSup, radInf);
//     realRat_div_ui(ExRad, ExRad, 2);
// //     realRat_div(ExRad, ExRad, cacheCauchy_isoRatioref(cacheCau));
//     realRat_div(ratio, Rad, ExRad);
//     realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
//     slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
//     /* */
//     nbCenters = CCLUSTER_MAX( 12, nbCenters );
//     
//     /* isoRatio = (4/5)*(1+ExRad/Rad) = (4/5)*(1+1/ratio) */
//     realRat_inv(isoRatio, ratio);
//     realRat_add_si( isoRatio, isoRatio, 1);
//     realRat_mul_si( isoRatio, isoRatio, 4);
//     realRat_div_ui( isoRatio, isoRatio, 5);
//     
// //     cacheCauchy_set_bounds( cacheCau, Rad, CCLUSTER_DEFAULT_PREC);
// //     cacheCauchy_set_bounds( cacheCau, radInf, CCLUSTER_DEFAULT_PREC);
//     
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("#---cauchy deterministic counting test: \n");
//         printf("#------ isoRatio for counting : "); realRat_print( isoRatio); printf("\n");
//         printf("#------ isoRatio for exclusion: "); realRat_print( cacheCauchy_isoRatioref(cacheCau) ); printf("\n");
//         printf("#------ nb of discs on contour for certification: %d \n", (int) nbCenters);
//         printf("#------ radius of discs on contour for certification: "); realRat_print( ExRad ); printf("\n");
// //         printf("#------ number of roots in disc: %ld\n", nbOfRoots);
//     }
//     
//     /* try to prove that D(c,Rad) is isoRatio-isolated */
//     cauchyTest_res resEx;
//     resEx.nbOfSol = 0;
//     resEx.appPrec = res.appPrec;
//     for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
//         resEx = cauchyTest_deterministic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cache, cacheCau, resEx.appPrec, CAUCHYTEST_INCOUNTIN, meta, depth );
//     }
//     
//     if (resEx.nbOfSol != 0) {
//         res.nbOfSol = -1;
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#------ certification failed!\n");
//         }
//     } else {
//     /* D(c,Rad) is isoRatio-isolated */
//         compDsk_t Delta;
//         compDsk_init(Delta);
//         compDsk_set_compRat_realRat(Delta, center, Rad);
//         slong nbOfRootsC = cauchyTest_computeS0compDsk( isoRatio, Delta, cache, cacheCau, meta, depth);
//         compDsk_clear(Delta);
//         res.nbOfSol = nbOfRootsC;
//     }
//     
//     realRat_clear(isoRatio);
//     realRat_clear(Rad);
//     realRat_clear(ExRad);
//     realRat_clear(ratio);
//     realApp_clear(nbCentersApp);
//     
//     if (metadatas_haveToCount(meta)) {
//         metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
//         metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
//     }
//     
//     return res;
// }

