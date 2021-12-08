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

void cauchyTest_fmatheta( realRat_t res, const realRat_t a, const realRat_t theta ) {
    
    realRat_t t1, t2;
    realRat_init(t1);
    realRat_init(t2);
    
    /* res = (5/4)theta */
    realRat_set_si( res, 5, 4 );
    realRat_mul   ( res, res, theta );
    realRat_one(t1);
    realRat_one(t2);
    /* t1 = (1-(5/4)theta) */
    realRat_sub(t1, t1, res);
    /* t2 = (1+(5/4)theta) */
    realRat_add(t2, t2, res);
    /* res = a^2*t1 + t2 */
    realRat_mul(res, a, a);
    realRat_mul(res, res, t1);
    realRat_add(res, res, t2);
    /* res = (1/(2a))*(a^2*t1 + t2) */
    realRat_div(res, res, a);
    realRat_div_ui(res, res, 2);
    
    realRat_clear(t1);
    realRat_clear(t2);
}

void cauchyTest_fpatheta( realRat_t res, const realRat_t a, const realRat_t theta ){
    
    realRat_t t1, t2;
    realRat_init(t1);
    realRat_init(t2);
    
    /* res = (5/4)theta */
    realRat_set_si( res, 5, 4 );
    realRat_mul   ( res, res, theta );
    realRat_one(t1);
    realRat_one(t2);
    /* t1 = (1+(5/4)theta) */
    realRat_add(t1, t1, res);
    /* t2 = (1-(5/4)theta) */
    realRat_sub(t2, t2, res);
    /* res = a^2*t1 + t2 */
    realRat_mul(res, a, a);
    realRat_mul(res, res, t1);
    realRat_add(res, res, t2);
    /* res = (1/(2a))*(a^2*t1 + t2) */
    realRat_div(res, res, a);
    realRat_div_ui(res, res, 2);
    
    realRat_clear(t1);
    realRat_clear(t2);
}

cauchyTest_res cauchyTest_rootFreeAnnulus_verification( const compDsk_t Delta,
                                                      slong nbOfRoots,
                                                      const realRat_t a,
                                                      cacheCauchy_t cacheCau,
                                                      slong prec,
                                                      metadatas_t meta, int depth){
    
    int level =4;
    cauchyTest_res res;
    res.appPrec = prec;
    
    realRat_t radInf, radSup;
    realRat_init(radInf);
    realRat_init(radSup);
    
    realRat_set(radInf, compDsk_radiusref(Delta));
    realRat_div(radInf, radInf, a);
    realRat_set(radSup, compDsk_radiusref(Delta));
    realRat_mul(radSup, radSup, a);
    
    /* first call probabilistic test */
    if (realRat_cmp(a, cacheCauchy_isoRatioref(cacheCau))==0) {
        res = cauchyTest_probabilistic_counting( Delta, cacheCau, res.appPrec, meta, depth);
    }
    else {
        res = cauchyTest_probabilistic_counting_withIsoRatio( a, Delta, cacheCau, res.appPrec, meta, depth);
    }
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#------------------cauchyTest_deterministic_verification: res of proba counting is %d\n", res.nbOfSol);
    }

    if ( res.nbOfSol != nbOfRoots )
        res.nbOfSol = -1;
    else {
        res = cauchyTest_rootFreeAnnulus( compDsk_centerref(Delta), radInf, radSup,
                                          cacheCau, prec, meta, depth);
        if (res.nbOfSol == 0)
           res.nbOfSol = nbOfRoots;
        else
           res.nbOfSol = -1;
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------------cauchyTest_deterministic_verification: res of verification is %d\n", res.nbOfSol);
        }
    }
    
    realRat_clear(radInf);
    realRat_clear(radSup);
    
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
                                           cacheCauchy_t cacheCau,
                                           slong prec,
                                           metadatas_t meta, int depth){
    
    int level =4;
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
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---------cauchy RootFreeAnnulus: \n");
        printf("#------------nb of discs on contour for certification: %d \n", (int) nbCenters);
        printf("#------------radius of discs on contour for certification: 5/4*"); realRat_print( ExRad ); printf("\n");
        printf("#------------isoRatio for exclusion: "); realRat_print( cacheCauchy_isoRatioref(cacheCau) ); printf("\n");
    }
    
    /* try to prove that A(c, Rad-ExRad, Rad+ExRad) contains no root */
    /*ExRad = (5/4)*ExRad */
    realRat_mul_si( ExRad, ExRad, 5);
    realRat_div_ui( ExRad, ExRad, 4);
    
    cauchyTest_res resEx;
    resEx.nbOfSol = 0;
    resEx.appPrec = res.appPrec;
    for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
        resEx = cauchyTest_probabilistic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cacheCau, resEx.appPrec, meta, depth );
    }
    
    if (resEx.nbOfSol != 0) {
        res.nbOfSol = -1;
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------certification failed!\n");
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

/* version for any disk */
// cauchyTest_res cauchyTest_rootFreeAnnulus_counting( const compDsk_t Delta,
//                                                   const realRat_t a,
//                                                   cacheCauchy_t cacheCau,
//                                                   slong prec,
//                                                   metadatas_t meta, int depth){
//     
//     realRat_t radInf, radSup;
//     realRat_init(radInf);
//     realRat_init(radSup);
//     
//     realRat_set(radInf, compDsk_radiusref(Delta));
//     realRat_div(radInf, radInf, a);
//     realRat_set(radSup, compDsk_radiusref(Delta));
//     realRat_mul(radSup, radSup, a);
//     
//     cauchyTest_res res = cauchyTest_rootFreeAnnulus( compDsk_centerref(Delta), radInf, radSup,
//                                                      cacheCau, prec, meta, depth);
//     
//     if (res.nbOfSol==0) {
//         /* Delta is a-isolated */
//         res.nbOfSol = cauchyTest_computeS0compDsk( a, Delta, cacheCau, meta, depth);
//     }
//     
//     realRat_clear(radInf);
//     realRat_clear(radSup);
//     
//     return res;
//                 
// }

