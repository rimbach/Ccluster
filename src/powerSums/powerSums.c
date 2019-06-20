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

#include "powerSums/powerSums.h"

slong powerSums_getNbOfPointsForCounting( const realRat_t wantedPrec, slong degree, const realRat_t isoRatio ){
    
    slong res;
    
    realApp_t wP, iR, den;
    realRat_t iR_inv;
    realApp_init(wP);
    realApp_init(iR);
    realApp_init(den);
    realRat_init(iR_inv);
    
    realRat_inv(iR_inv, isoRatio);
    realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC );
    realApp_add_si( den, wP, degree, CCLUSTER_DEFAULT_PREC );
    realApp_div( wP, wP, den, CCLUSTER_DEFAULT_PREC );
    realApp_log( wP, wP, CCLUSTER_DEFAULT_PREC );
    
    realApp_set_realRat( iR, iR_inv, CCLUSTER_DEFAULT_PREC );
    realApp_log( iR, iR, CCLUSTER_DEFAULT_PREC );
    
    realApp_div( wP, wP, iR, CCLUSTER_DEFAULT_PREC );
    res = realApp_ceil_si( wP, CCLUSTER_DEFAULT_PREC );
    
    realApp_clear(wP);
    realApp_clear(iR);
    realApp_clear(den);
    realRat_clear(iR_inv);
    
    return res;
}

void powerSums_getEvaluationPoints( compApp_ptr points, 
                                    compApp_ptr pointsShifted,
                                    const compRat_t center,
                                    const realRat_t radius,
                                    slong nbPoints,
                                    slong prec ) {
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
    compApp_set_compRat(c, center, prec);
    for(slong i=0; i<nbPoints; i++) {
        realRat_set_si(argu, 2*i, nbPoints);
        compApp_set_realRat(a, argu, prec);
        acb_exp_pi_i( points + i, a, prec);
        compApp_mul_realRat_in_place(points + i, radius, prec);
        compApp_add( pointsShifted + i, c, points + i, prec);
    }
    
    compApp_clear(c);
    compApp_clear(a);
    realRat_clear(argu);
}

void powerSums_evaluateAtPoints( compApp_ptr f_val,
                                 compApp_ptr fder_val,
                                 const compApp_ptr points,
                                 slong nbPoints,
                                 cacheApp_t cache,
                                 slong prec ){
    
    compApp_poly_ptr app = cacheApp_getApproximation ( cache, prec );
    for (slong i=0; i<nbPoints; i++)
        compApp_poly_evaluate2_horner(f_val + i, fder_val + i, app, points + i, prec);
}

void powerSums_computeS0_fromVals( compApp_t s0, 
                                   compApp_ptr points,
                                   compApp_ptr f_val,
                                   compApp_ptr fder_val,
                                   slong nbPoints,
                                   slong prec ){
    
    compApp_t temp;
    compApp_init(temp);
    compApp_div(s0, fder_val + 0, f_val + 0, prec);
    compApp_mul(s0, s0, points + 0, prec);
    for (slong i = 1; i<nbPoints; i++) {
        compApp_div(temp, fder_val + i, f_val + i, prec);
        compApp_mul(temp, temp, points + i, prec);
        compApp_add(s0, s0, temp, prec);
    }
    compApp_div_si(s0, s0, nbPoints, prec);
                
    compApp_clear(temp);
}

void powerSums_computeS0_prec(     compApp_t s0, 
                                   compApp_ptr points,
                                   compApp_ptr pointsShifted,
                                   compApp_ptr f_val,
                                   compApp_ptr fder_val,
                                   const compRat_t center,
                                   const realRat_t radius,
                                   cacheApp_t cache,
                                   slong nbPoints,
                                   slong prec ){
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, prec);
    powerSums_evaluateAtPoints( f_val, fder_val, pointsShifted, nbPoints, cache, prec);
    powerSums_computeS0_fromVals( s0, points, f_val, fder_val, nbPoints, prec );
}

powerSums_res powerSums_countingTest( const compRat_t center,
                                      const realRat_t radius,
                                      cacheApp_t cache,
                                      const realRat_t isoRatio,
                                      slong prec ){
    powerSums_res res;
    res.appPrec = prec;
    
    compApp_t s0;
    realApp_t radRe, radIm, wP;
    realRat_t wantedPrec;
    compApp_ptr points;
    compApp_ptr pointsShifted;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    
    realRat_init(wantedPrec);
    compApp_init(s0);
    realApp_init(radRe);
    realApp_init(radIm);
    realApp_init(wP);
    
    realRat_set_si(wantedPrec, 1, 4);
    realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
    
    slong degree = cacheApp_getDegree (cache);
    slong nbPoints = powerSums_getNbOfPointsForCounting( wantedPrec, degree, isoRatio );
//     printf(" nb evaluation points: %d \n", (int) nbPoints);
    
    points =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fvals =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdervals =      (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
    }
    
    powerSums_computeS0_prec( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec );
    realApp_get_rad_realApp( radRe, compApp_realref(s0) );
    realApp_get_rad_realApp( radIm, compApp_imagref(s0) );
//     printf(" s0 at prec %d: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
    while ( (!compApp_is_finite(s0)) || (!realApp_lt( radRe, wP )) || (!realApp_lt( radIm, wP )) ) {
        res.appPrec = 2*res.appPrec;
        powerSums_computeS0_prec( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec );
        realApp_get_rad_realApp( radRe, compApp_realref(s0) );
        realApp_get_rad_realApp( radIm, compApp_imagref(s0) );
//         printf(" s0 at prec %d: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
    }
    
    realApp_add_error( compApp_realref(s0), wP );
    realApp_add_error( compApp_imagref(s0), wP );
//     printf(" s0 at prec %d with 1/4 error: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
    
    slong nbOfSol = -1;
    int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
    int containsZero = realApp_contains_zero( compApp_imagref(s0) );
    if (! ( unique && containsZero) ){
        res.nbOfSol = -1;
    }
    else {
        res.nbOfSol = (int) nbOfSol;
    }
    
    for (int i=0; i<nbPoints; i++){
        compApp_clear( points +i );
        compApp_clear( pointsShifted +i );
        compApp_clear( fvals +i );
        compApp_clear( fdervals +i );
    }
    
    ccluster_free(points);
    ccluster_free(pointsShifted);
    ccluster_free(fvals);
    ccluster_free(fdervals);
    
    realRat_clear(wantedPrec);
    compApp_clear(s0);
    realApp_clear(radRe);
    realApp_clear(radIm);
    realApp_clear(wP);
    
    return res;
}
