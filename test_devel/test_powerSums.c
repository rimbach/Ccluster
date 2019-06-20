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

#include <stdio.h>
#include "numbers/realApp.h"
#include "polynomials/realRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "caches/cacheApp.h"
#include "powerSums/powerSums.h"

int main() {
    
    /* test ceiling */
    realApp_t a;
    realApp_init(a);
    realApp_set_d(a, 3.2);
    slong ca = realApp_ceil_si(a, 53);
    printf(" ceil: %d \n", (int) ca);
    realApp_clear(a);
    
    /* test nbOfPoints for eval */
    realRat_t wP, iR;
    realRat_init(wP);
    realRat_init(iR);
    realRat_set_si(wP, 1, 4);
    realRat_set_si(iR, 5, 4);
    
    slong nb = powerSums_getNbOfPointsForCounting( wP, 128, iR );
    printf(" nb: %d \n", (int) nb);
    
    /*test checkAccuracy */
    realApp_t bre, bim;
    compApp_t b;
    realApp_init(bre);
    realApp_init(bim);
    compApp_init(b);
    realRat_set_si(wP, 1, 10);
    
    realApp_set_realRat(bre, wP, 53);
    realApp_set_realRat(bim, wP, 53);
    compApp_set_real_realApp(b, bre);
    compApp_set_imag_realApp(b, bim);
    
    printf( " check accuracy 32: %d\n", compApp_checkAccuracy( b, 32));
    printf( " check accuracy 51: %d\n", compApp_checkAccuracy( b, 51));
    printf( " check accuracy 52: %d\n", compApp_checkAccuracy( b, 52));
    printf( " check accuracy 53: %d\n", compApp_checkAccuracy( b, 53));
    
    realApp_clear(bre);
    realApp_clear(bim);
    compApp_clear(b);
    
    realRat_clear(wP);
    realRat_clear(iR);
    
    /*test powerSums_getEvaluationPoints*/
    
    slong nbPoints = 3;
    compApp_ptr points;
    compApp_ptr pointsShifted;
    
    points = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
    }
    
    compRat_t center;
    realRat_t radius;
    compRat_init(center);
    realRat_init(radius);
    compRat_set_sisi(center, 1,1,1,1);
    realRat_set_si(radius, 1,1);
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, 53);
    
    for (int i=0; i<nbPoints; i++){
        printf("point %d:         ", i); compApp_printd(points + i, 10); printf("\n");
        printf("point %d shifted: ", i); compApp_printd(pointsShifted + i, 10); printf("\n");
        printf("\n");
    }
    
    for (int i=0; i<nbPoints; i++){
        compApp_clear( points +i );
        compApp_clear( pointsShifted +i );
    }
    ccluster_free(points);
    ccluster_free(pointsShifted);
    
    /*test powerSums_evaluateAtPoints*/
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    bernoulli_polynomial(pbern , 128);
    cacheApp_t cache;
    cacheApp_init_realRat_poly(cache, pbern);
    
    slong prec = 212;
    nbPoints = 3;
    compRat_set_sisi(center, 1,1,1,1);
    realRat_set_si(radius, 1,32);
    
    compApp_ptr fvals;
    compApp_ptr fdervals;
    
    points = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fvals = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdervals = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
    }
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, prec);
    powerSums_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, prec);
    
    for (int i=0; i<nbPoints; i++){
        printf("point %d shifted: ", i); compApp_printd(pointsShifted + i, 10); printf("\n");
        printf("fval            : "   ); compApp_printd(fvals + i, 10); printf("\n");
        printf("fderval         : "   ); compApp_printd(fdervals + i, 10); printf("\n");
        printf("\n");
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
    
    /*test powerSums_computeS0_fromVals*/
    
    slong degree = 128;
    prec = 106;
    realRat_t isoRatio, wantedPrec;
    compApp_t s0;
    compApp_init(s0);
    realRat_init(isoRatio);
    realRat_init(wantedPrec);
    realRat_set_si(isoRatio, 5, 4);
    realRat_set_si(wantedPrec, 1, 4);
    
    nbPoints = powerSums_getNbOfPointsForCounting( wantedPrec, degree, isoRatio );
    printf(" nbPoints: %d \n", (int) nbPoints);
    
    bernoulli_polynomial(pbern , 128);
    cacheApp_init_realRat_poly(cache, pbern);
    
    compRat_set_sisi(center, 0,1,0,1);
    realRat_set_si(radius, 18,1);
    
    points = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fvals = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdervals = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
    }
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, prec);
    powerSums_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, prec);
    powerSums_computeS0_fromVals( s0, points, fvals, fdervals, nbPoints, prec );
    
    printf(" s0 at prec %d: ", (int) prec); compApp_printd(s0, 20); printf("\n");
    
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
    
    
    compApp_clear(s0);
    realRat_clear(wantedPrec);
    
    printf("\n\n");
    
    /*test powerSums_counting test*/
    
    degree = 128;
    prec = 53;
    bernoulli_polynomial(pbern , 128);
    cacheApp_clear(cache);
    cacheApp_init_realRat_poly(cache, pbern);
    
    realRat_set_si(isoRatio, 2, 1);
    compRat_set_sisi(center, 0,1,1,1);
    realRat_set_si(radius, 1,1);
    powerSums_res res = powerSums_countingTest(center, radius, cache, isoRatio, prec );
    printf("--- test for disk centered in "); compRat_print( center);
    printf(" with radius "); realRat_print( radius );
    printf(" isolation ratio "); realRat_print( isoRatio ); printf(":\n");
    printf("------nbOfSols: %d, precision to decide: %d \n", res.nbOfSol, (int) res.appPrec );
    
    realRat_set_si(isoRatio, 5, 4);
    compRat_set_sisi(center, 0,1,0,1);
    realRat_set_si(radius, 18,1);
    res = powerSums_countingTest(center, radius, cache, isoRatio, prec );
    printf("--- test for disk centered in "); compRat_print( center);
    printf(" with radius "); realRat_print( radius );
    printf(" isolation ratio "); realRat_print( isoRatio ); printf(":\n");
    printf("------nbOfSols: %d, precision to decide: %d \n", res.nbOfSol, (int) res.appPrec );
    
    realRat_set_si(isoRatio, 2, 1);
    compRat_set_str(center, "-2437278117108061801009625","302231454903657293676544","303638123570424451456025","151115727451828646838272",10);
    realRat_set_str(radius, "525","604462909807314587353088",10);
    res = powerSums_countingTest(center, radius, cache, isoRatio, prec );
    printf("--- test for disk centered in "); compRat_print( center);
    printf(" with radius "); realRat_print( radius );
    printf(" isolation ratio "); realRat_print( isoRatio ); printf(":\n");
    printf("------nbOfSols: %d, precision to decide: %d \n", res.nbOfSol, (int) res.appPrec );
    
    
    compRat_clear(center);
    realRat_clear(radius);
    realRat_clear(isoRatio);
    realRat_poly_clear(pbern);
    cacheApp_clear(cache);
    
    return 0;
    
}
