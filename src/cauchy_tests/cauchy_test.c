/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "cauchy_tests/cauchy_test.h"

void cauchyTest_getEvaluationPoints( compApp_ptr points, 
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

void cauchyTest_evaluateAtPoints( compApp_ptr f_val,
                                 compApp_ptr fder_val,
                                 const compApp_ptr points,
                                 slong nbPoints,
                                 cacheApp_t cache,
                                 slong prec,
                                 metadatas_t meta){
    
    if (metadatas_pwSumref(meta)->evalPoly == NULL) {
        compApp_poly_ptr app = cacheApp_getApproximation ( cache, prec );
        for (slong i=0; i<nbPoints; i++)
            compApp_poly_evaluate2_rectangular(f_val + i, fder_val + i, app, points + i, prec);
//             compApp_poly_evaluate2_horner(f_val + i, fder_val + i, app, points + i, prec);
    }
    else {
        for (slong i=0; i<nbPoints; i++)
            metadatas_pwSumref(meta)->evalPoly( f_val+i, fder_val + i, points+i, prec);
    }
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeFdivs_fromVals(const realApp_t lb,
                                     const realApp_t ub,
                                     compApp_ptr fvals,
                                     compApp_ptr fdervals,
                                     compApp_ptr fdivs,
                                     slong nbPoints,
                                     slong prec,
                                     metadatas_t meta){
    
    int res=1;
    
    realApp_t modulus;
    realApp_init(modulus);
    
    /* compute fdivs; check if fvals contains zero and 
     *                has modulus less than lower bound*/
    for (slong i = 0; (i<nbPoints) && (res==1) ; i++) {
        compApp_abs(modulus, fvals +i, prec);
        if (compApp_contains_zero( fvals +i )){
            /* compute modulus of fvals +i */
            if (realApp_lt( modulus, lb )){
                res=-2;
            }
            else {
                res=-1;
            }
        }
        else if (realApp_lt( modulus, lb )) {
            res = -2;
        }
        else if (!realApp_ge( modulus, lb )){
            res = -1;
        }
        compApp_div(fdivs +i, fdervals + i, fvals + i, prec);
        /* check if the ratio is less than the upper bound */
        compApp_abs(modulus, fdivs +i, prec);
        if (realApp_gt( modulus, ub )) {
            res = -2;
        }
        
    }

    realApp_clear(modulus);
    
    return res;
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeS0Approx_fromVals(compApp_t ps,
                                        const compRat_t center,
                                        const realRat_t radius,
                                        const realApp_t lb,
                                        const realApp_t ub,
                                        const realApp_t wP,
                                        compApp_ptr points,
//                                         compApp_ptr pointsShifted,
                                        compApp_ptr fvals,
                                        compApp_ptr fdervals,
                                        compApp_ptr fdivs,
                                        slong nbPoints,
                                        slong rotation,
//                                         slong nbPowerSums,
//                                         cacheApp_t cache,
                                        slong prec,
                                        metadatas_t meta){
    
    int res=1;
    
//     printf("cauchyTest_computeS0Approx_fromVals, rotation: %ld\n", rotation);
//     printf("cauchyTest_computeS0Approx_fromVals, (nbPoints-rotation) modulo nbPoints: %ld\n", (nbPoints- rotation)%nbPoints);
//     if (res==1){
        
        realApp_t radRe, radIm;
        realApp_init(radRe);
        realApp_init(radIm);
        
        /* compute s0 */
        compApp_mul(ps, fdivs + (nbPoints - rotation)%nbPoints, points + nbPoints%nbPoints, prec);
        for (slong i = 1; i<nbPoints; i++)
            compApp_addmul(ps, fdivs + (nbPoints + i - rotation)%nbPoints, points + (nbPoints + i)%nbPoints, prec);  
        compApp_div_si(ps, ps, nbPoints, prec);
        /* no need to scale by the radius because points are already */
//         compApp_mul_realRat_in_place(ps, radius, CCLUSTER_DEFAULT_PREC);
        
        /* check if precision is OK */
        realApp_get_rad_realApp( radRe, compApp_realref(ps) );
        realApp_get_rad_realApp( radIm, compApp_imagref(ps) );
        res = res && (realApp_lt( radRe, wP )) && (realApp_lt( radIm, wP ));
        res = ((res==1)? 1:-1);
        
//         if (metadatas_getVerbo(meta)>3){
//             printf("--- %d-th power sum approximation: ", (int) j);
//             compApp_printd( ps+j, 10 ); printf("\n");
//             printf("--- errors: ");
//             realApp_printd(radRe, 5); printf(", "); realApp_printd(radIm, 5); printf("\n");
//             printf("--- comparaison: %d\n", realApp_lt( radRe, wP ) && realApp_lt( radIm, wP ) );
//             printf("--- res: %d\n", res );
//         }
        
        realApp_clear(radRe);
        realApp_clear(radIm);
        
//     }
    
    return res;
}

cauchyTest_res cauchyTest_computeS0Approx(compApp_t ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          compApp_ptr points,
                                          compApp_ptr pointsShifted,
                                          compApp_ptr fvals,
                                          compApp_ptr fdervals,
                                          compApp_ptr fdivs,
                                          slong nbPoints,
                                          slong rotation,
                                          int *alreadyEvaluated,
                                          realApp_t lb,
                                          realApp_t ub,
                                          realApp_t wP,
                                          cacheApp_t cache,
                                          slong prec,
                                          metadatas_t meta, int depth ){
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
    
    if (*alreadyEvaluated == 0) {
        
        while (res.nbOfSol==-1){
            /* compute points and evals at prec res.appPrec*/
            cauchyTest_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, res.appPrec);
            clock_t start = clock();
            cauchyTest_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, res.appPrec, meta);
            if (metadatas_haveToCount(meta))
                metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start) );   
            res.nbOfSol = cauchyTest_computeFdivs_fromVals(lb, ub, fvals, fdervals, fdivs, nbPoints, res.appPrec, meta);
            if ( res.nbOfSol ==-1 )
                res.appPrec = 2*res.appPrec;
        }
        *alreadyEvaluated = 1;
    }
    
    if (res.nbOfSol == -2)
        return res;
    
    /* compute approximation of Power sums */
    res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, center, radius, lb, ub, wP, points, fvals, fdervals, fdivs, nbPoints, rotation, res.appPrec, meta);
    
    while ( res.nbOfSol ==-1 ) {
        res.appPrec = 2*res.appPrec;
        
        cauchyTest_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, res.appPrec);
        clock_t start2 = clock();
        cauchyTest_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, res.appPrec, meta);
        if (metadatas_haveToCount(meta))
            metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start2) );
        /* compute approximation of Power sums */
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(lb, ub, fvals, fdervals, fdivs, nbPoints, res.appPrec, meta);
        res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, center, radius, lb, ub, wP, points, fvals, fdervals, fdivs, nbPoints, rotation, res.appPrec, meta);
    }
    
    return res;
}

cauchyTest_res cauchyTest_exclusionTest( const compRat_t center,
                                          const realRat_t radius,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
                                          slong prec,
                                          metadatas_t meta, int depth){
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("---cauchy exclusion test: \n");
//         printf("------ isoRatio: "); realRat_print(metadatas_getIsoRatio(meta)); printf("\n");
        printf("------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
    }
    
    /* compute q1 = log_theta (2*degree +1) */
    realApp_t liR;
    realApp_init(liR);
    realApp_set_realRat(liR, metadatas_getIsoRatio(meta), CCLUSTER_DEFAULT_PREC);
    realApp_log(liR, liR, CCLUSTER_DEFAULT_PREC );
    
    realApp_t q1App;
    realApp_init(q1App);
    realApp_set_si(q1App, 4*cacheApp_getDegree(cache)+1);
    realApp_log(q1App, q1App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q1App, q1App, liR, CCLUSTER_DEFAULT_PREC );
    slong q1 = realApp_ceil_si( q1App, CCLUSTER_DEFAULT_PREC ) +1;
    
    /* compute q2 = max ( log_theta (2*degree*q1 +1), degree +1 ) s.t. q2 multiple of q1 */
    realApp_t q2App;
    realApp_init(q2App);
    realApp_set_si(q2App, 4*cacheApp_getDegree(cache)*q1+1);
    realApp_log(q2App, q2App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q2App, q2App, liR, CCLUSTER_DEFAULT_PREC );
    slong q2 = realApp_ceil_si( q2App, CCLUSTER_DEFAULT_PREC ) +1;
    q2 = CCLUSTER_MAX(q2, cacheApp_getDegree(cache) +1);
    slong quo = ((slong) q2/q1) +1;
    q2 = q1*quo;
    
    /* compute error wP = (d*isoRatio^(-q1))/(1-isoRatio^(-q1)) */
    realApp_t wP, tempApp;
    realApp_init(wP);
    realApp_init(tempApp);
    realApp_set_realRat(wP, metadatas_getIsoRatio(meta), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP, wP, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP, wP, q1, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP, wP, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP, wP, cacheApp_getDegree(cache), CCLUSTER_DEFAULT_PREC);
    
    /* compute error wP2 = (d*isoRatio^(-q2))/(1-isoRatio^(-q2)) */
    realApp_t wP2;
    realApp_init(wP2);
    realApp_set_realRat(wP2, metadatas_getIsoRatio(meta), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP2, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP2, wP2, q2, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP2, wP2, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP2, wP2, cacheApp_getDegree(cache), CCLUSTER_DEFAULT_PREC);
    
    if (metadatas_getVerbo(meta)>=3) {
//         printf("------ number of eval points q1: %ld\n", q1);
        printf("------ number of eval points q1: %ld\n", cacheCauchy_nbEvalCoref(cacheCau));
//         printf("------ number of eval points q2: %ld\n", q2);
        printf("------ number of eval points q2: %ld\n", cacheCauchy_nbEvalCeref(cacheCau));
//         printf("------ error wP1: "); realApp_printd(wP, 10); printf("\n");
        printf("------ error wP1: "); realApp_printd(cacheCauchy_wanErrCoref(cacheCau), 10); printf("\n");
//         printf("------ error wP2: "); realApp_printd(wP2, 10); printf("\n");
        printf("------ error wP2: "); realApp_printd(cacheCauchy_wanErrCeref(cacheCau), 10); printf("\n");
    }
    
    /* evaluation points */
    compApp_ptr points;
    compApp_ptr pointsShifted;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    compApp_ptr fdivs;
    
    points =        (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    fvals =         (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    fdervals =      (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    fdivs =         (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    
    for (int i=0; i<q1; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
        compApp_init( fdivs +i );
    }
    
    /* compute lower bound: (radius^d*(isoRatio-1)^d/isoRatio^d */
    realRat_t lb;
    realRat_init(lb);
    realRat_add_si(lb, metadatas_getIsoRatio(meta), -1);
    realRat_mul(lb, lb, radius);
    realRat_div(lb, lb, metadatas_getIsoRatio(meta));
    realRat_pow_si(lb, lb, cacheApp_getDegree(cache));
    
    /* compute upper bound: (d*(isoRatio+1)/r*(isoRatio-1) */
    realRat_t ub, temp;
    realRat_init(ub);
    realRat_init(temp);
    realRat_add_si(ub, metadatas_getIsoRatio(meta), +1);
    realRat_mul_si(ub, ub, cacheApp_getDegree(cache));
    realRat_add_si(temp, metadatas_getIsoRatio(meta), -1);
    realRat_mul(temp, temp, radius);
    realRat_div(ub, ub, temp);
    
    realApp_t lbApp, ubApp;
    realApp_init(lbApp);
    realApp_init(ubApp);
    realApp_set_realRat(lbApp, lb, CCLUSTER_DEFAULT_PREC);
    realApp_set_realRat(ubApp, ub, CCLUSTER_DEFAULT_PREC);
    
    int alreadyEvaluated = 0;
    compApp_t s0;
    compApp_init(s0);
    res = cauchyTest_computeS0Approx(s0, center, radius, 
                                    points, pointsShifted, fvals, fdervals, fdivs, 
                                    q1, 0, &alreadyEvaluated, lbApp, ubApp, wP, cache, prec, meta, depth );
    
    res.nbOfSol = ( res.nbOfSol==-2? 1:0 );
    
    if (res.nbOfSol==0) {
        realApp_add_error( compApp_realref(s0), wP );
        realApp_add_error( compApp_imagref(s0), wP );
        res.nbOfSol = ( compApp_contains_zero(s0)==1? 0:1 );
    }
    
    if (metadatas_getVerbo(meta)>=3) {
//         printf("---cauchy exclusion test: \n");
        printf("------ res: %i\n", res.nbOfSol);
//         if (res.nbOfSol==1){
            printf("------ s0: "); compApp_printd(s0,10); printf("\n");
//         }
    }
    
    for (int i=0; i<q1; i++){
        compApp_clear( points +i );
        compApp_clear( pointsShifted +i );
        compApp_clear( fvals +i );
        compApp_clear( fdervals +i );
        compApp_clear( fdivs +i );
    }
    
    ccluster_free(points);
    ccluster_free(pointsShifted);
    ccluster_free(fvals);
    ccluster_free(fdervals);
    ccluster_free(fdivs);
    
    if (res.nbOfSol==0) {
        
        
        points =        (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
        pointsShifted = (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
        fvals =         (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
        fdervals =      (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
        fdivs =         (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    
        for (int i=0; i<q2; i++){
            compApp_init( points +i );
            compApp_init( pointsShifted +i );
            compApp_init( fvals +i );
            compApp_init( fdervals +i );
            compApp_init( fdivs +i );
        } 
        
        alreadyEvaluated = 0;
        slong g = 0;
        while ((g<q1)&&(res.nbOfSol==0)) {
            
            res = cauchyTest_computeS0Approx(s0, center, radius, 
                                    points, pointsShifted, fvals, fdervals, fdivs, 
                                    q2, quo*g, &alreadyEvaluated, lbApp, ubApp, wP2, cache, res.appPrec, meta, depth );
            
            res.nbOfSol = ( res.nbOfSol==-2? 1:0 );
            
            if (res.nbOfSol==0) {
                realApp_add_error( compApp_realref(s0), wP2 );
                realApp_add_error( compApp_imagref(s0), wP2 );
                res.nbOfSol = ( compApp_contains_zero(s0)==1? 0:1 );
            }
            
            g = g+1;
        }
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("------ res: %i\n", res.nbOfSol);
            printf("------ s0: "); compApp_printd(s0,10); printf("\n");
//         }
        }
        
        for (int i=0; i<q2; i++){
            compApp_clear( points +i );
            compApp_clear( pointsShifted +i );
            compApp_clear( fvals +i );
            compApp_clear( fdervals +i );
            compApp_clear( fdivs +i );
        }  
        
        ccluster_free(points);
        ccluster_free(pointsShifted);
        ccluster_free(fvals);
        ccluster_free(fdervals);
        ccluster_free(fdivs);
        
    }
    
    
    
    realApp_clear(liR);
    realApp_clear(q1App);
    realApp_clear(q2App);
    
    realApp_clear(lbApp);
    realApp_clear(ubApp);
    
    realRat_clear(ub);
    realRat_clear(lb);
    realRat_clear(temp);
    
    realApp_clear(wP);
    realApp_clear(wP2);
    realApp_clear(tempApp);
    
    compApp_clear(s0);
    
    if (metadatas_haveToCount(meta))
        metadatas_add_time_PSTests(meta, (double) (clock() - start));
    
    return res;
}
