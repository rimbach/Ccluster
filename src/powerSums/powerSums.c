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
#include <time.h>

void pwSuDatas_set( pwSuDatas_t p, 
                    void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                    slong degree,
                    slong iRnum, slong iRden,
                    slong nbPws,
                    int verb ) {
    
    p->evalPoly = evalFast;
    pwSuDatas_set_isolaRatio_si( p, iRnum, iRden );
    pwSuDatas_set_nbPwSuComp( p, nbPws );
    slong nbP = powerSums_getNbOfPointsForCounting( pwSuDatas_wantedPrec_ptr(p), 
                                                    degree, 
                                                    pwSuDatas_isolaRatio_ptr(p) )
                + pwSuDatas_nbPwSuComp(p) -1;
    pwSuDatas_set_nbPntsEval(p, nbP);
    if (verb>=2) {
        
        printf("nb of power sums computed: %d\n", (int) pwSuDatas_nbPwSuComp(p) );
        printf("iso ratio used for tests:  "); realRat_print( pwSuDatas_isolaRatio_ptr(p) ); printf("\n");
        printf("nb points for eval:        %d\n", (int) pwSuDatas_nbPntsEval(p) );
        
//         realRat_t errorNum, errorDen;
//         realRat_init(errorNum);
//         realRat_init(errorDen);
//         realRat_inv(errorDen, isoRatio);
//         realRat_pow_si(errorDen, errorDen, nbP);
//         realRat_add_si(errorDen, errorDen, -1);
//         realRat_neg(errorDen, errorDen);
//         for (int h = 0; h<nbP/2; h++){
//             realRat_inv(errorNum, isoRatio);
//             realRat_pow_si(errorNum, errorNum, nbP-h);
//             realRat_mul_si(errorNum, errorNum, cacheApp_getDegree(cache));
//             realRat_div( errorNum, errorNum, errorDen );
//             printf("error for h=%d: ", h); realRat_print( errorNum ); printf("\n");
//         }
//         realRat_clear(errorNum);
//         realRat_clear(errorDen);
    }
}

void metadatas_set_pwSuDatas( metadatas_t meta, 
                              void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                              slong degree,
                              slong iRnum, slong iRden,
                              slong nbPws,
                              int verb ) {
    pwSuDatas_set( metadatas_pwSumref(meta), evalFast, degree, iRnum, iRden, nbPws, verb );
}

/* case where only 0-th power sum is computed */
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
                                 slong prec,
                                 metadatas_t meta){
    
    if (metadatas_pwSumref(meta)->evalPoly == NULL) {
        compApp_poly_ptr app = cacheApp_getApproximation ( cache, prec );
        for (slong i=0; i<nbPoints; i++)
            compApp_poly_evaluate2_rectangular(f_val + i, fder_val + i, app, points + i, prec);
    }
    else {
        for (slong i=0; i<nbPoints; i++)
            metadatas_pwSumref(meta)->evalPoly( f_val+i, fder_val + i, points+i, prec);
    }
}

void powerSums_evaluateAtPoints_fast( compApp_ptr f_val,
                                 compApp_ptr fder_val,
                                 const compApp_ptr points,
                                 slong nbPoints,
                                 cacheApp_t cache,
                                 slong prec ){
    compApp_poly_t appDer;
    compApp_poly_init(appDer);
    
    compApp_poly_ptr app = cacheApp_getApproximation ( cache, prec );
    compApp_poly_derivative( appDer, app, prec);
    
    acb_ptr * tree;
    tree = _acb_poly_tree_alloc(nbPoints);
    _acb_poly_tree_build(tree, points, nbPoints, prec);
        
//     acb_poly_evaluate_vec_fast(f_val, app, points, nbPoints, prec);
    _acb_poly_evaluate_vec_fast_precomp(f_val, app->coeffs, app->length, tree, nbPoints, prec);
//     acb_poly_evaluate_vec_fast(fder_val, appDer, points, nbPoints, prec);
    _acb_poly_evaluate_vec_fast_precomp(fder_val, appDer->coeffs, appDer->length, tree, nbPoints, prec);
        
    _acb_poly_tree_free(tree, nbPoints);
    compApp_poly_clear(appDer);
}

// void powerSums_computefdivs_fromVals( compApp_ptr fdivs,
//                                       compApp_ptr f_val,
//                                       compApp_ptr fder_val,
//                                       slong nbPoints,
//                                       slong prec ){
//     
//     for (slong i = 0; i<nbPoints; i++) {
//         compApp_div(fdivs +i, fder_val + i, f_val + i, prec);
//     }
//     
// }

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

void powerSums_computeS0_prec_without_meta(     compApp_t s0, 
                                   compApp_ptr points,
                                   compApp_ptr pointsShifted,
                                   compApp_ptr f_val,
                                   compApp_ptr fder_val,
                                   const compRat_t center,
                                   const realRat_t radius,
                                   cacheApp_t cache,
                                   slong nbPoints,
                                   slong prec){
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, prec);
    powerSums_evaluateAtPoints( f_val, fder_val, pointsShifted, nbPoints, cache, prec, NULL);
//     powerSums_evaluateAtPoints_fast( f_val, fder_val, pointsShifted, nbPoints, cache, prec);
    powerSums_computeS0_fromVals( s0, points, f_val, fder_val, nbPoints, prec );
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
                                   slong prec,
                                   metadatas_t meta, int depth ){
    
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, prec);
    clock_t start = clock();
    
    powerSums_evaluateAtPoints( f_val, fder_val, pointsShifted, nbPoints, cache, prec, meta);
//     powerSums_evaluateAtPoints_fast( f_val, fder_val, pointsShifted, nbPoints, cache, prec);
    if (metadatas_haveToCount(meta))
            metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start) );
    
    powerSums_computeS0_fromVals( s0, points, f_val, fder_val, nbPoints, prec );
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int powerSums_computePsApprox_fromVals(compApp_ptr ps,
                                        const compRat_t center,
                                        const realRat_t radius,
                                        const realRat_t lowerBound,
                                        const realRat_t upperBound,
                                        compApp_ptr points,
//                                         compApp_ptr pointsShifted,
                                        compApp_ptr fvals,
                                        compApp_ptr fdervals,
                                        compApp_ptr fdivs,
                                        slong nbPoints,
                                        slong nbPowerSums,
//                                         cacheApp_t cache,
                                        slong prec,
                                        metadatas_t meta){
    
    int res=1;
    
    realApp_t lb, ub, modulus;
    realApp_init(lb);
    realApp_init(ub);
    realApp_init(modulus);
    realApp_set_realRat(lb, lowerBound, prec);
    realApp_set_realRat(ub, upperBound, prec);
    
    /* compute fdivs; check if fvals contains zero and 
     *                has modulus less than lower bound*/
    for (slong i = 0; (i<nbPoints) && (res==1) ; i++) {
        compApp_abs(modulus, fvals +i, prec);
        if (compApp_contains_zero( fvals +i )){
            /* compute modulus of fvals +i */
            if (realApp_lt( modulus, lb )){
//                 if (metadatas_getVerbo(meta)>3){
//                     printf("---------------------------------------\n");
//                     printf("---test for disk centered in "); compRat_print(center); printf("\n");
//                     printf("---with radius "); realRat_print(radius); printf("\n");
//                     printf("---fvals + %d contains zero:", (int) i);
//                     compApp_printd(fvals + i, 10); printf("\n");
//                     printf("---lower bound on fvals + %d, assuming f is monic: ", (int) i);
//                     realApp_printd(lb, 10); printf("\n");
//                     printf("---|fvals + %d| is less that lower bound:", (int) i);
//                     printf("%d\n",realApp_lt( modulus, lb ));
//                     printf("---------------------------------------\n");
//                 }
                res=-2;
            }
            else {
//                 if (metadatas_getVerbo(meta)>3) {
//                     printf("---------------------------------------\n");
//                     printf("---test for disk centered in "); compRat_print(center); printf("\n");
//                     printf("---with radius "); realRat_print(radius); printf("\n");
//                     printf("---fvals + %d contains zero:", (int) i);
//                     compApp_printd(fvals + i, 10); printf("\n");
//                     printf("---lower bound on fvals + %d, assuming f is monic: ", (int) i);
//                     realApp_printd(lb, 10); printf("\n");
//                     printf("---|fvals + %d| is less that lower bound:", (int) i);
//                     printf("%d\n",realApp_lt( modulus, lb ));
//                     printf("---------------------------------------\n");
//                 }
                res=-1;
            }
        }
        else if (realApp_lt( modulus, lb )) {
//             if (metadatas_getVerbo(meta)>3) {
//                 printf("---------------------------------------\n");
//                 printf("---test for disk centered in "); compRat_print(center); printf("\n");
//                 printf("---with radius "); realRat_print(radius); printf("\n");
//                 printf("---fvals + %d:", (int) i);
//                 compApp_printd(fvals + i, 10); printf("\n");
//                 printf("---lower bound on fvals + %d, assuming f is monic: ", (int) i);
//                 realApp_printd(lb, 10); printf("\n");
//                 printf("---Disk can not have isolation ratio >="); 
//                 realRat_print( metadatas_getIsoRatio(meta) ); printf("\n");
//                 printf("---------------------------------------\n");
//             }
            res = -2;
        }
        else if (!realApp_ge( modulus, lb )){
//             if (metadatas_getVerbo(meta)>3) {
//                 printf("---------------------------------------\n");
//                 printf("---test for disk centered in "); compRat_print(center); printf("\n");
//                 printf("---with radius "); realRat_print(radius); printf("\n");
//                 printf("---prec: %ld\n ", prec);
//                 printf("---fvals + %d:", (int) i);
//                 compApp_printd(fvals + i, 10); printf("\n");
//                 printf("---lower bound on fvals + %d, assuming f is monic: ", (int) i);
//                 realApp_printd(lb, 10); printf("\n");
//                 printf("---Disk may have isolation ratio <"); 
//                 realRat_print( metadatas_getIsoRatio(meta) ); printf("\n");
//                 printf("---------------------------------------\n");
//             }
            res = -1;
        }
        compApp_div(fdivs +i, fdervals + i, fvals + i, prec);
        /* check if the ratio is less than the upper bound */
        compApp_abs(modulus, fdivs +i, prec);
        if (realApp_gt( modulus, ub )) {
//             if (metadatas_getVerbo(meta)>=2) {
//                 printf("---------------------------------------\n");
//                 printf("---test for disk centered in "); compRat_print(center); printf("\n");
//                 printf("---with radius "); realRat_print(radius); printf("\n");
//                 printf("---fdivs + %d:", (int) i);
//                 compApp_printd(fdivs + i, 10); printf("\n");
//                 printf("---upper bound on fdivs + %d: ", (int) i);
//                 realApp_printd(ub, 50); printf("\n");
//                 printf("---Disk can not have isolation ratio >="); 
//                 realRat_print( metadatas_getIsoRatio(meta) ); printf("\n");
//                 printf("---------------------------------------\n");
//             }
            res = -2;
        }
        
    }
    
    if (res==1){
        
        realApp_t radRe, radIm, wP;
        realRat_t wantedPrec;
        realRat_init(wantedPrec);
        realApp_init(radRe);
        realApp_init(radIm);
        realApp_init(wP);
        realRat_set_si(wantedPrec, 1, 4);
        realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
        
        /* compute powerSums */
        for (slong j = 0; j<nbPowerSums; j++)
            compApp_mul(ps+j, fdivs + 0, points + 0, prec);
        for (slong i = 1; i<nbPoints; i++)
            for (slong j = 0; j<nbPowerSums; j++)
                compApp_addmul(ps+j , fdivs + i, points + ((j+1)*i)%nbPoints, prec);  
        for (slong j = 0; j<nbPowerSums; j++)
            compApp_div_si(ps+j, ps+j, nbPoints, prec);
        /* check if precision is OK */
        
        for (slong j = 0; j<nbPowerSums; j++){
            realApp_get_rad_realApp( radRe, compApp_realref(ps+j) );
            realApp_get_rad_realApp( radIm, compApp_imagref(ps+j) );
            res = res && (realApp_lt( radRe, wP )) && (realApp_lt( radIm, wP ));
            res = ((res==1)? 1:-1);
            
//             if (metadatas_getVerbo(meta)>3){
//                 printf("--- %d-th power sum approximation: ", (int) j);
//                 compApp_printd( ps+j, 10 ); printf("\n");
//                 printf("--- errors: ");
//                 realApp_printd(radRe, 5); printf(", "); realApp_printd(radIm, 5); printf("\n");
//                 printf("--- comparaison: %d\n", realApp_lt( radRe, wP ) && realApp_lt( radIm, wP ) );
//                 printf("--- res: %d\n", res );
//             }
        }
        
        realRat_clear(wantedPrec);
        realApp_clear(wP);
        realApp_clear(radRe);
        realApp_clear(radIm);
        
    }
    
    realApp_clear(lb);
    realApp_clear(ub);
    realApp_clear(modulus);
    
    return res;
}

powerSums_res powerSums_computePsApprox(compApp_ptr ps,
                                        const compRat_t center,
                                        const realRat_t radius,
                                        compApp_ptr points,
                                        compApp_ptr pointsShifted,
                                        compApp_ptr fvals,
                                        compApp_ptr fdervals,
                                        compApp_ptr fdivs,
                                        slong nbPoints,
                                        slong nbPowerSums,
                                        cacheApp_t cache,
                                        slong prec,
                                        metadatas_t meta, int depth ){
    powerSums_res res;
    res.appPrec = prec;
    
    /* compute lower bound: (radius^d*(isoRatio-1)^d/isoRatio^d */
    realRat_t lb;
    realRat_init(lb);
    realRat_add_si(lb, metadatas_getIsoRatio(meta), -1);
    realRat_mul(lb, lb, radius);
    realRat_div(lb, lb, metadatas_getIsoRatio(meta));
    realRat_pow_si(lb, lb, cacheApp_getDegree(cache));
    
    /* compute upper bound: (d*(isoRatio-1)/r*(isoRatio-1) */
    realRat_t ub, temp;
    realRat_init(ub);
    realRat_init(temp);
    realRat_add_si(ub, metadatas_getIsoRatio(meta), +1);
    realRat_mul_si(ub, ub, cacheApp_getDegree(cache));
    realRat_add_si(temp, metadatas_getIsoRatio(meta), -1);
    realRat_mul(temp, temp, radius);
    realRat_div(ub, ub, temp);
    
    /* compute points and evals at prec res.appPrec*/
    powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, res.appPrec);
    clock_t start = clock();
    powerSums_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, res.appPrec, meta);
    if (metadatas_haveToCount(meta))
            metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start) );   

    /* compute approximation of Power sums */
    res.nbOfSol = powerSums_computePsApprox_fromVals(ps, center, radius, lb, ub, points, fvals, fdervals, fdivs, nbPoints, nbPowerSums, res.appPrec, meta);
    
    while ( res.nbOfSol ==-1 ) {
        res.appPrec = 2*res.appPrec;
        
        powerSums_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, res.appPrec);
        clock_t start2 = clock();
        powerSums_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, res.appPrec, meta);
        if (metadatas_haveToCount(meta))
            metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start2) );
        /* compute approximation of Power sums */
        res.nbOfSol = powerSums_computePsApprox_fromVals(ps, center, radius, lb, ub, points, fvals, fdervals, fdivs, nbPoints, nbPowerSums, res.appPrec, meta);
    }
    
    realRat_clear(lb);
    realRat_clear(ub);
    realRat_clear(temp);
    return res;
}
    
                                
powerSums_res powerSums_discardingTest( const compRat_t center,
                                        const realRat_t radius,
                                        cacheApp_t cache,
                                        slong nbPoints,
                                        slong nbPowerSums,
                                        slong prec,
                                        metadatas_t meta, int depth){
    clock_t start = clock();
    
    powerSums_res res;
    res.appPrec = prec;
    
    realApp_t wP;
    realRat_t wantedPrec;
    compApp_ptr points;
    compApp_ptr pointsShifted;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    compApp_ptr fdivs;
    compApp_ptr ps;
    
    realRat_init(wantedPrec);
//     realApp_init(radRe);
//     realApp_init(radIm);
    realApp_init(wP);
    
    realRat_set_si(wantedPrec, 1, 4);
    realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
    
    points =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fvals =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdervals =      (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdivs =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    
    ps =        (compApp_ptr) ccluster_malloc( nbPowerSums*sizeof(compApp) );
    
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
        compApp_init( fdivs +i );
    }
    
    for (int j=0; j<nbPowerSums; j++)
        compApp_init( ps +j );
    
    res = powerSums_computePsApprox(ps, center, radius, 
                                    points, pointsShifted, fvals, fdervals, fdivs, 
                                    nbPoints, nbPowerSums, cache, prec, meta, depth );
    
    /* here res.nbOfSol is either -2 and we return -1 */
    /*                    or 1 and we continue */
    res.nbOfSol = ( res.nbOfSol==-2? -1:1 );
    
    int j=0;
    while ( (j<nbPowerSums) && (res.nbOfSol==1) ) {
//     for (int j=0; j<nbPowerSums; j++){
        
        /* get error to which j-th power sum is known */
        /* assume it is (1/4)*2^{ -((nbPowerSums - 1) - j)} */
        realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
        realApp_mul_2exp_si(wP, wP, j+1-nbPowerSums);
//         if (metadatas_getVerbo(meta)>3) {
//             printf("%d-th power sum, error: ", j); realApp_printd(wP, 20); printf("\n"); }
        /* add error to j-th power sum */
//         if (metadatas_getVerbo(meta)>3) {
//             printf("%d-th power sum approx, prec %d: ", j, (int) res.appPrec); compApp_printd(ps+j, 20); printf("\n"); }
        realApp_add_error( compApp_realref(ps+j), wP );
        realApp_add_error( compApp_imagref(ps+j), wP );
//         if (metadatas_getVerbo(meta)>3) {
//             printf("%d-th power sum approx, prec %d: ", j, (int) res.appPrec); compApp_printd(ps+j, 20); printf("\n"); }
        res.nbOfSol = compApp_contains_zero( ps + j ); /* contains at most one integer */
//         if (metadatas_getVerbo(meta)>3) {
//             printf("%d-th power sum contains zero: %d\n", j, res.nbOfSol);  }
        res.nbOfSol = ( (res.nbOfSol==0) ? -1:1 );
        j++;
    }
    
    /* here if (res.nbOfSol==1) then j first power sums contain 0 */
    /* conclude 0 root in the disk and returns 0 */
    res.nbOfSol = ( (res.nbOfSol==1) ? 0:-1 ); 
    
    for (int i=0; i<nbPoints; i++){
        compApp_clear( points +i );
        compApp_clear( pointsShifted +i );
        compApp_clear( fvals +i );
        compApp_clear( fdervals +i );
        compApp_clear( fdivs +i );
    }
    
    for (int j=0; j<nbPowerSums; j++)
        compApp_clear( ps +j );
    
    ccluster_free(points);
    ccluster_free(pointsShifted);
    ccluster_free(fvals);
    ccluster_free(fdervals);
    ccluster_free(fdivs);
    ccluster_free(ps);
    
    realRat_clear(wantedPrec);
//     realApp_clear(radRe);
//     realApp_clear(radIm);
    realApp_clear(wP);
    
    if (metadatas_haveToCount(meta))
        metadatas_add_time_PSTests(meta, (double) (clock() - start));
    
    return res;
}


powerSums_res powerSums_countingTest( const compRat_t center,
                                        const realRat_t radius,
                                        cacheApp_t cache,
                                        slong nbPoints,
                                        int isIsolated,
                                        slong prec,
                                        metadatas_t meta, int depth){
    clock_t start = clock();
    
    powerSums_res res;
    res.appPrec = prec;
    
    realApp_t wP;
    realRat_t wantedPrec;
    compApp_ptr points;
    compApp_ptr pointsShifted;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    compApp_ptr fdivs;
    compApp_ptr ps;
    
    realRat_init(wantedPrec);
//     realApp_init(radRe);
//     realApp_init(radIm);
    realApp_init(wP);
    
    realRat_set_si(wantedPrec, 1, 4);
    realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
    
    points =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fvals =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdervals =      (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    fdivs =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
    
    ps =        (compApp_ptr) ccluster_malloc( sizeof(compApp) );
    
    for (int i=0; i<nbPoints; i++){
        compApp_init( points +i );
        compApp_init( pointsShifted +i );
        compApp_init( fvals +i );
        compApp_init( fdervals +i );
        compApp_init( fdivs +i );
    }
    
    compApp_init( ps );
    
    res = powerSums_computePsApprox(ps, center, radius, 
                                    points, pointsShifted, fvals, fdervals, fdivs, 
                                    nbPoints, 1, cache, prec, meta, depth );
    
    /* here res.nbOfSol is either -2 and we return -1 */
    /*                    or 1 and we continue */
    res.nbOfSol = ( res.nbOfSol==-2? -1:1 );
    
    if (res.nbOfSol==1) {
        /* get error to which j-th power sum is known */
        /* assume it is (1/4) */
        realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
        /* add error to 0-th power sum */
        realApp_add_error( compApp_realref(ps), wP );
        realApp_add_error( compApp_imagref(ps), wP );
        
        slong nbOfSol = -1;
        int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(ps) );
        int containsZero = realApp_contains_zero( compApp_imagref(ps) );
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
    
    
    for (int i=0; i<nbPoints; i++){
        compApp_clear( points +i );
        compApp_clear( pointsShifted +i );
        compApp_clear( fvals +i );
        compApp_clear( fdervals +i );
        compApp_clear( fdivs +i );
    }
    
    compApp_clear( ps);
    
    ccluster_free(points);
    ccluster_free(pointsShifted);
    ccluster_free(fvals);
    ccluster_free(fdervals);
    ccluster_free(fdivs);
    ccluster_free(ps);
    
    realRat_clear(wantedPrec);
//     realApp_clear(radRe);
//     realApp_clear(radIm);
    realApp_clear(wP);
    
    if (metadatas_haveToCount(meta))
        metadatas_add_time_PSTests(meta, (double) (clock() - start));
    
    return res;
}
                                        
    
// powerSums_res powerSums_countingTest( const compRat_t center,
//                                       const realRat_t radius,
//                                       cacheApp_t cache,
//                                       slong nbPoints,
//                                       int isIsolated,
//                                       slong prec,
//                                       metadatas_t meta, int depth){
//     
//     clock_t start = clock();
//     
//     powerSums_res res;
//     res.appPrec = prec;
//     
//     compApp_t s0;
//     realApp_t radRe, radIm, wP;
//     realRat_t wantedPrec;
//     compApp_ptr points;
//     compApp_ptr pointsShifted;
//     compApp_ptr fvals;
//     compApp_ptr fdervals;
// //     compApp_ptr fdivs;
//     
//     realRat_init(wantedPrec);
//     compApp_init(s0);
//     realApp_init(radRe);
//     realApp_init(radIm);
//     realApp_init(wP);
//     
//     realRat_set_si(wantedPrec, 1, 4);
//     realApp_set_realRat( wP, wantedPrec, CCLUSTER_DEFAULT_PREC);
//     
// //     slong degree = cacheApp_getDegree (cache);
// //     slong nbPoints = powerSums_getNbOfPointsForCounting( wantedPrec, degree, isoRatio );
// //     printf(" nb evaluation points: %d \n", (int) nbPoints);
//     
//     points =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     fvals =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     fdervals =      (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
// //     fdivs =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     
//     for (int i=0; i<nbPoints; i++){
//         compApp_init( points +i );
//         compApp_init( pointsShifted +i );
//         compApp_init( fvals +i );
//         compApp_init( fdervals +i );
// //         compApp_init( fdivs +i );
//     }
//     
//     powerSums_computeS0_prec( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec, meta, depth );
//     realApp_get_rad_realApp( radRe, compApp_realref(s0) );
//     realApp_get_rad_realApp( radIm, compApp_imagref(s0) );
//     
//     res.nbOfSol = -1;
//     
//     while ( (res.nbOfSol==-1)&&( (!compApp_is_finite(s0)) || (!realApp_lt( radRe, wP )) || (!realApp_lt( radIm, wP )) )) {
//         res.appPrec = 2*res.appPrec;
//         powerSums_computeS0_prec( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec, meta, depth );
//         realApp_get_rad_realApp( radRe, compApp_realref(s0) );
//         realApp_get_rad_realApp( radIm, compApp_imagref(s0) );
// //            printf(" s0 at prec %d: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
//         
//         if ( (isIsolated==0) && (!compApp_is_finite(s0)) ){
//         /* check if one fval contains zero and has prec >= 53*/  
//             for (int i=0; i<nbPoints; i++) {
//                 if ( ( compApp_contains_zero( fvals +i ) ) && ( compApp_checkAccuracy( fvals +i, 53) )){
//                     res.nbOfSol = -2;
//                     i = nbPoints;
//                     printf("--- ps counting test returns -2; prec: %d \n", (int) res.appPrec);
//                     printf("------ fvals + %d: ", i); compApp_printd(fvals +i, 20); printf("\n");
//                 }
//             }
//         }
//     }
//     
// //     printf(" s0 at prec %d: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
//     
//     if (res.nbOfSol==-1) {
//         if (isIsolated==0) {
// #ifdef CCLUSTER_STATS_PS_MACIS
//             /* adding too much error -> find the closest integer */
//             realApp_add( wP, wP, wP, CCLUSTER_DEFAULT_PREC);
// #endif
//         }
//         realApp_add_error( compApp_realref(s0), wP );
//         realApp_add_error( compApp_imagref(s0), wP );
// //          printf(" s0 at prec %d with 1/4 error: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
//     
//         slong nbOfSol = -1;
//         int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
//         int containsZero = realApp_contains_zero( compApp_imagref(s0) );
//         if (! ( unique && containsZero) ){
//             res.nbOfSol = -1;
// //         printf("--- ps counting test returns -1; prec: %d \n", (int) res.appPrec);
// //         printf("------ s0 with 1/2 error: "); compApp_printd(s0, 20); printf("\n");
// //         printf("------ unique: %d\n", unique);
// //         printf("------ contains 0: %d\n", containsZero);
//         }
//         else {
//             res.nbOfSol = (int) nbOfSol;
//         }
//     }
//     
//     for (int i=0; i<nbPoints; i++){
//         compApp_clear( points +i );
//         compApp_clear( pointsShifted +i );
//         compApp_clear( fvals +i );
//         compApp_clear( fdervals +i );
// //         compApp_clear( fdivs +i );
//     }
//     
//     ccluster_free(points);
//     ccluster_free(pointsShifted);
//     ccluster_free(fvals);
//     ccluster_free(fdervals);
// //     ccluster_free(fdivs);
//     
//     realRat_clear(wantedPrec);
//     compApp_clear(s0);
//     realApp_clear(radRe);
//     realApp_clear(radIm);
//     realApp_clear(wP);
//     
//     if (metadatas_haveToCount(meta))
//         metadatas_add_time_PSTests(meta, (double) (clock() - start));
//     
//     return res;
// }

powerSums_res powerSums_countingTest_with_isoRatio( const compRat_t center,
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
    
    powerSums_computeS0_prec_without_meta( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec );
    realApp_get_rad_realApp( radRe, compApp_realref(s0) );
    realApp_get_rad_realApp( radIm, compApp_imagref(s0) );
//     printf(" s0 at prec %d: ", (int) res.appPrec); compApp_printd(s0, 20); printf("\n");
    while ( (!compApp_is_finite(s0)) || (!realApp_lt( radRe, wP )) || (!realApp_lt( radIm, wP )) ) {
        res.appPrec = 2*res.appPrec;
        powerSums_computeS0_prec_without_meta( s0, points, pointsShifted, fvals, fdervals, center, radius, cache, nbPoints, res.appPrec );
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
