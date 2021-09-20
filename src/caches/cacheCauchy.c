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

#include <stdlib.h>
#include "cacheCauchy.h"

/* compute q = ceil ( log_isoRatio (4*degree +1) ) + nbPs */
slong cacheCauchy_get_NbOfEvalPoints( slong degree, const realRat_t isoRatio, slong nbPs, slong prec ) {
    
    realApp_t liR;
    realApp_init(liR);
    realApp_t qApp;
    realApp_init(qApp);
    
    realApp_set_realRat(liR, isoRatio, prec);
    realApp_log(liR, liR, prec );
    realApp_set_si(qApp, 4*degree+1);
    realApp_log(qApp, qApp, prec );
    realApp_div(qApp, qApp, liR, prec );
    slong q = realApp_ceil_si( qApp, prec ) + nbPs;
    
    realApp_clear(liR);
    realApp_clear(qApp);
    
    return q;
}

/* compute q2 = log_isoRatio (2*degree*q1 +1) s.t. q2 multiple of q1 */
slong cacheCauchy_get_NbOfEvalPoints_cert( slong degree, slong q1, const realRat_t isoRatio, slong prec ) {
    
    realApp_t liR;
    realApp_init(liR);
    realApp_t qApp;
    realApp_init(qApp);
    
    realApp_set_realRat(liR, isoRatio, prec);
    realApp_log(liR, liR, prec );
    realApp_set_si(qApp, 4*degree*q1+1);
    realApp_log(qApp, qApp, prec );
    realApp_div(qApp, qApp, liR, prec );
    slong q2 = realApp_ceil_si( qApp, prec ) + 1;
    
    realApp_clear(liR);
    realApp_clear(qApp);
    
    return q2;
}

/* compute error wP = (d*isoRatio^(-q))/(1-isoRatio^(-q)) */
void cacheCauchy_wantedErrorOnS0 (realApp_t wP, slong degree, slong q, const realRat_t isoRatio, slong prec ){
    
    realApp_t tempApp;
    realApp_init(tempApp);
    
    realApp_set_realRat(wP, isoRatio, prec);
    realApp_inv(wP, wP, prec);
    realApp_pow_ui(wP, wP, q, prec);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP, prec);
    realApp_div(wP, wP, tempApp, prec);
    realApp_mul_si(wP, wP, degree, prec);
    
    realApp_clear(tempApp);
}

/* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
void cacheCauchy_lowerBoundUnit( realRat_t lowerBoundUnit, slong degree, const realRat_t isoRatio ) {
    realRat_add_si(lowerBoundUnit, isoRatio, -1);
    realRat_div(lowerBoundUnit, lowerBoundUnit, isoRatio);
    realRat_pow_si(lowerBoundUnit, lowerBoundUnit, degree);
}

/* compute upper bound unit: (d*(isoRatio+1)/(isoRatio-1) */
void cacheCauchy_upperBoundUnit( realRat_t upperBoundUnit, slong degree, const realRat_t isoRatio ) {
    realRat_t temp;
    realRat_init(temp);
    
    realRat_add_si(upperBoundUnit, isoRatio, +1);
    realRat_mul_si(upperBoundUnit, upperBoundUnit, degree);
    realRat_add_si(temp, isoRatio, -1);
    realRat_div(upperBoundUnit, upperBoundUnit, temp);
    
    realRat_clear(temp);
}

void cacheCauchy_lBoundApp( realApp_t lbApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec){
    
    realRat_t lb, rd;
    realRat_init(lb);
    realRat_init(rd);
    
    cacheCauchy_lowerBoundUnit( lb, degree, isoRatio );
    realRat_pow_si(rd, radius, degree);
    realRat_mul( lb, lb, rd );
    
    /* approximate */
    realApp_set_realRat(lbApp, lb, prec);
    
    realRat_clear(lb);
    realRat_clear(rd);
}

void cacheCauchy_uBoundApp( realApp_t ubApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec){
    realRat_t ub;
    realRat_init(ub);
    cacheCauchy_upperBoundUnit( ub, degree, isoRatio );
    realRat_div(ub, ub, radius);
    realApp_set_realRat(ubApp, ub, prec);
    realRat_clear(ub);
}

void cacheCauchy_init_sparseEval ( cacheCauchy_t cache, 
                                   cacheApp_t cachePol ) {
    
    compApp_poly_ptr app = cacheApp_getApproximation ( cachePol, CCLUSTER_DEFAULT_PREC );
    cacheCauchy_inNZCref(cache) = (slong *) ccluster_malloc ( (cacheCauchy_degreeref(cache) + 1)*sizeof(slong) );
    
    cacheCauchy_nbNZCref(cache)=0;
    
    for (slong ind=0; ind <= cacheCauchy_degreeref(cache); ind++) {
        if ( ! compApp_is_zero( (app->coeffs) + ind ) ) {
            cacheCauchy_inNZCref(cache)[cacheCauchy_nbNZCref(cache)]=ind;
            cacheCauchy_nbNZCref(cache)+=1;
        }
    }
    
}

void cacheCauchy_sparseEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, cacheApp_t cachePol,
                              const compApp_t point, slong prec) {
    
    compApp_poly_ptr app    = cacheApp_getApproximation ( cachePol, prec );
    
    slong nbNZC = cacheCauchy_nbNZCref(cache);
    slong * inNZC = cacheCauchy_inNZCref(cache);
    
    compApp_set( fval, (app->coeffs) + 0);
    compApp_zero( fderval );
    
    compApp_t x, xp, mon;
    compApp_init(x);
    compApp_init(xp);
    compApp_init(mon);
    
    slong ind=0;
    if (inNZC[ind] == 0)
        ind++;
    compApp_pow_si( x, point, inNZC[ind], prec );
    compApp_pow_si( xp, point, inNZC[ind]-1, prec );
    
    while ( ind < nbNZC ){
        compApp_mul(mon, (app->coeffs) + inNZC[ind], x, prec);
        compApp_add(fval, fval, mon, prec);
        
        compApp_mul(mon, (app->coeffs) + inNZC[ind], xp, prec);
        compApp_mul_si( mon, mon, inNZC[ind], prec);
        compApp_add(fderval, fderval, mon, prec);
        
        ind++;
        
        if (ind < nbNZC) {
            compApp_pow_si(mon, point, inNZC[ind] - inNZC[ind-1], prec);
            compApp_mul( x, x, mon, prec);
            compApp_mul( xp, xp, mon, prec);
        }
    }
    
//     printf("\n");
    
    compApp_clear(x);
    compApp_clear(xp);
    compApp_clear(mon);
    
}

void cacheCauchy_rectangularEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, cacheApp_t cachePol,
                                   const compApp_t point, slong prec) {
    
    compApp_poly_ptr app    = cacheApp_getApproximation ( cachePol, prec );
    compApp_poly_evaluate2_rectangular(fval, fderval, app, point, prec);
    
}

void cacheCauchy_init ( cacheCauchy_t cache, 
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        const realRat_t isoRatio,
                        slong nbPows,
                        const metadatas_t meta
                      ){
    
    int level = 3;
    
    cacheCauchy_evalFastref(cache) = evalFast;
    cacheCauchy_degreeref(cache) = degree;
    
    /* for sparse evaluation */
    cacheCauchy_nbNZCref(cache) = 0;
    cacheCauchy_inNZCref(cache) = NULL;
    cacheCauchy_choiceref(cache) = -1;
    
    realRat_init(cacheCauchy_isoRatioref(cache));
    realRat_set(cacheCauchy_isoRatioref(cache), isoRatio);
    
    realApp_init(cacheCauchy_precfdivref(cache));
    realApp_set_d(cacheCauchy_precfdivref(cache), .25);
    
    cacheCauchy_nbPwSuExref(cache) = nbPows;
    
//     realApp_init(cacheCauchy_wanErrExref(cache));
    cacheCauchy_wanErrExref(cache) = (realApp_ptr) ccluster_malloc( cacheCauchy_nbPwSuExref(cache)*sizeof(realApp) );
    for (int i=0; i< cacheCauchy_nbPwSuExref(cache); i++)
        realApp_init( cacheCauchy_wanErrExref(cache) + i );
    
    
    realApp_init(cacheCauchy_wanErrCeref(cache));
    
    realRat_init(cacheCauchy_lBoundUnref(cache));
    realRat_init(cacheCauchy_uBoundUnref(cache));
    
    realRat_init(cacheCauchy_curRadiuref(cache));
    realApp_init(cacheCauchy_lBoundApref(cache));
    realApp_init(cacheCauchy_uBoundApref(cache));
    
    /* compute q1 = ceil ( log_isoRatio (4*degree +1) ) + nbPs*/
    slong q1 = cacheCauchy_get_NbOfEvalPoints( degree, cacheCauchy_isoRatioref(cache), cacheCauchy_nbPwSuExref(cache), CCLUSTER_DEFAULT_PREC );
    cacheCauchy_nbEvalExref(cache) = q1;
    
    /* compute nbEvalCe = max ( log_isoRatio (2*degree*nbEvalCo +1), degree +1 ) s.t. nbEvalCe multiple of nbEvalCo */
    slong q2 = cacheCauchy_get_NbOfEvalPoints_cert( degree, q1, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC );
//     printf("q2: %ld\n");
    q2 = CCLUSTER_MAX(q2, degree +1);
//     printf("q2: %ld\n");
    slong quo = ((slong) q2/q1) +1;
    q2 = q1*quo;
    cacheCauchy_quotientref(cache) = quo;
    cacheCauchy_nbEvalCeref(cache) = q2;
    
    /* compute error cacheCauchy_wanErrExref(cache)[i] = (d*isoRatio^(-q1+i))/(1-isoRatio^(-q1)) */
    realApp_ptr wP = cacheCauchy_wanErrExref(cache) + 0;
    cacheCauchy_wantedErrorOnS0 (wP, degree, q1, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC );
    
    for (int i = 1; i<cacheCauchy_nbPwSuExref(cache); i++) {
        wP = cacheCauchy_wanErrExref(cache) + i;
        realApp_mul_realRat(wP, cacheCauchy_wanErrExref(cache) + (i-1), cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    }
        
    
    /* compute error cacheCauchy_wanErrCeref(cache) = (d*isoRatio^(-q2))/(1-isoRatio^(-q2)) */
    realApp_ptr wP2 = cacheCauchy_wanErrCeref(cache);
    cacheCauchy_wantedErrorOnS0 (wP2, degree, q2, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC );
    
    /* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
    realRat_ptr lb = cacheCauchy_lBoundUnref(cache);
    cacheCauchy_lowerBoundUnit( lb, degree, cacheCauchy_isoRatioref(cache) );
    
    /* compute upper bound unit: (d*(isoRatio+1)/(isoRatio-1) */
    realRat_ptr ub = cacheCauchy_uBoundUnref(cache);
    cacheCauchy_upperBoundUnit( ub, degree, cacheCauchy_isoRatioref(cache) );
    
    /* evaluation points for uncertified*/
    cacheCauchy_precEvalExref(cache) = 0;
    compApp_ptr pointsEx        = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    compApp_ptr pointsShiftedEx = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    compApp_ptr fvalsEx         = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    compApp_ptr fdervalsEx      = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    compApp_ptr fdivsEx         = (compApp_ptr) ccluster_malloc( q1*sizeof(compApp) );
    
    for (int i=0; i<q1; i++){
        compApp_init( pointsEx +i );
        compApp_init( pointsShiftedEx +i );
        compApp_init( fvalsEx +i );
        compApp_init( fdervalsEx +i );
        compApp_init( fdivsEx +i );
    }
    
    cacheCauchy_pointsExref(cache)        = pointsEx       ; 
    cacheCauchy_pointsShiftedExref(cache) = pointsShiftedEx; 
    cacheCauchy_fvalsExref(cache)         = fvalsEx        ; 
    cacheCauchy_fdervalsExref(cache)      = fdervalsEx     ; 
    cacheCauchy_fdivsExref(cache)         = fdivsEx        ; 
    
    /* evaluation points for certified*/
    compApp_ptr pointsCe =        (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    compApp_ptr pointsShiftedCe = (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    compApp_ptr fvalsCe =         (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    compApp_ptr fdervalsCe =      (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    compApp_ptr fdivsCe =         (compApp_ptr) ccluster_malloc( q2*sizeof(compApp) );
    
    for (int i=0; i<q2; i++){
        compApp_init( pointsCe +i );
        compApp_init( pointsShiftedCe +i );
        compApp_init( fvalsCe +i );
        compApp_init( fdervalsCe +i );
        compApp_init( fdivsCe +i );
    }
    
    cacheCauchy_pointsCeref(cache)        = pointsCe       ; 
    cacheCauchy_pointsShiftedCeref(cache) = pointsShiftedCe; 
    cacheCauchy_fvalsCeref(cache)         = fvalsCe        ; 
    cacheCauchy_fdervalsCeref(cache)      = fdervalsCe     ; 
    cacheCauchy_fdivsCeref(cache)         = fdivsCe        ;
    
    /*initialize shifted poly */
//     compApp_poly_init2( cacheCauchy_shiftedPolyref(cache), q2);
//     compApp_poly_init2( cacheCauchy_shiftedPolyDerref(cache), q2);
//     /* initialize cache for mod pols */
//     compRat_poly_init(cache->_RP);
//     compRat_poly_init(cache->_RPp);
//     compRat_poly_init(cache->_QP);
//     compRat_poly_init(cache->_QPp);
//     compRat_poly_zero(cache->_RP);
//     compRat_poly_zero(cache->_RPp);
//     compRat_poly_zero(cache->_QP);
//     compRat_poly_zero(cache->_QPp);
//     
//     cache->_sizeCache            = 0;
//     cache->_allocsizeCache       = CACHE_DEFAULT_SIZE;
//     cache->_cacheRP              = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cache->_cacheRPp             = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cache->_cacheQP              = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cache->_cacheQPp             = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsizeCache) * sizeof(compApp_poly_t) );
    
    /* initialize compression */
//        cache->_nbEvalComp = 0;
       
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cacheCauchy: \n");
        printf("#------ assumed isolation ratio                : "); realRat_print(cacheCauchy_isoRatioref(cache)); printf("\n");
        printf("#------ number of power sums  for exclusion    : %ld\n", cacheCauchy_nbPwSuExref(cache));
        printf("#------ number of eval points for exclusion    : %ld\n", cacheCauchy_nbEvalExref(cache));
        printf("#------ number of eval points for certification: %ld\n", cacheCauchy_nbEvalCeref(cache));
        printf("#------ error on si* for exclusion             : "); realApp_printd(cacheCauchy_precfdivref(cache), 10); printf("\n");
        for (int i=0; i<cacheCauchy_nbPwSuExref(cache); i++) {
            printf("#------ error on s%d for exclusion              : ", i); realApp_printd(cacheCauchy_wanErrExref(cache) +i, 10); printf("\n");
        }
        printf("#------ error on s0 for certification          : "); realApp_printd(cacheCauchy_wanErrCeref(cache), 10); printf("\n");
    }
}

void cacheCauchy_get_lBoundApp( realApp_t lbApp, cacheCauchy_t cache, const realRat_t radius, slong prec){
    /* set lb to radius^d */
    realRat_t lb;
    realRat_init(lb);
    realRat_set(lb, radius );
    realRat_pow_si(lb, lb, cacheCauchy_degreeref(cache));
    /* mult by lBoundUn =  (isoRatio-1)^d/isoRatio^d*/
    realRat_mul(lb, lb, cacheCauchy_lBoundUnref(cache));
    /* approximate */
    realApp_set_realRat(lbApp, lb, prec);
    
    realRat_clear(lb);
}

void cacheCauchy_get_uBoundApp( realApp_t ubApp, cacheCauchy_t cache, const realRat_t radius, slong prec){
    realRat_t ub;
    realRat_init(ub);
    realRat_div(ub, cacheCauchy_uBoundUnref(cache), radius);
    realApp_set_realRat(ubApp, ub, prec);
    realRat_clear(ub);
}

void cacheCauchy_set_bounds( cacheCauchy_t cache, const realRat_t radius, slong prec ) {
    
    cacheCauchy_precEvalExref(cache) = 0;
    realRat_set( cacheCauchy_curRadiuref(cache), radius );
        
    /* set temp to radius^d */
    realRat_t temp;
    realRat_init(temp);
    realRat_set(temp, radius );
    realRat_pow_si(temp, temp, cacheCauchy_degreeref(cache));
    /* mult by lBoundUn =  (isoRatio-1)^d/isoRatio^d*/
    realRat_mul(temp, temp, cacheCauchy_lBoundUnref(cache));
    /* approximate */
    realApp_set_realRat(cacheCauchy_lBoundApref(cache), temp, prec);
    
    realRat_div(temp, cacheCauchy_uBoundUnref(cache), radius);
    realApp_set_realRat(cacheCauchy_uBoundApref(cache), temp, prec);
    
    cacheCauchy_precBounref(cache) = prec;
    
    realRat_clear(temp);
    
}

// void cacheCauchy_init_comp ( cacheCauchy_t cache, const realRat_t isoRatio, const compDsk_t Delta, slong nbOfRoots, const realRat_t eps ){
//     
//     realRat_t error;
//     realRat_init  (error);
//     realRat_div_ui(error, eps, 2);
//     /* set cache->_errorComp */
//     realApp_init       ( cacheCauchy_errorCompref(cache)                              );
//     realApp_set_realRat( cacheCauchy_errorCompref(cache), error, CCLUSTER_DEFAULT_PREC);
//     /* set cache->_nbEvalComp */
//     realApp_t qApp, logIsoRatio;
//     realApp_init(qApp);
//     realApp_init(logIsoRatio);
//     realApp_set_realRat(logIsoRatio, isoRatio,                                               CCLUSTER_DEFAULT_PREC);
//     realApp_log        (logIsoRatio, logIsoRatio,                                            CCLUSTER_DEFAULT_PREC);
//     realApp_set_realRat(qApp,        compDsk_radiusref(Delta),                               CCLUSTER_DEFAULT_PREC);
//     realApp_mul_realRat(qApp,        qApp,                     isoRatio,                     CCLUSTER_DEFAULT_PREC);
//     realApp_mul_si     (qApp,        qApp,                     cacheCauchy_degreeref(cache), CCLUSTER_DEFAULT_PREC);
//     /* qApp = (rd*isoRatio)/error */
//     realApp_div_realRat(qApp,        qApp,                     error,                        CCLUSTER_DEFAULT_PREC);
//     /* qApp = ( 1 + (rd*isoRatio)/error ) */
//     realApp_add_si     (qApp,        qApp,                     1,                            CCLUSTER_DEFAULT_PREC);
//     /* qApp = log_{e}       ( 1 + (rd*isoRatio)/error ) */
//     realApp_log        (qApp,        qApp,                                                   CCLUSTER_DEFAULT_PREC);
//     /* qApp = log_{isoRatio}( 1 + (rd*isoRatio)/error ) */
//     realApp_div        (qApp,        qApp,                     logIsoRatio,                  CCLUSTER_DEFAULT_PREC); 
//     /* qApp = ceil( log_{isoRatio}(1 + (rd*isoRatio)/error ) ) + 1 */
//     slong q = realApp_ceil_si( qApp, CCLUSTER_DEFAULT_PREC ) + 1;
//     cacheCauchy_nbEvalCompref(cache) = q;
// //     printf(" Number of eval Points: %ld\n", cacheCauchy_nbEvalCompref(cache) );
//     
//     compApp_ptr pointsComp        = (compApp_ptr) ccluster_malloc( q*sizeof(compApp) );
//     compApp_ptr pointsShiftedComp = (compApp_ptr) ccluster_malloc( q*sizeof(compApp) );
//     compApp_ptr fvalsComp         = (compApp_ptr) ccluster_malloc( q*sizeof(compApp) );
//     compApp_ptr fdervalsComp      = (compApp_ptr) ccluster_malloc( q*sizeof(compApp) );
//     compApp_ptr fdivsComp         = (compApp_ptr) ccluster_malloc( q*sizeof(compApp) );
//     
//     for (slong i=0; i<q; i++){
//         compApp_init( pointsComp +i );
//         compApp_init( pointsShiftedComp +i );
//         compApp_init( fvalsComp +i );
//         compApp_init( fdervalsComp +i );
//         compApp_init( fdivsComp +i );
//     }
//     
//     cacheCauchy_pointsCompref(cache)        = pointsComp       ; 
//     cacheCauchy_pointsShiftedCompref(cache) = pointsShiftedComp; 
//     cacheCauchy_fvalsCompref(cache)         = fvalsComp        ; 
//     cacheCauchy_fdervalsCompref(cache)      = fdervalsComp     ; 
//     cacheCauchy_fdivsCompref(cache)         = fdivsComp        ;
//     
//     realRat_clear(error);
//     realApp_clear(qApp);
//     realApp_clear(logIsoRatio);
// }
// 
// void cacheCauchy_clear_comp ( cacheCauchy_t cache ){
//     realApp_clear( cacheCauchy_errorCompref(cache) );
//     
//     compApp_ptr pointsComp        = cacheCauchy_pointsCompref(cache)       ;
//     compApp_ptr pointsShiftedComp = cacheCauchy_pointsShiftedCompref(cache);
//     compApp_ptr fvalsComp         = cacheCauchy_fvalsCompref(cache)        ;
//     compApp_ptr fdervalsComp      = cacheCauchy_fdervalsCompref(cache)        ;
//     compApp_ptr fdivsComp         = cacheCauchy_fdivsCompref(cache)           ;
//     slong q = cacheCauchy_nbEvalCompref(cache);
//     
//     for (slong i=0; i<q; i++){
//         compApp_clear( pointsComp +i );
//         compApp_clear( pointsShiftedComp +i );
//         compApp_clear( fvalsComp +i );
//         compApp_clear( fdervalsComp +i );
//         compApp_clear( fdivsComp +i );
//     }
//     
//     ccluster_free(pointsComp        );
//     ccluster_free(pointsShiftedComp );
//     ccluster_free(fvalsComp         );
//     ccluster_free(fdervalsComp      );
//     ccluster_free(fdivsComp         );
//     
//     cacheCauchy_nbEvalCompref(cache) = 0;
//     
// }

// //requires: prec is 2^n*CCLUSTER_DEFAULT_PREC
// compApp_poly_ptr cacheCauchy_getApproximation_RP      ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec ) {
//     
//     if (cache->_from_poly==0)
//         return cacheApp_getApproximation ( cache, prec );
//     
//     //get index in cache
//     slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
//     int index = 0;
//     while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
// //     printf("index: %d\n", index); 
//     
//     if (index < cacheCau->_sizeCache)
//         return (cacheCau->_cacheRP)[index];
//     
//     if (cacheCau->_sizeCache == 0) {
//         /* let q=_nbEvalEx */
//         /*compute p modulo x^q -1 and p' modulo x^q -1*/
//         realRat_poly_t div;
//         /* set div to x^q -1 */
//         realRat_poly_init(div);
//         realRat_poly_one(div);
// //         realRat_poly_shift_left(div, div, cacheCau->_nbEvalEx);
//         realRat_poly_shift_left(div, div, (slong) (cacheCau->_degree)/2 );
// //         realRat_poly_set_coeff_si_ui (div, 0, -1, 1);
//         
// //         printf("x^q - 1: ");
// //         realRat_poly_print_pretty( div, "x" );
// //         printf("\n");
//         
// //         compRat_poly_divrem(cacheCau->_QP, cacheCau->_Pmod, cache->_poly, div);
//         realRat_poly_divrem(compRat_poly_realref(cacheCau->_QP), compRat_poly_realref(cacheCau->_RP), compRat_poly_realref(cache->_poly), div);
//         
// //         printf("\n P = (x^q - 1)Q + R: degree of Q: %ld, degree of R: %ld\n", compRat_poly_degree(cacheCau->_QP), compRat_poly_degree(cacheCau->_RP));
// //         realRat_poly_print_pretty( compRat_poly_realref(cacheCau->_QP), "x" );
// //         printf("\n");
// //         realRat_poly_print_pretty( compRat_poly_realref(cacheCau->_RP), "x" );
// //         realRat_poly_t temp;
// //         realRat_poly_init(temp);
// //         realRat_poly_mul(temp, compRat_poly_realref(cacheCau->_QP), div);
// //         realRat_poly_add(temp, temp, compRat_poly_realref(cacheCau->_RP) );
// //         printf("\n");
// //         realRat_poly_print_pretty( temp, "x" );
// //         printf("\n");
// //         realRat_poly_print_pretty( compRat_poly_realref(cache->_poly), "x" );
// //         realRat_poly_clear(temp);
// //         printf("\n");
//         
//         compRat_poly_t Pp;
//         compRat_poly_init(Pp);
//         compRat_poly_derivative(Pp, cache->_poly);
//         realRat_poly_divrem(compRat_poly_realref(cacheCau->_QPp), compRat_poly_realref(cacheCau->_RPp), compRat_poly_realref(Pp), div);
//         
// //         printf("\n p' / (x^q - 1): ");
// //         realRat_poly_print_pretty( compRat_poly_realref(cacheCau->_QPp), "x" );
// //         printf("\n");
// //         realRat_poly_print_pretty( compRat_poly_realref(cacheCau->_RPp), "x" );
// // //         realRat_poly_t temp;
// //         realRat_poly_init(temp);
// //         realRat_poly_mul(temp, compRat_poly_realref(cacheCau->_QPp), div);
// //         realRat_poly_add(temp, temp, compRat_poly_realref(cacheCau->_RPp) );
// //         printf("\n");
// //         realRat_poly_print_pretty( temp, "x" );
// //         printf("\n");
// //         realRat_poly_print_pretty( compRat_poly_realref(Pp), "x" );
// //         realRat_poly_clear(temp);
// //         printf("\n");
//         
//         compRat_poly_clear(Pp);
//         realRat_poly_clear(div);
//     }
//     
//     if (index < cacheCau->_allocsizeCache) {
//         while (index >= cacheCau->_sizeCache){
//             compApp_poly_init(cacheCau->_cacheRP[cacheCau->_sizeCache]);
//             compApp_poly_init(cacheCau->_cacheRPp[cacheCau->_sizeCache]);
//             compApp_poly_init(cacheCau->_cacheQP[cacheCau->_sizeCache]);
//             compApp_poly_init(cacheCau->_cacheQPp[cacheCau->_sizeCache]);
//             slong nprec = (0x1<<(cacheCau->_sizeCache))*CCLUSTER_DEFAULT_PREC;
//             compApp_poly_set_compRat_poly(cacheCau->_cacheRP[cacheCau->_sizeCache], cacheCau->_RP, nprec);
//             compApp_poly_set_compRat_poly(cacheCau->_cacheRPp[cacheCau->_sizeCache], cacheCau->_RPp, nprec);
//             compApp_poly_set_compRat_poly(cacheCau->_cacheQP[cacheCau->_sizeCache], cacheCau->_QP, nprec);
//             compApp_poly_set_compRat_poly(cacheCau->_cacheQPp[cacheCau->_sizeCache], cacheCau->_QPp, nprec);
//             
//             cacheCau->_sizeCache +=1;
//         }
//         return (cacheCau->_cacheRP)[index];
//     }
//     
//     while (index >= cacheCau->_allocsizeCache) 
//         cacheCau->_allocsizeCache += CACHE_DEFAULT_SIZE;
//     
//     cacheCau->_cacheRP  = (compApp_poly_t *) ccluster_realloc (cacheCau->_cacheRP, (cacheCau->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cacheCau->_cacheRPp = (compApp_poly_t *) ccluster_realloc (cacheCau->_cacheRPp, (cacheCau->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cacheCau->_cacheQP  = (compApp_poly_t *) ccluster_realloc (cacheCau->_cacheQP, (cacheCau->_allocsizeCache) * sizeof(compApp_poly_t) );
//     cacheCau->_cacheQPp = (compApp_poly_t *) ccluster_realloc (cacheCau->_cacheQPp, (cacheCau->_allocsizeCache) * sizeof(compApp_poly_t) );
//    
//     while (index >= cacheCau->_sizeCache){
//         
//         compApp_poly_init(cacheCau->_cacheRP[cacheCau->_sizeCache]);
//         compApp_poly_init(cacheCau->_cacheRPp[cacheCau->_sizeCache]);
//         compApp_poly_init(cacheCau->_cacheQP[cacheCau->_sizeCache]);
//         compApp_poly_init(cacheCau->_cacheQPp[cacheCau->_sizeCache]);
//         slong nprec = (0x1<<(cacheCau->_sizeCache))*CCLUSTER_DEFAULT_PREC;
//         compApp_poly_set_compRat_poly(cacheCau->_cacheRP[cacheCau->_sizeCache], cacheCau->_RP, nprec);
//         compApp_poly_set_compRat_poly(cacheCau->_cacheRPp[cacheCau->_sizeCache], cacheCau->_RPp, nprec);
//         compApp_poly_set_compRat_poly(cacheCau->_cacheQP[cacheCau->_sizeCache], cacheCau->_QP, nprec);
//         compApp_poly_set_compRat_poly(cacheCau->_cacheQPp[cacheCau->_sizeCache], cacheCau->_QPp, nprec);
//             
//         cacheCau->_sizeCache +=1;
//     }
//     
//     return (cacheCau->_cacheRP)[index];
//     
// }
// 
// compApp_poly_ptr cacheCauchy_getApproximation_RPp     ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec ){
//     
//     if (cache->_from_poly==0)
//         return NULL;
//     
//     //get index in cache
//     slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
//     int index = 0;
//     while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
//     
//     cacheCauchy_getApproximation_RP      ( cacheCau, cache, prec );
//     return (cacheCau->_cacheRPp)[index];
//     
// }
// 
// compApp_poly_ptr cacheCauchy_getApproximation_QP     ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec ){
//     
//     if (cache->_from_poly==0)
//         return NULL;
//     
//     //get index in cache
//     slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
//     int index = 0;
//     while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
//     
//     cacheCauchy_getApproximation_RP      ( cacheCau, cache, prec );
//     return (cacheCau->_cacheQP)[index];
//     
// }
// 
// compApp_poly_ptr cacheCauchy_getApproximation_QPp     ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec ){
//     
//     if (cache->_from_poly==0)
//         return NULL;
//     
//     //get index in cache
//     slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
//     int index = 0;
//     while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
//     
//     cacheCauchy_getApproximation_RP      ( cacheCau, cache, prec );
//     return (cacheCau->_cacheQPp)[index];
//     
// }
// 
void cacheCauchy_clear ( cacheCauchy_t cache ){
    
    if (cacheCauchy_inNZCref(cache))
        ccluster_free( cacheCauchy_inNZCref(cache) );
    realRat_clear(cacheCauchy_isoRatioref(cache));
    realApp_clear(cacheCauchy_precfdivref(cache));
    
//     realApp_clear(cacheCauchy_wanErrExref(cache));
    for (int i=0; i< cacheCauchy_nbPwSuExref(cache); i++)
        realApp_clear( cacheCauchy_wanErrExref(cache) + i );
    ccluster_free( cacheCauchy_wanErrExref(cache) );
    
    realApp_clear(cacheCauchy_wanErrCeref(cache));
    
    realRat_clear(cacheCauchy_lBoundUnref(cache));
    realRat_clear(cacheCauchy_uBoundUnref(cache));
    
    realRat_clear(cacheCauchy_curRadiuref(cache));
    realApp_clear(cacheCauchy_lBoundApref(cache));
    realApp_clear(cacheCauchy_uBoundApref(cache));
    
    /* evaluation points for uncertified*/
    compApp_ptr pointsEx        = cacheCauchy_pointsExref(cache)       ;
    compApp_ptr pointsShiftedEx = cacheCauchy_pointsShiftedExref(cache);
    compApp_ptr fvalsEx         = cacheCauchy_fvalsExref(cache)        ;
    compApp_ptr fdervalsEx      = cacheCauchy_fdervalsExref(cache)        ;
    compApp_ptr fdivsEx         = cacheCauchy_fdivsExref(cache)           ;
    slong q1 = cacheCauchy_nbEvalExref(cache);
//     printf("q1: %ld\n", q1);
    
    for (int i=0; i<q1; i++){
        compApp_clear( pointsEx +i );
        compApp_clear( pointsShiftedEx +i );
        compApp_clear( fvalsEx +i );
        compApp_clear( fdervalsEx +i );
        compApp_clear( fdivsEx +i );
    }
    
    ccluster_free(pointsEx        );
    ccluster_free(pointsShiftedEx );
    ccluster_free(fvalsEx         );
    ccluster_free(fdervalsEx      );
    ccluster_free(fdivsEx         );
    
    /* evaluation points for certified*/
    compApp_ptr pointsCe        = cacheCauchy_pointsCeref(cache)       ;
    compApp_ptr pointsShiftedCe = cacheCauchy_pointsShiftedCeref(cache);
    compApp_ptr fvalsCe         = cacheCauchy_fvalsCeref(cache)        ;
    compApp_ptr fdervalsCe      = cacheCauchy_fdervalsCeref(cache)        ;
    compApp_ptr fdivsCe         = cacheCauchy_fdivsCeref(cache)           ;
    slong q2 = cacheCauchy_nbEvalCeref(cache);
//     printf("q2: %ld\n", q2);
    
    for (int i=0; i<q2; i++){
        compApp_clear( pointsCe +i );
        compApp_clear( pointsShiftedCe +i );
        compApp_clear( fvalsCe +i );
        compApp_clear( fdervalsCe +i );
        compApp_clear( fdivsCe +i );
    }
    
    ccluster_free(pointsCe        );
    ccluster_free(pointsShiftedCe );
    ccluster_free(fvalsCe         );
    ccluster_free(fdervalsCe      );
    ccluster_free(fdivsCe         );
    
    /*initialize shifted poly */
//     compApp_poly_clear( cacheCauchy_shiftedPolyref(cache));
//     compApp_poly_clear( cacheCauchy_shiftedPolyDerref(cache));
    
//     /* free the caches */
//     compRat_poly_clear(cache->_RP);
//     compRat_poly_clear(cache->_RPp);
//     compRat_poly_clear(cache->_QP);
//     compRat_poly_clear(cache->_QPp);
//     for (int i=0; i<cache->_sizeCache; i++) {
//         compApp_poly_clear( (cache->_cacheRP)[i] );
//         compApp_poly_clear( (cache->_cacheRPp)[i] );
//         compApp_poly_clear( (cache->_cacheQP)[i] );
//         compApp_poly_clear( (cache->_cacheQPp)[i] );
//     }
//     ccluster_free(cache->_cacheRP);
//     ccluster_free(cache->_cacheRPp);
//     ccluster_free(cache->_cacheQP);
//     ccluster_free(cache->_cacheQPp);
}
