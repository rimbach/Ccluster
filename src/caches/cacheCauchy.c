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

void cacheCauchy_init ( cacheCauchy_t cache, 
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        const realRat_t isoRatio,
                        slong nbPows,
                        const metadatas_t meta
                      ){
    
    cacheCauchy_evalFastref(cache) = evalFast;
    cacheCauchy_degreeref(cache) = degree;
    
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
    q2 = CCLUSTER_MAX(q2, degree +1);
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
    
    if (metadatas_getVerbo(meta)>=2) {
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

void cacheCauchy_clear ( cacheCauchy_t cache ){
    
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
}
