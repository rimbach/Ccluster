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
    
    /* compute nbEvalEx = ceil ( log_isoRatio (2*degree +1) ) + 1*/
    realApp_t liR;
    realApp_init(liR);
    realApp_set_realRat(liR, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_log(liR, liR, CCLUSTER_DEFAULT_PREC );
    realApp_t q1App;
    realApp_init(q1App);
    realApp_set_si(q1App, 4*degree+1);
    realApp_log(q1App, q1App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q1App, q1App, liR, CCLUSTER_DEFAULT_PREC );
    slong q1 = realApp_ceil_si( q1App, CCLUSTER_DEFAULT_PREC ) + cacheCauchy_nbPwSuExref(cache);
    cacheCauchy_nbEvalExref(cache) = q1;
    
    /* compute nbEvalCe = max ( log_isoRatio (2*degree*nbEvalCo +1), degree +1 ) s.t. nbEvalCe multiple of nbEvalCo */
    realApp_t q2App;
    realApp_init(q2App);
    realApp_set_si(q2App, 4*degree*q1+1);
    realApp_log(q2App, q2App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q2App, q2App, liR, CCLUSTER_DEFAULT_PREC );
    slong q2 = realApp_ceil_si( q2App, CCLUSTER_DEFAULT_PREC ) +1;
    q2 = CCLUSTER_MAX(q2, degree +1);
    slong quo = ((slong) q2/q1) +1;
    q2 = q1*quo;
    cacheCauchy_quotientref(cache) = quo;
    cacheCauchy_nbEvalCeref(cache) = q2;
    
    /* compute error cacheCauchy_wanErrExref(cache)[i] = (d*isoRatio^(-q1+i))/(1-isoRatio^(-q1)) */
    realApp_ptr wP = cacheCauchy_wanErrExref(cache) + 0;
    realApp_t tempApp;
    realApp_init(tempApp);
    realApp_set_realRat(wP, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP, wP, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP, wP, q1, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP, wP, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP, wP, degree, CCLUSTER_DEFAULT_PREC);
    for (int i = 1; i<cacheCauchy_nbPwSuExref(cache); i++) {
        wP = cacheCauchy_wanErrExref(cache) + i;
        realApp_mul_realRat(wP, cacheCauchy_wanErrExref(cache) + (i-1), cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    }
        
    
    /* compute error cacheCauchy_wanErrCeref(cache) = (d*isoRatio^(-q2))/(1-isoRatio^(-q2)) */
    realApp_ptr wP2 = cacheCauchy_wanErrCeref(cache);
    realApp_set_realRat(wP2, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP2, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP2, wP2, q2, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP2, wP2, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP2, wP2, degree, CCLUSTER_DEFAULT_PREC);
    
    /* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
    realRat_ptr lb = cacheCauchy_lBoundUnref(cache);
    realRat_add_si(lb, cacheCauchy_isoRatioref(cache), -1);
    realRat_div(lb, lb, cacheCauchy_isoRatioref(cache));
    realRat_pow_si(lb, lb, degree);
    
    /* compute upper bound unit: (d*(isoRatio+1)/(isoRatio-1) */
    realRat_ptr ub = cacheCauchy_uBoundUnref(cache);
    realRat_t temp;
    realRat_init(temp);
    realRat_add_si(ub, cacheCauchy_isoRatioref(cache), +1);
    realRat_mul_si(ub, ub, degree);
    realRat_add_si(temp, cacheCauchy_isoRatioref(cache), -1);
    realRat_div(ub, ub, temp);
    
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
    
    realApp_clear(liR);
    realApp_clear(q1App);
    realApp_clear(q2App);
    realApp_clear(tempApp);
    
    realRat_clear(temp);
    
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
