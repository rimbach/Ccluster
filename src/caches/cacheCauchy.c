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

void cacheCauchy_eval( compApp_ptr fvals, compApp_ptr fdervals, compApp_ptr points, slong nbPoints, 
                       cacheCauchy_t cache, slong prec ){
    
    if (cacheCauchy_evalFastref(cache)) {
        
        for (slong i = 0; i< nbPoints; i++)
            (cacheCauchy_evalFastref(cache))( fvals + i, fdervals + i, points + i, prec);
        
    } else {
        
        if (cacheCauchy_choiceref(cache) == 1) {
            
            for (slong i = 0; i< nbPoints; i++)
                cacheCauchy_sparseEval2 ( fvals + i, fdervals + i, cache, points + i, prec);
        
        } else if ( (cacheCauchy_choiceref(cache) == -1 )&&(nbPoints > 1) ){
            
            clock_t start2 = clock();
            for (slong i=0; i<nbPoints; i++)
                cacheCauchy_rectangularEval ( fvals + i, fdervals + i, cache, points+ i, prec);
            double timeInRectangular = (double) (clock() - start2);
            start2 = clock();
            for (slong i=0; i<nbPoints; i++)
                cacheCauchy_sparseEval2 ( fvals + i, fdervals + i, cache, points + i, prec);
            double timeInSparse = (double) (clock() - start2);
            if (timeInRectangular <= timeInSparse)
                cacheCauchy_choiceref(cache) = 0;
            else
                cacheCauchy_choiceref(cache) = 1;
        
        } else {
            for (slong i=0; i<nbPoints; i++)
                cacheCauchy_rectangularEval ( fvals + i, fdervals + i, cache, points + i, prec);
        
        }
    }
    
}

void cacheCauchy_sparseEval2 ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
                              const compApp_t point, slong prec) {
    
    compApp_poly_ptr app    = cacheApp_getApproximation ( cacheCauchy_cacheAppref(cache), prec );
    
//     compApp_poly_printd(app, 10);
//     printf("\n");
    
    slong nbNZC = cacheCauchy_nbNZCref(cache);
    slong * inNZC = cacheCauchy_inNZCref(cache);
    
    int log2deg = (int) ceil(log2( (inNZC[nbNZC-1]) - 1 ));
//     printf("deg: %ld, log2deg: %d\n", inNZC[nbNZC-1], log2deg);
    compApp_ptr pows = (compApp_ptr) ccluster_malloc (log2deg*sizeof(compApp));
    compApp_t powp, powd;
    compApp_init(powp);
    compApp_init(powd);
    
    compApp_zero( fval );
    compApp_zero( fderval );
    
    compApp_init(pows+0);
    compApp_set(pows+0, point);
    for (int i=1; i<log2deg; i++){
        compApp_init(pows+i);
        compApp_sqr(pows+i, pows+(i-1), prec); 
    }
    
    for (slong i=0; i<nbNZC; i++){
//         printf("i: %ld, inNZC[i]: %ld\n", i, inNZC[i]);
        if ( inNZC[i]==0 ){
            compApp_zero(powd);
            compApp_one (powp);
        } else if ( inNZC[i]==1 ){ 
            compApp_one(powd);
            compApp_set (powp, point);
        } else if ( inNZC[i]==2 ){
            compApp_set (powd, pows + 0);
            compApp_set (powp, pows + 1);
        } else {
            compApp_one (powd);
            slong pow = inNZC[i] - 1;
            int indmax = (int) ceil(log2( pow ));
//             printf("pow: %ld, ind: %d\n", pow, indmax);
            for(int ind = 0; ind<indmax; ind++){
                if (pow%2)
                    compApp_mul( powd, powd, pows + ind, prec );
                pow=pow>>1;
            }
            compApp_mul(powp, powd, point, prec);
        }
        compApp_mul_si(powd, powd, inNZC[i], prec);
        compApp_addmul(fval,    powp, (app->coeffs) + inNZC[i], prec);
        compApp_addmul(fderval, powd, (app->coeffs) + inNZC[i], prec);
    }
    
    for (int i=0; i<log2deg; i++)
        compApp_clear(pows+i);
    
    ccluster_free(pows);
    compApp_clear(powp);
    compApp_clear(powd);   
}

// void cacheCauchy_sparseEval2 ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
//                               const compApp_t point, slong prec) {
//     
//     compApp_poly_ptr app    = cacheApp_getApproximation ( cacheCauchy_cacheAppref(cache), prec );
//     
//     slong nbNZC = cacheCauchy_nbNZCref(cache);
//     slong * inNZC = cacheCauchy_inNZCref(cache);
//     
// //     printf("nbNZC: %ld\n", nbNZC);
// //     printf("inNZC: [ ");
// //     for (slong i=0; i<nbNZC; i++)
// //         printf("%ld, ", inNZC[i]);
// //     printf(" ]\n");
// //     compApp_poly_printd(app, 10);
// //     printf("\n");
//     
//     compApp_t xp, mon;
//     compApp_init(xp);
//     compApp_init(mon);
//      
//     compApp_set( fval, (app->coeffs) + inNZC[nbNZC-1] );
//     compApp_mul_si( fderval, fval, inNZC[nbNZC-1], prec );
//     
// //     printf("fval   : "); compApp_printd(fval,    10); printf("\n");
// //     printf("fderval: "); compApp_printd(fderval, 10); printf("\n");
//     
//     for (slong ind = nbNZC-2; ind >=0; ind--){
// //         printf("ind: %ld, inNZC[%ld]: %ld, pow: %ld,\n", ind, ind, inNZC[ind], inNZC[ind+1]-inNZC[ind]);
//         if ((ind==0)&&(inNZC[0]==0)){
//             compApp_pow_si( xp, point, inNZC[ind+1]-inNZC[ind] -1 , prec );
//             compApp_mul( fderval, fderval, xp, prec);
//             compApp_mul(xp, xp, point, prec);
//             compApp_mul(fval, fval, xp, prec);
//         } else {
//             compApp_pow_si( xp, point, inNZC[ind+1]-inNZC[ind] , prec );
//             compApp_mul( fderval, fderval, xp, prec);
//             compApp_mul( fval, fval, xp, prec);
//         }
//         compApp_mul_si( mon, (app->coeffs) + inNZC[ind], inNZC[ind], prec);
//         compApp_add( fval, fval, (app->coeffs) + inNZC[ind], prec);
//         compApp_add( fderval, fderval, mon, prec);
//     }
//     
//     if (inNZC[0] > 0) {
//         compApp_pow_si( xp, point, inNZC[0] - 1 , prec );
//         compApp_mul( fderval, fderval, xp, prec );
//         compApp_mul(xp, xp, point, prec);
//         compApp_mul(fval, fval, xp, prec);
// //     compApp_mul_si( mon, (app->coeffs) + inNZC[0], inNZC[0], prec);
// //     compApp_add( fval, fval, (app->coeffs) + inNZC[0], prec);
// //     compApp_add( fderval, fderval, mon, prec);
//     }
//     
//     compApp_clear(xp);
//     compApp_clear(mon);
// }

void cacheCauchy_sparseEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
                              const compApp_t point, slong prec) {
    
    compApp_poly_ptr app    = cacheApp_getApproximation ( cacheCauchy_cacheAppref(cache), prec );
    
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
//     compApp_pow_si( x, point, inNZC[ind], prec );
    compApp_pow_si( xp, point, inNZC[ind]-1, prec );
    compApp_mul(x, xp, point, prec);
    
    while ( ind < nbNZC ){
        compApp_mul(mon, (app->coeffs) + inNZC[ind], x, prec);
        compApp_add(fval, fval, mon, prec);
//         compApp_addmul(fval, (app->coeffs) + inNZC[ind], x, prec);
        
        compApp_mul(mon, (app->coeffs) + inNZC[ind], xp, prec);
        compApp_mul_si( mon, mon, inNZC[ind], prec);
        compApp_add(fderval, fderval, mon, prec);
        
        ind++;
        
        if (ind < nbNZC) {
            /* write (inNZC[ind]-1) = q*(inNZC[ind-1]-1) + r */
//             slong q = (inNZC[ind]-1) / (inNZC[ind-1]-1);
//             slong r = (inNZC[ind]-1) % (inNZC[ind-1]-1);
//             compApp_pow_si( mon, point, r, prec);
//             compApp_pow_si( xp, xp, q, prec );
//             compApp_mul(xp, xp, mon, prec);
//             compApp_mul(x, xp, point, prec);
            
            /* old */
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

void cacheCauchy_rectangularEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, const compApp_t point, slong prec) {
    
    compApp_poly_ptr app    = cacheApp_getApproximation ( cacheCauchy_cacheAppref(cache), prec );
    compApp_poly_evaluate2_rectangular(fval, fderval, app, point, prec);
    
}

void cacheCauchy_init ( cacheCauchy_t cache, cacheApp_t cacheA,
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        const realRat_t isoRatio,
                        slong nbPows,
                        const metadatas_t meta
                      ){
    
    int level = 3;
    
    cacheCauchy_cacheAppref(cache) = (cacheApp_ptr) cacheA;
        
    cacheCauchy_evalFastref(cache) = evalFast;
    cacheCauchy_degreeref(cache) = degree;
    
    /* for sparse evaluation */
    cacheCauchy_nbNZCref(cache) = 0;
    cacheCauchy_inNZCref(cache) = NULL;
    cacheCauchy_choiceref(cache) = -1;
    if (evalFast==NULL) {
        compApp_poly_ptr app = cacheApp_getApproximation ( cacheCauchy_cacheAppref(cache), CCLUSTER_DEFAULT_PREC );
        cacheCauchy_inNZCref(cache) = (slong *) ccluster_malloc ( (cacheCauchy_degreeref(cache) + 1)*sizeof(slong) );
    
        cacheCauchy_nbNZCref(cache)=0;
    
        for (slong ind=0; ind <= cacheCauchy_degreeref(cache); ind++) {
            if ( ! compApp_is_zero( (app->coeffs) + ind ) ) {
                cacheCauchy_inNZCref(cache)[cacheCauchy_nbNZCref(cache)]=ind;
                cacheCauchy_nbNZCref(cache)+=1;
            }
        }
    }
    
    realRat_init(cacheCauchy_isoRatioref(cache));
    realRat_set(cacheCauchy_isoRatioref(cache), isoRatio);
    
    realApp_init(cacheCauchy_precfdivref(cache));
    realApp_set_d(cacheCauchy_precfdivref(cache), .25);
    
    cacheCauchy_nbPwSuExref(cache) = nbPows;
    
    cacheCauchy_wanErrExref(cache) = (realApp_ptr) ccluster_malloc( cacheCauchy_nbPwSuExref(cache)*sizeof(realApp) );
    for (int i=0; i< cacheCauchy_nbPwSuExref(cache); i++)
        realApp_init( cacheCauchy_wanErrExref(cache) + i );
    
    realRat_init(cacheCauchy_lBoundUnref(cache));
    realRat_init(cacheCauchy_uBoundUnref(cache));
    
    realRat_init(cacheCauchy_curRadiuref(cache));
    realRat_zero(cacheCauchy_curRadiuref(cache));
    realApp_init(cacheCauchy_lBoundApref(cache));
    realApp_init(cacheCauchy_uBoundApref(cache));
    
    
    /* compute q1 = ceil ( log_isoRatio (4*degree +1) ) + nbPs*/
    slong q1 = cacheCauchy_get_NbOfEvalPoints( degree, cacheCauchy_isoRatioref(cache), cacheCauchy_nbPwSuExref(cache), CCLUSTER_DEFAULT_PREC );
    cacheCauchy_nbEvalExref(cache) = q1;
    
    /* compute error cacheCauchy_wanErrExref(cache)[i] = (d*isoRatio^(-q1+i))/(1-isoRatio^(-q1)) */
    realApp_ptr wP = cacheCauchy_wanErrExref(cache) + 0;
    cacheCauchy_wantedErrorOnS0 (wP, degree, q1, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC );
    
    for (int i = 1; i<cacheCauchy_nbPwSuExref(cache); i++) {
        wP = cacheCauchy_wanErrExref(cache) + i;
        realApp_mul_realRat(wP, cacheCauchy_wanErrExref(cache) + (i-1), cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    }
    
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
       
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---cacheCauchy: \n");
        printf("#------ assumed isolation ratio                : "); realRat_print(cacheCauchy_isoRatioref(cache)); printf("\n");
        printf("#------ number of power sums  for exclusion    : %ld\n", cacheCauchy_nbPwSuExref(cache));
        printf("#------ number of eval points for exclusion    : %ld\n", cacheCauchy_nbEvalExref(cache));
        printf("#------ error on si* for exclusion             : "); realApp_printd(cacheCauchy_precfdivref(cache), 10); printf("\n");
        for (int i=0; i<cacheCauchy_nbPwSuExref(cache); i++) {
            printf("#------ error on s%d for exclusion              : ", i); realApp_printd(cacheCauchy_wanErrExref(cache) +i, 10); printf("\n");
        }
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

#ifdef CCLUSTER_TIMINGS    
    clock_t start = clock();
#endif 
    
    if ( !(realRat_cmp(radius, cacheCauchy_curRadiuref(cache))==0) ){
//     cacheCauchy_precEvalExref(cache) = 0;
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

#ifdef CCLUSTER_TIMINGS    
    time_in_cacheCauchy_set_bounds += (double) (clock() - start);
#endif
    
}

void cacheCauchy_clear ( cacheCauchy_t cache ){
    
    if (cacheCauchy_inNZCref(cache))
        ccluster_free( cacheCauchy_inNZCref(cache) );
    
    realRat_clear(cacheCauchy_isoRatioref(cache));
    realApp_clear(cacheCauchy_precfdivref(cache));
    
    for (int i=0; i< cacheCauchy_nbPwSuExref(cache); i++)
        realApp_clear( cacheCauchy_wanErrExref(cache) + i );
    ccluster_free( cacheCauchy_wanErrExref(cache) );
    
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
    
}
