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

#include "cauchy_tests/cauchy_tests.h"

void cauchyTest_computePointPointShifted( compApp_t point, compApp_t pointShifted, 
                                          const compApp_t center, slong q, slong i, const realRat_t radius,
                                          slong prec ) {

    realApp_set_si     ( compApp_realref(point), 2*i );
    realApp_div_ui     ( compApp_realref(point), compApp_realref(point), q, prec );
    realApp_zero       ( compApp_imagref(point) );
    compApp_exp_pi_i   ( point,                  point, prec);
    compApp_mul_realRat( pointShifted,           point, radius, prec );
    compApp_add        ( pointShifted,           pointShifted, center, prec );
    
}

int  cauchyTest_compute_fdiv_checkPrecAndBounds( compApp_t fdiv, 
                                                 const compApp_t fval, 
                                                 const compApp_t fderval,
                                                 const realApp_t lbound,
                                                 const realApp_t ubound,
                                                 slong prec ) {
    
    int res = 1;
    
    realApp_t modulus;
    realApp_init(modulus);
    compApp_abs(modulus, fval, prec);
    
    if (realApp_lt( modulus, lbound )){
        res=-2;
    } else if (compApp_contains_zero( fval )){
        res=-1;
    } else {
        compApp_div(fdiv, fderval, fval, prec);
        compApp_abs(modulus, fdiv, prec);
        if (realApp_gt( modulus, ubound )) {
            res = -2;
        }
    }
       
    realApp_clear(modulus);
    return res;
}

/* if (radius2==NULL), computes: 
 * c             = approx( center )
 * otherwise, computes:
 * c             = approx( center + radius2 * exp( 2*pi*vindex/vangle ) )
 * then:
 * points        = [     radius*exp(2*pi*0/nbPoints), ...,     radius*exp( 2*pi*(nbPoints-1)/nbPoints ) ]
 * pointsShifted = [ c + radius*exp(2*pi*0/nbPoints), ..., c + radius*exp( 2*pi*(nbPoints-1)/nbPoints ) ]
 */

void cauchyTest_getEvaluationPoints_fpri( const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,
                                          slong vindex,
                                          cacheCauchy_t cacheCau ) {
#ifdef CCLUSTER_TIMINGS    
    clock_t start = clock();
#endif
    
    slong nbPoints         = cacheCauchy_nbEvalExref(cacheCau);
    compApp_ptr points          = cacheCauchy_pointsExref(cacheCau);
    fpci_ptr fpci_points        = cacheCauchy_fpci_pointsExref(cacheCau);
//     compApp_ptr pointsShifted   = cacheCauchy_pointsShiftedExref(cacheCau);
    fpci_ptr fpci_pointsShifted = cacheCauchy_fpci_pointsShiftedExref(cacheCau);
    
    compApp_t c, a;
    realApp_t r;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realApp_init(r);
    realRat_init(argu);

    fpci_t cfpci;
    fpci_init(cfpci);
    fpri_t rfpri;
    fpri_init(rfpri);
    
    realApp_set_realRat(r, radius, CCLUSTER_DEFAULT_PREC);
    fpri_set_arb(rfpri, r);
    
#ifdef CCLUSTER_TIMINGS    
    clock_t start2 = clock();
#endif    
    /* recompute roots of unity only if necessary */
    if ( cacheCauchy_fpci_precEvalExref(cacheCau) == 0 ){
        
        if ( cacheCauchy_precEvalExref(cacheCau) == 0 ){
            for(slong i=0; i<nbPoints; i++) {
                realRat_set_si(argu, -2*i, nbPoints);
                compApp_set_realRat(a, argu, CCLUSTER_DEFAULT_PREC);
                acb_exp_pi_i( points + i, a, CCLUSTER_DEFAULT_PREC);
//                 fpci_set_acb( fpci_points + i, points + i );
//                 printf(" fpci_points + i: "); fpci_print( fpci_points + i ); printf("\n");
//                 printf("  acb_points + i: "); acb_printd( points + i, 16 ); printf("\n\n");
            }
            cacheCauchy_precEvalExref(cacheCau) = CCLUSTER_DEFAULT_PREC;
        }
        for(slong i=0; i<nbPoints; i++) {
            fpci_set_acb( fpci_points + i, points + i );
//             printf(" fpci_points + i: "); fpci_print( fpci_points + i ); printf("\n");
        }
        cacheCauchy_fpci_precEvalExref(cacheCau) = CCLUSTER_FPRI_PREC;
    }
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_rootsOfUnits += (double) (clock() - start2);
#endif
    if (radius2==NULL)
        compApp_set_compRat(c, center, CCLUSTER_DEFAULT_PREC);
    else {
    /* compute approximation of the center */
        compApp_set_compRat         (c,    center,   2*CCLUSTER_DEFAULT_PREC);
        realRat_set_si              (argu, 2*vindex, vangle);
        compApp_set_realRat         (a,    argu,     2*CCLUSTER_DEFAULT_PREC);
        acb_exp_pi_i                (a,    a,        2*CCLUSTER_DEFAULT_PREC);
        compApp_mul_realRat_in_place(a,    radius2,  2*CCLUSTER_DEFAULT_PREC);
        compApp_add                 (c,    c, a,     2*CCLUSTER_DEFAULT_PREC);
    }
    
    fpci_set_acb(cfpci, c);
    
#ifdef CCLUSTER_TIMINGS    
    start2 = clock();
#endif
    for(slong i=0; i<nbPoints; i++) {
//         printf(" rfpri: "); fpri_print( rfpri ); printf("\n");
//         printf(" r    : ");  arb_printd( r, 16 ); printf("\n");
//         printf(" fpci_points + i: "); fpci_print( fpci_points + i ); printf("\n");
//         printf("  acb_points + i: "); acb_printd( points + i, 16 ); printf("\n\n");
        _fpri_mul( fpci_realref(fpci_pointsShifted + i), fpci_realref(fpci_points + i), rfpri );
        _fpri_mul( fpci_imagref(fpci_pointsShifted + i), fpci_imagref(fpci_points + i), rfpri );
//         printf(" fpci_pointsShifted + i: "); fpci_print( fpci_pointsShifted + i ); printf("\n");
        fpci_add(  fpci_pointsShifted + i, fpci_pointsShifted + i, cfpci);
        
//         printf(" fpci_pointsShifted + i: "); fpci_print( fpci_pointsShifted + i ); printf("\n");
//         
//         compApp_set( pointsShifted + i, c);
//         compApp_addmul_realApp(pointsShifted + i, points + i, r, CCLUSTER_DEFAULT_PREC);
//         printf("  acb_pointsShifted + i: "); acb_printd( pointsShifted + i, 16 ); printf("\n\n");
    }
//     printf("\n\n");
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_shift_points += (double) (clock() - start2);
#endif    
    compApp_clear(c);
    compApp_clear(a);
    realApp_clear(r);
    realRat_clear(argu);
    
    fpci_clear(cfpci);
    fpri_clear(rfpri);
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_getEvaluationPoints += (double) (clock() - start);
#endif
    
}

void cauchyTest_getEvaluationPoints( const compRat_t center,
                                     const realRat_t radius,
                                     const realRat_t radius2,
                                     slong vangle,
                                     slong vindex,
                                     cacheCauchy_t cacheCau,
                                     slong prec ) {
    
    if (prec==CCLUSTER_FPRI_PREC) {
        cauchyTest_getEvaluationPoints_fpri( center, radius, radius2, vangle, vindex, cacheCau );
        return;
    }
    
#ifdef CCLUSTER_TIMINGS    
    clock_t start = clock();
#endif   
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    compApp_ptr points = cacheCauchy_pointsExref(cacheCau);
    compApp_ptr pointsShifted = cacheCauchy_pointsShiftedExref(cacheCau);
    
    compApp_t c, a;
    realApp_t r;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realApp_init(r);
    realRat_init(argu);

    realApp_set_realRat(r, radius, prec);
    
#ifdef CCLUSTER_TIMINGS    
    clock_t start2 = clock();
#endif    
    /* recompute roots of unity only if necessary */
    if ( prec > cacheCauchy_precEvalExref(cacheCau) ){
        
        for(slong i=0; i<nbPoints; i++) {
            realRat_set_si(argu, -2*i, nbPoints);
            compApp_set_realRat(a, argu, prec);
            acb_exp_pi_i( points + i, a, prec);
        }
        cacheCauchy_precEvalExref(cacheCau) = prec;
    }
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_rootsOfUnits += (double) (clock() - start2);
#endif
    
    if (radius2==NULL)
        compApp_set_compRat(c, center, prec);
    else {
    /* compute approximation of the center */
        compApp_set_compRat         (c,    center,   2*prec);
        realRat_set_si              (argu, 2*vindex, vangle);
        compApp_set_realRat         (a,    argu,     2*prec);
        acb_exp_pi_i                (a,    a,        2*prec);
        compApp_mul_realRat_in_place(a,    radius2,  2*prec);
        compApp_add                 (c,    c, a,     2*prec);
    }
#ifdef CCLUSTER_TIMINGS    
    start2 = clock();
#endif
    for(slong i=0; i<nbPoints; i++) {
//         compApp_mul_realRat(pointsShifted + i, points + i, radius, prec);
//         compApp_mul_realApp(pointsShifted + i, points + i, r, prec);
//         compApp_add( pointsShifted + i, c, pointsShifted + i, prec);
        compApp_set( pointsShifted + i, c);
        compApp_addmul_realApp(pointsShifted + i, points + i, r, prec);
    }
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_shift_points += (double) (clock() - start2);
#endif    
    compApp_clear(c);
    compApp_clear(a);
    realApp_clear(r);
    realRat_clear(argu);
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_getEvaluationPoints += (double) (clock() - start);
#endif
}

void cauchyTest_evaluateAtPoints_fpri( cacheCauchy_t cacheCau,
                                       int inCounting,
                                       metadatas_t meta, int depth){
    
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    fpci_ptr fpci_pointsShifted = cacheCauchy_fpci_pointsShiftedExref(cacheCau);
    fpci_ptr fpci_fvals         = cacheCauchy_fpci_fvalsExref(cacheCau);
    fpci_ptr fpci_fdervals      = cacheCauchy_fpci_fdervalsExref(cacheCau);
    
    clock_t start = clock();
    
    cacheCauchy_eval_fpri( fpci_fvals, fpci_fdervals, fpci_pointsShifted, nbPoints, cacheCau, meta);
    
//     for (slong i=0; i<nbPoints; i++) {
//         fpci_get_acb( cacheCauchy_fvalsExref(cacheCau) + i, fpci_fvals + i );
//         fpci_get_acb( cacheCauchy_fdervalsExref(cacheCau) + i, fpci_fdervals + i );
//     }
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoEv(meta, (double) (clock() - start));
            metadatas_add_CauchyCoEvals(meta, depth, nbPoints);
        } else {
            metadatas_add_time_CauExEv(meta, (double) (clock() - start));
            metadatas_add_CauchyExEvals(meta, depth, nbPoints);
        }
    }
}

void cauchyTest_evaluateAtPoints( cacheCauchy_t cacheCau,
                                  slong prec,
                                  int inCounting,
                                  metadatas_t meta, int depth){
    
    if (prec==CCLUSTER_FPRI_PREC) {
        cauchyTest_evaluateAtPoints_fpri( cacheCau, inCounting, meta, depth);
        return;
    }
    
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    compApp_ptr pointsShifted = cacheCauchy_pointsShiftedExref(cacheCau);
    compApp_ptr fvals = cacheCauchy_fvalsExref(cacheCau);
    compApp_ptr fdervals = cacheCauchy_fdervalsExref(cacheCau);
    
    clock_t start = clock();
    
    cacheCauchy_eval( fvals, fdervals, pointsShifted, nbPoints, cacheCau, prec, meta);
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoEv(meta, (double) (clock() - start));
            metadatas_add_CauchyCoEvals(meta, depth, nbPoints);
        } else {
            metadatas_add_time_CauExEv(meta, (double) (clock() - start));
            metadatas_add_CauchyExEvals(meta, depth, nbPoints);
        }
    }
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeFdivs_fromVals_fpri(cacheCauchy_t cacheCau, 
                                          int inCounting,
                                          metadatas_t meta){
    
    int res=1;
    
    slong nbPoints    = cacheCauchy_nbEvalExref(cacheCau);
    fpci_ptr fvals    = cacheCauchy_fpci_fvalsExref(cacheCau);
    fpci_ptr fdervals = cacheCauchy_fpci_fdervalsExref(cacheCau);
    fpci_ptr fdivs    = cacheCauchy_fpci_fdivsExref(cacheCau);
    
    clock_t start = clock();
    
    fpri_t modulus;
    fpri_init(modulus);
    
    fpri_ptr lb = cacheCauchy_fpri_lBoundApref(cacheCau);
    fpri_ptr ub = cacheCauchy_fpri_uBoundApref(cacheCau);
    
    /* compute fdivs; check if fvals contains zero and 
     *                has modulus less than lower bound*/
    for (slong i = 0; (i<nbPoints) && (res==1) ; i++) {
        fpci_sqrabs(modulus, fvals +i);
        if (fpci_contains_zero( fvals +i )){
            /* compute modulus of fvals +i */
            if (fpri_lt( modulus, lb )){
                res=-2;
            }
            else {
                res=-1;
            }
        }
        else if (fpri_lt( modulus, lb )) {
            res = -2;
        }
        else if (!fpri_ge( modulus, lb )){
            res = -1;
        }
        fpci_div(fdivs +i, fdervals + i, fvals + i);
        /* check if the ratio is less than the upper bound */
       fpci_sqrabs(modulus, fdivs +i);
        if (fpri_gt( modulus, ub )) {
            res = -2;
        }
        
        fpci_get_acb( cacheCauchy_fdivsExref(cacheCau) + i, fdivs + i );
    }

    fpri_clear(modulus);
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoDS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExDS(meta, (double) (clock() - start));
        }
    }
    
    return res;
}

int cauchyTest_computeFdivs_fromVals(cacheCauchy_t cacheCau,
                                     slong prec,
                                     int inCounting,
                                     metadatas_t meta){
    
    if (prec == CCLUSTER_FPRI_PREC) {
        return cauchyTest_computeFdivs_fromVals_fpri(cacheCau, inCounting, meta);
    }
    
    int res=1;
    
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    compApp_ptr fvals = cacheCauchy_fvalsExref(cacheCau);
    compApp_ptr fdervals = cacheCauchy_fdervalsExref(cacheCau);
    compApp_ptr fdivs = cacheCauchy_fdivsExref(cacheCau);
    
    clock_t start = clock();
    
    realApp_t modulus;
    realApp_init(modulus);
    
    realApp_ptr lb = cacheCauchy_lBoundApref(cacheCau);
    realApp_ptr ub = cacheCauchy_uBoundApref(cacheCau);
    
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
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoDS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExDS(meta, (double) (clock() - start));
        }
    }
    
    return res;
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeS0Approx_fromVals(compApp_t ps,
                                        const realRat_t radius,
                                        cacheCauchy_t cacheCau,
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta){
    
    int res=1;
    
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    realApp_ptr wP = cacheCauchy_precfdivref(cacheCau);
    compApp_ptr points = cacheCauchy_pointsExref(cacheCau);
    compApp_ptr fdivs = cacheCauchy_fdivsExref(cacheCau);
        
    clock_t start = clock();
        
    realApp_t radRe, radIm;
    realApp_init(radRe);
    realApp_init(radIm);
        
    /* compute s0 */
    compApp_mul(ps, fdivs + 0, points + 0, prec);
    for (slong i = 1; i<nbPoints; i++)
        compApp_addmul(ps, fdivs + (nbPoints + i)%nbPoints, points + (nbPoints + i)%nbPoints, prec);  
    compApp_div_si(ps, ps, nbPoints, prec);
    /* scale by the radius */
    compApp_mul_realRat_in_place(ps, radius, prec);
        
    /* check if precision is OK */
    realApp_get_rad_realApp( radRe, compApp_realref(ps) );
    realApp_get_rad_realApp( radIm, compApp_imagref(ps) );
    res = res && (realApp_lt( radRe, wP )) 
              && (realApp_lt( radIm, wP ));
    res = ((res==1)? 1:-1);
    
    realApp_clear(radRe);
    realApp_clear(radIm);
        
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }
    
    return res;
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeSsApprox_fromVals_fpri(compApp_ptr ps,
                                             const realRat_t radius,
                                             cacheCauchy_t cacheCau,
                                             int inCounting,
                                             metadatas_t meta){
#ifdef CCLUSTER_TIMINGS    
    clock_t start2 = clock();
#endif   
    
    int res=1;
    
    slong nbPoints     = cacheCauchy_nbEvalExref(cacheCau);
    slong nbPowerSums  = cacheCauchy_nbPwSuExref(cacheCau);
    realApp_ptr    wP  = cacheCauchy_precfdivref(cacheCau);
    fpci_ptr fpci_points = cacheCauchy_fpci_pointsExref(cacheCau);
    fpci_ptr fpci_fdivs  = cacheCauchy_fpci_fdivsExref(cacheCau);
//     compApp_ptr points = cacheCauchy_pointsExref(cacheCau);
//     compApp_ptr fdivs  = cacheCauchy_fdivsExref(cacheCau);
    
    clock_t start = clock();
    
    realApp_t radRe, radIm;
    realApp_init(radRe);
    realApp_init(radIm);
    fpci_t temp;
    fpci_init(temp);
    fpri_t r;
    fpri_init(r);
    realApp_set_realRat(radRe, radius, CCLUSTER_DEFAULT_PREC);
    fpri_set_arb(r, radRe);
    
//     printf("-----------------------------------------------\n");
    fpci_ptr psfpci = (fpci_ptr) ccluster_malloc ( nbPowerSums * sizeof(fpci_struct) );
    for (slong j = 0; j<nbPowerSums; j++)
        fpci_init(psfpci+j);
    
    /* compute powerSums */
    for (slong j = 0; j<nbPowerSums; j++) {
        fpci_mul(psfpci+j, fpci_fdivs + 0, fpci_points + 0);
//         compApp_mul(ps+j, fdivs + 0, points + 0, CCLUSTER_DEFAULT_PREC);
//         printf("psfpci+%ld: ", j); fpri_print(fpci_realref(psfpci+j)); printf("\n");
//         printf("psfpci+%ld: ", j); fpri_print(fpci_imagref(psfpci+j)); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_realref(ps+j),16); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_imagref(ps+j),16); printf("\n");
    }
//     printf("\n");
    for (slong i = 1; i<nbPoints; i++) {
//         printf("&&& i: %ld &&& \n", i);
        for (slong j = 0; j<nbPowerSums; j++) {
            fpci_mul(temp, fpci_fdivs + i, fpci_points + ((j+1)*i)%nbPoints);
            fpci_add(psfpci+j, psfpci+j, temp);
//             compApp_addmul(ps+j , fdivs + i, points + ((j+1)*i)%nbPoints, CCLUSTER_DEFAULT_PREC);
//             printf("psfpci+%ld: ", j); fpri_print(fpci_realref(psfpci+j)); printf("\n");
//             printf("psfpci+%ld: ", j); fpri_print(fpci_imagref(psfpci+j)); printf("\n");
//             printf("    ps+%ld: ", j); arb_printd(acb_realref(ps+j),16); printf("\n");
//             printf("    ps+%ld: ", j); arb_printd(acb_imagref(ps+j),16); printf("\n");
        }
//         printf("\n");
    }
    for (slong j = 0; j<nbPowerSums; j++) {
//         compApp_div_si(ps+j, ps+j, nbPoints, CCLUSTER_DEFAULT_PREC);
        fpci_div_si(psfpci+j, psfpci+j, nbPoints);
//         printf("psfpci+%ld: ", j); fpri_print(fpci_realref(psfpci+j)); printf("\n");
//         printf("psfpci+%ld: ", j); fpri_print(fpci_imagref(psfpci+j)); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_realref(ps+j),16); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_imagref(ps+j),16); printf("\n");
        
//         compApp_mul_realRat_in_place(ps+j, radius, CCLUSTER_DEFAULT_PREC);
        fpri_mul(fpci_realref(psfpci+j), fpci_realref(psfpci+j), r);
        fpri_mul(fpci_imagref(psfpci+j), fpci_imagref(psfpci+j), r);
//         printf("psfpci+%ld: ", j); fpri_print(fpci_realref(psfpci+j)); printf("\n");
//         printf("psfpci+%ld: ", j); fpri_print(fpci_imagref(psfpci+j)); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_realref(ps+j),16); printf("\n");
//         printf("    ps+%ld: ", j); arb_printd(acb_imagref(ps+j),16); printf("\n");
//         printf("\n");
        fpci_get_acb(ps+j, psfpci+j);
//         if (fpci_get_acb(ps+j, psfpci+j) == 0)
//             res = -1;
    }
    
    /* check if precision is OK */
    if (res==1)
        for (slong j = 0; j<nbPowerSums; j++){
            realApp_get_rad_realApp( radRe, compApp_realref(ps+j) );
            realApp_get_rad_realApp( radIm, compApp_imagref(ps+j) );
            res = res && (realApp_lt( radRe, wP )) && (realApp_lt( radIm, wP ));
            res = ((res==1)? 1:-1);
        }
        
    realApp_clear(radRe);
    realApp_clear(radIm);
    fpci_clear(temp);
    fpri_clear(r);
    for (slong j = 0; j<nbPowerSums; j++)
        fpci_clear(psfpci+j);
    ccluster_free(psfpci);
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }

#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_computeSsApprox_fromVals += (double) (clock() - start2);
#endif
    
    return res;
}

int cauchyTest_computeSsApprox_fromVals(compApp_ptr ps,
                                        const realRat_t radius,
                                        cacheCauchy_t cacheCau,
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta){
#ifdef CCLUSTER_TIMINGS    
    clock_t start2 = clock();
#endif 
    
    if (prec==CCLUSTER_FPRI_PREC) {
        return cauchyTest_computeSsApprox_fromVals_fpri(ps, radius, cacheCau, inCounting, meta);
    }
    
    int res=1;
    
    slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    slong nbPowerSums = cacheCauchy_nbPwSuExref(cacheCau);
    realApp_ptr wP = cacheCauchy_precfdivref(cacheCau);
    compApp_ptr points = cacheCauchy_pointsExref(cacheCau);
    compApp_ptr fdivs = cacheCauchy_fdivsExref(cacheCau);
    
    clock_t start = clock();
    
    realApp_t radRe, radIm;
    realApp_init(radRe);
    realApp_init(radIm);
    
    /* compute powerSums */
    for (slong j = 0; j<nbPowerSums; j++)
        compApp_mul(ps+j, fdivs + 0, points + 0, prec);
    for (slong i = 1; i<nbPoints; i++)
        for (slong j = 0; j<nbPowerSums; j++)
            compApp_addmul(ps+j , fdivs + i, points + ((j+1)*i)%nbPoints, prec);   
    for (slong j = 0; j<nbPowerSums; j++) {
        compApp_div_si(ps+j, ps+j, nbPoints, prec);
        compApp_mul_realRat_in_place(ps+j, radius, prec);
    }
    
    /* check if precision is OK */
    for (slong j = 0; j<nbPowerSums; j++){
        realApp_get_rad_realApp( radRe, compApp_realref(ps+j) );
        realApp_get_rad_realApp( radIm, compApp_imagref(ps+j) );
        res = res && (realApp_lt( radRe, wP )) && (realApp_lt( radIm, wP ));
        res = ((res==1)? 1:-1);
    }
        
    realApp_clear(radRe);
    realApp_clear(radIm);
        
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }

#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_computeSsApprox_fromVals += (double) (clock() - start2);
#endif
    
    return res;
}

cauchyTest_res cauchyTest_computeS0Approx(compApp_t ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex, 
                                          int *alreadyEvaluated,
                                          cacheCauchy_t cacheCau,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth ){
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
    
    if (*alreadyEvaluated == 0) {
        
        while (res.nbOfSol==-1){
            cauchyTest_getEvaluationPoints( center, radius,
                                                radius2, vangle, vindex,
                                                cacheCau, res.appPrec);
                
            cauchyTest_evaluateAtPoints( cacheCau, res.appPrec, inCounting, meta, depth);
        
            res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, res.appPrec, inCounting, meta);
            if ( res.nbOfSol ==-1 )
                res.appPrec = 2*res.appPrec;
        }
        *alreadyEvaluated = 1;
    }
    
    if (res.nbOfSol == -2)
        return res;
    
    /* compute approximation of Power sums */
    res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, radius, cacheCau, res.appPrec, inCounting, meta);
    
    while ( res.nbOfSol ==-1 ) {
        res.appPrec = 2*res.appPrec;
        cauchyTest_getEvaluationPoints( center, radius, 
                                            radius2, vangle, vindex,
                                            cacheCau, res.appPrec);
            
        cauchyTest_evaluateAtPoints( cacheCau, res.appPrec, inCounting, meta, depth);
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, res.appPrec, inCounting, meta);
         
        /* compute approximation of Power sums */
        res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, radius, cacheCau, res.appPrec, inCounting, meta);
    }
    
    return res;
}

cauchyTest_res cauchyTest_computeSsApprox(compApp_ptr ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,
                                          cacheCauchy_t cacheCau,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth ){
    
//     slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
#ifdef CCLUSTER_TIMINGS    
    clock_t start = clock();
#endif  
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
        
    while (res.nbOfSol==-1){
        /* compute points and evals at prec res.appPrec*/
        cauchyTest_getEvaluationPoints( center, radius,
                                    radius2, vangle, vindex,
                                    cacheCau, res.appPrec);
        
        cauchyTest_evaluateAtPoints( cacheCau, res.appPrec, inCounting, meta, depth);
        
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, res.appPrec, inCounting, meta);
        if ( res.nbOfSol ==-1 )
//             res.appPrec = 2*res.appPrec;
            res.appPrec = ( res.appPrec==CCLUSTER_FPRI_PREC? CCLUSTER_DEFAULT_PREC : 2*res.appPrec);
    }
    
    if (res.nbOfSol == -2)
        return res;
    
    /* compute approximation of Power sums */
    res.nbOfSol = cauchyTest_computeSsApprox_fromVals(ps, radius, cacheCau, res.appPrec, inCounting, meta);
    
    while ( res.nbOfSol ==-1 ) {
//         res.appPrec = 2*res.appPrec;
        res.appPrec =  ( res.appPrec==CCLUSTER_FPRI_PREC? CCLUSTER_DEFAULT_PREC : 2*res.appPrec);
        
        cauchyTest_getEvaluationPoints( center, radius, 
                                        radius2, vangle, vindex,
                                        cacheCau, res.appPrec);
        
        cauchyTest_evaluateAtPoints( cacheCau, res.appPrec, inCounting, meta, depth);
        /* compute approximation of Power sums */
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, res.appPrec, inCounting, meta);
        res.nbOfSol = cauchyTest_computeSsApprox_fromVals(ps, radius, cacheCau, res.appPrec, inCounting, meta);
    }

#ifdef CCLUSTER_TIMINGS    
    time_in_cauchyTest_computeSsApprox += (double) (clock() - start);
#endif
    
    return res;
    
}

void cauchyTest_computeBounds ( realApp_t ubound, realApp_t lbound, 
                                const realRat_t theta, const realRat_t r, slong d, slong prec ) {
    
    realApp_set_realRat( ubound, theta,         prec);
    realApp_add_si     ( ubound, ubound, -1,    prec);
    realApp_mul_realRat( ubound, ubound, r,     prec);
    realApp_div_realRat( ubound, ubound, theta, prec);
    realApp_pow_ui     ( lbound, ubound, d,     prec);
    realApp_inv        ( ubound, ubound,        prec);
    realApp_mul_si     ( ubound, ubound, d,     prec);

} 

slong cauchyTest_computeS0compDsk( const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   cacheCauchy_t cacheCau,
                                   metadatas_t meta, int depth) {
    
    slong res = -1;
    
    clock_t start = clock();
    double evalTime=0;
    
    slong appPrec = CCLUSTER_DEFAULT_PREC;
    /* want |s0*-s0| less than 1/4 */
    /* => number of evaluation points = ceil ( log_isoRatio (4*degree +1) ) + 1*/
    slong q = cacheCauchy_get_NbOfEvalPoints( cacheCauchy_degreeref(cacheCau), isoRatio, 1, CCLUSTER_DEFAULT_PREC );
    
    /* want w(s0*)<1/4 */
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_d(errAp, 0.25);
    
    compApp_t point       ;
    compApp_t pointShifted;
    compApp_t fval        ;
    compApp_t fderval     ;
    compApp_t fdiv        ;
    compApp_t s0star      ;
    
    compApp_init( point        );
    compApp_init( pointShifted );
    compApp_init( fval         );
    compApp_init( fderval      );
    compApp_init( fdiv         );
    compApp_init( s0star       );
    
    realApp_t radRe;
    realApp_t radIm;
    realApp_init(radRe);
    realApp_init(radIm);
    
    compApp_t c;
    
    compApp_init(c);
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=3)
            printf("#------ appPrec: %ld\n", appPrec);
        
        enoughPrec = 1;
        compApp_zero(s0star);
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        for(slong i=0; i<q; i++) {
            cauchyTest_computePointPointShifted( point, pointShifted, c, q, i, compDsk_radiusref(Delta), appPrec );
            
            start = clock();
            cacheCauchy_eval( fval, fderval, pointShifted, 1, cacheCau, appPrec, meta );
            evalTime += (double) (clock() - start);
            
            if (compApp_contains_zero( fval )){
                if (metadatas_getVerbo(meta)>=3)
                    printf("#------ fval contains 0, i: %ld\n", i);
                appPrec = 2*appPrec;
                enoughPrec = 0;
                break;
            }
            
            compApp_div(fdiv, fderval, fval, appPrec);
            compApp_mul(fdiv, fdiv, point, appPrec);
            compApp_mul_realRat(fdiv, fdiv, compDsk_radiusref(Delta), appPrec);
            compApp_div_si( fdiv, fdiv, q, appPrec);
            compApp_add(s0star, s0star, fdiv, appPrec);
            
            /* check error of s0star */
            realApp_get_rad_realApp( radRe, compApp_realref(s0star) );
            realApp_get_rad_realApp( radIm, compApp_imagref(s0star) );
//             realApp_mul_realRat( radRe, radRe, compDsk_radiusref(Delta), appPrec );
//             realApp_mul_realRat( radIm, radIm, compDsk_radiusref(Delta), appPrec );
            if ( realApp_ge( radRe, errAp ) || realApp_ge( radIm, errAp ) ){
                if (metadatas_getVerbo(meta)>=3)
                    printf("#------ s0star not precise enough, i: %ld\n", i);
                appPrec = 2*appPrec;
                enoughPrec = 0;
                break;
            }
                
        }
        if (metadatas_getVerbo(meta)>=3)
            if (enoughPrec == 1)
                printf("#------ s0star precise enough\n");
        
    }
    
    /* add error */
    realApp_add_error( compApp_realref(s0star), errAp );
    realApp_add_error( compApp_imagref(s0star), errAp );
    
    /* s0star contains a unique integer */
    realApp_get_unique_si( &res, compApp_realref(s0star) );
    
//     realApp_clear(liR);
//     realApp_clear(qApp);
    realApp_clear(errAp);
    compApp_clear( point        );
    compApp_clear( pointShifted );
    compApp_clear( fval         );
    compApp_clear( fderval      );
    compApp_clear( fdiv         );
    compApp_clear( s0star       );
    realApp_clear(radRe);
    realApp_clear(radIm);
    compApp_clear(c);
    
//     if (metadatas_haveToCount(meta)) {
//         metadatas_add_time_CauCoED(meta, evalTime);
//         metadatas_add_CauchyCoEvalsD(meta, depth, q);
//     }
    
    return res;
}

/* computes e = m*theta^(-h) + (d-m)*theta^(h)             */
/*          q = log ( e / eps ) / log ( theta ) + 1 */
void  cauchyTest_compute_required_numOfPoints ( realApp_t e, realApp_t q,
                                                const realApp_t theta_h, 
                                                const realApp_t eps,
                                                const realApp_t logtheta, 
                                                slong m, slong d,
                                                slong prec ) {
    
    realApp_mul_si     (e, theta_h,    d-m,               prec);
    realApp_inv        (q, theta_h,                       prec);
    realApp_mul_si     (q, q,          m,                 prec);
    realApp_add        (e, e,          q,                 prec);
    
    realApp_div        (q, e,          eps,               prec);
    realApp_add_si     (q, q,          1,                 prec);
    realApp_log        (q, q,                             prec);
    realApp_div        (q, q,          logtheta,          prec);
    
}

/* want to compute sN, s2N, ..., sgN, ..., shN with max error eps */
slong  cauchyTest_compute_q_and_errors ( realApp_ptr errors,
                                         const realRat_t theta, 
                                         const realRat_t eps, 
                                         slong m, slong d, slong h, slong N,
                                         slong prec ) {
    
    realApp_t epsApp, qApp, thetaTotheN, thetaTothegN, logtheta;
    realApp_init(epsApp);
    realApp_init(qApp);
    realApp_init(thetaTotheN);
    realApp_init(thetaTothegN);
    realApp_init(logtheta);
    
    realApp_set_realRat(epsApp,       eps,            prec);
    realApp_set_realRat(logtheta,     theta,          prec);
    realApp_log        (logtheta,     logtheta,       prec);
    realApp_set_realRat(thetaTotheN,  theta,          prec);
    realApp_pow_ui     (thetaTotheN,  thetaTotheN, N, prec);
    realApp_set        (thetaTothegN, thetaTotheN);
    
    slong qmax = 0;
    
    for (slong g = 1; g <= h; g++) {
        cauchyTest_compute_required_numOfPoints ( errors + (g-1), qApp, thetaTothegN, epsApp, logtheta, m, d, prec);
        slong q = realApp_ceil_si( qApp, prec ) + 1;
        qmax = CCLUSTER_MAX( qmax, q );
        
        realApp_mul(thetaTothegN, thetaTothegN,  thetaTotheN,  prec);
    }
    
    qmax = CCLUSTER_MAX( qmax, h*N + 1 );
    
    realApp_set_realRat( thetaTotheN, theta, prec);
    realApp_pow_ui     ( thetaTotheN, thetaTotheN, qmax, prec);
    realApp_add_si     ( thetaTotheN, thetaTotheN, -1, prec);
    for (slong g = 1; g <= h; g++) {
        realApp_div( errors + (g-1), errors + (g-1), thetaTotheN, prec);
    }
    
    realApp_clear(epsApp);
    realApp_clear(qApp);
    realApp_clear(thetaTotheN);
    realApp_clear(thetaTothegN);
    realApp_clear(logtheta);
    
    return qmax;
}

slong cauchyTest_computeS1compDsk( compDsk_t res,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   metadatas_t meta, int depth) {
    
    int level = 4;
    slong appPrec = CCLUSTER_DEFAULT_PREC;
    
    realRat_t error;
    realRat_init(error);
    realRat_div_ui(error, eps, 2);
    realRat_div(error, error, compDsk_radiusref(Delta) );
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_realRat(errAp, error, CCLUSTER_DEFAULT_PREC);    /* errAp = contains error */
    
//     if (metadatas_getVerbo(meta)>=level) {
//         printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: nb of eval points: %ld\n", q);
//         printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: error on S1: ");
//         realApp_printd(errorS1, 10);
//         printf("\n");
//     }
    
    realApp_t qApp, logIsoRatio;
    realApp_init(qApp);
    realApp_init(logIsoRatio);
    realApp_set_realRat(logIsoRatio, isoRatio,                                               CCLUSTER_DEFAULT_PREC);
    realApp_log        (logIsoRatio, logIsoRatio,                                            CCLUSTER_DEFAULT_PREC);
//     realApp_set_realRat(qApp,        compDsk_radiusref(Delta),                               CCLUSTER_DEFAULT_PREC);
    realApp_one        (qApp);
    realApp_mul_realRat(qApp,        qApp,                     isoRatio,                     CCLUSTER_DEFAULT_PREC);
    realApp_mul_si     (qApp,        qApp,                     cacheCauchy_degreeref(cacheCau), CCLUSTER_DEFAULT_PREC);
    realApp_div_realRat(qApp,        qApp,                     error,                        CCLUSTER_DEFAULT_PREC);
    /* qApp = (rd*isoRatio)/error */
    realApp_add_si     (qApp,        qApp,                     1,                            CCLUSTER_DEFAULT_PREC);
    /* qApp = ( 1 + (rd*isoRatio)/error ) */
    realApp_log        (qApp,        qApp,                                                   CCLUSTER_DEFAULT_PREC);
    /* qApp = log_{e}       ( 1 + (rd*isoRatio)/error ) */
    realApp_div        (qApp,        qApp,                     logIsoRatio,                  CCLUSTER_DEFAULT_PREC); 
    /* qApp = log_{isoRatio}( 1 + (rd*isoRatio)/error ) */
    slong q = realApp_ceil_si( qApp, CCLUSTER_DEFAULT_PREC ) + 1;
    /* qApp = ceil( log_{isoRatio}(1 + (rd*isoRatio)/error ) ) + 1 */
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: nb of eval points: %ld\n", q);
    }
    
    compApp_t point       ;
    compApp_t pointShifted;
    compApp_t fval        ;
    compApp_t fderval     ;
    compApp_t fdiv        ;
    compApp_t s1          ;
    
    compApp_init( point        );
    compApp_init( pointShifted );
    compApp_init( fval         );
    compApp_init( fderval      );
    compApp_init( fdiv         );
    compApp_init( s1           );
    
    realApp_t radRe;
    realApp_t radIm;
    realApp_init(radRe);
    realApp_init(radIm);
    
    compApp_t c;
    compApp_init(c);
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: precision: %ld\n",  appPrec);
        }
        
        enoughPrec = 1;
        compApp_zero(s1);
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        for(slong i=0; i<q; i++) {
            cauchyTest_computePointPointShifted( point, pointShifted, c, q, i, compDsk_radiusref(Delta), appPrec );
            cacheCauchy_eval( fval, fderval, pointShifted, 1, cacheCau, appPrec, meta );
            
            if (compApp_contains_zero( fval )){
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: fval contains 0, i: %ld\n", i);
                }
                appPrec = 2*appPrec;
                enoughPrec = 0;
                break;
            }
            
            compApp_div(fdiv, fderval, fval, appPrec);
            compApp_pow_si(point, point, 2, appPrec);
            compApp_mul(fdiv, fdiv, point, appPrec);
            compApp_mul_realRat(fdiv, fdiv, compDsk_radiusref(Delta), appPrec);
            compApp_mul_realRat(fdiv, fdiv, compDsk_radiusref(Delta), appPrec);
            compApp_div_si( fdiv, fdiv, q, appPrec);
            compApp_add(s1, s1, fdiv, appPrec);
            
            /* check error of s1 */
            realApp_get_rad_realApp( radRe, compApp_realref(s1) );
            realApp_get_rad_realApp( radIm, compApp_imagref(s1) );
//             realApp_mul_realRat( radRe, radRe, compDsk_radiusref(Delta), appPrec );
//             realApp_mul_realRat( radIm, radIm, compDsk_radiusref(Delta), appPrec );
            if ( realApp_ge( radRe, errAp ) || realApp_ge( radIm, errAp ) ){
                if (metadatas_getVerbo(meta)>=level)
                    printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: s1 not precise enough, i: %ld\n", i);
                appPrec = 2*appPrec;
                enoughPrec = 0;
                break;
            }
                
        }
        if (metadatas_getVerbo(meta)>=level)
            if (enoughPrec == 1)
                printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: s1 precise enough\n");
            
        
    }
//     if (metadatas_getVerbo(meta)>=3) {
//         compApp_t temp;
//         compApp_init(temp);
//         compApp_div_realRat(temp, s1, compDsk_radiusref(Delta), appPrec);
//         printf(" s1: "); compApp_printd(temp, 10); printf("\n");
//         compApp_clear(temp);
//     }
        
    /* the radius of the ball s1 is less than errAp => it is less than error = eps/2 */
    /* s1 approximates rs1(p_\Delta,0,1) with |s1 - rs1(p_\Delta,0,1)|\leq errAp */
    /* thus the disk center(s1), 2*error contains rs1(p_\Delta,0,1) */
    compRat_t mc;
    compRat_init(mc);
    compRat_set(mc, compDsk_centerref(Delta));
    realRat_mul_si( compRat_realref(mc), compRat_realref(mc), nbOfRoots);
    realRat_mul_si( compRat_imagref(mc), compRat_imagref(mc), nbOfRoots);
    compApp_get_compRat( compDsk_centerref(res), s1 );
    compRat_add(compDsk_centerref(res), compDsk_centerref(res), mc);
//     realRat_set_si(compDsk_radiusref(res), 2, 1);
//     realRat_pow_si(compDsk_radiusref(res), compDsk_radiusref(res), -prec);
    realRat_set(compDsk_radiusref(res), eps);
    
    realRat_clear(error);
    realApp_clear(errAp);
    realApp_clear(qApp);
    realApp_clear(logIsoRatio);
    
    compApp_clear( point       );
    compApp_clear( pointShifted);
    compApp_clear( fval        );
    compApp_clear( fderval     );
    compApp_clear( fdiv        );
    compApp_clear( s1           );
    realApp_clear(radRe);
    realApp_clear(radIm);
    
    compApp_clear(c);
    compRat_clear(mc);
    
    return appPrec;
}

/* Assume Delta is isoRatio-isolated and contains m roots     */
/* Let h be a non-negative integer                            */
/* Assume SS is initialized to contain at least h numbers */
/* Computes h first power sums  s1, s2, ..., sh with error less than eps   */
cauchyTest_res cauchyTest_computeSScompDsk( compApp_ptr SS,
                                            const realRat_t isoRatio,
                                            const compDsk_t Delta,
                                            slong m,
                                            slong h,
                                            cacheCauchy_t cacheCau,
                                            const realRat_t eps,
                                            slong prec,
                                            metadatas_t meta, int depth){
    
    int level = 4;
    slong appPrec = CCLUSTER_DEFAULT_PREC;
    slong degree  = cacheCauchy_degreeref(cacheCau);
    /* number of power sums to compute s1, ..., sh: h */

    /* compute the required number of points */
    realRat_t error;
    realRat_init(error);
    realRat_div_ui(error, eps, 2);
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_realRat(errAp, error, CCLUSTER_DEFAULT_PREC);    /* errAp = contains error */
    
    /* initialize the vector of errors on the sjs */
    realApp_ptr errorsSS = (realApp_ptr) ccluster_malloc ( h*sizeof(realApp));
    for (slong j = 0; j < h; j++)
        realApp_init( errorsSS + j );
    
    slong qmax = cauchyTest_compute_q_and_errors ( errorsSS, isoRatio, error, m, degree, h, 1, CCLUSTER_DEFAULT_PREC );
    
    realApp_t ubound, lbound;
    realApp_init( ubound );
    realApp_init( lbound );
    cauchyTest_computeBounds ( ubound, lbound, isoRatio, compDsk_radiusref(Delta), degree, CCLUSTER_DEFAULT_PREC );
    if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: lb:");
            realApp_printd( lbound, 10);
            printf("\n");
            printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: ub:");
            realApp_printd( ubound, 10);
            printf("\n");
    }
    /* compute the power sums */
    compApp_t fval, fderval, fdiv;
    realApp_t radRe, radIm;
    compApp_t c;
    
    compApp_init(fval         );
    compApp_init(fderval      );
    compApp_init(fdiv         );
    realApp_init(radRe);
    realApp_init(radIm);
    compApp_init(c);
    
    int enoughPrec = -1;
    /* compute the points */
    while (enoughPrec == -1) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: precision: %ld\n",  appPrec);
        }
        
        enoughPrec = 1;
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        compApp_ptr points = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        compApp_ptr pointsShifted = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        
        for(slong i=0; i<qmax; i++) {
            compApp_init (points + i);
            compApp_init (pointsShifted + i);
            cauchyTest_computePointPointShifted( points + i, pointsShifted + i, c, qmax, i, compDsk_radiusref(Delta), appPrec );
        }
        
        /* initialize power sums */
        for (slong j = 0; j<h ; j++)
            compApp_zero( SS + j );
        
        for(slong i=0; (i<qmax) && (enoughPrec==1) ; i++) {
            
            cacheCauchy_eval( fval, fderval, pointsShifted + i, 1, cacheCau, appPrec, meta );
            enoughPrec = cauchyTest_compute_fdiv_checkPrecAndBounds( fdiv, fval, fderval, lbound, ubound, appPrec );
            if (enoughPrec==1) {
                for (slong j = 0; j< h; j++) 
                        compApp_addmul(SS+j , fdiv, points + ((j+2)*i)%qmax, appPrec);
            } else {
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: return of compute fdiv for %ld-th point: %d\n", i, enoughPrec);
                }
            }
            
        }
        
        if (enoughPrec==1) {
            for (slong j = 0; j<h ; j++) {
                compApp_div_si(SS+j, SS+j, qmax, appPrec);
                compApp_mul_realRat(SS+j, SS+j, compDsk_radiusref(Delta), appPrec);
            }
            /* no need to scale by the radius because points are already */
            
            /* check if precisions are ok */
            for (slong j = 0; j<h ; j++) {
                realApp_get_rad_realApp( radRe, compApp_realref(SS+j) );
                realApp_get_rad_realApp( radIm, compApp_imagref(SS+j) );
                if ( (realApp_ge( radRe, errAp )) || (realApp_ge( radIm, errAp )) ) {
                    enoughPrec=-1;
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: s%ld not precise enough:",j);
                        compApp_printd( SS+j, 10);
                        printf("\n");
                    }
                } else {
                    realApp_add_error( compApp_realref(SS+j), errorsSS + j );
                    realApp_add_error( compApp_imagref(SS+j), errorsSS + j );
                }
            }
        }
        
        if (enoughPrec==-1) {
            appPrec = 2*appPrec;
        } else {
            if (metadatas_getVerbo(meta)>=level) {
                printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: SS precise enough\n");
            }
        }
        
        for ( slong i =0; i<qmax; i++ ) {
            compApp_clear(points + i);
            compApp_clear(pointsShifted + i);
        }
        
        ccluster_free(points);
        ccluster_free(pointsShifted);
        
    }
    
    
    for (slong i = 0; i <h ; i++)
        realApp_clear( errorsSS + i );
    ccluster_free(errorsSS);
    
    realRat_clear(error);
    realApp_clear(errAp);
//     realApp_clear(qApp);
//     realApp_clear(IsoRatioTotheN);
//     realApp_clear(IsoRatioInvTotheN);
//     realApp_clear(logIsoRatio);
    
    realApp_clear( ubound );
    realApp_clear( lbound );
    
    compApp_clear( fval         );
    compApp_clear( fderval      );
    compApp_clear( fdiv         );
    realApp_clear(radRe);
    realApp_clear(radIm);
    compApp_clear(c);
    
    cauchyTest_res res;
    res.appPrec = appPrec;
    res.nbOfSol = enoughPrec;
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: power sums:\n");
        for (slong j=0; j<h ; j++) {
            printf("#--------------- %ld-th: ", j); compApp_printd( SS+j, 10); printf("\n");
        }
    }
    
    return res;
}

/* Assume Delta is isoRatio-isolated and contains m roots           */
/* Let N be a non-negative integer                                  */
/* Assume SgN is initialized to contain at least h numbers          */
/* Computes power sums sN, s2N, ..., shN with error less than eps   */
cauchyTest_res cauchyTest_computeSgNcompDsk( compApp_ptr SgN,
                                             const realRat_t isoRatio,
                                             const compDsk_t Delta,
                                             slong m,
                                             ulong N,
                                             slong h,
                                             cacheCauchy_t cacheCau,
                                             const realRat_t eps,
                                             slong prec,
                                             metadatas_t meta, int depth){
    
    int level = 3;
    slong appPrec = prec;
    slong degree  = cacheCauchy_degreeref(cacheCau);
    /* compute power sums sN, ..., shN: h */

    /* compute the required number of points */
    realRat_t error;
    realRat_init(error);
    realRat_div_ui(error, eps, 2);
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_realRat(errAp, error, CCLUSTER_DEFAULT_PREC);    /* errAp = contains error */
    
    realApp_ptr errorsSgN = (realApp_ptr) ccluster_malloc ( m*sizeof(realApp));
    for (slong g = 1; g <= h; g++)
        realApp_init( errorsSgN + (g-1) );
    
    slong qmax = cauchyTest_compute_q_and_errors ( errorsSgN, isoRatio, error, m, degree, h, N, CCLUSTER_DEFAULT_PREC );
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: number of eval points: %ld\n", qmax);
    }
    
    realApp_t ubound, lbound;
    realApp_init( ubound );
    realApp_init( lbound );
    cauchyTest_computeBounds ( ubound, lbound, isoRatio, compDsk_radiusref(Delta), degree, CCLUSTER_DEFAULT_PREC );
//     if (metadatas_getVerbo(meta)>=level) {
//             printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: lb:");
//             realApp_printd( lbound, 10);
//             printf("\n");
//             printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: ub:");
//             realApp_printd( ubound, 10);
//             printf("\n");
//     }
    /* compute the power sums */
    compApp_t fval, fderval, fdiv, temp;
    realApp_t radRe, radIm;
    compApp_t c;
    
    compApp_init(fval         );
    compApp_init(fderval      );
    compApp_init(fdiv         );
    realApp_init(radRe);
    realApp_init(radIm);
    compApp_init(temp);
    compApp_init(c);
    
    int enoughPrec = -1;
    /* compute the points */
    while (enoughPrec == -1) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: precision: %ld\n",  appPrec);
        }
        
        enoughPrec = 1;
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        compApp_ptr points = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        compApp_ptr pointsShifted = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        
        for(slong i=0; i<qmax; i++) {
            compApp_init (points + i);
            compApp_init (pointsShifted + i);
            cauchyTest_computePointPointShifted( points + i, pointsShifted + i, c, qmax, i, compDsk_radiusref(Delta), appPrec );
        }
        
        /* initialize power sums */
        for (slong g = 1; g <= h ; g++)
            compApp_zero( SgN + (g-1) );
        
        for(slong i=0; (i<qmax) && (enoughPrec==1) ; i++) {
            cacheCauchy_eval( fval, fderval, pointsShifted + i, 1, cacheCau, appPrec, meta );
            enoughPrec = cauchyTest_compute_fdiv_checkPrecAndBounds( fdiv, fval, fderval, lbound, ubound, appPrec );
            if (!(enoughPrec==1)) {
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: return of compute fdiv for %ld-th point: %d\n", i, enoughPrec);
                }
            } else {
                for (slong g = 1; g <= h && (enoughPrec==1); g++) {
//                         compApp_addmul(SgN+(g-1) , fdiv, points + ((g*N+1)*i)%qmax, appPrec); 
                        compApp_mul( temp, fdiv, points + ((g*N+1)*i)%qmax, appPrec);
                        compApp_div_si(temp, temp, qmax, appPrec);
                        compApp_mul_realRat(temp, temp, compDsk_radiusref(Delta), appPrec);
                        compApp_add( SgN+(g-1), SgN+(g-1), temp, appPrec ); 
                        /* check precision */
                        realApp_get_rad_realApp( radRe, compApp_realref(SgN + (g-1)) );
                        realApp_get_rad_realApp( radIm, compApp_imagref(SgN + (g-1)) );
                        if ( (realApp_ge( radRe, errAp )) || (realApp_ge( radIm, errAp )) ) {
                            enoughPrec=-1;
                            if (metadatas_getVerbo(meta)>=level) {
                                printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: s%ld not precise enough after %ld evals:",(g*N), i);
                                compApp_printd( SgN + (g-1), 10);
                                printf("\n");
                            }
                        }
                    }
            }
        }
        
        if (enoughPrec==1) {
            
            /* add errors */
            for (slong g = 1; g<=h ; g++) {
                realApp_add_error( compApp_realref(SgN + (g-1)), errorsSgN + (g-1) );
                realApp_add_error( compApp_imagref(SgN + (g-1)), errorsSgN + (g-1) );
            }
        }
        
        if (enoughPrec==-1) {
            appPrec = 2*appPrec;
        } else {
            if (metadatas_getVerbo(meta)>=level) {
                printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: SgN precise enough\n");
            }
        }
        
        for ( slong i =0; i<qmax; i++ ) {
            compApp_clear(points + i);
            compApp_clear(pointsShifted + i);
        }
        
        ccluster_free(points);
        ccluster_free(pointsShifted);
        
    }
    
    
    for (slong g = 1; g <= h ; g++)
        realApp_clear( errorsSgN + (g-1) );
    ccluster_free(errorsSgN);
    
    realRat_clear(error);
    realApp_clear(errAp);
    
    realApp_clear( ubound );
    realApp_clear( lbound );
    
    compApp_clear( fval         );
    compApp_clear( fderval      );
    compApp_clear( fdiv         );
    realApp_clear(radRe);
    realApp_clear(radIm);
    compApp_clear(temp);
    compApp_clear(c);
    
    cauchyTest_res res;
    res.appPrec = appPrec;
    res.nbOfSol = enoughPrec;
    
//     if (metadatas_getVerbo(meta)>=level) {
//         printf("#------------cauchy_tests.c: cauchyTest_computeSgNcompDsk: power sums:\n");
//         for (slong g=1; g<=m ; g++) {
//             printf("#--------------- %ld-th: ", g*N); compApp_printd( SgN+(g-1), 10); printf("\n");
//         }
//     }
    
    return res;
}

/* DEPRECATED */
// void cauchyTest_shiftFFT(                 const compRat_t center,
//                                           const realRat_t radius,
//                                           const realRat_t radius2,
//                                           slong vangle,           
//                                           slong vindex,           
//                                           cacheApp_t cache,
//                                           cacheCauchy_t cacheCau,
//                                           int certified,
//                                           slong prec,
//                                           int inCounting,
//                                           metadatas_t meta, int depth ) {
//     
//     // //             printf("here\n");
//     
//                 clock_t start = clock();
//                 
//                 slong nbPoints = cacheCauchy_nbEvalCeref(cacheCau);
//                 
//                 compApp_poly_t shiftedPoly;
//                 compApp_poly_init2( shiftedPoly, nbPoints );
//                 /* compute approximation of the pol */
//                 compApp_poly_set(shiftedPoly, cacheApp_getApproximation ( cache, prec ) );
//                 
//                 compApp_t c;
//                 compApp_init(c);
//                 realRat_t argu;
//                 realRat_init(argu);
//                 compApp_t a;
//                 compApp_init(a);
//                 /* compute approximation of the center */
//                 if (radius2==NULL)
//                     compApp_set_compRat(c, center, 2*prec);
//                 else {
//                     compApp_set_compRat         (c,    center,   2*prec);
//                     realRat_set_si              (argu, 2*vindex, vangle);
//                     compApp_set_realRat         (a,    argu,     2*prec);
//                     acb_exp_pi_i                (a,    a,        2*prec);
//                     compApp_mul_realRat_in_place(a,    radius2,  2*prec);
//                     compApp_add                 (c,    c, a,     2*prec);
//                 }
//                 /* shift in the center */
//                 clock_t start2 = clock();
//                 _acb_poly_taylor_shift_convolution(shiftedPoly->coeffs, c, shiftedPoly->length, prec);
//                 compApp_poly_scale_realRat_in_place( shiftedPoly->coeffs, radius, shiftedPoly->length, prec );
//                 if (metadatas_haveToCount(meta))
//                     metadatas_add_time_Taylors(meta, (double) (clock() - start2) );
//                 
//                 /* compute fft */
//                 acb_dft_pre_t t;
//                 acb_dft_precomp_init(t, nbPoints, prec);
//                 acb_dft_precomp(cacheCauchy_fvalsCeref(cacheCau), shiftedPoly->coeffs, t, prec);
//                 /* derivative */
//                 compApp_poly_derivative( shiftedPoly, shiftedPoly, prec);
//                 slong index = 0;
//                 while (index < shiftedPoly->length){
//                     compApp_div_realRat( (shiftedPoly->coeffs) + index, (shiftedPoly->coeffs) + index, radius, prec );
//                     index++;
//                 }
//                 compApp_zero((shiftedPoly->coeffs) + (shiftedPoly->length));
//                 acb_dft_precomp(cacheCauchy_fdervalsCeref(cacheCau), shiftedPoly->coeffs, t, prec);
//                 
//                 compApp_poly_clear(shiftedPoly);
//                 compApp_clear(c);
//                 realRat_clear(argu);
//                 compApp_clear(a);
//                 acb_dft_precomp_clear(t);
//                 
//     if (metadatas_haveToCount(meta)) {
//         slong nbEvals = nbPoints;
//         if (inCounting == CAUCHYTEST_INCOUNTIN) {
//             if (certified == CAUCHYTEST_CERTIFIED) {
//                 metadatas_add_time_CauCoED(meta, (double) (clock() - start));
//                 metadatas_add_CauchyCoEvalsD(meta, depth, nbEvals);
//             } else {
//                 metadatas_add_time_CauCoEP(meta, (double) (clock() - start));
//                 metadatas_add_CauchyCoEvalsP(meta, depth, nbEvals);
//             }
//         } else {
//             if (certified == CAUCHYTEST_CERTIFIED) {
//                 metadatas_add_time_CauExED(meta, (double) (clock() - start));
//                 metadatas_add_CauchyExEvalsD(meta, depth, nbEvals);
//             } else {
//                 metadatas_add_time_CauExEP(meta, (double) (clock() - start));
//                 metadatas_add_CauchyExEvalsP(meta, depth, nbEvals);
//             }
//         }
//     }
//     
// }
