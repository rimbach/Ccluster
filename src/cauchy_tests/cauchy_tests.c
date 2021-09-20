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

/* if (radius2==NULL), computes: 
 * c             = approx( center )
 * otherwise, computes:
 * c             = approx( center + radius2 * exp( 2*pi*vindex/vangle ) )
 * then:
 * points        = [     radius*exp(2*pi*0/nbPoints), ...,     radius*exp( 2*pi*(nbPoints-1)/nbPoints ) ]
 * pointsShifted = [ c + radius*exp(2*pi*0/nbPoints), ..., c + radius*exp( 2*pi*(nbPoints-1)/nbPoints ) ]
 */

void cauchyTest_getEvaluationPoints( const compRat_t center,
                                     const realRat_t radius,
                                     const realRat_t radius2,
                                     slong vangle,
                                     slong vindex,
                                     cacheCauchy_t cacheCau,
                                     int certified,
                                     slong prec ) {
    
    slong nbPoints;
    compApp_ptr points;
    compApp_ptr pointsShifted;
    
    if (certified == CAUCHYTEST_CERTIFIED) {
        nbPoints      = cacheCauchy_nbEvalCeref(cacheCau);
        points        = cacheCauchy_pointsCeref(cacheCau);
        pointsShifted = cacheCauchy_pointsShiftedCeref(cacheCau);
    } else {
        nbPoints      = cacheCauchy_nbEvalExref(cacheCau);
        points        = cacheCauchy_pointsExref(cacheCau);
        pointsShifted = cacheCauchy_pointsShiftedExref(cacheCau);
    }
    
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
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
    for(slong i=0; i<nbPoints; i++) {
        realRat_set_si(argu, -2*i, nbPoints);
        compApp_set_realRat(a, argu, prec);
        acb_exp_pi_i( points + i, a, prec);
        compApp_mul_realRat_in_place(points + i, radius, prec);
        compApp_add( pointsShifted + i, c, points + i, prec);
    }
    
    compApp_clear(c);
    compApp_clear(a);
    realRat_clear(argu);
    
}

void cauchyTest_evaluateAtPoints( cacheApp_t cache,
                                  cacheCauchy_t cacheCau,
                                  int certified,
                                  slong prec,
                                  int inCounting,
                                  metadatas_t meta, int depth){
    
    slong nbPoints;
    compApp_ptr pointsShifted;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    int reUseEvals = 0;
    slong quo = 0;
    
    if (certified == CAUCHYTEST_CERTIFIED) {
        nbPoints      = cacheCauchy_nbEvalCeref(cacheCau);
        pointsShifted = cacheCauchy_pointsShiftedCeref(cacheCau);
        fvals         = cacheCauchy_fvalsCeref(cacheCau);
        fdervals      = cacheCauchy_fdervalsCeref(cacheCau);
        reUseEvals    = (cacheCauchy_precEvalExref(cacheCau) == prec);
        quo           = cacheCauchy_quotientref(cacheCau);
    } else {
        nbPoints      = cacheCauchy_nbEvalExref(cacheCau);
        pointsShifted = cacheCauchy_pointsShiftedExref(cacheCau);
        fvals         = cacheCauchy_fvalsExref(cacheCau);
        fdervals      = cacheCauchy_fdervalsExref(cacheCau);
    }
    
    clock_t start = clock();
    
    if (cacheCauchy_evalFastref(cacheCau) == NULL) {
        
        if (cacheCauchy_choiceref(cacheCau) == -1) {
            
            clock_t start2 = clock();
            for (slong i=0; i<nbPoints; i++)
                cacheCauchy_rectangularEval ( fvals + i, fdervals + i, cacheCau, cache, pointsShifted + i, prec);
            double timeInRectangular = (double) (clock() - start2);
            start2 = clock();
            for (slong i=0; i<nbPoints; i++)
                cacheCauchy_sparseEval ( fvals + i, fdervals + i, cacheCau, cache, pointsShifted + i, prec);
            double timeInSparse = (double) (clock() - start2);
            if (timeInRectangular <= timeInSparse)
                cacheCauchy_choiceref(cacheCau) = 0;
            else
                cacheCauchy_choiceref(cacheCau) = 1;
        }
        
//         compApp_poly_ptr app    = cacheApp_getApproximation ( cache, prec );
//      
        for (slong i=0; i<nbPoints; i++)
            if ( ( reUseEvals ) && ( (i%quo)==0 ) ) {
                compApp_set( fvals+i, cacheCauchy_fvalsExref(cacheCau) + i/quo );
                compApp_set( fdervals+i, cacheCauchy_fdervalsExref(cacheCau) + i/quo );
            }
            else {
                if (cacheCauchy_choiceref(cacheCau) == 1)
                    cacheCauchy_sparseEval ( fvals + i, fdervals + i, cacheCau, cache, pointsShifted + i, prec);
                else
                    cacheCauchy_rectangularEval ( fvals + i, fdervals + i, cacheCau, cache, pointsShifted + i, prec);
                
            }
    }
    else {
        for (slong i=0; i<nbPoints; i++) {
            if ( ( reUseEvals ) && ( (i%quo)==0 ) ) {
                compApp_set( fvals+i, cacheCauchy_fvalsExref(cacheCau) + i/quo );
                compApp_set( fdervals+i, cacheCauchy_fdervalsExref(cacheCau) + i/quo );
            }
            else {
                (cacheCauchy_evalFastref(cacheCau))( fvals+i, fdervals + i, pointsShifted+i, prec);
            }
        }
    }
    
    if (metadatas_haveToCount(meta)) {
        slong nbEvals = nbPoints - reUseEvals*cacheCauchy_nbEvalExref(cacheCau);
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            if (certified == CAUCHYTEST_CERTIFIED) {
                metadatas_add_time_CauCoED(meta, (double) (clock() - start));
                metadatas_add_CauchyCoEvalsD(meta, depth, nbEvals);
            } else {
                metadatas_add_time_CauCoEP(meta, (double) (clock() - start));
                metadatas_add_CauchyCoEvalsP(meta, depth, nbEvals);
            }
        } else {
            if (certified == CAUCHYTEST_CERTIFIED) {
                metadatas_add_time_CauExED(meta, (double) (clock() - start));
                metadatas_add_CauchyExEvalsD(meta, depth, nbEvals);
            } else {
                metadatas_add_time_CauExEP(meta, (double) (clock() - start));
                metadatas_add_CauchyExEvalsP(meta, depth, nbEvals);
            }
        }
    }
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeFdivs_fromVals(cacheCauchy_t cacheCau,
                                     int certified,
                                     slong prec,
                                     int inCounting,
                                     metadatas_t meta){
    
    int res=1;
    
    slong nbPoints;
    compApp_ptr fvals;
    compApp_ptr fdervals;
    compApp_ptr fdivs;
    
    if (certified == CAUCHYTEST_CERTIFIED) {
        nbPoints      = cacheCauchy_nbEvalCeref(cacheCau);
        fvals         = cacheCauchy_fvalsCeref(cacheCau);
        fdervals      = cacheCauchy_fdervalsCeref(cacheCau);
        fdivs         = cacheCauchy_fdivsCeref(cacheCau);
    } else {
        nbPoints      = cacheCauchy_nbEvalExref(cacheCau);
        fvals         = cacheCauchy_fvalsExref(cacheCau);
        fdervals      = cacheCauchy_fdervalsExref(cacheCau);
        fdivs         = cacheCauchy_fdivsExref(cacheCau);
    }
    
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
                                        cacheCauchy_t cacheCau,
                                        int certified,
                                        slong rotation,
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta){
    
    int res=1;
    
    slong nbPoints;
    realApp_ptr wP = cacheCauchy_precfdivref(cacheCau);
    compApp_ptr points;
    compApp_ptr fdivs;
    
    if (certified == CAUCHYTEST_CERTIFIED) {
        nbPoints      = cacheCauchy_nbEvalCeref(cacheCau);
//         wP            = cacheCauchy_wanErrCeref(cacheCau);
        points        = cacheCauchy_pointsCeref(cacheCau);
        fdivs         = cacheCauchy_fdivsCeref(cacheCau);
    } else {
        nbPoints      = cacheCauchy_nbEvalExref(cacheCau);
//         wP            = cacheCauchy_wanErrExref(cacheCau) + 0;
        points        = cacheCauchy_pointsExref(cacheCau);
        fdivs         = cacheCauchy_fdivsExref(cacheCau);
    }
    
//     printf("cauchyTest_computeS0Approx_fromVals, rotation: %ld\n", rotation);
//     printf("cauchyTest_computeS0Approx_fromVals, (nbPoints-rotation) modulo nbPoints: %ld\n", (nbPoints- rotation)%nbPoints);
//     if (res==1){
        
        clock_t start = clock();
        
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
        res = res && (realApp_lt( radRe, wP )) 
                  && (realApp_lt( radIm, wP ));
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
        
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }
    
    return res;
}

void cauchyTest_computeShApprox_fromVals(compApp_t sh,
                                         slong h,
                                         cacheCauchy_t cacheCau,
                                         slong prec,
                                         int inCounting,
                                         metadatas_t meta){
    
    slong nbPoints = cacheCauchy_nbEvalCeref(cacheCau);
    compApp_ptr points = cacheCauchy_pointsCeref(cacheCau);
    compApp_ptr fdivs  = cacheCauchy_fdivsCeref(cacheCau);
    
    clock_t start = clock();
    
    /* compute h-th powerSum */
    compApp_mul(sh, fdivs + 0, points + 0, prec);
    for (slong i = 1; i<nbPoints; i++)
        compApp_addmul(sh , fdivs + i, points + ((h+1)*i)%nbPoints, prec);
    compApp_div_si(sh, sh, nbPoints, prec);
    /* no need to scale by the radius because points are already */
    
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }
    
}

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeSsApprox_fromVals(compApp_ptr ps,
                                        cacheCauchy_t cacheCau,
//                                         int certified, /* should be 0 */
//                                         slong rotation,/* should be 1 */
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta){
    
    int res=1;
    
    slong nbPoints;
    slong nbPowerSums = cacheCauchy_nbPwSuExref(cacheCau);
    realApp_ptr wP = cacheCauchy_precfdivref(cacheCau);
    compApp_ptr points;
    compApp_ptr fdivs;
    
//     if (certified == CAUCHYTEST_CERTIFIED) {
//         nbPoints      = cacheCauchy_nbEvalCeref(cacheCau);
// //         wP            = cacheCauchy_wanErrCeref(cacheCau);
//         points        = cacheCauchy_pointsCeref(cacheCau);
//         fdivs         = cacheCauchy_fdivsCeref(cacheCau);
//     } else {
        nbPoints      = cacheCauchy_nbEvalExref(cacheCau);
//         wP            = cacheCauchy_wanErrExref(cacheCau) + 0;
        points        = cacheCauchy_pointsExref(cacheCau);
        fdivs         = cacheCauchy_fdivsExref(cacheCau);
//     }
    
//     printf("cauchyTest_computeS0Approx_fromVals, rotation: %ld\n", rotation);
//     printf("cauchyTest_computeS0Approx_fromVals, (nbPoints-rotation) modulo nbPoints: %ld\n", (nbPoints- rotation)%nbPoints);
//     if (res==1){
        
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
        for (slong j = 0; j<nbPowerSums; j++)
            compApp_div_si(ps+j, ps+j, nbPoints, prec);
        /* no need to scale by the radius because points are already */
//         compApp_mul_realRat_in_place(ps, radius, CCLUSTER_DEFAULT_PREC);
        
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
        
        realApp_clear(radRe);
        realApp_clear(radIm);
        
//     }
    if (metadatas_haveToCount(meta)) {
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            metadatas_add_time_CauCoCS(meta, (double) (clock() - start));
        } else {
            metadatas_add_time_CauExCS(meta, (double) (clock() - start));
        }
    }
    
    return res;
}

void cauchyTest_shiftFFT(                 const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,           
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          int certified,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth ) {
    
    // //             printf("here\n");
    
                clock_t start = clock();
                
                slong nbPoints = cacheCauchy_nbEvalCeref(cacheCau);
                
                compApp_poly_t shiftedPoly;
                compApp_poly_init2( shiftedPoly, nbPoints );
                /* compute approximation of the pol */
                compApp_poly_set(shiftedPoly, cacheApp_getApproximation ( cache, prec ) );
                
                compApp_t c;
                compApp_init(c);
                realRat_t argu;
                realRat_init(argu);
                compApp_t a;
                compApp_init(a);
                /* compute approximation of the center */
                if (radius2==NULL)
                    compApp_set_compRat(c, center, 2*prec);
                else {
                    compApp_set_compRat         (c,    center,   2*prec);
                    realRat_set_si              (argu, 2*vindex, vangle);
                    compApp_set_realRat         (a,    argu,     2*prec);
                    acb_exp_pi_i                (a,    a,        2*prec);
                    compApp_mul_realRat_in_place(a,    radius2,  2*prec);
                    compApp_add                 (c,    c, a,     2*prec);
                }
                /* shift in the center */
                clock_t start2 = clock();
                _acb_poly_taylor_shift_convolution(shiftedPoly->coeffs, c, shiftedPoly->length, prec);
                compApp_poly_scale_realRat_in_place( shiftedPoly->coeffs, radius, shiftedPoly->length, prec );
                if (metadatas_haveToCount(meta))
                    metadatas_add_time_Taylors(meta, (double) (clock() - start2) );
                
                /* compute fft */
                acb_dft_pre_t t;
                acb_dft_precomp_init(t, nbPoints, prec);
                acb_dft_precomp(cacheCauchy_fvalsCeref(cacheCau), shiftedPoly->coeffs, t, prec);
                /* derivative */
                compApp_poly_derivative( shiftedPoly, shiftedPoly, prec);
                slong index = 0;
                while (index < shiftedPoly->length){
                    compApp_div_realRat( (shiftedPoly->coeffs) + index, (shiftedPoly->coeffs) + index, radius, prec );
                    index++;
                }
                compApp_zero((shiftedPoly->coeffs) + (shiftedPoly->length));
                acb_dft_precomp(cacheCauchy_fdervalsCeref(cacheCau), shiftedPoly->coeffs, t, prec);
                
                compApp_poly_clear(shiftedPoly);
                compApp_clear(c);
                realRat_clear(argu);
                compApp_clear(a);
                acb_dft_precomp_clear(t);
                
    if (metadatas_haveToCount(meta)) {
        slong nbEvals = nbPoints;
        if (inCounting == CAUCHYTEST_INCOUNTIN) {
            if (certified == CAUCHYTEST_CERTIFIED) {
                metadatas_add_time_CauCoED(meta, (double) (clock() - start));
                metadatas_add_CauchyCoEvalsD(meta, depth, nbEvals);
            } else {
                metadatas_add_time_CauCoEP(meta, (double) (clock() - start));
                metadatas_add_CauchyCoEvalsP(meta, depth, nbEvals);
            }
        } else {
            if (certified == CAUCHYTEST_CERTIFIED) {
                metadatas_add_time_CauExED(meta, (double) (clock() - start));
                metadatas_add_CauchyExEvalsD(meta, depth, nbEvals);
            } else {
                metadatas_add_time_CauExEP(meta, (double) (clock() - start));
                metadatas_add_CauchyExEvalsP(meta, depth, nbEvals);
            }
        }
    }
    
}

cauchyTest_res cauchyTest_computeS0Approx(compApp_t ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,           
                                          slong rotation,
                                          int *alreadyEvaluated,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          int certified,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth ){
    
//     slong nbPoints = (certified == CAUCHYTEST_CERTIFIED ? cacheCauchy_nbEvalCeref(cacheCau) : cacheCauchy_nbEvalExref(cacheCau) );
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
    
    if (*alreadyEvaluated == 0) {
        
        while (res.nbOfSol==-1){
            cauchyTest_getEvaluationPoints( center, radius,
                                                radius2, vangle, vindex,
                                                cacheCau, certified, res.appPrec);
                
            cauchyTest_evaluateAtPoints( cache, cacheCau, certified, res.appPrec, inCounting, meta, depth);
        
            res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, certified, res.appPrec, inCounting, meta);
            if ( res.nbOfSol ==-1 )
                res.appPrec = 2*res.appPrec;
        }
        *alreadyEvaluated = 1;
    }
    
    if (res.nbOfSol == -2)
        return res;
    
    /* compute approximation of Power sums */
    res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, cacheCau, certified, rotation, res.appPrec, inCounting, meta);
    
    while ( res.nbOfSol ==-1 ) {
        res.appPrec = 2*res.appPrec;
        cauchyTest_getEvaluationPoints( center, radius, 
                                            radius2, vangle, vindex,
                                            cacheCau, certified, res.appPrec);
            
        cauchyTest_evaluateAtPoints( cache, cacheCau, certified, res.appPrec, inCounting, meta, depth);
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, certified, res.appPrec, inCounting, meta);
         
        /* compute approximation of Power sums */
        res.nbOfSol = cauchyTest_computeS0Approx_fromVals(ps, cacheCau, certified, rotation, res.appPrec, inCounting, meta);
    }
    
    return res;
}

cauchyTest_res cauchyTest_computeSsApprox(compApp_ptr ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,           
//                                           slong rotation,
//                                           int *alreadyEvaluated,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          int certified, /* should be CAUCHYTEST_UNCERTIFI */
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth ){
    
//     slong nbPoints = cacheCauchy_nbEvalExref(cacheCau);
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
        
    while (res.nbOfSol==-1){
        /* compute points and evals at prec res.appPrec*/
        cauchyTest_getEvaluationPoints( center, radius,
                                    radius2, vangle, vindex,
                                    cacheCau, certified, res.appPrec);
        
        cauchyTest_evaluateAtPoints( cache, cacheCau, certified, res.appPrec, inCounting, meta, depth);
        
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, certified, res.appPrec, inCounting, meta);
        if ( res.nbOfSol ==-1 )
            res.appPrec = 2*res.appPrec;
    }
    
    if (res.nbOfSol == -2)
        return res;
    
    /* compute approximation of Power sums */
    res.nbOfSol = cauchyTest_computeSsApprox_fromVals(ps, cacheCau, res.appPrec, inCounting, meta);
    
    while ( res.nbOfSol ==-1 ) {
        res.appPrec = 2*res.appPrec;
        
        cauchyTest_getEvaluationPoints( center, radius, 
                                        radius2, vangle, vindex,
                                        cacheCau, certified, res.appPrec);
        
        cauchyTest_evaluateAtPoints( cache, cacheCau, certified, res.appPrec, inCounting, meta, depth);
        /* compute approximation of Power sums */
        res.nbOfSol = cauchyTest_computeFdivs_fromVals(cacheCau, certified, res.appPrec, inCounting, meta);
        res.nbOfSol = cauchyTest_computeSsApprox_fromVals(ps, cacheCau, res.appPrec, inCounting, meta);
    }
    
    cacheCauchy_precEvalExref(cacheCau) = res.appPrec;
    
    return res;
    
}

slong cauchyTest_computeS0compDsk( const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   cacheApp_t cache,
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
    
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
    compApp_poly_ptr app = NULL;
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=3)
            printf("#------ appPrec: %ld\n", appPrec);
        
        enoughPrec = 1;
        compApp_zero(s0star);
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        if (cacheCauchy_evalFastref(cacheCau) == NULL)
            app = cacheApp_getApproximation ( cache, appPrec );
        
        for(slong i=0; i<q; i++) {
            realRat_set_si(argu, 2*i, q);
            compApp_set_realRat(a, argu, appPrec);
            acb_exp_pi_i( point, a, appPrec);
            compApp_mul_realRat(pointShifted, point, compDsk_radiusref(Delta), appPrec);
            compApp_add( pointShifted, c, pointShifted, appPrec);
            
            start = clock();
            if (cacheCauchy_evalFastref(cacheCau) == NULL)
                compApp_poly_evaluate2_rectangular(fval, fderval, app, pointShifted, appPrec);
            else
                (cacheCauchy_evalFastref(cacheCau))( fval, fderval, pointShifted, appPrec);
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
    compApp_clear(a);
    realRat_clear(argu);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoED(meta, evalTime);
        metadatas_add_CauchyCoEvalsD(meta, depth, q);
    }
    
    return res;
}

slong cauchyTest_computeS1compDsk( compDsk_t res,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   metadatas_t meta, int depth) {
    
    int level = 5;
    slong appPrec = CCLUSTER_DEFAULT_PREC;
    
    realRat_t error;
    realRat_init(error);
//     realRat_set_si(error, 2, 1);
//     realRat_pow_si(error, error, -prec-1); /* error = (1/2)*2^{-prec} */
    realRat_div_ui(error, eps, 2);
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_realRat(errAp, error, CCLUSTER_DEFAULT_PREC);    /* errAp = contains error */
    
    realApp_t qApp, logIsoRatio;
    realApp_init(qApp);
    realApp_init(logIsoRatio);
    realApp_set_realRat(logIsoRatio, isoRatio,                                               CCLUSTER_DEFAULT_PREC);
    realApp_log        (logIsoRatio, logIsoRatio,                                            CCLUSTER_DEFAULT_PREC);
    realApp_set_realRat(qApp,        compDsk_radiusref(Delta),                               CCLUSTER_DEFAULT_PREC);
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
    
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
    compApp_poly_ptr app = NULL;
    
    int enoughPrec = 0;
    while (enoughPrec == 0) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: precision: %ld\n",  appPrec);
        }
        
        enoughPrec = 1;
        compApp_zero(s1);
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        if (cacheCauchy_evalFastref(cacheCau) == NULL)
            app = cacheApp_getApproximation ( cache, appPrec );
        
        for(slong i=0; i<q; i++) {
            realRat_set_si(argu, 2*i, q);
            compApp_set_realRat(a, argu, appPrec);
            acb_exp_pi_i( point, a, appPrec);
            compApp_mul_realRat(pointShifted, point, compDsk_radiusref(Delta), appPrec);
            compApp_add( pointShifted, c, pointShifted, appPrec);
            
            if (cacheCauchy_evalFastref(cacheCau) == NULL)
                compApp_poly_evaluate2_rectangular(fval, fderval, app, pointShifted, appPrec);
            else
                (cacheCauchy_evalFastref(cacheCau))( fval, fderval, pointShifted, appPrec);
            
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
    /* s1 approximates rs1(p_\Delta,0,1) with |s1 - rs1(p_\Delta,0,1)|\geq errAp */
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
    compApp_clear(a);
    realRat_clear(argu);
    compRat_clear(mc);
    
    return appPrec;
}

/* Assume Delta is isoRatio-isolated and contains nbOfRoots roots */
/* Assume SS is initialized to contain at least nbOfRoots + 1 numbers */
/* Computes nbOfRoots + 1 power sums with error less than eps     */
cauchyTest_res cauchyTest_computeSScompDsk( compApp_ptr SS,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   slong prec,
                                   metadatas_t meta, int depth){
    
    int level = 3;
    slong appPrec = CCLUSTER_DEFAULT_PREC;
    
    /* number of power sums to compute */
    slong n = nbOfRoots + 1;
    /* compute the required number of points */
    realRat_t error;
    realRat_init(error);
    realRat_div_ui(error, eps, 2);
    realApp_t errAp;
    realApp_init(errAp);
    realApp_set_realRat(errAp, error, CCLUSTER_DEFAULT_PREC);    /* errAp = contains error */
    
    realApp_ptr errorsSS = (realApp_ptr) ccluster_malloc ( n*sizeof(realApp));
    for (slong i = 0; i < n; i++)
        realApp_init( errorsSS + i );
    
    realApp_t qApp, IsoRatioTotheN, IsoRatioInvTotheN, logIsoRatio;
    realApp_init(qApp);
    realApp_init(IsoRatioTotheN);
    realApp_init(IsoRatioInvTotheN);
    realApp_init(logIsoRatio);
    realApp_set_realRat(logIsoRatio, isoRatio,                                               CCLUSTER_DEFAULT_PREC);
    realApp_log        (logIsoRatio, logIsoRatio,                                            CCLUSTER_DEFAULT_PREC);
    realApp_one(IsoRatioTotheN);
    slong qmax = 0;
    
    for (slong i = 0; i < n; i++) {
        realApp_mul_si     (errorsSS + i, IsoRatioTotheN, cacheCauchy_degreeref(cacheCau)-nbOfRoots, CCLUSTER_DEFAULT_PREC);
        realApp_inv        (IsoRatioInvTotheN, IsoRatioTotheN,                               CCLUSTER_DEFAULT_PREC);
        realApp_mul_si     (IsoRatioInvTotheN, IsoRatioInvTotheN, nbOfRoots,                         CCLUSTER_DEFAULT_PREC);
        realApp_add        (errorsSS + i, errorsSS + i, IsoRatioInvTotheN,                   CCLUSTER_DEFAULT_PREC);
        
        realApp_div_realRat(qApp,         errorsSS + i,             error,                   CCLUSTER_DEFAULT_PREC);
        realApp_add_si     (qApp,         qApp,                     1,                       CCLUSTER_DEFAULT_PREC);
        realApp_log        (qApp,         qApp,                                              CCLUSTER_DEFAULT_PREC);
        realApp_div        (qApp,         qApp,                     logIsoRatio,             CCLUSTER_DEFAULT_PREC);
        slong q = realApp_ceil_si( qApp, CCLUSTER_DEFAULT_PREC ) + 1;
        qmax = CCLUSTER_MAX( qmax, q );
        
        realApp_mul_realRat(IsoRatioTotheN,    IsoRatioTotheN,     isoRatio,                 CCLUSTER_DEFAULT_PREC);
    }
    
//     if (metadatas_getVerbo(meta)>=level) {
//         printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: nb of eval points: %ld\n", qmax);
//     }
    
    qmax = CCLUSTER_MAX( qmax, n );
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: nb of eval points: %ld\n", qmax);
        printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: required max erro: ");
        realApp_printd( errAp, 10);
        printf("\n");
    }
    
    realApp_set_realRat( IsoRatioTotheN, isoRatio,       CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui     ( IsoRatioTotheN, IsoRatioTotheN, qmax, CCLUSTER_DEFAULT_PREC);
    realApp_add_si     ( IsoRatioTotheN, IsoRatioTotheN, -1,   CCLUSTER_DEFAULT_PREC);
    for (slong i = 0; i < n; i++) {
        realApp_div( errorsSS + i, errorsSS + i, IsoRatioTotheN, CCLUSTER_DEFAULT_PREC);
//         if (metadatas_getVerbo(meta)>=level) {
//             printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: %ld-th error:", i);
//             realApp_printd( errorsSS + i, 10);
//             printf("\n");
//         }
    }
    
    realApp_t ubound, lbound;
    realApp_init( ubound );
    realApp_init( lbound );
    realApp_set_realRat( ubound, isoRatio,   CCLUSTER_DEFAULT_PREC);
    realApp_add_si     ( ubound, ubound, -1, CCLUSTER_DEFAULT_PREC);
    realApp_mul_realRat( ubound, ubound, compDsk_radiusref(Delta), CCLUSTER_DEFAULT_PREC);
    realApp_div_realRat( ubound, ubound, isoRatio,                 CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui     ( lbound, ubound, cacheCauchy_degreeref(cacheCau), CCLUSTER_DEFAULT_PREC);
    realApp_inv        ( ubound, ubound,                           CCLUSTER_DEFAULT_PREC);
    realApp_mul_si     ( ubound, ubound, cacheCauchy_degreeref(cacheCau), CCLUSTER_DEFAULT_PREC);
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
    realApp_t modulus, radRe, radIm;
    compApp_t c, a;
    realRat_t argu;
    
    compApp_init( fval         );
    compApp_init( fderval      );
    compApp_init( fdiv         );
    realApp_init( modulus);
    realApp_init(radRe);
    realApp_init(radIm);
    compApp_init(c);
    compApp_init(a);
    realRat_init(argu);
    
    compApp_poly_ptr app = NULL;
    
    int enoughPrec = -1;
    /* compute the points */
    while (enoughPrec == -1) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: precision: %ld\n",  appPrec);
        }
        
        enoughPrec = 1;
        compApp_set_compRat(c, compDsk_centerref(Delta), appPrec);
        
        if (cacheCauchy_evalFastref(cacheCau) == NULL)
            app = cacheApp_getApproximation ( cache, appPrec );
        
        compApp_ptr points = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        compApp_ptr pointsShifted = (compApp_ptr) ccluster_malloc ( qmax*sizeof(compApp) );
        
        for(slong i=0; i<qmax; i++) {
            compApp_init (points + i);
            compApp_init (pointsShifted + i);
            
            realRat_set_si(argu, 2*i, qmax);
            compApp_set_realRat(a, argu, appPrec);
            acb_exp_pi_i( points + i, a, appPrec);
            compApp_mul_realRat(pointsShifted + i, points + i, compDsk_radiusref(Delta), appPrec);
            compApp_add( pointsShifted + i, c, pointsShifted + i, appPrec);
        }
        
        /* initialize power sums */
        for (slong j = 0; j< n; j++)
            compApp_zero( SS + j );
        
        for(slong i=0; (i<qmax) && (enoughPrec==1) ; i++) {
            
            if (cacheCauchy_evalFastref(cacheCau) == NULL)
                compApp_poly_evaluate2_rectangular(fval, fderval, app, pointsShifted + i, appPrec);
            else
                (cacheCauchy_evalFastref(cacheCau))( fval, fderval, pointsShifted + i, appPrec);
            
            /* check if fval has modulus less that lbound */
            compApp_abs(modulus, fval, appPrec);
            if (realApp_lt( modulus, lbound )){
                enoughPrec=-2;
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: fval has modulus less than lbound, %ld\n", i);
                }
            } else if (compApp_contains_zero( fval )){
                enoughPrec=-1;
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: fval contains 0, i: %ld\n", i);
                }
            } else {
                /* compute fdiv */
                compApp_div(fdiv, fderval, fval, appPrec);
                compApp_abs(modulus, fdiv, appPrec);
                /* check if the ratio is less than the upper bound */
                if (realApp_gt( modulus, ubound )) {
                    enoughPrec = -2;
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------------cauchy_tests.c: cauchyTest_computeS1compDsk: ratio has modulus more than ubound, %ld\n", i);
                    }
                } else {
                    for (slong j = 0; j< n; j++) 
                        compApp_addmul(SS+j , fdiv, points + ((j+1)*i)%qmax, appPrec); 
                }
            }
            
        }
        
        if (enoughPrec==1) {
            for (slong j = 0; j< n; j++) {
                compApp_div_si(SS+j, SS+j, qmax, appPrec);
                compApp_mul_realRat(SS+j, SS+j, compDsk_radiusref(Delta), appPrec);
            }
            /* no need to scale by the radius because points are already */
            
            /* check if precisions are ok */
            for (slong j = 0; j< n; j++) {
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
    
    
    for (slong i = 0; i < n; i++)
        realApp_clear( errorsSS + i );
    ccluster_free(errorsSS);
    
    realRat_clear(error);
    realApp_clear(errAp);
    realApp_clear(qApp);
    realApp_clear(IsoRatioTotheN);
    realApp_clear(IsoRatioInvTotheN);
    realApp_clear(logIsoRatio);
    
    realApp_clear( ubound );
    realApp_clear( lbound );
    
    compApp_clear( fval         );
    compApp_clear( fderval      );
    compApp_clear( fdiv         );
    realApp_clear( modulus);
    realApp_clear(radRe);
    realApp_clear(radIm);
    compApp_clear(c);
    compApp_clear(a);
    realRat_clear(argu);
    
    cauchyTest_res res;
    res.appPrec = appPrec;
    res.nbOfSol = enoughPrec;
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------cauchy_tests.c: cauchyTest_computeSScompDsk: power sums:\n");
        for (slong j=0; j<n; j++) {
            printf("#--------------- %ld-th: ", j); compApp_printd( SS+j, 10); printf("\n");
        }
    }
    
    return res;
}
