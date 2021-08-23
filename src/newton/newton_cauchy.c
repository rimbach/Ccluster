/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "newton/newton.h"

newton_res newton_cauchy_first_condition( compApp_t fcenter, 
                                          compApp_t fpcenter, 
                                          cacheApp_t cache, 
                                          cacheCauchy_t cacheCau, 
                                          const compRat_t c, 
                                          const realRat_t d, 
                                          slong prec, 
                                          metadatas_t meta){
    
    int soft_comp_res = 0;
    newton_res res;
    if (metadatas_usePredictPrec(meta))
        res.appPrec = prec;
    else
        res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    compApp_t center;
    compApp_init(center);
    
    realApp_t diam, fcenterabs, fpcenterabs;
    realApp_init(diam);
    realApp_init(fcenterabs);
    realApp_init(fpcenterabs);
    
    compApp_poly_t pApprox;
    compApp_poly_init(pApprox);
    
    compApp_set_compRat( center, c, res.appPrec );
    realApp_set_realRat( diam, d, res.appPrec);
    if (cacheCauchy_evalFastref(cacheCau) == NULL) {
        tstar_getApproximation( pApprox, cache, res.appPrec, meta);
        compApp_poly_evaluate2(fcenter, fpcenter, pApprox, center, res.appPrec);
    } else {
//         printf("\n\n  prec: %ld\n", res.appPrec);
//         printf("     point: "); compApp_printd(center, 20); printf("\n");
        (cacheCauchy_evalFastref(cacheCau))(fcenter, fpcenter, center, res.appPrec);
//         printf("  f(point): "); compApp_printd(fcenter, 20); printf("\n");
//         printf(" f'(point): "); compApp_printd(fpcenter, 20); printf("\n");
    }
    compApp_abs(fcenterabs, fcenter, res.appPrec);
    compApp_abs(fpcenterabs, fpcenter, res.appPrec);
    realApp_mul(fpcenterabs, fpcenterabs, diam, res.appPrec);
    soft_comp_res = realApp_soft_compare(fpcenterabs, fcenterabs, res.appPrec);

    while( soft_comp_res == -2 ){
        res.appPrec *=2;
        compApp_set_compRat( center, c, res.appPrec );
        realApp_set_realRat( diam, d, res.appPrec);
        if (cacheCauchy_evalFastref(cacheCau) == NULL) {
            tstar_getApproximation( pApprox, cache, res.appPrec, meta);
            compApp_poly_evaluate2(fcenter, fpcenter, pApprox, center, res.appPrec);
        } else {
            (cacheCauchy_evalFastref(cacheCau))(fcenter, fpcenter, center, res.appPrec);
        }
        compApp_abs(fcenterabs, fcenter, res.appPrec);
        compApp_abs(fpcenterabs, fpcenter, res.appPrec);
        realApp_mul(fpcenterabs, fpcenterabs, diam, res.appPrec);
        soft_comp_res = realApp_soft_compare(fpcenterabs, fcenterabs, res.appPrec);  
    }
    
    compApp_clear(center);
    realApp_clear(diam);
    realApp_clear(fcenterabs);
    realApp_clear(fpcenterabs);
    compApp_poly_clear(pApprox);
    
    if (soft_comp_res==0)
        res.nflag = 0;
    else
        res.nflag = 1;
    return res;
}

newton_res newton_cauchy_iteration( compApp_t iteration, 
                                    cacheApp_t cache, 
                                    cacheCauchy_t cacheCau,
                                    const connCmp_t CC, 
                                    const compRat_t c, 
                                    compApp_t fcenter, 
                                    compApp_t fpcenter,
                                    slong prec, metadatas_t meta){
    
    newton_res res;
    res.nflag=1;
    if (metadatas_usePredictPrec(meta))
        res.appPrec = prec;
    else
        res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    realRat_t errorBound;
    compApp_t center;
    realApp_t iterationError, errorBoundApp;
    compApp_poly_t pApprox;
    realRat_init(errorBound);
    compApp_init(center);
    realApp_init(iterationError);
    realApp_init(errorBoundApp);
    compApp_poly_init(pApprox);
    
    realRat_set_si(errorBound, 1, 64);
    realRat_div_fmpz(errorBound, errorBound, connCmp_nwSpdref(CC));
    realRat_mul(errorBound, errorBound, connCmp_widthref(CC));
    
    compApp_set_compRat( center, c, res.appPrec );
    compApp_div(iteration, fcenter, fpcenter, res.appPrec);
    compApp_mul_si(iteration, iteration, connCmp_nSolsref(CC), res.appPrec);
    compApp_sub(iteration, center, iteration, res.appPrec);
    compApp_abs(iterationError, iteration, res.appPrec);
    realApp_get_rad_realApp(iterationError, iterationError);
    realApp_set_realRat(errorBoundApp, errorBound, res.appPrec);
    
    while ( (realApp_is_finite(iterationError)==0)||(realApp_ge(iterationError, errorBoundApp)==1)) {
        res.appPrec *=2;
        compApp_set_compRat( center, c, res.appPrec );
        if (cacheCauchy_evalFastref(cacheCau) == NULL) {
            tstar_getApproximation( pApprox, cache, res.appPrec, meta);
            compApp_poly_evaluate2(fcenter, fpcenter, pApprox, center, res.appPrec);
        } else {
            (cacheCauchy_evalFastref(cacheCau))(fcenter, fpcenter, center, res.appPrec);
        }
        compApp_div(iteration, fcenter, fpcenter, res.appPrec);
        compApp_mul_si(iteration, iteration, connCmp_nSolsref(CC), res.appPrec);
        compApp_sub(iteration, center, iteration, res.appPrec);
        compApp_abs(iterationError, iteration, res.appPrec);
        realApp_get_rad_realApp(iterationError, iterationError);
        realApp_set_realRat(errorBoundApp, errorBound, res.appPrec);
    }
    
    compApp_clear(center);
    realRat_clear(errorBound);
    realApp_clear(iterationError);
    realApp_clear(errorBoundApp);
    compApp_poly_clear(pApprox);
    return res;
}

/* interval newton */
/* performs a newton test for the compBox b contained in the compDisk d;
 * returns 1 if interval newton certifies the existence of a solution in b;
 *         0 otherwise */
newton_res newton_cauchy_interval(  compDsk_t d, 
                                    cacheApp_t cache, 
                                    cacheCauchy_t cacheCau,
                                    slong prec, 
                                    metadatas_t meta) {
    
    newton_res res;
    res.nflag=0;
    if (metadatas_usePredictPrec(meta))
        res.appPrec = prec;
    else
        res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    compApp_t cBall, ball, fcBall, fpBall, dummy;
    realApp_t bRe, bIm, error;
    realRat_t nwidth;
    compApp_poly_t pApprox, ppApprox;
    
    compApp_init(cBall);
    compApp_init(ball);
    compApp_init(fcBall);
    compApp_init(fpBall);
    compApp_init(dummy);
    realApp_init(bRe);
    realApp_init(bIm);
    realApp_init(error);
    realRat_init(nwidth);
    compApp_poly_init(pApprox);
    compApp_poly_init(ppApprox);
    
    realApp_set_realRat( bRe, compRat_realref( compDsk_centerref(d) ), res.appPrec );
    realApp_set_realRat( bIm, compRat_imagref( compDsk_centerref(d) ), res.appPrec );
    realRat_set_si( nwidth, 2,3);
    realRat_mul( nwidth, nwidth, compDsk_radiusref(d));
    realApp_set_realRat( error, nwidth, res.appPrec );
    arb_add_error(bRe, error);
    arb_add_error(bIm, error);
    compApp_set_real_realApp(ball, bRe);
    compApp_set_imag_realApp(ball, bIm);
    compApp_set_compRat(cBall, compDsk_centerref(d), res.appPrec );
    
    if (cacheCauchy_evalFastref(cacheCau) == NULL) {
        tstar_getApproximation( pApprox, cache, res.appPrec, meta);
        compApp_poly_derivative(ppApprox, pApprox, res.appPrec);
        compApp_poly_evaluate(fpBall, ppApprox, ball, res.appPrec);
    } else {
        (cacheCauchy_evalFastref(cacheCau))(dummy, fpBall, ball, res.appPrec);
    }
    
    if (compApp_contains_zero(fpBall)){
        res.nflag=0; /* do nothing */
//         printf("fpBall contains 0: %d, fpBall: ", compApp_contains_zero(fpBall)); compApp_print(fpBall); printf("\n");
    }
    else {
        if (cacheCauchy_evalFastref(cacheCau) == NULL) {
            compApp_poly_evaluate(fcBall, pApprox, cBall, res.appPrec);
        } else {
            (cacheCauchy_evalFastref(cacheCau))(fcBall, dummy, cBall, res.appPrec);
        }
        compApp_div( fcBall, fcBall, fpBall, res.appPrec);
        compApp_sub( fcBall, cBall, fcBall, res.appPrec);
        if (compApp_contains(ball, fcBall)) {
            res.nflag=1;
        }
    }
    
//     printf("result of interval newton: "); compApp_print(fcBall); printf("\n");
//     printf("initial box              : "); compApp_print(ball); printf("\n");
//     printf("interval newton result: %d \n", res.nflag);
    
    
    compApp_clear(cBall);
    compApp_clear(ball);
    compApp_clear(fcBall);
    compApp_clear(fpBall);
    compApp_clear(dummy);
    realApp_clear(bRe);
    realApp_clear(bIm);
    realApp_clear(error);
    realRat_clear(nwidth);
    compApp_poly_clear(pApprox);
    compApp_poly_clear(ppApprox);
    
    return res;
}

newton_res newton_cauchy_newton_connCmp( connCmp_t nCC,
                                         connCmp_t CC,
                                         cacheApp_t cache,
                                         cacheCauchy_t cacheCau,
                                         const compRat_t c,
                                         slong prec, 
                                         metadatas_t meta) {
    
    int level = 4;
//     printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: begin --------------------------- \n");
    newton_res res;
    if (metadatas_usePredictPrec(meta))
        res.appPrec = prec;
    else
        res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    realRat_t fourcc, two, nwidth;
    compDsk_t ndisk;
    compApp_t fcenter, fpcenter, iteration;
    compBox_list_ptr ltemp;
    compBox_ptr btemp;
    realRat_init(fourcc);
    realRat_init(two);
    realRat_init(nwidth);
    compDsk_init(ndisk);
    compApp_init(fcenter);
    compApp_init(fpcenter);
    compApp_init(iteration);
    
    realRat_set_si(two,2,1);
    connCmp_diameter(fourcc, CC);
    realRat_mul(fourcc, fourcc, two);
    
//     slong precsave = prec;
    
    res = newton_cauchy_first_condition( fcenter, fpcenter, cache, cacheCau, c, fourcc, res.appPrec, meta);
//     printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: firstCondOk (nflag: %d) prec avant: %ld, prec apres: %ld \n", res.nflag, prec, res.appPrec);
    
    if (res.nflag){
//         precsave = res.appPrec;
        res = newton_cauchy_iteration( iteration, cache, cacheCau, CC, c, fcenter, fpcenter, res.appPrec, meta);
        compApp_get_compRat( compDsk_centerref(ndisk), iteration);
        realRat_set_si(compDsk_radiusref(ndisk),1,8);
        realRat_div_fmpz(compDsk_radiusref(ndisk), compDsk_radiusref(ndisk), connCmp_nwSpdref(CC));
        realRat_mul(compDsk_radiusref(ndisk), compDsk_radiusref(ndisk), connCmp_widthref(CC));
        res.nflag = connCmp_intersection_has_non_empty_interior_compDsk(CC, ndisk);
//         res.appPrec = precsave;
    }
    
//     printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: iterationOK (nflag: %d) prec avant: %d, prec apres: %d \n", res.nflag, prec, res.appPrec);
    
    if (res.nflag) {
//         precsave = res.appPrec;
//         if (res.appPrec == prec) res.appPrec *=2;
        
        /*improvement: if one solution, try to validate with interval newton*/
        /* to avoid a costly taylor-shift in the tstar tes                  */ 
        newton_res nres;
        nres.nflag = 0;
        if (connCmp_nSolsref(CC)==1) {
            if (metadatas_getVerbo(meta)>=level)
                printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: the CC contains one solution: try to validate with interval newton\n");
            nres = newton_cauchy_interval( ndisk, cache, cacheCau, res.appPrec, meta);
            if (metadatas_getVerbo(meta)>level) {
                printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: nres.nflag: %d\n", nres.nflag);
//                 newton_res nrest;
//                 nrest = newton_interval( ndisk, cache, res.appPrec, meta);
//                 printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: nrest.nflag: %d\n", nrest.nflag);
            }
        }
        
        if (nres.nflag==0) {
//             printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: failed... use cauchy deterministic counting test\n");
            slong depth = connCmp_getDepth(CC, metadatas_initBref(meta));
            
            res.nflag = 0;
            cauchyTest_res cres;
            if ( connCmp_nSolsref(CC) <= 3 ) {
                /* try the combinatorial version of Cauchy counting test */
                realRat_mul_si(compDsk_radiusref(ndisk), compDsk_radiusref(ndisk), 2);
                cres = cauchyTest_deterministic_counting_combinatorial( compDsk_centerref(ndisk),
                                                                        compDsk_radiusref(ndisk),
                                                                        connCmp_nSolsref(CC),
                                                                        cache, cacheCau, res.appPrec,
                                                                        meta, depth);
                res.appPrec = cres.appPrec;
                res.nflag = (cres.nbOfSol == connCmp_nSolsref(CC));
                if (! res.nflag )
                    realRat_div_ui(compDsk_radiusref(ndisk), compDsk_radiusref(ndisk), 2);
                if (metadatas_getVerbo(meta)>level) {
                    printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: cauchy test comb res %d\n", cres.nbOfSol);
                }
            } else {
//             if (res.nflag==0){
                /* version with exclusion of discs on contour */
                cres = cauchyTest_deterministic_verification( ndisk, connCmp_nSolsref(CC), cacheCauchy_isoRatioref(cacheCau),
                                                               cache, cacheCau, res.appPrec, meta, depth);
                res.appPrec = cres.appPrec;
//                 cres = cauchyTest_probabilistic_counting( ndisk, cache, cacheCau, res.appPrec, meta, depth);
//                 res.appPrec = cres.appPrec;
//                 if (cres.nbOfSol == connCmp_nSolsref(CC)) {
//                     cres = cauchyTest_deterministic_counting( ndisk, cacheCauchy_isoRatioref(cacheCau), cache, cacheCau, res.appPrec, meta, depth);
//                     res.appPrec = cres.appPrec;
//                 }
                res.nflag = (cres.nbOfSol == connCmp_nSolsref(CC));
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: cauchy test res %d\n", cres.nbOfSol);
                    /* version with Tstar test: only for debug*/
                    tstar_res tres = tstar_interface( cache, ndisk, connCmp_nSolsref(CC), 0, 1, res.appPrec, depth, NULL, meta);
                    printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: tstar  test res %d\n\n", tres.nbOfSol);
                }
            }
        }
//         printf("number of sols in ndisk from pellet test: %d\n\n", tres.nbOfSol);
        /* end improvement */
        
    }
    
//     printf("Newton: tstarOk     (nflag: %d) prec avant: %d, prec apres: %d \n", res.nflag, prec, res.appPrec);
    if (res.nflag) {
        realRat_set_si(nwidth, 1,2);
        realRat_div_fmpz(nwidth, nwidth, connCmp_nwSpdref(CC));
        realRat_mul(nwidth, nwidth, connCmp_widthref(CC));
        
        ltemp = connCmp_boxesref(CC);
        /*old version with successive quadrisections and intersections with ndisk*/
        /*while (realRat_cmp( compBox_bwidthref(compBox_list_first(ltemp)), nwidth)>0){ */
        /*    btemp = compBox_list_pop(ltemp);                                          */
        /*    subdBox_quadrisect_intersect_compDsk(ltemp, btemp, ndisk);                */
        /*    compBox_clear(btemp);                                                     */
        /*    ccluster_free(btemp); */                                 /*comment it for julia... */
        /*}                                                                             */
        /*new version: we directly compute the boxes of width nwidth intersectiong ndisk*/
        compBox_list_t ltemp2;
        compBox_list_init(ltemp2);
        
        while (compBox_list_get_size(ltemp)>0){
            btemp = compBox_list_pop(ltemp);
            subdBox_quadrisect_with_compDsk( ltemp2, btemp, ndisk, nwidth);
            compBox_clear(btemp);
            ccluster_free(btemp); /*comment it for julia...*/
        }
        compBox_list_swap(ltemp, ltemp2);
        compBox_list_clear(ltemp2);
        
        /*printf("length of list: %d\n", compBox_list_get_size(ltemp));*/
        /*compBox_list_print(ltemp); printf("\n");                     */
        
        
        /* printf("Newton: bisectOk (nflag: %d) (nbboxes: %d)--------------------------- \n", res.nflag, compBox_list_get_size(ltemp));*/
        
        btemp = compBox_list_pop(ltemp);
        realRat_set(connCmp_widthref(nCC), compBox_bwidthref(btemp));
        connCmp_insert_compBox(nCC, btemp);
        while (!compBox_list_is_empty(ltemp))
            connCmp_insert_compBox(nCC, compBox_list_pop(ltemp));
        connCmp_nSols(nCC) = connCmp_nSols(CC);
        fmpz_set(connCmp_nwSpdref(nCC), connCmp_nwSpdref(CC));
        /* test */
        connCmp_isSep(nCC) = connCmp_isSep(CC);
        /* end test */
        /*connCmp_appPrref(nCC) = res.appPrec;*/ /*adjust the precision in the main loop*/
        
    }
    /* printf("Newton: new CC Ok (nflag: %d)--------------------------- \n", res.nflag);*/
    
    realRat_clear(fourcc);
    realRat_clear(two);
    realRat_clear(nwidth);
    compDsk_clear(ndisk);
    compApp_clear(fcenter);
    compApp_clear(fpcenter);
    compApp_clear(iteration);
    
//     chronos_toc_Newtons(metadatas_chronref(meta));
    
//     printf("#newton_cauchy.c, newton_cauchy_newton_connCmp: end --------------------------- \n");
    
    return res;
}

