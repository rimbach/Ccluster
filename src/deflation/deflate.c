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

#include "deflation/deflate.h"

void deflate_taylor_shift_interval_inplace( realApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        realApp_poly_taylorShift_interval_in_place( res, compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d), prec );

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefTayl(meta, (double) (clock() - start) );
}

void deflate_real_graeffe_iterations_inplace( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        for(int i = 0; i < N; i++)
            realApp_poly_oneGraeffeIteration_in_place( res, prec );
        
        if (metadatas_haveToCount(meta))
            metadatas_add_time_Graeffe(meta, (double) (clock() - start) );
}

void deflate_real_graeffe_iterations_abs_two_first_coeffs( realApp_t coeff0, realApp_t coeff1, const realApp_poly_t pApprox, int N, slong prec, metadatas_t meta){
    realApp_poly_t p1, p2;
    realApp_poly_init2(p1, realApp_poly_degree(pApprox)+1);
    realApp_poly_init2(p2, realApp_poly_degree(pApprox)+1);
    realApp_poly_set(p1, pApprox);
    slong bound = 0x1 << N; /* assume 2^N fits in an slong: N<32 for 32bits machine... */
    for( int i = 0; i < N; i++) {
        bound = bound >> 1;
//         printf("bound: %d\n", (int) bound);
        realApp_poly_oneGraeffeIteration_first_n_coeff( p2, p1, CCLUSTER_MIN(realApp_poly_degree(pApprox), bound), realApp_poly_degree(pApprox)+1, prec);
        realApp_poly_swap(p1,p2);
    }
    
    realApp_abs( coeff0, realApp_poly_getCoeff(p1, 0));
    realApp_abs( coeff1, realApp_poly_getCoeff(p1, 1));
    
    realApp_poly_clear(p1);
    realApp_poly_clear(p2);
}

void deflate_derivative_inplace( realApp_poly_t res, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        realApp_poly_derivative(res, res, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefDeri(meta, (double) (clock() - start) );
}

void deflate_evaluate(realApp_t y, const realApp_poly_t f, const realApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
//         arb_poly_evaluate_horner(y, f, x, prec);
        realApp_poly_evaluate(y, f, x, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefEval(meta, (double) (clock() - start) );
}

void realApp_poly_oneGraeffeIteration_lastTerms(realApp_poly_t ls, const realApp_poly_t f, slong split, slong prec, metadatas_t meta ){
    
    clock_t start = clock();
    
    realApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    realApp_poly_t fe, fo;
    realApp_poly_init2(fe, len2);
    realApp_poly_init2(fo, len2);
    realApp_ptr feptr = fe->coeffs;
    realApp_ptr foptr = fo->coeffs;
    
    for (i = split; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) 
            realApp_set( feptr + quo, fptr+i);
        else
            realApp_set( foptr + quo, fptr+i);
    }
    realApp_poly_set_length(fe, len2);
    realApp_poly_set_length(fo, len2);
    
    realApp_poly_t fes, fos;
    realApp_poly_init2(fes, len1);
    realApp_poly_init2(fos, len1);
    realApp_poly_mullow( fes, fe, fe, len1, prec);
    realApp_poly_mullow( fos, fo, fo, len1, prec);
    realApp_poly_shift_left( fos, fos, 1 );
    realApp_poly_sub(ls, fes, fos, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(ls, ls);
    
    realApp_poly_clear(fe);
    realApp_poly_clear(fo);
    realApp_poly_clear(fes);
    realApp_poly_clear(fos);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}

void realApp_poly_oneGraeffeIteration_with_lastTerms_inPlace(realApp_poly_t f, 
                                                              const realApp_poly_t lastTerms, 
                                                              const realRat_t radius, 
                                                              slong split, 
                                                              slong prec,
                                                              metadatas_t meta){
    
    clock_t start = clock();
    /* scale lastTerms */
    realRat_t radS;
    realRat_init(radS);
    realRat_pow_si(radS, radius, 2);
    
    realApp_poly_t lsScaled;
    realApp_poly_init(lsScaled);
    realApp_poly_set(lsScaled, lastTerms);
    realApp_poly_scale_realRat_in_place( lsScaled->coeffs, radS, lsScaled->length, prec);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
    start = clock();
    
    realApp_poly_t f1, f2;
    realApp_poly_init(f1);
    realApp_poly_init(f2);
    
    realApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    realApp_poly_t fe1, fe2, fo1, fo2;
    realApp_poly_init2(fe1, len2);
    realApp_poly_init2(fe2, len2);
    realApp_poly_init2(fo1, len2);
    realApp_poly_init2(fo2, len2);
    realApp_ptr fe1ptr = fe1->coeffs;
    realApp_ptr fe2ptr = fe2->coeffs;
    realApp_ptr fo1ptr = fo1->coeffs;
    realApp_ptr fo2ptr = fo2->coeffs;
    
    for (i = 0; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) {
            if (i<=split)
                realApp_set( fe1ptr + quo, fptr+i);
            else
                realApp_set( fe2ptr + quo, fptr+i);
        }
        else {
            if (i<=split)
                realApp_set( fo1ptr + quo, fptr+i);
            else
                realApp_set( fo2ptr + quo, fptr+i);
        }
    }
    realApp_poly_set_length(fe1, len2);
    realApp_poly_set_length(fe2, len2);
    realApp_poly_set_length(fo1, len2);
    realApp_poly_set_length(fo2, len2);
    
    realApp_poly_t fes, fos;
    realApp_poly_init2(fes, len1);
    realApp_poly_init2(fos, len1);
    /* compute f1 */
    realApp_poly_mullow( fes, fe1, fe1, len1, prec);
    realApp_poly_mullow( fos, fo1, fo1, len1, prec);
    realApp_poly_shift_left( fos, fos, 1 );
    realApp_poly_sub(f1, fes, fos, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(f1, f1);
    /* compute f2 */
    realApp_poly_mullow( fes, fe1, fe2, len1, prec);
    realApp_poly_mullow( fos, fo1, fo2, len1, prec);
    realApp_poly_shift_left( fos, fos, 1 );
    realApp_poly_sub(f2, fes, fos, prec);
    realApp_poly_add(f2, f2, f2, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(f1, f1);
    
    realApp_poly_add(f, f1, f2, prec);
    realApp_poly_add(f, f, lsScaled, prec);
    if ((len1%2)==0)
        realApp_poly_neg(f, f);
    
    realRat_clear(radS);
    realApp_poly_clear(f1);
    realApp_poly_clear(f2);
    realApp_poly_clear(lsScaled);
    realApp_poly_clear(fe1);
    realApp_poly_clear(fe2);
    realApp_poly_clear(fo1);
    realApp_poly_clear(fo2);
    realApp_poly_clear(fes);
    realApp_poly_clear(fos);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}

void deflate_connCmp_init  (connCmp_t x){
    realApp_poly_init( connCmp_defPoref(x) );
    realApp_poly_init( connCmp_defFGref(x) );
}

void deflate_connCmp_clear (connCmp_t x){
    connCmp_isDefref(x) = 0;
    connCmp_degDeref(x) = 0;
    connCmp_isDFGref(x) = 0;
    realApp_poly_clear( connCmp_defPoref(x) );
    realApp_poly_clear( connCmp_defFGref(x) );
}



void deflate_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta ){
    
    connCmp_isDefref(x) = 1;
    connCmp_degDeref(x) = nbSols;
    /* compute the interval taylor shift at any point in the disk */
    realApp_poly_set(connCmp_defPoref(x), cacheApp_getApproximation_real ( cache, prec ));
    deflate_taylor_shift_interval_inplace( connCmp_defPoref(x), disk, prec, meta);
    
}

void deflate_copy( connCmp_t dest, const connCmp_t src ){
    if (connCmp_isDefref(src) != 0) {
        connCmp_isDefref(dest) = connCmp_isDefref(src);
        connCmp_degDeref(dest) = connCmp_degDeref(src);
        connCmp_isDFGref(dest) = connCmp_isDFGref(src);
        realApp_poly_set( connCmp_defPoref(dest), connCmp_defPoref(src) );
        realApp_poly_set( connCmp_defFGref(dest), connCmp_defFGref(src) );
    }
}

void deflate_compute_trailing_coeffs(realApp_ptr coeffs, const connCmp_t x, cacheApp_t cache, const compDsk_t d, slong prec, metadatas_t meta){
    
    int nbCoeffs = connCmp_degDeref(x) +1;
    /* pol */
    realApp_poly_t fapprox;
    realApp_poly_init(fapprox);
    realApp_poly_set(fapprox, cacheApp_getApproximation_real ( cache, prec ));
    
    
    realApp_t center, coeff;
    realRat_t factor;
    realApp_init(center);
    realApp_init(coeff);
    realRat_init(factor);
    
    realRat_set_si(factor, 1, 1);
    realApp_set_realRat(center, compRat_realref(compDsk_centerref(d)), prec);
    
    for (int index=0; index<nbCoeffs; index ++) {
        deflate_evaluate(coeffs+index, fapprox, center, prec, meta);
        realApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
        
        if (index<nbCoeffs){
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
            deflate_derivative_inplace( fapprox, prec, meta);
        }
        
    }
    
    realApp_poly_clear(fapprox);
    realApp_clear(center);
    realApp_clear(coeff);
    realRat_clear(factor);
}

void deflate_compute_leading_coeffs(realApp_ptr coeffs, const connCmp_t x, const compDsk_t d, slong prec, metadatas_t meta){
    
    realApp_t factor, temp;
    
    realApp_init(factor);
    realApp_init(temp);
    
    clock_t start = clock();
    
    realApp_set_realRat(temp, compDsk_radiusref(d), prec);
    realApp_pow_ui( factor, temp, connCmp_degDeref(x) + 1, prec);
    for (int index = connCmp_degDeref(x) + 1; index < connCmp_defPoref(x)->length; index ++ ) {
        
        realApp_mul( coeffs + index, connCmp_defPoref(x)->coeffs + index, factor, prec);
        realApp_mul(factor, factor, temp, prec);
    }
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
    realApp_clear(factor);
    realApp_clear(temp);
    
}

tstar_res deflate_tstar_test( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta) {
    
    clock_t start = clock();
    
    tstar_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    
    slong deg = cacheApp_getDegree(cache);
    int N = 0;
    int restemp = 0;
    int nbGraeffe = 0;
    int iteration = 0;
    
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,deg+1);
    pApprox->length = deg+1;
    realApp_t sum;
    realApp_init(sum);
    
    realApp_t coeff0, coeff1, coeffn; /* for anticipate */
    int anticipate_already_applied = 0;
    N = (int) 4+ceil(log2(1+log2(deg)));
    
    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
    deflate_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
    /*compute the trailing coefficients */
    deflate_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
    
    while( (iteration <= N)&&(restemp==0) ){
        
        if (iteration >= 1) {
            
//             if (iteration==1){
//                 
//                 if (connCmp_isDFGref(CC) == 0) {
//                     connCmp_isDFGref(CC) =1;
//                     realApp_poly_init2(connCmp_defFGref(CC) , deg+1);
//                     realApp_poly_oneGraeffeIteration_lastTerms(connCmp_defFGref(CC), connCmp_defPoref(CC), connCmp_degDeref(CC)+1, res.appPrec, meta );
//                 }
//                 realApp_poly_oneGraeffeIteration_with_lastTerms_inPlace(pApprox, connCmp_defFGref(CC), compDsk_radiusref(d),
//                                                                              connCmp_degDeref(CC)+1, res.appPrec, meta);
//             }
//             else
                deflate_real_graeffe_iterations_inplace( pApprox, 1, res.appPrec, meta);
            nbGraeffe +=1;
        }
        
        realApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
        res.nbOfSol = -1;
        while( (res.nbOfSol < connCmp_degDeref(CC))&&(restemp==0) ){
            res.nbOfSol += 1;
            restemp = realApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
            if (metadatas_getVerbo(meta)>=3)
                printf("Pellet test, discard: %d, %d-th coeff: %d, %i-th Graeffe it\n", discard, res.nbOfSol, restemp, iteration);
            if (restemp==-2){
                if (res.appPrec == prec) {
                    res.appPrec *=2;
                    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
                    deflate_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
                    /*compute the trailing coefficients */
                    deflate_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
                    deflate_real_graeffe_iterations_inplace( pApprox, iteration, res.appPrec, meta);
                    realApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
                    restemp = realApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
                    if (metadatas_getVerbo(meta)>=3)
                        printf("Pellet test, discard: %d, %d-th coeff: %d, %i-th Graeffe it\n", discard, res.nbOfSol, restemp, iteration);
                    
                }
            }

            if ( (discard) && (metadatas_useAnticipate(meta)) && (anticipate_already_applied==0) && (restemp == 0) ) {
                    int test_anticipate = ((0x1 << (N-iteration)) <= (realApp_poly_degree(pApprox)/4));
                    if (test_anticipate) {
                        clock_t start2 = clock();
                
                        anticipate_already_applied = 1;
                        realApp_init(coeff0);
                        realApp_init(coeff1);
                        realApp_init(coeffn);
                        
                        deflate_real_graeffe_iterations_abs_two_first_coeffs( coeff0, coeff1, pApprox, N-iteration, res.appPrec, meta);
                        realApp_abs( coeffn, realApp_poly_getCoeff(pApprox, realApp_poly_degree(pApprox)) );
                        realApp_pow_ui( coeffn, coeffn, (ulong)(0x1 << (N-iteration)), res.appPrec);
                        realApp_add(coeffn, coeffn, coeff1, res.appPrec);
                        restemp = realApp_soft_compare( coeff0, coeffn, res.appPrec );
                        
                        if (restemp==0){
//                             printf("la, nbMSol = %d, %i-th Graeffe\n", max_nb_sols, iteration);
                            restemp = -1;
                        }
                        else
                            restemp = 0;
                
                        realApp_clear(coeff0);
                        realApp_clear(coeff1);
                        realApp_clear(coeffn);

                        if (metadatas_haveToCount(meta)) 
                            metadatas_add_time_Anticip(meta, (double) (clock() - start2) );
                    }
            }
        
        }
        iteration +=1;
            
    }
    
    realApp_poly_clear(pApprox);
    realApp_clear(sum);
    
    if ((restemp==0)||(restemp==-2)) res.nbOfSol = -2;
    if ((restemp==-1)) res.nbOfSol = -1;
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefTsta(meta, (double) (clock() - start) );
    
    return res;
}

/*DEPRECATED*/
// void deflate2_set_clearance (connCmp_t x, const compDsk_t disk){
//     compDsk_set(connCmp_defClref(x), disk);
//     connCmp_defiCref(x) = 1;
// }
// 
// void deflate2_copy_clearance (connCmp_t dest, const connCmp_t src){
//     compDsk_set(connCmp_defClref(dest), connCmp_defClref(src));
//     connCmp_defiCref(dest) = connCmp_defiCref(src);
// }
// 
// int deflate2_isDef_clearance (const connCmp_t x){
//     return (connCmp_defiCref(x)==1);
// }
// 
// 
// void deflate2_connCmp_init  (connCmp_t x){
//     compDsk_init     (connCmp_defDsref(x));
//     compApp_poly_init( connCmp_defPfref(x) );
//     fmpz_init        (connCmp_defIRref(x));
// }
// 
// void deflate2_connCmp_clear (connCmp_t x){
//     compDsk_clear     (connCmp_defDsref(x));
//     compApp_poly_clear( connCmp_defPfref(x) );
//     fmpz_clear        (connCmp_defIRref(x));
//     connCmp_defiDref(x)=0;
// }
// 
// void deflate2_set (connCmp_t x, const compDsk_t disk, const fmpz_t isoRatio){
// //     connCmp_degDeref(x) = nbSols;
//     compDsk_set(connCmp_defDsref(x), disk);
//     fmpz_set   (connCmp_defIRref(x), isoRatio);
//     connCmp_defiDref(x)=1;
// }
// 
// void deflate2_copy (connCmp_t dest, const connCmp_t src){
//     compDsk_set(connCmp_defDsref(dest), connCmp_defDsref(src));
//     compApp_poly_set( connCmp_defPfref(dest), connCmp_defPfref(src));
//     fmpz_set   (connCmp_defIRref(dest), connCmp_defIRref(src));
//     connCmp_defiDref(dest)=connCmp_defiDref(src);
//     connCmp_defPrref(dest)=connCmp_defPrref(src);
// }
// 
// slong deflate2_get_Nb_EvalPoints( const fmpz_t isoRatio, slong degree, slong nbPs, slong prec, metadatas_t meta ) {
//     
//     realApp_t liR, wP, nbP;
//     realApp_init(liR);
//     realApp_init(wP);
//     realApp_init(nbP);
//     /* log(1/isoRatio) */
//     realApp_set_fmpz(liR, isoRatio, prec );
//     realApp_log( liR, liR, prec );
//     realApp_neg(liR, liR);
//     /* wanted precision */
//     realApp_set_d(wP, 0.5 );
//     realApp_pow_ui(wP, wP, (ulong) prec, prec);
//     realApp_log( nbP, wP, prec );
//     realApp_div( nbP, nbP, liR, prec);
//     realApp_add_si(wP, wP, degree, prec);
//     realApp_log( wP, wP, prec );
//     realApp_div( wP, wP, liR, prec);
//     realApp_sub( nbP, nbP, wP, prec);
//     slong q = realApp_ceil_si(nbP, prec);
//     q = q + nbPs;
//     
//     printf("Nb of required evaluation points (old): %ld\n", q);
//     
//     realApp_set_fmpz(liR, isoRatio, prec );
//     realApp_log_base_ui( liR, liR, 2, prec );
//     realApp_set_si(nbP, degree);
//     realApp_log_base_ui( nbP, nbP, 2, prec );
//     realApp_add_si(nbP, nbP, prec+1, prec);
//     realApp_div(nbP, nbP, liR, prec);
//     realApp_add_si(nbP, nbP, nbPs, prec);
//     q = realApp_ceil_si(nbP, prec);
//     printf("Nb of required evaluation points (new): %ld\n", q);
//     
//     realApp_clear(liR);
//     realApp_clear(wP);
//     realApp_clear(nbP);
//     
//     return q;
// }
// 
// void deflate2_compute_factor( connCmp_t x, cacheApp_t cache, slong prec, metadatas_t meta ){
//     
//     prec = prec;
//     /* number of evaluation points */
//     slong nbPoints = deflate2_get_Nb_EvalPoints( connCmp_defIRref(x), cacheApp_getDegree(cache), connCmp_nSolsref(x), prec, meta );
//     slong nbPowerS = connCmp_nSolsref(x);
//     
//     compApp_ptr points;
//     compApp_ptr pointsShifted;
//     compApp_ptr fvals;
//     compApp_ptr fdervals;
//     compApp_ptr fdivs;
//     compApp_ptr ps;
//     
//     points =        (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     pointsShifted = (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     fvals =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     fdervals =      (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     fdivs =         (compApp_ptr) ccluster_malloc( nbPoints*sizeof(compApp) );
//     
//     ps =            (compApp_ptr) ccluster_malloc( nbPowerS*sizeof(compApp) );
//     
//     for (slong i=0; i<nbPoints; i++){
//         compApp_init( points +i );
//         compApp_init( pointsShifted +i );
//         compApp_init( fvals +i );
//         compApp_init( fdervals +i );
//         compApp_init( fdivs +i );
//     }
//     
//     for (slong j=0; j<nbPowerS; j++)
//         compApp_init( ps +j );
//     
//     slong precRes = deflate_computePsApprox(ps, compDsk_centerref(connCmp_defDsref(x)), compDsk_radiusref(connCmp_defDsref(x)),
//                                             points, pointsShifted, fvals, fdervals, fdivs, nbPoints, nbPowerS, cache, prec, meta);
//     
//     if (metadatas_getVerbo(meta)>=3){
//         printf("asked prec for ps. approx: %ld, out prec: %ld\n", prec, precRes);
//         for (slong j = 1; j<=nbPowerS; j++) {
//             printf("--- %d-th power sum approximation: ", (int) j);
//             compApp_printd( ps+j-1, 10 ); printf("\n");
//         }
//     }
//     
//     /* compute error */
//     realApp_t error, iir, iirqph, iirqmh, denom;
//     realApp_init(error);
//     realApp_init(iir);
//     realApp_init(iirqph);
//     realApp_init(iirqmh);
//     realApp_init(denom);
//     realApp_set_fmpz(iir, connCmp_defIRref(x), prec);
//     realApp_inv(iir, iir, prec);
//     realApp_pow_ui(iirqph, iir, nbPoints, prec);
//     realApp_set(iirqmh, iirqph);
//     realApp_one(denom);
//     realApp_sub(denom, denom, iirqph, prec);
//     realApp_mul_si(iirqph, iirqph, connCmp_nSolsref(x), prec);
//     realApp_mul_si(iirqmh, iirqmh, cacheApp_getDegree(cache) - connCmp_nSolsref(x), prec);
//     realApp_div(iirqph, iirqph, denom, prec);
//     realApp_div(iirqmh, iirqmh, denom, prec);
//     
//     for (slong j = 1; j<=nbPowerS; j++) {
//         realApp_mul(iirqph, iirqph, iir, prec);
//         realApp_div(iirqmh, iirqmh, iir, prec);
//         realApp_add(error, iirqph, iirqmh, prec);
//         realApp_add_error( compApp_realref(ps+j-1), error );
//         realApp_add_error( compApp_imagref(ps+j-1), error );
//         if (metadatas_getVerbo(meta)>=3){
//             printf("error of %ld-th power sum: ", j);
//             realApp_printd(error, 10);
//             printf("\n");
//         }
//     }
//     
//     if (metadatas_getVerbo(meta)>=3){
//         printf("with error: \n");
//         for (slong j = 1; j<=nbPowerS; j++) {
//             printf("--- %d-th power sum approximation: ", (int) j);
//             compApp_printd( ps+j-1, 10 ); printf("\n");
//         }
//     }
//     
//     /* compute factor with newton identities */
//     compApp_poly_init2( connCmp_defPfref(x), connCmp_nSolsref(x) + 1 );
//     connCmp_defPfref(x)->length = connCmp_nSolsref(x) + 1;
//     int index = connCmp_nSolsref(x);
//     compApp_one( connCmp_defPfref(x)->coeffs + index );
//     for (slong j = 1; j<=nbPowerS; j++) {
//         index = index - 1;
//         compApp_set( connCmp_defPfref(x)->coeffs + index, ps+j-1 );
//         for (slong i = 1; i < j; i++) {
// //             compApp_addmul( connCmp_defPfref(x)->coeffs + index, ps+j-1-i, connCmp_defPfref(x)->coeffs + index+i, prec );
//             compApp_addmul( connCmp_defPfref(x)->coeffs + index, ps+j-1-i, connCmp_defPfref(x)->coeffs + connCmp_nSolsref(x)-i, prec );
//         }
//         compApp_div_si(connCmp_defPfref(x)->coeffs + index, connCmp_defPfref(x)->coeffs + index, j, prec);
//         compApp_neg( connCmp_defPfref(x)->coeffs + index, connCmp_defPfref(x)->coeffs + index );
//     }
//     
//     if (metadatas_getVerbo(meta)>=3){
//         printf("factor: \n");
//         compApp_poly_printd( connCmp_defPfref(x), 10 );
//         printf("\n");
//         slong acc = compApp_poly_get_relOne_accuracy_min( connCmp_defPfref(x));
//         printf("accuracy: %ld\n", acc);
//     }
//     
//     /* shift to 0+0i with rad 1 */
//     realRat_t invR;
//     realRat_init(invR);
//     realRat_inv(invR, compDsk_radiusref(connCmp_defDsref(x)));
//     compApp_poly_scale_realRat_in_place( connCmp_defPfref(x)->coeffs, invR, connCmp_defPfref(x)->length, prec );
//     compRat_t mcenter;
//     compRat_init(mcenter);
//     compRat_neg(mcenter, compDsk_centerref(connCmp_defDsref(x)));
// //     compRat_set(mcenter, compDsk_centerref(connCmp_defDsref(x)));
//     compApp_poly_taylorShift_in_place_noscale( connCmp_defPfref(x), mcenter, prec );
//     
//     compRat_clear(mcenter);
//     realRat_clear(invR);
//     connCmp_defiDref(x)=1;
//     
// //     compApp_t numb;
// //     compApp_init(numb);
// //     compApp_set(numb, connCmp_defPfref(x)->coeffs + (connCmp_defPfref(x)->length -1));
// // //     realApp_zero(compApp_imagref(numb));
// //     compApp_one(connCmp_defPfref(x)->coeffs + (connCmp_defPfref(x)->length -1));
// //     for (slong index=0; index < (connCmp_defPfref(x)->length-1); index++){
// // //         realApp_zero(compApp_imagref(connCmp_defPfref(x)->coeffs + index));
// //         compApp_div(connCmp_defPfref(x)->coeffs + index, connCmp_defPfref(x)->coeffs + index, numb, prec);
// //         
// //     }
// //     
// //     compApp_clear(numb);
//     
//     if (metadatas_getVerbo(meta)>=3){
//         printf("factor: \n");
//         compApp_poly_printd( connCmp_defPfref(x), 100 );
//         printf("\n");
//         slong acc = compApp_poly_get_relOne_accuracy_min( connCmp_defPfref(x));
//         printf("accuracy: %ld\n", acc);
//     }
//     
//     realApp_clear(error);
//     realApp_clear(iir);
//     realApp_clear(iirqph);
//     realApp_clear(iirqmh);
//     realApp_clear(denom);
//     
//     for (int i=0; i<nbPoints; i++){
//         compApp_clear( points +i );
//         compApp_clear( pointsShifted +i );
//         compApp_clear( fvals +i );
//         compApp_clear( fdervals +i );
//         compApp_clear( fdivs +i );
//     }
//     
//     for (int j=0; j<nbPowerS; j++)
//         compApp_clear( ps +j );
//     
//     ccluster_free(points);
//     ccluster_free(pointsShifted);
//     ccluster_free(fvals);
//     ccluster_free(fdervals);
//     ccluster_free(fdivs);
//     ccluster_free(ps);
//     
// }
// 
// tstar_res deflate2_tstar_test( const connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, slong prec, metadatas_t meta) {
//     
//     tstar_res res;
//     res.nbOfSol = -1;
//     res.appPrec = prec;
//     
// //     realApp_ptr coeffs;
//     realApp_t sum;
//     realApp_init(sum);
//     
//     compApp_poly_t pshift;
//     compApp_poly_init(pshift);
//     
// //     compRat_t shiftedCenter;
// //     compRat_init(shiftedCenter);
// //     compRat_set(shiftedCenter, compDsk_centerref(d));
// //     realRat_sub( compRat_realref(shiftedCenter), compRat_realref(shiftedCenter), compRat_realref(compDsk_centerref(connCmp_defDsref(CC))) );
// //     realRat_sub( compRat_imagref(shiftedCenter), compRat_imagref(shiftedCenter), compRat_imagref(compDsk_centerref(connCmp_defDsref(CC))) );
//     
//     compApp_poly_taylorShift( pshift, connCmp_defPfref(CC), 
//                            compDsk_centerref(d), 
//                            compDsk_radiusref(d), 
//                            prec );
//     
//     printf("shifted factor: \n");
//     compApp_poly_printd( pshift, 10 );
//     printf("\n");
//     slong acc = compApp_poly_get_relOne_accuracy_min( pshift);
//     printf("accuracy: %ld\n", acc);
//     
//     compApp_poly_sum_abs_coeffs( sum, pshift, prec );
//     
//     int ind = 0;
//     int resPellet=0;
//     while ( ( ind <= connCmp_nSolsref(CC) ) 
//          && ( resPellet != 1 ) ) {
//         resPellet = compApp_poly_TkGtilda_with_sum( pshift, sum, ind, prec);
//         printf("Pellet test, %d-th coeff: %d\n", ind, resPellet);
// // //         if (resPellet == -2) {
// //             printf("---%d-th coeff: ", ind);
// //             realApp_printd(coeffs + ind, res.appPrec);
// //             printf("\n---      sum  : ");
// //             realApp_printd(sum, res.appPrec);
// //             printf("\n");
// // //         }
//         ind++;
//     }
//     
//     if (resPellet==1)
//         res.nbOfSol = ind -1;
//     else 
//         res.nbOfSol = -2;
//     
//     compApp_poly_clear(pshift);
//     realApp_clear(sum);
// 
//     return res;
// }
// 
// void deflate_getEvaluationPoints( compApp_ptr points, 
//                                     compApp_ptr pointsShifted,
//                                     const compRat_t center,
//                                     const realRat_t radius,
//                                     slong nbPoints,
//                                     slong prec ) {
//     compApp_t c, a;
//     realRat_t argu;
//     
//     compApp_init(c);
//     compApp_init(a);
//     realRat_init(argu);
//     
//     compApp_set_compRat(c, center, prec);
//     for(slong i=0; i<nbPoints; i++) {
//         realRat_set_si(argu, 2*i, nbPoints);
//         compApp_set_realRat(a, argu, prec);
//         acb_exp_pi_i( points + i, a, prec);
// //         compApp_mul_realRat_in_place(points + i, radius, prec);
// //         compApp_add( pointsShifted + i, c, points + i, prec);
//         compApp_mul_realRat(pointsShifted + i, points + i, radius, prec);
//         compApp_add( pointsShifted + i, c, pointsShifted + i, prec);
//     }
//     
//     compApp_clear(c);
//     compApp_clear(a);
//     realRat_clear(argu);
// }
// 
// void deflate_evaluateAtPoints( compApp_ptr f_val,
//                                  compApp_ptr fder_val,
//                                  const compApp_ptr points,
//                                  slong nbPoints,
//                                  cacheApp_t cache,
//                                  slong prec,
//                                  metadatas_t meta){
//     
//     if (metadatas_pwSumref(meta)->evalPoly == NULL) {
//         compApp_poly_ptr app = cacheApp_getApproximation ( cache, prec );
//         for (slong i=0; i<nbPoints; i++)
//             compApp_poly_evaluate2_rectangular(f_val + i, fder_val + i, app, points + i, prec);
//     }
//     else {
//         for (slong i=0; i<nbPoints; i++)
//             metadatas_pwSumref(meta)->evalPoly( f_val+i, fder_val + i, points+i, prec);
//     }
// }
// 
// /* returns -1: should increase precision
//  *          1: OK! */
// int deflate_computePsApprox_fromVals(compApp_ptr ps,
//                                         const realRat_t radius,
//                                         compApp_ptr points,
//                                         compApp_ptr fvals,
//                                         compApp_ptr fdervals,
//                                         compApp_ptr fdivs,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
//                                         slong prec,
//                                         metadatas_t meta){
//     
//     int res=1;
//     
//     /* compute fdivs; check if fvals contains zero*/
//     for (slong i = 0; (i<nbPoints) && (res==1) ; i++) {
//         if (compApp_contains_zero( fvals +i )){
//             res = -1;
//         }
//         compApp_div(fdivs +i, fdervals + i, fvals + i, prec);
//     }
//     
//     if (res==1){
//         
//         realRat_t factor;
//         realRat_init(factor);
//         realRat_set(factor, radius);
// //         realRat_set_si(factor, 1, 1);
//         realRat_div_ui( factor, factor, nbPoints );
//         
//         
//         /* compute powerSums */
//         for (slong j = 1; j<=nbPowerSums; j++)
//             compApp_mul(ps+(j-1), fdivs + 0, points + 0, prec);
//         for (slong i = 1; i<nbPoints; i++)
//             for (slong j = 1; j<=nbPowerSums; j++) /* from first power sum */
//                 compApp_addmul(ps+(j-1) , fdivs + i, points + (((j+1)*i)%nbPoints), prec);  
//         for (slong j = 1; j<=nbPowerSums; j++) {
// //             realRat_mul(factor, factor, radius);
//             compApp_mul_realRat(ps+(j-1), ps+(j-1), factor, prec);
//         }
//         
//         
//         realRat_clear(factor);
//         
//     }
//     
//     return res;
// }
// 
// slong deflate_computePsApprox(compApp_ptr ps,
//                                         const compRat_t center,
//                                         const realRat_t radius,
//                                         compApp_ptr points,
//                                         compApp_ptr pointsShifted,
//                                         compApp_ptr fvals,
//                                         compApp_ptr fdervals,
//                                         compApp_ptr fdivs,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
//                                         cacheApp_t cache,
//                                         slong prec,
//                                         metadatas_t meta){
//     slong precRes = prec;
//     int resTemp = -1;
//     
//     /* compute points and evals at prec res.appPrec*/
//     deflate_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, precRes);
// //     clock_t start = clock();
//     deflate_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, precRes, meta);
// /*    if (metadatas_haveToCount(meta))
//             metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start) );  */ 
// 
//     /* compute approximation of Power sums */
//     resTemp = deflate_computePsApprox_fromVals(ps, radius, points, fvals, fdervals, fdivs, nbPoints, nbPowerSums, precRes, meta);
//     
//     while ( resTemp ==-1 ) {
//         precRes = 2*precRes;
//         
//         deflate_getEvaluationPoints( points, pointsShifted, center, radius, nbPoints, precRes);
// //         clock_t start2 = clock();
//         deflate_evaluateAtPoints( fvals, fdervals, pointsShifted, nbPoints, cache, precRes, meta);
// //         if (metadatas_haveToCount(meta))
// //             metadatas_add_Evals( meta, depth, nbPoints, (double) (clock() - start2) );
//         /* compute approximation of Power sums */
//         resTemp = deflate_computePsApprox_fromVals(ps, radius, points, fvals, fdervals, fdivs, nbPoints, nbPowerSums, precRes, meta);
//     }
//     
//     return precRes;
// }
