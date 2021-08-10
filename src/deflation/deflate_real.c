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

void deflate_real_taylor_shift_interval_inplace( realApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta){
    
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
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}

void deflate_real_derivative_inplace( realApp_poly_t res, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        realApp_poly_derivative(res, res, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefDeri(meta, (double) (clock() - start) );
}

void deflate_real_evaluate(realApp_t y, const realApp_poly_t f, const realApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
//         arb_poly_evaluate_horner(y, f, x, prec);
        realApp_poly_evaluate(y, f, x, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefEval(meta, (double) (clock() - start) );
}

void deflate_real_evaluate2(realApp_t y, realApp_t z, const realApp_poly_t f, const realApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
//         arb_poly_evaluate_horner(y, f, x, prec);
        realApp_poly_evaluate2(y, z, f, x, prec);

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

void deflate_real_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta ){
    
    connCmp_isDefref(x) = 1;
    connCmp_degDeref(x) = nbSols;
    /* compute the interval taylor shift at any point in the disk */
    realApp_poly_set(connCmp_defPoRref(x), cacheApp_getApproximation_real ( cache, prec ));
    deflate_real_taylor_shift_interval_inplace( connCmp_defPoRref(x), disk, prec, meta);
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("deflate.c: deflate_tstar_test: Interval polynomial: \n");
//         realApp_poly_printd(connCmp_defPoref(x), 10);
//         printf("\n\n");
//     }
    /* after the taylor shift: all leading coeffs (index > nbSols) that contain zero */
    /* are set to a ball centered in zero containing the coeff */
    realApp_ptr coeffs = connCmp_defPoRref(x)->coeffs;
    for ( slong index = connCmp_degDeref(x) +1; index < connCmp_defPoRref(x)->length; index ++){
        if ( realApp_contains_zero( coeffs + index ) )
            realApp_center_in_zero( coeffs + index );
    }
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("deflate.c: deflate_tstar_test: Interval polynomial after re-center: \n");
//         realApp_poly_printd(connCmp_defPoRref(x), 10);
//         printf("\n\n");
//     }
}

void deflate_real_compute_trailing_coeffs(realApp_ptr coeffs, const connCmp_t x, cacheApp_t cache, const compDsk_t d, slong prec, metadatas_t meta){
    
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
        
        if ((index+1)<nbCoeffs) {
            deflate_real_evaluate2(coeffs+index, coeffs+(index+1), fapprox, center, prec, meta);
            realApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
            deflate_real_derivative_inplace( fapprox, prec, meta);
            index = index + 1;
            realApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
            
        } else {
            deflate_real_evaluate(coeffs+index, fapprox, center, prec, meta);
            realApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
        }
        
        if (index<nbCoeffs){
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
            deflate_real_derivative_inplace( fapprox, prec, meta);
        }
        
    }
    
    realApp_poly_clear(fapprox);
    realApp_clear(center);
    realApp_clear(coeff);
    realRat_clear(factor);
}

void deflate_real_compute_leading_coeffs(realApp_ptr coeffs, const connCmp_t x, const compDsk_t d, slong prec, metadatas_t meta){
    
    realApp_t factor, temp;
    
    realApp_init(factor);
    realApp_init(temp);
    
    clock_t start = clock();
    
    realApp_set_realRat(temp, compDsk_radiusref(d), prec);
    realApp_pow_ui( factor, temp, connCmp_degDeref(x) + 1, prec);
    for (int index = connCmp_degDeref(x) + 1; index < connCmp_defPoRref(x)->length; index ++ ) {
        
        realApp_mul( coeffs + index, connCmp_defPoRref(x)->coeffs + index, factor, prec);
        realApp_mul(factor, factor, temp, prec);
    }
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
    realApp_clear(factor);
    realApp_clear(temp);
    
}

tstar_res deflate_real_tstar_test( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta) {
    
    clock_t start = clock();
    
    if (metadatas_getVerbo(meta)>=4) {
        compApp_t c;
        realApp_t r;
        compApp_init(c);
        realApp_init(r);
        compApp_set_compRat(c, compDsk_centerref(d), CCLUSTER_DEFAULT_PREC);
        realApp_set_realRat(r, compDsk_radiusref(d), CCLUSTER_DEFAULT_PREC);
        printf("\n#deflate.c: deflate_tstar_test: begin\n");
        printf("#---current disk: center: ");
        compApp_printd(c, 10);
        printf(" radius: ");
        realApp_printd(r, 10);
        printf("\n");
        compApp_clear(c);
        realApp_clear(r);
    }
    
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
    
//     realApp_t coeff0, coeff1, coeffn; /* for anticipate */
//     int anticipate_already_applied = 0;
    N = (int) 4+ceil(log2(1+log2(deg)));
    
    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
    deflate_real_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
    /*compute the leading coefficients */
    deflate_real_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
    
    while( (iteration <= N)&&(restemp==0) ){
        
        if (iteration >= 1) {
            
            if (iteration==1){
                
                if (connCmp_isDFGref(CC) == 0) {
                    connCmp_isDFGref(CC) =1;
                    realApp_poly_init2(connCmp_defFGRref(CC) , deg+1);
                    realApp_poly_oneGraeffeIteration_lastTerms(connCmp_defFGRref(CC), connCmp_defPoRref(CC), connCmp_degDeref(CC)+1, res.appPrec, meta );
                }
                realApp_poly_oneGraeffeIteration_with_lastTerms_inPlace(pApprox, connCmp_defFGRref(CC), compDsk_radiusref(d),
                                                                             connCmp_degDeref(CC)+1, res.appPrec, meta);
            }
            else
                deflate_real_graeffe_iterations_inplace( pApprox, 1, res.appPrec, meta);
            nbGraeffe +=1;
        }
        
        realApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
        int cond = 0;
        if (discard) {
            res.nbOfSol = -1;
            cond = (res.nbOfSol < connCmp_degDeref(CC));
        } else {
            res.nbOfSol = connCmp_degDeref(CC) +1;
            cond = (res.nbOfSol >= 1);
        }
        while( cond &&(restemp==0) ){
            if (discard) {
                res.nbOfSol += 1;
                cond = (res.nbOfSol < connCmp_degDeref(CC));
            } else {
                res.nbOfSol -= 1;
                cond = (res.nbOfSol >= 1);
            }
            restemp = realApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
            if (metadatas_getVerbo(meta)>=4)
                printf("#deflate.c: deflate_tstar_test: discard: %d, %d-th coeff: %d, %i-th Graeffe it, prec: %ld\n", discard, res.nbOfSol, restemp, nbGraeffe, res.appPrec);
            if (restemp==-2){
                if (res.appPrec == prec) {
                    res.appPrec *=2;
                    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
                    deflate_real_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
                    /*compute the trailing coefficients */
                    deflate_real_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
                    deflate_real_graeffe_iterations_inplace( pApprox, iteration, res.appPrec, meta);
                    realApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
                    restemp = realApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
                    if (metadatas_getVerbo(meta)>=4)
                        printf("#deflate.c: deflate_tstar_test: discard: %d, %d-th coeff: %d, %i-th Graeffe it, prec: %ld\n", discard, res.nbOfSol, restemp, nbGraeffe, res.appPrec);
                    
                }
            }
        
        }
        iteration +=1;
            
    }
    
    connCmp_reu_set_real( CC, compRat_realref( compDsk_centerref( d ) ), compDsk_radiusref( d ),
                                  nbGraeffe, res.appPrec, pApprox );
    
    realApp_poly_clear(pApprox);
    realApp_clear(sum);
    
    if ((restemp==0)||(restemp==-2)) res.nbOfSol = -2;
    if ((restemp==-1)) res.nbOfSol = -1;
    
    if (metadatas_haveToCount(meta))
//         if (discard)
            metadatas_add_time_DefTsta(meta, (double) (clock() - start) );
    
    if (metadatas_getVerbo(meta)>=3){
        printf("#deflate.c: deflate_tstar_test: discard: %d, res: %d, %i Graeffe its, prec: %ld\n", discard, res.nbOfSol, nbGraeffe, res.appPrec);
//         printf("#deflate.c: deflate_tstar_test: end\n\n");
    }
    
    return res;
}

tstar_res deflate_real_tstar_test_rescale( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta){
        return deflate_tstar_test_rescale( CC, cache, d, max_nb_sols, discard, prec, meta);
}

/* DEPRECATED */
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
