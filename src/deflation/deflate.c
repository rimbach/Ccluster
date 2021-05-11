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

// void deflate_connCmp_init  (connCmp_t x){
//     realApp_poly_init( connCmp_defPoRref(x) );
//     realApp_poly_init( connCmp_defFGRref(x) );
//     compApp_poly_init( connCmp_defPoCref(x) );
//     compApp_poly_init( connCmp_defFGCref(x) );
// }
// 
// void deflate_connCmp_clear (connCmp_t x){
//     connCmp_isDefref(x) = 0;
//     connCmp_degDeref(x) = 0;
//     connCmp_isDFGref(x) = 0;
//     realApp_poly_clear( connCmp_defPoRref(x) );
//     realApp_poly_clear( connCmp_defFGRref(x) );
//     compApp_poly_clear( connCmp_defPoCref(x) );
//     compApp_poly_clear( connCmp_defFGCref(x) );
// }

// void deflate_copy( connCmp_t dest, const connCmp_t src ){
//     if (connCmp_isDefref(src) != 0) {
//         connCmp_isDefref(dest) = connCmp_isDefref(src);
//         connCmp_degDeref(dest) = connCmp_degDeref(src);
//         connCmp_isDFGref(dest) = connCmp_isDFGref(src);
//         realApp_poly_set( connCmp_defPoRref(dest), connCmp_defPoRref(src) );
//         realApp_poly_set( connCmp_defFGRref(dest), connCmp_defFGRref(src) );
//         compApp_poly_set( connCmp_defPoCref(dest), connCmp_defPoCref(src) );
//         compApp_poly_set( connCmp_defFGCref(dest), connCmp_defFGCref(src) );
//     }
// }

void deflate_taylor_shift_interval_inplace( compApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_taylorShift_interval_in_place( res, compDsk_centerref(d), compDsk_radiusref(d), prec );

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefTayl(meta, (double) (clock() - start) );
}

void deflate_graeffe_iterations_inplace( compApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        for(int i = 0; i < N; i++)
            compApp_poly_oneGraeffeIteration_in_place( res, prec );
        
        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}

void deflate_derivative_inplace( compApp_poly_t res, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_derivative(res, res, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefDeri(meta, (double) (clock() - start) );
}

void deflate_evaluate(compApp_t y, const compApp_poly_t f, const compApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
        compApp_poly_evaluate(y, f, x, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefEval(meta, (double) (clock() - start) );
}

void deflate_evaluate2(compApp_t y, compApp_t z, const compApp_poly_t f, const compApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
//         arb_poly_evaluate_horner(y, f, x, prec);
        compApp_poly_evaluate2(y, z, f, x, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefEval(meta, (double) (clock() - start) );
}

void compApp_poly_oneGraeffeIteration_lastTerms(compApp_poly_t ls, const compApp_poly_t f, slong split, slong prec, metadatas_t meta ){
    
    clock_t start = clock();
    
    compApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    compApp_poly_t fe, fo;
    compApp_poly_init2(fe, len2);
    compApp_poly_init2(fo, len2);
    compApp_ptr feptr = fe->coeffs;
    compApp_ptr foptr = fo->coeffs;
    
    for (i = split; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) 
            compApp_set( feptr + quo, fptr+i);
        else
            compApp_set( foptr + quo, fptr+i);
    }
    compApp_poly_set_length(fe, len2);
    compApp_poly_set_length(fo, len2);
    
    compApp_poly_t fes, fos;
    compApp_poly_init2(fes, len1);
    compApp_poly_init2(fos, len1);
    compApp_poly_mullow( fes, fe, fe, len1, prec);
    compApp_poly_mullow( fos, fo, fo, len1, prec);
    compApp_poly_shift_left( fos, fos, 1 );
    compApp_poly_sub(ls, fes, fos, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(ls, ls);
    
    compApp_poly_clear(fe);
    compApp_poly_clear(fo);
    compApp_poly_clear(fes);
    compApp_poly_clear(fos);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}

void compApp_poly_oneGraeffeIteration_with_lastTerms_inPlace(compApp_poly_t f, 
                                                              const compApp_poly_t lastTerms, 
                                                              const realRat_t radius, 
                                                              slong split, 
                                                              slong prec,
                                                              metadatas_t meta){
    
    clock_t start = clock();
    /* scale lastTerms */
    realRat_t radS;
    realRat_init(radS);
    realRat_pow_si(radS, radius, 2);
    
    compApp_poly_t lsScaled;
    compApp_poly_init(lsScaled);
    compApp_poly_set(lsScaled, lastTerms);
    compApp_poly_scale_realRat_in_place( lsScaled->coeffs, radS, lsScaled->length, prec);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
    start = clock();
    
    compApp_poly_t f1, f2;
    compApp_poly_init(f1);
    compApp_poly_init(f2);
    
    compApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    compApp_poly_t fe1, fe2, fo1, fo2;
    compApp_poly_init2(fe1, len2);
    compApp_poly_init2(fe2, len2);
    compApp_poly_init2(fo1, len2);
    compApp_poly_init2(fo2, len2);
    compApp_ptr fe1ptr = fe1->coeffs;
    compApp_ptr fe2ptr = fe2->coeffs;
    compApp_ptr fo1ptr = fo1->coeffs;
    compApp_ptr fo2ptr = fo2->coeffs;
    
    for (i = 0; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) {
            if (i<=split)
                compApp_set( fe1ptr + quo, fptr+i);
            else
                compApp_set( fe2ptr + quo, fptr+i);
        }
        else {
            if (i<=split)
                compApp_set( fo1ptr + quo, fptr+i);
            else
                compApp_set( fo2ptr + quo, fptr+i);
        }
    }
    compApp_poly_set_length(fe1, len2);
    compApp_poly_set_length(fe2, len2);
    compApp_poly_set_length(fo1, len2);
    compApp_poly_set_length(fo2, len2);
    
    compApp_poly_t fes, fos;
    compApp_poly_init2(fes, len1);
    compApp_poly_init2(fos, len1);
    /* compute f1 */
    compApp_poly_mullow( fes, fe1, fe1, len1, prec);
    compApp_poly_mullow( fos, fo1, fo1, len1, prec);
    compApp_poly_shift_left( fos, fos, 1 );
    compApp_poly_sub(f1, fes, fos, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(f1, f1);
    /* compute f2 */
    compApp_poly_mullow( fes, fe1, fe2, len1, prec);
    compApp_poly_mullow( fos, fo1, fo2, len1, prec);
    compApp_poly_shift_left( fos, fos, 1 );
    compApp_poly_sub(f2, fes, fos, prec);
    compApp_poly_add(f2, f2, f2, prec);
//     if ((len1%2)==0)
//         realApp_poly_neg(f1, f1);
    
    compApp_poly_add(f, f1, f2, prec);
    compApp_poly_add(f, f, lsScaled, prec);
    if ((len1%2)==0)
        compApp_poly_neg(f, f);
    
    realRat_clear(radS);
    compApp_poly_clear(f1);
    compApp_poly_clear(f2);
    compApp_poly_clear(lsScaled);
    compApp_poly_clear(fe1);
    compApp_poly_clear(fe2);
    compApp_poly_clear(fo1);
    compApp_poly_clear(fo2);
    compApp_poly_clear(fes);
    compApp_poly_clear(fos);
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefGrae(meta, (double) (clock() - start) );
}




void deflate_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta ){
    
    connCmp_isDefref(x) = 1;
    connCmp_degDeref(x) = nbSols;
    /* compute the interval taylor shift at any point in the disk */
    compApp_poly_set(connCmp_defPoCref(x), cacheApp_getApproximation ( cache, prec ));
    deflate_taylor_shift_interval_inplace( connCmp_defPoCref(x), disk, prec, meta);
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("deflate.c: deflate_tstar_test: Interval polynomial: \n");
//         realApp_poly_printd(connCmp_defPoref(x), 10);
//         printf("\n\n");
//     }
    /* after the taylor shift: all leading coeffs (index > nbSols) that contain zero */
    /* are set to a ball centered in zero containing the coeff */
    /*TODO*/
//     compApp_ptr coeffs = connCmp_defPoCref(x)->coeffs;
//     for ( slong index = connCmp_degDeref(x) +1; index < connCmp_defPoCref(x)->length; index ++){
//         if ( compApp_contains_zero( coeffs + index ) )
//             compApp_center_in_zero( coeffs + index );
//     }
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("deflate.c: deflate_tstar_test: Interval polynomial after re-center: \n");
//         realApp_poly_printd(connCmp_defPoref(x), 10);
//         printf("\n\n");
//     }
}



void deflate_compute_trailing_coeffs(compApp_ptr coeffs, const connCmp_t x, cacheApp_t cache, const compDsk_t d, slong prec, metadatas_t meta){
    
    int nbCoeffs = connCmp_degDeref(x) +1;
    /* pol */
    compApp_poly_t fapprox;
    compApp_poly_init(fapprox);
    compApp_poly_set(fapprox, cacheApp_getApproximation( cache, prec ));
    
    
    compApp_t center, coeff;
    realRat_t factor;
    compApp_init(center);
    compApp_init(coeff);
    realRat_init(factor);
    
    realRat_set_si(factor, 1, 1);
    compApp_set_compRat(center, compDsk_centerref(d), prec);
    
    for (int index=0; index<nbCoeffs; index ++) {
        
        if ((index+1)<nbCoeffs) {
            deflate_evaluate2(coeffs+index, coeffs+(index+1), fapprox, center, prec, meta);
            compApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
            deflate_derivative_inplace( fapprox, prec, meta);
            index = index + 1;
            compApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
            
        } else {
            deflate_evaluate(coeffs+index, fapprox, center, prec, meta);
            compApp_mul_realRat(coeffs+index, coeffs+index, factor, prec);
        }
        
        if (index<nbCoeffs){
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
            deflate_derivative_inplace( fapprox, prec, meta);
        }
        
    }
    
    compApp_poly_clear(fapprox);
    compApp_clear(center);
    compApp_clear(coeff);
    realRat_clear(factor);
}

void deflate_compute_leading_coeffs(compApp_ptr coeffs, const connCmp_t x, const compDsk_t d, slong prec, metadatas_t meta){
    
    realApp_t factor, temp;
    
    realApp_init(factor);
    realApp_init(temp);
    
    clock_t start = clock();
    
    realApp_set_realRat(temp, compDsk_radiusref(d), prec);
    realApp_pow_ui( factor, temp, connCmp_degDeref(x) + 1, prec);
    for (int index = connCmp_degDeref(x) + 1; index < connCmp_defPoCref(x)->length; index ++ ) {
        
        compApp_mul_realApp( coeffs + index, connCmp_defPoCref(x)->coeffs + index, factor, prec);
        realApp_mul(factor, factor, temp, prec);
    }
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
    realApp_clear(factor);
    realApp_clear(temp);
    
}

tstar_res deflate_tstar_test( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta) {
    
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
    
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,deg+1);
    pApprox->length = deg+1;
    realApp_t sum;
    realApp_init(sum);
    
//     realApp_t coeff0, coeff1, coeffn; /* for anticipate */
//     int anticipate_already_applied = 0;
    N = (int) 4+ceil(log2(1+log2(deg)));
    
    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
    deflate_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
    /*compute the leading coefficients */
    deflate_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
    
    while( (iteration <= N)&&(restemp==0) ){
        
        if (iteration >= 1) {
            
            if (iteration==1){
                
                if (connCmp_isDFGref(CC) == 0) {
                    connCmp_isDFGref(CC) =1;
                    compApp_poly_init2(connCmp_defFGCref(CC) , deg+1);
                    compApp_poly_oneGraeffeIteration_lastTerms(connCmp_defFGCref(CC), connCmp_defPoCref(CC), connCmp_degDeref(CC)+1, res.appPrec, meta );
                }
                compApp_poly_oneGraeffeIteration_with_lastTerms_inPlace(pApprox, connCmp_defFGCref(CC), compDsk_radiusref(d),
                                                                             connCmp_degDeref(CC)+1, res.appPrec, meta);
            }
            else
                deflate_graeffe_iterations_inplace( pApprox, 1, res.appPrec, meta);
            nbGraeffe +=1;
        }
        
        compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
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
            restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
            if (metadatas_getVerbo(meta)>=4)
                printf("#deflate.c: deflate_tstar_test: discard: %d, %d-th coeff: %d, %i-th Graeffe it, prec: %ld\n", discard, res.nbOfSol, restemp, nbGraeffe, res.appPrec);
            if (restemp==-2){
                if (res.appPrec == prec) {
                    res.appPrec *=2;
                    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
                    deflate_compute_trailing_coeffs(pApprox->coeffs, CC, cache, d, res.appPrec, meta);
                    /*compute the trailing coefficients */
                    deflate_compute_leading_coeffs(pApprox->coeffs, CC, d, res.appPrec, meta);
                    deflate_graeffe_iterations_inplace( pApprox, iteration, res.appPrec, meta);
                    compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
                    restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
                    if (metadatas_getVerbo(meta)>=4)
                        printf("#deflate.c: deflate_tstar_test: discard: %d, %d-th coeff: %d, %i-th Graeffe it, prec: %ld\n", discard, res.nbOfSol, restemp, nbGraeffe, res.appPrec);
                    
                }
            }
        
        }
        iteration +=1;
            
    }
    
    connCmp_reu_set_comp( CC, compDsk_centerref( d ), compDsk_radiusref( d ),
                                  nbGraeffe, res.appPrec, pApprox );
    
    compApp_poly_clear(pApprox);
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

tstar_res deflate_tstar_test_rescale( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta) {
    
    clock_t start = clock();
    
    if (metadatas_getVerbo(meta)>=4) {
        compApp_t c;
        realApp_t r;
        compApp_init(c);
        realApp_init(r);
        compApp_set_compRat(c, compDsk_centerref(d), CCLUSTER_DEFAULT_PREC);
        realApp_set_realRat(r, compDsk_radiusref(d), CCLUSTER_DEFAULT_PREC);
        printf("\n#deflate.c: deflate_tstar_test_rescale: begin\n");
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
// //     int N = 0;
    int restemp = 0;
    int nbGraeffe = 0;
//     int iteration = 0;
    
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,deg+1);
    pApprox->length = deg+1;
    realApp_t sum;
    realApp_init(sum);
    
//     N = (int) 4+ceil(log2(1+log2(deg)));
    
    compApp_poly_set( pApprox, connCmp_reuPoref(CC) );
    res.appPrec = connCmp_reuPrref(CC);
    realRat_t ratio;
    realRat_init(ratio);
    realRat_set(ratio, compDsk_radiusref( d ));
    realRat_div(ratio, ratio, connCmp_reuRaref(CC));
    slong pow = 1 >> connCmp_reuNgref(CC);
    realRat_pow_si (ratio, ratio, pow);
    compApp_poly_scale_realRat_in_place( pApprox->coeffs, ratio, pApprox->length, res.appPrec );
    realRat_clear(ratio);
    
    compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
    
    while( (res.nbOfSol < max_nb_sols)&&(restemp==0)&&(res.nbOfSol<deg) ){
            res.nbOfSol += 1;
            
            restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
            
            if ( (restemp == -2)||(restemp == -1) )
                restemp = 0;
            
    }
    
//     if ((restemp==0)||(restemp==-1)||(restemp==-2)) res.nbOfSol = -1;
    if ((restemp==0)||(restemp==-2)) res.nbOfSol = -2;
    if ((restemp==-1)) res.nbOfSol = -1;
    
    
    compApp_poly_clear(pApprox);
    realApp_clear(sum);
    
    
    
    if (metadatas_haveToCount(meta))
//         if (discard)
            metadatas_add_time_DefTsta(meta, (double) (clock() - start) );
    
    if (metadatas_getVerbo(meta)>=4){
        printf("#deflate.c: deflate_tstar_test_rescale: discard: %d, res: %d, %i Graeffe its, prec: %ld\n", discard, res.nbOfSol, nbGraeffe, res.appPrec);
        printf("#deflate.c: deflate_tstar_test_rescale: end\n\n");
    }
    
    return res;
}

