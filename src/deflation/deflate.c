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

void deflate_connCmp_init  (connCmp_t x){
    compApp_poly_init( connCmp_defPoref(x) );
    realApp_init     ( connCmp_sumAbref(x) );
}

void deflate_connCmp_clear (connCmp_t x){
    compApp_poly_clear( connCmp_defPoref(x) );
    realApp_clear     ( connCmp_sumAbref(x) );
    connCmp_isDefref(x) = 0;
    connCmp_degDeref(x) = 0;
}

void deflate_taylor_shift_interval_inplace( compApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_taylorShift_interval_in_place( res, compDsk_centerref(d), compDsk_radiusref(d), prec );

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefTayl(meta, (double) (clock() - start) );
}

void deflate_derivative_inplace( compApp_poly_t res, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_derivative(res, res, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefDeri(meta, (double) (clock() - start) );
}

void deflate_evaluate(compApp_t y, compApp_poly_t f, const compApp_t x, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        
        compApp_poly_evaluate(y, f, x, prec);

        if (metadatas_haveToCount(meta))
            metadatas_add_time_DefEval(meta, (double) (clock() - start) );
}

void deflate_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta ){
    
    connCmp_isDefref(x) = 1;
    connCmp_degDeref(x) = nbSols;
    /* compute the interval taylor shift at any point in the disk */
    compApp_poly_set(connCmp_defPoref(x), cacheApp_getApproximation ( cache, prec ));
//     compApp_poly_taylorShift_interval_in_place( connCmp_defPoref(x), compDsk_centerref(disk), compDsk_radiusref(disk), prec );
    deflate_taylor_shift_interval_inplace( connCmp_defPoref(x), disk, prec, meta);
    
    /* compute the sum of the absolute values of the degree - nbSols leading coeffs */
    realApp_t abs;
    realApp_init(abs);
    compApp_abs( connCmp_sumAbref(x), 
                 connCmp_defPoref(x)->coeffs + nbSols + 1, prec );
    for (slong index = connCmp_nSolsref(x) + 2; index < connCmp_defPoref(x)->length; index++ ){
        compApp_abs( abs, connCmp_defPoref(x)->coeffs + index, prec );
        realApp_add( connCmp_sumAbref(x), connCmp_sumAbref(x), abs, prec );
    }
    
    realApp_clear(abs);
    
}

void deflate_copy( connCmp_t dest, const connCmp_t src ){
    if (connCmp_isDefref(src) != 0) {
        connCmp_isDefref(dest) = connCmp_isDefref(src);
        connCmp_degDeref(dest) = connCmp_degDeref(src);
        compApp_poly_set( connCmp_defPoref(dest), connCmp_defPoref(src) );
        realApp_set     ( connCmp_sumAbref(dest), connCmp_sumAbref(src) );
    }
}

void _deflate_compute_trailing_coeffs(realApp_ptr coeffs, cacheApp_t cache, const compDsk_t d, int nbCoeffs, slong prec, metadatas_t meta){
    
    /* pol */
    compApp_poly_t fapprox;
    compApp_poly_init(fapprox);
    compApp_poly_set(fapprox, cacheApp_getApproximation ( cache, prec ));
    
    compApp_t center, coeff;
    realRat_t factor;
    compApp_init(center);
    compApp_init(coeff);
    realRat_init(factor);
    
    realRat_set_si(factor, 1, 1);
    compApp_set_compRat(center, compDsk_centerref(d), prec);
    
    for (int index=0; index<nbCoeffs; index ++) {
        realApp_init( coeffs + index );
//         compApp_poly_evaluate(coeff, fapprox, center, prec);
        deflate_evaluate(coeff, fapprox, center, prec, meta);
        compApp_mul_realRat(coeff, coeff, factor, prec);
        compApp_abs( coeffs + index, coeff, prec);
        if (index<nbCoeffs){
            realRat_mul(factor, factor, compDsk_radiusref(d));
            realRat_div_ui(factor, factor, (ulong) (index+1));
//             compApp_poly_derivative(fapprox, fapprox, prec);
            deflate_derivative_inplace( fapprox, prec, meta);
        }
    }
    
    compApp_poly_clear(fapprox);
    compApp_clear(center);
    compApp_clear(coeff);
    realRat_clear(factor);
}

void _deflate_compute_sum_leading(realApp_t sum, const connCmp_t CC, const compDsk_t d, slong prec, metadatas_t meta){
    realRat_t factor;
    realRat_init(factor);
    
    realRat_set(factor, compDsk_radiusref(d));
    realRat_pow_si( factor, factor, connCmp_degDeref(CC) + 1);
    realApp_mul_realRat( sum, connCmp_sumAbref(CC), factor, prec );
    
    realRat_clear(factor);
}

void _deflate_compute_sum_leading_with_scaling(realApp_t sum, const connCmp_t CC, const compDsk_t d, slong prec, metadatas_t meta){
//     realRat_t factor;
    realApp_t factor, temp;
    compApp_t coeff;
    realApp_t abs;
    
//     realRat_init(factor);
    realApp_init(factor);
    realApp_init(temp);
    compApp_init(coeff);
    realApp_init(abs);
    
    clock_t start = clock();
    
//     realRat_set(factor, compDsk_radiusref(d));
//     realRat_pow_si( factor, factor, connCmp_degDeref(CC) + 1);
//     compApp_mul_realRat( coeff, connCmp_defPoref(CC)->coeffs + connCmp_degDeref(CC) + 1, factor, prec );
//     compApp_abs( sum, coeff, prec );
    realApp_set_realRat(temp, compDsk_radiusref(d), prec);
    realApp_pow_ui( factor, temp, connCmp_degDeref(CC) + 1, prec);
    compApp_abs( sum, connCmp_defPoref(CC)->coeffs + connCmp_degDeref(CC) + 1, prec);
    realApp_mul( sum, sum, factor, prec);
    for (int index = connCmp_degDeref(CC) + 2; index < connCmp_defPoref(CC)->length; index ++ ) {
//         realRat_mul(factor, factor, compDsk_radiusref(d));
//         compApp_mul_realRat( coeff, connCmp_defPoref(CC)->coeffs + index, factor, prec );
//         compApp_abs( abs, coeff, prec );
//         realApp_add(sum, sum, abs, prec);
        realApp_mul(factor, factor, temp, prec);
        compApp_abs( abs, connCmp_defPoref(CC)->coeffs + connCmp_degDeref(CC) + 1, prec);
        realApp_mul( abs, abs, factor, prec );
        realApp_add(sum, sum, abs, prec);
    }
    
    if (metadatas_haveToCount(meta))
            metadatas_add_time_DefScal(meta, (double) (clock() - start) );
    
//     realRat_clear(factor);
    realApp_clear(factor);
    realApp_clear(temp);
    compApp_clear(coeff);
    realApp_clear(abs);
}

tstar_res deflate_tstar_test( const connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, slong prec, metadatas_t meta) {
    
    tstar_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    
    realApp_ptr coeffs;
    realApp_t sum;
    realApp_init(sum);
    
    coeffs = (realApp_ptr) ccluster_malloc( (connCmp_degDeref(CC) +1)*sizeof(realApp) );
    for (int index=0; index<=connCmp_degDeref(CC); index ++)
        realApp_init( coeffs + index );
    
    /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
    _deflate_compute_trailing_coeffs(coeffs, cache, d, connCmp_degDeref(CC) +1, res.appPrec, meta);
    
    /*compute the sum of deg leading coefficients */
    /*scale the sum*/
    _deflate_compute_sum_leading(sum, CC, d, res.appPrec, meta);
//     _deflate_compute_sum_leading_with_scaling(sum, CC, d, res.appPrec, meta);
    for (int index = 1; index<=connCmp_degDeref(CC); index ++)
        realApp_add( sum, sum, coeffs + index, res.appPrec );
    
    int ind = 0;
    int resPellet=0;
    while ( ( ind <= connCmp_degDeref(CC) ) 
         && ( resPellet != 1 ) ) {
        resPellet = realApp_soft_compare(coeffs + ind, sum, res.appPrec);
        printf("Pellet test, %d-th coeff: %d\n", ind, resPellet);
//         if (resPellet == -2) {
//             printf("---%d-th coeff: ", ind);
//             realApp_printd(coeffs + ind, res.appPrec);
//             printf("\n---      sum  : ");
//             realApp_printd(sum, res.appPrec);
//             printf("\n");
//         }
        if (ind<connCmp_degDeref(CC)){
            realApp_add( sum, sum, coeffs + ind     , res.appPrec );
            realApp_sub( sum, sum, coeffs + (ind +1), res.appPrec );
        }
        ind++;
    }
    
    if (resPellet!=1) {
        /* try to double working precision */
        res.appPrec = 2*res.appPrec;
        /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
        _deflate_compute_trailing_coeffs(coeffs, cache, d, connCmp_degDeref(CC) +1, res.appPrec, meta);
        /*compute the sum of deg leading coefficients */
        /*scale the sum*/
//         _deflate_compute_sum_leading(sum, CC, d, res.appPrec, meta);
        _deflate_compute_sum_leading_with_scaling(sum, CC, d, res.appPrec, meta);
        for (int index = 1; index<=connCmp_degDeref(CC); index ++)
            realApp_add( sum, sum, coeffs + index, res.appPrec );
        
        ind = 0;
        resPellet=0;
        while ( ( ind <= connCmp_degDeref(CC) ) 
            && ( resPellet != 1 ) ) {
            resPellet = realApp_soft_compare(coeffs + ind, sum, res.appPrec);
            printf("Pellet test, %d-th coeff: %d\n", ind, resPellet);
//             if (resPellet == -2) {
//                 printf("---%d-th coeff: ", ind);
//                 realApp_printd(coeffs + ind, res.appPrec);
//                 printf("\n---      sum  : ");
//                 realApp_printd(sum, res.appPrec);
//                 printf("\n");
//             }
            if (ind<connCmp_degDeref(CC)){
                realApp_add( sum, sum, coeffs + ind     , res.appPrec );
                realApp_sub( sum, sum, coeffs + (ind +1), res.appPrec );
            }
            ind++;
        }
        
    }
//         else { /* try to decrease the precision */
//             if (res.appPrec > CCLUSTER_DEFAULT_PREC)
//                 res.appPrec = res.appPrec/2;
//         }
    
//     if (resPellet==-2) {
//         /* try to double working precision */
//         res.appPrec = 2*res.appPrec;
//         /* compute the connCmp_degDeref(CC) +1  trailing coefficients */
//         _deflate_compute_trailing_coeffs(coeffs, cache, d, connCmp_degDeref(CC) +1, res.appPrec, meta);
//         /*compute the sum of deg leading coefficients with better scaling */
//         /*scale the sum*/
//         _deflate_compute_sum_leading_with_scaling(sum, CC, d, res.appPrec, meta);
//         for (int index = 1; index<=connCmp_degDeref(CC); index ++)
//             realApp_add( sum, sum, coeffs + index, res.appPrec );
//         
//         ind = 0;
//         resPellet=0;
//         while ( ( ind <= connCmp_degDeref(CC) ) 
//             && ( resPellet != 1 ) ) {
//             resPellet = realApp_soft_compare(coeffs + ind, sum, res.appPrec);
//             printf("Pellet test, %d-th coeff: %d\n", ind, resPellet);
// //             if (resPellet == -2) {
// //                 printf("---%d-th coeff: ", ind);
// //                 realApp_printd(coeffs + ind, res.appPrec);
// //                 printf("\n---      sum  : ");
// //                 realApp_printd(sum, res.appPrec);
// //                 printf("\n");
// //             }
//             if (ind<connCmp_degDeref(CC)){
//                 realApp_add( sum, sum, coeffs + ind     , res.appPrec );
//                 realApp_sub( sum, sum, coeffs + (ind +1), res.appPrec );
//             }
//             ind++;
//         }
//         
//     }
    
    if (resPellet==1)
        res.nbOfSol = ind -1;
    else 
        res.nbOfSol = -2;
    
    realApp_clear(sum);
    return res;
}
