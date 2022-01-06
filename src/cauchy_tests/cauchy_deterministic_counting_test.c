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

#include "cauchy_tests/cauchy_tests.h"

cauchyTest_res cauchyTest_deterministic_counting_combinatorial( const compRat_t center,
                                                                realRat_t radius,
                                                                slong nbOfRoots,
                                                                cacheCauchy_t cacheCau,
                                                                slong prec,
                                                                metadatas_t meta, int depth){
    
    
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    realApp_ptr wP  = cacheCauchy_wanErrExref(cacheCau);
    
    realRat_t gamma, rad;
    realRat_init(gamma);
    realRat_init(rad);
    realRat_mul(gamma, cacheCauchy_isoRatioref(cacheCau), cacheCauchy_isoRatioref(cacheCau));
    fmpz_add_si( realRat_numref(gamma), realRat_numref(gamma), 1);
    realRat_set(rad, radius);
    
    if (metadatas_getVerbo(meta)>=4) {
        printf("#---cauchy deterministic counting test in Newton iteration, combinatorial version: \n");
        printf("#------ isoRatio: "); realRat_print(cacheCauchy_isoRatioref(cacheCau)); printf("\n");
        printf("#------ max number of roots in disc: %ld\n", nbOfRoots);
        printf("#------ gamma: "); realRat_print(gamma); printf("\n");
    }
    
    slong i=0;
    res.nbOfSol = nbOfRoots;
    compApp_t s0;
    compApp_init(s0);
    
    while ((i < nbOfRoots+1) && (res.nbOfSol == nbOfRoots) ) {
        
        cacheCauchy_set_bounds( cacheCau, rad, CCLUSTER_DEFAULT_PREC);
        int alreadyEvaluated = 0;
        /* apply counting test to D(center, rad) */
        res = cauchyTest_computeS0Approx(s0, center, rad,
                                         NULL, 0, 0, &alreadyEvaluated, cacheCau, prec, CAUCHYTEST_INCOUNTIN, meta, depth );
    
        res.nbOfSol = ( res.nbOfSol==-2? -1:0 );
        if (res.nbOfSol==0) {
            realApp_add_error( compApp_realref(s0), wP + 0 );
            realApp_add_error( compApp_imagref(s0), wP + 0 );
            slong nbOfSol = -1;
            int unique = realApp_get_unique_si( &nbOfSol, compApp_realref(s0) );
            int containsZero = realApp_contains_zero( compApp_imagref(s0) );
            if (! ( unique && containsZero) ){
                res.nbOfSol = -1;
            }
            else {
                res.nbOfSol = (int) nbOfSol;
            }
        }
        
        if (metadatas_getVerbo(meta)>=4) {
            if ((res.nbOfSol == -1)||(res.nbOfSol != nbOfRoots)){
                printf("# ------%ld-th iteration, certification failed; res.nbOfSol: %d\n", i, res.nbOfSol);
            }
        }
        
        i++;
        realRat_div( rad, rad, gamma );
    }
    
    realRat_clear(gamma);
    realRat_clear(rad);
    compApp_clear(s0);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
        metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
    }
    
    return res;
    
}

slong cauchyTest_getNbEvals_counting_combinatorial_with_rinfrsup(  const compRat_t c,
                                                                   const realRat_t rinf,
                                                                   const realRat_t rsup,
                                                                   slong m,
                                                                   cacheCauchy_t cacheCau){
    
    slong res = 0;
    
    realRat_t shift, ri, thetai;
    realRat_init(shift);
    realRat_init(ri);
    realRat_init(thetai);
    /* set shift to (rsup - rinf)/(m+1) */
    realRat_sub(shift, rsup, rinf);
    realRat_div_ui(shift, shift, m+1);
    
    for (slong i = 0; i<m; i++) {
        /* set ri to rinf + (i+1)*shift */
        realRat_mul_si(ri, shift, i+1);
        realRat_add( ri, ri, rinf );
        /*set thetai to (ri + shift/2)/ri */
        realRat_div_ui(thetai, shift, 2);
        realRat_add(thetai, thetai, ri);
        realRat_div(thetai, thetai, ri);
        
        res = res + cacheCauchy_get_NbOfEvalPoints( cacheCauchy_degreeref(cacheCau), thetai, 1, CCLUSTER_DEFAULT_PREC );
        
    }
    
    realRat_clear(shift);
    realRat_clear(ri);
    realRat_clear(thetai);
    
    return res;
}

slong cauchyTest_getNbEvals_counting_combinatorial_with_isoRatio( const compRat_t c,
                                                                  const realRat_t r,
                                                                  const realRat_t theta,
                                                                  slong m,
                                                                  cacheCauchy_t cacheCau){
    
    realRat_t rinf, rsup;
    realRat_init(rinf);
    realRat_init(rsup);
    realRat_mul(rsup, r, theta);
    realRat_div(rinf, r, theta);
    
    slong res = cauchyTest_getNbEvals_counting_combinatorial_with_rinfrsup( c, rinf, rsup, m, cacheCau);
    
    realRat_clear(rinf);
    realRat_clear(rsup);
    
    return res;

}

cauchyTest_res cauchyTest_deterministic_counting_combinatorial_with_rinfrsup( const compRat_t c,
                                                                              const realRat_t rinf,
                                                                              const realRat_t rsup,
                                                                              slong m,
                                                                              cacheCauchy_t cacheCau,
                                                                              slong prec,
                                                                              metadatas_t meta, int depth){
    
    int level = 4;
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("------validation of a cluster with %ld roots\n", m);
        printf("---------rinf: "); realRat_print(rinf); printf("\n");
        printf("---------rsup: "); realRat_print(rsup); printf("\n");
    }
        
    realRat_t shift, ri, thetai;
    realRat_init(shift);
    realRat_init(ri);
    realRat_init(thetai);
    
    /* set shift to (rsup - rinf)/(m+1) */
    realRat_sub(shift, rsup, rinf);
    realRat_div_ui(shift, shift, m+1);
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = m;
    
    compDsk_t Delta;
    compDsk_init(Delta);
    compRat_set( compDsk_centerref(Delta), c);
    
    for (slong i = 0; (i<m)&&(res.nbOfSol==m); i++) {
        /* set ri to rinf + (i+1)*shift */
        realRat_mul_si(ri, shift, i+1);
        realRat_add( ri, ri, rinf );
        /*set thetai to (ri + shift/2)/ri */
        realRat_div_ui(thetai, shift, 2);
        realRat_add(thetai, thetai, ri);
        realRat_div(thetai, thetai, ri);
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("---------i=%ld; isoRatio: ", i); realRat_print(thetai);
            printf("\n");
        }
        
        /* call */
        realRat_set(compDsk_radiusref(Delta), ri);
        
        res = cauchyTest_probabilistic_counting_withIsoRatio( thetai, Delta, cacheCau, res.appPrec, meta, depth);
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("---------i=%ld; res.nbOfSol: %d\n", i, res.nbOfSol); 
        }
    }
    
    if (res.nbOfSol != m)
        res.nbOfSol=-1;
    
    realRat_clear(shift);
    realRat_clear(ri);
    realRat_clear(thetai);
    compDsk_clear(Delta);
    
    return res;
}

cauchyTest_res cauchyTest_deterministic_counting_combinatorial_with_isoRatio( const compRat_t c,
                                                                              const realRat_t r,
                                                                              const realRat_t theta,
                                                                              slong m,
                                                                              cacheCauchy_t cacheCau,
                                                                              slong prec,
                                                                              metadatas_t meta, int depth){
    
    realRat_t rinf, rsup;
    realRat_init(rinf);
    realRat_init(rsup);
    realRat_mul(rsup, r, theta);
    realRat_div(rinf, r, theta);
    
    cauchyTest_res res = cauchyTest_deterministic_counting_combinatorial_with_rinfrsup( c, rinf, rsup, m, cacheCau, prec,
                                                                                        meta, depth);
    
    realRat_clear(rinf);
    realRat_clear(rsup);
    
    return res;

}

int cauchyTest_shift( compApp_poly_t dest, const compRat_t c, const realRat_t r, const realRat_t theta, 
                               cacheCauchy_t cacheCau, slong prec, metadatas_t meta ){
    
    int level = 4;
    slong d = cacheCauchy_degreeref(cacheCau);
    /* number q of evaluation points */
    slong q = d+1;
    
    /* roots of unity and roots of unity shifted */
    compApp_ptr roots        = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
    compApp_ptr rootsShifted = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
    /* values of p(c+rx) at rootsShifted */ 
    compApp_ptr pvals        = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
    
    
    compApp_t center;
    compApp_init(center);
    compApp_set_compRat(center, c, prec);
    
    /* initialize and compute roots of unity */
    for(slong i=0; i<q; i++) {
        compApp_init (roots + i);
        compApp_init (rootsShifted + i);
        compApp_init (pvals + i);
        cauchyTest_computePointPointShifted( roots + i, rootsShifted + i, center, q, -i, r, prec );
    }
    
    realApp_t ubound, lbound;
    realApp_init( ubound );
    realApp_init( lbound );
    cauchyTest_computeBounds ( ubound, lbound, theta, r, d, CCLUSTER_DEFAULT_PREC );
    
    compApp_t pdiv, pder;
    compApp_init(pdiv);
    compApp_init(pder);
    
    int res = 1;
    slong nbEvals = 0;
    clock_t start = clock();
    for(nbEvals=0; (nbEvals<q) && (res==1); nbEvals++) {
        cacheCauchy_eval( pvals + nbEvals, pder, rootsShifted + nbEvals, 1, cacheCau, prec, meta );
        res = cauchyTest_compute_fdiv_checkPrecAndBounds( pdiv, pvals + nbEvals, pder, lbound, ubound, prec );
    }
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CertPEv(meta, (double) (clock() - start));
        metadatas_add_nbEvalsInPellet ( meta, nbEvals );
    }
    
    /* compute the coefficients of p(c+rx) */
    compApp_poly_init2(dest, q);
    
    start = clock();
//     acb_poly_interpolate_fast( dest, roots, pvals, q, prec );
    acb_dft_inverse(dest->coeffs, pvals, q, prec);
    dest->length = q;
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CertPIn(meta, (double) (clock() - start));
    }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("interpolated pol: ");
        compApp_poly_printd(dest, 10);
        printf("\n");
    
        
//         compApp_poly_t dest2;
//         compApp_poly_init2(dest2, q);
//         
//         compApp_poly_set(dest2, cacheApp_getApproximation ( cache, prec ));
//         compApp_poly_taylorShift_in_place( dest2, c, r, prec );
//         printf("shifted pol: ");
//         compApp_poly_printd(dest2, 10);
//         printf("\n");
//         compApp_poly_clear(dest2);
    }
    
    for(slong i=0; i<q; i++){
        compApp_clear (roots + i);
        compApp_clear (rootsShifted + i);
        compApp_clear (pvals + i);
    }
    ccluster_free(roots);
    ccluster_free(rootsShifted);
    ccluster_free(pvals);
    realApp_clear( ubound );
    realApp_clear( lbound );
    compApp_clear(pdiv);
    compApp_clear(pder);
    compApp_clear(center);
    
    return res;
}

cauchyTest_res cauchyTest_Pellet_counting( const compRat_t c, const realRat_t r, const realRat_t theta, slong m,
                                                                              cacheCauchy_t cacheCau,
                                                                              slong prec,
                                                                              metadatas_t meta, int depth){
    
    int level = 2;
    clock_t start = clock();
    
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = -1;
    
    int N = (int) 4+ceil(log2(1+log2(cacheCauchy_degreeref(cacheCau))));
    int iteration = 0;
    
    compApp_poly_t pShifted;
    compApp_poly_init(pShifted);
    realApp_t sum;
    realApp_init(sum);
    
    /* compute shift with evaluation-interpolation */
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---Pellet counting: m: %ld initial prec: %ld\n", m, res.appPrec);
    }
    res.nbOfSol = cauchyTest_shift( pShifted, c, r, theta, cacheCau, res.appPrec, meta );
    while ( res.nbOfSol == -1 ) {
        res.appPrec = 2*res.appPrec;
        res.nbOfSol = cauchyTest_shift( pShifted, c, r, theta, cacheCau, res.appPrec, meta );
    }
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---Pellet counting: res, prec after first shift: %d, %ld\n", res.nbOfSol, res.appPrec);
    }
    
    if (res.nbOfSol == 1) {
        
        res.nbOfSol = -1;
        while( (iteration <= N)&&(res.nbOfSol ==-1) ) {
            
            if (metadatas_getVerbo(meta)>=level) {
                printf("#------iteration: %d\n", iteration);
            }
    
            /* do iteration DLG iterations */
            clock_t start2 = clock();
            for (int i = 0; i<iteration; i++)
                compApp_poly_oneGraeffeIteration_in_place( pShifted, res.appPrec);
            if ((iteration>0)&&(metadatas_haveToCount(meta))) {
                metadatas_add_time_CertPGr(meta, (double) (clock() - start2));
                metadatas_add_nbGraeffeInPellet ( meta, iteration );
            }
    
            /* compare abs value of m-th coeff to the sum of abs value of other coeffs */
            compApp_poly_sum_abs_coeffs( sum, pShifted, res.appPrec );
            res.nbOfSol  = compApp_poly_TkGtilda_with_sum( pShifted, sum, m, res.appPrec);
            
            if (metadatas_getVerbo(meta)>=level) {
                printf("#------res.nbOfSol: %d\n", res.nbOfSol);
            }
            
            
            while( res.nbOfSol == -2 ){
                res.appPrec *=2;
                res.nbOfSol = cauchyTest_shift( pShifted, c, r, theta, cacheCau, res.appPrec, meta );
                for (int i = 0; i<iteration; i++)
                    compApp_poly_oneGraeffeIteration_in_place( pShifted, res.appPrec);
                compApp_poly_sum_abs_coeffs( sum, pShifted, res.appPrec );
                res.nbOfSol  = compApp_poly_TkGtilda_with_sum( pShifted, sum, m, res.appPrec);
                if (metadatas_getVerbo(meta)>=level) {
                    printf("#------res.appPrec: %ld, res.nbOfSol: %d\n", res.appPrec, res.nbOfSol);
                }
            }
            
            if (( res.nbOfSol == -1 ) || ( res.nbOfSol == 0 )) {
                res.nbOfSol = -1;
                iteration++;
            }
        }
    }
    
    if (res.nbOfSol == 1)
        res.nbOfSol = m;
    else
        res.nbOfSol = -2;
    
    realApp_clear(sum);
    compApp_poly_clear(pShifted);
    
    if (metadatas_haveToCount(meta)) {
        metadatas_add_time_CertPel(meta, (double) (clock() - start));
        metadatas_add_nbCertifiedWithPellet(meta, 1);
    }
    
    return res;
}



// int cauchyTest_shift_and_DLGs( compApp_poly_t dest, const compRat_t c, const realRat_t r, const realRat_t theta, cacheApp_t cache, 
//                                cacheCauchy_t cacheCau, slong prec, metadatas_t meta ){
//     
//     slong d = cacheCauchy_degreeref(cacheCau);
//     /* number N of DLG iterations */
// //     int N = (int) 4+ceil(log2(1+log2(d)));
//     int N = 1;
//     slong twototheN = 1<<N;
//     slong twototheNm1 = 1<<(N-1);
//     /* number q of evaluation points: smallest odd number >=d+1 */
//     /* so that q and 2^N are coprime */
//     slong q = ( ((d+1)%2)? d+1: d+2 );
//     
//     printf("q: %ld, N: %d, twototheN: %ld\n", q, N, twototheN);
//     
//     /* roots of unity and roots of unity shifted */
//     compApp_ptr roots        = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
//     compApp_ptr rootsShifted = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
//     /* values of p(c+rx) and (p(c+rx))^[N] at rootsShifted */ 
//     compApp_ptr pvals        = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
//     compApp_ptr pmvals       = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
//     compApp_ptr pNvals       = (compApp_ptr) ccluster_malloc ( q*sizeof(compApp) );
//     
//     compApp_t center;
//     compApp_init(center);
//     compApp_set_compRat(center, c, prec);
//     
//     /* initialize and compute roots of unity */
//     for(slong i=0; i<q; i++) {
//         compApp_init (roots + i);
//         compApp_init (rootsShifted + i);
//         compApp_init (pvals + i);
//         compApp_init (pmvals + i);
//         compApp_init (pNvals + i);
//         cauchyTest_computePointPointShifted( roots + i, rootsShifted + i, center, q, i, r, prec );
//     }
//     
//     realApp_t ubound, lbound;
//     realApp_init( ubound );
//     realApp_init( lbound );
//     cauchyTest_computeBounds ( ubound, lbound, theta, r, d, CCLUSTER_DEFAULT_PREC );
//     
//     compApp_t pdiv;
//     compApp_init(pdiv);
//     
//     int res = 1;
//     for(slong i=0; (i<q) && (res==1); i++) {
//         cauchyTest_eval ( pvals + i, pNvals + i, rootsShifted + i, cache, cacheCau, prec );
//         res = cauchyTest_compute_fdiv_checkPrecAndBounds( pdiv, pvals + i, pNvals + i, lbound, ubound, prec );
//         compApp_neg( pdiv, roots+i );
//         compApp_mul_realRat( pdiv, pdiv, r, prec );
//         compApp_add(pdiv, pdiv, center, prec);
//         cauchyTest_eval ( pmvals + i, pNvals + i, pdiv, cache, cacheCau, prec );
//         res = cauchyTest_compute_fdiv_checkPrecAndBounds( pdiv, pmvals + i, pNvals + i, lbound, ubound, prec );
//     }
//     
//     /* compute the coefficients of p(c+rx) */
//     compApp_poly_init2(dest, q);
//     acb_poly_interpolate_fast( dest, roots, pvals, q, prec );
//     printf("interpolated pol: ");
//     compApp_poly_printd(dest, 10);
//     printf("\n");
//     
//     compApp_poly_set(dest, cacheApp_getApproximation ( cache, prec ));
//     compApp_poly_taylorShift_in_place( dest, c, r, prec );
//     printf("shifted pol: ");
//     compApp_poly_printd(dest, 10);
//     printf("\n");
//     
//     /* compute the values of (p(c+rx))^[N] at roots of unity */
//     for(slong i=0; (i<q); i++) {
//         slong idest = (i*twototheN)%q;
//         printf("i: %ld, idest: %ld\n", i, idest);
//         compApp_mul( pNvals + idest, pvals + i, pmvals + i, prec );
//         compApp_pow_si( pNvals + idest, pNvals + idest, twototheNm1, prec );
//         /* mult by (-1)^(d) */
//         if (d%2==1)
//             compApp_neg( pNvals + idest, pNvals + idest );
//     }
//     
//     /* compute the coefficients of (p(c+rx))^[N] */
//     compApp_poly_init2(dest, q);
//     acb_poly_interpolate_fast( dest, roots, pNvals, q, prec );
//     printf("interpolated pol: ");
//     compApp_poly_printd(dest, 10);
//     printf("\n");
//     
//     compApp_poly_set(dest, cacheApp_getApproximation ( cache, prec ));
//     compApp_poly_taylorShift_in_place( dest, c, r, prec );
//     for(int i = 0; i < N; i++)       
//             compApp_poly_oneGraeffeIteration_in_place( dest, prec );
//     printf("shifted and iterated pol: ");
//     compApp_poly_printd(dest, 10);
//     printf("\n");
//     
//     for(slong i=0; i<q; i++){
//         compApp_clear (roots + i);
//         compApp_clear (rootsShifted + i);
//         compApp_clear (pvals + i);
//         compApp_clear (pmvals + i);
//         compApp_clear (pNvals + i);
//     }
//     ccluster_free(roots);
//     ccluster_free(rootsShifted);
//     realApp_clear( ubound );
//     realApp_clear( lbound );
//     compApp_clear(pdiv);
//     ccluster_free(pvals);
//     ccluster_free(pmvals);
//     ccluster_free(pNvals);
//     compApp_clear(center);
//     
//     return res;
// }

// void cauchyTest_schroeders_iteration_inplace ( compApp_t x, slong m, cacheApp_t cache, cacheCauchy_t cacheCau, slong prec){
//     compApp_t fval, fderval;
//     compApp_init(fval);
//     compApp_init(fderval);
//     
//     cauchyTest_eval ( fval, fderval, x, cache, cacheCau, prec );
//     compApp_div(fval, fval, fderval, prec);
//     compApp_mul_si(fval, fval, m, prec);
//     compApp_sub(x, x, fval, prec);
//     
//     compApp_clear(fval);
//     compApp_clear(fderval);
// }
// 
// cauchyTest_res cauchyTest_schroeders_counting( const compRat_t c, const realRat_t r, const realRat_t theta, slong m,
//                                                                               cacheApp_t cache,
//                                                                               cacheCauchy_t cacheCau,
//                                                                               slong prec,
//                                                                               metadatas_t meta, int depth){
//     cauchyTest_res res;
//     res.appPrec = prec;
//     res.nbOfSol = m;
//     
//     int level = 2;
//     slong d, q;
//     d = cacheCauchy_degreeref(cacheCau);
//     q = d+1;
//     realApp_t tau;
//     realApp_init(tau);
//     realApp_set_si(tau, q);
//     realApp_sqrt(tau, tau, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_si(tau, tau, (2*d + 2*m + 1)*d, CCLUSTER_DEFAULT_PREC);
//     realApp_inv(tau, tau, CCLUSTER_DEFAULT_PREC);
//     
//     if (metadatas_getVerbo(meta)>=level) {
//             compApp_t center;
//             realApp_t radius;
//             compApp_init(center);
//             realApp_init(radius);
//             compApp_set_compRat(center, c, CCLUSTER_DEFAULT_PREC);
//             realApp_set_realRat(radius, r, CCLUSTER_DEFAULT_PREC);
//             printf("#---------Schroeder's counting for a cluster with %ld roots: \n", m);
//             printf("#------------tau   : "); realApp_printd(tau, 10); printf("\n");
//             printf("#------------center: "); compApp_printd(center, 10); printf("\n");
//             printf("#------------radius: "); realApp_printd(radius, 10); printf("\n");
//             realApp_clear(radius);
//             compApp_clear(center);
//     }
//     
//     realApp_t ubound, lbound, absuval, absfderval;
//     realApp_init( ubound );
//     realApp_init( lbound );
//     realApp_init( absuval);
//     realApp_init(absfderval);
//     
//     compApp_t center, point, pointShifted, fval, fderval, fdiv, uval;
//     compApp_init(center);
//     compApp_init(point); 
//     compApp_init(pointShifted);
//     compApp_init(fval);
//     compApp_init(fderval);
//     compApp_init(fdiv);
//     compApp_init(uval);
//     
//     cauchyTest_computeBounds ( ubound, lbound, theta, r, d, CCLUSTER_DEFAULT_PREC );
//     
//     int nbSchroedersIterations = 3;
//     compRat_t initPoint;
//     compRat_init(initPoint);
//     compRat_set(initPoint, c);
//     realRat_add( compRat_realref(initPoint), compRat_realref(initPoint), r );
//     compApp_set_compRat(center, initPoint, res.appPrec);
//     /* apply Schroeder's iterations to center to re-center the disc */
//     for (int it = 0; it < nbSchroedersIterations; it++){
//         cauchyTest_schroeders_iteration_inplace (center, m, cache, cacheCau, 2*prec);
//     }
//     compRat_clear(initPoint);
//     
//     if (metadatas_getVerbo(meta)>=level) {
//         printf("#------------center after %d shroeder's iterations: ", nbSchroedersIterations); compApp_printd(center, 10); printf("\n");
//     }
//     
//     for(slong i=0; (i<q)&&(res.nbOfSol==m); i++) {
//         
//         int enoughPrec = -1;
//         while ((enoughPrec==-1)&&(res.nbOfSol==m)) {
//             
//             compApp_set_compRat(center, c, res.appPrec);
//             
//             cauchyTest_computePointPointShifted( point, pointShifted, center, q, i, r, res.appPrec );
//             cauchyTest_eval ( fval, fderval, pointShifted, cache, cacheCau, res.appPrec );
//             enoughPrec = cauchyTest_compute_fdiv_checkPrecAndBounds( fdiv, fval, fderval, lbound, ubound, res.appPrec );
//             
//             if (enoughPrec==1) {
//                 
//                 /* (p(c+rx))' = rp'(c+rx) */
//                 compApp_mul_realRat( fderval, fderval, r, res.appPrec );
//                 /* |u(x)| = |m*p(c+rx) - x*(p(c+rx))'| */
//                 compApp_mul( uval, point, fderval, res.appPrec );
//                 compApp_mul_si(fval, fval, m, res.appPrec);
//                 compApp_sub(uval, fval, uval, res.appPrec);
//                 compApp_abs(absuval, uval, res.appPrec);
//                 /* tau*|(p(c+rx))'| */
//                 compApp_abs(absfderval, fderval, res.appPrec);
//                 realApp_mul(absfderval, absfderval, tau, res.appPrec);
//                 
//                 /* compute shroeder's iteration */
//                 compApp_div(uval, uval, fderval, res.appPrec);
//                 printf("i: %ld, prec: %ld\n", i, res.appPrec);
//                 printf("absuval: "); realApp_printd(absuval, 10); printf("\n");
//                 printf("schroeder's first iteration: "); compApp_printd(uval, 10); printf("\n");
//                 printf("absfderval: "); realApp_printd(absfderval, 10); printf("\n");
//                 printf("\n");
//                 
//                 if ( realApp_gt( absuval, absfderval ) ) {
//                     printf("here: %d\n", realApp_gt( absuval, absfderval ));
//                     res.nbOfSol=-1;
//                 } else if ( !realApp_lt(absuval, absfderval) ){
//                     enoughPrec=-1;
//                 }
//             } 
//             
//             if (enoughPrec == -1) {
//                 res.appPrec = 2*res.appPrec;
//             } else if (enoughPrec==-2) { // the disk is not theta-isolated
//                 res.nbOfSol = -2;
//             }
//             
//         }
//         
//     }
//     
//     compApp_clear(uval);
//     compApp_clear(fdiv);
//     compApp_clear(fval);
//     compApp_clear(fderval);
//     compApp_clear(center);
//     compApp_clear(point); 
//     compApp_clear(pointShifted);
//     realApp_clear( absuval);
//     realApp_clear(absfderval);
//     realApp_clear( ubound );
//     realApp_clear( lbound );
//     realApp_clear(tau);
//     
//     return res;
// }

/* OLD version: */
// cauchyTest_res cauchyTest_deterministic_counting( const compDsk_t Delta,
//                                                   cacheApp_t cache,
//                                                   cacheCauchy_t cacheCau,
//                                                   slong prec,
//                                                   metadatas_t meta, int depth){
//     
//     realRat_t radInf, radSup;
//     realRat_init(radInf);
//     realRat_init(radSup);
//     
//     realRat_set(radInf, compDsk_radiusref(Delta));
//     realRat_mul_si(radInf, radInf, 3);
//     realRat_div_ui(radInf, radInf, 4);
//     realRat_set(radSup, compDsk_radiusref(Delta));
//     realRat_mul_si(radSup, radSup, 2);
//     
//     cauchyTest_res res = cauchyTest_deterministic_counting_with_radInf_radSup( compDsk_centerref(Delta), radInf, radSup, 
//                                                                                cache, cacheCau, prec,
//                                                                                meta, depth);
//     
//     realRat_clear(radInf);
//     realRat_clear(radSup);
//     
//     return res;
//                 
// }
// 
// cauchyTest_res cauchyTest_deterministic_counting_with_radInf_radSup( const compRat_t center,
//                                                                      const realRat_t radInf,  
//                                                                      const realRat_t radSup,
//                                                                      cacheApp_t cache,
//                                                                      cacheCauchy_t cacheCau,
//                                                                      slong prec,
//                                                                      metadatas_t meta, int depth){
//     clock_t start = clock();
//     
//     cauchyTest_res res;
//     res.appPrec = prec;
//     
//     realRat_t Rad, ExRad, ratio, isoRatio;
//     realApp_t nbCentersApp;
//     realRat_init(Rad);
//     realRat_init(ExRad);
//     realRat_init(ratio);
//     realApp_init(nbCentersApp);
//     realRat_init(isoRatio);
//     
//     realRat_add(Rad, radInf, radSup);
//     realRat_div_ui(Rad, Rad, 2);
//     realRat_sub(ExRad, radSup, radInf);
//     realRat_div_ui(ExRad, ExRad, 2);
// //     realRat_div(ExRad, ExRad, cacheCauchy_isoRatioref(cacheCau));
//     realRat_div(ratio, Rad, ExRad);
//     realApp_pi(nbCentersApp, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_si(nbCentersApp, nbCentersApp, 2, CCLUSTER_DEFAULT_PREC);
//     realApp_mul_realRat(nbCentersApp, nbCentersApp, ratio, CCLUSTER_DEFAULT_PREC);
//     slong nbCenters = realApp_ceil_si(nbCentersApp, CCLUSTER_DEFAULT_PREC);
//     /* */
//     nbCenters = CCLUSTER_MAX( 12, nbCenters );
//     
//     /* isoRatio = (4/5)*(1+ExRad/Rad) = (4/5)*(1+1/ratio) */
//     realRat_inv(isoRatio, ratio);
//     realRat_add_si( isoRatio, isoRatio, 1);
//     realRat_mul_si( isoRatio, isoRatio, 4);
//     realRat_div_ui( isoRatio, isoRatio, 5);
//     
// //     cacheCauchy_set_bounds( cacheCau, Rad, CCLUSTER_DEFAULT_PREC);
// //     cacheCauchy_set_bounds( cacheCau, radInf, CCLUSTER_DEFAULT_PREC);
//     
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("#---cauchy deterministic counting test: \n");
//         printf("#------ isoRatio for counting : "); realRat_print( isoRatio); printf("\n");
//         printf("#------ isoRatio for exclusion: "); realRat_print( cacheCauchy_isoRatioref(cacheCau) ); printf("\n");
//         printf("#------ nb of discs on contour for certification: %d \n", (int) nbCenters);
//         printf("#------ radius of discs on contour for certification: "); realRat_print( ExRad ); printf("\n");
// //         printf("#------ number of roots in disc: %ld\n", nbOfRoots);
//     }
//     
//     /* try to prove that D(c,Rad) is isoRatio-isolated */
//     cauchyTest_res resEx;
//     resEx.nbOfSol = 0;
//     resEx.appPrec = res.appPrec;
//     for (int vindex = 0; vindex < nbCenters && (resEx.nbOfSol==0) ; vindex++){
//         resEx = cauchyTest_deterministic_exclusion_test( center, ExRad, Rad, nbCenters, vindex, cache, cacheCau, resEx.appPrec, CAUCHYTEST_INCOUNTIN, meta, depth );
//     }
//     
//     if (resEx.nbOfSol != 0) {
//         res.nbOfSol = -1;
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#------ certification failed!\n");
//         }
//     } else {
//     /* D(c,Rad) is isoRatio-isolated */
//         compDsk_t Delta;
//         compDsk_init(Delta);
//         compDsk_set_compRat_realRat(Delta, center, Rad);
//         slong nbOfRootsC = cauchyTest_computeS0compDsk( isoRatio, Delta, cache, cacheCau, meta, depth);
//         compDsk_clear(Delta);
//         res.nbOfSol = nbOfRootsC;
//     }
//     
//     realRat_clear(isoRatio);
//     realRat_clear(Rad);
//     realRat_clear(ExRad);
//     realRat_clear(ratio);
//     realApp_clear(nbCentersApp);
//     
//     if (metadatas_haveToCount(meta)) {
//         metadatas_add_time_CauCoTo(meta, (double) (clock() - start));
//         metadatas_add_CauchyCoTest(meta, depth, res.appPrec);
//     }
//     
//     return res;
// }

