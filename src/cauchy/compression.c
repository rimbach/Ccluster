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

#include "cauchy/cauchy.h"

/* let res = D(c',r'), Delta=D(c,r) */
/* set r' st |c-c'|+r/theta <= r' <= 5/4 r */
void cauchy_setRadSup( compDsk_t res, const compDsk_t Delta, const realRat_t theta ){
    
    compRat_t temp;
    compRat_init(temp);
    compRat_sub(temp, compDsk_centerref(Delta), compDsk_centerref(res) );
    
    compApp_t tempApp;
    compApp_init(tempApp);
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    compApp_set_compRat(tempApp, temp, prec);
    realApp_t distApp;
    realApp_init(distApp);
    compApp_abs( distApp, tempApp, prec );
    
    realRat_t dist, rOtheta, ubound;
    realRat_init(dist);
    realRat_init(rOtheta);
    realRat_init(ubound);
    
    realApp_get_realRat(dist, distApp);
    realRat_div( rOtheta, compDsk_radiusref(Delta), theta );
    realRat_set_si(ubound, 5,4);
    realRat_mul(ubound, ubound, compDsk_radiusref(Delta) );
    
    realRat_add(compDsk_radiusref(res), dist, rOtheta );
    
    while( realRat_cmp( compDsk_radiusref(res), ubound ) >0 ) {
        prec = 2*prec;
        compApp_set_compRat(tempApp, temp, prec);
        compApp_abs( distApp, tempApp, prec );
        realApp_get_realRat(dist, distApp);
        realRat_add(compDsk_radiusref(res), dist, rOtheta );
    }
    
    compRat_clear(temp);
    compApp_clear(tempApp);
    realApp_clear(distApp);
    realRat_clear(dist);
    realRat_init(rOtheta);
    realRat_init(ubound);
}

/* Let Delta = D(c,r) have m roots and be theta-isolated */
// slong cauchy_rootBound_deflation( realApp_t lbound,
//                                   realApp_t ubound,
//                                   const compDsk_t Delta,
//                                   const realRat_t theta,
//                                   slong m,
//                                   const compRat_t cprime,
//                                   cacheApp_t cache,
//                                   cacheCauchy_t cacheCau,
//                                   slong prec, metadatas_t meta, slong depth) {
//     
//     int level = 3;
//     
//     /* number of DLG iterations to apply */
//     int N = (int) ceil(log2(1+log2(m+1)));
//     int n = 0x1<<N;
//     if (metadatas_getVerbo(meta)>=level) {
//             printf("#---------compression.c: cauchy_rootBound_deflation; N: %d\n", N);
//     }
//     /* want Nth DLG iteration of factor with precision about prec2 */
//     /* compute power sums at precision 2^N*prec2 */
//     slong prec2 = CCLUSTER_DEFAULT_PREC;
// //     prec2 = N*m;
//     prec2 = prec2*N;
//     realRat_t eps;
//     realRat_init(eps);
//     realRat_set_si(eps, 1, 2);
//     realRat_pow_si(eps, eps, prec2);
//     prec2 = CCLUSTER_MAX(prec2, CCLUSTER_DEFAULT_PREC);
//     
//     compApp_ptr Shs = (compApp_ptr) ccluster_malloc ( m*sizeof(compApp) );
//     for (slong i=0; i<m; i++)
//         compApp_init( Shs + i );
//         
//     /* compute m power sums s1, s2, ..., sm of roots of p(c + rz) in D(0,1) at precision eps */
//     cauchyTest_res resT = cauchyTest_computeSgNcompDsk( Shs, theta, Delta, m, 1, m, cache, cacheCau, eps, prec2, meta, depth);
//     
//     realRat_clear(eps);
//     
//     if ( resT.nbOfSol == -2 ) {
//         if (metadatas_getVerbo(meta)>=level) {
//             printf("#---------compression.c: cauchy_rootBound_deflation; input disk is not theta-isolated\n");
//         }
//         realRat_clear(eps);
//         for (slong i=0; i<m; i++)
//             compApp_clear( Shs + i );
//         ccluster_free(Shs);
//         return -2;
//     }
//     
//     /* get monic polynomial f_{c,r}(z)=f(c+rz) of degree m+1 having power sums s1, s2, ..., sm in D(0,1)*/
//     compApp_poly_t factor;
//     compApp_poly_init2(factor, m+1);
//     /* apply Newton Identities to get monic pol*/
//     compApp_ptr coeffs = factor->coeffs;
//     factor->length = m+1;
//     compApp_one( coeffs + m );
//     for (slong j = m-1; j>=0; j-- ){
//         compApp_set( coeffs + j, Shs + (m-j) -1 );
//         for (slong k = 1; k< (m-j); k++)
//             compApp_addmul( coeffs + j,  Shs + k -1, coeffs + (j+k), resT.appPrec );
//         compApp_div_si( coeffs + j, coeffs + j, -(m-j), resT.appPrec  );
//     }
//     
//     for (slong i=0; i<m; i++)
//         compApp_clear( Shs + i );
//     ccluster_free(Shs);
//     
//     /* scale by 1/r => get monic f_{c,1}(z) = f(c + z) */
//     realRat_t rinv;
//     realRat_init(rinv);
//     realRat_inv( rinv, compDsk_radiusref(Delta) );
//     compApp_poly_scale_realRat_in_place_monic( factor->coeffs, rinv, factor->length, resT.appPrec );
//     realRat_clear(rinv);
//     
//     /* shift in cprime - c => get monic f_{cprime,1}(z) = f(cprime + z) */
//     compRat_t nc;
//     compRat_init(nc);
//     compRat_sub(nc, cprime, compDsk_centerref(Delta) );
//     compApp_poly_taylorShift_in_place_noscale( factor, nc, resT.appPrec );
//     compRat_clear(nc);
//     
//     /* apply N = ceil ( log_2 ( 1 + log_2 (m+1) ) ) DLG iterations */
//     for (int i=0; i<N; i++)
//         compApp_poly_oneGraeffeIteration_in_place( factor, resT.appPrec);
//     /* factor is still monic */
//     /* compute root bound for root of factor with greatest modulus */
//     compApp_poly_monic_bound_r1( lbound, ubound, factor, resT.appPrec);
//     realApp_pos_root_ui(lbound,  lbound, n, prec);
//     realApp_pos_root_ui(ubound,  ubound, n, prec);
//     if (metadatas_getVerbo(meta)>=level) {
//             printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; computed bounds for r1:\n");
//             realApp_printd( lbound, 10 ); printf("\n");
//             realApp_printd( ubound, 10 ); printf("\n");
//             printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; precision required: %ld\n", resT.appPrec);
//     }
//     
//     compApp_poly_clear(factor);
//     return resT.appPrec;
// }

//      (1/2m)^(1/2^g) >= 2/3
// <=>  1/2m >= (2/3)^(2^g)
// <=>  2m   <= (3/2)^(2^g)
// <=>  log(2m) <= 2^g log(3/2)
// <=>  loglog(2m) <= g log(3/2)
// <=>  g >= loglog(2m)/log(3/2)
/* Let Delta = D(c,r) have m roots and be theta-isolated */
slong cauchy_rootBound_deflation( realApp_t lbound,
                                   realApp_t ubound,
                                   const realRat_t relativeError,
                                   const compDsk_t Delta,
                                   const realRat_t theta,
                                   slong m,
                                   const compRat_t cprime,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   slong prec, metadatas_t meta, slong depth) {
    
    int level = 3;
    
    /* number of DLG iterations to apply */
    int N = (int) ceil( log2(1+log2(m+1))/log2(1.75) );
    int n = 0x1<<N;
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---------compression.c: cauchy_rootBound_deflation; N: %d\n", N);
    }
//     /* want Nth DLG iteration of factor with precision about prec2 */
//     /* compute power sums at precision 2^N*prec2 */
    slong prec2 = CCLUSTER_DEFAULT_PREC;
    prec2 = N*m;
    prec2 = CCLUSTER_MAX(prec2, CCLUSTER_DEFAULT_PREC);
//     prec2 = prec2*N;
    realRat_t eps, epsp;
    realRat_init(eps);
    realRat_init(epsp);
    realRat_set_si(eps, 1, 2);
    realRat_pow_si(eps, eps, prec2);
    realRat_set(epsp, eps);
//     prec2 = CCLUSTER_MAX(prec2, CCLUSTER_DEFAULT_PREC);
//     
    compApp_ptr Shs = (compApp_ptr) ccluster_malloc ( m*sizeof(compApp) );
    for (slong i=0; i<m; i++)
        compApp_init( Shs + i );
    
    int enoughPrec = -1;
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = enoughPrec;
    
    compApp_poly_t factor;
    compApp_poly_init2(factor, m+1);
    realRat_t rinv;
    realRat_init(rinv);
    realRat_inv( rinv, compDsk_radiusref(Delta) );
    compRat_t nc;
    compRat_init(nc);
    compRat_sub(nc, cprime, compDsk_centerref(Delta) );
    compApp_t temp;
    compApp_init(temp);
    realApp_t mid, rad;
    realApp_init(mid);
    realApp_init(rad);
    
    while ( enoughPrec == -1 ) {
        
        /* compute m power sums s1, s2, ..., sm of roots of p(c + rz) in D(0,1) at precision eps */
        res = cauchyTest_computeSgNcompDsk( Shs, theta, Delta, m, 1, m, cache, cacheCau, epsp, prec2, meta, depth);
    
        if ( res.nbOfSol == -2 ) {
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------compression.c: cauchy_rootBound_deflation; input disk is not theta-isolated\n");
            }
            break;
        }
        
        /* get monic polynomial f_{c,r}(z)=f(c+rz) of degree m+1 having power sums s1, s2, ..., sm in D(0,1)*/
        /* apply Newton Identities to get monic pol*/
        compApp_ptr coeffs = factor->coeffs;
        factor->length = m+1;
        compApp_one( coeffs + m );
        for (slong j = m-1; j>=0; j-- ){
            compApp_set( coeffs + j, Shs + (m-j) -1 );
            for (slong k = 1; k< (m-j); k++)
                compApp_addmul( coeffs + j,  Shs + k -1, coeffs + (j+k), res.appPrec );
            compApp_div_si( coeffs + j, coeffs + j, -(m-j), res.appPrec  );
        }
        /* scale by 1/r => get monic f_{c,1}(z) = f(c + z) */
        compApp_poly_scale_realRat_in_place_monic( factor->coeffs, rinv, factor->length, res.appPrec );
        /* shift in cprime - c => get monic f_{cprime,1}(z) = f(cprime + z) */
        compApp_poly_taylorShift_in_place_noscale( factor, nc, res.appPrec );
        /* apply N = ceil ( log_2 ( log_2 (5) ) ) DLG iterations */
        for (int i=0; i<N; i++)
            compApp_poly_oneGraeffeIteration_in_place( factor, res.appPrec);
        
        /* factor is still monic */
        
        /* compute root bound for root of factor with greatest modulus */
        compApp_poly_monic_bound_r1( lbound, ubound, factor, res.appPrec);
        realApp_pos_root_ui(lbound,  lbound, n, prec);
        realApp_pos_root_ui(ubound,  ubound, n, prec);
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------compression.c: cauchy_rootBound_deflation; computed bounds for r1:\n");
            realApp_printd( lbound, 10 ); printf("\n");
            realApp_printd( ubound, 10 ); printf("\n");
            printf("#---------compression.c: cauchy_rootBound_deflation; precision required: %ld\n", res.appPrec);
        }
        
        realApp_mul_realRat(compApp_realref(temp), lbound, relativeError, res.appPrec);
        /* check if lbound*relativeError >= ubound */
        if ( realApp_lt( compApp_realref(temp), ubound ) == 1 ) {
            enoughPrec=-1;
            realRat_pow_si(epsp, epsp, 2);
            prec2 = 2*prec2;
        } else 
            enoughPrec=1;
    }
        
    realRat_clear(eps);    
    for (slong i=0; i<m; i++)
        compApp_clear( Shs + i ); 
    ccluster_free(Shs);
    compApp_poly_clear(factor);
    realRat_clear(rinv);
    compRat_clear(nc);
    compApp_clear(temp);
    realApp_clear(mid);
    realApp_clear(rad);

    return res.appPrec;

}

/* Let Delta = D(c,r) have m roots and be theta-isolated */
slong cauchy_rootBound_deflation_Turan( realRat_t lbound,
                                        realRat_t ubound,
                                        const realRat_t relativeError,
                                        const compDsk_t Delta,
                                        const realRat_t theta,
                                        slong m,
                                        const compRat_t cprime,
                                        const realRat_t eps,
                                        cacheApp_t cache,
                                        cacheCauchy_t cacheCau,
                                        slong prec, metadatas_t meta, slong depth) {
    
    int level = 3;
    
    /* number of DLG iterations to apply */
    int N = (int) ceil(log2(log2(5))) + 1;
    int n = 0x1<<N;
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---------compression.c: cauchy_rootBound_deflation2; N: %d\n", N);
    }
//     /* want Nth DLG iteration of factor with precision about prec2 */
//     /* compute power sums at precision 2^N*prec2 */
    slong prec2 = CCLUSTER_DEFAULT_PREC;
    prec2 = N*m;
    prec2 = CCLUSTER_MAX(prec2, CCLUSTER_DEFAULT_PREC);
//     prec2 = prec2*N;
//     realRat_t eps, epsp;
    realRat_t epsp;
//     realRat_init(eps);
    realRat_init(epsp);
    realRat_set_si(epsp, 1, 2);
    realRat_pow_si(epsp, epsp, prec2);
//     realRat_set(epsp, eps);
//     prec2 = CCLUSTER_MAX(prec2, CCLUSTER_DEFAULT_PREC);
//     
    compApp_ptr Shs = (compApp_ptr) ccluster_malloc ( m*sizeof(compApp) );
    for (slong i=0; i<m; i++)
        compApp_init( Shs + i );
    
    int enoughPrec = -1;
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = enoughPrec;
    
    compApp_poly_t factor;
    compApp_poly_init2(factor, m+1);
    realRat_t rinv;
    realRat_init(rinv);
    realRat_inv( rinv, compDsk_radiusref(Delta) );
    compRat_t nc;
    compRat_init(nc);
    compRat_sub(nc, cprime, compDsk_centerref(Delta) );
    compApp_t temp;
    compApp_init(temp);
    realApp_t r;
    realApp_init(r);
    realRat_t t, exp;
    realRat_init(t);
    realRat_init(exp);
//     realApp_t mid, rad;
//     realApp_init(mid);
//     realApp_init(rad);
    
    while ( enoughPrec == -1 ) {
        
        /* compute m power sums s1, s2, ..., sm of roots of p(c + rz) in D(0,1) at precision eps */
        res = cauchyTest_computeSgNcompDsk( Shs, theta, Delta, m, 1, m, cache, cacheCau, epsp, prec2, meta, depth);
    
        if ( res.nbOfSol == -2 ) {
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------compression.c: cauchy_rootBound_deflation2; input disk is not theta-isolated\n");
            }
            break;
        }
        
        /* get monic polynomial f_{c,r}(z)=f(c+rz) of degree m+1 having power sums s1, s2, ..., sm in D(0,1)*/
        /* apply Newton Identities to get monic pol*/
        compApp_ptr coeffs = factor->coeffs;
        factor->length = m+1;
        compApp_one( coeffs + m );
        for (slong j = m-1; j>=0; j-- ){
            compApp_set( coeffs + j, Shs + (m-j) -1 );
            for (slong k = 1; k< (m-j); k++)
                compApp_addmul( coeffs + j,  Shs + k -1, coeffs + (j+k), res.appPrec );
            compApp_div_si( coeffs + j, coeffs + j, -(m-j), res.appPrec  );
        }
        /* scale by 1/r => get monic f_{c,1}(z) = f(c + z) */
        compApp_poly_scale_realRat_in_place_monic( factor->coeffs, rinv, factor->length, res.appPrec );
        /* shift in cprime - c => get monic f_{cprime,1}(z) = f(cprime + z) */
        compApp_poly_taylorShift_in_place_noscale( factor, nc, res.appPrec );
        /* apply N = ceil ( log_2 ( log_2 (5) ) ) DLG iterations */
        for (int i=0; i<N; i++)
            compApp_poly_oneGraeffeIteration_in_place( factor, res.appPrec);
        /* factor is still monic */
        /* apply newton identities to get the power sums s1, s2, ..., sm of the factor */  
        for (slong i=1; i<=m; i++) {
            compApp_mul_si( Shs + (i-1), (factor->coeffs) + (m-i), -i, res.appPrec );
            for (slong j=1; j<= i-1; j++) {
                compApp_mul( temp, Shs + (j-1), coeffs + (m-i)+j, res.appPrec );
                compApp_sub( Shs + (i-1), Shs + (i-1), temp, res.appPrec );
            }
        }
        realApp_zero(r);
        /* compute i-th root */
        for (slong i=1; i<=m; i++) {
            compApp_abs(compApp_realref(Shs + (i-1)), Shs + (i-1), res.appPrec);
            realApp_div_ui( compApp_realref(Shs + (i-1)), compApp_realref(Shs + (i-1)), m, res.appPrec );
            realApp_pos_root_ui(compApp_realref(Shs + (i-1)),  compApp_realref(Shs + (i-1)), (ulong) i, res.appPrec);
            realApp_zero( compApp_imagref(Shs + (i-1)) );
            realApp_max(r, r, compApp_realref(Shs + (i-1)), res.appPrec);
        }
        realApp_pos_root_ui(r,  r, n, res.appPrec);
        
        /* compute lbound */
        realRat_set_si(exp, 1, 1);
        realRat_set_si(lbound, 1, 1);
        arb_get_interval_fmpz_2exp(realRat_numref(lbound),realRat_numref(ubound),realRat_numref(exp),r);
        realRat_set_si(t, 2, 1);
        fmpq_pow_fmpz(t, t, realRat_numref(exp));
        realRat_mul(lbound, lbound, t);
        /* compute ubound */
        realRat_set_si(exp, 1, 1);
        realRat_set_si(ubound, 1, 1);
        realApp_set_si( compApp_realref(temp), 5 );
        realApp_root_ui(compApp_realref(temp), compApp_realref(temp), n, res.appPrec);
        realApp_mul(r, r, compApp_realref(temp), res.appPrec);
        arb_get_interval_fmpz_2exp(realRat_numref(t),realRat_numref(ubound),realRat_numref(exp),r);
        realRat_set_si(t, 2, 1);
        fmpq_pow_fmpz(t, t, realRat_numref(exp));
        realRat_mul(ubound, ubound, t);
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------compression.c: cauchy_rootBound_deflation_Turan; computed bounds for r1 as rationals:\n");
            realRat_print(lbound); printf("\n");
            realRat_print(ubound); printf("\n");
//             printf("#---------compression.c: cauchy_rootBound_deflation_Turan; computed bounds for r1 as reals:\n");
//             realApp_t lb, ub;
//             realApp_init(lb);
//             realApp_init(ub);
//             realApp_set_realRat(lb, lbound, res.appPrec);
//             realApp_set_realRat(ub, ubound, res.appPrec);
//             realApp_printd(lb, 10); printf("\n");
//             realApp_printd(ub, 10); printf("\n");
//             realApp_clear(lb);
//             realApp_clear(ub);
        }
        /* check if lbound*relativeError >= ubound */
        realRat_mul(t, lbound, relativeError);
        if ( (realRat_cmp(ubound, eps)<=0) ||
             (realRat_cmp(t, ubound)>=0) ){
            enoughPrec=1;
        } else {
            enoughPrec=-1;
            realRat_pow_si(epsp, epsp, 2);
            prec2 = 2*prec2;
        }
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------compression.c: cauchy_rootBound_deflation_Turan; check relative error: %d\n", (realRat_cmp(t, ubound)>=0));
        }
        
//         /* compute lbound and ubound */
//         realApp_set_si(ubound, 5);
//         realApp_root_ui(ubound, ubound, n, res.appPrec);
//         realApp_get_mid_realApp(mid, lbound);
//         realApp_get_rad_realApp(rad, lbound);
//         realApp_sub(lbound, mid, rad, res.appPrec);
//         realApp_add(mid, mid, rad, res.appPrec);
//         realApp_mul(ubound, mid, ubound, res.appPrec);
//         if (metadatas_getVerbo(meta)>=level) {
//             printf("#---------compression.c: cauchy_rootBound_deflation2; computed bounds for r1:\n");
//             realApp_printd( lbound, 10 ); printf("\n");
//             realApp_printd( ubound, 10 ); printf("\n");
//             printf("#---------compression.c: cauchy_rootBound_deflation_Turan; precision required: %ld\n", res.appPrec);
//         }
//         realApp_mul_realRat(compApp_realref(temp), lbound, relativeError, res.appPrec);
//         
//         /* check if lbound*relativeError >= ubound */
//         if ( realApp_lt( compApp_realref(temp), ubound ) == 1 ) {
//             
//         } else 
//             enoughPrec=1;
    }
        
    realRat_clear(epsp);    
    for (slong i=0; i<m; i++)
        compApp_clear( Shs + i ); 
    ccluster_free(Shs);
    compApp_poly_clear(factor);
    realRat_clear(rinv);
    compRat_clear(nc);
    compApp_clear(temp);
//     realApp_clear(mid);
//     realApp_clear(rad);
    realApp_clear(r);
    realRat_clear(t);
    realRat_clear(exp);
    
    return res.appPrec;

}

/* assume Delta = D(c,r) contains m and has isolation ratio theta >=2 */
/* computes a disk res = D(c',r') such that*/
/* Delta and res contain the same roots */
/* either r' <= eps */
/*     or res is m/((2m-2)*theta) rigid */
slong cauchy_compressionIntoRigidDisk( compDsk_t res, const compDsk_t Delta, slong m, const realRat_t theta, const realRat_t eps,
                                       cacheApp_t cache,
                                       cacheCauchy_t cacheCau,
                                       slong prec, metadatas_t meta, slong depth) {
    
    int level = 4;
    clock_t start=clock();
    
    realRat_t epsp;
    realRat_init(epsp);
    realRat_mul_si(epsp, eps, 2);
    /* if 2*eps \geq r, return the disk D(c,r/2) */
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; eps: ");
            realRat_print(eps);
            printf(", 2*eps: ");
            realRat_print(epsp);
            printf(", r:");
            realRat_print(compDsk_radiusref(Delta));
            printf("\n");
    }
        
    if ( (realRat_cmp( epsp, compDsk_radiusref(Delta) )>=0) ){
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; 2*eps >= r; \n");
        }
        
        compDsk_set(res, Delta );
        realRat_div_ui( compDsk_radiusref(res), compDsk_radiusref(res), 2);
        
        realRat_clear(epsp);
        return prec;
    }
    /* here eps < r/2 */ 
    /* set epsp = m*eps/theta */
    realRat_div(epsp, eps, theta);
    realRat_mul_si(epsp, epsp, m);
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; epsp: ");
            realRat_print(epsp);
            printf("\n");
    }
    
    clock_t start2=clock();
    /* compute a disk of radius less than m*eps/theta containing s1(p,Delta) */
    slong appPrec = cauchyTest_computeS1compDsk( res, theta, Delta, m, cache, cacheCau, epsp, meta, depth );
    metadatas_add_time_CompCen(meta, (double) (clock() - start2));
    
    /* if m=1 then res is equivalent to Delta and has radius less than eps/2 */
    /* return the disk res */
    if (m==1) {
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; m=1; use s1 \n");
        }
        
        metadatas_addComp_nb_1( meta, 1);
        metadatas_add_time_CompTot(meta, (double) (clock() - start));
    
        realRat_clear(epsp);
        return appPrec;
    }
    
    realRat_t relativeError;
    realRat_init(relativeError);
    realRat_set_si(relativeError, 19, 10);
    
    /* set epsp = eps/theta */
    realRat_div(epsp, eps, theta);
    
    /* set c' = center(s1)/m */
    compRat_div_ui( compDsk_centerref(res), compDsk_centerref(res), (ulong) m );
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; call root radii \n");
    }
    
#ifdef DEFLATION_TURAN
    
    realRat_t lbound, ubound;
    realRat_init(lbound);
    realRat_init(ubound);
    
    start2=clock();
    cauchy_rootBound_deflation_Turan( lbound, ubound, relativeError, Delta, theta, m, compDsk_centerref(res), epsp, cache, cacheCau, appPrec, meta, depth);
    metadatas_add_time_CompRR2(meta, (double) (clock() - start2));
    
    realRat_set( compDsk_radiusref(res), ubound );
    if (metadatas_getVerbo(meta)>=3) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; time spent in Turan RR: %f \n", (double) (clock() - start2)/CLOCKS_PER_SEC);
    }
    
    realRat_clear(lbound);
    realRat_clear(ubound);
    
#else
    
    /* set u st |c-c'|+r/theta <= u <= 5/4 r */
    cauchy_setRadSup( res, Delta, theta );
    /* compute r such that 1 <= r / r_{d+1-m} (s1/m, p) <= relativeError */ 
    realRat_t radInf;
    realRat_init(radInf);
    realRat_zero(radInf);
    
    /* Certified RR algo */
    start2=clock();
    cauchyRootRadii_root_radius( compDsk_centerref(res),
                                 radInf,        /* radInf = 0 < r_{d+1-m}(center, p) */
                                 compDsk_radiusref(res),        /* radSup > r_{d+1-m}(center, p) */
                                 relativeError, /* want relativeError*radInf >= radSup */ 
                                 epsp,
                                 theta,   /*isolation ratio of the disk in which is computed rr */ 
                                 m,
                                 cacheCau, cache, meta );
     metadatas_add_time_CompRR1(meta, (double) (clock() - start2));
    
     if (metadatas_getVerbo(meta)>=3) {
            printf("#---------cauchy.c: cauchy_compressionIntoRigidDisk; time spent in RR with double exponential sieve: %f \n", (double) (clock() - start2)/CLOCKS_PER_SEC);
    }
    
    realRat_clear(radInf);
    
#endif
    
    realRat_clear(epsp);
    realRat_clear(relativeError);
    
    
    metadatas_addComp_nb_p( meta, 1);
    metadatas_add_time_CompTot(meta, (double) (clock() - start));
    
    return appPrec;
}

connCmp_ptr cauchy_actualizeCCafterCompression( connCmp_ptr CC, const compDsk_t Delta, slong appPrec, metadatas_t meta ){
    
    int level = 4;
    
    realRat_t r;
    realRat_init(r);
    realRat_set(r, compDsk_radiusref(Delta) );
    realRat_mul_si(r, r, 2);
    
    connCmp_ptr nCC;
    nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init(nCC);
    
    compBox_list_ptr ltemp;
    ltemp = connCmp_boxesref(CC);
    compBox_ptr btemp;

    while (realRat_cmp( compBox_bwidthref(compBox_list_first(ltemp)), r)>0){ 
            btemp = compBox_list_pop(ltemp);                            
            subdBox_quadrisect_intersect_compDsk(ltemp, btemp, Delta);  
            compBox_clear(btemp);                                       
            ccluster_free(btemp);
        }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---------size of ltemp after bisection: %d\n", compBox_list_get_size(ltemp));
    }
    
    btemp = compBox_list_pop(ltemp);
    realRat_set(connCmp_widthref(nCC), compBox_bwidthref(btemp));
    connCmp_insert_compBox(nCC, btemp);
    while (!compBox_list_is_empty(ltemp))
        connCmp_insert_compBox(nCC, compBox_list_pop(ltemp));
    connCmp_nSols(nCC) = connCmp_nSols(CC);
    connCmp_isSep(nCC) = connCmp_isSep(CC);
    connCmp_isSepCertref(nCC) = connCmp_isSepCertref(CC);
    fmpz_set(connCmp_nwSpdref(nCC), connCmp_nwSpdref(CC));
    connCmp_appPrref(nCC) = appPrec;
    
    connCmp_clear(CC);
    ccluster_free(CC);
    
    realRat_clear(r);
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("#---------connCmp after compression: ");
        connCmp_print(nCC);
        printf("\n");
    }
    
    return nCC;
}
