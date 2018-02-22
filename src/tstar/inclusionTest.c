/* ************************************************************************** */
/*  Copyright (C) 2018 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "tstar/tstar.h"

#ifdef CCLUSTER_EXPERIMENTAL

void compApp_poly_abs( arb_poly_t res, const compApp_poly_t p, slong prec){
    compApp_srcptr pptr = p->coeffs;
    arb_ptr resptr = res->coeffs;
    const slong len = p->length;
    slong i;
    for (i = 0; i < len; i++)
        compApp_abs(resptr +i, pptr + i, prec);
    _arb_poly_set_length(res, len);
}

int tstar_inclusion_radius_th9_wn( cacheApp_t cache,
                          const realApp_t abspc,
                          const compApp_t nc,
                          const realApp_t nrad,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {
    
    int k=1;
    int res = 0;
    slong deg = cacheApp_getDegree(cache);
    slong nprec = prec;
    
    compApp_poly_t pk;
    compApp_t pkc;
    realApp_t abspkc;
    realRat_t facbin;
    fmpz_t fac;
//     realApp_t current_radius;
    realApp_t inclusion_radius;
    
    /*initialize variables*/
    compApp_poly_init2(pk, deg+1);
//     compApp_init(c);
    compApp_init(pkc);
    realApp_init(abspkc);
    realRat_init(facbin);
//     realRat_init(radtothek);
    fmpz_init(fac);
//     realApp_init(current_radius);
    realApp_init(inclusion_radius);
    
//     compApp_set_compRat(c, compDsk_centerref(d), nprec);
    tstar_getDerivative( pk, cache, nprec, (slong) 1, meta);
    tstar_evaluate_horner(pkc, pk, nc, nprec, meta, depth);
    compApp_abs(abspkc, pkc, nprec);
    realApp_get_mid_realApp(abspkc, abspkc);
    
    /* compute deg(p)!(bin(deg(p),k)*|p(c)|/|pk(c)| */
    fmpz_fac_ui(fac, (ulong) k);
    realRat_set_si(facbin, 1,1);
    fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
    fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
    realApp_div(inclusion_radius, abspc, abspkc, nprec);
    realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
    /* compute the k-th root */
//     realApp_root_ui( inclusion_radius, inclusion_radius, (ulong) k, nprec );
    /* test if inclusion_radius <= nrad */
    res = realApp_soft_compare( nrad, inclusion_radius, nprec );
    
    printf(" k: %d, inclusion radius: ", (int) k);
    realApp_printd(inclusion_radius, 20);
    printf("\n");
        
    while ((res !=1)&&(k<deg)) {
//     while ((res !=1)&&(k<1)) {
        k = k+1;
        /*compute pk*/
        tstar_getDerivative( pk, cache, nprec, (slong) k, meta);
        /* evaluate ||pk(c)| */
//         compApp_poly_evaluate_horner(pkc, pk, c, nprec);
        tstar_evaluate_horner(pkc, pk, nc, nprec, meta, depth);
        compApp_abs(abspkc, pkc, nprec);
        realApp_get_mid_realApp(abspkc, abspkc);
        
        /* compute the inclusion radius to the k*/
        fmpz_fac_ui(fac, (ulong) k);
        realRat_set_si(facbin, 1,1);
        fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
        fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
        realApp_div(inclusion_radius, abspc, abspkc, nprec);
        realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
        /* compute the k-th root */
        realApp_root_ui( inclusion_radius, inclusion_radius, (ulong) k, nprec );
        /* test if inclusion_radius <= nrad */
        res = realApp_soft_compare( nrad, inclusion_radius, nprec );
        printf(" k: %d, inclusion radius: ", (int) k);
        realApp_printd(inclusion_radius, 20);
        printf("\n");
    }
    printf("  :  , current radius:   ");
    realApp_printd(nrad, 20);
    printf("\n");
    
    if (res==1) {
        printf("---depth: %d, prec %d, success of th9, k: %d \n", depth, (int) nprec, k);
//         printf(" inclusion radius: "); 
//         realApp_printd(inclusion_radius, 20);
        
    }
    else {
//         printf("---depth: %d, prec %d, fail    of th9, k: %d \n", depth, (int) nprec, k);
//         printf(" inclusion radius: "); 
//         realApp_printd(inclusion_radius, 20);
//         printf(", current radius: ");
//         realApp_printd(current_radius, 20);
//         printf("\n");
        res = 0;
    }
    
    /*clear variables*/
    compApp_poly_clear(pk);
//     compApp_clear(c);
    compApp_clear(pkc);
    realApp_clear(abspkc);
    realRat_clear(facbin);
//     realRat_clear(radtothek);
    fmpz_clear(fac);
//     realApp_clear(current_radius);
    realApp_clear(inclusion_radius);
    
    return res;
    
}

int tstar_inclusion_radius_th10_wn( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {

     slong deg = cacheApp_getDegree(cache);
     compApp_poly_t pk;
     arb_poly_t sk;
     compApp_t c, pkc, den;
     realApp_t absc, sc, skc, abspkc, res, current_radius ;
     realRat_t epsilon, facbin, radtothek, factor;
     fmpz_t fac;
     slong k=0;
     int stop = 0;
     slong nprec = prec;
     
 //     printf("center: "); compRat_print(c); printf("\n");
     /*initialize variables*/
     compApp_poly_init2(pk, deg + 1);
     arb_poly_init2(sk, deg + 1);
     compApp_init(c);
     compApp_init(pkc);
     compApp_init(den);
     realApp_init(absc);
     realApp_init(sc);
     realApp_init(skc);
     realApp_init(abspkc);
     realApp_init(res);
     realApp_init(current_radius);
     realRat_init(epsilon);
     realRat_init(facbin);
     realRat_init(radtothek);
     realRat_init(factor);
     fmpz_init(fac);
     
     /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+2)) */
    realRat_set_si(epsilon, 2,1);
    realRat_pow_si(epsilon, epsilon, -nprec+1);
    realRat_mul_si(epsilon, epsilon, 3*( 4*(deg) +2 ) );
    /* compute 2*eps */
    realRat_set_si(factor, 2,1);
    realRat_mul(epsilon, factor, epsilon);
    
    /* compute 2epsilon s(|c|) */
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    compApp_abs(absc, c, nprec);
    tstar_getApproximation( pk, cache, nprec, meta);
    compApp_poly_abs( sk, pk, nprec); 
    arb_poly_evaluate_horner(sc, sk, absc, nprec);
    realApp_mul_realRat( sc, sc, epsilon, nprec );
    realApp_get_mid_realApp(sc, sc);
    
    stop = 0;
    while ( (stop==0) && (k < cacheApp_getDegree(cache))){
        k = k+1;
        /* compute p^(k)(c) */
        tstar_getDerivative( pk, cache, nprec, (slong) k, meta);
        tstar_evaluate_horner(pkc, pk, c, nprec, meta, depth);
        /* compute 2epsilon s^(k)(|c|) */
        arb_poly_derivative(sk,sk,nprec);
        arb_poly_evaluate_horner(skc, sk, absc, nprec);
        realApp_mul_realRat( skc, skc, epsilon, nprec );
        realApp_get_mid_realApp(skc, skc);
        /* test if |P^(k)(c)|>2epsilon s^(k)(|c|) */
        compApp_abs(abspkc, pkc, nprec);
        realApp_get_mid_realApp(abspkc, abspkc);
        stop = realApp_soft_compare( abspkc, skc, nprec );
        
        if (stop==1) {
            /*compute the radius*/ 
            /* compute |P^(k)(c) - 2epsilon s^(k)(|c|)| in abspkc */
            compApp_zero(den);
            compApp_set_real_realApp(den, skc);
            compApp_sub(den, pkc, den, nprec);
            compApp_abs(abspkc, den, nprec);
            realApp_get_mid_realApp(abspkc, abspkc);
            /* compute 2epsilon s(|c|) / |P^(k)(c) - 2epsilon s^(k)(|c|)| in res */
            realApp_div(res, sc, abspkc, nprec);
            /* compute the inclusion radius ^k in res */
            fmpz_fac_ui(fac, (ulong) k);
            realRat_set_si(facbin, 1,1);
            fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
            fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
            realApp_mul_realRat( res, res, facbin, nprec );
            /* compute the current radius to the k in current_radius */
            realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
            realApp_set_realRat(current_radius, radtothek, nprec);
            stop = realApp_soft_compare( current_radius, res, nprec );
            printf(" k: %d, inclusion radius: ", (int) k);
            realApp_printd(res, 20);
            printf("\n");
            printf("  :  , current radius:   ");
            realApp_printd(current_radius, 20);
            printf("\n");
        }
        
        else stop = 0;
        
    }
     
     
     
     if (stop==1) {
        printf("---depth: %d, prec %d, success of th10, k: %d \n", depth, (int) nprec, (int) k);
//         printf(" inclusion radius: "); 
//         realApp_printd(res, 20);
     }
     else {
        printf("---depth: %d, prec %d, fail    of th10, k: %d \n", depth, (int) nprec, (int) k);
//         printf(" inclusion radius: "); 
//         realApp_printd(res, 20);
         
         stop = 0;
     }
     
     /*clear variables*/
     compApp_poly_clear(pk);
     arb_poly_clear(sk);
     compApp_clear(c);
     compApp_clear(pkc);
     compApp_clear(den);
     realApp_clear(absc);
     realApp_clear(sc);
     realApp_clear(skc);
     realApp_clear(abspkc);
     realApp_clear(res);
     realApp_clear(current_radius);
     realRat_clear(epsilon);
     realRat_clear(facbin);
     realRat_clear(radtothek);
     realRat_clear(factor);
     fmpz_clear(fac);
     
     return stop;

}
    
int tstar_inclusion_test_wn( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {
    int choice = 0;
    int res = 0;
    slong deg = cacheApp_getDegree(cache);
    slong nprec = prec;
    compApp_poly_t p;
    arb_poly_t s;
    compApp_t c, nc, pc, ppc;
    realApp_t absc, sabsc, abspc, nrad;
    realRat_t factor, epsilon;
    
    /*initialize variables*/
    compApp_poly_init2(p, deg+1);
    arb_poly_init2(s, deg+1);
    compApp_init(c);
    compApp_init(nc);
    compApp_init(pc);
    compApp_init(ppc);
    realApp_init(absc);
    realApp_init(sabsc);
    realApp_init(abspc);
    realApp_init(nrad);
    realRat_init(factor);
    realRat_init(epsilon);
    
    /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+2)) */
    realRat_set_si(epsilon, 2,1);
    realRat_pow_si(epsilon, epsilon, -nprec+1);
    realRat_mul_si(epsilon, epsilon, 3*( 4*(deg) +2 ) );
    /* compute (2/3)*eps */
    realRat_set_si(factor, 2,3);
    realRat_mul(epsilon, factor, epsilon);
    
    /* apply a newton step */
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    tstar_getApproximation( p, cache, nprec, meta);
    tstar_evaluate(pc, p, c, nprec, meta, depth);
//     printf("value before newton iteration: "); compApp_printd(pc, 20); printf("\n");
    tstar_getDerivative( p, cache, nprec, 1, meta);
    tstar_evaluate(ppc, p, c, nprec, meta, depth);
    compApp_div(ppc, pc, ppc, nprec);
    compApp_sub(nc, c, ppc, nprec);
    
//     tstar_getApproximation( p, cache, nprec, meta);
//     tstar_evaluate(pc, p, nc, nprec, meta, depth);
//     tstar_getDerivative( p, cache, nprec, 1, meta);
//     tstar_evaluate(ppc, p, nc, nprec, meta, depth);
//     compApp_div(ppc, pc, ppc, nprec);
//     compApp_sub(nc, nc, ppc, nprec);
    
    /* check if the newton iterate is in d */
    compApp_sub(c, c, nc, nprec);
    compApp_abs(absc, c, nprec);
    realApp_set_realRat(sabsc, compDsk_radiusref(d), nprec);
    res = realApp_soft_compare( sabsc, absc, nprec );
    
    if (res == 1) { /* the newton iterate is in d */
        /* set c to the newton iterate and nrad to the new radius*/
        compApp_set(c, nc);
        realApp_sub(nrad, sabsc, absc, nprec);
        
        tstar_getApproximation( p, cache, nprec, meta);
        /*construct S from P */
        compApp_poly_abs( s, p, nprec);
        /*construct approximation of |c|*/
        compApp_abs(absc, c, nprec);
        /*evaluate polynomials p and s in c and |c|*/
        tstar_evaluate_horner(pc, p, c, nprec, meta, depth);
//         printf("value after  newton iteration: "); compApp_printd(pc, 20); printf("\n");
        arb_poly_evaluate_horner(sabsc, s, absc, nprec);
        compApp_abs(abspc, pc, nprec);
        realApp_mul_realRat(sabsc, sabsc, epsilon, nprec);
        realApp_get_mid_realApp(abspc, abspc);
        choice = realApp_soft_compare( abspc, sabsc, nprec );
        
        if (choice==1) {
//         printf("---depth: %d, prec before: %d, prec after: %d, inclusion radius: apply thm 9\n", depth, (int) prec, (int) nprec);
//         printf("   res of thm 9: %d\n", tstar_inclusion_radius_th9(cache, d, prec, depth, meta) );
            res = tstar_inclusion_radius_th9_wn(cache, abspc, nc, nrad, d, nprec, depth, meta);
//             printf("---depth: %d, prec: %d, inclusion radius wn: apply thm  9, success: %d\n", depth, (int) nprec, res); 
        }
        if (choice==0) {
//             res = tstar_inclusion_radius_th10_wn(cache, d, nprec, depth, meta);
            printf("---depth: %d, prec: %d, inclusion radius wn: apply thm 10, success: %d\n", depth, (int) nprec, 0);  
        }
    }
    else {
        res = 0;
//         printf("---depth: %d, prec: %d, inclusion radius wn: fail\n", depth, (int) nprec); 
    }
    
    
    
    /*clear variables*/
    compApp_poly_clear(p);
    arb_poly_clear(s);
    compApp_clear(c);
    compApp_clear(nc);
    compApp_clear(pc);
    compApp_clear(ppc);
    realApp_clear(absc);
    realApp_clear(sabsc);
    realApp_clear(abspc);
    realApp_clear(nrad);
    realRat_clear(factor);
    realRat_clear(epsilon);
    
    return res;
}

int tstar_inclusion_radius_th9( cacheApp_t cache,
                          const realApp_t abspc,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {
    
    int k=1;
    int res = 0;
    slong deg = cacheApp_getDegree(cache);
    slong nprec = prec;
    
    compApp_poly_t pk;
    compApp_t c, pkc;
    realApp_t abspkc;
    realRat_t facbin, radtothek;
    fmpz_t fac;
    realApp_t current_radius;
    realApp_t inclusion_radius;
    
    /*initialize variables*/
    compApp_poly_init2(pk, deg+1);
    compApp_init(c);
    compApp_init(pkc);
    realApp_init(abspkc);
    realRat_init(facbin);
    realRat_init(radtothek);
    fmpz_init(fac);
    realApp_init(current_radius);
    realApp_init(inclusion_radius);
    
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    tstar_getDerivative( pk, cache, nprec, (slong) 1, meta);
//     compApp_poly_evaluate_horner(pkc, pk, c, nprec);
    tstar_evaluate_horner(pkc, pk, c, nprec, meta, depth);
    compApp_abs(abspkc, pkc, nprec);
    realApp_get_mid_realApp(abspkc, abspkc);
    
    /* compute deg(p)*|p(c)|/|pk(c)| */
    fmpz_fac_ui(fac, (ulong) k);
    realRat_set_si(facbin, 1,1);
    fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
    fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
    realApp_div(inclusion_radius, abspc, abspkc, nprec);
    realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
    /* test if inclusion_radius <= radius of d^k */
    realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
    realApp_set_realRat(current_radius, radtothek, nprec);
    res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
    
    while ((res !=1)&&(k<deg)) {
        k = k+1;
        /*compute pk*/
        tstar_getDerivative( pk, cache, nprec, (slong) k, meta);
        /* evaluate ||pk(c)| */
//         compApp_poly_evaluate_horner(pkc, pk, c, nprec);
        tstar_evaluate_horner(pkc, pk, c, nprec, meta, depth);
        compApp_abs(abspkc, pkc, nprec);
        realApp_get_mid_realApp(abspkc, abspkc);
        
        /* compute the inclusion radius to the k*/
        fmpz_fac_ui(fac, (ulong) k);
        realRat_set_si(facbin, 1,1);
        fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
        fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
        realApp_div(inclusion_radius, abspc, abspkc, nprec);
        realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
        realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
        realApp_set_realRat(current_radius, radtothek, nprec);
        res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
//         printf(" k: %d, inclusion radius: ", (int) k);
//         realApp_printd(inclusion_radius, 20);
//         printf("\n");
//         printf("  :  , current radius:   ");
//         realApp_printd(current_radius, 20);
//         printf("\n");
    }
    
    if (res==1) {
//         printf("---depth: %d, prec %d, success of th9, k: %d \n", depth, (int) nprec, k);
//         printf(" inclusion radius: "); 
//         realApp_printd(inclusion_radius, 20);
        
    }
    else {
//         printf("---depth: %d, prec %d, fail    of th9, k: %d \n", depth, (int) nprec, k);
//         printf(" inclusion radius: "); 
//         realApp_printd(inclusion_radius, 20);
//         printf(", current radius: ");
//         realApp_printd(current_radius, 20);
//         printf("\n");
        res = 0;
    }
    
    /*clear variables*/
    compApp_poly_clear(pk);
    compApp_clear(c);
    compApp_clear(pkc);
    realApp_clear(abspkc);
    realRat_clear(facbin);
    realRat_clear(radtothek);
    fmpz_clear(fac);
    realApp_clear(current_radius);
    realApp_clear(inclusion_radius);
    
    return res;
    
}

int tstar_inclusion_radius_th10( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {

     slong deg = cacheApp_getDegree(cache);
     compApp_poly_t pk;
     arb_poly_t sk;
     compApp_t c, pkc, den;
     realApp_t absc, sc, skc, abspkc, res, current_radius ;
     realRat_t epsilon, facbin, radtothek, factor;
     fmpz_t fac;
     slong k=0;
     int stop = 0;
     slong nprec = prec;
     
 //     printf("center: "); compRat_print(c); printf("\n");
     /*initialize variables*/
     compApp_poly_init2(pk, deg + 1);
     arb_poly_init2(sk, deg + 1);
     compApp_init(c);
     compApp_init(pkc);
     compApp_init(den);
     realApp_init(absc);
     realApp_init(sc);
     realApp_init(skc);
     realApp_init(abspkc);
     realApp_init(res);
     realApp_init(current_radius);
     realRat_init(epsilon);
     realRat_init(facbin);
     realRat_init(radtothek);
     realRat_init(factor);
     fmpz_init(fac);
     
     /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+2)) */
    realRat_set_si(epsilon, 2,1);
    realRat_pow_si(epsilon, epsilon, -nprec+1);
    realRat_mul_si(epsilon, epsilon, 3*( 4*(deg) +2 ) );
    /* compute 2*eps */
    realRat_set_si(factor, 2,1);
    realRat_mul(epsilon, factor, epsilon);
    
    /* compute 2epsilon s(|c|) */
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    compApp_abs(absc, c, nprec);
    tstar_getApproximation( pk, cache, nprec, meta);
    compApp_poly_abs( sk, pk, nprec); 
    arb_poly_evaluate_horner(sc, sk, absc, nprec);
    realApp_mul_realRat( sc, sc, epsilon, nprec );
    realApp_get_mid_realApp(sc, sc);
    
    stop = 0;
    while ( (stop==0) && (k < cacheApp_getDegree(cache))){
        k = k+1;
        /* compute p^(k)(c) */
        tstar_getDerivative( pk, cache, nprec, (slong) k, meta);
        tstar_evaluate_horner(pkc, pk, c, nprec, meta, depth);
        /* compute 2epsilon s^(k)(|c|) */
        arb_poly_derivative(sk,sk,nprec);
        arb_poly_evaluate_horner(skc, sk, absc, nprec);
        realApp_mul_realRat( skc, skc, epsilon, nprec );
        realApp_get_mid_realApp(skc, skc);
        /* test if |P^(k)(c)|>2epsilon s^(k)(|c|) */
        compApp_abs(abspkc, pkc, nprec);
        realApp_get_mid_realApp(abspkc, abspkc);
        stop = realApp_soft_compare( abspkc, skc, nprec );
        
        if (stop==1) {
            /*compute the radius*/ 
            /* compute |P^(k)(c) - 2epsilon s^(k)(|c|)| in abspkc */
            compApp_zero(den);
            compApp_set_real_realApp(den, skc);
            compApp_sub(den, pkc, den, nprec);
            compApp_abs(abspkc, den, nprec);
            realApp_get_mid_realApp(abspkc, abspkc);
            /* compute 2epsilon s(|c|) / |P^(k)(c) - 2epsilon s^(k)(|c|)| in res */
            realApp_div(res, sc, abspkc, nprec);
            /* compute the inclusion radius ^k in res */
            fmpz_fac_ui(fac, (ulong) k);
            realRat_set_si(facbin, 1,1);
            fmpz_bin_uiui(realRat_numref(facbin), (ulong) deg, (ulong) k);
            fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
            realApp_mul_realRat( res, res, facbin, nprec );
            /* compute the current radius to the k in current_radius */
            realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
            realApp_set_realRat(current_radius, radtothek, nprec);
            stop = realApp_soft_compare( current_radius, res, nprec );
            printf(" k: %d, inclusion radius: ", (int) k);
            realApp_printd(res, 20);
            printf("\n");
            printf("  :  , current radius:   ");
            realApp_printd(current_radius, 20);
            printf("\n");
        }
        
        else stop = 0;
        
    }
     
     
     
     if (stop==1) {
        printf("---depth: %d, prec %d, success of th10, k: %d \n", depth, (int) nprec, (int) k);
//         printf(" inclusion radius: "); 
//         realApp_printd(res, 20);
     }
     else {
        printf("---depth: %d, prec %d, fail    of th10, k: %d \n", depth, (int) nprec, (int) k);
//         printf(" inclusion radius: "); 
//         realApp_printd(res, 20);
         
         stop = 0;
     }
     
     /*clear variables*/
     compApp_poly_clear(pk);
     arb_poly_clear(sk);
     compApp_clear(c);
     compApp_clear(pkc);
     compApp_clear(den);
     realApp_clear(absc);
     realApp_clear(sc);
     realApp_clear(skc);
     realApp_clear(abspkc);
     realApp_clear(res);
     realApp_clear(current_radius);
     realRat_clear(epsilon);
     realRat_clear(facbin);
     realRat_clear(radtothek);
     realRat_clear(factor);
     fmpz_clear(fac);
     
     return stop;

}
    
int tstar_inclusion_test( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {
    int choice = 0;
    int res = 0;
    slong deg = cacheApp_getDegree(cache);
    slong nprec = prec;
    compApp_poly_t p;
    arb_poly_t s;
    compApp_t c, pc;
    realApp_t absc, sabsc, abspc;
    realRat_t factor, epsilon;
    
    /*initialize variables*/
    compApp_poly_init2(p, deg+1);
    arb_poly_init2(s, deg+1);
    compApp_init(c);
    compApp_init(pc);
    realApp_init(absc);
    realApp_init(sabsc);
    realApp_init(abspc);
    realRat_init(factor);
    realRat_init(epsilon);
    
    /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+2)) */
    realRat_set_si(epsilon, 2,1);
    realRat_pow_si(epsilon, epsilon, -nprec+1);
    realRat_mul_si(epsilon, epsilon, 3*( 4*(deg) +2 ) );
    /* compute (2/3)*eps */
    realRat_set_si(factor, 2,3);
    realRat_mul(epsilon, factor, epsilon);
    
    tstar_getApproximation( p, cache, nprec, meta);
    /*construct S from P */
    compApp_poly_abs( s, p, nprec);
    /*construct approximations of c and |c|*/
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    compApp_abs(absc, c, nprec);
    /*evaluate polynomials p and s in c and |c|*/
//     compApp_poly_evaluate_horner(pc, p, c, nprec);
    tstar_evaluate_horner(pc, p, c, nprec, meta, depth);
    arb_poly_evaluate_horner(sabsc, s, absc, nprec);
    compApp_abs(abspc, pc, nprec);
    realApp_mul_realRat(sabsc, sabsc, epsilon, nprec);
//     printf("|p(c)|:         "); realApp_printd(abspc,50); printf("\n");
    realApp_get_mid_realApp(abspc, abspc);
//     printf("|p(c)|:         "); realApp_printd(abspc,50); printf("\n");
//     printf("2/3*eps*s(|c|): "); realApp_printd(sabsc,50); printf("\n");
    choice = realApp_soft_compare( abspc, sabsc, nprec );
    
//     while( choice == -2 ){
//         nprec *=2;
//         
//         /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+1)) */
//         realRat_set_si(epsilon, 2,1);
//         realRat_pow_si(epsilon, epsilon, -nprec+1);
//         realRat_mul_si(epsilon, epsilon, 3*( 4*(cacheApp_getDegree(cache)) +1 ) );
//         /* compute (2/3)*eps */
//         realRat_set_si(factor, 2,3);
//         realRat_mul(epsilon, factor, epsilon);
//     
//         tstar_getApproximation( p, cache, nprec, meta);
//         /*construct S from P */
//         compApp_poly_abs( s, p, nprec);
//         /*construct approximations of c and |c|*/
//         compApp_set_compRat(c, compDsk_centerref(d), nprec);
//         compApp_abs(absc, c, nprec);
//         /*evaluate polynomials p and s in c and |c|*/
//         compApp_poly_evaluate(pc, p, c, nprec);
//         arb_poly_evaluate_horner(sabsc, s, absc, nprec);
//         compApp_abs(abspc, pc, nprec);
//         realApp_get_mid_realApp(abspc, abspc);
//         realApp_mul_realRat(sabsc, sabsc, epsilon, nprec);
//         choice = realApp_soft_compare( abspc, sabsc, nprec );
//         
//     }
    
    if (choice==1) {
//         printf("---depth: %d, prec before: %d, prec after: %d, inclusion radius: apply thm 9\n", depth, (int) prec, (int) nprec);
//         printf("   res of thm 9: %d\n", tstar_inclusion_radius_th9(cache, d, prec, depth, meta) );
        res = tstar_inclusion_radius_th9(cache, abspc, d, nprec, depth, meta);
//         printf("---depth: %d, prec: %d, inclusion radius: apply thm  9, success: %d\n", depth, (int) nprec, res); 
        
    }
    if (choice==0) {
        res = tstar_inclusion_radius_th10(cache, d, nprec, depth, meta);
        printf("---depth: %d, prec: %d, inclusion radius: apply thm 10, success: %d\n", depth, (int) nprec, res);  
//         
//         printf("--------------------------------------------\n");
    }
    /*clear variables*/
    compApp_poly_clear(p);
    arb_poly_clear(s);
    compApp_clear(c);
    compApp_clear(pc);
    realApp_clear(absc);
    realApp_clear(sabsc);
    realApp_clear(abspc);
    realRat_clear(factor);
    realRat_clear(epsilon);
    
    return res;
}


// int tstar_inclusion_radius( realApp_t res, /*the inclusion radius*/
//                             cacheApp_t cache,
//                             const compRat_t c,
//                             slong prec,
//                             metadatas_t meta) {
//     
// }

// int tstar_inclusion_test( cacheApp_t cache,
//                           const compDsk_t d,
//                           slong prec,
//                           int depth, /* just for display*/
//                           metadatas_t meta) {
//     int res=0;
//     realApp_t radComp, rad;
//     realApp_init(radComp);
//     realApp_init(rad);
//     
//     int k = tstar_inclusion_radius( radComp, cache, compDsk_centerref(d), prec, meta);
//     
//     if (k>0){
//         realApp_set_realRat( rad, compDsk_radiusref(d), prec);
//         
//         res = realApp_le(radComp, rad);
//         if (res) {
//             printf("--depth: %d, succes of inclusion test, k: %d\n", depth, k);
//             printf("  r:      "); realApp_printd(radComp, 50); printf("\n"); 
//             printf("  radius: "); realApp_printd(rad, 50);     printf("\n"); 
//             printf("  radius: "); realRat_print(compDsk_radiusref(d));     printf("\n");
//         }
//     }
//     realApp_clear(radComp);
//     realApp_clear(rad);
//     return res;
// }
#endif