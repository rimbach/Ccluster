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

void compApp_poly_abs( arb_poly_t res, const compApp_poly_t p, slong prec){
    compApp_srcptr pptr = p->coeffs;
    arb_ptr resptr = res->coeffs;
    const slong len = p->length;
    slong i;
    for (i = 0; i < len; i++)
        compApp_abs(resptr +i, pptr + i, prec);
    _arb_poly_set_length(res, len);
}

int tstar_inclusion_radius_th9( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta) {
    
    int k=1;
    int res = 0;
    slong nprec = prec;
    
    compApp_poly_t p, pk;
    compApp_t c, pc, pkc;
    realApp_t abspc, abspkc;
    realRat_t facbin, radtothek;
    fmpz_t fac;
    realApp_t current_radius;
    realApp_t inclusion_radius;
    
    /*initialize variables*/
    compApp_poly_init2(p, cacheApp_getDegree(cache)+1);
    compApp_poly_init2(pk, cacheApp_getDegree(cache)+1);
    compApp_init(c);
    compApp_init(pc);
    compApp_init(pkc);
    realApp_init(abspc);
    realApp_init(abspkc);
    realRat_init(facbin);
    realRat_init(radtothek);
    fmpz_init(fac);
    realApp_init(current_radius);
    realApp_init(inclusion_radius);
//     realApp_init(temp);
    
    /* construct p and its first derivative pk */
    tstar_getApproximation( p, cache, nprec, meta);
    tstar_getApproximation( pk, cache, nprec, meta);
    compApp_poly_derivative(pk,pk,nprec);
    /* evaluate |p(c)| and |pk(c)| */
    compApp_set_compRat(c, compDsk_centerref(d), nprec);
    compApp_poly_evaluate(pc, p, c, nprec);
    compApp_poly_evaluate(pkc, pk, c, nprec);
    compApp_abs(abspc, pc, nprec);
    compApp_abs(abspkc, pkc, nprec);
    realApp_get_mid_realApp(abspc, abspc);
    realApp_get_mid_realApp(abspkc, abspkc);
    /* compute deg(p)*|p(c)|/|pk(c)| */
    fmpz_fac_ui(fac, (ulong) k);
    realRat_set_si(facbin, 1,1);
    fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
//     printf("------k: %d, factor: "); realRat_print(facbin); printf("\n");
    fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
    realApp_div(inclusion_radius, abspc, abspkc, nprec);
    realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
    /* test if inclusion_radius <= radius of d^k */
    realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
    realApp_set_realRat(current_radius, radtothek, nprec);
    res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
    
    while (res == -2) {
        nprec *=2;
        tstar_getApproximation( p, cache, nprec, meta);
        tstar_getApproximation( pk, cache, nprec, meta);
        compApp_poly_derivative(pk,pk,nprec);
        /* evaluate |p(c)| and |pk(c)| */
        compApp_set_compRat(c, compDsk_centerref(d), nprec);
        compApp_poly_evaluate(pc, p, c, nprec);
        compApp_poly_evaluate(pkc, pk, c, nprec);
        compApp_abs(abspc, pc, nprec);
        compApp_abs(abspkc, pkc, nprec);
        /* compute deg(p)*|p(c)|/|pk(c)| */
        fmpz_fac_ui(fac, (ulong) k);
        realRat_set_si(facbin, 1,1);
        fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
    //     printf("------k: %d, factor: "); realRat_print(facbin); printf("\n");
        fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
        realApp_div(inclusion_radius, abspc, abspkc, nprec);
        realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
        /* test if inclusion_radius <= radius of d^k */
        realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
        realApp_set_realRat(current_radius, radtothek, nprec);
        res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
    }
    
    while ((res !=1)&&(k<cacheApp_getDegree(cache))) {
        k = k+1;
        /*compute pk*/
        compApp_poly_derivative(pk,pk,nprec);
        /* evaluate ||pk(c)| */
        compApp_poly_evaluate(pkc, pk, c, nprec);
        compApp_abs(abspkc, pkc, nprec);
        /* compute the inclusion radius to the k*/
        fmpz_fac_ui(fac, (ulong) k);
        realRat_set_si(facbin, 1,1);
        fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
        fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
        realApp_div(inclusion_radius, abspc, abspkc, nprec);
        realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
        realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
        realApp_set_realRat(current_radius, radtothek, nprec);
        res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
        
        while (res == -2) {
            nprec *=2;
            tstar_getApproximation( p, cache, nprec, meta);
            tstar_getApproximation( pk, cache, nprec, meta);
            slong i;
            for(i=0;i<k;i++)
                compApp_poly_derivative(pk,pk,nprec);
            /* evaluate |p(c)| and |pk(c)| */
            compApp_set_compRat(c, compDsk_centerref(d), nprec);
            compApp_poly_evaluate(pc, p, c, nprec);
            compApp_poly_evaluate(pkc, pk, c, nprec);
            compApp_abs(abspc, pc, nprec);
            compApp_abs(abspkc, pkc, nprec);
            /* compute deg(p)*|p(c)|/|pk(c)| */
            fmpz_fac_ui(fac, (ulong) k);
            realRat_set_si(facbin, 1,1);
            fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
        //     printf("------k: %d, factor: "); realRat_print(facbin); printf("\n");
            fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
            realApp_div(inclusion_radius, abspc, abspkc, nprec);
            realApp_mul_realRat(inclusion_radius, inclusion_radius, facbin, nprec);
            /* test if inclusion_radius <= radius of d^k */
            realRat_pow_si(radtothek, compDsk_radiusref(d), (slong) k);
            realApp_set_realRat(current_radius, radtothek, nprec);
            res = realApp_soft_compare( current_radius, inclusion_radius, nprec );
        }
    }
    
    if (res==1) {
        printf("---depth: %d, prec before: %d, prec after: %d, success of th9. \n", depth, prec, nprec);
        printf("   k: %d, inclusion radius: ", k); 
        realApp_printd(inclusion_radius, 20);
        printf(", current radius: ");
        realApp_printd(current_radius, 20);
        printf("\n");
    }
    
    /*clear variables*/
    compApp_poly_clear(p);
    compApp_poly_clear(pk);
    compApp_clear(c);
    compApp_clear(pc);
    compApp_clear(pkc);
    realApp_clear(abspc);
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

     compApp_poly_t p, pk;
     arb_poly_t s, sk;
     compApp_t c, pc, pkc, den;
     realApp_t absc, sc, skc, abspkc, res ;
     realRat_t epsilon, facbin, tothek;
     fmpz_t fac;
     slong k=0;
     int stop = 0;
     slong nprec = prec;
     
 //     printf("center: "); compRat_print(c); printf("\n");
     /*initialize variables*/
     compApp_poly_init2(p, cacheApp_getDegree(cache)+1);
     compApp_poly_init2(pk, cacheApp_getDegree(cache)+1);
     arb_poly_init2(s, cacheApp_getDegree(cache)+1);
     arb_poly_init2(sk, cacheApp_getDegree(cache)+1);
     compApp_init(c);
     compApp_init(pc);
     compApp_init(pkc);
     compApp_init(den);
     realApp_init(absc);
     realApp_init(sc);
     realApp_init(skc);
     realApp_init(abspkc);
     realApp_init(res);
     realRat_init(epsilon);
     realRat_init(facbin);
     realRat_init(tothek);
     fmpz_init(fac);
     
     /*compute 2*epsilon*/
     realRat_set_si(epsilon, 2,1);
     realRat_pow_si(epsilon, epsilon, -nprec+1);
     realRat_mul_si(epsilon, epsilon, 3*( 4*(cacheApp_getDegree(cache)) +1 ) ); /*test*/
     realRat_mul_si(epsilon, epsilon, 2);
     
     tstar_getApproximation( pk, cache, nprec, meta);
     compApp_poly_abs( sk, pk, nprec);
     compApp_set_compRat(c, compDsk_centerref(d), nprec);
     compApp_abs(absc, c, nprec);
     stop = 0;
     /*find appropriated k*/
     while ( (stop==0) && (k < cacheApp_getDegree(cache))){
         k = k+1;
         /*compute the derivatives of p2prec and sprec*/
         compApp_poly_derivative(pk,pk,nprec);
         arb_poly_derivative(sk,sk,nprec);
         /*evaluate the derivatives*/
         compApp_poly_evaluate(pkc, pk, c, nprec);
         arb_poly_evaluate_horner(skc, sk, absc, nprec);
         compApp_abs(abspkc, pkc, nprec);
         realApp_get_mid_realApp(abspkc, abspkc);
         realApp_mul_realRat( skc, skc, epsilon, nprec );
         stop = realApp_soft_compare( abspkc, skc, nprec );
         
         while (stop==-2) {
             nprec *=2;
             tstar_getApproximation( pk, cache, nprec, meta);
             compApp_poly_abs( sk, pk, nprec);
             slong i;
             for(i=0;i<k;i++) {
                compApp_poly_derivative(pk,pk,nprec);
                arb_poly_derivative(sk,sk,nprec);
             }
             compApp_set_compRat(c, compDsk_centerref(d), nprec);
             compApp_abs(absc, c, nprec);
             compApp_poly_evaluate(pkc, pk, c, nprec);
             arb_poly_evaluate_horner(skc, sk, absc, nprec);
             compApp_abs(abspkc, pkc, nprec);
             realApp_get_mid_realApp(abspkc, abspkc);
             realRat_set_si(epsilon, 2,1);
             realRat_pow_si(epsilon, epsilon, -nprec+1);
             realRat_mul_si(epsilon, epsilon, 3*( 4*(cacheApp_getDegree(cache)) +1 ) ); /*test*/
             realRat_mul_si(epsilon, epsilon, 2);
             realApp_mul_realRat( skc, skc, epsilon, nprec );
             stop = realApp_soft_compare( abspkc, skc, nprec );
         }
         if (stop==-1) stop=0;
     }
     
     /*compute the radius*/
     tstar_getApproximation( p, cache, nprec, meta);
     /*construct S from P */
     compApp_poly_abs( s, p, nprec);
//      compApp_poly_evaluate(pc, p, c, nprec);
     arb_poly_evaluate_horner(sc, s, absc, nprec);
     /*compute 2*epsilon*s(|c|) */
     realApp_mul_realRat( sc, sc, epsilon, nprec );
     compApp_zero(den);
     compApp_set_real_realApp(den, skc);
     realApp_get_mid_realApp(compApp_realref(pkc), compApp_realref(pkc));
     realApp_get_mid_realApp(compApp_imagref(pkc), compApp_imagref(pkc));
     compApp_sub(den, pkc, den, nprec);
     compApp_abs(abspkc, den, nprec);
     realApp_div(res, sc, abspkc, nprec);
     fmpz_fac_ui(fac, (ulong) k);
     realRat_set_si(facbin, 1,1);
     fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
     fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
     realApp_mul_realRat( res, res, facbin, nprec );
     realApp_root_ui(res, res, (ulong) k, nprec );
     
//      realRat_pow_si(tothek, compDsk_radiusref(d), (slong) k);
//      realApp_set_realRat(sc, tothek, nprec);
     realApp_set_realRat(sc, compDsk_radiusref(d), nprec);
     
     stop = realApp_soft_compare( sc, res, nprec );
     
     if (stop==1) {
         printf("---depth: %d, prec before: %d, prec after: %d, success of th10\n", depth, prec, nprec);
         printf("k:              %d\n", k);
//          compApp_abs(abspkc, pc, nprec);
//          printf("|p(c)|:             "); realApp_printd(abspkc, 50); printf("\n");
//          arb_poly_evaluate_horner(abspkc, s, absc, nprec);
//          printf("s(|c|):             "); realApp_printd(abspkc, 50); printf("\n");
         printf("res to k:           "); realApp_printd(res, 50); printf("\n");
         printf("actual radius to k: "); realApp_printd( sc, 50 ); printf("\n");
     }
     else {
         printf("---depth: %d, prec before: %d, prec after: %d, failure of th10, stop: %d\n\n", depth, prec, nprec, stop);
         printf("k:             %d\n", k);
         printf("res          : "); realApp_printd(res, 50); printf("\n");
         printf("actual radius: "); realApp_printd( sc, 50 ); printf("\n");
         
//          while (k < cacheApp_getDegree(cache)) {
//              k=k+1;
//              compApp_poly_derivative(pk,pk,nprec);
//              arb_poly_derivative(sk,sk,nprec);
//              /*evaluate the derivatives*/
//              compApp_poly_evaluate(pkc, pk, c, nprec);
//              arb_poly_evaluate_horner(skc, sk, absc, nprec);
//              compApp_abs(abspkc, pkc, nprec);
//              realApp_get_mid_realApp(abspkc, abspkc);
//              realApp_mul_realRat( skc, skc, epsilon, nprec );
//              stop = realApp_soft_compare( abspkc, skc, nprec );
//              if (stop){
//                 arb_poly_evaluate_horner(sc, s, absc, nprec); 
//                 realApp_mul_realRat( sc, sc, epsilon, nprec );
//                 compApp_zero(den);
//                 compApp_set_real_realApp(den, skc);
//                 realApp_get_mid_realApp(compApp_realref(pkc), compApp_realref(pkc));
//                 realApp_get_mid_realApp(compApp_imagref(pkc), compApp_imagref(pkc));
//                 compApp_sub(den, pkc, den, nprec);
//                 compApp_abs(abspkc, den, nprec);
//                 realApp_div(res, sc, abspkc, nprec);
//                 fmpz_fac_ui(fac, (ulong) k);
//                 realRat_set_si(facbin, 1,1);
//                 fmpz_bin_uiui(realRat_numref(facbin), (ulong) cacheApp_getDegree(cache), (ulong) k);
//                 fmpz_mul(realRat_numref(facbin), realRat_numref(facbin), fac);
//                 realApp_mul_realRat( res, res, facbin, nprec );
//                 realApp_root_ui(res, res, (ulong) k, nprec );
//                 printf("k:             %d\n", k);
//                 printf("res          : "); realApp_printd(res, 50); printf("\n");
//              }
//          }
         
         stop = 0;
     }
     
     /*clear variables*/
     compApp_poly_clear(p);
     compApp_poly_clear(pk);
     arb_poly_clear(s);
     arb_poly_clear(sk);
     compApp_clear(c);
     compApp_clear(pc);
     compApp_clear(pkc);
     compApp_clear(den);
     realApp_clear(absc);
     realApp_clear(sc);
     realApp_clear(skc);
     realApp_clear(abspkc);
     realApp_clear(res);
     realRat_clear(epsilon);
     realRat_clear(facbin);
     realRat_clear(tothek);
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
    slong nprec = prec;
    compApp_poly_t p;
    arb_poly_t s;
    compApp_t c, pc;
    realApp_t absc, sabsc, abspc;
    realRat_t factor, epsilon;
    
    /*initialize variables*/
    compApp_poly_init2(p, cacheApp_getDegree(cache)+1);
    arb_poly_init2(s, cacheApp_getDegree(cache)+1);
    compApp_init(c);
    compApp_init(pc);
    realApp_init(absc);
    realApp_init(sabsc);
    realApp_init(abspc);
    realRat_init(factor);
    realRat_init(epsilon);
    
    /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+1)) */
    realRat_set_si(epsilon, 2,1);
    realRat_pow_si(epsilon, epsilon, -nprec+1);
    realRat_mul_si(epsilon, epsilon, 3*( 4*(cacheApp_getDegree(cache)) +1 ) );
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
    compApp_poly_evaluate(pc, p, c, nprec);
    arb_poly_evaluate_horner(sabsc, s, absc, nprec);
    compApp_abs(abspc, pc, nprec);
    realApp_mul_realRat(sabsc, sabsc, epsilon, nprec);
//     printf("|p(c)|:         "); realApp_printd(abspc,50); printf("\n");
    realApp_get_mid_realApp(abspc, abspc);
//     printf("|p(c)|:         "); realApp_printd(abspc,50); printf("\n");
//     printf("2/3*eps*s(|c|): "); realApp_printd(sabsc,50); printf("\n");
    choice = realApp_soft_compare( abspc, sabsc, nprec );
    
    while( choice == -2 ){
        nprec *=2;
        
        /* compute eps = (2^(-prec + 1))*3*((4*deg(p)+1)) */
        realRat_set_si(epsilon, 2,1);
        realRat_pow_si(epsilon, epsilon, -nprec+1);
        realRat_mul_si(epsilon, epsilon, 3*( 4*(cacheApp_getDegree(cache)) +1 ) );
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
        compApp_poly_evaluate(pc, p, c, nprec);
        arb_poly_evaluate_horner(sabsc, s, absc, nprec);
        compApp_abs(abspc, pc, nprec);
        realApp_get_mid_realApp(abspc, abspc);
        realApp_mul_realRat(sabsc, sabsc, epsilon, nprec);
        choice = realApp_soft_compare( abspc, sabsc, nprec );
        
    }
    
    if (choice==1) {
//         printf("---depth: %d, prec before: %d, prec after: %d, inclusion radius: apply thm 9\n", depth, prec, nprec);
//         printf("   res of thm 9: %d\n", tstar_inclusion_radius_th9(cache, d, prec, depth, meta) );
        res = tstar_inclusion_radius_th9(cache, d, nprec, depth, meta);
    }
    if (choice==0) {
        printf("---depth: %d, prec before: %d, prec after: %d, inclusion radius: apply thm 10\n", depth, prec, nprec);   
        res = tstar_inclusion_radius_th10(cache, d, nprec, depth, meta);
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
