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

#include "cauchy/cauchy.h"

#ifdef CCLUSTER_TIMINGS
double time_in_cauchy_discard_compBox_list;
double time_in_cauchy_bisect_connCmp;
double time_in_cauchy_MaxNumberOfRootsInCC;
double time_in_cauchy_compression;
double time_in_cauchy_newton;
double time_in_cauchy_connCmp_is_separated;
double time_in_cauchy_connCmp_is_separated_certified;
double time_in_cauchy_main_loop;
double time_in_proba_counting;
double time_in_real_coeffs;
double time_in_certified;
double time_in_insert_sorted;
int    nbPos;
#endif

int   schroeders_inclusion_test_get_rho( realApp_t rho, const compDsk_t Delta, cacheCauchy_t cacheCau, slong prec, metadatas_t meta){
    /* let Delta = D(c, r)                                */
    /*     pc(x) = p(c+rx)                                */
    /*     t(x)  = pc_reverse(x) = x^d * p(c + r/x)       */
    /* if returns 1 then rho > r1(0, t)                 */
    /* if returns 0 then r_d(0, pc) <= rb with b \simeq 2^-53 */
//     int level = 2;
    int res = 0;
    slong s = -CCLUSTER_DEFAULT_PREC;
//     slong s = -30;
    compApp_t cApp, pval, pderval;
    compApp_init(cApp);
    compApp_init(pval);
    compApp_init(pderval);
    compApp_set_compRat(cApp, compDsk_centerref(Delta), prec);
    mag_add_ui_2exp_si( arb_radref(compApp_realref(cApp)), arb_radref(compApp_realref(cApp)), 1, s );
    mag_mul_fmpz( arb_radref(compApp_realref(cApp)), arb_radref(compApp_realref(cApp)), fmpq_numref( compDsk_radiusref(Delta) ) );
    mag_div_fmpz( arb_radref(compApp_realref(cApp)), arb_radref(compApp_realref(cApp)), fmpq_denref( compDsk_radiusref(Delta) ) );
    mag_add_ui_2exp_si( arb_radref(compApp_imagref(cApp)), arb_radref(compApp_realref(cApp)), 1, s );
    mag_mul_fmpz( arb_radref(compApp_imagref(cApp)), arb_radref(compApp_imagref(cApp)), fmpq_numref( compDsk_radiusref(Delta) ) );
    mag_div_fmpz( arb_radref(compApp_imagref(cApp)), arb_radref(compApp_imagref(cApp)), fmpq_denref( compDsk_radiusref(Delta) ) );
    mag_max( arb_radref(compApp_realref(cApp)), arb_radref(compApp_realref(cApp)), arb_radref(compApp_imagref(cApp)) );
    mag_set( arb_radref(compApp_imagref(cApp)), arb_radref(compApp_realref(cApp)) );
    /* cApp = c +/- rb with b > 2^-s */
    /* eval p at cApp */
    cacheCauchy_eval( pval, pderval, cApp, 1, cacheCau, prec, meta );
    /* check if pval contains 0 */
    if ( compApp_contains_zero( pval ) ) {
        /* r_d(0, pc) <= rb : do not want to apply the test because delta probably contains a root */
        res = 0;
    } else {
        res = 1;
        arf_t rho2;
        arf_init(rho2);
        realApp_get_rad_realApp( rho, compApp_realref(cApp) );
        /* rho < r_d(0, t) */
        realApp_inv( rho, rho, prec );
        arb_get_ubound_arf(rho2, rho, prec);
        arb_set_arf(rho, rho2);
        /* rho > r_1(0, t) */
        arf_clear(rho2);
    }
    
    compApp_clear(cApp);
    compApp_clear(pval);
    compApp_clear(pderval);
    return res;
}

int   schroeders_inclusion_test( compDsk_t Delta, cacheCauchy_t cacheCau, slong prec, metadatas_t meta){
    /* assume d>=2 */
    /* let Delta = D(c, r)                                */
    /*     pc(x) = p(c+rx)                                */
    /*     t(x)  = pc_reverse(x) = x^d * p(c + r/x)       */
    /* then t'(x) = dx^(d-1) * p(c + r/x) - rx^(d-2) * p'(c + r/x) */
    int level = 2;
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---schroeders_inclusion_test: begin \n");
    }
    int res=0;
    slong d = cacheCauchy_degreeref(cacheCau);
    int l=0;
    slong q = 0x1;
    while (q<=d) {
        l++;
        q = q << 1;
    }
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---schroeders_inclusion_test: d: %ld, l:%d, q:%ld \n", d, l, q);
    }
    realApp_t rho;
    realApp_init(rho);
    res = schroeders_inclusion_test_get_rho( rho, Delta, cacheCau, prec, meta);
//     if (res==1) {
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---                               rho: "); realApp_printd(rho, 10); printf("\n");
        }
        realApp_mul_si( rho, rho, 2, prec );
        compApp_t x, y, c, pval, pderval, t, u;
        compApp_init(x);
        compApp_init(y);
        compApp_init(c);
        compApp_init(pval);
        compApp_init(pderval);
        compApp_init(t);
        compApp_init(u);
        realApp_t abs, one;
        realApp_init(abs);
        realApp_init(one);
        realApp_one(one);
        compApp_set_compRat(c, compDsk_centerref(Delta), prec);
        res = 0;
        for (slong g = 0; (g < q) && (res==0) ; g++){
            compApp_zero(x);
            /* compute x = 2\rho*e^(2*i*pi g / q) */
            realApp_set_si     ( compApp_realref(x), 2*g );
            realApp_div_ui     ( compApp_realref(x), compApp_realref(x), q, prec );
            compApp_exp_pi_i   ( x,                  x, prec);
//             compApp_mul_realApp( x,                  x, rho, prec );
            compApp_mul_si( x,                  x, 2, prec );
            /* compute y = c + r/x */
            compApp_inv( y, x, prec );
            compApp_mul_realRat( y, y, compDsk_radiusref(Delta), prec );
//             compApp_mul_si( y,                  y, 4, prec );
            compApp_add(y, c, y, prec);                                    /* y = c + r/x */
            
//             compApp_one(y);
//             compApp_mul_si( y,                  y, 40, prec );
            /* compute t(x)  = x^d * p(c + r/x) */
            /*     and t'(x) = dx^(d-1) * p(c + r/x) - rx^(d-2) * p'(c + r/x) */
            /* evaluate p and p' at y = c+r/x */
            cacheCauchy_eval( pval, pderval, y, 1, cacheCau, prec, meta );
            compApp_pow_si( t, x, d-2, prec );                             /* t <- x^(d-2) */
            compApp_mul(pderval, pderval, t, prec);                        /* pderval <- x^(d-2)*p'(c + r/x) */
            compApp_mul_realRat(pderval, pderval, compDsk_radiusref(Delta), prec ); /* pderval <- r*x^(d-2)*p'(c + r/x) */
            compApp_mul(t, t, x, prec );                                   /* t <- x^(d-1) */
            compApp_mul(u, t, pval, prec);                                 /* u <- x^(d-1)*p(c + r/x) */
            compApp_mul_si(u, u, d, prec );                                /* u <- d*x^(d-1)*p(c + r/x) */       
            compApp_sub(pderval, u, pderval, prec);                        /* pderval <- u - pderval */
            compApp_mul(t, t, x, prec );                                   /* t <- x^d */
            compApp_mul(pval, pval, t, prec );                             /* pval <- x^d*p(c + r/x) */
            /* compute one schroeder iteration: c = x - d*t(x)/t'(x) */
            compApp_mul_si(pval, pval, d, prec);
            compApp_div( pval, pval, pderval, prec );
            compApp_sub( y, x, pval, prec );
            compApp_abs( abs, y, prec);
            if (metadatas_getVerbo(meta)>=level) {
//                 printf("#---                               x  : "); compApp_printd(x, 10); printf("\n");
//                 printf("#---                               z  : "); compApp_printd(y, 10); printf("\n");
//                 printf("#---                              |z| : "); realApp_printd(abs, 10); printf("\n");
                printf("#---                     g = %ld |z|>1? : %d, |z|: ", g, realApp_gt(abs, one)); realApp_printd(abs, 10); printf("\n");
//                 realApp_set_realRat( compApp_realref(x), compDsk_radiusref(Delta), prec);
//                 printf("#---                               r  : "); realApp_printd(compApp_realref(x), 10); printf("\n");
            }
            res = realApp_gt(abs, one);
//             /* compute t2(x)  = p(c + r/x) */
//             /*     and t2'(x) = - rx^(-2) * p'(c + r/x) */
//             /* evaluate p and p' at y = c+r/x */
//             cacheCauchy_eval( pval, pderval, y, 1, cacheCau, prec, meta );
//             compApp_pow_si( t, x, -2, prec );                             /* t <- x^(-2) */
//             compApp_mul(pderval, pderval, t, prec);                        /* pderval <- x^(-2)*p'(c + r/x) */
//             compApp_mul_realRat(pderval, pderval, compDsk_radiusref(Delta), prec ); /* pderval <- r*x^(-2)*p'(c + r/x) */
//             compApp_neg(pderval, pderval); /* pderval <- -r*x^(-2)*p'(c + r/x) */
//             /* compute one schroeder iteration: y = x - d*t(x)/t'(x) */
//             compApp_mul_si(pval, pval, d, prec);
//             compApp_div( pval, pval, pderval, prec );
//             compApp_sub( y, x, pval, prec );
//             compApp_abs( abs, y, prec);
//             if (metadatas_getVerbo(meta)>=level) {
//                 printf("#---                               x  : "); compApp_printd(x, 10); printf("\n");
//                 printf("#---                               z  : "); compApp_printd(y, 10); printf("\n");
//                 printf("#---                              |z| : "); realApp_printd(abs, 10); printf("\n");
//                 printf("#---                           |z|>1? : %d\n", realApp_gt(abs, one));
//             }
//             res = realApp_gt(abs, one);
        }
        realApp_clear(one);
        realApp_clear(abs);
        compApp_clear(u);
        compApp_clear(t);
        compApp_clear(pval);
        compApp_clear(pderval);
        compApp_clear(c);
        compApp_clear(y);
        compApp_clear(x);
//     }
    realApp_clear(rho);
    if (metadatas_getVerbo(meta)>=level) {
            printf("#---schroeders_inclusion_test: end; res: %d \n", res);
    }
    return res;
}

// void  cauchy_num_lth_der( compApp_t res, const compApp_t point, int l, const realApp_t eps, 
//                           cacheCauchy_t cacheCau, slong prec, metadatas_t meta ) {
//         if (l==1) {
//             compApp_t dummy;
//             compApp_init(dummy);
//             cacheCauchy_eval( dummy, res, point, 1, cacheCau, prec, meta );
//             compApp_clear(dummy);
//         } else {
//             compApp_t pointPeps;
//             compApp_init(pointPeps);
//             realApp_add( compApp_realref(pointPeps), compApp_realref(point), eps, CCLUSTER_DEFAULT_PREC);
//             realApp_add( compApp_imagref(pointPeps), compApp_imagref(point), eps, CCLUSTER_DEFAULT_PREC);
//             cauchy_num_lth_der( res, pointPeps, l-1, eps, cacheCau, prec, meta );
//             cauchy_num_lth_der( pointPeps, point, l-1, eps, cacheCau, prec, meta );
//             compApp_sub(res, res, pointPeps, prec);
//             compApp_div_realApp(res, res, eps, prec);
//             compApp_clear(pointPeps);
//         }
// }

int   cauchy_inclusionTests( const compDsk_t Delta, cacheCauchy_t cacheCau, slong prec, metadatas_t meta ) {
    
    int res = 0;
    compApp_t c, cPe, cP2e;
    compApp_t Pc, Ppc, PpcPe, PpcP2e, Pppc, Ppppc, eps;
    compApp_init(c);
    compApp_init(cPe);
    compApp_init(cP2e);
    compApp_init(Pc);
    compApp_init(Ppc);
    compApp_init(PpcPe);
    compApp_init(PpcP2e);
    compApp_init(Pppc);
    compApp_init(Ppppc);
    compApp_init(eps);
    
    realApp_t r, rd;
    realApp_init(r);
    realApp_init(rd);
    compApp_set_compRat(c, compDsk_centerref(Delta), prec);
    realApp_set_realRat(r, compDsk_radiusref(Delta), prec);
    realApp_div_ui(compApp_realref(eps), r, 1000000, prec);
    realApp_div_ui(compApp_imagref(eps), r, 1000000, prec);
    
    slong d = cacheCauchy_degreeref(cacheCau);
    /* evaluate P and Pp at c */
    cacheCauchy_eval( Pc, Ppc, c, 1, cacheCau, prec, meta );
    /* compute rd = d|Pc/Ppc| */
    compApp_div(cPe, Pc, Ppc, prec);
    compApp_abs(rd, cPe, prec);
    realApp_mul_si(rd, rd, d, prec);
    res = realApp_le(rd, r);
    
    if (res){
            printf("INCLUSION TEST 1 IN EXCLUSION SUCCEEDS\n");
        }
        
    if (res==0) {
        /* compute Pppc numerically */
        compApp_add( cPe, c, eps, prec);
        /* evaluate Pp at c+eps */
        cacheCauchy_eval( Pppc, PpcPe, cPe, 1, cacheCau, prec, meta );
        /* compute Pppc = (Pp(c+eps) - Pp(c))/eps */
        compApp_sub( Pppc, PpcPe, Ppc, prec);
        compApp_div( Pppc, Pppc, eps, prec );
        /* compute 1 / (p'(x)/p(x))' evaluated at c := 1/ ( (p''(c)p(c) -p'(c)p'(c))/(p(c)*p(c)) ) */
        compApp_mul( c,     Pppc, Pc, prec);/* c     <- p''(c)p(c) */
        compApp_sqr( Ppppc, Ppc, prec );    /* Ppppc <- (p'(c)^2) */
        compApp_sub( c, c, Ppppc, prec);
        compApp_sqr( Ppppc,     Pc,  prec );/* Ppppc <- (p(c)^2) */
        compApp_div( c, Ppppc, c, prec );
        compApp_abs(rd, c, prec);
        realApp_mul_si(rd, rd, d, prec);
        realApp_sqrt(rd, rd, prec);
        res = realApp_le(rd, r);
        
        if (res){
            printf("INCLUSION TEST 2 IN EXCLUSION SUCCEEDS\n");
        }
    }
    
    if (res==0) {
        /* compute Ppppc numerically */
        compApp_add( cP2e, cPe, eps, prec);
        /* evaluate Pp at c+2eps */
        cacheCauchy_eval( Ppppc, PpcP2e, cP2e, 1, cacheCau, prec, meta );
        /* compute Ppppc = (PpcP2e - 2*PpcPe + Ppc)/eps^2 */
        compApp_mul_si(PpcPe, PpcPe, 2, prec);
        compApp_sub(   Ppppc, PpcP2e, PpcPe, prec );
        compApp_add(   Ppppc, Ppppc, Ppc, prec );
        compApp_sqr(   eps,   eps,        prec );
        compApp_div(   Ppppc, Ppppc, eps, prec ); /* Ppppc <- p'''(c) */
        compApp_div(   Ppppc, Ppppc, Pc,  prec ); /* Ppppc <- p'''(c)/p(c) */
        compApp_sqr(   c,     Pc,         prec ); /* c     <- p(c)^2 */
        compApp_mul_si(Pppc,  Pppc,  3,   prec ); /* Pppc  <- 3p''(c) */
        compApp_mul(   Pppc,  Pppc,  Ppc, prec ); /* Pppc  <- 3p''(c)p'(c) */
        compApp_div(   Pppc,  Pppc,  c,   prec ); /* Pppc  <- 3p''(c)p'(c) / p(c)^2 */
        compApp_pow_si(Ppc,   Ppc,   3,   prec ); /* Ppc   <- p'(c)^3 */
        compApp_mul_si(Ppc,   Ppc,   2,   prec ); /* Ppc   <- 2p'(c)^3 */
        compApp_pow_si(Pc,    Pc,    3,   prec ); /* Pc    <- p(c)^3 */
        compApp_div   (Ppc,   Ppc,   Pc,  prec ); /* Ppc   <- 2p'(c)^3/p(c)^3 */
        compApp_sub   (c,     Ppppc, Pppc,prec ); /*   c   <- */
        compApp_add   (c,     c,     Ppc, prec );
        compApp_abs   (rd,    c,          prec );
        realApp_inv   (rd,    rd,         prec );
        realApp_mul_si(rd,    rd,    d,   prec );
        realApp_root_ui(rd,   rd,    4,   prec );
        
        res = realApp_le(rd, r);
        
        if (res){
            printf("INCLUSION TEST 3 IN EXCLUSION SUCCEEDS\n");
        }
        
    }
    
    realApp_clear(r);
    realApp_clear(rd);
    compApp_clear(eps);
    compApp_clear(c);
    compApp_clear(cPe);
    compApp_clear(cP2e);
    compApp_clear(Pc);
    compApp_clear(Ppc);
    compApp_clear(PpcPe);
    compApp_clear(PpcP2e);
    compApp_clear(Pppc);
    compApp_clear(Ppppc);
    
    return res;
}

slong cauchy_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
                                     cacheCauchy_t cacheCau,
                                     slong nbMSols, 
                                     slong prec, metadatas_t meta){
#ifdef CCLUSTER_TIMINGS    
    clock_t start=clock();
#endif
    
//     int level = 4;
    
    /* for test */
//     nbMSols = -1;

//     tstar_res res;
//     res.appPrec = prec;
    cauchyTest_res res;
    if ( (metadatas_usefpri(meta)) && (cacheCauchy_evalFastref(cacheCau)==NULL) )
        res.appPrec = CCLUSTER_FPRI_PREC;
    else
        res.appPrec = CCLUSTER_DEFAULT_PREC;
        
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
//     compRat_t compDist;
//     compRat_init(compDist);
    
    compBox_ptr btemp;
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
        /* Real Coeffs */
        if (( metadatas_useRealCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
//             printf("ici\n");
//             if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
//             } else {
//                 compBox_list_push(bDiscarded, btemp);
//             }
            continue;
        }

// #ifdef CAUCHY_INCLUSION
//         int successTest = cauchy_inclusionTests( bdisk, cacheCau, 8*CCLUSTER_DEFAULT_PREC, meta );
// #endif

        /* exclusion test */

        res = cauchyTest_probabilistic_exclusion_test( compDsk_centerref(bdisk), compDsk_radiusref(bdisk), 
                                                            NULL, 0,0,
                                                            cacheCau, res.appPrec, meta, depth);
//         if (metadatas_getVerbo(meta)>=level) {
//             printf("#---cauchy exclusion test: res: %d\n", res.nbOfSol );
//         }
// #ifdef CAUCHY_INCLUSION
//         if (successTest){
//             printf("                 res: %d\n", res.nbOfSol);
//         }
// #endif  
        if (res.nbOfSol>0){
            printf("#---cauchy exclusion test: res: %d\n", res.nbOfSol );
            int successTest3 = schroeders_inclusion_test( bdisk, cacheCau, 8*res.appPrec, meta);
            printf("#---Schroeder's inclusion test, res: %d\n", successTest3 );
        }
#ifdef CCLUSTER_TIMINGS
        if (res.nbOfSol>0){
            nbPos++;
        }
#endif
        if (res.nbOfSol==0) {
            if (metadatas_haveToCount(meta)){
                metadatas_add_discarded( meta, depth);
            }
            if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
            } else {
                compBox_list_push(bDiscarded, btemp);
            }
        }
        
        else{
                if (res.nbOfSol>0) {
                    btemp->nbMSol = res.nbOfSol;
                }
                
                /* there is at least a root in 2bdisk */
//                 realRat_mul_si( compDsk_radiusref(bdisk), compDsk_radiusref(bdisk), 2);
//                 compBox_ptr btemp2;
//                 compBox_list_t ltemp2;
//                 compBox_list_init(ltemp2);
//                 while (!compBox_list_is_empty(boxes)){
//                     btemp2 = compBox_list_pop(boxes);
//                     /* check if btemp2 is 4-adjacent to btemp */
//                     compRat_comp_distance( compDist, compBox_centerref(btemp), compBox_centerref(btemp2) );
//                     if ( ( ( realRat_cmp( compRat_realref(compDist), compBox_bwidthref(btemp) ) ==0 ) &&
//                            ( realRat_is_zero( compRat_imagref(compDist) ) ) ) ||
//                          ( ( realRat_cmp( compRat_imagref(compDist), compBox_bwidthref(btemp) ) ==0 ) &&
//                            ( realRat_is_zero( compRat_realref(compDist) ) ) )  ) {
//                         compBox_list_push(ltemp, btemp2);
//                     } else {
//                         compBox_list_push(ltemp2, btemp2);
//                     }
//                 }
//                 compBox_list_swap(boxes, ltemp2);
//                 compBox_list_clear(ltemp2);
                
                if (nbMSols==1) {
                    /* there is at least a root in 2bdisk */
                    realRat_mul_si( compDsk_radiusref(bdisk), compDsk_radiusref(bdisk), 2);
                    compBox_ptr btemp2;
                    compBox_list_t ltemp2;
                    compBox_list_init(ltemp2);
                    while (!compBox_list_is_empty(boxes)){
                        btemp2 = compBox_list_pop(boxes);
                        /* test if btemp2 intersects 2bdisk */
                        if (compBox_intersection_is_not_empty_compDsk ( btemp2, bdisk)){
                            /* keep btemp2 */
                            compBox_list_push(ltemp2, btemp2);
                        } else {
//                             printf("##############################ICI##########################\n");
                            /* discard btemp2 */
                            if (metadatas_haveToCount(meta)){
                                metadatas_add_discarded( meta, depth);
                            }
                            if (metadatas_getDrSub(meta)==0){
                                compBox_clear(btemp2);
                                ccluster_free(btemp2);
                            } else {
                                compBox_list_push(bDiscarded, btemp2);
                            }
                        }
                    }
                    compBox_list_swap(boxes, ltemp2);
                    compBox_list_clear(ltemp2);
                }
                compBox_list_push(ltemp, btemp);
        }
    }
    
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
//     compRat_clear(compDist);
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_discard_compBox_list += (double) (clock() - start);
#endif    
    return res.appPrec;
}

void cauchy_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              compBox_list_t bDiscarded,
                              cacheCauchy_t cacheCau,
                              slong nbMsols, 
                              metadatas_t meta, 
                              slong nbThreads){
#ifdef CCLUSTER_TIMINGS     
    clock_t start=clock();
#endif
    
    slong prec = connCmp_appPr(cc);
    compBox_list_t subBoxes;
    connCmp_list_t ltemp;
    compBox_list_init(subBoxes);
    connCmp_list_init(ltemp);
    
    compBox_ptr btemp;
    connCmp_ptr ctemp;
    
    /* RealCoeffs */
    int cc_contains_real_line = 0;
    /* Check if cc contains the real line */
    if ( (metadatas_useRealCoeffs(meta)) && (!connCmp_is_imaginary_positive(cc)) )
        cc_contains_real_line = 1;
    /* end RealCoeffs */
    
    while (!connCmp_is_empty(cc)) {
        btemp = connCmp_pop(cc);
        subdBox_quadrisect( subBoxes, btemp );
        compBox_clear(btemp);
        ccluster_free(btemp);
    }
    
    prec = cauchy_discard_compBox_list( subBoxes, bDiscarded, cacheCau, nbMsols, prec, meta);
    
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    
    int specialFlag = 1;
    if (connCmp_list_get_size(ltemp) == 1)
        specialFlag = 0;
    
    /* RealCoeffs */
    if ( (metadatas_useRealCoeffs(meta)) && (connCmp_list_get_size(ltemp) == 1) && (cc_contains_real_line == 1) ){
        ctemp = connCmp_list_first(ltemp);
        /* test if cc has been separated from real case;
         in which case reset everything*/
        if ( connCmp_is_imaginary_positive(ctemp) )
            specialFlag = 1;
    }
    /* end RealCoeffs */
    
    slong nprec; 
    if (prec == connCmp_appPrref(cc)) {
        nprec = CCLUSTER_MAX(prec/2,CCLUSTER_DEFAULT_PREC);
//         printf("decrease precision\n");
    }
    else 
        nprec = prec;
    
    
    while (!connCmp_list_is_empty(ltemp)){
        ctemp = connCmp_list_pop(ltemp);
        
        if (connCmp_intersection_is_not_empty(ctemp, metadatas_initBref(meta))){
            connCmp_appPrref(ctemp) = nprec;
            if (specialFlag)
                connCmp_initiali_nwSpd(ctemp);
            else {
                connCmp_initiali_nwSpd_connCmp(ctemp, cc);
                connCmp_decrease_nwSpd(ctemp);
                /* copy the number of sols */
                connCmp_nSolsref(ctemp) = connCmp_nSolsref(cc);
                connCmp_isSep(ctemp) = connCmp_isSep(cc);
                connCmp_isSepCertref(ctemp) = connCmp_isSepCertref(cc);
                connCmp_isRigref(ctemp) = connCmp_isRigref(cc);
                compDsk_set( connCmp_grSepDiskref(ctemp), connCmp_grSepDiskref(cc) );
            }
            connCmp_list_push(dest, ctemp);
        }
        else {
            connCmp_appPrref(ctemp) = prec;
            connCmp_list_push(discardedCcs, ctemp);
        }
    }
    
    compBox_list_clear(subBoxes);
    connCmp_list_clear(ltemp);

#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_bisect_connCmp += (double) (clock() - start);
#endif
    
}

void cauchy_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t qPrepLoop, 
                         connCmp_list_t discardedCcs,
                         cacheCauchy_t cacheCau,
                         metadatas_t meta) {
    
    connCmp_list_t ltemp;
    realRat_t halfwidth, diam;
    connCmp_list_init(ltemp);
    realRat_init(halfwidth);
    realRat_init(diam);
    
    realRat_set_si(halfwidth, 1, 2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(metadatas_initBref(meta)) );
    
    connCmp_ptr ctemp;
    
    while (!connCmp_list_is_empty(qPrepLoop)) {
        
        ctemp = connCmp_list_pop(qPrepLoop);
        connCmp_diameter(diam, ctemp);
        
        if ( connCmp_is_confined(ctemp, metadatas_initBref(meta)) && (realRat_cmp(diam, halfwidth)<0) )
            connCmp_list_insert_sorted(qMainLoop, ctemp);
        else {
            cauchy_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cacheCau, -1, meta, metadatas_useNBThreads(meta));
            
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ctemp);
            ccluster_free(ctemp);
        }        
    }
    
    connCmp_list_clear(ltemp);
    realRat_clear(halfwidth);
    realRat_clear(diam);
}

/* try to upper bound the number of roots in ccur */
/* at least one sol per connected comp in qMainLoop, */
/* 2 if use real coeffs and the CC does not contain the real line */
slong cauchy_MaxNumberOfRootsInCC( const connCmp_t ccur,
                                   connCmp_list_t qMainLoop,
                                   slong nbSolsInQResults,
                                   cacheCauchy_t cacheCau, 
                                   metadatas_t meta ){
    
#ifdef CCLUSTER_TIMINGS    
    clock_t start = clock();
#endif
    
    slong nbMaxSol=connCmp_nSolsref(ccur);
    if (nbMaxSol==-1) {
        nbMaxSol=0;
        connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
        while (it!=connCmp_list_end()){
            
            if (connCmp_nSolsref( connCmp_list_elmt( it ) ) > -1) {
                nbMaxSol += connCmp_nSolsref( connCmp_list_elmt( it ) ) ;
                if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(connCmp_list_elmt( it ))) )
                    nbMaxSol += connCmp_nSolsref( connCmp_list_elmt( it ) ) ; 
            } else {
                nbMaxSol += 1 ;
                if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(connCmp_list_elmt( it ))) )
                    nbMaxSol += 1 ; 
            }
                
            it = connCmp_list_next(it);
        }
        nbMaxSol += nbSolsInQResults;
//         printf("nbMaxSol: %ld\n", nbMaxSol);
        nbMaxSol = cacheCauchy_degreeref(cacheCau) - nbMaxSol;
        if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(ccur)) )
            nbMaxSol = nbMaxSol/2;
    }
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_MaxNumberOfRootsInCC += (double) (clock() - start);
#endif
    
    return nbMaxSol;
}

/* let ccDisk = D(c,r) be at least 2-isolated and contain m roots and ccur*/
/* let epsp = max( eps, r/( (1 + 3d/m)^2 ) ) */
/* compute a CC with containing disk D(c',r') s.t   */
/* either D(c',r') is rigid with radius r' > epsp */
/* or     D(c',r')           has radius r' <= epsp */
/* in the latter case D(c', (1 + 3d/m)r') is (1 + 3d/m)-isolated */
connCmp_ptr cauchy_compression( connCmp_ptr ccur,
                                int *rigidFlag,
                                int *widthFlag,
                                int *isolaFlag,
                                const compDsk_t ccDisk,
                                slong m,
                                const realRat_t eps,
                                cacheCauchy_t cacheCau,
                                slong prec,
                                metadatas_t meta,
                                slong depth ) {

#ifdef CCLUSTER_TIMINGS
    clock_t start = clock();
#endif    
    
    int level = 3;
    
    /* let epsp = max( eps, r/( (1 + 3d/m)^2 ) ) */
    realRat_t epsp;
    realRat_init(epsp);
    realRat_set_si( epsp, (slong) connCmp_nSolsref(ccur), (slong) connCmp_nSolsref(ccur) + 3*cacheCauchy_degreeref(cacheCau) );
    realRat_mul( epsp, epsp, epsp );
    realRat_mul( epsp, epsp, compDsk_radiusref(ccDisk) );
    /* case where eps = +inf; eps = 1/0 */
    if (realRat_is_den_zero(eps)==0)
        realRat_max_2_realRat(epsp, eps);
    
    if (!(metadatas_useNewton(meta))) {
        if (realRat_is_den_zero(eps)==0)
            realRat_set(epsp, eps);
        else 
            realRat_one(epsp);
    }
    
    /* let epspp = (2/9) epsp */
    /* compute a disk D(c'',r'') s.t. */
    /* either D(c'',r'') is rigid with radius r'' > epspp = (2/9) epsp */
    /*     or D(c'',r'')           has radius r'' <= epspp = (2/9) epsp */
    realRat_t epspp;
    realRat_init(epspp);
    realRat_mul_si(epspp, epsp, 2);
    realRat_div_ui(epspp, epspp, 9);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("\n#------Test compression into rigid disc for a CC with depth %ld with %d roots \n", depth, connCmp_nSolsref(ccur));
        printf("#---------Delta: "); compDsk_print( ccDisk ); printf("\n");
        printf("#---------Required radius of containing disk of CC: "); realRat_print(epsp); printf("\n");
        printf("#---------Required radius of disk:                  "); realRat_print(epspp); printf("\n");
        printf("#---------Epsilon                :                  "); realRat_print(eps); printf("\n");
    }
    
    compDsk_t res;
    compDsk_init(res);
    /* ccDisk is at least 2 isolated */
    realRat_t theta;
    realRat_init(theta);
    realRat_set_si(theta, 2, 1);
    
    slong precres = cauchy_compressionIntoRigidDisk( res, ccDisk, m, theta, epspp, cacheCau, prec, meta, depth);
    
//     compDsk_t twoCCDisk;
//     compDsk_init(twoCCDisk);
//     compDsk_set(twoCCDisk, res);
//     realRat_mul_si( compDsk_radiusref(twoCCDisk), compDsk_radiusref(res), 2);
//     cauchyTest_res resCauchy = cauchyTest_probabilistic_counting( twoCCDisk, cacheCau, prec, meta, depth);
//     printf("m: %ld, resCauchy.nbOfSol: %d\n", m, resCauchy.nbOfSol);
//     if ( resCauchy.nbOfSol == 0 ) {
//         printf("m: %ld, resCauchy.nbOfSol: %d\n", m, resCauchy.nbOfSol);
//         printf("#FAILURE: A DISK is not certified to contain at least one root \n");
//     }
//     compDsk_clear(twoCCDisk);
//     printf("#---------Precision after compression: %ld\n", precres);
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---------Precision after compression: %ld\n", precres);
        printf("#---------res: "); compDsk_print( res ); printf("\n");
        compApp_t c;
        compApp_init(c);
        compApp_set_compRat(c, compDsk_centerref(res), precres);
        printf("#---------app: center: "); compApp_printd( c, 10 );
        realApp_set_realRat(compApp_realref(c), compDsk_radiusref(res), precres);
        printf("# radius: "); realApp_printd( compApp_realref(c), 10 ); printf("\n");
        compApp_clear(c);
    }
    
    /* actualize connected component: */
    if ( realRat_cmp( compDsk_radiusref(res), epspp) > 0) { /* if r'' > epspp = (2/9)epsp */
        /* cover D(c'',r'') with at most 9 sub-boxes of the subdivision tree of width >r'' and <= 2r'' */
        /* resulting CC has width at most 6*r'' */
        /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*6*r'' = (9/2)r'' */
        ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
        *rigidFlag = 1;
        connCmp_isRigref(ccur) = 1;
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------rigid disk\n");
        }
    } else {
        if ( realRat_cmp( epsp, eps) == 0) { /* r'' <= (2/9)eps */
            /* cover D(c'',r'') with at most 9 sub-boxes of the subdivision tree of width >r'' and <= 2r'' */
            /* resulting CC has width at most 6*r'' <= 6*(2/9)eps */
            /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*6*(2/9)eps = eps */
            ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
            *widthFlag = 1;
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------disk with size smaller than epsilon\n");
            }
    
        } else {
            /* the Disk D(c'', (1 + 3d/m)r'') is (1 + 3d/m) isolated */
            realRat_mul_si( compDsk_radiusref(res), compDsk_radiusref(res), 
                            connCmp_nSolsref(ccur) + 3*cacheCauchy_degreeref(cacheCau) );
            realRat_div_ui( compDsk_radiusref(res), compDsk_radiusref(res), 
                            connCmp_nSolsref(ccur) );
            /* cover D(c'',(1 + 3d/m)r'') with at most 9 sub-boxes of the subdivision tree */
            /* of size (1 + 3d/m)r''=(1/3)(1 + 3d/m)epsp */
            /* resulting CC has width at most 3*(1/3)(1 + 3d/m)*epsp = (1 + 3d/m)*epsp */
            /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*(1 + 3d/m)*epsp */
            /* and is (1 + 3d/m) isolated */
            ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
            *isolaFlag = 1;
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------isolated disk\n");
            }
        }
        
    }
    
    realRat_clear(theta);
    compDsk_clear(res);
    realRat_clear(epsp);
    realRat_clear(epspp);
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_compression += (double) (clock() - start);
#endif
    
    return ccur;
}

newton_res cauchy_newton ( connCmp_ptr * ccur,
                           const compBox_t componentBox,
                           cacheCauchy_t cacheCau,
                           slong prec,
                           metadatas_t meta){

#ifdef CCLUSTER_TIMINGS
    clock_t start = clock();
#endif 
    
    newton_res resNewton;
    
    compRat_t initPoint;
    compRat_init (initPoint);
//     if (connCmp_nSols(*ccur)==1) 
//         compRat_set(initPoint, compBox_centerref(componentBox));
//     else
        connCmp_find_point_outside_connCmp( initPoint, *ccur, metadatas_initBref(meta) );
    
    connCmp_ptr nCC;
    nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init(nCC);
    resNewton = newton_cauchy_newton_connCmp( nCC, *ccur, cacheCau, initPoint, prec, meta);
            
    compRat_clear(initPoint);
    
    if (resNewton.nflag) {
        connCmp_clear(*ccur);
        ccluster_free(*ccur);
        *ccur = nCC;
        connCmp_increase_nwSpd(*ccur);
        connCmp_newSuref(*ccur) = 1;
        connCmp_appPrref(*ccur) = resNewton.appPrec;
        
        /* it was separated; 8 time the new cc disk is included in four times the old one */
//         connCmp_isSepCertref(*ccur) = 1; -> can not say about separated from 2*qMainLoop!
    
    }
    else {
        connCmp_newSuref(*ccur) = 0;
        connCmp_clear(nCC);
        ccluster_free(nCC);
    }
       
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_newton += (double) (clock() - start);
#endif
    
    return resNewton;
}

void cauchy_conjugate ( connCmp_ptr * ccurConj, int * pushConjugFlag, int * separationFlag, 
                        connCmp_ptr ccur, metadatas_t meta ){
    
    if (connCmp_is_imaginary_positive(ccur)) {
        * pushConjugFlag = 1;
        /*compute the complex conjugate*/
        * ccurConj = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
        connCmp_init( *ccurConj );
        connCmp_set_conjugate(*ccurConj, ccur);
        
        /* test if initial box is symetric relatively to real axe */
        if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
            /* test if the cc intersects initial box */
            if ( connCmp_intersection_is_not_empty(*ccurConj, metadatas_initBref(meta)) ) {
                /* test if the cc is confined */
                if (!connCmp_is_confined(*ccurConj, metadatas_initBref(meta))) {
                    *pushConjugFlag = 0;
                    *separationFlag = 0; 
                  /* delete ccurConj*/
                    connCmp_clear(*ccurConj);
                    ccluster_free(*ccurConj);
                }
            }
            else {
                *pushConjugFlag = 0;
                /* delete ccurConj*/
                connCmp_clear(*ccurConj);
                ccluster_free(*ccurConj);
            }
        } 
    } else {
        /* test if initial box is symetric relatively to real axe */
        if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
            /* test if the cc is confined and intersects initial box */
            if (! ( connCmp_is_confined(ccur, metadatas_initBref(meta)) 
                 && connCmp_intersection_is_not_empty(ccur, metadatas_initBref(meta)) ) ){
                /* bisect ccur until this hold */
                *separationFlag = 0;
            }
        }
    }
}

connCmp_ptr cauchy_connCmp_set_conjugate_closure(connCmp_ptr cc, metadatas_t meta){
    
    connCmp_ptr ccurConjClo = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
    connCmp_init( ccurConjClo );
    connCmp_set_conjugate_closure(ccurConjClo, cc, metadatas_initBref(meta));
    connCmp_clear(cc);
    ccluster_free(cc);
    return ccurConjClo;
                
}

int  cauchy_compDsk_is_separated( const compDsk_t d, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs ){
    int res = 1;
    connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    it = connCmp_list_begin(discardedCcs);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    return res;
}

int cauchy_connCmp_is_separated( connCmp_t cc, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, metadatas_t meta ) {
    
#ifdef CCLUSTER_TIMINGS
    clock_t start = clock();
#endif
    
    if (connCmp_isSepref(cc) == 1)
        return 1;
    
    compBox_t componentBox;
    compDsk_t ccDisk;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    connCmp_componentBox(componentBox, cc, metadatas_initBref(meta));
    compBox_get_containing_dsk(ccDisk, componentBox);
    realRat_mul_si( compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 4);
    
    int res = cauchy_compDsk_is_separated( ccDisk, qMainLoop, discardedCcs );
    
    if (res && metadatas_useRealCoeffs(meta) && connCmp_is_imaginary_positive(cc) ) {
        /* check if the cc is separated from its complex conjugate */
        realRat_neg( compRat_imagref(compDsk_centerref(ccDisk)), compRat_imagref(compDsk_centerref(ccDisk)) );
        res = res&&(!connCmp_intersection_is_not_empty_compDsk( cc, ccDisk));
        /* in order to save the containing disk later */
        realRat_neg( compRat_imagref(compDsk_centerref(ccDisk)), compRat_imagref(compDsk_centerref(ccDisk)) );
    }
    
    if ( (!connCmp_isSepref(cc)) && (res==1) ) {
        compBox_get_containing_dsk(ccDisk, componentBox);
        realRat_mul_si( compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 2);
        compDsk_set( connCmp_grSepDiskref(cc) , ccDisk );
    }
    
    connCmp_isSepref(cc) = res;
        
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_connCmp_is_separated += (double) (clock() - start);
#endif
    
    return res;
}

// int cauchy_connCmp_is_separated_certified( int *separated_certified, connCmp_t cc, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, metadatas_t meta ) {
int cauchy_connCmp_is_separated_certified( connCmp_t cc, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, metadatas_t meta ) {
    
    if (connCmp_isSepCertref(cc) == 1)
        return 1;
    
#ifdef CCLUSTER_TIMINGS
    clock_t start = clock();
#endif
    
    compBox_t componentBox;
    compDsk_t sixCCDisk, twoCCDisk, otherCCDisk;
    compBox_init(componentBox);
    compDsk_init(sixCCDisk);
    compDsk_init(twoCCDisk);
    compDsk_init(otherCCDisk);
    connCmp_componentBox(componentBox, cc, metadatas_initBref(meta));
    compBox_get_containing_dsk(twoCCDisk, componentBox);
    realRat_mul_si( compDsk_radiusref(twoCCDisk), compDsk_radiusref(twoCCDisk), 2);
    compDsk_set(sixCCDisk, twoCCDisk);
    realRat_mul_si( compDsk_radiusref(sixCCDisk), compDsk_radiusref(twoCCDisk), 3);
    
//     printf(" cauchy_connCmp_is_separated_certified: \n");
//     printf(" --- greatest separated containing disk       : "); compDsk_print( connCmp_grSepDiskref(cc) ); printf("\n");
//     printf(" --- six time containing disk                 : "); compDsk_print( sixCCDisk ); printf("\n");
//     printf(" --- is sixCCDisk in connCmp_grSepDiskref(cc) : %d\n", compDsk_is_stricly_in(sixCCDisk, connCmp_grSepDiskref(cc) ) );
    
    int res = 1;
    
    if (!compDsk_is_stricly_in(sixCCDisk, connCmp_grSepDiskref(cc) )) {
//     *separated_certified = 1;
//     *sep_certified = 1;
        connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
        while ( res && (it!=connCmp_list_end()) ) {
            connCmp_componentBox(componentBox, connCmp_list_elmt(it), metadatas_initBref(meta));
            compBox_get_containing_dsk(otherCCDisk, componentBox);
            realRat_mul_si( compDsk_radiusref(otherCCDisk), compDsk_radiusref(otherCCDisk), 2);
            res = res && (! compDsk_intersect_compDsk(sixCCDisk, otherCCDisk) );
            /* test */
    //         realRat_mul_si( compDsk_radiusref(otherCCDisk), compDsk_radiusref(otherCCDisk), 3);
    //         res = res && (! compDsk_intersect_compDsk(twoCCDisk, otherCCDisk) );
            /* end test */
            it = connCmp_list_next(it);
        }
        it = connCmp_list_begin(discardedCcs);
        while ( res && (it!=connCmp_list_end()) ) {
            connCmp_componentBox(componentBox, connCmp_list_elmt(it), metadatas_initBref(meta));
            compBox_get_containing_dsk(otherCCDisk, componentBox);
            realRat_mul_si( compDsk_radiusref(otherCCDisk), compDsk_radiusref(otherCCDisk), 2);
            res = res && (! compDsk_intersect_compDsk(sixCCDisk, otherCCDisk) );
            it = connCmp_list_next(it);
        }
        
        if (res && metadatas_useRealCoeffs(meta) && connCmp_is_imaginary_positive(cc) ) {
            /* check if the cc is separated from its complex conjugate */
            connCmp_componentBox(componentBox, cc, metadatas_initBref(meta));
            compBox_get_containing_dsk(otherCCDisk, componentBox);
            realRat_mul_si( compDsk_radiusref(otherCCDisk), compDsk_radiusref(otherCCDisk), 2);
            
            realRat_neg( compRat_imagref(compDsk_centerref(sixCCDisk)), compRat_imagref(compDsk_centerref(sixCCDisk)) );
            res = res && (! compDsk_intersect_compDsk(sixCCDisk, otherCCDisk) );
        }
    
    }
    
//     connCmp_isSepCertref(cc) = res && (*separated_certified);
    connCmp_isSepCertref(cc) = res;
        
    compBox_clear(componentBox);
    compDsk_clear(sixCCDisk);
    compDsk_clear(twoCCDisk);
    compDsk_clear(otherCCDisk);
    
#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_connCmp_is_separated_certified += (double) (clock() - start);
#endif
    
//     if ( (res==0) && compDsk_is_stricly_in(sixCCDisk, connCmp_grSepDiskref(cc) ) )
//         printf(" cauchy_connCmp_is_separated_certified ERROR    : %d\n", res);
    
    return res;
}

int cauchy_main_loop( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps,
                         cacheCauchy_t cacheCau,
                         int certified,
                         metadatas_t meta){
 
#ifdef CCLUSTER_TIMINGS
    clock_t startt = clock();
#endif
    
    /* for prints */
    int level = 3;
    
    int separationFlag;
//     int separationFlag_certified;
//     int do_not_bisect_flag=0;
    int widthFlag;
    int compactFlag;
    int rigidFlag; /* the disk is m/( (2m-2)theta )-rigid */
    int isolaFlag; /* the disk is (1 + 3d/m)-isolated */
    
    slong prec, depth;
//     tstar_res resTstar;
    newton_res resNewton;
    slong nbSolsInQResults = 0;
    
    compBox_t componentBox;
    compDsk_t ccDisk, twoCCDisk, fourCCDisk;
    realRat_t two, three, four, mp2, threeWidth, neps;
    connCmp_list_t ltemp;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    compDsk_init(twoCCDisk);
    compDsk_init(fourCCDisk);
    realRat_init(two);
    realRat_init(three);
    realRat_init(four);
    realRat_init(mp2);
    realRat_init(neps);
    realRat_init(threeWidth);
    connCmp_list_init(ltemp);
    
    connCmp_ptr ccur;
    
    clock_t start=clock();
    
    /* Real Coeff */
    connCmp_ptr ccurConj;
    ccurConj = NULL;
    int pushConjugFlag = 0;
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    realRat_set_si(two, 2, 1);
    realRat_set(neps, eps);
    if (certified) 
        realRat_div_ui(neps, eps, 2);
    
    int failure = 0;
    
    connCmp_list_iterator it = connCmp_list_end();
    
    while ( (failure==0)&&(!connCmp_list_is_empty(qMainLoop)) ) {
        
        //         if (metadatas_getVerbo(meta)>0) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }
        widthFlag      = 0;
        compactFlag    = 0;
        rigidFlag      = 0;
        isolaFlag      = 0;
        resNewton.nflag = 0;
//         do_not_bisect_flag = 0;
//         separationFlag_certified = 0;
        
        if ( (!(it==NULL))&&(it==connCmp_list_begin(qMainLoop))) {
            it=NULL;
        }
        
        ccur = connCmp_list_pop(qMainLoop);
        
        /* try to upper bound the number of roots in ccur */
        /* at least one sol per connected comp in qMainLoop, */
        /* 2 if use real coeffs and the CC does not contain the real line */
        slong nbMaxSol= cauchy_MaxNumberOfRootsInCC( ccur, qMainLoop, nbSolsInQResults, cacheCau,  meta );
        
#ifdef CCLUSTER_TIMINGS
        clock_t startt3 = clock();
#endif
        /* Real Coeff */
        pushConjugFlag = 0;
        if (metadatas_useRealCoeffs(meta)) {
            /* if the component contains the real line in its interior */
            if (!connCmp_is_imaginary_positive(ccur))
                ccur = cauchy_connCmp_set_conjugate_closure( ccur, meta);
        }
#ifdef CCLUSTER_TIMINGS    
        time_in_real_coeffs += (double) (clock() - startt3);
#endif        
        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
        compBox_get_containing_dsk(ccDisk, componentBox);
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        separationFlag = cauchy_connCmp_is_separated( ccur, qMainLoop, discardedCcs, meta );
        
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), neps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        rigidFlag      = 0;
        isolaFlag      = 0;
        
//         printf("#------number of boxes in ccur:          %d\n", connCmp_nb_boxes(ccur));
//         printf("#------connCmp_nSolsref(ccur):           %d\n", connCmp_nSolsref(ccur));
//         printf("#------max number of roots in ccur:      %ld\n\n", nbMaxSol); 
            
        if (metadatas_getVerbo(meta)>=level) {
            printf("\n#---depth: %d\n", (int) depth);
//             printf("------component Box:"); compBox_print(componentBox); printf("\n");
            printf("#------length of working queue :         %d\n", connCmp_list_get_size(qMainLoop));
            printf("#------length of results queue :         %d\n", connCmp_list_get_size(qResults));
            printf("#------number of boxes in ccur:          %d\n", connCmp_nb_boxes(ccur));
            printf("#------connCmp_nSolsref(ccur):           %d\n", connCmp_nSolsref(ccur)); 
            if (connCmp_nSolsref(ccur) == -1)
                printf("#------max number of roots in ccur:      %ld\n", nbMaxSol); 
            printf("#------separation Flag: %d\n", separationFlag);
            printf("#------widthFlag: %d\n", widthFlag); 
            printf("#------compactFlag: %d\n", compactFlag);
            printf("#------rigidFlag: %d\n", rigidFlag); 
            printf("#------isolaFlag: %d\n", isolaFlag); 
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
            if (connCmp_nSolsref(ccur)==-1){

                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------run Cauchy probabilistic counter:\n");
                    }
                    
#ifdef CCLUSTER_TIMINGS
                    clock_t startt2 = clock();
#endif
                    compDsk_inflate_realRat(twoCCDisk, ccDisk, two);
                    cauchyTest_res resCauchy = cauchyTest_probabilistic_counting( twoCCDisk, cacheCau, prec, meta, depth);
                    connCmp_nSolsref(ccur) = resCauchy.nbOfSol;
                    prec = resCauchy.appPrec;
                    slong m = connCmp_nSols(ccur);
                    
//                     if (metadatas_getVerbo(meta)>=level) {
//                         tstar_res resTstar = tstar_interface( cacheCauchy_cacheAppref(cacheCau), ccDisk, cacheApp_getDegree(cacheCauchy_cacheAppref(cacheCau)), 0, 0, prec, depth, NULL, meta);
//                         printf("#------resTstar: %d\n", (int) resTstar.nbOfSol);
//                     }
                    
#ifdef CCLUSTER_TIMINGS    
                    time_in_proba_counting += (double) (clock() - startt2);
#endif
    
                    if (resCauchy.nbOfSol == -1) {
                        failure = 1;
                        continue;
                        printf("#FAILURE: A DISC IS NOT 2-ISOLATED \n");
                    }                   
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------nb sols: %d\n", (int) connCmp_nSolsref(ccur));
                    }
                    
                    if (metadatas_useCompression(meta)){
                        clock_t start2 = clock();
                        
                        if (metadatas_getVerbo(meta)>=level) {
                            printf("\n#---Compression into rigid disc for a CC with depth %ld with %d roots \n", depth, connCmp_nSolsref(ccur));
                            printf("#------Disk: "); compDsk_print(twoCCDisk); printf("\n");
                        }
                        
                        realRat_t widthSave;
                        realRat_init(widthSave);
                        realRat_set(widthSave, connCmp_widthref(ccur));
                        
                        ccur = cauchy_compression( ccur, &rigidFlag, &widthFlag, &isolaFlag, 
                                                   twoCCDisk, m, neps, cacheCau, prec, meta, depth );
                        
                        if (metadatas_getVerbo(meta)>=level) {
                            printf("#------time spent in compression: %f\n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
                            printf("\n");
                        }
                        
                        metadatas_add_time_CompTot(meta, (double) (clock() - start2));
                        
                        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
                        compBox_get_containing_dsk(ccDisk, componentBox);
                        realRat_mul(threeWidth, three, connCmp_widthref(ccur)); /* not used???*/
                        
                        /* test if component has been reduced */
                        if ( realRat_cmp( connCmp_widthref(ccur), widthSave ) < 0 ) {
#ifdef CCLUSTER_TIMINGS
                            startt3 = clock();
#endif            
//                          connCmp_list_insert_sorted(qMainLoop, ccur);
                            it = connCmp_list_insert_sorted_from_end(qMainLoop, ccur, it);
#ifdef CCLUSTER_TIMINGS    
                            time_in_insert_sorted += (double) (clock() - startt3);
#endif
                            continue;
                        }
                        
                    }
                    
                
            }
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
               ( (!widthFlag)||( connCmp_nSols(ccur)== cacheCauchy_degreeref(cacheCau) ) )
               && ( (!metadatas_useCompression(meta))||(connCmp_isRigref(ccur)==0)) ) ) {
            
            if (metadatas_getVerbo(meta)>=level)
                printf("#------run Newton:\n");
        
            start = clock();
        
            resNewton = cauchy_newton ( &ccur, componentBox, cacheCau, prec, meta);
            
            if (metadatas_getVerbo(meta)>=level)
                printf("#---------res_newton: %d \n", resNewton.nflag);
            
            if (metadatas_haveToCount(meta)){
                metadatas_add_Newton   ( meta, depth, resNewton.nflag, (double) (clock() - start) );
            }
            
            if (resNewton.nflag) {
#ifdef CCLUSTER_TIMINGS
                startt3 = clock();
#endif            
//              connCmp_list_insert_sorted(qMainLoop, ccur);
                it = connCmp_list_insert_sorted_from_end(qMainLoop, ccur, it);
#ifdef CCLUSTER_TIMINGS    
                time_in_insert_sorted += (double) (clock() - startt3);
#endif
                continue;                
            }
                
        }
        
        if (certified && (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag) {
            separationFlag = cauchy_connCmp_is_separated_certified( ccur, qMainLoop, discardedCcs, meta );
//             separationFlag = cauchy_connCmp_is_separated_certified( &separationFlag_certified, ccur, qMainLoop, discardedCcs, meta );
            if (metadatas_getVerbo(meta)>=level)
                printf("#------is separated Certified: %d\n", separationFlag );
//             if (metadatas_useNewton(meta) && certified) {
//                 separationFlag = separationFlag && separationFlag_certified;
//                 separationFlag_certified = separationFlag;
//                 do_not_bisect_flag = 0;
//             } else {
//                 if ( certified && separationFlag && !separationFlag_certified )
//                     do_not_bisect_flag = 1;
//             }
        }
        
#ifdef CCLUSTER_TIMINGS
        startt3 = clock();
#endif
        /* Real Coeff */
        if ( metadatas_useRealCoeffs(meta)
            && ( ( (connCmp_nSols(ccur)>0) && separationFlag 
//             && ((!certified)||separationFlag_certified)
            && widthFlag && compactFlag ) ) ) {
            
            cauchy_conjugate ( &ccurConj, &pushConjugFlag, &separationFlag, ccur, meta );

        }
#ifdef CCLUSTER_TIMINGS    
        time_in_real_coeffs += (double) (clock() - startt3);
#endif 
        
        if ( (connCmp_nSols(ccur)>0) && separationFlag 
//             && ((!certified)||separationFlag_certified)
            && widthFlag && compactFlag && (connCmp_nSols(ccur)<cacheCauchy_degreeref(cacheCau))) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );

#ifdef CCLUSTER_TIMINGS
        startt3 = clock();
#endif            
            if (certified) { 
                slong m = connCmp_nSols(ccur);
                connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
                compBox_get_containing_dsk(ccDisk, componentBox);
                        
                compDsk_set(twoCCDisk, ccDisk);
                realRat_mul_si(compDsk_radiusref(twoCCDisk), compDsk_radiusref(ccDisk), 2);
                
                cauchyTest_res resCert; 
                resCert.nbOfSol = 1;
                if (m>1) {
                    slong q = cauchyTest_getNbEvals_counting_combinatorial_with_isoRatio(compDsk_centerref(ccDisk), compDsk_radiusref(twoCCDisk), two, m, cacheCau);
                    if (metadatas_getVerbo(meta)>=2) {
                        printf("#precision before Pellet counting: %ld\n", prec);
                        printf("#number of eval points with Pellet: %ld\n", cacheCauchy_degreeref(cacheCau)+1);
                        printf("#number of eval points with Combin: %ld\n", q);
                    }
                    resCert = cauchyTest_Pellet_counting( compDsk_centerref(ccDisk), compDsk_radiusref(twoCCDisk), two, m, cacheCau, prec, 
                                                            meta, depth);
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#result: %d, precision: %ld\n", resCert.nbOfSol, resCert.appPrec);
                    }
                } 
                else if (metadatas_useCompression(meta) ){ /* m=1 */
                    clock_t start2 = clock();
//                     realRat_mul_si(compDsk_radiusref(twoCCDisk), compDsk_radiusref(ccDisk), 1);
                    /* check that twoCCDisk contains at least a root with proba root counter */
                    cauchyTest_res resCauchy = cauchyTest_probabilistic_counting( ccDisk, cacheCau, prec, meta, depth);
                    if ( resCauchy.nbOfSol == 0 ) {
                        resCert.nbOfSol = -1;
                        failure = 1;
//                         printf("m: %ld, resCauchy.nbOfSol: %d\n", m, resCauchy.nbOfSol);
                        printf("#FAILURE: CERTIFICATION FAILED: A DISK is not certified to contain at least one root \n");
                        continue;
                    }
                    metadatas_add_time_CompTot(meta, (double) (clock() - start2));
                }
                
                if (resCert.nbOfSol <= -1) {
                                failure = 1;
                                printf("#FAILURE: CERTIFICATION FAILED: A DISK is NOT 2-ISOLATED \n");
                                continue;
                            }
                            
                /* else inflate the cc */
        //         realRat_set_si(mp2, m+2, 1);
                metadatas_add_nbClusterCertified( meta, 1);
                connCmp_infate_realRat_inPlace(ccur,  two, metadatas_initBref(meta) );
            }
#ifdef CCLUSTER_TIMINGS    
        time_in_certified += (double) (clock() - startt3);
#endif
            
            connCmp_list_push(qResults, ccur);
            nbSolsInQResults += connCmp_nSols(ccur);
            
            if (metadatas_getVerbo(meta)>=level)
                printf("------validated with %d roots\n", connCmp_nSols(ccur));
            
//             printf("metadatas_useRealCoeffs(meta): %d, pushConjugFlag: %d\n", metadatas_useRealCoeffs(meta), pushConjugFlag);
            /* Real Coeff */
            if ((metadatas_useRealCoeffs(meta))&&(pushConjugFlag)){
                /*compute the complex conjugate*/
                metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
                connCmp_list_push(qResults, ccurConj);
                nbSolsInQResults += connCmp_nSols(ccur);
                if (certified) {
                    metadatas_add_nbClusterCertified( meta, 1);
                }
            }
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
            
            if (metadatas_getVerbo(meta)>=level)
                printf("------push in working queue after newton success %d roots\n", connCmp_nSols(ccur));
#ifdef CCLUSTER_TIMINGS
            startt3 = clock();
#endif            
//             connCmp_list_insert_sorted(qMainLoop, ccur);
            it = connCmp_list_insert_sorted_from_end(qMainLoop, ccur, it);
#ifdef CCLUSTER_TIMINGS    
            time_in_insert_sorted += (double) (clock() - startt3);
#endif
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            connCmp_decrease_nwSpd(ccur);

            if (metadatas_getVerbo(meta)>=level)
                printf("------push in working queue after newton fail %d roots\n", connCmp_nSols(ccur));
#ifdef CCLUSTER_TIMINGS
            startt3 = clock();
#endif            
//             connCmp_list_insert_sorted(qMainLoop, ccur);
            it = connCmp_list_insert_sorted_from_end(qMainLoop, ccur, it);
#ifdef CCLUSTER_TIMINGS    
            time_in_insert_sorted += (double) (clock() - startt3);
#endif
        }
//         else if (do_not_bisect_flag ){
//             it = connCmp_list_insert_sorted_from_end(qMainLoop, ccur, it);
//         }
        else {

            if (metadatas_getVerbo(meta)>=level)
                printf("------bisect and push in working queue\n");
            
            cauchy_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cacheCau, nbMaxSol, meta,1);
            while (!connCmp_list_is_empty(ltemp)) {
#ifdef CCLUSTER_TIMINGS
                startt3 = clock();
#endif
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
                it = connCmp_list_insert_sorted_from_end(qMainLoop, connCmp_list_pop(ltemp), it);
#ifdef CCLUSTER_TIMINGS    
                time_in_insert_sorted += (double) (clock() - startt3);
#endif
            }
            connCmp_clear(ccur);
            ccluster_free(ccur);
        }
    }
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    compDsk_clear(twoCCDisk);
    realRat_clear(two);
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(mp2);
    realRat_clear(neps);
    realRat_clear(threeWidth);
    connCmp_list_clear(ltemp);
    
//     printf("time in UnionFind: %f\n", timeInUnionFind/CLOCKS_PER_SEC);

#ifdef CCLUSTER_TIMINGS    
    time_in_cauchy_main_loop += (double) (clock() - startt);
#endif
    
    return failure;
}

int cauchy_algo_global( connCmp_list_t qResults, 
                           compBox_list_t bDiscarded,
                           const realRat_t eps, 
                           cacheCauchy_t cacheCau,
                           int certified,
                           metadatas_t meta){
    
#ifdef CCLUSTER_TIMINGS
    time_in_cauchy_discard_compBox_list   = 0.0;
    time_in_cauchy_bisect_connCmp         = 0.0;
    time_in_cauchy_MaxNumberOfRootsInCC   = 0.0;
    time_in_cauchy_compression            = 0.0;
    time_in_cauchy_newton                 = 0.0;
    time_in_cauchy_connCmp_is_separated           = 0.0;
    time_in_cauchy_connCmp_is_separated_certified = 0.0;
    time_in_cauchy_main_loop                      = 0.0;
    time_in_proba_counting                        = 0.0;
    time_in_real_coeffs                           = 0.0;
    time_in_certified                             = 0.0;
    time_in_insert_sorted                         = 0.0;
    
    time_in_cauchyTest_getEvaluationPoints =0.0;
    time_in_cauchyTest_rootsOfUnits        =0.0; 
    time_in_cauchyTest_shift_points        =0.0; 
    time_in_cauchyTest_computeSsApprox        =0.0; 
    time_in_cauchyTest_computeSsApprox_fromVals        =0.0; 
    time_in_cacheCauchy_set_bounds        =0.0; 
    time_in_cacheCauchy_eval              =0.0;
    time_in_cacheCauchy_eval_powering     =0.0;
    nbPos = 0;
#endif
    
    clock_t start = clock();
    int level = 3;
    
//     realRat_t factor;
//     realRat_init(factor);
//     realRat_set_si(factor, 5, 4);
    
    compBox_ptr initialBox;
    initialBox = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    
    /* automaticly set initialBox */
    /* with evaluation function */
    cauchyRootRadii_root_bound( compBox_bwidthref(initialBox), cacheCau, meta );
    if (metadatas_getVerbo(meta)>=level) {
        printf("#root bound with eval function: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
//     compBox_set_si(initialBox, 0,1,0,1,1024,1);
    
    metadatas_setInitBox(meta, initialBox);
    compBox_nbMSolref(initialBox) = cacheCauchy_degreeref(cacheCau);
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, initialBox);
    
    connCmp_list_t qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
//     connCmp_list_init(qPrepLoop);
    connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qMainLoop, initialCC);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster preploop: \n");
//     ccluster_prep_loop( qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster mainloop: \n");
    int failure = cauchy_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cacheCau, certified, meta);
    
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
#ifdef CCLUSTER_TIMINGS
//     printf("time_in_cauchy_discard_compBox_list          : %f\n", time_in_cauchy_discard_compBox_list           /CLOCKS_PER_SEC);
    printf("time_in_cauchy_bisect_connCmp                : %f\n", time_in_cauchy_bisect_connCmp                 /CLOCKS_PER_SEC);
    printf("time_in_cauchy_MaxNumberOfRootsInCC          : %f\n", time_in_cauchy_MaxNumberOfRootsInCC           /CLOCKS_PER_SEC);
    printf("time_in_proba_counting                       : %f\n", time_in_proba_counting                        /CLOCKS_PER_SEC);
    printf("time_in_cauchy_compression                   : %f\n", time_in_cauchy_compression                    /CLOCKS_PER_SEC);
    printf("time_in_cauchy_newton                        : %f\n", time_in_cauchy_newton                         /CLOCKS_PER_SEC);
    printf("time_in_cauchy_connCmp_is_separated          : %f\n", time_in_cauchy_connCmp_is_separated           /CLOCKS_PER_SEC);
    printf("time_in_cauchy_connCmp_is_separated_certified: %f\n", time_in_cauchy_connCmp_is_separated_certified /CLOCKS_PER_SEC);
    printf("time_in_real_coeffs                          : %f\n", time_in_real_coeffs                           /CLOCKS_PER_SEC);
    printf("time_in_certified                            : %f\n", time_in_certified                             /CLOCKS_PER_SEC);
    printf("time_in_insert_sorted                        : %f\n", time_in_insert_sorted                         /CLOCKS_PER_SEC);
    printf("                                             ---------\n");
    double total = time_in_cauchy_bisect_connCmp + time_in_cauchy_MaxNumberOfRootsInCC + time_in_proba_counting
                 + time_in_cauchy_compression    + time_in_cauchy_newton + time_in_cauchy_connCmp_is_separated
                 + time_in_cauchy_connCmp_is_separated_certified + time_in_real_coeffs + time_in_certified
                 + time_in_insert_sorted; 
    printf("                                               %f\n", total                                         /CLOCKS_PER_SEC);
    printf("time_in_cauchy_main_loop                     : %f\n", time_in_cauchy_main_loop                      /CLOCKS_PER_SEC);
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
    printf("time_in_cauchyTest_getEvaluationPoints       : %f\n", time_in_cauchyTest_getEvaluationPoints        /CLOCKS_PER_SEC);
    printf("time_in_cauchyTest_rootsOfUnits              : %f\n", time_in_cauchyTest_rootsOfUnits               /CLOCKS_PER_SEC);
    printf("time_in_cauchyTest_shift_points              : %f\n", time_in_cauchyTest_shift_points               /CLOCKS_PER_SEC);
    printf("time_in_cauchyTest_computeSsApprox           : %f\n", time_in_cauchyTest_computeSsApprox            /CLOCKS_PER_SEC);
    printf("time_in_cauchyTest_computeSsApprox_fromVals  : %f\n", time_in_cauchyTest_computeSsApprox_fromVals   /CLOCKS_PER_SEC);
    printf("time_in_cacheCauchy_set_bounds               : %f\n", time_in_cacheCauchy_set_bounds                /CLOCKS_PER_SEC);
    printf("time_in_cacheCauchy_eval                     : %f\n", time_in_cacheCauchy_eval                 /CLOCKS_PER_SEC);
    printf("time_in_cacheCauchy_eval_powering            : %f\n", time_in_cacheCauchy_eval_powering                /CLOCKS_PER_SEC);
    printf("number of exclusion test returning >0        :%d \n", nbPos);
#endif
    
    return failure;
}

int metadatas_cauchy_fprint(FILE * file, metadatas_t meta, const realRat_t eps, cacheCauchy_t cacheCau, int certified){
    int r=1;
    char temp[1000];
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, "# -------------------Cauchy Root Finder: ------------------------------\n");
    r = fprintf(file, "# -------------------Input:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    slong degree = cacheCauchy_degreeref(cacheCau);
    cacheApp_ptr cache = cacheCauchy_cacheAppref(cacheCau);
    if (cacheCauchy_evalFastref(cacheCau) == NULL) {
        if (cache->_from_poly) {
            slong bitsize = cacheApp_getBitsize (cache);
            r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10ld|\n", "exact coefficients", degree, bitsize);
        } else
            r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10c|\n", "approximated coefficients", degree, '?');
    }
    else
        r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10c|\n", "evaluation procedure", degree, '?');
    }
    compBox_sprint_for_stat( temp, metadatas_initBref(meta) );
    r = fprintf(file, "#|box:%-65s\n", temp);
    if (realRat_is_den_zero( eps ))
        r = fprintf(file, "#|eps: %-64s|\n", "+inf");
    else {
        realRat_sprint_for_stat( temp, eps );
        r = fprintf(file, "#|eps: %-64s|\n", temp);
    }
    int len = 0;
    //TODO find a better way for this...
    if ( metadatas_useNewton(meta) &&
         metadatas_useRealCoeffs(meta) &&
         metadatas_useCompression(meta)
       ) len += sprintf( temp + len, " default");
    else {    
        if (metadatas_useNewton(meta)) len += sprintf( temp + len, " newton");
        if (metadatas_useRealCoeffs(meta)) len += sprintf( temp + len, " realCoeffs");
        if (metadatas_useCompression(meta)) len += sprintf( temp + len, " compression");
    }
    if (metadatas_usePowerSums(meta)) len += sprintf( temp + len, " + powerSums");
    if (metadatas_forTests(meta)) len += sprintf( temp + len, " + test");
    if (metadatas_usefpri(meta)) len += sprintf( temp + len, " + fpri");
#ifdef CCLUSTER_HAVE_PTHREAD
    if (metadatas_useNBThreads(meta)>1) len += sprintf( temp + len, " %d threads", metadatas_useNBThreads(meta));
#endif
    if (metadatas_stratref(meta)->_additionalFlags !=0) 
        len += sprintf(temp +len, " %d", metadatas_stratref(meta)->_additionalFlags);
    r = fprintf(file, "#|strat:%-63s|\n", temp);
    
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, "# -------------------Cauchy exclusion tests----------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number ET:",                       metadatas_getNbCauchyExTests(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in tests  ET:",          metadatas_get_time_CauExTo(meta),    " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals            ET:",         metadatas_getNbCauchyExEvals(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals            ET:",         metadatas_get_time_CauExEv(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing divs     ET:",         metadatas_get_time_CauExDS(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing s0s      ET:",         metadatas_get_time_CauExCS(meta),    " " );
    r = fprintf(file, "# -------------------Cauchy counting tests----------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number CT:",                       metadatas_getNbCauchyCoTests(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in tests  CT:",          metadatas_get_time_CauCoTo(meta),    " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals            CT:",         metadatas_getNbCauchyCoEvals(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals            CT:",         metadatas_get_time_CauCoEv(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing divs     CT:",         metadatas_get_time_CauCoDS(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing s0s      ET:",         metadatas_get_time_CauCoCS(meta),    " " );
//     r = fprintf(file, "# ---------------------------------------------------------------------\n");
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time in shift for FFT        :",         metadatas_get_time_Taylors(meta),    " " );
    if (metadatas_useCompression(meta)) {
    r = fprintf(file, "# -------------------Compression into rigid discs---------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number for clus of  1 root:",       metadatas_getComp_nb_1(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number for clus of >1 root:",       metadatas_getComp_nb_p(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in comp s1/m:",                metadatas_get_time_CompCen(meta),    " " );
#ifdef DEFLATION_TURAN
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR with Turan's theorem:", metadatas_get_time_CompRR2(meta),    " " );
#else
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR with doub. exp. sieve:", metadatas_get_time_CompRR1(meta),    " " );
#endif
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR algo 3:",                metadatas_get_time_CompRR3(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in compression:",        metadatas_get_time_CompTot(meta),    " " );
    }
    if (metadatas_useNewton(meta)){
    r = fprintf(file, "# -------------------Newton Iterations---------------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number NE:",                     metadatas_getNbNewton(meta),         " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of fails:",                    metadatas_getNbFailingNewton(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in newton:",         metadatas_get_time_Newtons(meta),    " " );
    }
    if (certified) {
    r = fprintf(file, "# -------------------Final Certification-------------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of certified clusters:",       metadatas_get_nbClusterCertified(meta),         " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of Pellet's counting tests:",  metadatas_get_nbCertifiedWithPellet(meta),      " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time in Pellet:",               metadatas_get_time_CertPel(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "      time in Evaluations:",          metadatas_get_time_CertPEv(meta),    " " );
    r = fprintf(file, "#|%-39s %14ld %14s|\n", "        nb of Evaluations:",          metadatas_get_nbEvalsInPellet(meta),      " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "      time in Interpolations:",       metadatas_get_time_CertPIn(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "      time in DLG iterations:",       metadatas_get_time_CertPGr(meta),    " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "        nb of DLG iterations:",       metadatas_get_nbGraeffeInPellet(meta),      " " );
    }
    r = fprintf(file, "# -------------------Other---------------------------------------------\n");
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in getApproximation:",           metadatas_get_time_Approxi(meta),    " " );
    double timeInEval = metadatas_get_time_CauExEv(meta) + metadatas_get_time_CauCoEv(meta) + metadatas_get_time_CertPEv(meta);
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in Evaluation:",                 timeInEval,    " " );
    int    nbOfEvals  = metadatas_getNbCauchyExEvals(meta) + metadatas_getNbCauchyCoEvals(meta) + metadatas_get_nbEvalsInPellet(meta);
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number of evaluations:",        nbOfEvals,    " " );
    r = fprintf(file, "# -------------------Precision-----------------------------------------\n");
    r = metadatas_boxes_by_prec_fprint ( file, meta );
     
    }
   
    r = fprintf(file, "# -------------------Output:   ----------------------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of clusters:",                 metadatas_getNbValidated(meta),      " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of solutions:",                metadatas_getNbSolutions(meta),      " " );
    r = fprintf(file, "# -------------------Stats:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, "#|%-39s %14d %14s|\n", "tree depth:",                         metadatas_getDepth(meta),            " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "tree size:",                          metadatas_getNbExplored(meta),       " " );
    }
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time:",                         metadatas_get_time_CclusAl(meta),    " " );
    r = fprintf(file, "# ---------------------------------------------------------------------\n");
    }
    return r;
}
