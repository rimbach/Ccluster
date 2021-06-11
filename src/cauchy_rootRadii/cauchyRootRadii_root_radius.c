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

#include "cauchy_rootRadii/cauchy_rootRadii.h"

/* returns 1 if ( radSup <= eps ) or ( radSup > radInf*relativeError ) */
/*         0 otherwise                                                 */
int _cauchyRootRadii_stoppingCriterion ( const realRat_t radInf, const realRat_t radSup, const realRat_t relativeError, const realRat_t eps ) {
    
    int stop = ( realRat_cmp( radSup, eps ) <= 0 );
    
    if (stop==0) {
        realRat_t temp;
        realRat_init(temp);
        realRat_mul(temp, relativeError, radInf);
        stop = ( realRat_cmp( radSup, temp ) <= 0 );
        realRat_clear(temp);
    }
    
    return stop;
}

/* compute a rational t so that */
/* sqrt(radInf*radSup) <= t < (relativeError)^(1/4)*sqrt(radInf*radSup) */
/* returns required intermediate precision */
slong _cauchyRootRadii_findRational ( realRat_t t, const realRat_t radInf, const realRat_t radSup, const realRat_t relativeError, slong prec ){
    
    
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational begin\n");
//     }
        
    realApp_t tApp;
    realApp_init(tApp);
    
    fmpz_t a,b,exp;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(exp);
    
    realRat_t temp1, temp2;
    realRat_init(temp1);
    realRat_init(temp2);
    
    realRat_mul( t, radInf, radSup );
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
    realApp_set_realRat( tApp, t, prec );
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    realApp_sqrt(tApp, tApp, prec);
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    arb_get_interval_fmpz_2exp(a, b, exp, tApp);
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational a: "); fmpz_print(a); printf("\n");
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational b: "); fmpz_print(b); printf("\n");
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational exp: "); fmpz_print(exp); printf("\n");
    realRat_set_si( t, 2, 1);
    realRat_pow_fmpz( t, t, exp);
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
    realRat_mul_fmpz( t, t, b );
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
    realApp_set_realRat( tApp, t, prec );
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    
    /* check t^4 <= relativeError*(radInf*radSup)^2 */
    realRat_pow_si(temp1, t, 4);
    realRat_mul(temp2, radInf, radSup);
    realRat_pow_si(temp2, temp2, 2);
    realRat_mul(temp2, temp2, relativeError);
    
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational temp1: "); realRat_print(temp1); printf("\n");
//     printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational temp2: "); realRat_print(temp2); printf("\n");
    
    while (realRat_cmp(temp1, temp2)>0) {
        
        prec = 2*prec;
        
        realRat_mul( t, radInf, radSup );
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
        realApp_set_realRat( tApp, t, prec );
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
        realApp_sqrt(tApp, tApp, prec);
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
        arb_get_interval_fmpz_2exp(a, b, exp, tApp);
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational a: "); fmpz_print(a); printf("\n");
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational b: "); fmpz_print(b); printf("\n");
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational exp: "); fmpz_print(exp); printf("\n");
        realRat_set_si( t, 2, 1);
        realRat_pow_fmpz( t, t, exp);
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
        realRat_mul_fmpz( t, t, b );
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational t: "); realRat_print(t); printf("\n");
        
        /* check t^4 <= relativeError*(radInf*radSup)^2 */
        realRat_pow_si(temp1, t, 4);
        realRat_mul(temp2, radInf, radSup);
        realRat_pow_si(temp2, temp2, 2);
        realRat_mul(temp2, temp2, relativeError);
        
    }
        
    realApp_clear(tApp);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(exp);
    realRat_clear(temp1);
    realRat_clear(temp2);
    
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("# --------- cauchyRootRadii_root_radius.c: _cauchyRootRadii_findRational end\n");
//     }
    
    return prec;
}

/* Assume radInf < r_{d+1-m}(center, p) < radSup */
/* Modifies radInf, radSup so that: */
/*   either radInf < r_{d+1-m}(center, p) < radSup <= relativeError*radInf */
/*       or radInf < r_{d+1-m}(center, p) < radSup <= eps */
/* assume relativeError >= 19/10 */
/* so that when 1 < a <= 12/11, fmatheta(a,4/3) > relativeError^(1/4) */
void cauchyRootRadii_root_radius( const compRat_t center,
                                  realRat_t radInf,        /* radInf < r_{d+1-m}(center, p) */
                                  realRat_t radSup,        /* radSup > r_{d+1-m}(center, p) */
                                  const realRat_t relativeError, /* want relativeError*radInf >= radSup */ 
                                  const realRat_t eps,
                                  const realRat_t theta,   /*isolation ratio of the disk in which is computed rr */ 
                                  slong nbOfRoots,
                                  cacheCauchy_t cacheCau,
                                  cacheApp_t cache,
                                  metadatas_t meta ){
    
    cauchyTest_res cres;
    cres.appPrec = CCLUSTER_DEFAULT_PREC;
    slong precForT = CCLUSTER_DEFAULT_PREC;
    
    compDsk_t Delta;
    compDsk_init(Delta);
    compRat_set( compDsk_centerref(Delta), center );
    
    realRat_t epsprime, a, fma;
    realRat_init(a);
    realRat_init(fma);
    realRat_init(epsprime);
    
    realRat_t t;
    realRat_init(t);
    
    realRat_div(epsprime, eps, theta);
    
    realRat_set_si(a, 12, 11);
    cauchyTest_fmatheta(fma, a, cacheCauchy_isoRatioref( cacheCau ) );
    
    int stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps );
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("# ------ cauchyRootRadii_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
        printf("# ------ cauchyRootRadii_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
        printf("# ------ cauchyRootRadii_root_radius.c: stop is %d \n", stop);
    }
        
    if ( (stop==0) && ( realRat_is_zero( radInf ) ) ) {
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("# ------ cauchyRootRadii_root_radius.c: radInf is 0 \n");
        }
    
        /* Apply deterministic root counter to D(center, epsprime)*/
        realRat_set( compDsk_radiusref(Delta), epsprime );
        cres = cauchyTest_deterministic_counting( Delta, a, cache, cacheCau, cres.appPrec, meta, 0);
        if (cres.nbOfSol == -1) { /* A(center, fma*epsprime, fpa*epsprime) contains a root */
                                  /* r_{d+1-m}(center, p) > fma*epsprime */
            realRat_mul( radInf, fma, epsprime );
        } else if ( cres.nbOfSol < nbOfRoots ) {
                                  /* r_{d+1-m}(center, p) > epsprime */
            realRat_set( radInf, epsprime );
        } else {                  /* cres.nbOfSol = nbOfRoots */
                                  /* r_{d+1-m}(center, p) < epsprime < eps */
            realRat_set( radSup, eps );
        }
        
        
    }
    
    stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps );
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("# ------ cauchyRootRadii_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
        printf("# ------ cauchyRootRadii_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
        printf("# ------ cauchyRootRadii_root_radius.c: stop is %d \n", stop);
    }
    
    while( stop == 0 ) {
    
        /* compute a rational t so that */
        /* sqrt(radInf*radSup) <= t < (relativeError)^(1/4)*sqrt(radInf*radSup) */
        precForT = _cauchyRootRadii_findRational ( t, radInf, radSup, relativeError, precForT );
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("# ------ cauchyRootRadii_root_radius.c: t is "); realRat_print(t); printf("\n");
        }
        
        /* Apply deterministic root counter to D(center, t)*/
        cres = cauchyTest_deterministic_counting( Delta, a, cache, cacheCau, cres.appPrec, meta, 0);
        if (cres.nbOfSol == -1) { /* A(center, fma*t, fpa*t) contains a root */
                                  /* r_{d+1-m}(center, p) > fma*t */
            realRat_mul( radInf, fma, t );
        } else if ( cres.nbOfSol < nbOfRoots ) {
                                  /* r_{d+1-m}(center, p) > t */
            realRat_set( radInf, t );
        } else {                  /* cres.nbOfSol = nbOfRoots */
                                  /* r_{d+1-m}(center, p) < t */
            realRat_set( radSup, t );
        }
        
        stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps );
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("# ------ cauchyRootRadii_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
            printf("# ------ cauchyRootRadii_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
            printf("# ------ cauchyRootRadii_root_radius.c: stop is %d \n", stop);
        }
    }
    
    compDsk_clear(Delta);
    realRat_clear(a);
    realRat_clear(fma);
    realRat_clear(epsprime);
    realRat_clear(t);
}

// void cauchyRootRadii_root_radius( const compRat_t center,
//                                   realRat_t radInf, /* radInf < r_{d+1-m}(center, p) */
//                                   realRat_t radSup, /* radSup > r_{d+1-m}(center, p) */
//                                   realRat_t relativeError, /* want relativeError*radInf >= radSup */ 
//                                   slong nbOfRoots,
//                                   cacheCauchy_t cacheCau,
//                                   cacheApp_t cache,
//                                   metadatas_t meta ){
//     
//     /* 2^logRelError <= relativeError */
//     ulong logRelError = fmpz_flog_ui( realRat_numref( relativeError ), 2 ) - fmpz_clog_ui( realRat_denref( relativeError ), 2 );
//     slong precForMiddle = 2*logRelError;
//     
//     cauchyTest_res cres;
//     cres.appPrec = CCLUSTER_DEFAULT_PREC;
//     
//     int stop = 0;
//     realRat_t ratio, middle;
//     realRat_init(ratio);
//     realRat_init(middle);
//     
//     realApp_t mid;
//     realApp_init(mid);
//     
//     fmpz_t a,b,exp;
//     fmpz_init(a);
//     fmpz_init(b);
//     fmpz_init(exp);
//     
//     compDsk_t Delta;
//     compDsk_init(Delta);
//     compRat_set( compDsk_centerref(Delta), center );
//     
//     realRat_div(ratio, radSup, radInf);
//     stop = ( realRat_cmp( ratio, relativeError ) <= 0);
//     
//     
//     while (stop==0) {
//         
//         realRat_mul( middle, radInf, radSup );
//         realApp_set_realRat( mid, middle, precForMiddle );
//         realApp_sqrt(mid, mid, precForMiddle);
//         arb_get_interval_fmpz_2exp(a, b, exp, mid);
//         realRat_set_fmpz(middle, a);
//         fmpz_set_si(b, 2);
//         fmpz_pow_fmpz(b, exp);
//         fmpq_mul_fmpz(middle, middle, b);
//         
//         /* call probabilistic root counter */
//         realRat_set( compDsk_radiusref(Delta), middle );
//         cres = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//         if (cres.nbOfSol == -1) { /* there are roots in A(c, middle/(theta), middle*theta ) with theta = 4/3 */
//                                   /* r_{d+1-m}(center, p) > middle/(theta) */
//             realRat_mul( radInf, middle, cacheCauchy_isoRatioref( cacheCau ) );
//         } else if (cres.nbOfSol != nbOfRoots) { /* either the result is correct and r_{d+1-m}(center, p) > middle */
//                                                /* or     the result is not correct and there are roots in A(c, middle/(theta), middle*theta ) */
//                                                /* r_{d+1-m}(center, p) > middle/(theta) */
//             realRat_mul( radInf, middle, cacheCauchy_isoRatioref( cacheCau ) );
//         } else { /* need to check the result with deterministic counting */
//             cres = cauchyTest_deterministic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//             if (cres.nbOfSol == -1) { /* there are roots in A(c, 13/24*middle, 53/24*middle); */
//                                       /* r_{d+1-m}(center, p) > 13/24*middle */
//                 realRat_mul_si(radInf, middle, 13);
//                 realRat_div_ui(radInf, radInf, 24);
//         } else if (cres.nbOfSol != nbOfRoots) { /* necessarily cres.nbOfSol < nbOfRoots because the result is correct */
//                                                 /* r_{d+1-m}(center, p) > middle */
//         
//                 realRat_set( radInf, middle );
//         } else { /* cres.nbOfSol == nbOfRoots */
//                  /* r_{d+1-m}(center, p) < middle */
//                 realRat_set( radSup, middle ); 
//         }
//         /* check stop criterion */
//         realRat_div(ratio, radSup, radInf);
//         stop = ( realRat_cmp( ratio, relativeError ) <= 0);
//     }
//     
//     realRat_clear(ratio);
//     realRat_clear(middle);
//     realApp_clear(mid);
//     fmpz_clear(a);
//     fmpz_clear(b);
//     fmpz_clear(exp);
// }
