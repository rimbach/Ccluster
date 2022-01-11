/* ************************************************************************** */
/*  Copyright (C) 2022 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "pellet_rootRadii/pellet_rootRadii.h"

/* returns 1 if ( radSup <= eps ) or ( radSup > radInf*relativeError ) */
/*         0 otherwise                                                 */
int _pellet_RootRadii_stoppingCriterion ( const realRat_t radInf, const realRat_t radSup, const realRat_t relativeError, const realRat_t eps, const metadatas_t meta ) {
    
    int level = 5;
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------pellet_root_radius.c: radInf        is "); realRat_print(radInf); printf("\n");
        printf("#                                  radSup        is "); realRat_print(radSup); printf("\n");
        printf("#                                  eps           is "); realRat_print(eps   ); printf("\n");
        printf("#                                  relativeError is "); realRat_print(relativeError); printf("\n");
    }
    
    int stop = ( realRat_cmp( radSup, eps ) <= 0 );
    
    if (stop==0) {
        realRat_t temp;
        realRat_init(temp);
        realRat_mul(temp, relativeError, radInf);
        stop = ( realRat_cmp( radSup, temp ) <= 0 );
        realRat_clear(temp);
    }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#                                  stop          is %d \n", stop);
    }
    
    return stop;
}

/* compute a rational t so that */
/* sqrt(radInf*radSup) <= t < (relativeError)^(1/4)*sqrt(radInf*radSup) */
/* returns required intermediate precision */
slong _pellet_RootRadii_findRational ( realRat_t t, const realRat_t radInf, const realRat_t radSup, const realRat_t relativeError, slong prec ){
    
    
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("# ---------pellet_root_radius.c: _pellet_RootRadii_findRational begin\n");
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
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
    realApp_set_realRat( tApp, t, prec );
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    realApp_sqrt(tApp, tApp, prec);
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    arb_get_interval_fmpz_2exp(a, b, exp, tApp);
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational a: "); fmpz_print(a); printf("\n");
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational b: "); fmpz_print(b); printf("\n");
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational exp: "); fmpz_print(exp); printf("\n");
    realRat_set_si( t, 2, 1);
    realRat_pow_fmpz( t, t, exp);
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
    realRat_mul_fmpz( t, t, b );
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
    realApp_set_realRat( tApp, t, prec );
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
    
    /* check t^4 <= relativeError*(radInf*radSup)^2 */
    realRat_pow_si(temp1, t, 4);
    realRat_mul(temp2, radInf, radSup);
    realRat_pow_si(temp2, temp2, 2);
    realRat_mul(temp2, temp2, relativeError);
    
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational temp1: "); realRat_print(temp1); printf("\n");
//     printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational temp2: "); realRat_print(temp2); printf("\n");
    
    while (realRat_cmp(temp1, temp2)>0) {
        
        prec = 2*prec;
        
        realRat_mul( t, radInf, radSup );
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
        realApp_set_realRat( tApp, t, prec );
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
        realApp_sqrt(tApp, tApp, prec);
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational tApp: "); realApp_printd(tApp, 53); printf("\n");
        arb_get_interval_fmpz_2exp(a, b, exp, tApp);
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational a: "); fmpz_print(a); printf("\n");
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational b: "); fmpz_print(b); printf("\n");
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational exp: "); fmpz_print(exp); printf("\n");
        realRat_set_si( t, 2, 1);
        realRat_pow_fmpz( t, t, exp);
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
        realRat_mul_fmpz( t, t, b );
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational t: "); realRat_print(t); printf("\n");
        
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
//         printf("# --------- pellet_root_radius.c: _pellet_RootRadii_findRational end\n");
//     }
    
    return prec;
}

/* Assume radInf < r_{d+1-m}(center, p) < radSup */
/* Modifies radInf, radSup so that: */
/*   either radInf < r_{d+1-m}(center, p) < radSup <= relativeError*radInf */
/*       or radInf < r_{d+1-m}(center, p) < radSup <= eps */
/* assume relativeError >= 19/10 */
/* the Pellet test succeeds if there are no root in A(c, 0.94, 4/3) */
/* thus 0.94 > relativeError^(-1/4) ... also works when relativeError >=19/10 */
void pellet_root_radius( const compRat_t center,
                         realRat_t radInf,        /* radInf < r_{d+1-m}(center, p) */
                         realRat_t radSup,        /* radSup > r_{d+1-m}(center, p) */
                         const realRat_t relativeError, /* want relativeError*radInf >= radSup */ 
                         const realRat_t eps,
                         const realRat_t theta,   /*isolation ratio of the disk in which is computed rr */ 
                         slong nbOfRoots,
                         cacheApp_t cache,
                         metadatas_t meta,
                         slong prec){
    
    int level = 4;
    slong precForT = CCLUSTER_DEFAULT_PREC;
    
    compDsk_t Delta;
    compDsk_init(Delta);
    compRat_set( compDsk_centerref(Delta), center );
    
    realRat_t epsprime, rho1;
    realRat_init(rho1);
    realRat_init(epsprime);
    
    realRat_t t;
    realRat_init(t);
    
    realRat_div(epsprime, eps, theta);
    
    realRat_set_si(rho1, 94, 100);
    
    int stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps, meta );
    
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("#------------pellet_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
//         printf("#------------pellet_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
//         printf("#------------pellet_root_radius.c: stop is %d \n", stop);
//     }
    
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,deg+1);

    
    if ( (stop==0) && ( realRat_is_zero( radInf ) ) ) {
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------------pellet_root_radius.c: radInf is 0 \n");
        }
 
        /* Apply deterministic root counter to D(center, epsprime)*/
        realRat_set( compDsk_radiusref(Delta), epsprime );
        cres = cauchyTest_rootFreeAnnulus_verification( Delta, nbOfRoots, a,
                                                      cacheCau, cres.appPrec, meta, 0);
        if (metadatas_getVerbo(meta)>=level)
                printf("#---------------pellet_root_radius.c: return of deterministic verification: cres.nbOfSol: %d \n", cres.nbOfSol);
        
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
        
        stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps, meta );
    
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#---------------pellet_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
//             printf("#---------------pellet_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
//             printf("#---------------pellet_root_radius.c: stop is %d \n", stop);
//         }
    
    }
    
    while( stop == 0 ) {
    
        /* compute a rational t so that */
        /* sqrt(radInf*radSup) <= t < (relativeError)^(1/4)*sqrt(radInf*radSup) */
        precForT = _cauchyRootRadii_findRational ( t, radInf, radSup, relativeError, precForT );
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------pellet_root_radius.c: t is "); realRat_print(t); printf("\n");
        }
        realRat_set( compDsk_radiusref(Delta), t );
        
        /* Apply deterministic verification to D(center, t)*/
        cres = cauchyTest_rootFreeAnnulus_verification( Delta, nbOfRoots, a,
                                                      cacheCau, cres.appPrec, meta, 0);
        if (metadatas_getVerbo(meta)>=level)
                printf("#------------pellet_root_radius.c: return of deterministic verification: cres.nbOfSol: %d \n", cres.nbOfSol);
        
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
        
        stop = _cauchyRootRadii_stoppingCriterion ( radInf, radSup, relativeError, eps, meta );
        
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#------------pellet_root_radius.c: radInf is "); realRat_print(radInf); printf("\n");
//             printf("#------------pellet_root_radius.c: radSup is "); realRat_print(radSup); printf("\n");
//             printf("#------------pellet_root_radius.c: stop is %d \n", stop);
//         }
    }
    
    compApp_poly_clear(pApprox);
    
    compDsk_clear(Delta);
    realRat_clear(a);
    realRat_clear(fma);
    realRat_clear(epsprime);
    realRat_clear(t);
}
