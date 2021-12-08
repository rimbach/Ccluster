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

#include "turan_rootRadii/turan_rootRadii.h"

slong turanRootRadii_root_bound( realApp_t rm,
                                 const compDsk_t Delta,
                                 slong           m,
                                 const realRat_t theta,
                                 ulong           N,
                                 const realRat_t eps,
                                 cacheCauchy_t cacheCau,
                                 slong prec,
                                 metadatas_t meta, int depth) {
    
    int level = 3;
    
    realRat_t epsprime, eps2;
    realRat_init(epsprime);
    realRat_init(eps2);
    
    /* set relative error */
    realApp_t delta, deltam1, temp, epsApp;
    realApp_init(delta);
    realApp_init(deltam1);
    realApp_init(temp);
    realApp_init(epsApp);
    
    realApp_set_si(delta, 5);
    realApp_root_ui(delta, delta, N, prec);
    realApp_sub_si (deltam1, delta, 1, prec);
    /* automatically set eps2 */
    realRat_set_si(eps2, 1, 2);
    realApp_set_realRat(temp, eps2, prec);
    while ( realApp_gt( temp, deltam1 ) ) {
        realRat_pow_si(eps2, eps2, 2);
        realApp_set_realRat(temp, eps2, prec);
    }
//     realRat_set_si(eps2, 1, 2);
//     realRat_pow_si(eps2, eps2, 30);
    realApp_set_realRat(epsApp, eps2, prec);
    realRat_set_si(epsprime, 1, 2);
    realRat_pow_si(epsprime, epsprime, 53);
    
    
    
    /* initialized power sums */
    compApp_ptr sgNs = (compApp_ptr) ccluster_malloc ( m*sizeof(compApp) );
    for (slong g=1; g<=m; g++)
        compApp_init( sgNs + (g-1) );
    
    realApp_t absSgN, mid, rad;
    realApp_init(absSgN);
    realApp_init(mid);
    realApp_init(rad);
    
    realApp_zero(rm);
    int enoughPrec = -1;
    cauchyTest_res res;
    res.appPrec = prec;
    res.nbOfSol = enoughPrec;
    
    while ( enoughPrec == -1 ) {
        realRat_pow_si(epsprime, epsprime, 2);
        res = cauchyTest_computeSgNcompDsk( sgNs, theta, Delta, m, N, m, cacheCau, epsprime, res.appPrec, meta, depth);
        
        if (res.nbOfSol==-2)
            break;
        
        /* compute max_j( |sj/m| to the (1/j) ) */
        for (slong g=1; g<=m; g++) {
            compApp_abs(    absSgN, sgNs + (g-1), res.appPrec );
            realApp_div_si( absSgN, absSgN, m,    res.appPrec ); 
            if (realApp_contains_zero(absSgN)) {
                realApp_get_mid_realApp(mid, absSgN);
                realApp_get_rad_realApp(rad, absSgN);
                realApp_add(absSgN, mid, rad, res.appPrec);
                realApp_root_ui(absSgN,  absSgN, (ulong) g*N, res.appPrec);
                realApp_zero(mid);
                realApp_union(absSgN, mid, absSgN, res.appPrec);
            } else
                realApp_root_ui(absSgN,  absSgN, (ulong) g*N, res.appPrec);
            
            realApp_max(rm, rm, absSgN, res.appPrec);
            
//             if (metadatas_getVerbo(meta)>=level) {
//                 printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: abs of (s%ld/%ld)^(1/%ld): ", g*N, m, g*N);
//                 realApp_printd( absSgN, 10);
//                 printf("\n");
//             }
            
        }
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: max: ");
            realApp_printd( rm, 10);
            printf("\n");
            printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: wanted error: ");
            realApp_printd( epsApp, 10);
            printf("\n");
//             printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: enoughPrec: %d \n", enoughPrec);
        }
        
        realApp_get_rad_realApp( rad, rm );
        if ( realApp_ge( rad, epsApp ) >=1 ) {
            enoughPrec=-1;
            realRat_pow_si(epsprime, epsprime, 2);
            realApp_zero(rm);
        } else
            enoughPrec = 1;
    }
    
    res.nbOfSol = enoughPrec;
    
    if (! (res.nbOfSol==-2)) {
        realApp_get_mid_realApp(mid, rm);
        realApp_get_rad_realApp(rad, rm);
//         printf("mid: "); realApp_printd(mid, 10); printf("\n");
//         printf("rad: "); realApp_printd(rad, 10); printf("\n");
        realApp_add(mid, mid, rad, res.appPrec);
//         printf("mid: "); realApp_printd(mid, 10); printf("\n");
        realApp_mul(mid, mid, delta, res.appPrec);
//         printf("mid: "); realApp_printd(mid, 10); printf("\n");
        realApp_union(rm, mid, rm, res.appPrec);
        realApp_mul_realRat(rm, rm, compDsk_radiusref(Delta), res.appPrec);
    }
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: res.nbOfSol: %d \n", res.nbOfSol);
        printf("#------------turanRootRadii_root_bound.c: turanRootRadii_root_bound: max: ");
        realApp_printd( rm, 10);
        printf("\n");
    }
    
    for (slong g=1; g<=m; g++)
        compApp_clear( sgNs + (g-1) );
    ccluster_free(sgNs);
    
    realRat_clear(epsprime);
    realApp_clear(absSgN);
    realApp_clear(mid);
    realApp_clear(rad);
    
    return res.appPrec;
}
