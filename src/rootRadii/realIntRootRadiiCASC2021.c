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

#include "realIntRootRadii.h"

/* returns 0 if after some iterations, the relative accuracy is less than 1 */
int realIntRootRadiiCASC2021_Ngraeffe_iterations_inplace_real( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        int ret = 1;
        slong lastAcc=prec;
        slong curAcc;
        
        for(int i = 0; (i < N) && (ret == 1); i++) {
            curAcc = realApp_poly_get_relOne_accuracy_min(res);
            printf("#i = %d, Working precision: %ld, max relative acc: %ld, min relative acc: %ld\n", i,
                    prec, realApp_poly_get_relOne_accuracy_max(res), 
                    realApp_poly_get_relOne_accuracy_min(res)
                  );
//             if ( ( ( i == (int) N/4 ) || ( i == (int) N/2 ) || ( i == (int) 3*N/4 ) ) && (prec > CCLUSTER_DEFAULT_PREC) ) {
            if ( ( curAcc < prec/2 ) && ( lastAcc < prec/2 ) && (prec > CCLUSTER_DEFAULT_PREC) ) {
//                 printf("old prec: %ld, ", prec);
                prec = prec/2;
//                 printf("new prec: %ld \n", prec);
            }
//             if ( (lastAcc< -prec) && (curAcc<lastAcc) )
//                 ret = 0;

            if ( curAcc < 1 )
                ret = 0;
            else
                realApp_poly_oneGraeffeIteration_in_place( res, prec );
            lastAcc=curAcc;
        }
        
        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Graeffe(meta, (double) (end - start) );
            metadatas_add_time_RRGraef(meta, (double) (end - start) );
        }
        
        return ret;
        
}

/* returns 0 if after some iterations, the relative accuracy is less than 1 */
int realIntRootRadiiCASC2021_Ngraeffe_iterations_inplace_comp( compApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        int ret = 1;
        slong lastAcc=prec;
        slong curAcc;
        
        for(int i = 0; i < N && (ret == 1); i++) {
            curAcc = compApp_poly_get_relOne_accuracy_min(res);
//             printf("#Working precision: %ld, max relative acc: %ld, min relative acc: %ld\n",
//                     prec, compApp_poly_get_relOne_accuracy_max(res), 
//                     compApp_poly_get_relOne_accuracy_min(res)
//                   );
//             if ( ( ( i == (int) N/4 ) || ( i == (int) N/2 ) || ( i == (int) 3*N/4 ) )  && (prec > CCLUSTER_DEFAULT_PREC) ) {
            if ( ( curAcc < prec/2 ) && ( lastAcc < prec/2 ) && (prec > CCLUSTER_DEFAULT_PREC) ) {
//                 printf("old prec: %ld, ", prec);
                prec = prec/2;
//                 printf("new prec: %ld \n", prec);
            }
            if ( (lastAcc< -prec) && (curAcc<lastAcc) )
                ret = 0;
            else
                compApp_poly_oneGraeffeIteration_in_place( res, prec );
            
            lastAcc=curAcc;
        }
        
        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Graeffe(meta, (double) (end - start) );
            metadatas_add_time_RRGraef(meta, (double) (end - start) );
        }
        
        return ret;
}

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* returns 0 if needs more precision on the coeffs */
/* otherwise returns the length of the convex hull */
slong realIntRootRadiiCASC2021_convexHull( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec ){
    
    slong res = 0;
    /* push two first points */
    convexHull[res]=0;
    res++;
    convexHull[res]=1;
    res++;
    /* loop on other points */
    for (slong i = 2; i<len; i++){
        int liesBelow = 1;
        while ((res >= 2) && (liesBelow==1) ) {
            liesBelow = realIntRootRadii_liesBelow( convexHull[res-2], abscoeffs + convexHull[res-2], 
                                                 convexHull[res-1], abscoeffs + convexHull[res-1],
                                                 i,     abscoeffs + i, 
                                                 prec);
            if (liesBelow == 1)
                res--;
            if (liesBelow == -1) {
                res = 0;
                return res;
            }
        }
//         printf("#realIntRootRadii.c , realIntRootRadii_convexHull: res: %ld\n", res);
        convexHull[res] = i;
        res++;
    }
    
    return res;
}

slong realIntRootRadiiCASC2021_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                          slong centerRe,
                                          cacheApp_t cache,        /* polynomial */
//                                        const realRat_t delta,
                                          slong prec,
                                          metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
    int N = metadatas_getNbGIt(meta);
    ulong pow = 0x1<<N;
    realApp_t relError, relErrorInv;
    realApp_init(relError);
    realApp_init(relErrorInv);
    realApp_set_si(relError, 2*degree);
    realApp_root_ui(relError, relError, pow, prec);
    realApp_inv(relErrorInv, relError, prec);
    
//     printf("#realIntRootRadiiCASC2021.c; realIntRootRadii_rootRadiiCASC : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
        printf("#---realIntRootRadiiCASC2021.c; realIntRootRadii_rootRadii : prec: %ld \n", nprec);
        
        realIntRootRadii_getApproximation_real( pApprox, cache, nprec, meta );
        if (centerRe != 0)
            realIntRootRadii_taylor_shift_inplace_real( pApprox, centerRe, nprec, meta);
        int enoughRelacc = realIntRootRadiiCASC2021_Ngraeffe_iterations_inplace_real( pApprox, N, nprec, meta);
//         printf("---realIntRootRadiiCASC2021_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute abs of coeffs */
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++) {
                realApp_abs( (pApprox->coeffs)+i, (pApprox->coeffs)+i );
            }
            /* compute convex hull */
            lenCh = realIntRootRadiiCASC2021_convexHull( convexHull, (pApprox->coeffs), degree+1, nprec );
        }
//         printf("---realIntRootRadiiCASC2021_rootRadii: lenCh: %ld\n", lenCh);
        
        if (lenCh==0) { /* double precision */
            nprec = 2*nprec;
            nbRep++;
        }
    }
    
    if (metadatas_haveToCount(meta)) {
        if (centerRe != 0) {
                (metadatas_countref(meta))[0].RR_nbTaylors += 1;
                (metadatas_countref(meta))[0].RR_nbTaylorsRepeted += nbRep;
        }
        (metadatas_countref(meta))[0].RR_nbGraeffe += N;
        (metadatas_countref(meta))[0].RR_nbGraeffeRepeted += N*nbRep;
//         (metadatas_countref(meta))[0].RR_prec      = nprec;
//         (metadatas_countref(meta))[0].RR_predPrec      = prec;
    }
    
    if (metadatas_getVerbo(meta)>=3){
        printf("# Convex hull: %ld vertices: ", lenCh );
        for (slong ind = 0; ind < lenCh; ind++)
            printf("%ld, ", convexHull[ind]);
        printf("\n");
    }
    
    /* create list of annulii */
    compAnn_ptr cur;
    prec = CCLUSTER_DEFAULT_PREC;
    
    slong left = convexHull[0];
    for (slong ind = 1; ind < lenCh; ind++){
        
        /* create annulus */
        cur = ( compAnn_ptr ) ccluster_malloc (sizeof(compAnn));
        compAnn_init(cur);
        
        slong right = convexHull[ind];
        slong shift = right - left;
        compAnn_indMaxref(cur) = degree + 1 - (left +1);
        compAnn_indMinref(cur) = degree + 1 - (right);
        compAnn_centerReref(cur) = centerRe;
        compAnn_centerImref(cur) = 0;
        if ( realApp_contains_zero( (pApprox->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
            realApp_zero( compAnn_radSupref(cur) );
        } else {
            realApp_div( compAnn_radInfref(cur), (pApprox->coeffs) + right, (pApprox->coeffs) + left, prec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, prec );
            realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), prec );
//             ulong pow = 0x1<<N;
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, prec );
//             realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, nprec );
//             realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, nprec );
            realApp_mul( compAnn_radSupref(cur), compAnn_radInfref(cur), relError, prec );
            realApp_mul( compAnn_radInfref(cur), compAnn_radInfref(cur), relErrorInv, prec );
        }
        
        left = convexHull[ind];
//         compAnn_printd(cur, 10); printf("\n");
        compAnn_list_push(annulii, cur);
    }
    
    realApp_clear(relError);
    realApp_clear(relErrorInv);
    
    realApp_poly_clear(pApprox);
//     realRat_clear(oneplusdelta);
//     realRat_clear(oneplusdeltainv);
    ccluster_free(convexHull);
    
    return nprec;
}

slong realIntRootRadiiCASC2021_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
                                                     slong centerIm,
                                                     cacheApp_t cache,        /* polynomial */
//                                                   const realRat_t delta,
                                                     slong prec,
                                                     metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
    realRat_t oneplusdelta, oneplusdeltainv;
    
    realRat_init(oneplusdelta);
    realRat_init(oneplusdeltainv);
    
    realRat_add_si(oneplusdelta, metadatas_getRelPr(meta), 1);
    realRat_inv( oneplusdeltainv, oneplusdelta );
    
    int N = metadatas_getNbGIt(meta);
//     printf("#realIntRootRadiiCASC2021.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,degree+1);
    realApp_poly_t pSquares;
    realApp_poly_init2(pSquares,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
//         printf("#---realIntRootRadiiCASC2021.c; realIntRootRadii_rootRadii_imagCenter : prec: %ld \n", nprec);
        
        realIntRootRadii_getApproximation_comp( pApprox, cache, nprec, meta );
        if (centerIm != 0)
            realIntRootRadii_taylor_shift_inplace_comp( pApprox, 0, centerIm, nprec, meta);
        int enoughRelacc = realIntRootRadiiCASC2021_Ngraeffe_iterations_inplace_comp( pApprox, N, nprec, meta);
//         printf("---realIntRootRadiiCASC2021_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute sum of squares of real and imaginary parts of coeffs */
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++){
                realApp_sqr( compApp_realref((pApprox->coeffs)+i), compApp_realref((pApprox->coeffs)+i), nprec );
                realApp_sqr( compApp_imagref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec );
                realApp_add( (pSquares->coeffs)+i, compApp_realref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec);
            }
            /* compute convex hull */
            lenCh = realIntRootRadiiCASC2021_convexHull( convexHull, (pSquares->coeffs), degree+1, nprec );
//             printf("---realIntRootRadii_rootRadii: lenCh: %ld\n", lenCh);
        }
        
        if (lenCh==0){ /* double precision */
            nprec = 2*nprec;
            nbRep ++;
        }
    }
    
    if (metadatas_haveToCount(meta)) {
        if (centerIm != 0) {
                (metadatas_countref(meta))[0].RR_nbTaylors += 1;
                (metadatas_countref(meta))[0].RR_nbTaylorsRepeted += nbRep;
        }
        (metadatas_countref(meta))[0].RR_nbGraeffe += N;
        (metadatas_countref(meta))[0].RR_nbGraeffeRepeted += N*nbRep;
//         (metadatas_countref(meta))[0].RR_prec      = nprec;
    }
    
    for(slong i = 0; i <= degree; i++)
        realApp_sqrt( (pSquares->coeffs)+i, (pSquares->coeffs)+i, nprec);
    
//     printf("# Convex hull: %ld vertices: ", lenCh );
//     for (slong ind = 0; ind < lenCh; ind++)
//         printf("%ld, ", convexHull[ind]);
//     printf("\n");
    
    /* create list of annulii */
    compAnn_ptr cur;
    
    slong left = convexHull[0];
    for (slong ind = 1; ind < lenCh; ind++){
        
        /* create annulus */
        cur = ( compAnn_ptr ) ccluster_malloc (sizeof(compAnn));
        compAnn_init(cur);
        
        slong right = convexHull[ind];
        slong shift = right - left;
        compAnn_indMaxref(cur) = degree + 1 - (left +1);
        compAnn_indMinref(cur) = degree + 1 - (right);
        compAnn_centerReref(cur) = 0;
        compAnn_centerImref(cur) = centerIm;
        realApp_div( compAnn_radInfref(cur), (pSquares->coeffs) + right, (pSquares->coeffs) + left, nprec );
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, nprec );
        realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), nprec );
        ulong pow = 0x1<<N;
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, nprec );
        realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, nprec );
        realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, nprec );
        
        left = convexHull[ind];
//         compAnn_printd(cur, 10); printf("\n");
        compAnn_list_push(annulii, cur);
    }
    
    compApp_poly_clear(pApprox);
    realApp_poly_clear(pSquares);
    realRat_clear(oneplusdelta);
    realRat_clear(oneplusdeltainv);
    ccluster_free(convexHull);
    
    return nprec;
}

