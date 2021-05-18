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

#include "realIntRootRadii.h"

void realIntRootRadii_getApproximation_real( realApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta){
        clock_t start = clock();

        realApp_poly_set(res, cacheApp_getApproximation_real ( cache, prec ));

        if (metadatas_haveToCount(meta))
            metadatas_add_time_Approxi(meta, (double) (clock() - start) );
        
}

void realIntRootRadii_getApproximation_comp( compApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta){
        clock_t start = clock();

        compApp_poly_set(res, cacheApp_getApproximation ( cache, prec ));

        if (metadatas_haveToCount(meta))
            metadatas_add_time_Approxi(meta, (double) (clock() - start) );
        
}

void realIntRootRadii_taylor_shift_inplace_real( realApp_poly_t res, slong centerRe, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        realApp_poly_taylorShift_in_place_slong( res, 
                                           centerRe, 
                                           prec );

        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Taylors(meta, (double) (end - start) );
            metadatas_add_time_RRTaylo(meta, (double) (end - start) );
        }
}

void realIntRootRadii_taylor_shift_inplace_comp( compApp_poly_t res, slong centerRe, slong centerIm, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_taylorShift_in_place_slong( res, 
                                           centerRe,
                                           centerIm,
                                           prec );

        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Taylors(meta, (double) (end - start) );
            metadatas_add_time_RRTaylo(meta, (double) (end - start) );
        }
}

/* returns 0 if after some iterations, the relative accuracy is less than 1 */
int realIntRootRadii_Ngraeffe_iterations_inplace_real( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        int ret = 1;
        slong lastAcc=prec;
        slong curAcc;
        for(int i = 0; (i < N) && (ret == 1); i++) {
            curAcc = realApp_poly_get_relOne_accuracy_min(res);
//             printf("#Working precision: %ld, max relative acc: %ld, min relative acc: %ld\n",
//                     prec, realApp_poly_get_relOne_accuracy_max(res), 
//                     realApp_poly_get_relOne_accuracy_min(res)
//                   );
//             printf("#Working precision: %ld, max absolute acc: %ld, min absolute acc: %ld\n",
//                     prec, realApp_poly_get_absolute_accuracy_max(res), 
//                     realApp_poly_get_absolute_accuracy_min(res)
//                   );
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
int realIntRootRadii_Ngraeffe_iterations_inplace_comp( compApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        int ret = 1;
        slong lastAcc=prec;
        slong curAcc;
        
        for(int i = 0; i < N && (ret == 1); i++) {
            curAcc = compApp_poly_get_relOne_accuracy_min(res);
//             printf("Working precision: %ld, max relative acc: %ld, min relative acc: %ld\n",
//                     prec, compApp_poly_get_relOne_accuracy_max(res), compApp_poly_get_relOne_accuracy_min(res));
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

/* assume i<j<k */
/* assume absPi=|pi|, absPj=|pj|, absPk=|pk| are approximations of integers */
/* decide if [j,log|pj|] lies below the line passing trough [i,log|pi|] and [k,log|pk|]*/
/* returns 1 if yes */
/*         0 if no  */
/*        -1 if it can not be decided */
int realIntRootRadii_liesBelow( slong i, const realApp_t absPi,
                             slong j, const realApp_t absPj,
                             slong k, const realApp_t absPk,
                             slong prec ){
    /* decide if [j,log|pj|] lies below the line passing trough [i,log|pi|] and [k,log|pk|]  */
    /*  <=> decide if (log|pj| - log|pi|)/(j-i) <= (log|pk| - log|pi|)/(k-i)                 */
    /*  <=> decide if (|pj|/|pi|)^(1/(j-i)) <= (|pk|/|pi|)^(1/(k-i))                         */
    /*  <=> decide if (|pj|/|pi|)^(k-i) <= (|pk|/|pi|)^(j-i)                                 */
    /*  <=> decide if |pj|^(k-i) * |pi|^(j-i) <= |pk|^(j-i) * |pi|^(k-i)                     */
    /*  <=> decide if |pj|^(k-i) * |pi|^(j-i) - |pk|^(j-i) * |pi|^(k-i) <= 0                 */
    /*                left side of inequality is integer                                     */
    
//     printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow: i: %ld, j:%ld, k:%ld \n", i, j, k);
    int res=-1;
    /* check if absPi is 0 */
    /* in which case j lies ABOVE the line passing trough [i,log|pi|] and [k,log|pk|]*/
    if ( realApp_contains_zero(absPi) ){
        realApp_t width, half;
        realApp_init( width );
        realApp_init( half );
        realApp_get_rad_realApp(width, absPi);
        realApp_set_d(half, 0.5);
        if (realApp_lt(width, half)==1) {
            /* with strictly less that 1, contains as unique integer 0*/
            res = 0;
//             printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow: absPi is zero \n");
        }
        else {/* otherwise can not say if it is zero or not */
            res = -1;
//             printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow: absPi may be zero, but need more prec \n");
        }
        
        realApp_clear( width );
        realApp_clear( half );
        return res;
    }
    
    ulong kmi = k-i;
    ulong jmi = j-i;
    realApp_t leftSide, rightSide, temp;
    realApp_init(leftSide);
    realApp_init(rightSide);
    realApp_init(temp);
    
    realApp_pow_ui(leftSide,  absPj,     kmi,       prec);
    realApp_pow_ui(temp,      absPi,     jmi,       prec);
    realApp_mul   (leftSide,  leftSide,  temp,      prec);
    realApp_pow_ui(rightSide, absPk,     jmi,       prec);
    realApp_pow_ui(temp,      absPi,     kmi,       prec);
    realApp_mul   (rightSide, rightSide, temp,      prec);
    realApp_sub   (leftSide,  leftSide,  rightSide, prec);
    
//     printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow: leftSide: ");
//     realApp_printd(leftSide, 10);
//     printf("\n");
    
    realApp_zero  (temp);
    if (realApp_lt(leftSide, temp)==1)
        res = 1;
    else if (realApp_gt(leftSide, temp)==1)
        res = 0;
    else {
             realApp_get_rad_realApp(leftSide, leftSide);
             realApp_one(temp);
             realApp_div_ui(temp, temp, 2, prec);
             if (realApp_lt(leftSide, temp)==1) 
                 /* with strictly less that 1, contains as unique integer 0*/
                 res = 1;
             else /* otherwise can not say if it is zero or not */
                 res = -1;
    }
    
    realApp_clear(leftSide);
    realApp_clear(rightSide);
    realApp_init(temp);

//     printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow: leq 0: %d \n", res);
    
    return res;
}

/* assume i<j<k */
/* assume logAbsPi=log|pi|, logAbsPj=log|pj|, logAbsPk=log|pk| have absolute error less than 1/2 */
/* decide if [j,log|pj|] lies below the line passing trough [i,log|pi|] and [k,log|pk|]*/
/* returns 1 if yes */
/*         0 if no  */
/*        -1 if it can not be decided => should never happen*/
int realIntRootRadii_liesBelow_2( slong i, const realApp_t logAbsPi,
                             slong j, const realApp_t logAbsPj,
                             slong k, const realApp_t logAbsPk,
                             slong prec ){
    /* decide if [j,log|pj|] lies below the line passing trough [i,log|pi|] and [k,log|pk|]  */
    /*  <=> decide if (log|pj| - log|pi|)/(j-i) <= (log|pk| - log|pi|)/(k-i)                 */
    /*  <=> decide if (log|pj| - log|pi|)*(k-i) - (log|pk| - log|pi|)*(j-i) <= 0             */
    
//     printf( "#realIntRootRadii.c , realIntRootRadii_liesBelow_2: i: %ld, j:%ld, k:%ld \n", i, j, k);
    int res=-1;
    /* check if absPi is <0 */
    /* in which case j lies ABOVE the line passing trough [i,log|pi|] and [k,log|pk|]*/
    if ( realApp_is_negative(logAbsPi) || realApp_is_negative(logAbsPj) || realApp_is_negative(logAbsPk) ){
        return 0;
    }
    
    slong tempPrec = prec;
    int isExact  = 0;
    int isStrLe  = 0;
    int isStrGr  = 0;
    int isEqZer  = 0;
    int decided  = (isStrGr || isStrLe || (isExact && isEqZer));
    
    realApp_t leftSide, rightSide;
    realApp_init(leftSide);
    realApp_init(rightSide);
    
    while (!decided){
        realApp_sub(leftSide, logAbsPj, logAbsPi, tempPrec);
        realApp_mul_si(leftSide, leftSide, k-i, tempPrec);
        realApp_sub(rightSide, logAbsPk, logAbsPi, tempPrec);
        realApp_mul_si(rightSide, rightSide, j-i, tempPrec);
        realApp_sub(leftSide, leftSide, rightSide, tempPrec);
        isExact = realApp_is_exact(leftSide);
        isStrLe = realApp_is_negative(leftSide);
        isStrGr = realApp_is_positive(leftSide);
        isEqZer = realApp_is_zero(leftSide);
        decided  = (isStrGr || isStrLe || (isExact && isEqZer));
        if (!decided)
            tempPrec = 2*tempPrec;
    }
    
    
    realApp_clear(leftSide);
    realApp_clear(rightSide);
    
    if (isStrGr)
        res = 0;
    if (isStrLe || (isExact && isEqZer))
        res = 1;
    
    return res;
}

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* returns 0 if needs more precision on the coeffs */
/* otherwise returns the length of the convex hull */
slong realIntRootRadii_convexHull( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec ){
    
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
                res=0;
                return res;
            }
        }
//         printf("#realIntRootRadii.c , realIntRootRadii_convexHull: res: %ld\n", res);
        convexHull[res] = i;
        res++;
    }
    
    return res;
}

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* returns 0 if needs more precision on the coeffs */
/* otherwise returns the length of the convex hull */
slong realIntRootRadii_convexHull_2( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec ){
    
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
            liesBelow = realIntRootRadii_liesBelow_2( convexHull[res-2], abscoeffs + convexHull[res-2], 
                                                 convexHull[res-1], abscoeffs + convexHull[res-1],
                                                 i,     abscoeffs + i, 
                                                 prec);
            if (liesBelow == 1)
                res--;
            if (liesBelow == -1) {
                res=0;
                return res;
            }
        }
//         printf("#realIntRootRadii.c , realIntRootRadii_convexHull: res: %ld\n", res);
        convexHull[res] = i;
        res++;
    }
    
    return res;
}

slong realIntRootRadii_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                                  cacheApp_t cache,        /* polynomial */
//                                   const realRat_t delta,
                                  slong prec,
                                  metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
//     realRat_t oneplusdelta, oneplusdeltainv;
//     
//     realRat_init(oneplusdelta);
//     realRat_init(oneplusdeltainv);
//     
//     realRat_add_si(oneplusdelta, delta, 1);
//     realRat_inv( oneplusdeltainv, oneplusdelta );
//     
//     double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
//     log2_1pdelta = log2_1pdelta / log(2);
//     int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    
    int N = metadatas_getNbGIt(meta);
    ulong pow = 0x1<<N;
    /* test */
    realApp_t relError, relErrorInv;
    realApp_init(relError);
    realApp_init(relErrorInv);
    realApp_set_si(relError, 2*degree);
    realApp_root_ui(relError, relError, pow, prec);
    realApp_inv(relErrorInv, relError, prec);
    /* fin test*/
    
//     printf("#realIntRootRadii.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
        realIntRootRadii_getApproximation_real( pApprox, cache, nprec, meta );
//         realApp_poly_set(pApprox, cacheApp_getApproximation_real ( cache, nprec ));
        if (centerRe != 0)
            realIntRootRadii_taylor_shift_inplace_real( pApprox, centerRe, nprec, meta);
        int enoughRelacc = realIntRootRadii_Ngraeffe_iterations_inplace_real( pApprox, N, nprec, meta);
//         printf("realIntRootRadii_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute abs of coeffs */
        
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++) {
                realApp_abs( (pApprox->coeffs)+i, (pApprox->coeffs)+i );
            }
            /* compute convex hull */
            lenCh = realIntRootRadii_convexHull( convexHull, (pApprox->coeffs), degree+1, nprec );
        }
        
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
        (metadatas_countref(meta))[0].RR_prec      = nprec;
        (metadatas_countref(meta))[0].RR_predPrec      = prec;
    }
    
//     if (metadatas_getVerbo(meta)>=3){
//         printf("# Convex hull: %ld vertices: ", lenCh );
//         for (slong ind = 0; ind < lenCh; ind++)
//             printf("%ld, ", convexHull[ind]);
//         printf("\n");
//     }
    
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
        compAnn_centerReref(cur) = centerRe;
        compAnn_centerImref(cur) = 0;
        if ( realApp_contains_zero( (pApprox->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
            realApp_zero( compAnn_radSupref(cur) );
        } else {
            realApp_div( compAnn_radInfref(cur), (pApprox->coeffs) + right, (pApprox->coeffs) + left, nprec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, nprec );
            realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), nprec );
//             ulong pow = 0x1<<N;
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, nprec );
//             realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, nprec );
//             realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, nprec );
            realApp_mul( compAnn_radSupref(cur), compAnn_radInfref(cur), relError, nprec );
            realApp_mul( compAnn_radInfref(cur), compAnn_radInfref(cur), relErrorInv, nprec );
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

slong realIntRootRadii_rootRadii_2( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                                  cacheApp_t cache,        /* polynomial */
//                                   const realRat_t delta,
                                  slong prec,
                                  metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
//     realRat_t oneplusdelta, oneplusdeltainv;
//     
//     realRat_init(oneplusdelta);
//     realRat_init(oneplusdeltainv);
//     
//     realRat_add_si(oneplusdelta, delta, 1);
//     realRat_inv( oneplusdeltainv, oneplusdelta );
//     
//     double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
//     log2_1pdelta = log2_1pdelta / log(2);
//     int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    
    int N = metadatas_getNbGIt(meta);
    ulong pow = 0x1<<N;
    /* test */
    realApp_t relError, relErrorInv;
    realApp_init(relError);
    realApp_init(relErrorInv);
    realApp_set_si(relError, 2*degree);
    realApp_root_ui(relError, relError, pow, prec);
    realApp_inv(relErrorInv, relError, prec);
    /* fin test*/
    
//     printf("#realIntRootRadii.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    realApp_poly_t pApprox, logPApprox;
    realApp_poly_init2(pApprox,degree+1);
    realApp_poly_init2(logPApprox,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
        realIntRootRadii_getApproximation_real( pApprox, cache, nprec, meta );
//         realApp_poly_set(pApprox, cacheApp_getApproximation_real ( cache, nprec ));
        if (centerRe != 0)
            realIntRootRadii_taylor_shift_inplace_real( pApprox, centerRe, nprec, meta);
        int enoughRelacc = realIntRootRadii_Ngraeffe_iterations_inplace_real( pApprox, N, nprec, meta);
//         printf("realIntRootRadii_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute abs of coeffs */
        
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++) {
                realApp_abs( (pApprox->coeffs)+i, (pApprox->coeffs)+i );
                realApp_log_base_ui((logPApprox->coeffs)+i, (pApprox->coeffs)+i, 2, nprec);
                realApp_set_rad_zero((logPApprox->coeffs)+i);
//                 printf("realIntRootRadii_rootRadii: %d-th abs of coeff: ", i);
//                 realApp_printd((pApprox->coeffs)+i, 10);
//                 printf("\nrealIntRootRadii_rootRadii: %d-th log abs of coeff: ", i);
//                 realApp_printd((logPApprox->coeffs)+i, 10);
//                 printf("\n");
            }
            /* compute convex hull */
            lenCh = realIntRootRadii_convexHull_2( convexHull, (logPApprox->coeffs), degree+1, nprec );
        }
        
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
        (metadatas_countref(meta))[0].RR_prec      = nprec;
        (metadatas_countref(meta))[0].RR_predPrec      = prec;
    }
    
    if (metadatas_getVerbo(meta)>=3){
        printf("# Convex hull: %ld vertices: ", lenCh );
        for (slong ind = 0; ind < lenCh; ind++)
            printf("%ld, ", convexHull[ind]);
        printf("\n");
    }
    
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
        compAnn_centerReref(cur) = centerRe;
        compAnn_centerImref(cur) = 0;
        if ( realApp_contains_zero( (pApprox->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
            realApp_zero( compAnn_radSupref(cur) );
        } else {
            realApp_div( compAnn_radInfref(cur), (pApprox->coeffs) + right, (pApprox->coeffs) + left, nprec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, nprec );
            realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), nprec );
//             ulong pow = 0x1<<N;
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, nprec );
//             realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, nprec );
//             realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, nprec );
            realApp_mul( compAnn_radSupref(cur), compAnn_radInfref(cur), relError, nprec );
            realApp_mul( compAnn_radInfref(cur), compAnn_radInfref(cur), relErrorInv, nprec );
        }
        
        left = convexHull[ind];
//         compAnn_printd(cur, 10); printf("\n");
        compAnn_list_push(annulii, cur);
    }
    
    realApp_clear(relError);
    realApp_clear(relErrorInv);
    
    realApp_poly_clear(pApprox);
    realApp_poly_clear(logPApprox);
//     realRat_clear(oneplusdelta);
//     realRat_clear(oneplusdeltainv);
    ccluster_free(convexHull);
    
    return nprec;
}

slong realIntRootRadii_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
                                             slong centerIm,
                                             cacheApp_t cache,        /* polynomial */
//                                              const realRat_t delta,
                                             slong prec,
                                             metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
    realRat_t oneplusdelta, oneplusdeltainv;
    
    realRat_init(oneplusdelta);
    realRat_init(oneplusdeltainv);
    
    realRat_add_si(oneplusdelta, metadatas_getRelPr(meta), 1);
    realRat_inv( oneplusdeltainv, oneplusdelta );
    
//     double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
//     log2_1pdelta = log2_1pdelta / log(2);
//     int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    
    int N = metadatas_getNbGIt(meta);
//     printf("#realIntRootRadii.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,degree+1);
    realApp_poly_t pSquares;
    realApp_poly_init2(pSquares,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
        realIntRootRadii_getApproximation_comp( pApprox, cache, nprec, meta );
//         compApp_poly_set(pApprox, cacheApp_getApproximation ( cache, nprec ));
        if (centerIm != 0)
            realIntRootRadii_taylor_shift_inplace_comp( pApprox, 0, centerIm, nprec, meta);
        int enoughRelacc = realIntRootRadii_Ngraeffe_iterations_inplace_comp( pApprox, N, nprec, meta);
//         for(int i = 0; i < N; i++)
//             compApp_poly_oneGraeffeIteration_in_place( pApprox, nprec );
        /* compute sum of squares of real and imaginary parts of coeffs */
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++){
                realApp_sqr( compApp_realref((pApprox->coeffs)+i), compApp_realref((pApprox->coeffs)+i), nprec );
                realApp_sqr( compApp_imagref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec );
                realApp_add( (pSquares->coeffs)+i, compApp_realref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec);
            }
            /* compute convex hull */
            lenCh = realIntRootRadii_convexHull( convexHull, (pSquares->coeffs), degree+1, nprec );
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
        (metadatas_countref(meta))[0].RR_prec      = nprec;
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

void realIntRootRadii_connectedComponents( compAnn_list_t annulii, slong prec ){
    
    compAnn_ptr cur, curnext;
    compAnn_list_iterator it, itnext;
    
    /* group annulii into connected components */
    
    it = compAnn_list_begin(annulii);
    cur = compAnn_list_elmt( it ); /* annulii contains at least one element */
    
    itnext = compAnn_list_next(it);
    
    while (itnext!=compAnn_list_end() ) {
        curnext = compAnn_list_elmt( itnext );
        
        if (! (realApp_lt(compAnn_radSupref(cur), compAnn_radInfref(curnext))==1)) {
            /* merge the two annulii into cur */
            compAnn_indMinref(cur) = compAnn_indMinref(curnext);
            realApp_set( compAnn_radSupref(cur), compAnn_radSupref(curnext) );
            /*remove curnext from the list*/
            curnext = compAnn_list_remove_at(annulii, it);
            /*delete curnext*/
            compAnn_clear( curnext );
            ccluster_free( curnext );
            itnext = compAnn_list_next(it);
        } else {
            it = compAnn_list_next(it);
            cur = compAnn_list_elmt( it );
            itnext = compAnn_list_next(itnext);
        }
        
    }
    
}

void realIntRootRadii_containsRealRoot( compAnn_list_t annulii, cacheApp_t cache, slong prec ){
    
    compAnn_ptr cur;
    compAnn_list_iterator it;
    
    slong degree = cacheApp_getDegree(cache);
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    //     realIntRootRadii_getApproximation( pApprox, cache, prec, meta );
    realApp_poly_set(pApprox, cacheApp_getApproximation_real ( cache, prec ));
    
    realApp_t centerLeft, centerLeftVal, centerRight, centerRightVal;
    realApp_init(centerLeft);
    realApp_init(centerLeftVal);
    realApp_init(centerRight);
    realApp_init(centerRightVal);
    
//     realApp_poly_t pApproxDer;
//     realApp_poly_init2(pApproxDer,degree);
//     realApp_poly_derivative(pApproxDer, pApprox);
    realApp_t interval;
    realApp_init(interval);
    
    it = compAnn_list_begin(annulii);
    while (it!=compAnn_list_end() ) {
        cur = compAnn_list_elmt( it );
        /* check if is is 0 */
        if ( realApp_is_zero(compAnn_radInfref(cur)) && realApp_is_zero(compAnn_radSupref(cur)) ){
            compAnn_rrInPoref(cur) = 1;
            compAnn_rrInNeref(cur) = 1;
        } else
        if ( compAnn_indMinref(cur) == compAnn_indMaxref(cur) ) { 
            /*contains a unique root => it is real (otherwise would contain 2)*/
            /* intersection with R+ */
            realApp_get_mid_realApp( centerLeft, compAnn_radInfref(cur) );
            realApp_get_mid_realApp( centerRight, compAnn_radSupref(cur) );
            realApp_poly_evaluate(centerLeftVal, pApprox, centerLeft, prec);
            realApp_poly_evaluate(centerRightVal, pApprox, centerRight, prec);
//             printf("Annulus %ld, %ld: \n", compAnn_indMinref(cur), compAnn_indMinref(cur) );
//             printf("--Value on Left : "); realApp_printd(centerLeftVal, 10); printf("\n");
//             printf("--Value on Right: "); realApp_printd(centerRightVal, 10); printf("\n");
            int prodSgn = realApp_sgn_nonzero(centerLeftVal) * realApp_sgn_nonzero(centerRightVal);
            if ( prodSgn == -1 ) { /*opposite sign, contains a unique real root*/
                compAnn_rrInPoref(cur) = 1;
                compAnn_rrInNeref(cur) = 0;
            } else if ( prodSgn == 1 ) { 
                /*same sign, contains no real root => intersection with R- contains a unique real root*/
                compAnn_rrInPoref(cur) = 0;
                compAnn_rrInNeref(cur) = 1;
            } else { /* ( prodSgn == 0 ) one evaluation contains zero -> try the same for negative segment*/
                realApp_neg( centerLeft, centerLeft );
                realApp_neg( centerRight, centerRight ); /* no need to swap */
                realApp_poly_evaluate(centerLeftVal, pApprox, centerLeft, prec);
                realApp_poly_evaluate(centerRightVal, pApprox, centerRight, prec);
                prodSgn = realApp_sgn_nonzero(centerLeftVal) * realApp_sgn_nonzero(centerRightVal);
                if ( prodSgn == -1 ) { /*opposite sign, contains a unique real root*/
                    compAnn_rrInPoref(cur) = 0;
                    compAnn_rrInNeref(cur) = 1;
                } else if ( prodSgn == 1 ) { 
                    /*same sign, contains no real root => intersection with R+ contains a unique real root*/
                    compAnn_rrInPoref(cur) = 1;
                    compAnn_rrInNeref(cur) = 0;
                }
                /* else can not decide */
            }
        } else if ( compAnn_indMinref(cur) == (compAnn_indMaxref(cur) -1) ) {
            /* contains 2 roots */
            /* intersection with R+ */
            realApp_get_mid_realApp( centerLeft, compAnn_radInfref(cur) );
            realApp_get_mid_realApp( centerRight, compAnn_radSupref(cur) );
            realApp_poly_evaluate(centerLeftVal, pApprox, centerLeft, prec);
            realApp_poly_evaluate(centerRightVal, pApprox, centerRight, prec);
//             printf("Annulus %ld, %ld: \n", compAnn_indMinref(cur), compAnn_indMinref(cur) );
//             printf("--Value on Left : "); realApp_printd(centerLeftVal, 10); printf("\n");
//             printf("--Value on Right: "); realApp_printd(centerRightVal, 10); printf("\n");
            int prodSgn = realApp_sgn_nonzero(centerLeftVal) * realApp_sgn_nonzero(centerRightVal);
            if ( prodSgn == -1 ) { /* opposite sign, contains a unique real root*/
                                   /* it can not be double otherwise signs would be the same */
                                   /* intersection with R- also contains a unique real root */
                compAnn_rrInPoref(cur) = 1;
                compAnn_rrInNeref(cur) = 1;
            }
        } else if ( compAnn_indMinref(cur) < (compAnn_indMaxref(cur) -1) ) {
//             printf("index min: %d, index max: %d \n", compAnn_indMinref(cur), compAnn_indMaxref(cur));
            /* intersection with R+ */
            realApp_get_mid_realApp( centerLeft, compAnn_radInfref(cur) );
            realApp_get_mid_realApp( centerRight, compAnn_radSupref(cur) );
            realApp_poly_evaluate(centerLeftVal, pApprox, centerLeft, prec);
            realApp_poly_evaluate(centerRightVal, pApprox, centerRight, prec);
//             printf("Annulus %ld, %ld: \n", compAnn_indMinref(cur), compAnn_indMinref(cur) );
//             printf("--Value on Left : "); realApp_printd(centerLeftVal, 10); printf("\n");
//             printf("--Value on Right: "); realApp_printd(centerRightVal, 10); printf("\n");
            int prodSgn = realApp_sgn_nonzero(centerLeftVal) * realApp_sgn_nonzero(centerRightVal);
            if ( prodSgn == -1 ) { /* opposite sign, contains at least one real root*/
                compAnn_rrInPoref(cur) = 2;
            }
            /* intersection with R- */
            realApp_neg( centerLeft, centerLeft );
            realApp_neg( centerRight, centerRight ); /* no need to swap */
            realApp_poly_evaluate(centerLeftVal, pApprox, centerLeft, prec);
            realApp_poly_evaluate(centerRightVal, pApprox, centerRight, prec);
//             printf("Annulus %ld, %ld: \n", compAnn_indMinref(cur), compAnn_indMinref(cur) );
//             printf("--Value on Left : "); realApp_printd(centerLeftVal, 10); printf("\n");
//             printf("--Value on Right: "); realApp_printd(centerRightVal, 10); printf("\n");
            prodSgn = realApp_sgn_nonzero(centerLeftVal) * realApp_sgn_nonzero(centerRightVal);
            if ( prodSgn == -1 ) { /* opposite sign, contains at least one real root*/
                compAnn_rrInNeref(cur) = 2;
            }
            
            /* else undetermined */
        }
        /* try interval evaluations */
//         else {
        if (compAnn_rrInPoref(cur) == -1 ) {
            realApp_get_mid_realApp( centerLeft, compAnn_radInfref(cur) );
            realApp_get_mid_realApp( centerRight, compAnn_radSupref(cur) );
            realApp_union(interval, centerLeft, centerRight, prec);
            realApp_poly_evaluate(interval, pApprox, interval, prec);
            if (realApp_contains_zero(interval)==0)
                compAnn_rrInPoref(cur) = 0;
        }
        if (compAnn_rrInNeref(cur) == -1 ) {
            realApp_neg( centerLeft, centerLeft );
            realApp_neg( centerRight, centerRight ); /* no need to swap */
            realApp_union(interval, centerLeft, centerRight, prec);
            realApp_poly_evaluate(interval, pApprox, interval, prec);
            if (realApp_contains_zero(interval)==0)
                compAnn_rrInNeref(cur) = 0;
        }
//         }
        it = compAnn_list_next(it);
    }
    
    realApp_poly_clear(pApprox);
//     realApp_poly_clear(pApproxDer);
    realApp_clear(interval);
    
    realApp_clear(centerLeft);
    realApp_clear(centerLeftVal);
    realApp_clear(centerRight);
    realApp_clear(centerRightVal);
    
}

void realIntRootRadii_bisect_connCmp( connCmp_list_t dest, 
                                      connCmp_t cc){
    
    compBox_list_t subBoxes;
    compBox_list_init(subBoxes);
    
    compBox_ptr btemp, bstemp;
    
    while (!connCmp_is_empty(cc)) {
        btemp = connCmp_pop(cc);
        /* bisect */
        subdBox_risolate_bisect( subBoxes, btemp );
        /* remove boxes that do not intersect an annulus */
        while (!compBox_list_is_empty(subBoxes)) {
            bstemp = compBox_list_pop(subBoxes);
            compBox_actualize_anulii_risolate( bstemp, btemp );
            
            int nbSol=-1;
            if ( compAnn_list_get_size(compBox_annuli0ref(bstemp)) == 0 )
                nbSol = 0;
            if ( compAnn_list_get_size(compBox_annuli0ref(bstemp)) == 1 ) {
                realApp_t center;
                realApp_t left, right, rad;
                realApp_init(center);
                realApp_init(left);
                realApp_init(right);
                realApp_init(rad);
                
                realApp_set_realRat( center, compRat_realref(compBox_centerref(bstemp)), CCLUSTER_DEFAULT_PREC );
                realApp_set_realRat( rad,    compBox_bwidthref(bstemp), CCLUSTER_DEFAULT_PREC );
                realApp_div_si     ( rad,    rad,                  2, CCLUSTER_DEFAULT_PREC );
                realApp_sub        ( left,   center,                 rad, CCLUSTER_DEFAULT_PREC );
                realApp_add        ( right,  center,                 rad, CCLUSTER_DEFAULT_PREC );
                
                compAnn_ptr ann = compAnn_list_first( compBox_annuli0ref(bstemp) ); /* there is only one element in the list */
                
                if ( ( (realApp_is_positive(left)==1)  && (compAnn_rrInPoref(ann) == 0) ) ||
                     ( (realApp_is_negative(right)==1)  && (compAnn_rrInNeref(ann) == 0) ) )
                    nbSol=0;
                
                realApp_clear(center);
                realApp_clear(left);
                realApp_clear(right);
                realApp_clear(rad);
            }
            
            if (nbSol==0) { /* delete bstemp */
                compBox_clear(bstemp);
                ccluster_free(bstemp);
            } else {
                connCmp_union_compBox( dest, bstemp);
            }
        }
        
        /* delete btemp */
        compBox_clear(btemp);
        ccluster_free(btemp);
        
        
    }
    
    compBox_list_clear(subBoxes);
}
