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
        }
        else {/* otherwise can not say if it is zero or not */
            res = -1;
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
             else {/* otherwise can not say if it is zero or not */
                 res = -1;
             }
    }
    
    realApp_clear(leftSide);
    realApp_clear(rightSide);
    realApp_init(temp);
    
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
                /* it is not possible to decide if (res-1, convexHull[res-1]) lies below     */
                /* the line passing trough (res-2, convexHull[res-2]) and (i, convexHull[i]) */
                /* try to figure out if there exist a k s.t. it is possible to decide that  */
                /* it lies below line passing trough (res-2, convexHull[res-2]) and (k, convexHull[k]) */
//                 printf("#%ld-th point lies below the line passing trough %ld, %ld ?\n", res-1, res-2, i);
                int liesBelow2 = -1;
                for (slong k = i+1; (k<len)&&(liesBelow2<1); k++) {
                    liesBelow2 = realIntRootRadii_liesBelow( convexHull[res-2], abscoeffs + convexHull[res-2], 
                                                 convexHull[res-1], abscoeffs + convexHull[res-1],
                                                 k,     abscoeffs + k, 
                                                 prec);
                        
                }
                if (liesBelow2 == 1) {
                    liesBelow = 1;
                    res--;
                } else {
                    res = 0;
                    return res;
                }
            }
        }
        convexHull[res] = i;
        res++;
    }
    
    return res;
}

int   realIntRootRadii_GraeffeAndCH_real ( slong convexHull[], slong * lenCH, slong * nprec, realApp_poly_t absCoeffs,
                                           realApp_poly_t pApprox, int N, slong prec, metadatas_t meta ) {
    
    int level = 5;
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#realIntRootRadii_GraeffeAndCH_real: begin; N = %d, precision: %ld\n", 
                                                                     N,            prec);
    }
        
    slong degree = realApp_poly_degree(pApprox);
     
    slong lastAcc = prec, curAcc;
    
    *lenCH = 1;
    int i = 1;
    for( ; ((i <= N) && (*lenCH)) ; i++) {
        
        curAcc = realApp_poly_get_relOne_accuracy_min(pApprox);
        if (metadatas_getVerbo(meta)>=level) {
            printf("#i = %d, working precision: %ld, last min relative acc: %ld, current min relative accuracy: %ld\n", 
                          i,                    prec,                      lastAcc,                            curAcc);
        }
        
        /* try to round to lower precision to save on bit operations */
        if ( ( curAcc < prec/2 ) && ( lastAcc < prec/2 ) && (prec > CCLUSTER_DEFAULT_PREC) ) {
            prec = prec/2;
            
            if (metadatas_getVerbo(meta)>=level) {
                printf("#i = %d, new working precision: %ld\n", 
                              i,                        prec  );
            }
        
        }
        
        /* apply one graeffe iteration to precision prec */
        clock_t start = clock();
        realApp_poly_oneGraeffeIteration_in_place( pApprox, prec );
        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Graeffe(meta, (double) (end - start) );
            metadatas_add_time_RRGraef(meta, (double) (end - start) );
        }
        
        /* if current minimum accuracy is less than 1,
         * try to compute the convexHull 
         * if N graeffe iterations have been computed,
         * try to compute the convexHull */
        if ( (curAcc <=1) || (i==N) ){
            /* first compute absolute values of coeffs */
            for(slong j = 0; j <= degree; j++)
                realApp_abs( (absCoeffs->coeffs)+j, (pApprox->coeffs)+j );
            /* then compute convex hull */
            *lenCH = realIntRootRadii_convexHull( convexHull, (absCoeffs->coeffs), degree + 1, prec );
            if (metadatas_getVerbo(meta)>=level) {
                printf("#i = %d, length of computed convex hull: %ld\n", 
                              i,                                 *lenCH);
            }
            if (*lenCH==0)
                i--;
        }
        
        /* update lastAcc */
        lastAcc=curAcc;
    }
    i=i-1;
    
    *nprec = prec;
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#realIntRootRadii_GraeffeAndCH_real: end; i = %d, lenCH = %ld, nprec = %ld\n",
                                                                   i,           *lenCH,      *nprec);
    }
    return i;
}

int   realIntRootRadii_GraeffeAndCH_comp ( slong convexHull[], slong * lenCH, slong * nprec, realApp_poly_t absCoeffs,
                                           compApp_poly_t pApprox, int N, slong prec, metadatas_t meta ) {
    
    int level = 5;
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#realIntRootRadii_GraeffeAndCH_comp: begin; N = %d, precision: %ld\n", 
                                                                     N,            prec);
    }
        
    slong degree = compApp_poly_degree(pApprox);
    compApp_poly_t pSquares;
    compApp_poly_init2(pSquares, degree +1);
    slong lastAcc = prec, curAcc;
    
    *lenCH = 1;
    int i = 1;
    for( ; ((i <= N) && (*lenCH)) ; i++) {
        
        curAcc = compApp_poly_get_relOne_accuracy_min(pApprox);
        if (metadatas_getVerbo(meta)>=level) {
            printf("#i = %d, working precision: %ld, last min relative acc: %ld, current min relative accuracy: %ld\n", 
                          i,                    prec,                      lastAcc,                            curAcc);
        }
        
        /* try to round to lower precision to save on bit operations */
        if ( ( curAcc < prec/2 ) && ( lastAcc < prec/2 ) && (prec > CCLUSTER_DEFAULT_PREC) ) {
            prec = prec/2;
            
            if (metadatas_getVerbo(meta)>=level) {
                printf("#i = %d, new working precision: %ld\n", 
                              i,                        prec  );
            }
        
        }
        
        /* apply one graeffe iteration to precision prec */
        clock_t start = clock();
        compApp_poly_oneGraeffeIteration_in_place( pApprox, prec );
        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Graeffe(meta, (double) (end - start) );
            metadatas_add_time_RRGraef(meta, (double) (end - start) );
        }
        
        /* if current minimum accuracy is less than 1,
         * try to compute the convexHull 
         * if N graeffe iterations have been computed,
         * try to compute the convexHull */
        if ( (curAcc <=1) || (i==N) ){
            /* first compute square of modulii of coeffs, which are integers */
            for(slong j = 0; j <= degree; j++) {
                realApp_sqr( compApp_realref((pSquares->coeffs)+j), compApp_realref((pApprox->coeffs)+j), prec );
                realApp_sqr( compApp_imagref((pSquares->coeffs)+j), compApp_imagref((pApprox->coeffs)+j), prec );
                realApp_add( (absCoeffs->coeffs)+j, compApp_realref((pSquares->coeffs)+j), compApp_imagref((pSquares->coeffs)+j), prec);
            }
            /* then compute convex hull */
            *lenCH = realIntRootRadii_convexHull( convexHull, (absCoeffs->coeffs), degree + 1, prec );
            if (metadatas_getVerbo(meta)>=level) {
                printf("#i = %d, length of computed convex hull: %ld\n", 
                              i,                                 *lenCH);
            }
            if (*lenCH==0)
                i--;
        }
        
        /* update lastAcc */
        lastAcc=curAcc;
    }
    i=i-1;
    
    *nprec = prec;
    
    if (metadatas_getVerbo(meta)>=level) {
            printf("#realIntRootRadii_GraeffeAndCH_comp: end; i = %d, lenCH = %ld, nprec = %ld\n",
                                                                   i,           *lenCH,      *nprec);
    }
    compApp_poly_clear(pSquares);
    
    return i;
}

slong realIntRootRadii_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                                  cacheApp_t cache,        /* polynomial */
                                  slong prec,
                                  metadatas_t meta ){
    
    int level = 4;
    slong degree = cacheApp_getDegree(cache);
    
    int N = metadatas_getNbGIt(meta);
    ulong pow = 0x1<<N;
    realApp_t relError, relErrorInv;
    realApp_init(relError);
    realApp_init(relErrorInv);
    realApp_set_si(relError, 2*degree);
    realApp_root_ui(relError, relError, pow, prec);
    realApp_inv(relErrorInv, relError, prec);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    slong nnprec = prec;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    realApp_poly_t absCoeffs;
    realApp_poly_init2(absCoeffs,degree+1);
    
    while ( lenCh == 0 ) {
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---realIntRootRadii.c; realIntRootRadii_rootRadii : center: %ld + 0i,  prec: %ld \n", 
                                                                                  centerRe,       nprec);
        }
        
        realIntRootRadii_getApproximation_real( pApprox, cache, nprec, meta );
        if (centerRe != 0){
            realIntRootRadii_taylor_shift_inplace_real( pApprox, centerRe, nprec, meta);
            if (metadatas_haveToCount(meta)) {
                if (nprec==prec)
                    (metadatas_countref(meta))[0].RR_nbTaylors += 1;
                else
                    (metadatas_countref(meta))[0].RR_nbTaylorsRepeted += 1;
            }
        }
        
        int res = realIntRootRadii_GraeffeAndCH_real ( convexHull, &lenCh, &nnprec, absCoeffs, pApprox, N, nprec, meta );
        if (res < N) { /* double precision */
            nprec = 2*nprec;
            lenCh = 0;
            
        }
        if (metadatas_haveToCount(meta)) {
            if (lenCh == 0)
                (metadatas_countref(meta))[0].RR_nbGraeffeRepeted += res;
            else
                (metadatas_countref(meta))[0].RR_nbGraeffe += N;
        }
    }
    
//     if (metadatas_getVerbo(meta)>=level){
//         printf("# Convex hull: %ld vertices: ", lenCh );
//         for (slong ind = 0; ind < lenCh; ind++)
//             printf("%ld, ", convexHull[ind]);
//         printf("\n");
//     }
    
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
        
        if ( realApp_contains_zero( (absCoeffs->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
            realApp_zero( compAnn_radSupref(cur) );
        } else {
            realApp_div( compAnn_radInfref(cur), (absCoeffs->coeffs) + right, (absCoeffs->coeffs) + left, prec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, prec );
            realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), prec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, prec );
            realApp_mul( compAnn_radSupref(cur), compAnn_radInfref(cur), relError, prec );
            realApp_mul( compAnn_radInfref(cur), compAnn_radInfref(cur), relErrorInv, prec );
        }
        
        left = convexHull[ind];
        compAnn_list_push(annulii, cur);
    }
    
    realApp_clear(relError);
    realApp_clear(relErrorInv);
    
    realApp_poly_clear(pApprox);
    realApp_poly_clear(absCoeffs);
    ccluster_free(convexHull);
    
    return nprec;
}

slong realIntRootRadii_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
                                             slong centerIm,
                                             cacheApp_t cache,        /* polynomial */
                                             slong prec,
                                             metadatas_t meta ){
    
    int level = 4;
    slong degree = cacheApp_getDegree(cache);
    
    int N = metadatas_getNbGIt(meta);
    ulong pow = 0x1<<N;
    realApp_t relError, relErrorInv;
    realApp_init(relError);
    realApp_init(relErrorInv);
    realApp_set_si(relError, 2*degree);
    realApp_root_ui(relError, relError, pow, prec);
    realApp_inv(relErrorInv, relError, prec);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    slong nnprec = prec;
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,degree+1);
    realApp_poly_t pSquares;
    realApp_poly_init2(pSquares,degree+1);
    
    while ( lenCh == 0 ) {
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---realIntRootRadii.c; realIntRootRadii_rootRadii : center: 0 + %ldi,  prec: %ld \n", 
                                                                                     centerIm,       nprec);
        }
        
        realIntRootRadii_getApproximation_comp( pApprox, cache, nprec, meta );
        if (centerIm != 0) {
            realIntRootRadii_taylor_shift_inplace_comp( pApprox, 0, centerIm, nprec, meta);
            if (metadatas_haveToCount(meta)) {
                if (nprec==prec)
                    (metadatas_countref(meta))[0].RR_nbTaylors += 1;
                else
                    (metadatas_countref(meta))[0].RR_nbTaylorsRepeted += 1;
            }
        }
        
        int res = realIntRootRadii_GraeffeAndCH_comp ( convexHull, &lenCh, &nnprec, pSquares, pApprox, N, nprec, meta );
        if (res < N) { /* double precision */
            nprec = 2*nprec;
            lenCh = 0;
            
        }
        if (metadatas_haveToCount(meta)) {
            if (lenCh == 0)
                (metadatas_countref(meta))[0].RR_nbGraeffeRepeted += res;
            else
                (metadatas_countref(meta))[0].RR_nbGraeffe += N;
        }
    }
    
    for(slong i = 0; i <= degree; i++)
        realApp_sqrt( (pSquares->coeffs)+i, (pSquares->coeffs)+i, nprec);
    
//     if (metadatas_getVerbo(meta)>=level){
//         printf("# Convex hull: %ld vertices: ", lenCh );
//         for (slong ind = 0; ind < lenCh; ind++)
//             printf("%ld, ", convexHull[ind]);
//         printf("\n");
//     }
    
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
        compAnn_centerReref(cur) = 0;
        compAnn_centerImref(cur) = centerIm;
        
        if ( realApp_contains_zero( (pSquares->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
            realApp_zero( compAnn_radSupref(cur) );
        } else {
            realApp_div( compAnn_radInfref(cur), (pSquares->coeffs) + right, (pSquares->coeffs) + left, nprec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, prec );
            realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), prec );
            realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, nprec );
            realApp_mul( compAnn_radSupref(cur), compAnn_radInfref(cur), relError, prec );
            realApp_mul( compAnn_radInfref(cur), compAnn_radInfref(cur), relErrorInv, prec );
        }
        left = convexHull[ind];
        compAnn_list_push(annulii, cur);
    }
    
    realApp_clear(relError);
    realApp_clear(relErrorInv);
    
    compApp_poly_clear(pApprox);
    realApp_poly_clear(pSquares);
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
    realApp_poly_set(pApprox, cacheApp_getApproximation_real ( cache, prec ));
    
    realApp_t centerLeft, centerLeftVal, centerRight, centerRightVal;
    realApp_init(centerLeft);
    realApp_init(centerLeftVal);
    realApp_init(centerRight);
    realApp_init(centerRightVal);
    
    realApp_t interval;
    realApp_init(interval);
    
    it = compAnn_list_begin(annulii);
    while (it!=compAnn_list_end() ) {
        cur = compAnn_list_elmt( it );
        /* check if it is 0 */
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
    realApp_clear(interval);
    
    realApp_clear(centerLeft);
    realApp_clear(centerLeftVal);
    realApp_clear(centerRight);
    realApp_clear(centerRightVal);
    
}

/* DEPRECATED VERSION */
/* returns 0 if after some iterations, the relative accuracy is less than 1 */
int realIntRootRadiiOLD_Ngraeffe_iterations_inplace_real( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
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
int realIntRootRadiiOLD_Ngraeffe_iterations_inplace_comp( compApp_poly_t res, int N, slong prec, metadatas_t meta){
    
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
slong realIntRootRadiiOLD_convexHull( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec ){
    
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

slong realIntRootRadiiOLD_rootRadii( compAnn_list_t annulii,  /* list of annulii */
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
    
//     printf("#realIntRootRadiiOLD.c; realIntRootRadii_rootRadiiOLD : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
        printf("#---realIntRootRadiiOLD.c; realIntRootRadii_rootRadii : prec: %ld \n", nprec);
        
        realIntRootRadii_getApproximation_real( pApprox, cache, nprec, meta );
        if (centerRe != 0)
            realIntRootRadii_taylor_shift_inplace_real( pApprox, centerRe, nprec, meta);
        int enoughRelacc = realIntRootRadiiOLD_Ngraeffe_iterations_inplace_real( pApprox, N, nprec, meta);
//         printf("---realIntRootRadiiOLD_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute abs of coeffs */
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++) {
                realApp_abs( (pApprox->coeffs)+i, (pApprox->coeffs)+i );
            }
            /* compute convex hull */
            lenCh = realIntRootRadiiOLD_convexHull( convexHull, (pApprox->coeffs), degree+1, nprec );
        }
//         printf("---realIntRootRadiiOLD_rootRadii: lenCh: %ld\n", lenCh);
        
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

slong realIntRootRadiiOLD_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
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
//     printf("#realIntRootRadiiOLD.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,degree+1);
    realApp_poly_t pSquares;
    realApp_poly_init2(pSquares,degree+1);
    
    int nbRep = 0;
    
    while ( lenCh == 0 ) {
        
//         printf("#---realIntRootRadiiOLD.c; realIntRootRadii_rootRadii_imagCenter : prec: %ld \n", nprec);
        
        realIntRootRadii_getApproximation_comp( pApprox, cache, nprec, meta );
        if (centerIm != 0)
            realIntRootRadii_taylor_shift_inplace_comp( pApprox, 0, centerIm, nprec, meta);
        int enoughRelacc = realIntRootRadiiOLD_Ngraeffe_iterations_inplace_comp( pApprox, N, nprec, meta);
//         printf("---realIntRootRadiiOLD_rootRadii: enoughRelacc: %d, prec: %ld\n", enoughRelacc, nprec);
        /* compute sum of squares of real and imaginary parts of coeffs */
        if (enoughRelacc==1) {
            for(slong i = 0; i <= degree; i++){
                realApp_sqr( compApp_realref((pApprox->coeffs)+i), compApp_realref((pApprox->coeffs)+i), nprec );
                realApp_sqr( compApp_imagref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec );
                realApp_add( (pSquares->coeffs)+i, compApp_realref((pApprox->coeffs)+i), compApp_imagref((pApprox->coeffs)+i), nprec);
            }
            /* compute convex hull */
            lenCh = realIntRootRadiiOLD_convexHull( convexHull, (pSquares->coeffs), degree+1, nprec );
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
