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

void realIntRootRadii_getApproximation( realApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta){
        clock_t start = clock();

        realApp_poly_set(res, cacheApp_getApproximation_real ( cache, prec ));

        if (metadatas_haveToCount(meta))
            metadatas_add_time_Approxi(meta, (double) (clock() - start) );
        
}

void realIntRootRadii_graeffe_iterations_inplace( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        for(int i = 0; i < N; i++)
            realApp_poly_oneGraeffeIteration_in_place( res, prec );
        
        if (metadatas_haveToCount(meta))
            metadatas_add_time_Graeffe(meta, (double) (clock() - start) );
}

void realIntRootRadii_taylor_shift_inplace( realApp_poly_t res, slong centerRe, slong prec){
    
//         clock_t start = clock();
        realApp_poly_taylorShift_in_place_slong( res, 
                                           centerRe, 
                                           prec );

//         if (metadatas_haveToCount(meta))
//             metadatas_add_time_Taylors(meta, (double) (clock() - start) );
}

/* assume i<j<k */
/* assume absPi=|pi|, absPj=|pj|, absPk=|pk| are approximations of integers */
/* decide if [j,log|pj|] lies below of on the line passing trough [i,log|pi|] and [k,log|pk|]*/
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
    ulong kmi = k-i;
    ulong jmi = j-i;
    realApp_t leftSide, rightSide, temp;
    realApp_init(leftSide);
    realApp_init(rightSide);
    realApp_init(temp);
    int res=-1;
    
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

slong realIntRootRadii_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                               cacheApp_t cache,        /* polynomial */
                               const realRat_t delta){
    
    slong degree = cacheApp_getDegree(cache);
    
    realRat_t oneplusdelta, oneplusdeltainv;
    
    realRat_init(oneplusdelta);
    realRat_init(oneplusdeltainv);
    
    realRat_add_si(oneplusdelta, delta, 1);
    realRat_inv( oneplusdeltainv, oneplusdelta );
    
    double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
    log2_1pdelta = log2_1pdelta / log(2);
    int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    
    printf("#realIntRootRadii.c; realIntRootRadii_rootRadii : number of Graeffe iterations: %d \n", N);
    
    slong lenCh = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    if (centerRe != 0)
        realIntRootRadii_taylor_shift_inplace( pApprox, centerRe, prec);
    
    while ( lenCh == 0 ) {
        
        //     realIntRootRadii_getApproximation( pApprox, cache, prec, meta );
        realApp_poly_set(pApprox, cacheApp_getApproximation_real ( cache, prec ));
        if (centerRe != 0)
            realIntRootRadii_taylor_shift_inplace( pApprox, centerRe, prec);
        //  realIntRootRadii_graeffe_iterations_inplace( pApprox, N, prec, meta);
        for(int i = 0; i < N; i++)
            realApp_poly_oneGraeffeIteration_in_place( pApprox, prec );
        /* compute abs of coeffs */
        for(slong i = 0; i <= degree; i++)
            realApp_abs( (pApprox->coeffs)+i, (pApprox->coeffs)+i );
        /* compute convex hull */
        lenCh = realIntRootRadii_convexHull( convexHull, (pApprox->coeffs), degree+1, prec );
        
        if (lenCh==0) /* double precision */
            prec = 2*prec;
    }
    
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
        compAnn_centerref(cur) = centerRe;
        realApp_div( compAnn_radInfref(cur), (pApprox->coeffs) + right, (pApprox->coeffs) + left, prec );
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, prec );
        realApp_inv( compAnn_radInfref(cur), compAnn_radInfref(cur), prec );
        ulong pow = 0x1<<N;
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, prec );
        realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, prec );
        realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, prec );
        
        left = convexHull[ind];
//         compAnn_printd(cur, 10); printf("\n");
        compAnn_list_push(annulii, cur);
    }
    
    realApp_poly_clear(pApprox);
    realRat_clear(oneplusdelta);
    realRat_clear(oneplusdeltainv);
    ccluster_free(convexHull);
    
    return prec;
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
        } else { /* try interval evaluations */
            realApp_get_mid_realApp( centerLeft, compAnn_radInfref(cur) );
            realApp_get_mid_realApp( centerRight, compAnn_radSupref(cur) );
            realApp_union(interval, centerLeft, centerRight, prec);
            realApp_poly_evaluate(interval, pApprox, interval, prec);
            if (realApp_contains_zero(interval)==0)
                compAnn_rrInPoref(cur) = 0;
            
            realApp_neg( centerLeft, centerLeft );
            realApp_neg( centerRight, centerRight ); /* no need to swap */
            realApp_union(interval, centerLeft, centerRight, prec);
            realApp_poly_evaluate(interval, pApprox, interval, prec);
            if (realApp_contains_zero(interval)==0)
                compAnn_rrInNeref(cur) = 0;
        }
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
            compBox_init_annuli(bstemp);
            compBox_actualize_anulii_risolate( bstemp, compBox_annuliref(btemp) );
            
            int nbSol=-1;
            if ( compAnn_list_get_size(compBox_annuliref(bstemp)) == 0 )
                nbSol = 0;
            if ( compAnn_list_get_size(compBox_annuliref(bstemp)) == 1 ) {
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
                
                compAnn_ptr ann = compAnn_list_first( compBox_annuliref(bstemp) ); /* there is only one element in the list */
                
                if ( ( (realApp_is_positive(left)==1)  && (compAnn_rrInPoref(ann) == 0) ) ||
                     ( (realApp_is_negative(right)==1)  && (compAnn_rrInNeref(ann) == 0) ) )
                    nbSol=0;
                
                realApp_clear(center);
                realApp_clear(left);
                realApp_clear(right);
                realApp_clear(rad);
            }
            
            if (nbSol==0) { /* delete bstemp */
                compBox_clear_annuli(bstemp);
                compBox_clear(bstemp);
                ccluster_free(bstemp);
            } else {
                connCmp_union_compBox( dest, bstemp);
            }
        }
        
        /* delete btemp */
        compBox_clear_annuli(btemp);
        compBox_clear(btemp);
        ccluster_free(btemp);
        
        
    }
    
    compBox_list_clear(subBoxes);
}
