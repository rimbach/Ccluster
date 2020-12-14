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

#include "realApp_rootRadii.h"

void realApp_rootRadii_getApproximation_real( realApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta){
        clock_t start = clock();

        realApp_poly_set(res, cacheApp_getApproximation_real ( cache, prec ));

        if (metadatas_haveToCount(meta))
            metadatas_add_time_Approxi(meta, (double) (clock() - start) );
        
}

void realApp_rootRadii_Ngraeffe_iterations_inplace_real( realApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        for(int i = 0; i < N; i++)
            realApp_poly_oneGraeffeIteration_in_place( res, prec );
        
        if (metadatas_haveToCount(meta)) {
            clock_t end = clock();
            metadatas_add_time_Graeffe(meta, (double) (end - start) );
            metadatas_add_time_RRGraef(meta, (double) (end - start) );
        }
        
}

/* assume i<j<k */
/* assume logAbsPi, logAbsPj, logAbsPk are exact and >=-prec */
/* returns 1 if [j, logAbsPj] lies strictly below or on the line */
/*                            passing trough [i,logAbsPi] and [k,logAbsPk]*/
/* returns 0 if it is strictly above */
int realApp_rootRadii_liesBelow( slong i, const realApp_t logAbsPi,
                                 slong j, const realApp_t logAbsPj,
                                 slong k, const realApp_t logAbsPk, slong prec ){
    
    if ( (realApp_is_finite(logAbsPj)==0)&&(realApp_is_negative(logAbsPj)==1) ) { /* in which case it is -inf */
        /* no point (u,-inf) lies above any line in the plane */
        return 1;
    } else if ( (realApp_is_finite(logAbsPk)==0)&&(realApp_is_negative(logAbsPk)==1) ) { /* in which case it is -inf */
        /* no point (u,-inf) lies above any line in the plane */
        return 0;
    }
    realApp_t slopeij, slopejk;
    realApp_init(slopeij);
    realApp_init(slopejk);
    int res = -1;
    slong nprec = CCLUSTER_DEFAULT_PREC;
    while (res==-1){
        realApp_sub(slopeij, logAbsPj, logAbsPi, nprec);
        realApp_mul_si( slopeij, slopeij, k-j, nprec);
        realApp_sub(slopejk, logAbsPk, logAbsPj, nprec);
        realApp_mul_si( slopejk, slopejk, j-i, nprec);
        if (realApp_gt( slopeij, slopejk ))
            res = 0;
        else if (realApp_lt( slopeij, slopejk ))
            res = 1;
        else if ( realApp_is_exact(slopeij) &&
                 realApp_is_exact(slopejk) &&
                 realApp_eq(slopeij, slopejk) )
            res = 1;
        else{
            nprec = 2*nprec;
            res = -1;
        }
    }
    
    return res;
    
    realApp_clear(slopeij);
    realApp_clear(slopejk);
}

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* assume logAbsCoeffs are len exact numbers >=-prec */
/* returns the length of the convex hull */
slong realApp_rootRadii_convexHull( slong * convexHull, const realApp_ptr logAbsCoeffs, slong len, slong prec ){
    
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
            liesBelow = realApp_rootRadii_liesBelow( convexHull[res-2], logAbsCoeffs + convexHull[res-2], 
                                                 convexHull[res-1], logAbsCoeffs + convexHull[res-1],
                                                 i,     logAbsCoeffs + i, 
                                                 prec);
            if (liesBelow == 1)
                res--;
//             if (liesBelow == -1) {
//                 res=0;
//                 return res;
//             }
        }
//         printf("#realIntRootRadii.c , realIntRootRadii_convexHull: res: %ld\n", res);
        convexHull[res] = i;
        res++;
    }
    
    return res;
}

/* returns the precision used to carry out root radii */
slong realApp_rootRadii_fromZero( compAnn_list_t annulii,  /* list of annulii */
                                  cacheApp_t cache,        /* polynomial */
                                  const realRat_t delta,
                                  slong prec,
                                  metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    
    realRat_t oneplusdelta, oneplusdeltainv;
    
    realRat_init(oneplusdelta);
    realRat_init(oneplusdeltainv);
    
    realRat_add_si(oneplusdelta, delta, 1);
    realRat_inv( oneplusdeltainv, oneplusdelta );
    
    double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
    log2_1pdelta = log2_1pdelta / log(2);
    int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    N = N+1; /* for the approximated algorithm */
//     N=0;
    
//     if (metadatas_getVerbo(meta)>=3)
        printf("#realApp_rootRadii.c; realApp_rootRadii_fromZero : number of Graeffe iterations: %d \n", N);
    slong recprec  = (slong) ceil( log2(2*degree) ) +1;
    slong nrecprec = recprec;
    //     if (metadatas_getVerbo(meta)>=3)
//         printf("#realApp_rootRadii.c; realApp_rootRadii_fromZero : required precision on coeffs for task AS': %ld \n", recprec2);
    realApp_t threshold, logthreshold;
    realApp_init(threshold);
    realApp_init(logthreshold);
    realApp_set_si(logthreshold, -(nrecprec-1));
    realApp_one(threshold);
    realApp_mul_2exp_si(threshold, threshold, -(nrecprec-1));
    
    int checkPrec = 0;
    slong * convexHull = (slong *) ccluster_malloc ( (degree+1)*sizeof(slong) );
    
    slong nprec = prec;
    realApp_poly_t pApprox;
    realApp_poly_init2(pApprox,degree+1);
    
    realApp_ptr    absCoeffs = (realApp_ptr) ccluster_malloc ( (degree+1)*sizeof(realApp) );
    realApp_ptr logAbsCoeffs = (realApp_ptr) ccluster_malloc ( (degree+1)*sizeof(realApp) );
    for (slong index = 0; index <= degree; index++) {
        realApp_init(    absCoeffs + index );
        realApp_init( logAbsCoeffs + index );
    }
    
    slong lenCh;
    realApp_t firstWidth;
    realApp_init( firstWidth );
    
    /* compute the logs of the abs of the coeffs at precision recprec */
    while ( checkPrec==0 ) {
        printf("actual working precision: %ld\n", nprec);
        checkPrec = 1;
        realApp_rootRadii_getApproximation_real( pApprox, cache, nprec, meta );
        realApp_rootRadii_Ngraeffe_iterations_inplace_real( pApprox, N, nprec, meta);
        for(slong index = 0; (index <= degree)&&(checkPrec); index++){
            realApp_abs( (pApprox->coeffs)+index, (pApprox->coeffs)+index );
            realApp_set( absCoeffs+index, (pApprox->coeffs)+index );
            if ( realApp_contains_zero (absCoeffs+index) ) {
                
//                 /* check if it is exactly zero */
//                 if (realApp_is_exact(absCoeffs+index)){
//                     printf("# abs of %ld-th coeff is exactly zero: ", index); realApp_printd(absCoeffs+index,10); printf("\n");
//                     /* do not check precision */
//                     realApp_neg_inf( logAbsCoeffs+index );
//                     printf("# set to infinity: "); realApp_printd(logAbsCoeffs+index,10); printf("\n");
// //                     printf("# test is_finite: %d, test is_negative: %d\n", realApp_is_finite( logAbsCoeffs+index ), realApp_is_negative(logAbsCoeffs+index) );
//                 } else {
                    printf("# abs of %ld-th coeff contains zero: ", index); realApp_printd(absCoeffs+index,10); printf("\n");
                    checkPrec = (checkPrec) && (realApp_check_absolute_accuracy(absCoeffs+index, nrecprec));
                    realApp_set(    absCoeffs+index, threshold );
                    realApp_set( logAbsCoeffs+index, logthreshold );
                    printf("# set to minus threshold: "); realApp_printd(logAbsCoeffs+index,10); printf("\n");
//                 }
            } else {
                realApp_log_base_ui(logAbsCoeffs+index, absCoeffs+index, 2, nprec);
                checkPrec = (checkPrec) && (realApp_check_absolute_accuracy(logAbsCoeffs+index, nrecprec));
                realApp_get_mid_realApp(absCoeffs+index, absCoeffs+index);
                realApp_get_mid_realApp(logAbsCoeffs+index, logAbsCoeffs+index);
                if (realApp_gt(logAbsCoeffs+index, logthreshold)!=1) {
                    printf("# log of abs of %ld-th coeff less than %ld: ", index, -(nrecprec-1)); realApp_printd(absCoeffs+index,10); printf("\n");
                    realApp_set(   absCoeffs+index,    threshold );
                    realApp_set(logAbsCoeffs+index, logthreshold );
                }
            }
        }
        if (checkPrec==0) {
            nprec = 2*nprec;
        } 
//         else {
//             realApp_one(threshold);
//             realApp_mul_2exp_si(threshold, threshold, -prec);
//     
//             lenCh = realApp_rootRadii_convexHull( convexHull, logAbsCoeffs, degree+1, (recprec-1) );
//             if ( realApp_contains_zero( (pApprox->coeffs) + 0 ) ) {
//                 realApp_div( firstWidth, absCoeffs + convexHull[0], absCoeffs + convexHull[1], nprec );
//                 realApp_root_ui( firstWidth, firstWidth, convexHull[1] - convexHull[0], nprec );
//                 ulong pow = 0x1<<N;
//                 printf("Width of first annulii: "); realApp_printd(firstWidth,10); printf("\n");
//                 realApp_root_ui( firstWidth, firstWidth, pow, nprec );
//                 printf("Width of first annulii: "); realApp_printd(firstWidth,10); printf("\n");
//                 if (realApp_gt(firstWidth, threshold)){
//                     printf("not enough precision\n");
//                     checkPrec = 0;
//                 } else {
//                     checkPrec = 1;
//                 }
//             } else {
//                 checkPrec = 1;
//             }
//         }
//         if (checkPrec==0) {
//             if (nrecprec<prec)
//                 nrecprec = 2*prec;
//             else 
//                 nrecprec = 2*nrecprec;
//             realApp_set_si(logthreshold, -(nrecprec-1));
//             realApp_one(threshold);
//             realApp_mul_2exp_si(threshold, threshold, -(nrecprec-1));
//             nprec = 2*nrecprec;
//         }
    }
//     realApp_clear( firstWidth );
//     if (metadatas_getVerbo(meta)>=3)
        printf("#realApp_rootRadii.c; realApp_rootRadii_fromZero : required precision: %ld \n", nprec);
    lenCh = realApp_rootRadii_convexHull( convexHull, logAbsCoeffs, degree+1, (recprec-1) );
    
    
//     if (metadatas_getVerbo(meta)>=3){
        printf("# Convex hull: %ld vertices: ", lenCh );
        for (slong ind = 0; ind < lenCh; ind++)
            printf("%ld, ", convexHull[ind]);
        printf("\n");
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
        compAnn_centerReref(cur) = 0;
        compAnn_centerImref(cur) = 0;
        
        realApp_div( compAnn_radInfref(cur), absCoeffs + left, absCoeffs + right, nprec );
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), shift, nprec );
        ulong pow = 0x1<<N;
        realApp_root_ui( compAnn_radInfref(cur), compAnn_radInfref(cur), pow, nprec );
        realApp_mul_realRat( compAnn_radSupref(cur), compAnn_radInfref(cur), oneplusdelta, nprec );
        realApp_mul_realRat_in_place( compAnn_radInfref(cur), oneplusdeltainv, nprec );
        
        if ( realApp_contains_zero( (pApprox->coeffs) + left ) ) {
            realApp_zero( compAnn_radInfref(cur) );
//             if (realApp_is_exact( (pApprox->coeffs) + left ) ){
//                 realApp_zero( compAnn_radSupref(cur) );
//             }
        }
        if ( realApp_contains_zero( (pApprox->coeffs) + right ) ) {
            realApp_pos_inf( compAnn_radInfref(cur) );
        }
        left = convexHull[ind];
//         compAnn_printd(cur, 10); printf("\n");
        compAnn_list_push(annulii, cur);
    }
    
    for (slong index = 0; index <= degree; index++){
        realApp_clear(    absCoeffs + index );
        realApp_clear( logAbsCoeffs + index );
    }
    ccluster_free(   absCoeffs);
    ccluster_free(logAbsCoeffs);
    
    realApp_clear(threshold);
    
    realApp_poly_clear(pApprox);
    realRat_clear(oneplusdelta);
    realRat_clear(oneplusdeltainv);
    ccluster_free(convexHull);
    
    return nprec;
    
}

void realApp_rootRadii_connectedComponents( compAnn_list_t annulii, slong prec ){
    
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
