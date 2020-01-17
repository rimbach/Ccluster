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
#include <time.h>

void tstar_getApproximation( compApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta){
        clock_t start = clock();

// #ifdef CCLUSTER_HAVE_PTHREAD
//         if (metadatas_useNBThreads(meta) >1)
//             cacheApp_lock(cache);
// #endif
        compApp_poly_set(res, cacheApp_getApproximation ( cache, prec ));
// #ifdef CCLUSTER_HAVE_PTHREAD
//         if (metadatas_useNBThreads(meta) >1)
//             cacheApp_unlock(cache);
// #endif
        if (metadatas_haveToCount(meta))
            metadatas_add_time_Approxi(meta, (double) (clock() - start) );
        
}

void tstar_taylor_shift_inplace( compApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        compApp_poly_taylorShift_in_place( res, 
                                           compDsk_centerref(d), 
                                           compDsk_radiusref(d), 
                                           prec );

        if (metadatas_haveToCount(meta))
            metadatas_add_time_Taylors(meta, (double) (clock() - start) );
}

void tstar_graeffe_iterations_inplace( compApp_poly_t res, int N, slong prec, metadatas_t meta){
    
        clock_t start = clock();
        for(int i = 0; i < N; i++)
            compApp_poly_oneGraeffeIteration_in_place( res, prec );
        
        if (metadatas_haveToCount(meta))
            metadatas_add_time_Graeffe(meta, (double) (clock() - start) );
}

void tstar_graeffe_iterations_abs_two_first_coeffs( realApp_t coeff0, realApp_t coeff1, const compApp_poly_t pApprox, int N, slong prec, metadatas_t meta){
    compApp_poly_t p1, p2;
    compApp_poly_init2(p1, compApp_poly_degree(pApprox)+1);
    compApp_poly_init2(p2, compApp_poly_degree(pApprox)+1);
    compApp_poly_set(p1, pApprox);
    slong bound = 0x1 << N; /* assume 2^N fits in an slong: N<32 for 32bits machine... */
    for( int i = 0; i < N; i++) {
        bound = bound >> 1;
//         printf("bound: %d\n", (int) bound);
        compApp_poly_oneGraeffeIteration_first_n_coeff( p2, p1, CCLUSTER_MIN(compApp_poly_degree(pApprox), bound), compApp_poly_degree(pApprox)+1, prec);
        compApp_poly_swap(p1,p2);
    }
    
    compApp_abs( coeff0, compApp_poly_getCoeff(p1, 0), prec);
    compApp_abs( coeff1, compApp_poly_getCoeff(p1, 1), prec);
    
    compApp_poly_clear(p1);
    compApp_poly_clear(p2);
}



tstar_res tstar_interface( cacheApp_t cache,
                           const compDsk_t d,
                           int max_nb_sols,    /*the maximum number of sols in the disk          */
                           int discard,        /*a flag saying if it is a discarding test or not */
                           int inNewton,      /*a flag saying if it is for newton refinement     */
                           slong prec,         /*the "default" arithmetic precision              */
                           int depth,          /*the depth for counter                           */
                           metadatas_t meta){
    slong nprec = CCLUSTER_DEFAULT_PREC;
    
    if (metadatas_usePredictPrec(meta))
        nprec = prec;
    
    if (metadatas_useTstarOptim(meta)) {
        if (discard&&CCLUSTER_V2(meta)){
            return tstar_optimized( cache, d, 0, discard, inNewton, nprec, depth, meta);
        }
        else {
            return tstar_optimized( cache, d, max_nb_sols, discard, inNewton, nprec, depth, meta);
        }
    }
    if (discard)
        return tstar_asInPaper( cache, d, 0, discard, inNewton, nprec, depth, meta);
    
    return tstar_asInPaper( cache, d, max_nb_sols, discard, inNewton, nprec, depth, meta);
    
}

tstar_res tstar_asInPaper( cacheApp_t cache,
                           const compDsk_t d,
                           int max_nb_sols,    /*the maximum number of sols in the disk,         */
                           int discard,        /*a flag saying if it is a discarding test or not */
                           int inNewton,      /*a flag saying if it is for newton refinement     */
                           slong prec,         /*the "default" arithmetic precision              */
                           int depth,          /*the depth for counter                           */
                           metadatas_t meta){
    
    clock_t start = clock();
    
    tstar_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    int restemp = 0;
    int nbTaylorsRepeted = 0;
    int nbGraeffeRepeted = 0;
    int N = 0;
    slong deg = cacheApp_getDegree(cache);
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox, deg+1);
    realApp_t sum;
    realApp_init(sum);
    N = (int) 4+ceil(log2(1+log2(deg)));
    
    tstar_getApproximation( pApprox, cache, res.appPrec, meta);
    tstar_taylor_shift_inplace( pApprox, d, res.appPrec, meta);
    tstar_graeffe_iterations_inplace( pApprox, N, res.appPrec, meta);
    compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
    
    while( (res.nbOfSol < max_nb_sols)&&(restemp==0) ){
        res.nbOfSol += 1;
        restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
        
        while( restemp == -2 ){
            res.appPrec *=2;
            tstar_getApproximation( pApprox, cache, res.appPrec, meta);
            tstar_taylor_shift_inplace( pApprox, d, res.appPrec, meta);
            tstar_graeffe_iterations_inplace( pApprox, N, res.appPrec, meta);
            compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
            restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
            nbTaylorsRepeted +=1;
            nbGraeffeRepeted +=N;
        }
        if (restemp == -1)
            restemp = 0;
        
    }
        
    compApp_poly_clear(pApprox);
    realApp_clear(sum);
    
    if (restemp==0) res.nbOfSol = -1;

    if (metadatas_haveToCount(meta))
        metadatas_add_Test( meta, depth, (restemp==1), discard, inNewton, 1, nbTaylorsRepeted, N, 
                            nbGraeffeRepeted, (int) res.appPrec, (double) (clock() - start) );
    return res;
}

tstar_res tstar_optimized( cacheApp_t cache,
                           const compDsk_t d,
                           int max_nb_sols,   /*the maximum number of sols in the disk          */
                           int discard,       /*a flag saying if it is a discarding test or not */
                           int inNewton,      /*a flag saying if it is for newton refinement     */
                           slong prec,        /*the "default" arithmetic precision              */
                           int depth,         /*the depth for counter                           */
                           metadatas_t meta){
    
    clock_t start = clock();
    
    tstar_res res;
    res.nbOfSol = -1;
    res.appPrec = prec;
    int restemp = 0;
    int TS_has_been_computed=0;
    int nbTaylorsRepeted = 0;
    int nbGraeffeRepeted = 0;
    int nbGraeffe = 0;
    int iteration = 0;
    int N = 0;
    slong deg = cacheApp_getDegree(cache);
    compApp_poly_t pApprox;
    compApp_poly_init2(pApprox,deg+1);
    realApp_t sum;
    realApp_init(sum);
    
    realApp_t coeff0, coeff1, coeffn; /* for anticipate */
    int anticipate_already_applied = 0;
    N = (int) 4+ceil(log2(1+log2(deg)));
    
#ifdef CCLUSTER_EXPERIMENTAL
    if ((discard)&&( CCLUSTER_INC_TEST(meta) )) {
        restemp = tstar_inclusion_test_wn( cache, d, res.appPrec, depth, meta);
        if (restemp==1) {
            compApp_poly_clear(pApprox);
            realApp_clear(sum);
            res.nbOfSol = -1;
            metadatas_add_Test( meta, depth, 0, discard, inNewton, 0, 0, 0, 0, (double) (clock() - start));
            return res;
        }
        restemp = 0;
    }
    
    if ( discard && (CCLUSTER_EXP_NUM_T0(meta)||CCLUSTER_EXP_NUM_T1(meta)) && (max_nb_sols<=1) ) {
        
        restemp = D0N1_test ( cache, pApprox, d, depth, res.appPrec, meta );
        
        if (restemp==0) {
            compApp_poly_clear(pApprox);
            realApp_clear(sum);
            res.nbOfSol = 0;
            metadatas_add_Test( meta, depth, 1, discard, inNewton, 0, 0, 0, 0, (double) (clock() - start));
//             printf("&&depth: %d, success of D0 discarding predicate\n", depth);
            return res;
        }
        
        if (restemp==1) {
            compApp_poly_clear(pApprox);
            realApp_clear(sum);
            res.nbOfSol = 1;
            metadatas_add_Test( meta, depth, 1, discard, inNewton, 0, 0, 0, 0, (double) (clock() - start));
//             printf("&&depth: %d, success of N1 non-discarding discarding predicate\n", depth);
            return res;
        }
    
        TS_has_been_computed = 1;
        restemp = 0;
    }
#endif

    if (TS_has_been_computed==0) {
        tstar_getApproximation( pApprox, cache, res.appPrec, meta);
        tstar_taylor_shift_inplace( pApprox, d, res.appPrec, meta);
    }
    
    
    if ( (discard)&&(metadatas_useAnticipate(meta)) ){
        realApp_init(coeff0);
        realApp_init(coeffn);
        
        compApp_abs( coeff0, compApp_poly_getCoeff(pApprox, 0), res.appPrec);
        compApp_abs( coeffn, compApp_poly_getCoeff(pApprox, compApp_poly_degree(pApprox)), res.appPrec);
        restemp = realApp_soft_compare( coeff0, coeffn, res.appPrec );
        while( restemp == -2 ){
            res.appPrec *=2;
            tstar_getApproximation( pApprox, cache, res.appPrec, meta);
            tstar_taylor_shift_inplace( pApprox, d, res.appPrec, meta);
            nbTaylorsRepeted +=1;
            compApp_abs( coeff0, compApp_poly_getCoeff(pApprox, 0), res.appPrec);
            compApp_abs( coeffn, compApp_poly_getCoeff(pApprox, compApp_poly_degree(pApprox)), res.appPrec);
            restemp = realApp_soft_compare( coeff0, coeffn, res.appPrec );
        }
            if (restemp==0)
                restemp = -1;
            else
                restemp = 0;
        
        realApp_clear(coeff0);
        realApp_clear(coeffn);
    }
    
    while( (iteration <= N)&&(restemp==0) ){
        
        if (iteration >= 1) {
            tstar_graeffe_iterations_inplace( pApprox, 1, res.appPrec, meta);
            nbGraeffe +=1;
        }
        
        compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
        
        res.nbOfSol = -1;
        while( (res.nbOfSol < max_nb_sols)&&(restemp==0)&&(res.nbOfSol<deg) ){
            res.nbOfSol += 1;
            
            restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
        
            while( restemp == -2 ){
                res.appPrec *=2;
                tstar_getApproximation( pApprox, cache, res.appPrec, meta);
                tstar_taylor_shift_inplace( pApprox, d, res.appPrec, meta);
                tstar_graeffe_iterations_inplace( pApprox, iteration, res.appPrec, meta);
                compApp_poly_sum_abs_coeffs( sum, pApprox, res.appPrec );
                restemp = compApp_poly_TkGtilda_with_sum( pApprox, sum, res.nbOfSol, res.appPrec);
                nbTaylorsRepeted +=1;
                nbGraeffeRepeted +=(iteration);
            }
        }
        
        if (restemp == -1)
            restemp = 0;
        
        if ( (discard) && (metadatas_useAnticipate(meta)) && (anticipate_already_applied==0) && (restemp == 0) ) {
           
            int test_anticipate = ((0x1 << (N-iteration)) <= (compApp_poly_degree(pApprox)/4));
            if (test_anticipate) {
                
                clock_t start2 = clock();
                
                anticipate_already_applied = 1;
                realApp_init(coeff0);
                realApp_init(coeff1);
                realApp_init(coeffn);
                
                tstar_graeffe_iterations_abs_two_first_coeffs( coeff0, coeff1, pApprox, N-iteration, res.appPrec, meta);
                compApp_abs( coeffn, compApp_poly_getCoeff(pApprox, compApp_poly_degree(pApprox)), res.appPrec);
                realApp_pow_ui( coeffn, coeffn, (ulong)(0x1 << (N-iteration)), res.appPrec);
                realApp_add(coeffn, coeffn, coeff1, res.appPrec);
                restemp = realApp_soft_compare( coeff0, coeffn, res.appPrec );
                
                if (restemp==0)
                    restemp = -1;
                else
                    restemp = 0;
        
                realApp_clear(coeff0);
                realApp_clear(coeff1);
                realApp_clear(coeffn);

                if (metadatas_haveToCount(meta)) 
                    metadatas_add_time_Anticip(meta, (double) (clock() - start2) );
           }
            
        }
        iteration +=1;
    }
    
    /*for test */
//     if ((metadatas_useCountSols(meta))&&(restemp>0)&&(discard)){
//         compApp_poly_set(cacheApp_workref(cache),pApprox); 
//         cacheApp_nbItref(cache) = nbGraeffe;
//     }
    /*end for test */
    
    compApp_poly_clear(pApprox);
    realApp_clear(sum);
    
    if ((restemp==0)||(restemp==-1)) res.nbOfSol = -1;
    if (metadatas_haveToCount(meta))
        metadatas_add_Test( meta, depth, (restemp==1), discard, inNewton, 1, nbTaylorsRepeted, nbGraeffe, 
                            nbGraeffeRepeted, (int) res.appPrec, (double) (clock() - start) );
        
//     if (discard)
//         printf(" prec for discarding test: %d\n", (int) res.appPrec );
//     else
//         printf(" --- prec for validating test: %d\n", (int) res.appPrec );
    return res;
    
}

/* EXPERIMENTAL */
#ifdef CCLUSTER_EXPERIMENTAL
void tstar_evaluate( compApp_t res, const compApp_poly_t p, const compApp_t point, slong prec, metadatas_t meta, int depth){
        chronos_tic_Evaluat(metadatas_chronref(meta));
        compApp_poly_evaluate(res, p, point, prec);
        chronos_toc_Evaluat(metadatas_chronref(meta));
        counters_add_Eval( metadatas_countref(meta), 1, depth );
}

void tstar_evaluate_horner( compApp_t res, const compApp_poly_t p, const compApp_t point, slong prec, metadatas_t meta, int depth){
        chronos_tic_Evaluat(metadatas_chronref(meta));
        compApp_poly_evaluate_horner(res, p, point, prec);
        chronos_toc_Evaluat(metadatas_chronref(meta));
        counters_add_Eval( metadatas_countref(meta), 1, depth );
}

void tstar_evaluate2( compApp_t res, compApp_t res2, const compApp_poly_t p, const compApp_t point, slong prec, metadatas_t meta, int depth){
        chronos_tic_Evaluat(metadatas_chronref(meta));
        compApp_poly_evaluate2(res, res2, p, point, prec);
        chronos_toc_Evaluat(metadatas_chronref(meta));
        counters_add_Eval( metadatas_countref(meta), 2, depth );
}

void tstar_getDerivative( compApp_poly_t res, cacheApp_t cache, slong prec, slong order, metadatas_t meta){
        chronos_tic_Derivat(metadatas_chronref(meta));
        compApp_poly_set(res, cacheApp_getDerivative ( cache, prec, order ));
        chronos_toc_Derivat(metadatas_chronref(meta));
}
#endif
/* DEPRECATED */
// tstar_res tstar_count_nb_Sols( cacheApp_t cache,
//                                const compDsk_t d,
//                                int nb_sols,   /*the number of sols in d          */
//                                slong prec,        /*the "default" arithmetic precision              */
//                                int depth,         /*the depth for counter                           */
//                                metadatas_t meta){
//     tstar_res res;
//     res.nbOfSol = -1;
//     res.appPrec = prec;
//     int N = (int) 5+ceil(log2(1+log2(cacheApp_getDegree(cache))));
//         
//     realRat_t factor;
//     realRat_init(factor);
//     realApp_t sum;
//     realApp_init(sum);
//     
//     realRat_set_si(factor, 2, 3);
// //     slong twoToTheN = 2^(cacheApp_nbItref(cache));
//     slong twoToTheN = 0x1<<(cacheApp_nbItref(cache));
//     realRat_pow_si(factor, factor, twoToTheN);
// //     printf("nb It: %d, twoToTheN: %d, factor:",cacheApp_nbItref(cache), twoToTheN); realRat_print(factor); printf("\n");
//     
//     compApp_poly_scale_realRat_in_place( cacheApp_workref(cache)->coeffs, factor, cacheApp_workref(cache)->length, res.appPrec );
//     
//     compApp_poly_sum_abs_coeffs( sum, cacheApp_workref(cache), res.appPrec );
//     
//     res.nbOfSol = compApp_poly_TkGtilda_with_sum( cacheApp_workref(cache), sum, nb_sols, res.appPrec);
//     
// //     if (((res.nbOfSol==-1)||(res.nbOfSol==0))&&(cacheApp_nbItref(cache)<N)){
// // //         printf("Ici\n");
// // //         tstar_graeffe_iterations_inplace( cacheApp_workref(cache), N-cacheApp_nbItref(cache), res.appPrec, meta);
// //         tstar_graeffe_iterations_inplace( cacheApp_workref(cache), 1, res.appPrec, meta);
// //         compApp_poly_sum_abs_coeffs( sum, cacheApp_workref(cache), res.appPrec );
// //         res.nbOfSol = compApp_poly_TkGtilda_with_sum( cacheApp_workref(cache), sum, nb_sols, res.appPrec);
// //     }
//     
//     realRat_clear(factor);
//     realApp_clear(sum);
// //     printf("res.nbOfSol: %d\n", res.nbOfSol);
//     return res;
// }

/*for julia*/
/*int tstar_res_getNbOfSol_forJulia( tstar_res r) { return r.nbOfSol; }   */
/*slong tstar_res_getAppPrec_forJulia( tstar_res r) { return r.appPrec; } */
