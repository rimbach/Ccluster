/* ************************************************************************** */
/*  Copyright (C) 2019 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "ccluster/ccluster.h"
    
void ccluster_interface_func( void(*func)(compApp_poly_t, slong), 
                              const compBox_t initialBox, 
                              const realRat_t eps, 
                              char * stratstr,
                              int nbThreads,
                              int output,
                              int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* no root radii */
    strategies_set_RootRadii( strat, 0 );
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
    if (metadatas_usePowerSums(meta)) {
        if (strat->_pwSuNbPs>1)
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
        else 
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    }
    
    /* check if input polynomial is close to zero */
    if (cacheApp_is_near_zero ( cache )) {
        printf("#ccluster_interface.c: ccluster_interface_func \n");
        printf("# input polynomial is nearly zero polynomial: can have infinite number of roots \n");
    }
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    if (cacheApp_getDegree(cache)>0)
        ccluster_algo( qRes, bDis, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 1);
    } else if (output==-3){
//         connCmp_list_gnuplot(stdout, qRes, meta, 1);
        connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
}

void ccluster_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int output,
                                     int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    compAnn_list_t qAnn;
    compAnn_list_t qAnn1;
    compAnn_list_t qAnn2;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automaticly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    
    /* check if input polynomial is close to zero */
    if (cacheApp_is_near_zero ( cache )) {
        printf("#ccluster_interface.c: ccluster_global_interface_func \n");
        printf("# input polynomial is nearly zero polynomial: can have infinite number of roots \n");
    }
    
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
            printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); 
            if (realRat_is_zero(compBox_bwidthref(initialBox))) {
                printf("; use 1 instead");
            }
            printf("\n");
    }
    
    if (realRat_is_zero(compBox_bwidthref(initialBox))) {
        realRat_set_si(compBox_bwidthref(initialBox), 1, 1);
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* no root radii */
    strategies_set_RootRadii( strat, 0 );
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    /* initialize power sums */
    if (metadatas_usePowerSums(meta)) {
        if (strat->_pwSuNbPs>1)
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
        else 
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_init(qAnn);
        compAnn_list_init(qAnn1);
        compAnn_list_init(qAnn2);
        if (cacheApp_getDegree(cache)>0)
            ccluster_algo_global_rootRadii( qRes, bDis, qAnn, qAnn1, qAnn2, initialBox, eps, cache, meta);
    }
    else
        if (cacheApp_getDegree(cache)>0)
            ccluster_algo_global( qRes, bDis, initialBox, eps, cache, meta);
    
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 0);
    } else if (output==-3){
        if (metadatas_useRootRadii(meta))
            connCmp_list_gnuplot_drawSubdiv_rootRadii(stdout, qRes, bDis, qAnn, qAnn1, qAnn2, meta);
        else
            connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_clear(qAnn);
        compAnn_list_clear(qAnn1);
        compAnn_list_clear(qAnn2);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
    compBox_clear(initialBox);
}

void ccluster_interface_realRat_poly( const realRat_poly_t poly,
                                      const compBox_t initialBox,
                                      const realRat_t eps, 
                                      char * stratstr,
                                      int nbThreads,
                                      int output,
                                      int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    /* set cache and strategy */
    cacheApp_init_realRat_poly ( cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#ccluster_interface.c: ccluster_interface_realRat_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        return;
    }
    strategies_set_str( strat, stratstr, nbThreads );
    /* no root radii */
    strategies_set_RootRadii( strat, 0 );
    
    /* automaticly set realCoeffs: */
    /*                             input poly IS real */
    /*                             input initial box is centered in 0 */
    strategies_set_realCoeffs(strat, 1);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
    if (metadatas_usePowerSums(meta)) {
        if (strat->_pwSuNbPs>1)
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
        else 
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    }
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    if (cacheApp_getDegree(cache)>0)
        ccluster_algo( qRes, bDis, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 1);
    } else if (output==-3){
//         connCmp_list_gnuplot(stdout, qRes, meta, 1);
        connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
}

void ccluster_global_interface_realRat_poly( const realRat_poly_t poly,
                                             const realRat_t eps, 
                                             char * stratstr,
                                             int nbThreads,
                                             int output,
                                             int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    compAnn_list_t qAnn;
    compAnn_list_t qAnn1;
    compAnn_list_t qAnn2;
    
    /* set cache and strategy */
    cacheApp_init_realRat_poly ( cache, poly);
//     cacheApp_canonicalise( cache );
    if (cacheApp_is_zero(cache)) {
        printf("#ccluster_interface.c: ccluster_global_interface_realRat_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        return;
    }
    
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs: */
    /*                             input poly IS real */
    /*                             input initial box is centered in 0 */
    strategies_set_realCoeffs(strat, 1);
    
    /* automaticaly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); 
        if (realRat_is_zero(compBox_bwidthref(initialBox))) {
            printf("; use 1 instead");
        }
        printf("\n");
    }
    if (realRat_is_zero(compBox_bwidthref(initialBox))) {
        realRat_set_si(compBox_bwidthref(initialBox), 1, 1);
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    /* initialize power sums */
    if (metadatas_usePowerSums(meta)) {
        if (strat->_pwSuNbPs>1)
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
        else 
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_init(qAnn);
        compAnn_list_init(qAnn1);
        compAnn_list_init(qAnn2);
        if (cacheApp_getDegree(cache)>0)
            ccluster_algo_global_rootRadii( qRes, bDis, qAnn, qAnn1, qAnn2, initialBox, eps, cache, meta);
    }
    else
        if (cacheApp_getDegree(cache)>0)
            ccluster_algo_global( qRes, bDis, initialBox, eps, cache, meta);
    
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 0);
    } else if (output==-3){
        if (metadatas_useRootRadii(meta))
            connCmp_list_gnuplot_drawSubdiv_rootRadii(stdout, qRes, bDis, qAnn, qAnn1, qAnn2, meta);
        else
            connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_clear(qAnn);
        compAnn_list_clear(qAnn1);
        compAnn_list_clear(qAnn2);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
    compBox_clear(initialBox);
}

/* version with function for fast evaluation */
void ccluster_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                   void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                                   const compBox_t initialBox, 
                                   const realRat_t eps, 
                                   char * stratstr,
                                   int nbThreads,
                                   int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
 
    connCmp_list_init(qRes);
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    if (cacheApp_getDegree(cache)>0)
        ccluster_algo( qRes, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
}

/* version with function for fast evaluation */
void ccluster_global_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                   void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong), 
                                   const realRat_t eps, 
                                   char * stratstr,
                                   int nbThreads,
                                   int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automaticly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); 
        if (realRat_is_zero(compBox_bwidthref(initialBox))) {
            printf("; use 1 instead");
        }
        printf("\n");
    }
    if (realRat_is_zero(compBox_bwidthref(initialBox))) {
        realRat_set_si(compBox_bwidthref(initialBox), 1, 1);
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
 
    connCmp_list_init(qRes);
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    if (cacheApp_getDegree(cache)>0)
        ccluster_algo_global( qRes, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    
    compBox_clear(initialBox);
}

// void ccluster_forJulia_func( connCmp_list_t qResults, 
//                              void(*func)(compApp_poly_t, slong), 
//                              const compBox_t initialBox, 
//                              const realRat_t eps, 
//                              char * stratstr,
//                              int nbThreads,
//                              int verb){
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
//     
//     cacheApp_t cache;
//     strategies_t strat;
//     metadatas_t meta;
//     
//     cacheApp_init(cache, func);
//     strategies_init(strat);
//     
//     strategies_set_str( strat, stratstr, nbThreads );
//     /* automatically set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
//     
//     metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//     
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
//     
//     metadatas_count(meta);
//     metadatas_fprint(stdout, meta, cache, eps);
//     if (verb>=3) {
//         connCmp_list_print_for_results(stdout, qResults, meta);
//     }
//     
//     cacheApp_clear(cache);
//     strategies_clear(strat);
//     metadatas_clear(meta);
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
// }

// void ccluster_global_forJulia_func( connCmp_list_t qResults, 
//                                     void(*func)(compApp_poly_t, slong),  
//                                     const realRat_t eps, 
//                                     char * stratstr,
//                                     int nbThreads,
//                                     int verb){
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
//     
//     cacheApp_t cache;
//     strategies_t strat;
//     metadatas_t meta;
//     
//     cacheApp_init(cache, func);
//     strategies_init(strat);
//     
//     /* automaticly set initialBox */
// //     compBox_t initialBox;
// //     compBox_init(initialBox);
//     compBox_set_si(initialBox, 0,1,0,1,0,1);
//     cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
//     if (verb>=3) {
//         printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
//     }
//     realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
//     
//     strategies_set_str( strat, stratstr, nbThreads );
//     /* automatically set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
//     
//     metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//     
//     ccluster_algo_global( qResults, initialBox, eps, cache, meta);
//     
//     metadatas_count(meta);
//     metadatas_fprint(stdout, meta, cache, eps);
//     if (verb>=3) {
//         connCmp_list_print_for_results(stdout, qResults, meta);
//     }
//     
//     cacheApp_clear(cache);
//     strategies_clear(strat);
//     metadatas_clear(meta);
// //     compBox_clear(initialBox);
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
// }

int ccluster_interface_poly( realRat_t * centerRe, realRat_t * centerIm, int * mults, 
                             const compRat_poly_t poly, 
                             const compBox_t initialBox, 
                             const realRat_t eps, 
                             int st, 
                             int verb){
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init_compRat_poly(cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#ccluster_interface.c: ccluster_interface_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        return 0;
    }
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4),0, (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), (st&( ((0x1<<10)-1)<<6 ))>>6, st>>17);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), (st&( ((0x1<<10)-1)<<7 ))>>7, st>>18);
    strategies_set_int ( strat, 
                         st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), (st&( ((0x1<<10)-1)<<7 ))>>7, st>>18);
    
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    connCmp_list_init(qRes);
    
    if (cacheApp_getDegree(cache)>0)
        ccluster_algo( qRes, NULL, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, cache, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
//         connCmp_list_print_for_results(stdout, qRes, 500, 40, meta);
    }
    
    /* feed the results */
//     int nbClus = connCmp_list_get_size(qRes);
    int nbClus = 0;
    compBox_t containingBox;
    compBox_init(containingBox);
    connCmp_list_iterator it = connCmp_list_begin(qRes);
    while (it!=connCmp_list_end() ) {
        
        connCmp_componentBox( containingBox, connCmp_list_elmt(it), metadatas_initBref(meta));
        realRat_set( centerRe[nbClus], compRat_realref(compBox_centerref(containingBox)) );
        realRat_set( centerIm[nbClus], compRat_imagref(compBox_centerref(containingBox)) );
        mults[nbClus] = connCmp_nSols(connCmp_list_elmt(it));
        
        it = connCmp_list_next(it);
        nbClus++;
    }
    compBox_clear(containingBox);
       
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    
    return nbClus;
}

int ccluster_interface_poly_real( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                  const realRat_poly_t poly, 
                                  const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                  const realRat_t eps, 
                                  int st, 
                                  int verb){
    
    /* initial Box */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_3realRat(initialBox, initialBox_cr, initialBox_ci, initialBox_wi);
    /* polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p,poly);
    
    /* call */
    int nbClus = ccluster_interface_poly( centerRe, centerIm, mults, p, initialBox, eps, st, verb);
    
    /* clear */
    compBox_clear(initialBox);
    compRat_poly_clear(p);
    
    return nbClus;
    
}

int ccluster_interface_poly_real_imag( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                       const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                       const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                       const realRat_t eps, 
                                       int st, 
                                       int verb){
    
    /* initial Box */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_3realRat(initialBox, initialBox_cr, initialBox_ci, initialBox_wi);
    /* polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    int nbClus = ccluster_interface_poly( centerRe, centerIm, mults, p, initialBox, eps, st, verb);
    
    /* clear */
    compBox_clear(initialBox);
    compRat_poly_clear(p);
    
    return nbClus;
    
}
