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

void ccluster_forJulia_func( connCmp_list_t qResults, 
                             void(*func)(compApp_poly_t, slong), 
                             const compBox_t initialBox, 
                             const realRat_t eps, 
                             char * stratstr,
                             int nbThreads,
                             int verb){
    
//     printf("ccluster_forJulia.c: ccluster_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
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
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );

    ccluster_algo( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_forJulia_func: end\n");
}

void ccluster_global_forJulia_func( connCmp_list_t qResults, 
                                    void(*func)(compApp_poly_t, slong), 
                                    compBox_t initialBox,
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb){
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automatically set initialBox */
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
    /* no root radii */
    strategies_set_RootRadii( strat, 0 );
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    ccluster_algo_global( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_func: end\n");
}

void ccluster_forJulia_compRat_poly( connCmp_list_t qResults, 
                                     const compRat_poly_t poly,
                                     const compBox_t initialBox, 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){
    
//     printf("ccluster_forJulia.c: ccluster_forJulia_compRat_poly: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init_compRat_poly(cache, poly);
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
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    ccluster_algo( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_forJulia_compRat_poly: end\n");
}

void ccluster_global_forJulia_compRat_poly( connCmp_list_t qResults, 
                                     const compRat_poly_t poly,
                                     compBox_t initialBox, 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){
 
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_compRat_poly: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init_compRat_poly(cache, poly);
    strategies_init(strat);
    
    /* automatically set initialBox */
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
    /* no root radii */
    strategies_set_RootRadii( strat, 0 );
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_compRat_poly: end\n");
}

void ccluster_forJulia_realRat_poly( connCmp_list_t qResults, 
                                    const realRat_poly_t poly, 
                                    const compBox_t initialBox, 
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb){
//     printf("ccluster_forJulia.c: ccluster_forJulia_realRat_poly: begin\n");
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p, poly);
    
    /* call */
    ccluster_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
//     printf("ccluster_forJulia.c: ccluster_forJulia_realRat_poly: end\n");
}

void ccluster_global_forJulia_realRat_poly( connCmp_list_t qResults, 
                                            const realRat_poly_t poly,  
                                            compBox_t initialBox,
                                            const realRat_t eps, 
                                            char * stratstr,
                                            int nbThreads,
                                            int verb){
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_realRat_poly: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    compBox_t initBox;
    compBox_init(initBox);
    
    compAnn_list_t qAnn;
    compAnn_list_t qAnn1;
    compAnn_list_t qAnn2;
    
    /* set cache and strategy */
    cacheApp_init_realRat_poly ( cache, poly);
    cacheApp_canonicalise( cache );
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs: */
    /*                             input poly IS real */
    /*                             input initial box is centered in 0 */
    strategies_set_realCoeffs(strat, 1);
    
    /* automatically set initialBox */
    compBox_set_si(initBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); 
        if (realRat_is_zero(compBox_bwidthref(initialBox))) {
            printf("; use 1 instead");
        }
        printf("\n");
    }
    if (realRat_is_zero(compBox_bwidthref(initBox))) {
        realRat_set_si(compBox_bwidthref(initBox), 1, 1);
    }
    realRat_mul_si(compBox_bwidthref(initBox), compBox_bwidthref(initBox), 2);
    
    metadatas_init(meta, initBox, strat, verb);
    
    /*copy initBox in initialBox */
    compBox_set(initialBox, initBox);
    
    
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
        ccluster_algo_global_rootRadii( qResults, NULL, qAnn, qAnn1, qAnn2, initBox, eps, cache, meta);
    }
    else
        ccluster_algo_global( qResults, NULL, initBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_clear(qAnn);
        compAnn_list_clear(qAnn1);
        compAnn_list_clear(qAnn2);
    }
    
    compBox_clear(initBox);
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_realRat_poly: end\n");
}

void ccluster_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                               const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                               const compBox_t initialBox, 
                                               const realRat_t eps, 
                                               char * stratstr,
                                               int nbThreads,
                                               int verb){
//     printf("ccluster_forJulia.c: ccluster_forJulia_realRat_poly_real_imag: begin\n");
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    ccluster_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
//     printf("ccluster_forJulia.c: ccluster_forJulia_realRat_poly_real_imag: end\n");
}

void ccluster_global_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                                      const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                                      compBox_t initialBox,
                                                      const realRat_t eps, 
                                                      char * stratstr,
                                                      int nbThreads,
                                                      int verb){
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_realRat_poly_real_imag: begin\n");
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    ccluster_global_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_realRat_poly_real_imag: end\n");
}

void ccluster_forJulia_refine( connCmp_list_t qResults,
                               connCmp_list_t qMainLoop,
                               void(*func)(compApp_poly_t, slong), 
                               const compBox_t initialBox, 
                               const realRat_t eps, 
                               char * stratstr,
                               int nbThreads,
                               int verb){
//     printf("ccluster_forJulia.c: ccluster_forJulia_refine: begin\n");
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    /* automaticly set realCoeffs: realCoeffs not implemented for refine */
    strategies_set_powerSums(strat, 0);
    strategies_set_realCoeffs(strat, 0);
    strategies_set_RootRadii( strat, 0 );
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_refine( qResults, qMainLoop, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
//     printf("ccluster_forJulia.c: ccluster_forJulia_refine: end\n");
}

/* res = 1: OK */
/* res = -1: lcf vanishes; may miss solutions with huge norm */
int ccluster_global_forJulia_forTcluster_func( connCmp_list_t qResults, 
                                                void(*func)(compApp_poly_t, slong), 
                                                compBox_t initialBox,
                                                const realRat_t eps, 
                                                char * stratstr,
                                                int nbThreads,
                                                int verb){
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_forTcluster_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automatically set initialBox */
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    int res = cacheApp_root_bound_unsure ( compBox_bwidthref(initialBox), cache );
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
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    ccluster_algo_global( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_global_forJulia_forTcluster_func: end\n");
    
    return res;
}

void ccluster_forJulia_draw( connCmp_list_t qResults, 
                            compBox_list_t qDiscarded, 
                            void(*func)(compApp_poly_t, slong), 
                            const compBox_t initialBox, 
                            const realRat_t eps, 
                            char * stratstr,
                            int nbThreads,
                            int verb){
    
    //     printf("ccluster_forJulia.c: ccluster_forJulia_draw: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
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
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo_draw( qResults, qDiscarded, initialBox, eps, cache, meta);
    ccluster_algo( qResults, qDiscarded, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster_forJulia.c: ccluster_forJulia_draw: end\n");
}

/*DEPRECATED:*/

void ccluster_interface_forJulia ( connCmp_list_t qResults, 
                                   void(*func)(compApp_poly_t, slong), 
                                   const compBox_t initialBox, 
                                   const realRat_t eps, 
                                   int st, 
                                   int verb){
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_interface_forJulia_draw( connCmp_list_t qResults, 
                                       compBox_list_t qDiscarded, 
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb){
    
//     printf("ccluster_interface.c: ccluster_interface_forJulia_draw \n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo_draw( qResults, qDiscarded, initialBox, eps, cache, meta);
    ccluster_algo( qResults, qDiscarded, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
}

void ccluster_interface_forJulia_realRat_poly( connCmp_list_t qResults, 
                                              const realRat_poly_t poly, 
                                              const compBox_t initialBox, 
                                              const realRat_t eps, 
                                              int st, 
                                              int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p, poly);
    
    /* call */
    ccluster_interface_forJulia_compRat_poly( qResults, p, initialBox, eps, st, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_interface_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                                         const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                                         const compBox_t initialBox, 
                                                         const realRat_t eps, 
                                                         int st, 
                                                         int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    ccluster_interface_forJulia_compRat_poly( qResults, p, initialBox, eps, st, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_interface_forJulia_compRat_poly( connCmp_list_t qResults, 
                                              const compRat_poly_t poly,
                                              const compBox_t initialBox, 
                                              const realRat_t eps, 
                                              int st, 
                                              int verb){
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init_compRat_poly(cache, poly);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_refine_forJulia( connCmp_list_t qResults,
                               connCmp_list_t qMainLoop,
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st&(0x1<<7), st>>8);
    
    /* automaticly set realCoeffs: realCoeffs not implemented for refine */
    strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_refine( qResults, qMainLoop, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
}
