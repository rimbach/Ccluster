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
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_global_forJulia_func( connCmp_list_t qResults, 
                                    void(*func)(compApp_poly_t, slong), 
                                    compBox_t initialBox,
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb){
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automatically set initialBox */
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo_global( qResults, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_forJulia_compRat_poly( connCmp_list_t qResults, 
                                     const compRat_poly_t poly,
                                     const compBox_t initialBox, 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init_compRat_poly(cache, poly);
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_global_forJulia_compRat_poly( connCmp_list_t qResults, 
                                     const compRat_poly_t poly,
                                     compBox_t initialBox, 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){
//     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init_compRat_poly(cache, poly);
    strategies_init(strat);
    
    /* automatically set initialBox */
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
}

void ccluster_forJulia_realRat_poly( connCmp_list_t qResults, 
                                    const realRat_poly_t poly, 
                                    const compBox_t initialBox, 
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p, poly);
    
    /* call */
    ccluster_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_global_forJulia_realRat_poly( connCmp_list_t qResults, 
                                            const realRat_poly_t poly,  
                                            compBox_t initialBox,
                                            const realRat_t eps, 
                                            char * stratstr,
                                            int nbThreads,
                                            int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p, poly);
    
    /* call */
    ccluster_global_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                               const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                               const compBox_t initialBox, 
                                               const realRat_t eps, 
                                               char * stratstr,
                                               int nbThreads,
                                               int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    ccluster_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_global_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                                      const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                                      compBox_t initialBox,
                                                      const realRat_t eps, 
                                                      char * stratstr,
                                                      int nbThreads,
                                                      int verb){
    /* set the polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    ccluster_global_forJulia_compRat_poly( qResults, p, initialBox, eps, stratstr, nbThreads, verb);
    
    /* clear */
    compRat_poly_clear(p);
}

void ccluster_forJulia_refine( connCmp_list_t qResults,
                               connCmp_list_t qMainLoop,
                               void(*func)(compApp_poly_t, slong), 
                               const compBox_t initialBox, 
                               const realRat_t eps, 
                               char * stratstr,
                               int nbThreads,
                               int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    /* automaticly set realCoeffs: realCoeffs not implemented for refine */
    strategies_set_realCoeffs(strat, 0);
    
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
}

void ccluster_forJulia_draw( connCmp_list_t qResults, 
                            compBox_list_t qDiscarded, 
                            void(*func)(compApp_poly_t, slong), 
                            const compBox_t initialBox, 
                            const realRat_t eps, 
                            char * stratstr,
                            int nbThreads,
                            int verb){
    
    //     printf("ccluster_interface.c: ccluster_interface_forJulia_draw \n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo_draw( qResults, qDiscarded, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
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
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, initialBox, eps, cache, meta);
    
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
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo_draw( qResults, qDiscarded, initialBox, eps, cache, meta);
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
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    ccluster_algo( qResults, initialBox, eps, cache, meta);
    
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
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    
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
