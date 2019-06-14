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

void ccluster_interface_func( void(*func)(compApp_poly_t, slong), const compBox_t initialBox, const realRat_t eps, int st, int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4),0, (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), (st&( ((0x1<<10)-1)<<6 ))>>6, st>>17);
    
    /* automaticly set realCoeffs */
//     printf("cacheApp_is_real(cache): %d\n", cacheApp_is_real(cache));
    if (cacheApp_is_real(cache)==0)
        strategies_set_realCoeffs(strat, 0);
//     printf("strategies_realCoeffs(strat): %d\n", strategies_realCoeffs(strat));
    
    metadatas_init(meta, initialBox, strat, verb);
    connCmp_list_init(qRes);
    
    ccluster_algo( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
//         connCmp_list_print_for_results(stdout, qRes, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
}

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
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4),0, (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), (st&( ((0x1<<10)-1)<<6 ))>>6, st>>17);
    metadatas_init(meta, initialBox, strat, verb);
    connCmp_list_init(qRes);
    
    ccluster_algo( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
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
    
        /* automaticly set realCoeffs */
//     printf("cacheApp_is_real(cache): %d\n", cacheApp_is_real(cache));
    if (cacheApp_is_real(cache)==0)
        strategies_set_realCoeffs(strat, 0);
//     printf("strategies_realCoeffs(strat): %d\n", strategies_realCoeffs(strat));
    
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
