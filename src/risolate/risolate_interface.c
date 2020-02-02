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

#include "ccluster/ccluster.h"
#include "risolate/risolate.h"

void risolate_interface_poly( const realRat_poly_t poly,
                              const compBox_t initialBox, 
                              const realRat_t eps, 
                              char * stratstr,
                              int nbThreads,
                              int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init_realRat_poly ( cache, poly);
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
//     /* automatically set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    connCmp_list_init(qRes);
    
    risolate_algo( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
}

void risolate_global_interface_poly( const realRat_poly_t poly,
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init_realRat_poly ( cache, poly);
    strategies_init(strat);
    
    /* automaticly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
//     /* automaticly set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
    
    connCmp_list_init(qRes);
    
    metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    risolate_algo_global( qRes, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_clear(initialBox);
}
