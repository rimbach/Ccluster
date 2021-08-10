/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "risolate/risolate.h"

void risolate_forJulia_realRat_poly( connCmp_list_t qResults, 
                                     const realRat_poly_t poly,  
                                     const compBox_t initialBox,
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int verb){
    
//     printf("risolate_forJulia.c: risolate_forJulia_realRat_poly: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    compBox_t initBox;
    compBox_init(initBox);
    compBox_set(initBox, initialBox);
    /* force imaginary part of initial box to be centered in 0 */
    realRat_set_si(compRat_imagref(compBox_centerref(initBox)),0,1);
    
    /* set cache and strategy */
    cacheApp_init_realRat_poly ( cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#risolate_forJulia.c: risolate_forJulia_realRat_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        compBox_clear(initBox);
        return;
    }
    
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* separation bound */
    realRat_t sepBound;
    realRat_init(sepBound);
    cacheApp_separation_bound ( sepBound, cache);
    if (verb>=3) {
        printf("separation bound: "); realRat_print(sepBound); printf("\n");
    }
    metadatas_setSepBound(meta, sepBound);
    
    if (cacheApp_getDegree(cache)>0)
        risolate_algo( qResults, NULL, initBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_risolate_fprint(stdout, meta, cache, eps);
    
    if (verb>=3) {
        connCmp_list_risolate_print_for_results(stdout, qResults, meta);
    }
    
    compBox_clear(initBox);
    realRat_clear(sepBound);
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("risolate_forJulia.c: risolate_forJulia_realRat_poly: end\n");
}

void risolate_global_forJulia_realRat_poly( connCmp_list_t qResults, 
                                            const realRat_poly_t poly,  
                                            compBox_t initialBox,
                                            const realRat_t eps, 
                                            char * stratstr,
                                            int nbThreads,
                                            int verb){
    
//     printf("risolate_forJulia.c: risolate_global_forJulia_realRat_poly: begin\n");
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    compAnn_list_t qAnn;
    
    /* set cache and strategy */
    cacheApp_init_realRat_poly ( cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#risolate_forJulia.c: risolate_global_forJulia_realRat_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        return;
    }
    
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
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
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* separation bound */
    realRat_t sepBound;
    realRat_init(sepBound);
    cacheApp_separation_bound ( sepBound, cache);
    if (verb>=3) {
        printf("separation bound: "); realRat_print(sepBound); printf("\n");
    }
    metadatas_setSepBound(meta, sepBound);
    
    if (metadatas_useRootRadii(meta)) {
        compAnn_list_init(qAnn);
        if (cacheApp_getDegree(cache)>0)
            risolate_algo_global_rootRadii( qResults, NULL, qAnn, initialBox, eps, cache, meta);
    }
    else
        if (cacheApp_getDegree(cache)>0)
            risolate_algo_global( qResults, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_risolate_fprint(stdout, meta, cache, eps);
    
    if (verb>=3) {
        connCmp_list_risolate_print_for_results(stdout, qResults, meta);
    }
    
    if (metadatas_useRootRadii(meta)){
        connCmp_list_iterator it = connCmp_list_begin(qResults);
        while (it != connCmp_list_end()){
            compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref(connCmp_list_elmt(it)) );
            while (itb != compBox_list_end() ) {
                compBox_clear_annulii( compBox_list_elmt(itb) );
                itb = compBox_list_next(itb);
            }
            it = connCmp_list_next(it);
        }
        compAnn_list_clear(qAnn);
    }
    
    realRat_clear(sepBound);
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    
//     printf("risolate_forJulia.c: risolate_global_forJulia_realRat_poly: end\n");
}
