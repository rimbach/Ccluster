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
                              int output,
                              int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    compBox_t initBox;
    compBox_init(initBox);
    compBox_set(initBox, initialBox);
    /* force imaginary part of initial box to be centered in 0 */
    realRat_set_si(compRat_imagref(compBox_centerref(initBox)),0,1);
    
    
    cacheApp_init_realRat_poly ( cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#risolate_interface.c: risolate_interface_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        compBox_clear(initBox);
        return;
    }
    
    cacheApp_canonicalise( cache );
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
    
    metadatas_init(meta, initBox, strat, verb);
    
    /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );

    /* separation bound */
    realRat_t sepBound;
    realRat_init(sepBound);
    cacheApp_separation_bound ( sepBound, cache);
    if (verb>=3) {
        printf("separation bound: "); realRat_print(sepBound); printf("\n");
    }
    metadatas_setSepBound(meta, sepBound);
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    if (cacheApp_getDegree(cache)>0)
        risolate_algo( qRes, bDis, initBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_risolate_fprint(stdout, meta, eps);
    
    if (output==-2) {
        risolate_connCmp_list_gnuplot(stdout, qRes, meta, 1);
    } else if (output==-3){
        risolate_connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
        connCmp_list_risolate_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    compBox_clear(initBox);
    realRat_clear(sepBound);
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
}

void risolate_global_interface_poly( const realRat_poly_t poly,
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
    
    cacheApp_init_realRat_poly ( cache, poly);
    if (cacheApp_is_zero(cache)) {
        printf("#risolate_interface.c: risolate_global_interface_poly \n");
        printf("# input polynomial is zero polynomial: infinite number of roots \n");
        printf("# exit\n");
        cacheApp_clear(cache);
        return;
    } 
    
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
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    /* separation bound */
    realRat_t sepBound;
    realRat_init(sepBound);
    cacheApp_separation_bound ( sepBound, cache);
    if (verb>=3) {
        printf("separation bound: "); realRat_print(sepBound); printf("\n");
    }
    metadatas_setSepBound(meta, sepBound);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    if (metadatas_useRootRadii(meta)) {
        compAnn_list_init(qAnn);
        if (cacheApp_getDegree(cache)>0)
            risolate_algo_global_rootRadii( qRes, bDis, qAnn, initialBox, eps, cache, meta);
    }
    else
        if (cacheApp_getDegree(cache)>0)
            risolate_algo_global( qRes, bDis, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_risolate_fprint(stdout, meta, eps);
    
    if (output==-2) {
        risolate_connCmp_list_gnuplot(stdout, qRes, meta, 1);
    } else if (output==-3){
        if (metadatas_useRootRadii(meta)) 
            risolate_connCmp_list_gnuplot_drawSubdiv_rootRadii(stdout, qRes, bDis, qAnn, meta);
        else
            risolate_connCmp_list_gnuplot_drawSubdiv(stdout, qRes, bDis, meta);
    } else if (output!=0) {
        connCmp_list_risolate_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    if (metadatas_useRootRadii(meta)){
        compAnn_list_clear(qAnn);
    }
    
    realRat_clear(sepBound);
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
    compBox_clear(initialBox);
}
