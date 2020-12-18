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

#include "cauchy/cauchy.h"

void cauchy_global_interface_func( void(*func)(compApp_poly_t, slong), 
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
//     if (metadatas_usePowerSums(meta)) {
//         if (strat->_pwSuNbPs>1)
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
//         else 
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
            metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 4, 3, 1, verb );
//     }
    
    cauchy_algo_global( qRes, bDis, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_cauchy_fprint(stdout, meta, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 0);
    } else if (output==-3){
//         connCmp_list_gnuplot(stdout, qRes, meta, 0);
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
    compBox_clear(initialBox);
}

/* version with function for fast evaluation */
void cauchy_global_interface_func_eval( void(*func)(compApp_poly_t, slong),
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
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
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
//     if (metadatas_usePowerSums(meta)) 
//         metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 2, 1, 1, verb );
//     metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 2, 1, 1, verb );
    metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 4, 3, 1, verb );
//     metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 16, 15, 1, verb );
    
    cauchy_algo_global( qRes, NULL, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_cauchy_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    
    compBox_clear(initialBox);
}
