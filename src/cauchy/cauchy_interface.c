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
                                     const realRat_t isoRatio,
                                     int nbPows,
                                     char * stratstr,
                                     int nbThreads,
                                     int output,
                                     int verb){

    cacheApp_t cache;
    cacheCauchy_t cacheCau;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    metadatas_init(meta, initialBox, strat, verb);
    
    cacheCauchy_init(cacheCau, NULL, cacheApp_getDegree(cache), isoRatio, (slong) nbPows, meta);
    
    /* automaticly set initialBox */
    /* with evaluation function */
    cauchyRootRadii_root_bound( compBox_bwidthref(initialBox), cacheCau, cache, meta );
    if (verb>=3) {
        printf("#root bound with eval function: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    /* with coefficients */
//     cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
//     if (verb>=3) {
//         printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
//     }
//     realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    metadatas_setInitBox(meta, initialBox);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(metadatas_stratref(meta), 0);
    
//     strategies_set_realCoeffs(metadatas_stratref(meta), 0);
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
        
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    /* initialize power sums */
//     if (metadatas_usePowerSums(meta)) {
//         if (strat->_pwSuNbPs>1)
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
//         else 
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 4, 3, 1, verb );
//     }
    
    cauchy_algo_global( qRes, bDis, initialBox, eps, cache, cacheCau, meta);
    
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
    cacheCauchy_clear(cacheCau);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
    compBox_clear(initialBox);
}

void cauchy_global_interface_realRat_poly( const realRat_poly_t poly,
                                           const realRat_t eps,
                                           const realRat_t isoRatio,
                                           int nbPows,
                                           char * stratstr,
                                           int nbThreads,
                                           int output,
                                           int verb){
    cacheApp_t cache;
    cacheCauchy_t cacheCau;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    cacheApp_init_realRat_poly(cache, poly);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    metadatas_init(meta, initialBox, strat, verb);
    
    cacheCauchy_init(cacheCau, NULL, cacheApp_getDegree(cache), isoRatio, (slong) nbPows, meta);
    
    /* automaticly set initialBox */
    /* with evaluation function */
    cauchyRootRadii_root_bound( compBox_bwidthref(initialBox), cacheCau, cache, meta );
    if (verb>=3) {
        printf("#root bound with eval function: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    /* with coefficients */
//     cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
//     if (verb>=3) {
//         printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
//     }
//     realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    metadatas_setInitBox(meta, initialBox);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(metadatas_stratref(meta), 0);
    
//     strategies_set_realCoeffs(metadatas_stratref(meta), 0);
    
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
        
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    /* initialize power sums */
//     if (metadatas_usePowerSums(meta)) {
//         if (strat->_pwSuNbPs>1)
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, strat->_pwSuNbPs, verb );
//         else 
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//             metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 4, 3, 1, verb );
//     }
    
    cauchy_algo_global( qRes, bDis, initialBox, eps, cache, cacheCau, meta);
    
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
    cacheCauchy_clear(cacheCau);
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
                                   const realRat_t isoRatio,
                                   int nbPows,
                                   char * stratstr,
                                   int nbThreads,
                                   int output,
                                   int verb){

    cacheApp_t cache;
    cacheCauchy_t cacheCau;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    compBox_list_t bDis;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    metadatas_init(meta, initialBox, strat, verb);
    
    cacheCauchy_init(cacheCau, evalFast, cacheApp_getDegree(cache), isoRatio, (slong) nbPows, meta);
    
    /* automaticly set initialBox */
    /* with evaluation function */
    cauchyRootRadii_root_bound( compBox_bwidthref(initialBox), cacheCau, cache, meta );
    if (verb>=3) {
        printf("root bound with eval function: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    /* with coefficients */
//     cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
//     if (verb>=3) {
//         printf("root bound with coefficients : "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
//     }
//     realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    metadatas_setInitBox(meta, initialBox);
    
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
 
    connCmp_list_init(qRes);
    compBox_list_init(bDis);
    
    if (output==-3) 
        metadatas_setDrSub(meta, 1);
    
    cauchy_algo_global( qRes, bDis, initialBox, eps, cache, cacheCau, meta);
    
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
    cacheCauchy_clear(cacheCau);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_list_clear(bDis);
    
    compBox_clear(initialBox);
}
