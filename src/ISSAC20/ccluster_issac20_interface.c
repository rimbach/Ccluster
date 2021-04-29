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

#include "ISSAC20/ccluster_issac20.h"

/* experimental version */
int ccluster_issac20_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                          const realRat_t eps, 
                                          char * stratstr,
                                          int nbThreads,
                                          int output,
                                          int verb){

    int res = 0;
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
    strategies_set_powerSums( strat, 1 );
    strategies_set_useNewton( strat, 0 );
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
//     metadatas_setNbPowerSums(meta, 2);
//     metadatas_setIsoRatio_si(meta, 2, 1);
// //     metadatas_setIsoRatio_si(meta, 4, 3);
// //     if ( metadatas_pwSuTest(meta) )
//     ccluster_initialize_pwSuTest(NULL, meta, cache, verb);
//     if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 4, 3, 3, verb );
    
    res = ccluster_issac20_algo_global( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    
    /*check if number of sols is equal to degree*/
    if ( metadatas_getNbSolutions(meta) != (int) cacheApp_getDegree(cache) ){
        res = 2;
        if (metadatas_getVerbo(meta)>=2) {
            printf("FAILURE: %d solutions found, degree of input pol: %ld\n", metadatas_getNbSolutions(meta), cacheApp_getDegree(cache));
        }
    }
    
    metadatas_issac20_fprint(stdout, res, meta, eps);
    
    if ((output==-1)||(output>0)) {
//         connCmp_list_print_for_results(stdout, qRes, meta);
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_clear(initialBox);
    
    return res;
}

int ccluster_issac20_global_interface_func_eval( void(*func)(compApp_poly_t, slong), 
                                              void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                                          const realRat_t eps, 
                                          char * stratstr,
                                          int nbThreads,
                                          int output,
                                          int verb){
    
    int res = 0;
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
    strategies_set_powerSums( strat, 1 );
    strategies_set_useNewton( strat, 0 );
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
//     metadatas_setNbPowerSums(meta, 2);
//     metadatas_setIsoRatio_si(meta, 2, 1);
// //     metadatas_setIsoRatio_si(meta, 4, 3);
// //     if ( metadatas_pwSuTest(meta) )
//     ccluster_initialize_pwSuTest(NULL, meta, cache, verb);
//     if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 4, 3, 3, verb );
    
    res = ccluster_issac20_algo_global( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    
    /*check if number of sols is equal to degree*/
    if ( metadatas_getNbSolutions(meta) != (int) cacheApp_getDegree(cache) ){
        res = 2;
        if (metadatas_getVerbo(meta)>=2) {
            printf("FAILURE: %d solutions found, degree of input pol: %ld\n", metadatas_getNbSolutions(meta), cacheApp_getDegree(cache));
        }
    }
    
    metadatas_issac20_fprint(stdout, res, meta, eps);
    
    if ((output==-1)||(output>0)) {
//         connCmp_list_print_for_results(stdout, qRes, meta);
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_clear(initialBox);
    
    return res;
    
}
