/* ************************************************************************** */
/*  Copyright (C) 2018 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <string.h>
#include <stdio.h>
#include "polynomials/compRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"
#include "ccluster/cclusterDAC.h"

#include "./parseArgs.h"

compRat_poly_t p_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

int main(int argc, char **argv){
    
    if (argc<6){
        printf("usage: %s degree ", argv[0]);
        printf("initial_box epsilon strategy verbosity\n");
        printf("initial_box: for instance 0,1,1,2,100,1 i.e. the square centered in 0/1 +i*(1/2) of width 100/1\n");
        printf("eps:         for instance1,100 (1/100) or -53 (1/2^(-53))\n");  
        printf("strategy:    default for default strategy\n");
        printf("             test for testing mode\n");
        printf("verbosity:   0: nothing\n");
        printf("             1: abstract of input and output\n");
        printf("             2: detailed reports concerning algorithm\n");
        printf("             3: same as 2 + prints the clusters\n");
        return -1;
    }
    
    int parse = 1;
    int degree;
    int st;
    int verbosity;
    int nbthreads = 1;
    int output = 16;
    
    compBox_t bInit;
    realRat_t eps;
    
    compBox_init(bInit);
    realRat_init(eps);
    
    parse = parse*scan_degree( argv[1], &degree);
    parse = parse*scan_initialBox( argv[2], bInit );
    parse = parse*scan_epsilon( argv[3], eps );
    parse = parse*scan_strategy(argv[4], &st );
    parse = parse*scan_verbosity(argv[5], &verbosity );
    
    if (argc>=7) {
        parse = parse*scan_nbthreads(argv[6], &nbthreads );
    }
    nbthreads = (nbthreads<<5);
    int add_temp = (st>>6)<<16;
    st = st&((0x1<<6)-1);
    st = st + nbthreads + add_temp;
    
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    compRat_poly_init(p_global);
        
    if (parse==1) {
        
        bernoulli_polynomial( pbern, degree);
        compRat_poly_set_realRat_poly(p_global,pbern);
        
//         ccluster_interface_func( getApprox, bInit, eps, st, verbosity);
        
        cacheApp_t cache;
        strategies_t strat;
        metadatas_t meta;
        connCmp_list_t qRes, qResTemp;
        connCmp_list_t qMainLoop, discardedCcs;
        connCmp_ptr ccur;
        
        cacheApp_init(cache, getApprox);
        strategies_init(strat);
        strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
        
        metadatas_init(meta, bInit, strat, verbosity);
        connCmp_list_init(qRes);
        connCmp_list_init(qResTemp);
        connCmp_list_init(qMainLoop);
        connCmp_list_init(discardedCcs);
        
        ccluster_DAC_first( qResTemp, qMainLoop, discardedCcs, bInit, eps, cache, meta);
        metadatas_count(meta);
        metadatas_fprint(stdout, meta, eps);
        
        while (!connCmp_list_is_empty(qMainLoop)) {
            ccur = connCmp_list_pop(qResTemp);
            connCmp_list_push(qRes, ccur);
            ccluster_DAC_next( qResTemp, qMainLoop, discardedCcs, eps, cache, meta);
            metadatas_count(meta);
            metadatas_fprint(stdout, meta, eps);
        }
        
        if (!connCmp_list_is_empty(qResTemp)){
            ccur = connCmp_list_pop(qResTemp);
            connCmp_list_push(qRes, ccur);
        }
        
        printf("number of clusters: %d\n", connCmp_list_get_size(qRes));
        if (verbosity>=3) {
            connCmp_list_print_for_results(stdout, qRes, meta);
        }
        
        cacheApp_clear(cache);
        strategies_clear(strat);
        metadatas_clear(meta);
        connCmp_list_clear(qRes);
        connCmp_list_clear(qResTemp);
        
        connCmp_list_clear(qMainLoop);
        connCmp_list_clear(discardedCcs);
    
//         ccluster_DAC_first_interface_func( getApprox, bInit, eps, st, verbosity);
    }
    
    realRat_poly_clear(pbern);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    return 0;
}
