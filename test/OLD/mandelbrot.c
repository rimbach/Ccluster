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

#include "./parseArgs.h"

compRat_poly_t p_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

int main(int argc, char **argv){
    
    if (argc<6){
        printf("usage: %s iteration ", argv[0]);
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
    char * st;
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
//     parse = parse*scan_strategy(argv[4], &st );
    st = argv[4];
    parse = parse*scan_verbosity(argv[5], &verbosity );
    
    if (argc>=7) {
        parse = parse*scan_nbthreads(argv[6], &nbthreads );
    } 
    
//     nbthreads = (nbthreads<<6);
//     int add_temp = (st>>7)<<17;
//     st = st&((0x1<<7)-1);
//     st = st + nbthreads + add_temp;
    
    realRat_poly_t pmand, pone, px;
    realRat_poly_init(pmand);
    realRat_poly_init(pone);
    realRat_poly_init(px);
    
    compRat_poly_init(p_global);
        
    if (parse==1) {
        
        realRat_poly_one(pmand);
        realRat_poly_one(pone);
        realRat_poly_zero(px);
        realRat_poly_set_coeff_si_ui(px, 1, 1, 1);
        
        for (int i = 1; i<=degree; i++) {
            realRat_poly_pow(pmand, pmand, 2);
            realRat_poly_mul(pmand, pmand, px);
            realRat_poly_add(pmand, pmand, pone);
        }
        
        compRat_poly_set_realRat_poly(p_global,pmand);
        
//         ccluster_interface_func( getApprox, bInit, eps, st, verbosity);
        ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, output, verbosity);
    }
    
    realRat_poly_clear(pmand);
    realRat_poly_clear(pone);
    realRat_poly_clear(px);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
