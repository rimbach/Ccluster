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
int d_global;
int b_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

void evaluateMignotteFast( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec ){
    
    compApp_t temp1, temp2;
    compApp_init(temp1);
    compApp_init(temp2);
    
    compApp_set_si(temp1, 2);
    compApp_pow_si(temp1, temp1, b_global, prec);
    compApp_mul(temp1, temp1, point, prec);
    compApp_sub_ui(temp1, temp1, 1, prec); /* 2^a*x-1 */
    
    compApp_pow_si(temp2, temp1, 2, prec);
    compApp_mul_si(temp2, temp2, 2, prec);
    compApp_pow_si(dest, point, d_global, prec);
    compApp_sub(dest, dest, temp2, prec);
    
    compApp_set_si(temp2, 2);
    compApp_pow_si(temp2, temp2, b_global + 2, prec);
    compApp_mul(temp2, temp2, temp1, prec);
    compApp_pow_si(dest2, point, d_global-1, prec);
    compApp_mul_si(dest2, dest2, d_global, prec);
    compApp_sub(dest2, dest2, temp2, prec);
    
    compApp_clear(temp1);
    compApp_clear(temp2);
    
    return;
}

int main(int argc, char **argv){
    
    if (argc<7){
        printf("usage: %s degree bitsize ", argv[0]);
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
    int bitsize;
//     int st;
    char * st;
    int verbosity;
    int nbthreads = 1;
    int output = 16;
    
    compBox_t bInit;
    realRat_t eps;
    
    compBox_init(bInit);
    realRat_init(eps);
    
    parse = parse*scan_degree( argv[1], &degree);
    parse = parse*scan_bitsize( argv[2], &bitsize);
    parse = parse*scan_initialBox( argv[3], bInit );
    parse = parse*scan_epsilon( argv[4], eps );
//     parse = parse*scan_strategy(argv[5], &st );
    st = argv[5];
    parse = parse*scan_verbosity(argv[6], &verbosity );
    
    if (argc>=8) {
        parse = parse*scan_nbthreads(argv[7], &nbthreads );
    }
//     nbthreads = (nbthreads<<6);
//     int add_temp = (st>>7)<<17;
//     st = st&((0x1<<7)-1);
//     st = st + nbthreads + add_temp;
    
    realRat_poly_t pmign;
    realRat_poly_init(pmign);
    compRat_poly_init(p_global);
    
    d_global = degree;
    b_global = bitsize;
    
    if (parse==1) {
        mignotte_polynomial(pmign, degree, bitsize);
        compRat_poly_set_realRat_poly(p_global,pmign);
        
//         printf("poly: "); compRat_poly_print(p_global); printf("\n");
    
//         ccluster_interface_func( getApprox, bInit, eps, st, verbosity);
//         ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, verbosity);
        ccluster_interface_func_eval( getApprox, evaluateMignotteFast, bInit, eps, st, nbthreads, verbosity);
//         ccluster_interface_funcPS( getApprox, NULL, bInit, eps, st, nbthreads, verbosity);
    }
    
    realRat_poly_clear(pmign);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
