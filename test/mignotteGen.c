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
    
    if (argc<8){
        printf("usage: %s degree  bitsize power ", argv[0]);
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
    int power;
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
    parse = parse*scan_power( argv[3], &power);
    parse = parse*scan_initialBox( argv[4], bInit );
    parse = parse*scan_epsilon( argv[5], eps );
//     parse = parse*scan_strategy(argv[6], &st );
    st = argv[6];
    parse = parse*scan_verbosity(argv[7], &verbosity );
    
    if (argc>=9) {
        parse = parse*scan_nbthreads(argv[8], &nbthreads );
    }
//     nbthreads = (nbthreads<<6);
//     int add_temp = (st>>7)<<17;
//     st = st&((0x1<<7)-1);
//     st = st + nbthreads + add_temp;
    
    realRat_poly_t pmign;
    realRat_poly_init(pmign);
    compRat_poly_init(p_global);
    
    if (parse==1) {
        mignotte_generalized(pmign, degree,(ulong) power, bitsize);
        compRat_poly_set_realRat_poly(p_global,pmign);
    
//         ccluster_interface_func( getApprox, bInit, eps, st, verbosity);
        ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, output, verbosity);
    }
    
    realRat_poly_clear(pmign);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
