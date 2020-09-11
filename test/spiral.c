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
#include "numbers/compApp.h"

#include "./parseArgs.h"

// compRat_poly_t p_global;
slong p_degree;

void getApprox_temp(compApp_poly_t dest, slong prec){
    
    realRat_t modu;
    realRat_t argu;
    compApp_t a_modu;
    compApp_t a_argu;
    compApp_t coeff;
    
    realRat_init(modu);
    realRat_init(argu);
    compApp_init(a_modu);
    compApp_init(a_argu);
    compApp_init(coeff);
    
    compApp_poly_t temp;
    compApp_poly_init2(temp,2);
    compApp_poly_set_coeff_si(temp, 1, 1);
    
    compApp_poly_one(dest);
    
    for(int i=1; i<=p_degree; i++){
        realRat_set_si(modu, -i, (ulong) p_degree);
        realRat_set_si(argu, 4*i, (ulong) p_degree);
        compApp_set_realRat( a_modu, modu, prec);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_poly_set_coeff_compApp(temp, 0, coeff);
        compApp_poly_mul(dest, dest, temp, prec);
    }
    
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_clear(coeff);
    compApp_poly_clear(temp);
}

void getApprox(compApp_poly_t dest, slong prec){
    
    slong prectemp = 2*prec;
    getApprox_temp( dest, prectemp );
    
    while (!compApp_poly_check_relOne_accuracy( dest, prec)){
        prectemp = 2*prectemp;
        getApprox_temp( dest, prectemp );
    }
    
    compApp_poly_set_round( dest, dest, prec);
}

int main(int argc, char **argv){
    
    if (argc<6){
        printf("usage: %s degree ", argv[0]);
        printf("initial_box epsilon strategy verbosity\n");
        printf("initial_box: global of a box,\n");
        printf("             for instance 0,1,1,2,100,1 i.e. the square centered in 0/1 +i*(1/2) of width 100/1\n");
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
    int global = 0;
    int output = 16;
    
    compBox_t bInit;
    realRat_t eps;
    
    compBox_init(bInit);
    realRat_init(eps);
    
    parse = parse*scan_degree( argv[1], &degree);
    global = scan_initialBox( argv[2], bInit );
    parse = parse*global;
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
    
    p_degree = (slong) degree;
    
    if (parse) {
    
        if (global==2)
            ccluster_global_interface_func( getApprox, eps, st, nbthreads, output, verbosity);
        else
            ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, output, verbosity);
    }
    
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
