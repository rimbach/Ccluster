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

#include <string.h>
#include <stdio.h>
#include "polynomials/compRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"

#include "./parseArgs.h"

compRat_poly_t p_global;
int d_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

void evaluateRunnelsFast( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec ){
    
    compApp_one(dest);
    compApp_zero(dest2);
    
    compApp_t destm1, destm2;
    compApp_init(destm1);
    compApp_init(destm2);
    
    compApp_t dest2m1, dest2m2, temp;
    compApp_init(dest2m1);
    compApp_init(dest2m2);
    compApp_init(temp);
    
    compApp_one(destm2);
    compApp_set(destm1, point);
    
    compApp_zero(dest2m2);
    compApp_one(dest2m1);
    
    
    for (int i = 2; i <= d_global; i++){
        
        /*compute value of derivative*/
        compApp_pow_si(temp, destm1, 1, prec);
        compApp_mul_si(temp, temp, 2, prec);
        compApp_mul( dest2, temp, dest2m1, prec); /* r*q_k(z)^{r-1}*q_k'(z) */
        
        compApp_pow_si(temp, destm2, 3, prec);
        compApp_mul_si(temp, temp, 4, prec);
        compApp_mul(temp, temp, point, prec);
        compApp_mul(temp, temp, dest2m2, prec); /* z*r^2*q_{k-1}(z)^{r^2-1}*q_{k-1}'(z) */
        compApp_pow_si(dest2m2, destm2,4,prec); /* q_{k-1}(z)^{r^2} */
        
        compApp_add(dest2, dest2, dest2m2, prec);
        compApp_add(dest2, dest2, temp, prec);
        
        compApp_set(dest2m2, dest2m1);
        compApp_set(dest2m1, dest2);
        
        
        /*compute value of poly*/
        compApp_pow_si(destm2, destm2, 4, prec);
        compApp_mul(destm2, destm2, point, prec);
        
        compApp_pow_si(dest, destm1, 2, prec);
        compApp_add(dest, dest, destm2, prec);
        
        compApp_set(destm2, destm1);
        compApp_set(destm1, dest);
    }
    
    compApp_clear(destm1);
    compApp_clear(destm2);
    
    compApp_clear(dest2m1);
    compApp_clear(dest2m2);
    compApp_clear(temp);
    
    return;
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
    
    realRat_poly_t prun, prunm1, prunm2, pone, px;
    realRat_poly_init(prun);
    realRat_poly_init(prunm1);
    realRat_poly_init(prunm2);
    realRat_poly_init(pone);
    realRat_poly_init(px);
    
    compRat_poly_init(p_global);
    d_global = degree;
        
    if (parse==1) {
        
        realRat_poly_zero(px);
        realRat_poly_set_coeff_si_ui(px, 1, 1, 1);
        
        realRat_poly_one(prunm2);
        realRat_poly_set(prunm1, px);
        realRat_poly_one(prun);
        
        for (int i = 2; i<=degree; i++) {
            
            realRat_poly_pow(prunm2, prunm2, 4);
            realRat_poly_mul(prunm2, prunm2, px);
            
            realRat_poly_pow(prun, prunm1, 2);
            realRat_poly_add(prun, prun, prunm2);
            
            realRat_poly_set(prunm2, prunm1);
            realRat_poly_set(prunm1, prun);
        }
        
        compRat_poly_set_realRat_poly(p_global,prun);
        
//         ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, verbosity);
        ccluster_interface_func_eval( getApprox, evaluateRunnelsFast, bInit, eps, st, nbthreads, verbosity);
    }
    
    realRat_poly_clear(prun);
    realRat_poly_clear(prunm1);
    realRat_poly_clear(prunm2);
    realRat_poly_clear(pone);
    realRat_poly_clear(px);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
