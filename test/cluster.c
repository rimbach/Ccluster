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

void clustersIterate( compApp_poly_ptr tabres, compApp_poly_ptr tabprec, int i, slong prec){
    // tabres is a table of 3^i compApp_poly
    // tabres is a table of 3^(i-1) compApp_poly
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
    int indexInTabRes = 0;
//     printf("pow(3,i-1): %d\n", (int) pow(3,i-1));
    for (int j = 0; j<((int) pow(3,i-1)); j++){
//         printf("%d\n", j);
//         realRat_set_si(modu, -1, (ulong) 0x1<<(4*(i-1)));
        realRat_set_si(modu, -1, (ulong) pow(4, 2*(i-1)));
        compApp_set_realRat( a_modu, modu, prec);
        
        realRat_set_si(argu, 2, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 4, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 6, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
    }
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_clear(coeff);    
}

void getApprox(compApp_poly_t dest, slong prec){
    
    compApp_poly_ptr tabprec;
    tabprec = (compApp_poly_ptr) malloc (1*sizeof(compApp_poly));
    compApp_poly_init2( tabprec, 2);
    compApp_poly_zero(tabprec);
    compApp_poly_set_coeff_si(tabprec, 1, 1);
    int degree = 3;
    
    for (int i=1; i<=p_degree; i++) {
        compApp_poly_ptr tabres;
        tabres = (compApp_poly_ptr) malloc (degree*sizeof(compApp_poly));
        for(int j = 0; j<degree; j++) compApp_poly_init2( tabres + j, 2);
        
        clustersIterate( tabres, tabprec, i, prec);
        for(int j = 0; j<((int) (degree/3)); j++) compApp_poly_clear( tabprec + j);
        free(tabprec);
        tabprec = tabres;
        degree = degree*3;
    }
    
    degree = (int) degree/3;
//     compApp_poly_t dest;
//     compApp_poly_init2(dest, degree+1);
    compApp_poly_one(dest);
    for(int j = 0; j<degree; j++) compApp_poly_mul(dest, dest, tabprec +j, prec);
    for(int j = 0; j<degree; j++) compApp_poly_clear( tabprec + j);
    free(tabprec);
}

int main(int argc, char **argv){
    
    if (argc<6){
        printf("usage: %s degree ", argv[0]);
        printf("initial_box epsilon strategy verbosity\n");
        printf("initial_box: for instance 0,1,1,2,100,1 i.e. the square centered in 0/1 +i*(1/2) of width 100/1\n");
        printf("eps:         for instance1,100 (1/100) or -53 (1/2^(-53))\n");  
        printf("strategy:    7: newton iterations, prediction of precision, optimized counting and discarding tests\n");
        printf("             15: same as 7 + eps is used only as an escape bound\n");
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
    
    p_degree = (slong) degree;
    
    if (parse==1)
        ccluster_interface_func( getApprox, bInit, eps, st, verbosity);
    
    realRat_clear(eps);
    compBox_clear(bInit);
    
    return 0;
}
