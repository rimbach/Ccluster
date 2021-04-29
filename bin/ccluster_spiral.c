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
// #include "polynomials/compRat_poly.h"
// #include "polynomials/compApp_poly.h"
// #include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"
// #include "numbers/compApp.h"

#include "./parseArgs.h"

/* degree of the spiral polynomial */
slong p_degree;
/* writes in dest a prec-bit approximation of the spiral polynomial of degree p_degree*/
void getApprox_Spiral(compApp_poly_t dest, slong prec);

int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s degree [OPTIONS]", argv[0]);
        printf("                                 \n");
        printf("      -d , --domain: the initial region of interest\n");
        printf("                     global [default] finds all the roots\n");
        printf("                     a box, for instance 0,1/2,100 i.e. the square centered in 0 +i*(1/2) of width 100\n");
        printf("                     if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("      -e , --epsilon: the size of output clusters\n");
        printf("                     +inf [default] output natural clusters wits less roots than degree of input polynomial\n");
        printf("                                                                 or size less than 2^(-53)\n");
        printf("                      a positive rational as 1 or 1/100 or a negative power of 2 as -53 for 2^(-53)\n");
        printf("      -o , --output: the way cluster are output; default is NO OUTPUT\n");
        printf("                     0: [default] NO OUTPUT\n");
        printf("                     d>0: d digit precision floating point numbers\n");
        printf("                     -1: rational numbers\n");
        printf("                     -2 or g or G: gnuplot output: can be piped to gnuplot \n");
        printf("                     -3 or gs or GS: gnuplot output with subdivision tree \n");
        printf("      -m, --mode: the version of the algorithm\n");
        printf("                     default value is \"default\"  \n");
        printf("      -v, --verbose: an integer for verbosity\n");
        printf("                     0: nothing\n");
        printf("                     1 [default]: abstract of input and output\n");
        printf("                     2: detailed reports concerning algorithm\n");
        printf("                     >=3: debugging mode\n");
        printf("      -j, --nbThreads: an positive integer for the number of threads\n");
        printf("                       1 [default]: one thread is used\n");
        printf("                       >1: no compatibility with -o -3 option\n");
        return -1;
    }
    
    if (argc<=2){ /* display usage */
        printf("usage: %s degree [OPTIONS]\n", argv[0]);
        printf("   or: %s to see options\n", argv[0]);
    }
    
    int parse = 1;
    int degree;
    char * st;
    int verbosity=1;
    int nbthreads = 1;
    int global = 2; /* by default, search all the roots */
    int infinity = 2;
    int output = 0;
    
    compBox_t bInit;
    realRat_t eps;
    
    char stDefault[] = "default"; 
    st = stDefault;
    
    compBox_init(bInit);
    realRat_init(eps);
    scan_epsilon( "+inf", eps );
    global = scan_initialBox( "global", bInit );
    parse = parse*scan_degree( argv[1], &degree);
    p_degree = (slong) degree;
    
    /* loop on arguments to figure out options */
    
    for (int arg = 2; arg< argc; arg++) {
        
        if ( (strcmp( argv[arg], "-v" ) == 0) || (strcmp( argv[arg], "--verbose" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_verbosity(argv[arg+1], &verbosity );
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-d" ) == 0) || (strcmp( argv[arg], "--domain" ) == 0) ) {
            if (argc>arg+1) {
                global = scan_initialBox( argv[arg+1], bInit );
                parse = parse*global;
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-e" ) == 0) || (strcmp( argv[arg], "--epsilon" ) == 0) ) {
            if (argc>arg+1) {
                infinity = scan_epsilon( argv[arg+1], eps );
                parse = parse*infinity;
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-o" ) == 0) || (strcmp( argv[arg], "--output" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_output(argv[arg+1], &output);
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-j" ) == 0)  ) {
            if (argc>arg+1) {
                parse = parse*scan_nbthreads(argv[arg+1], &nbthreads);
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-m" ) == 0) || (strcmp( argv[arg], "--mode" ) == 0) ) {
            if (argc>arg+1) {
                st = argv[arg+1];
                arg++;
            }
        }
        
    }
    
    
    
    if (parse) {
       if (global==2)
            ccluster_global_interface_func( getApprox_Spiral, eps, st, nbthreads, output, verbosity);
       else
            ccluster_interface_func( getApprox_Spiral, bInit, eps, st, nbthreads, output, verbosity);

        
    }
    
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}

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

void getApprox_Spiral(compApp_poly_t dest, slong prec){
    
    slong prectemp = 2*prec;
    getApprox_temp( dest, prectemp );
    
    while (!compApp_poly_check_relOne_accuracy( dest, prec)){
        prectemp = 2*prectemp;
        getApprox_temp( dest, prectemp );
    }
    
    compApp_poly_set_round( dest, dest, prec);
}
