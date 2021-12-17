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

#include "cauchy/cauchy.h"

#include "../bin/parseArgs.h"

compRat_poly_t p_global;
int d_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

void Runnels_polynomial( realRat_poly_t dest, int iterations);
void Runnels_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec );

int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s iteration [OPTIONS] ", argv[0]);
        printf("                                 \n");
//         printf("      -d , --domain: the initial region of interest\n");
//         printf("                     global [default] finds all the roots\n");
//         printf("                     a box, for instance 0,1/2,100 i.e. the square centered in 0 +i*(1/2) of width 100\n");
//         printf("                     if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("      -c, --certified: certify the output\n");
        printf("                       1 [default] output is certified\n");
        printf("                       0 output is not certified; whp correct!\n");
        printf("      -e , --epsilon: the size of output clusters\n");
        printf("                     +inf [default] output natural clusters with less roots than degree of input polynomial\n");
        printf("                                                                 or size less than 2^(-53)\n");
        printf("                     a positive rational as 1 or 1/100 or a negative power of 2 as -53 for 2^(-53)\n");
        printf("      -o , --output: the way cluster are output; default is NO OUTPUT\n");
        printf("                     0: [default] NO OUTPUT\n");
        printf("                     d>0: d digit precision floating point numbers\n");
        printf("                     -1: rational numbers\n");
        printf("                     -2 or g or G: gnuplot output: can be piped to gnuplot \n");
        printf("      -m, --mode: the version of the algorithm\n");
        printf("                     C1: subdivision without compression  \n");
        printf("                     C2 [default]: subdivision with compression  \n");
        printf("      -i, --isoRatio: the assumed isolation ratio\n");
        printf("                      a positive, >1, rational as 2 or 4/3 [default]\n");
        printf("      -n, --nbPows:  the number of power sums computed for uncertified exclusion\n");
        printf("                      a positive integer as 3 [default]\n");
        printf("      -v, --verbose: an integer for verbosity\n");
        printf("                     0: nothing\n");
        printf("                     1 [default]: abstract of input and output\n");
        printf("                     2: detailed reports concerning algorithm\n");
        printf("                     >=3: debugging mode\n");
        printf("      -f, --machineFloats: use or not machine numbers when possible\n");
        printf("                     0 [default] : no\n");
        printf("                     1: yes\n");
        return -1;
    }
    
    if (argc<=2){ /* display usage */
        printf("usage: %s iteration [OPTIONS]\n", argv[0]);
        printf("   or: %s to see options\n", argv[0]);
    }
    
    int parse = 1;
    int degree;
    char * st;
    int verbosity = 1;
    int nbthreads = 1;
    int global = 2; /* by default, search all the roots */
    int infinity = 2;
    int output = 0;
    int usefpri = 0;
    
    compBox_t bInit;
    realRat_t eps;
    realRat_t isoRatio;
    int nbPows = 3;
    int certified = 1;
    
    char stDefault[] = "C2"; 
    st = stDefault;
    
    compBox_init(bInit);
    realRat_init(eps);
    realRat_init(isoRatio);
    scan_epsilon( "+inf", eps );
    scan_epsilon( "4/3", isoRatio );
    global = scan_initialBox( "global", bInit );
    
    parse = parse*scan_degree( argv[1], &degree);
    
    /* loop on arguments to figure out options */
    for (int arg = 2; arg< argc; arg++) {
        
        if ( (strcmp( argv[arg], "-c" ) == 0) || (strcmp( argv[arg], "--certified" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &certified);
                certified = (certified == 0? 0 : 1);
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-v" ) == 0) || (strcmp( argv[arg], "--verbose" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_verbosity(argv[arg+1], &verbosity );
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-d" ) == 0) || (strcmp( argv[arg], "--domain" ) == 0) ) {
            if (argc>arg+1) {
                global = scan_initialBox( argv[arg+1], bInit );
                if (global==1) {
                    printf("option --domain %s is not valid; doing global\n", argv[arg+1]);
                    global =2;
                }
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
        
        if ( (strcmp( argv[arg], "-m" ) == 0) || (strcmp( argv[arg], "--mode" ) == 0) ) {
            if (argc>arg+1) {
//                 parse = parse*scan_strategy( argv[arg+1], st );
                st = argv[arg+1];
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-i" ) == 0) || (strcmp( argv[arg], "--isoRatio" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_isoRatio( argv[arg+1], isoRatio );
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-n" ) == 0) || (strcmp( argv[arg], "--nbPows" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_nbPows( argv[arg+1], &nbPows );
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-f" ) == 0) || (strcmp( argv[arg], "--machineFloats" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &usefpri);
                if ((usefpri <0) || (usefpri>1)) {
                    usefpri=0;
                }
                arg++;
            }
        }
        
    }
    
    realRat_poly_t prunn;
    realRat_poly_init(prunn);
    
    compRat_poly_init(p_global);
    d_global = degree;
    
    Runnels_polynomial( prunn, d_global);
    compRat_poly_set_realRat_poly(p_global,prunn);
    
    if (parse==1) {
        cauchy_global_interface_func_eval( getApprox, Runnels_evaluate, eps, isoRatio, nbPows, certified, usefpri, st, nbthreads, output, verbosity);
    }
    
    realRat_poly_clear(prunn);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    realRat_clear(isoRatio);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}

void Runnels_polynomial( realRat_poly_t prun, int iterations){
    
    realRat_poly_t prunm1, prunm2, pone, px;
    realRat_poly_init(prunm1);
    realRat_poly_init(prunm2);
    realRat_poly_init(pone);
    realRat_poly_init(px);
    
    realRat_poly_zero(px);
    realRat_poly_set_coeff_si_ui(px, 1, 1, 1);
    
    realRat_poly_one(prunm2);
    realRat_poly_set(prunm1, px);
    realRat_poly_one(prun);
    
    for (int i = 2; i<=iterations; i++) {
        
        realRat_poly_pow(prunm2, prunm2, 4);
        realRat_poly_mul(prunm2, prunm2, px);
        
        realRat_poly_pow(prun, prunm1, 2);
        realRat_poly_add(prun, prun, prunm2);
        
        realRat_poly_set(prunm2, prunm1);
        realRat_poly_set(prunm1, prun);
    }
     
    realRat_poly_clear(prunm1);
    realRat_poly_clear(prunm2);
    realRat_poly_clear(pone);
    realRat_poly_clear(px); 
}

void Runnels_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec ){
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
