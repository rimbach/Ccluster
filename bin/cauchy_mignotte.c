/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
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
int bitsize_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

void Mignotte_polynomial( realRat_poly_t dest, slong degree, slong bitsize);
void Mignotte_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec );

int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s degree bitsize [OPTIONS] ", argv[0]);
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
        return -1;
    }
    
    if (argc<=3){ /* display usage */
        printf("usage: %s degree bitsize [OPTIONS]\n", argv[0]);
        printf("   or: %s to see options\n", argv[0]);
    }
    
    int parse = 1;
    int degree;
    int bitsize;
    char * st;
    int verbosity = 1;
    int nbthreads = 1;
    int global = 2; /* by default, search all the roots */
    int infinity = 2;
    int output = 0;
    
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
    parse = parse*scan_degree( argv[2], &bitsize);
    
    /* loop on arguments to figure out options */
    for (int arg = 3; arg< argc; arg++) {
        
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
        
    }
    
    realRat_poly_t pmig;
    realRat_poly_init(pmig);
    
    compRat_poly_init(p_global);
    d_global = degree;
    bitsize_global = bitsize;
    
    Mignotte_polynomial( pmig, d_global, bitsize_global);
    compRat_poly_set_realRat_poly(p_global,pmig);
    
    if (parse==1) {
        cauchy_global_interface_func_eval( getApprox, Mignotte_evaluate, eps, isoRatio, nbPows, certified, 0, st, nbthreads, output, verbosity);
    }
    
    realRat_poly_clear(pmig);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    realRat_clear(isoRatio);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}

void Mignotte_polynomial( realRat_poly_t dest, slong degree, slong bitsize){
    
    mignotte_polynomial(dest, degree, bitsize);   
}

void Mignotte_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec ){
    
    compApp_t temp, temp2;
    compApp_init(temp);
    compApp_init(temp2);
    compApp_mul_2exp_si(temp, point, bitsize_global);
    /* temp = (x*2^bitsize -1) */
    compApp_add_si( temp, temp, -1, prec);
    
    /* temp2 = (x*2^bitsize -1)*2^(bitsize+2) */
    compApp_mul_2exp_si(temp2, temp, bitsize_global+2);
    /* temp = 2(x*2^bitsize -1)^2 */
    compApp_pow_si(temp, temp, 2, prec);
    compApp_mul_2exp_si(temp, temp, 1);
    
    /* dest2 = x^degree-1 */
    compApp_pow_si( dest2, point, d_global-1, prec);
    /* dest =  x^degree */
    compApp_mul( dest, dest2, point, prec);
    /* dest = x^degree - 2(x*2^bitsize -1)^2 */
    compApp_sub( dest, dest, temp, prec);
    /* dest = degree*x^(degree-1) - (x*2^bitsize -1)*2^(bitsize+2) */
    compApp_mul_si( dest2, dest2, d_global, prec);
    compApp_sub( dest2, dest2, temp2, prec);
    
    
    compApp_clear(temp);
    compApp_clear(temp2);
    
    /* check */
//     printf("----------------------------\n");
//     printf("   pval:"); compApp_printd(dest, 10); printf("\n");
//     printf("pderval:"); compApp_printd(dest2, 10); printf("\n");
//     compApp_poly_t p;
//     compApp_poly_init(p);
//     getApprox(p, prec);
//     compApp_poly_evaluate2_rectangular(dest, dest2, p, point, prec);
//     compApp_poly_clear(p);
//     printf("----------------------------\n");
//     printf("   pval:"); compApp_printd(dest, 10); printf("\n");
//     printf("pderval:"); compApp_printd(dest2, 10); printf("\n");
//     printf("----------------------------\n");
    return;
}
