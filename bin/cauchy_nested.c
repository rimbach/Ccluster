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
// #include "polynomials/compRat_poly.h"
// #include "polynomials/compApp_poly.h"
// #include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"
#include "cauchy/cauchy.h"
// #include "numbers/compApp.h"

#include "./parseArgs.h"

/* iteration */
slong p_iteration;
int   p_c;
int   p_a;
/* tables for writing intermediate results */
compApp_ptr p_dests;
compApp_ptr p_dest2s;
/* p_iteration first powers of a */
realRat_ptr aPows;

/* writes in dest a prec-bit approximation of the spiral polynomial of degree p_degree*/
void getApprox_Nested(compApp_poly_t dest, slong prec);

void _Nested_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec, int iteration, int nbRoots, int relativeWidth );

void Nested_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec ){
    
//     printf("Nested Evaluate, begin; prec: %ld\n", prec);
//     slong prectemp = 2*prec;
//     _Nested_evaluate( dest, dest2, point, prectemp, p_iteration, p_c, p_a );
//     printf("dest : "); compApp_printd(dest, 20); printf("\n"); 
//     printf("dest2: "); compApp_printd(dest2, 20); printf("\n"); 
//     while ( (!compApp_check_relOne_accuracy( dest, prec))||(!compApp_check_relOne_accuracy( dest2, prec)) ) {
//         prectemp = 2*prectemp;
//         _Nested_evaluate( dest, dest2, point, prectemp, p_iteration, p_c, p_a );
//         printf("---prectemp: %ld\n", prectemp);
//         printf("---dest : "); compApp_printd(dest, 20); printf("\n"); 
//         printf("---dest2: "); compApp_printd(dest2, 20); printf("\n");
//     }
//     compApp_set_round( dest, dest, prec);
//     compApp_set_round( dest2, dest2, prec);
    
    _Nested_evaluate( dest, dest2, point, prec, p_iteration, p_c, p_a );
}


int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s iteration [OPTIONS]", argv[0]);
        printf("                                 \n");
//         printf("      -d , --domain: the initial region of interest\n");
//         printf("                     global [default] finds all the roots\n");
//         printf("                     a box, for instance 0,1/2,100 i.e. the square centered in 0 +i*(1/2) of width 100\n");
//         printf("                     if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("      -c, --certified: certify the output\n");
        printf("                       1 [default] output is certified\n");
        printf("                       0 output is not certified; whp correct!\n");
        printf("      -e , --epsilon: the size of output clusters\n");
        printf("                     +inf [default] output natural clusters wits less roots than degree of input polynomial\n");
        printf("                                                                 or size less than 2^(-53)\n");
        printf("                     a positive rational as 1 or 1/100 or a negative power of 2 as -53 for 2^(-53)\n");
        printf("      -o , --output: the way cluster are output; default is NO OUTPUT\n");
        printf("                     0: [default] NO OUTPUT\n");
        printf("                     d>0: d digit precision floating point numbers\n");
        printf("                     -1: rational numbers\n");
        printf("                     -2 or g or G: gnuplot output: can be piped to gnuplot \n");
        printf("                     -3 or gs or GS: gnuplot output with subdivision tree \n");
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
        printf("       -c : nb of roots in the clusters, \n");
        printf("                       3: [default] or a positive integer            \n");
        printf("       -a : relative width of nested clusters, \n");
        printf("                       16: [default] or an integer >=2            \n");
        return -1;
    }
    
    if (argc<=2){ /* display usage */
        printf("usage: %s iteration [OPTIONS]\n", argv[0]);
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
    int c = 3;
    int a = 16;
    
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
        
        if ( (strcmp( argv[arg], "-c" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &c);
                if (c<=0){
                    printf("%s ERROR: NON-VALID c (should be >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-a" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &a);
                if (a<=1){
                    printf("%s ERROR: NON-VALID a (should be >=2) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
    }
    
    p_iteration = (slong) degree;
    p_c = c;
    p_a = a;
    /* initialize tables for intermediate results */
    p_dests  = (compApp_ptr) ccluster_malloc ( p_iteration*p_c*sizeof(compApp) );
    p_dest2s = (compApp_ptr) ccluster_malloc ( p_iteration*p_c*sizeof(compApp) );
    for (slong i=0; i<p_iteration*p_c; i++) {
        compApp_init(p_dests+i);
        compApp_init(p_dest2s+i);
    }
    aPows    = (realRat_ptr) ccluster_malloc ( (p_iteration+2)*sizeof(realRat) );
    for (slong i=0; i<=(p_iteration+1); i++) {
        realRat_init( aPows+i );
        if (i==0)
            realRat_one( aPows + i );
        if (i==1)
            realRat_set_si( aPows + i, (slong) p_a, 1);
        if (i==2)
            realRat_pow_si( aPows + i, aPows + (i-1), (slong) p_c);
        if (i>=2)
            realRat_mul( aPows+i, aPows + (i-1), aPows + 2 );
    }
    
    if (parse==1) {
//        cauchy_global_interface_func( getApprox_Nested, eps, isoRatio, nbPows, st, nbthreads, output, verbosity); 
       cauchy_global_interface_func_eval( getApprox_Nested, Nested_evaluate, eps, isoRatio, nbPows, certified, 0, st, nbthreads, output, verbosity);
    }
    
    realRat_clear(eps);
    compBox_clear(bInit);
    
    for (slong i=0; i<p_iteration*p_c; i++) {
        compApp_clear(p_dests+i);
        compApp_clear(p_dest2s+i);
    }
    ccluster_free ( p_dests  );
    ccluster_free ( p_dest2s );
    for (slong i=0; i<=(p_iteration+1); i++) {
        realRat_clear( aPows+i );
    }
    ccluster_free ( aPows );
    
    flint_cleanup();
    
    return 0;
}

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

void getApprox_temp(compApp_poly_t dest, slong prec){
    
    compApp_poly_ptr tabprec;
    tabprec = (compApp_poly_ptr) malloc (1*sizeof(compApp_poly));
    compApp_poly_init2( tabprec, 2);
    compApp_poly_zero(tabprec);
    compApp_poly_set_coeff_si(tabprec, 1, 1);
    int degree = 3;
    
    for (int i=1; i<=p_iteration; i++) {
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
    compApp_poly_one(dest);
    for(int j = 0; j<degree; j++) compApp_poly_mul(dest, dest, tabprec +j, prec);
    for(int j = 0; j<degree; j++) compApp_poly_clear( tabprec + j);
    free(tabprec);
}

void getApprox_Nested(compApp_poly_t dest, slong prec){
    
    slong prectemp = 2*prec;
    getApprox_temp( dest, prectemp );
    
    while (!compApp_poly_check_relOne_accuracy( dest, prec)){
        prectemp = 2*prectemp;
        getApprox_temp( dest, prectemp );
    }
    
    compApp_poly_set_round( dest, dest, prec);
}

void _Nested_evaluate( compApp_t dest, compApp_t dest2, const compApp_t point, slong prec, int n, int c, int a ){
    
    if (n==0) {
        compApp_set(dest, point);
        compApp_one(dest2);
        return;
    }
    
//     prec = 2*prec;
    
    slong degree = 1;
    for (int i = 1; i<= n; i++)
        degree = degree*c;
    
    realRat_ptr coeff = aPows + 1;
//     realRat_init(coeff);
//     realRat_set_si(coeff, (slong) a, 1);
    
//     compApp_ptr dests  = (compApp_ptr) ccluster_malloc ( c*sizeof(compApp) );
//     compApp_ptr dest2s = (compApp_ptr) ccluster_malloc ( c*sizeof(compApp) );
//     for (int i=0; i<c; i++) {
//         compApp_init(dests+i);
//         compApp_init(dest2s+i);
//     }
    compApp_ptr dests  = (compApp_ptr) (p_dests  + (n-1)*c);
    compApp_ptr dest2s = (compApp_ptr) (p_dest2s + (n-1)*c);
    
    realRat_t argu;
    compApp_t zeta, npoint;
    realRat_init(argu);
    compApp_init(zeta);
    compApp_init(npoint);
    
    for (int i=0; i<c; i++) {
        realRat_set_si( argu, 2*i, c);
        compApp_set_realRat( zeta, argu, prec);
        compApp_exp_pi_i( zeta, zeta, prec);
        compApp_sub(npoint, point, zeta, prec);
        compApp_mul_realRat(npoint, npoint, coeff, prec);
        _Nested_evaluate( dests+i, dest2s+i, npoint, prec, n-1, c, a );
    }
    
    compApp_one(dest);
    for (int i=0; i<c; i++) {
        compApp_mul( dest, dest, dests+i, prec );
    }
    
    compApp_zero(dest2);
    for (int i=0; i<c; i++) {
        compApp_one(zeta);
        for (int j=0; j<c; j++) {
            if (j==i)
                compApp_mul( zeta, zeta, dest2s+j, prec );
            else
                compApp_mul( zeta, zeta, dests+j, prec );
        }
        compApp_add( dest2, dest2, zeta, prec);
    }
    compApp_mul_realRat( dest2, dest2, coeff, prec );
    
//     realRat_pow_si(coeff, coeff, degree);
    
    compApp_div_realRat(dest, dest, aPows + (n + 1), prec);
    compApp_div_realRat(dest2, dest2, aPows + (n + 1), prec);
    
//     for (int i=0; i<c; i++) {
//         compApp_clear(dests+i);
//         compApp_clear(dest2s+i);
//     }
//     ccluster_free(dests);
//     ccluster_free(dest2s);
    
    realRat_clear(argu);
    compApp_clear(zeta);
    compApp_clear(npoint);
    realRat_clear(coeff);
}
