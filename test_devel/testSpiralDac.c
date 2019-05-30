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
#include <time.h>
#include "polynomials/compRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"
#include "ccluster/cclusterDAC.h"
#include "numbers/compApp.h"

#include "../test/parseArgs.h"

// compRat_poly_t p_global;
slong p_degree;

int checkAccuracy( compApp_poly_t p, slong prec ){
    for(int i=0; i<=p_degree; i++){
//         printf("degree: %d, accuracy: %d, prec: %d \n", i, (int) acb_rel_error_bits( compApp_poly_getCoeff(p, i)), (int) prec );
        if ( -acb_rel_error_bits( compApp_poly_getCoeff(p, i) ) < prec )
            return 0;
    }
    return 1;
}

void getApprox(compApp_poly_t dest, slong prec){
    
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

void getApprox2(compApp_poly_t dest, slong prec){
    
    slong prectemp = 2*prec;
    getApprox( dest, prectemp );
    
    while (!checkAccuracy( dest, prec)){
        prectemp = 2*prectemp;
        getApprox( dest, prectemp );
    }
    
    acb_poly_set_round( dest, dest, prec);
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
    nbthreads = (nbthreads<<6);
    int add_temp = (st>>7)<<17;
    st = st&((0x1<<7)-1);
    st = st + nbthreads + add_temp;
    
    p_degree = (slong) degree;
    
    clock_t tic = clock();
    
    if (parse==1)
        ccluster_interface_func( getApprox2, bInit, eps, st, verbosity);
    
    double ellapsedBFS = ((double) (clock() - tic))/ CLOCKS_PER_SEC;
    
    int paquets = (int) degree/8;
    
    connCmp_list_t qResults;
    connCmp_list_t qAllResults;
    connCmp_list_t qMainLoop;
    connCmp_list_t discardedCcs;
    
    connCmp_list_init( qResults    );
    connCmp_list_init( qAllResults    );
    connCmp_list_init( qMainLoop    );
    connCmp_list_init( discardedCcs    );
    
    tic = clock();
    
    ccluster_DAC_first_interface_forJulia( qResults, qAllResults, qMainLoop, discardedCcs,
                                            getApprox2, paquets, bInit, eps, 23, verbosity);
    while (!connCmp_list_is_empty(qResults) )
        connCmp_list_push(qAllResults, connCmp_list_pop(qResults));
    
    while (!connCmp_list_is_empty(qMainLoop) ) {
        ccluster_DAC_next_interface_forJulia( qResults, qAllResults, qMainLoop, discardedCcs,
                                            getApprox2, paquets, bInit, eps, 23, verbosity);
        while (!connCmp_list_is_empty(qResults) )
            connCmp_list_push(qAllResults, connCmp_list_pop(qResults));
    }
        
    double ellapsedDFS = ((double) (clock() - tic))/ CLOCKS_PER_SEC;
    printf("time BFS: %f\n", ellapsedBFS);
    printf("time DFS: %f\n", ellapsedDFS);
    
        
    connCmp_list_clear( qResults    );
    connCmp_list_clear( qAllResults    );
    connCmp_list_clear( qMainLoop    );
    connCmp_list_clear( discardedCcs    );
    
    
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
