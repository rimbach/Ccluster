/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpq_mpoly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include "polynomials/realRat_mpoly.h"
#include <time.h>

int realRat_poly_evaluate_all_but_one_fmpq( realRat_poly_t ev, const fmpq_mpoly_t A, fmpq *const *vals, const fmpq_mpoly_ctx_t ctx ){
    
    fmpq_mpoly_univar_t p_uni;
    fmpq_mpoly_univar_init( p_uni, ctx );
    slong nvars = fmpq_mpoly_ctx_nvars(ctx);
    
    fmpq_mpoly_to_univar(p_uni, A, nvars-1, ctx);
    
    slong l = fmpq_mpoly_univar_length(p_uni, ctx);
    slong d = fmpq_mpoly_univar_get_term_exp_si(p_uni, 0, ctx);
    
    realRat_poly_zero(ev);
    realRat_poly_fit_length(ev, d+1);
    ev->length = d+1;
    
    fmpq_t temp;
    fmpq_init(temp);
        
    fmpq_mpoly_t ptemp;
    fmpq_mpoly_init(ptemp, ctx);
    
    for (slong i=0; i<l; i++) {
        
        fmpq_mpoly_univar_get_term_coeff(ptemp, p_uni, i, ctx);
        slong deg = fmpq_mpoly_univar_get_term_exp_si(p_uni, i, ctx);
        fmpq_mpoly_evaluate_all_fmpq( temp, ptemp, vals, ctx );
        realRat_poly_set_coeff_realRat(ev, deg, temp);
    }
    
    fmpq_mpoly_clear(ptemp, ctx);
    fmpq_clear(temp);
    fmpq_mpoly_univar_clear( p_uni, ctx );
    
    return 1;
}

int main() {
    
    realRat_mpoly_t A;
    realRat_mpoly_init(A);
    
    realRat_mpoly_print(A);
    printf("\n");
    
    realRat_mpoly_clear(A);
    
    fmpq_mpoly_ctx_t context;
    fmpq_mpoly_t p;
    
    /* test univariate */
    const char *vars1[1] = {"x"};
    fmpq_mpoly_ctx_init( context, 1, ORD_LEX );
    fmpq_mpoly_init(p, context);
    
    char *pstr1 = "2*x^3 + 3*x^2 + 5"; 
    int res = fmpq_mpoly_set_str_pretty(p, pstr1, vars1, context);
    
    printf("successful parsing: %d\n", res);
    fmpq_mpoly_print_pretty(p, vars1, context);
    printf("\n\n");
    
    realRat_mpoly_init_fmpq_mpoly( A, p, context );
    realRat_mpoly_print(A);
    printf("\n");
    
    realRat_mpoly_clear(A);
    
    fmpq_mpoly_clear(p, context);
    fmpq_mpoly_ctx_clear( context );
    
    /* test bivariate */
    const char *vars2[2] = {"x", "y"};
    fmpq_mpoly_ctx_init( context, 2, ORD_LEX );
    fmpq_mpoly_init(p, context);
    
    char *pstr2 = "2*x^2*y^4 + 3*x*y^2 + x^2 + 3"; 
    int res2 = fmpq_mpoly_set_str_pretty(p, pstr2, vars2, context);
    
    printf("successful parsing: %d\n", res2);
    fmpq_mpoly_print_pretty(p, vars2, context);
    printf("\n\n");
    
    realRat_mpoly_init_fmpq_mpoly( A, p, context );
    realRat_mpoly_print(A);
    printf("\n");
    
    realRat_mpoly_clear(A);
    
    fmpq_mpoly_clear(p, context);
    fmpq_mpoly_ctx_clear( context );
    
        /* test trivariate */
    const char *vars3[3] = {"x", "y", "z"};
    fmpq_mpoly_ctx_init( context, 3, ORD_LEX );
    fmpq_mpoly_init(p, context);
    
//     char *pstr3 = "2*x^2*y^4 + 3*x*y^2 + x^2 + 3"; 
    char *pstr3 = "(1/27)*(2*x^2*y^4*z^2 + 3*x*y^3*z + 3*x*y^2 + x^2 + 3)";
    int res3 = fmpq_mpoly_set_str_pretty(p, pstr3, vars3, context);
    
    printf("successful parsing: %d\n", res3);
    fmpq_mpoly_print_pretty(p, vars3, context);
    printf("\n\n");
    
    realRat_mpoly_init_fmpq_mpoly( A, p, context );
    realRat_mpoly_print(A);
    printf("\n");
    printf("\n");
    
    realRat_mpoly_clear(A);
    
    fmpq_mpoly_clear(p, context);
    fmpq_mpoly_ctx_clear( context );
    
    /* test evaluate */
//     slong nbvars = 1;
//     ulong degrees[1] = {300};
//     slong nbvars = 2;
//     ulong degrees[2] = {100, 2};
//     slong nbvars = 3;
//     ulong degrees[3] = {100, 5, 4};
    slong nbvars = 4;
    ulong degrees[4] = {100, 10, 10, 10};
//     slong nbvars = 5;
//     ulong degrees[5] = {150,2,3,4,5};

    slong bitsize = 20;
    slong nbterms = (degrees[0]+1);
    for (int i = 1; i<nbvars; i++)
        nbterms = nbterms*(degrees[i]+1);
//     slong nbterms = 10;
    slong nbtests = 100;
    slong prec = 530;
    
    clock_t ti;
    double timeIn1=0;
    double timeIn2=0;
    
    fmpq_mpoly_ctx_init( context, nbvars, ORD_LEX );
    
    fmpq_mpoly_t prand;
    fmpq_mpoly_init(prand, context);
    
    flint_rand_t rt;
    flint_randinit( rt );
    fmpq_t val, val2;
    fmpq_init(val);
    fmpq_init(val2);
    
    compApp_t val3;
    compApp_init(val3);
    
    realRat_poly_t evpolyexact;
    realRat_poly_init(evpolyexact);
    compApp_poly_t evpolyapprox;
    compApp_poly_init( evpolyapprox);
    
    fmpq **vec2;
    vec2 = (fmpq **) flint_malloc ( nbvars * sizeof(fmpq) );
    for (int var = 0; var<nbvars; var++){
        vec2[var] = (fmpq *) flint_malloc ( sizeof(fmpq_t) );
        fmpq_init( vec2[var] );
        fmpq_randtest_not_zero( vec2[var], rt, 8);
    }
    
    compApp **vec3;
    vec3 = (compApp **) flint_malloc ( nbvars * sizeof(compApp) );
    for (int var = 0; var<nbvars; var++){
        vec3[var] = (compApp *) flint_malloc ( sizeof(compApp_t) );
        compApp_init( vec3[var] );
        compApp_set_realRat( vec3[var], vec2[var], prec );
    }
    
    slong * degsAsGet = (slong *) ccluster_malloc ( nbvars*sizeof(slong) );
    
    for (int test = 0; test<nbtests; test++){
            fmpq_mpoly_randtest_bounds(prand, rt, nbterms, bitsize, degrees, context);
            realRat_mpoly_init_fmpq_mpoly( A, prand, context );
//             fmpq_mpoly_print_pretty(prand, NULL, context);
//             printf("\n length: %ld\n", fmpq_mpoly_length(prand, context));
//             printf("\n\n");
            
            fmpq_mpoly_degrees_si ( degsAsGet, prand, context );
            printf("degrees: [ ");
            for (slong i = 0; i<nbvars; i++ ){
                printf(" %ld, ", degsAsGet[i] );
                degsAsGet[i]=0;
            }
            printf("]\n\n");
            
            realRat_mpoly_getDegrees( degsAsGet, A );
            printf("degrees: [ ");
            for (slong i = 0; i<nbvars; i++ ){
                printf(" %ld, ", degsAsGet[i] );
                degsAsGet[i]=0;
            }
            printf("]\n\n");
            
            ti = clock();
            fmpq_mpoly_evaluate_all_fmpq( val, prand, vec2, context );
//             printf("value of evaluation: ");
// //             fmpq_print(val);
//             compApp_set_realRat(val3, val, prec);
//             compApp_printd(val3, 10);
//             printf("\n\n");
            ti = clock() - ti;
            timeIn1 += ((float)ti)/CLOCKS_PER_SEC;
            
            
            compApp_zero(val3);
            ti = clock();
            realRat_mpoly_eval_compApp( val3, A, vec3, prec );
//             printf("value of evaluation: ");
//             compApp_printd(val3, 10);
//             printf("\n\n");
            ti = clock() - ti;
            timeIn2 += ((float)ti)/CLOCKS_PER_SEC;
            
            realRat_poly_evaluate_all_but_one_fmpq( evpolyexact, prand, vec2, context);
            
//             printf("exact polynomial: \n");
//             realRat_poly_print( evpolyexact );
//             printf("\n\n");
//             
//             printf("approximated polynomial: \n");
//             compApp_poly_set_realRat_poly(evpolyapprox, evpolyexact, prec);
//             compApp_poly_printd( evpolyapprox, 10);
//             printf("\n\n");
//             
//             printf("approximated polynomial 2: \n");
//             compApp_poly_ptr res = realRat_mpoly_eval_m1_compApp_poly_ptr( A , vec3, prec );
//             compApp_poly_printd( res, 10);
//             printf("\n\n");
            
            realRat_mpoly_clear(A);
            
    }
    
    printf(" time library   function: %f\n", timeIn1);
    printf(" time home made function: %f\n", timeIn2);
    
    for (int var = 0; var<nbvars; var++){
        fmpq_clear( vec2[var] );
        flint_free( vec2[var] );
    } 
    flint_free(vec2);
    
    for (int var = 0; var<nbvars; var++){
        compApp_clear( vec3[var] );
        flint_free( vec3[var] );
    } 
    flint_free(vec3);
    
    fmpq_clear(val);
    fmpq_clear(val2);
    compApp_clear(val3);
    flint_randclear(rt);
    
    fmpq_mpoly_clear(prand, context);
    fmpq_mpoly_ctx_clear( context );
    
    realRat_poly_clear(evpolyexact);
    compApp_poly_clear( evpolyapprox);
    
    ccluster_free(degsAsGet);
    
    flint_cleanup();
    
    return 0;
    
}
