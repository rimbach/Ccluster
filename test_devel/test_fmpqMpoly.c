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
#include <time.h>

double time_in_uni;

int my_fmpq_mpoly_evaluate_all_fmpq( fmpq_t ev, const fmpq_mpoly_t A, fmpq *const *vals, 
//                                      const char **x,
                                     const fmpq_mpoly_ctx_t ctx ){
    
    slong len = fmpq_mpoly_length(A, ctx);
    slong nvars = fmpq_mpoly_ctx_nvars(ctx);
//     printf("--evaluate: nb of terms of A: %ld\n", len);
//     printf("--number of variables in the context: %ld\n", nvars );
    
    fmpq_mpoly_t M;
    fmpq_mpoly_init(M, ctx);
    
//     ulong *exps;
//     exps = (ulong *) flint_malloc( nvars*sizeof(ulong) );
    
    fmpq_t coeff, temp;
    fmpq_init(coeff);
    fmpq_init(temp);
    
    fmpq_zero(ev);
    
    
    for( slong i = 0; i<len; i++){
//         fmpq_mpoly_get_term( M, A, i, ctx);
//         printf("--evaluate: %ldth monomial: ", i);
//         fmpq_mpoly_print_pretty(M, x, ctx);
//         printf("\n");
//         
//         fmpq_mpoly_get_term_exp_ui(exps, A, i, ctx);
//         printf("--                powers: ");
//         for (int j=0;j<nvars; j++) printf( "%s^%lu, ", x[j], exps[j] );
//         printf("\n");
//         
//         fmpq_mpoly_get_term_coeff_fmpq(coeff, A, i, ctx);
//         printf("--                coeff:  ");
//         fmpq_print(coeff);
//         printf("\n");
        
        fmpq_mpoly_get_term_coeff_fmpq(coeff, A, i, ctx);
        for (slong var=0; var<nvars; var++){
            fmpq_pow_si(temp, vals[var], (slong) fmpq_mpoly_get_term_var_exp_ui(A, i, var, ctx) );
            fmpq_mul(coeff, coeff, temp);
        }
        
        fmpq_add(ev, ev, coeff);
    }
    
    fmpq_clear(coeff);
    fmpq_clear(temp);
//     flint_free(exps);
    fmpq_mpoly_clear(M, ctx);
    
    return 1;
}

int my_fmpq_mpoly_evaluate_all_compApp( compApp_t ev, const fmpq_mpoly_t A, compApp *const *vals, 
//                                      const char **x,
                                     const fmpq_mpoly_ctx_t ctx, slong prec ){
    
    slong len = fmpq_mpoly_length(A, ctx);
    slong nvars = fmpq_mpoly_ctx_nvars(ctx);
    
    fmpq_mpoly_t M;
    fmpq_mpoly_init(M, ctx);
    
    fmpq_t coeff;
    compApp_t coeffApp, temp;
    fmpq_init(coeff);
    compApp_init(coeffApp);
    compApp_init(temp);
    
    compApp_zero(ev);
    
    
    for( slong i = 0; i<len; i++){
        
        fmpq_mpoly_get_term_coeff_fmpq(coeff, A, i, ctx);
        compApp_set_realRat( coeffApp, coeff, prec );
                             
        for (slong var=0; var<nvars; var++){
            compApp_pow_si(temp, vals[var], (slong) fmpq_mpoly_get_term_var_exp_ui(A, i, var, ctx), prec );
            compApp_mul(coeffApp, coeffApp, temp, prec);
        }
        
        compApp_add(ev, ev, coeffApp, prec);
    }
    
    fmpq_clear(coeff);
    compApp_clear(temp);
    compApp_clear(coeffApp);
//     flint_free(exps);
    fmpq_mpoly_clear(M, ctx);
    
    return 1;
}

int my_fmpq_mpoly_evaluate_all_compApp_Rec( compApp_t ev, const fmpq_mpoly_t A, compApp *const *vals, 
//                                               const char **x,
                                                 slong index,
                                                 const fmpq_mpoly_ctx_t ctx, slong prec ){
        
        fmpq_mpoly_univar_t p_uni;
        fmpq_mpoly_univar_init( p_uni, ctx );
        
        clock_t ti;
        ti = clock();
        
        fmpq_mpoly_to_univar(p_uni, A, index, ctx);
        
        ti = clock() - ti;
        time_in_uni += ((float)ti)/CLOCKS_PER_SEC;
            
        slong l = fmpq_mpoly_univar_length(p_uni, ctx);
        slong d = fmpq_mpoly_univar_get_term_exp_si(p_uni, 0, ctx);
        
        compApp_poly_t p;
        compApp_poly_init2(p, d+1);
        p->length = d+1;
        
        fmpq_t temp;
        fmpq_init(temp);
        
        fmpq_mpoly_t ptemp;
        fmpq_mpoly_init(ptemp, ctx);
        
        if (index == 0) {
            
//             printf("---index: %ld\n", index);
            
            for (slong i=0; i<l; i++) {
                
                fmpq_mpoly_univar_get_term_coeff(ptemp, p_uni, i, ctx);
//                 printf("--- index %ld --- %ld-th term:\n------",index,  i);
//                 fmpq_mpoly_print_pretty(ptemp, NULL, ctx);
//                 printf("\n");
                
                fmpq_mpoly_get_term_coeff_fmpq(temp, ptemp, 0, ctx);
                slong deg = fmpq_mpoly_univar_get_term_exp_si(p_uni, i, ctx);
                compApp_set_realRat( p->coeffs + deg, temp, prec );
//                 printf("------degree: %ld\n", deg); 
            }
//             printf("------univariate poly: \n------");
//             compApp_poly_printd(p, 5);
//             printf("\n");
            compApp_poly_evaluate( ev, p, vals[index], prec );
                                         
        } else {
            
//             printf("---index: %ld\n", index);
            
            for (slong i=0; i<l; i++) {
                
                fmpq_mpoly_univar_get_term_coeff(ptemp, p_uni, i, ctx);
                
//                 printf("--- index %ld --- %ld-th term:\n------",index,  i);
//                 fmpq_mpoly_print_pretty(ptemp, NULL, ctx);
//                 printf("\n");
                
                slong deg = fmpq_mpoly_univar_get_term_exp_si(p_uni, i, ctx);
                
                my_fmpq_mpoly_evaluate_all_compApp_Rec( p->coeffs + deg, ptemp, vals, 
//                                               x,
                                                 index -1,
                                                 ctx, prec );
                
                compApp_poly_evaluate( ev, p, vals[index], prec );
            }
            
        }
        
        fmpq_mpoly_clear(ptemp, ctx);
        fmpq_clear(temp);
        compApp_poly_clear(p);
        fmpq_mpoly_univar_clear( p_uni, ctx );
        
        return 1;
}

int main() {
    
    const char *vars[2] = {"x","y"};
    
    fmpq_mpoly_ctx_t context;
    fmpq_mpoly_ctx_init( context, 2, ORD_LEX );
    
    fmpq_mpoly_t p;
    fmpq_mpoly_init(p, context);
    
    char *pstr = "2*x^2*y^4 + 3*x*y^2 + x^2 + 3"; 
    int res = fmpq_mpoly_set_str_pretty(p, pstr, vars, context);
    
    printf("successful parsing: %d\n", res);
    fmpq_mpoly_print_pretty(p, vars, context);
    printf("\n\n");
    
    pstr = "x^2*(2*y^4 + 1) + 3*x*y^2 + 3"; 
    res = fmpq_mpoly_set_str_pretty(p, pstr, vars, context);
    printf("successful parsing: %d\n", res);
    
    fmpq_mpoly_print_pretty(p, vars, context);
    printf("\n\n");
    
    fmpq_t a1, a2, ev;
    fmpq_init(a1);
    fmpq_init(a2);
    fmpq_init(ev);
    
    fmpq_set_si(a1,1,2);
    fmpq_set_si(a2,1,3);
    
    fmpq* vec[2] = {a1,a2};
    fmpq_mpoly_evaluate_all_fmpq( ev, p, vec, context );
    printf("value of evaluation: ");
    fmpq_print(ev);
    printf("\n\n");
    
    fmpq_zero(ev);
    
//     my_fmpq_mpoly_evaluate_all_fmpq( ev, p, vec, vars, context );
    my_fmpq_mpoly_evaluate_all_fmpq( ev, p, vec, context );
    printf("value of evaluation: ");
    fmpq_print(ev);
    
    printf("\n");
    
    fmpq_clear(a1);
    fmpq_clear(a2);
    fmpq_clear(ev);
    
    fmpq_mpoly_clear(p, context);
    fmpq_mpoly_ctx_clear( context );
    
    
    printf("\n\n");
    
    slong nbvars = 5;
//     slong nbvars = 1;
    ulong degrees[5] = {150,2,2,2,2};
//     ulong degrees[1] = {300};
    slong bitsize = 20;
    slong nbterms = (degrees[0]+1);
    for (int i = 1; i<nbvars; i++)
        nbterms = nbterms*(degrees[i]+1);
//     slong nbterms = 100000;
    slong nbtests = 100;
    
    slong prec = 530;
    
    
    clock_t ti;
    double timeIn1=0;
    double timeIn2=0;
    double timeIn3=0;
    double timeIn4=0;
    
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
    
    time_in_uni = 0;
    
    for (int test = 0; test<nbtests; test++){
            fmpq_mpoly_randtest_bounds(prand, rt, nbterms, bitsize, degrees, context);
//          fmpq_mpoly_print_pretty(prand, NULL, context);
//             printf("\n length: %ld\n", fmpq_mpoly_length(prand, context));
//          printf("\n\n");
            
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
            my_fmpq_mpoly_evaluate_all_fmpq( val2, prand, vec2, context );
//             printf("value of evaluation: ");
// //             fmpq_print(val2);
//             compApp_set_realRat(val3, val, prec);
//             compApp_printd(val3, 10);
//             printf("\n\n");
            ti = clock() - ti;
            timeIn2 += ((float)ti)/CLOCKS_PER_SEC;
            
            compApp_zero(val3);
            ti = clock();
            my_fmpq_mpoly_evaluate_all_compApp( val3, prand, vec3, context, prec );
//             printf("value of evaluation: ");
//             compApp_printd(val3, 10);
//             printf("\n\n");
            ti = clock() - ti;
            timeIn3 += ((float)ti)/CLOCKS_PER_SEC;
            
            ti = clock();
            my_fmpq_mpoly_evaluate_all_compApp_Rec( val3, prand, vec3, nbvars-1, context, prec );
//             printf("value of evaluation: ");
//             compApp_printd(val3, 10);
//             printf("\n\n");
            ti = clock() - ti;
            timeIn4 += ((float)ti)/CLOCKS_PER_SEC;
            
    }
    
    printf(" time library   function: %f\n", timeIn1);
    printf(" time home made function: %f\n", timeIn2);
    printf(" time home made function: %f\n", timeIn3);
    printf(" time home made function: %f\n", timeIn4);
    printf(" time in converting     : %f\n", time_in_uni);
    
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
    
    return 0;
}
