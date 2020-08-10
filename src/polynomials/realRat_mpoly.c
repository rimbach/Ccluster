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

#include "polynomials/realRat_mpoly.h"

/* init to zero polynomial */
void realRat_mpoly_init( realRat_mpoly_t poly ) {
    realRat_mpoly_degref(poly) = 0;
    realRat_mpoly_dimref(poly) = 1;
    realRat_mpoly_lenref(poly) = 1;
    realRat_mpoly_degreesref(poly) = NULL;
    realRat_mpoly_monosref(poly) = NULL;
    realRat_poly_init(realRat_mpoly_monos_uniref(poly));
    compApp_poly_init(realRat_mpoly_forEvalref(poly));
}

void _realRat_mpoly_init_fmpq_mpoly( realRat_mpoly_t poly, slong nbvars, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx ){
    
//     const char *vars2[2] = {"x", "y"};
//     printf("--- init from fmpq_mpoly: nbvars: %ld\n", nbvars);
    /* get A as univariate polynomial */
    fmpq_mpoly_univar_t A_uni;
    fmpq_mpoly_univar_init( A_uni, ctx );
    fmpq_mpoly_to_univar(A_uni, A, nbvars-1, ctx);
       
    realRat_mpoly_degref(poly) = fmpq_mpoly_univar_get_term_exp_si(A_uni, 0, ctx);
    realRat_mpoly_dimref(poly) = nbvars;
    realRat_mpoly_lenref(poly) = fmpq_mpoly_univar_length(A_uni, ctx);
    
//     printf("--- init from fmpq_mpoly: deg: %ld, dim: %ld, len: %ld\n", poly->_deg, poly->_dim, poly->_len);
    
    fmpq_mpoly_t Acoeff;
    fmpq_mpoly_init(Acoeff, ctx);
    
    slong deg;
        
    if ( realRat_mpoly_dimref(poly) == 1 ){
        
        realRat_mpoly_degreesref(poly) = NULL;
        realRat_mpoly_monosref(poly) = NULL;
//         printf("--- init from fmpq_mpoly: dim %ld\n", realRat_mpoly_dimref(poly) );
        
        fmpq_t coeff;
        fmpq_init(coeff);
    
        realRat_poly_ptr p = realRat_mpoly_monos_uniref(poly);
        realRat_poly_init2( p , poly->_deg +1);
        p->length = realRat_mpoly_dimref(poly) + 1;
        
//         realRat_poly_print( realRat_mpoly_monos_uniref(poly) );
//         printf("\n");
        
        for (slong i=0; i<poly->_len; i++) {
            fmpq_mpoly_univar_get_term_coeff(Acoeff, A_uni, i, ctx);
            fmpq_mpoly_get_term_coeff_fmpq(coeff, Acoeff, 0, ctx);
            deg = fmpq_mpoly_univar_get_term_exp_si(A_uni, i, ctx);
            realRat_poly_set_coeff_realRat( p, deg, coeff );
        }
        
        fmpq_clear(coeff);
        
    } else {
        
        /* create and initialize degrees */
        realRat_mpoly_degreesref(poly) = (slong *) ccluster_malloc ( (realRat_mpoly_dimref(poly) - 1)*sizeof(slong) );
        for (slong v = 0; v < (realRat_mpoly_dimref(poly) - 1); v++)
            (realRat_mpoly_degreesref(poly))[v] = 0;
//         printf("--- init from fmpq_mpoly: dim %ld\n", realRat_mpoly_dimref(poly) );
        
        /* initialize table of monomials */
        realRat_mpoly_monosref(poly) = (realRat_mmono_ptr) ccluster_malloc ( realRat_mpoly_lenref(poly) * sizeof(struct realRat_mmono) );
        
        for (slong i=0; i < realRat_mpoly_lenref(poly); i++) {
            
//             printf("--- init from fmpq_mpoly: dim %ld, i: %ld\n", realRat_mpoly_dimref(poly), i );
            
            fmpq_mpoly_univar_get_term_coeff(Acoeff, A_uni, i, ctx);
            deg = fmpq_mpoly_univar_get_term_exp_si(A_uni, i, ctx);
            
//             printf("--- init from fmpq_mpoly: dim %ld, coeff: ", realRat_mpoly_dimref(poly));
//             fmpq_mpoly_print_pretty(Acoeff, vars2, ctx);
//             printf("\n");
            
            /* initialize monomial */
            realRat_mmono_ptr m = realRat_mpoly_monosref(poly) + i;
            realRat_mmono_dimref(m) = nbvars - 1;
            realRat_mmono_degref(m) = deg;
            realRat_mmono_coeffref(m) = (realRat_mpoly_ptr) ccluster_malloc(sizeof( struct realRat_mpoly));
            compApp_init( realRat_mmono_evalref(m) );
            /* recursive call to set coeff of monomial */
            _realRat_mpoly_init_fmpq_mpoly( realRat_mmono_coeffref(m), nbvars - 1, Acoeff, ctx );
//             printf("--- back to init from fmpq_mpoly: dim %ld\n", realRat_mpoly_dimref(poly) );
            
            /* actualize degrees */
            realRat_mpoly_ptr pm = realRat_mmono_coeffref(m);
            for (slong v = 1; v < (realRat_mpoly_dimref(poly) - 1); v++) {
                if ( (realRat_mpoly_degreesref( pm ))[v-1] > (realRat_mpoly_degreesref(poly))[v] )
                    (realRat_mpoly_degreesref(poly))[v] = (realRat_mpoly_degreesref( pm ))[v-1];
            }
            if ( realRat_mpoly_degref( pm ) > (realRat_mpoly_degreesref(poly))[0] )
                (realRat_mpoly_degreesref(poly))[0] = realRat_mpoly_degref( pm );
            
            compApp_init(realRat_mmono_evalref(m));
        }
        
    }
    
    compApp_poly_init2( realRat_mpoly_forEvalref(poly) , realRat_mpoly_dimref(poly) +1);
    
    fmpq_mpoly_clear(Acoeff, ctx);
    fmpq_mpoly_univar_clear( A_uni, ctx );
    
//     printf("--- end init from fmpq_mpoly: dim %ld\n", realRat_mpoly_dimref(poly) );
}

void realRat_mpoly_init_fmpq_mpoly( realRat_mpoly_t poly, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx ){
    slong nbvars = fmpq_mpoly_ctx_nvars(ctx);
    _realRat_mpoly_init_fmpq_mpoly( poly, nbvars, A, ctx );
}

void realRat_mpoly_eval_compApp( compApp_t ev, realRat_mpoly_t poly, compApp *const *vals, slong prec ){
    
    if (realRat_mpoly_dimref(poly)==1) {
        compApp_poly_set_realRat_poly(realRat_mpoly_forEvalref(poly), realRat_mpoly_monos_uniref(poly), prec);
    } else {
        /* compute univariate polynomial by recursive evaluation */
        for (slong i = 0; i<realRat_mpoly_lenref(poly); i++) {
            realRat_mmono_ptr m = realRat_mpoly_monosref(poly) + i;
            /* recursive call to evaluate the coeff */
            realRat_mpoly_eval_compApp( realRat_mmono_evalref(m), realRat_mmono_coeffref(m), vals, prec );
            compApp_poly_set_coeff_compApp( realRat_mpoly_forEvalref(poly), realRat_mmono_degref(m), realRat_mmono_evalref(m));
        }
    }
    
    compApp_poly_evaluate( ev, realRat_mpoly_forEvalref(poly), vals[realRat_mpoly_dimref(poly) -1], prec );
    
}

/* result is stored in realRat_mpoly_forEvalref(x) */
void _realRat_mpoly_eval_m1_compApp_poly_ptr( compApp_ptr ev, realRat_mpoly_t poly, slong stop, compApp *const *vals, slong prec ){
    
    if (realRat_mpoly_dimref(poly)==1) {
        compApp_poly_set_realRat_poly(realRat_mpoly_forEvalref(poly), realRat_mpoly_monos_uniref(poly), prec);
    } else {
        /* compute univariate polynomial by recursive evaluation */
        for (slong i = 0; i<realRat_mpoly_lenref(poly); i++) {
            realRat_mmono_ptr m = realRat_mpoly_monosref(poly) + i;
            /* recursive call to evaluate the coeff */
            realRat_mpoly_eval_compApp( realRat_mmono_evalref(m), realRat_mmono_coeffref(m), vals, prec );
            compApp_poly_set_coeff_compApp( realRat_mpoly_forEvalref(poly), realRat_mmono_degref(m), realRat_mmono_evalref(m));
        }
    }
    
    if (realRat_mpoly_dimref(poly) < stop) {
        compApp_poly_evaluate( ev, realRat_mpoly_forEvalref(poly), vals[realRat_mpoly_dimref(poly) -1], prec );
    }
    
}

compApp_poly_ptr realRat_mpoly_eval_m1_compApp_poly_ptr( realRat_mpoly_t poly, compApp *const *vals, slong prec ){
    slong nbvars = realRat_mpoly_dimref(poly);
    _realRat_mpoly_eval_m1_compApp_poly_ptr( NULL, poly, nbvars, vals, prec );
    return realRat_mpoly_forEvalref(poly);
}

void realRat_mpoly_getDegrees( slong * degrees, const realRat_mpoly_t poly ){
    slong nbvars = realRat_mpoly_dimref(poly);
    degrees[nbvars-1] = realRat_mpoly_degref(poly);
    for (slong v = nbvars-2; v>=0; v-- )
        degrees[v] = (realRat_mpoly_degreesref(poly))[nbvars-2-v];
}

void realRat_mpoly_print( realRat_mpoly_t poly ){
    char var[10];
    sprintf(var, "z%ld", realRat_mpoly_dimref(poly) );
    if (realRat_mpoly_dimref(poly)==1) {
        realRat_poly_print_pretty(realRat_mpoly_monos_uniref(poly), var); 
    } else {
        for (slong i=0; i < realRat_mpoly_lenref(poly); i++) {
            
            printf(" ( ");
            
            realRat_mmono_ptr m = realRat_mpoly_monosref(poly) + i;
            realRat_mpoly_print( realRat_mmono_coeffref(m) );
            
            printf(" )*%s^%ld", var, realRat_mmono_degref(m));
            if (i< realRat_mpoly_lenref(poly) -1)
                printf(" + ");
        }
    }
    
}

void realRat_mpoly_clear( realRat_mpoly_t poly ) {
    
//     printf("--- clear from fmpq_mpoly: dim: %ld, len: %ld\n", poly->_dim, poly->_len);
    
    if (realRat_mpoly_dimref(poly)==1){
        
        realRat_poly_clear(realRat_mpoly_monos_uniref(poly));
        
    } else {
        /* clear each monomial */
        for (slong i = 0; i<realRat_mpoly_lenref(poly); i++) {
//             printf("--- clear from fmpq_mpoly: dim: %ld, i: %ld\n", poly->_dim, i);
            realRat_mmono_ptr m = realRat_mpoly_monosref(poly) + i;
            compApp_clear(realRat_mmono_evalref(m));
            /* recursive call to clear the coeff */
            realRat_mpoly_clear( realRat_mmono_coeffref(m) );
//             printf("--- back to clear: dim %ld\n", realRat_mpoly_dimref(poly) );
            /* free the monomial */
            ccluster_free(realRat_mmono_coeffref(m));
//             printf("--- back to clear: dim %ld, after free monomial\n", realRat_mpoly_dimref(poly) );
        }
        /* clear the table of degrees */
        ccluster_free(realRat_mpoly_degreesref(poly));
        /* clear the table of monomials */
        ccluster_free(realRat_mpoly_monosref(poly));
    }
    
    compApp_poly_clear(realRat_mpoly_forEvalref(poly));
    
//     printf("--- end clear from fmpq_mpoly: dim: %ld, len: %ld\n", poly->_dim, poly->_len);
}
