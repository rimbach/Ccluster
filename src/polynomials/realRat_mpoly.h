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

#ifndef REALRAT_MPOLY_H
#define REALRAT_MPOLY_H

#ifdef POLYNOMIALS_INLINE_C
#define POLYNOMIALS_INLINE
#else
#define POLYNOMIALS_INLINE static __inline__
#endif

#include "flint/flint.h"
#include "flint/fmpq_mpoly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"

#ifdef __cplusplus
extern "C" {
#endif
    
struct realRat_mpoly;

/* a monomial of a multivariate polynomial        */
/* assume _dim > 1, represents:                   */
/* _coeff*z_dim^(_deg)                            */
/* where _coeff is in Q[z_1, z_2, ..., z_{_dim-1} */
/* _eval is to store evaluation of _coeff at a    */ 
/* {_dim-1}-dimensional vector of compApp         */
struct realRat_mmono {
    slong                 _deg;    /* degree of variable           */
    slong                 _dim;    /* number of variables in coeff */
    struct realRat_mpoly *_coeff;  /* coefficient: a mpoly */
    compApp               _eval;   /* compApp for result of evaluation */
};

typedef struct realRat_mmono * realRat_mmono_ptr;

#define realRat_mmono_degref(x)       (x->_deg)
#define realRat_mmono_dimref(x)       (x->_dim)
#define realRat_mmono_coeffref(x)     (x->_coeff)
#define realRat_mmono_evalref(x)      (&(x)->_eval)

/* a multivariate polynomial in Q[z_1,..., z_dim]               */
/* assume _dim>=1 (number of variables)                         */
/* _deg is the degree in z_dim                                  */
/* _len is the number of monomials with non-zero coefficients   */
/* if _dim = 1 :                                                */
/*     _monos is an uninitialized table of realRat_mmono        */
/*     _monos_uni is an initialized realRat_poly of degree _deg */
/*                and _length _deg+1                            */
/* else (_dim>1):                                               */
/*     _monos is a table of _len realRat_mmono                  */
/*            ordered by decreasing degree                      */
/*     _monos_uni is an uninitialized pointer on realRat_poly   */
/* _forEval is to store the evaluation of the polynomial at a   */
/*  {_dim-1}-dimensional vector of compApp                      */
struct realRat_mpoly {
    slong                  _deg;       /* degree of polynomial              */
    slong                  _dim;       /* number of variables in polynomial */
    slong                  _len;       /* number of terms                   */    
    slong                  *_degrees;  /* degrees in z_{_dim-1},...,z_1     */
    struct realRat_mmono   *_monos;     /* a table of monomials              */
    realRat_poly           _monos_uni; /* or a realRat_poly when _dim = 1   */
    compApp_poly           _forEval;   /* compApp_poly for evaluation       */
};

typedef struct realRat_mpoly realRat_mpoly_t[1];
typedef struct realRat_mpoly * realRat_mpoly_ptr;

#define realRat_mpoly_degref(x)       (x->_deg)
#define realRat_mpoly_dimref(x)       (x->_dim)
#define realRat_mpoly_lenref(x)       (x->_len)
#define realRat_mpoly_degreesref(x)   (x->_degrees)
#define realRat_mpoly_monosref(x)     (x->_monos)
#define realRat_mpoly_monos_uniref(x) (&(x)->_monos_uni)
#define realRat_mpoly_forEvalref(x)   (&(x)->_forEval)

/* init to zero polynomial */
void realRat_mpoly_init( realRat_mpoly_t poly );

/* init from fmpq_mpoly and context */
void _realRat_mpoly_init_fmpq_mpoly( realRat_mpoly_t poly, slong nbvars, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx );

void realRat_mpoly_init_fmpq_mpoly( realRat_mpoly_t poly, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx );

/* evaluate at a _dim-dimensional vector of compApp          */
/* result is stored in ev which is assumed to be initialized */
void realRat_mpoly_eval_compApp( compApp_t ev, realRat_mpoly_t poly , compApp *const *vals, slong prec );

/* evaluate at a {_dim-1}-dimensional vector of compApp */
/* result is stored in _forEval                         */
/* a pointer on result is returned                      */
void _realRat_mpoly_eval_m1_compApp_poly_ptr( compApp_ptr ev, realRat_mpoly_t poly, slong stop, compApp *const *vals, slong prec );
compApp_poly_ptr realRat_mpoly_eval_m1_compApp_poly_ptr( realRat_mpoly_t poly , compApp *const *vals, slong prec );

/* give degrees in z_1,...,z_dim in the table degrees, assumed to be initialized
 * at appropriated size */
void realRat_mpoly_getDegrees( slong * degrees, const realRat_mpoly_t poly );

void realRat_mpoly_print( realRat_mpoly_t poly );

void realRat_mpoly_clear( realRat_mpoly_t poly );
    

#ifdef __cplusplus
}
#endif

#endif
