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

#ifndef REALAPP_POLY_H
#define REALAPP_POLY_H

#ifdef POLYNOMIALS_INLINE_C
#define POLYNOMIALS_INLINE
#else
#define POLYNOMIALS_INLINE static __inline__
#endif

#include "arb_poly.h"
#include "base/base.h"
#include "numbers/realApp.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef arb_poly_struct realApp_poly;
typedef arb_poly_struct realApp_poly_t[1];
typedef arb_poly_struct * realApp_poly_ptr;

/* memory managment */
POLYNOMIALS_INLINE void realApp_poly_init      (realApp_poly_t poly)            { arb_poly_init      (poly);}
POLYNOMIALS_INLINE void realApp_poly_init2     (realApp_poly_t poly, slong len) { arb_poly_init2     (poly, len);}
POLYNOMIALS_INLINE void realApp_poly_clear     (realApp_poly_t poly)            { arb_poly_clear     (poly);}
POLYNOMIALS_INLINE void realApp_poly_fit_length(realApp_poly_t poly, slong len) { arb_poly_fit_length(poly, len); }
POLYNOMIALS_INLINE void realApp_poly_set_length(realApp_poly_t poly, slong len) { _arb_poly_set_length(poly, len); }

POLYNOMIALS_INLINE slong realApp_poly_degree(const realApp_poly_t poly) { return arb_poly_degree(poly); }
/* bounds are not checked */
POLYNOMIALS_INLINE realApp_srcptr realApp_poly_getCoeff(const realApp_poly_t poly, slong degree) { return poly->coeffs + degree; }
POLYNOMIALS_INLINE void  realApp_poly_swap( realApp_poly_t poly1, realApp_poly_t poly2 ) { arb_poly_swap(poly1, poly2); }

/* printing */
POLYNOMIALS_INLINE void realApp_poly_fprintd(FILE * file, const realApp_poly_t poly, slong digits) {
    arb_poly_fprintd(file, poly, digits);
}
POLYNOMIALS_INLINE void realApp_poly_printd(const realApp_poly_t poly, slong digits) {
    arb_poly_printd(poly, digits);
}

/* setting */
POLYNOMIALS_INLINE void realApp_poly_zero(realApp_poly_t poly) { arb_poly_zero(poly); }
POLYNOMIALS_INLINE void realApp_poly_one (realApp_poly_t poly) { arb_poly_one (poly); }
POLYNOMIALS_INLINE void realApp_poly_set (realApp_poly_t dest, const realApp_poly_t src) { arb_poly_set (dest, src); }
POLYNOMIALS_INLINE void realApp_poly_set_coeff_si (realApp_poly_t dest, slong n, slong x)           { arb_poly_set_coeff_si (dest, n, x); }
POLYNOMIALS_INLINE void realApp_poly_set_coeff_realApp(realApp_poly_t dest, slong n, const realApp_t x) { arb_poly_set_coeff_arb(dest, n, x); }

/* Arithmetic */
POLYNOMIALS_INLINE void realApp_poly_mul( realApp_poly_t res, const realApp_poly_t poly1, const realApp_poly_t poly2, slong prec) {
    arb_poly_mul (res, poly1, poly2, prec);
}
POLYNOMIALS_INLINE void realApp_poly_mullow( realApp_poly_t res, const realApp_poly_t poly1, const realApp_poly_t poly2, slong n, slong prec) {
    arb_poly_mullow (res, poly1, poly2, n, prec);
}
POLYNOMIALS_INLINE void realApp_poly_sub( realApp_poly_t res, const realApp_poly_t poly1, const realApp_poly_t poly2, slong prec) {
    arb_poly_sub (res, poly1, poly2, prec);
}
POLYNOMIALS_INLINE void realApp_poly_add( realApp_poly_t res, const realApp_poly_t poly1, const realApp_poly_t poly2, slong prec) {
    arb_poly_add (res, poly1, poly2, prec);
}
POLYNOMIALS_INLINE void realApp_poly_shift_left( realApp_poly_t res, const realApp_poly_t poly, slong n) {
    arb_poly_shift_left (res, poly, n);
}
POLYNOMIALS_INLINE void realApp_poly_neg( realApp_poly_t res, const realApp_poly_t poly) {
    arb_poly_neg (res, poly);
}

/* evaluation */
POLYNOMIALS_INLINE  void realApp_poly_evaluate(realApp_t y, const realApp_poly_t f, const realApp_t x, slong prec){
    arb_poly_evaluate_rectangular(y, f, x, prec);
}

// /*order one center evaluation */
// void realApp_poly_evaluate_order_one( realApp_t y, const realApp_poly_t f, const realApp_poly_t fder, const realApp_t x, slong prec);

/* sum absolute values of coeffs */
void realApp_poly_sum_abs_coeffs( realApp_t res, const realApp_poly_t f, slong prec );

/* accuracy */
int realApp_poly_checkAccuracy( const realApp_poly_t poly, slong prec);

slong realApp_poly_getAccuracy_min( const realApp_poly_t poly);
slong realApp_poly_getAccuracy_max( const realApp_poly_t poly);

/* derivation */
POLYNOMIALS_INLINE  void realApp_poly_derivative(realApp_poly_t fp, const realApp_poly_t f, slong prec){
    arb_poly_derivative(fp, f, prec);
}

/* Graeffe iterations */
void realApp_poly_oneGraeffeIteration_coeff( realApp_ptr coeff, realApp_srcptr f, slong index, slong len, slong prec);
/* requires: f is a partially computed polynomial of degree len-1 */
/*           its coeffs from 0 to n + delta exist, where delta = min(n, len-1-n) */
void realApp_poly_oneGraeffeIteration_first_n_coeff( realApp_poly_t res, const realApp_poly_t f, slong n, slong len, slong prec);
/* requires: f is initialized */
void realApp_poly_oneGraeffeIteration_in_place( realApp_poly_t f, slong prec );

#ifdef __cplusplus
}
#endif

#endif
