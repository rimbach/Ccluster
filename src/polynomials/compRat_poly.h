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

#ifndef COMPRAT_POLY_H
#define COMPRAT_POLY_H

#ifdef POLYNOMIALS_INLINE_C
#define POLYNOMIALS_INLINE
#else
#define POLYNOMIALS_INLINE static __inline__
#endif

#include "base/base.h"
#include "polynomials/realRat_poly.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    realRat_poly real;
    realRat_poly imag;
} compRat_poly;

typedef compRat_poly compRat_poly_t[1];
typedef compRat_poly * compRat_poly_ptr;

#define compRat_poly_realref(x) (&(x)->real)
#define compRat_poly_imagref(x) (&(x)->imag)

/* memory managment */
POLYNOMIALS_INLINE void compRat_poly_init( compRat_poly_t x ){
    realRat_poly_init( compRat_poly_realref(x) );
    realRat_poly_init( compRat_poly_imagref(x) );
}

POLYNOMIALS_INLINE void compRat_poly_init2( compRat_poly_t x, slong size ){
    realRat_poly_init2( compRat_poly_realref(x), size );
    realRat_poly_init2( compRat_poly_imagref(x), size );
}

POLYNOMIALS_INLINE void compRat_poly_fit_length( compRat_poly_t x, slong size ){
    realRat_poly_fit_length( compRat_poly_realref(x), size );
    realRat_poly_fit_length( compRat_poly_imagref(x), size );
}

POLYNOMIALS_INLINE void compRat_poly_clear( compRat_poly_t x ){
    realRat_poly_clear( compRat_poly_realref(x) );
    realRat_poly_clear( compRat_poly_imagref(x) );
}

/* members access */
POLYNOMIALS_INLINE realRat_poly_ptr compRat_poly_real_ptr(compRat_poly_t z) { return compRat_poly_realref(z); }
POLYNOMIALS_INLINE realRat_poly_ptr compRat_poly_imag_ptr(compRat_poly_t z) { return compRat_poly_imagref(z); }

POLYNOMIALS_INLINE void compRat_poly_get_real(realRat_poly_t re, const compRat_poly_t z) {
    realRat_poly_set(re, compRat_poly_realref(z));
}

POLYNOMIALS_INLINE void compRat_poly_get_imag(realRat_poly_t im, const compRat_poly_t z) {
    realRat_poly_set(im, compRat_poly_imagref(z));
}

/* canonicalise */
POLYNOMIALS_INLINE void compRat_poly_canonicalise(compRat_poly_t poly) { 
    realRat_poly_canonicalise(compRat_poly_realref(poly)); 
    realRat_poly_canonicalise(compRat_poly_imagref(poly));
}
/* multiple with same roots and integer coeffs */
POLYNOMIALS_INLINE void compRat_poly_integers(compRat_poly_t poly) { 
    realRat_poly_integers(compRat_poly_realref(poly)); 
    realRat_poly_integers(compRat_poly_imagref(poly));
}

/* testing */
POLYNOMIALS_INLINE int  compRat_poly_is_zero(const compRat_poly_t poly) { 
    return  realRat_poly_is_zero(compRat_poly_realref(poly)) 
           && realRat_poly_is_zero(compRat_poly_imagref(poly)); 
}

/* assume real and imag parts of poly are canonical */
POLYNOMIALS_INLINE slong compRat_poly_degree(const compRat_poly_t poly) {
    return CCLUSTER_MAX(realRat_poly_degree(compRat_poly_realref(poly)),realRat_poly_degree(compRat_poly_imagref(poly)));
}

/* setting */
POLYNOMIALS_INLINE void compRat_poly_set(compRat_poly_t z, const compRat_poly_t x) {
    realRat_poly_set(compRat_poly_realref(z), compRat_poly_realref(x));
    realRat_poly_set(compRat_poly_imagref(z), compRat_poly_imagref(x));
}

POLYNOMIALS_INLINE void compRat_poly_zero(compRat_poly_t z) {
    realRat_poly_zero(compRat_poly_realref(z));
    realRat_poly_zero(compRat_poly_imagref(z));
}

POLYNOMIALS_INLINE void compRat_poly_set_realRat_poly(compRat_poly_t z, const realRat_poly_t x) {
    realRat_poly_set(compRat_poly_realref(z), x);
    realRat_poly_zero(compRat_poly_imagref(z));
}

POLYNOMIALS_INLINE void compRat_poly_set2_realRat_poly
(compRat_poly_t z, const realRat_poly_t x, const realRat_poly_t y) {
    realRat_poly_set(compRat_poly_realref(z), x);
    realRat_poly_set(compRat_poly_imagref(z), y);
}

/* Finds the quotient Q and remainder R of the Euclidean division of poly1 by poly2.*/
POLYNOMIALS_INLINE void compRat_poly_divrem( compRat_poly_t Q, compRat_poly_t R, const compRat_poly_t poly1, const realRat_poly_t poly2) {
    realRat_poly_divrem( compRat_poly_realref(Q), compRat_poly_realref(R), compRat_poly_realref(poly1), poly2);
    realRat_poly_divrem( compRat_poly_imagref(Q), compRat_poly_imagref(R), compRat_poly_imagref(poly1), poly2);
}
/* Finds the remainder R of the Euclidean division of poly1 by poly2.*/
POLYNOMIALS_INLINE void compRat_poly_rem( compRat_poly_t R, const compRat_poly_t poly1, const realRat_poly_t poly2) {
    realRat_poly_rem( compRat_poly_realref(R), compRat_poly_realref(poly1), poly2);
    realRat_poly_rem( compRat_poly_imagref(R), compRat_poly_imagref(poly1), poly2);
}
POLYNOMIALS_INLINE void compRat_poly_derivative( compRat_poly_t dest, const compRat_poly_t src) {
    realRat_poly_derivative( compRat_poly_realref(dest), compRat_poly_realref(src));
    realRat_poly_derivative( compRat_poly_imagref(dest), compRat_poly_imagref(src));
}

/* printing */
int compRat_poly_fprint(FILE * file, const compRat_poly_t z); 

POLYNOMIALS_INLINE int compRat_poly_print(const compRat_poly_t x) {
    return compRat_poly_fprint(stdout, x);
}

int compRat_poly_fprint_pretty(FILE * file, const compRat_poly_t poly, const char * var);

POLYNOMIALS_INLINE int compRat_poly_print_pretty(const compRat_poly_t poly, const char * var) {
    return compRat_poly_fprint_pretty(stdout, poly, var);
}

#ifdef __cplusplus
}
#endif

#endif
