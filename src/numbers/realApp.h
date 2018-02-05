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

#ifndef REALAPP_H
#define REALAPP_H

#include "arb.h"
#include "fmpq.h"

typedef arb_struct realApp;
typedef realApp realApp_t[1];
typedef realApp * realApp_ptr;

/* memory managment */
static __inline__ void realApp_init (realApp_t x) { arb_init (x); }
static __inline__ void realApp_clear(realApp_t x) { arb_clear(x); }

/* setting */
static __inline__ void realApp_zero      (realApp_t x                                ) { arb_zero     (x); }
static __inline__ void realApp_one       (realApp_t x                                ) { arb_one      (x); }
static __inline__ void realApp_set       (realApp_t y, const realApp_t x             ) { arb_set      (y, x); }
static __inline__ void realApp_set_fmpq  (realApp_t y, const fmpq_t    x, slong prec ) { arb_set_fmpq (y, x, prec); }

/* comparisons */
static __inline__ int realApp_eq(const realApp_t x, const realApp_t y) { return arb_eq(x,y); }
static __inline__ int realApp_ne(const realApp_t x, const realApp_t y) { return arb_ne(x,y); }
static __inline__ int realApp_lt(const realApp_t x, const realApp_t y) { return arb_lt(x,y); }
static __inline__ int realApp_le(const realApp_t x, const realApp_t y) { return arb_le(x,y); }
static __inline__ int realApp_gt(const realApp_t x, const realApp_t y) { return arb_gt(x,y); }
static __inline__ int realApp_ge(const realApp_t x, const realApp_t y) { return arb_ge(x,y); }

/* ball operations */
static __inline__ int  realApp_is_finite(const realApp_t x) { return arb_is_finite( x ); }
static __inline__ void realApp_get_rad_realApp(realApp_t z, const realApp_t x) { arb_get_rad_arb( z, x ); }
static __inline__ void realApp_get_mid_realApp(realApp_t z, const realApp_t x) { arb_get_mid_arb( z, x ); }
/* interval operations */
static __inline__ int  realApp_intersection(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { 
    return arb_intersection( z, x, y, prec ); 
}

/* arithmetic operations */
static __inline__ void realApp_add(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_add(z, x, y, prec); }
static __inline__ void realApp_sub(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_sub(z, x, y, prec); }
static __inline__ void realApp_mul(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_mul(z, x, y, prec); }
static __inline__ void realApp_div(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_div(z, x, y, prec); }
static __inline__ void realApp_pow_ui(realApp_t y, const realApp_t x, ulong e, slong prec) { arb_pow_ui(y, x, e, prec); }
static __inline__ void realApp_root_ui(realApp_t y, const realApp_t x, ulong e, slong prec) { arb_root_ui(y, x, e, prec); }

/* printing */
static __inline__ void realApp_fprint (FILE * file, const realApp_t x)                           { arb_fprint (file, x               ); }
static __inline__ void realApp_fprintd(FILE * file, const realApp_t x, slong digits)             { arb_fprintd(file, x, digits       ); }
static __inline__ void realApp_fprintn(FILE * file, const realApp_t x, slong digits, ulong flags){ arb_fprintn(file, x, digits, flags); }  

static __inline__ void realApp_print (const realApp_t x)                           { arb_print (x               ); }
static __inline__ void realApp_printd(const realApp_t x, slong digits)             { arb_printd(x, digits       ); }
static __inline__ void realApp_printn(const realApp_t x, slong digits, ulong flags){ arb_printn(x, digits, flags); }

#endif