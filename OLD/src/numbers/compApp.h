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

#ifndef COMPAPP_H
#define COMPAPP_H

#include "acb.h"
#include "realApp.h"

#include "fmpq.h"

typedef acb_struct compApp;
typedef compApp compApp_t[1];
typedef compApp * compApp_ptr;
typedef const compApp * compApp_srcptr;

#define compApp_realref(x) (&(x)->real)
#define compApp_imagref(x) (&(x)->imag)

/* memory managment */
static __inline__ void compApp_init (compApp_t x) { acb_init (x); }
static __inline__ void compApp_clear(compApp_t x) { acb_clear(x); }

/* members access */
static __inline__ realApp_ptr compApp_real_ptr(compApp_t z                    ) { return acb_realref(z); }
static __inline__ realApp_ptr compApp_imag_ptr(compApp_t z                    ) { return acb_imagref(z); }
static __inline__ void        compApp_get_real(realApp_t re, const compApp_t z) { realApp_set(re, acb_realref(z));}
static __inline__ void        compApp_get_imag(realApp_t im, const compApp_t z) { realApp_set(im, acb_imagref(z));}

/* setting */
static __inline__ void compApp_swap(compApp_t z, compApp_t x      ) { acb_swap(z, x); }
static __inline__ void compApp_zero(compApp_t z                   ) { acb_zero(z); }
static __inline__ void compApp_one (compApp_t z                   ) { acb_one (z); }
static __inline__ void compApp_onei(compApp_t z                   ) { acb_onei(z); }
static __inline__ void compApp_set (compApp_t z, const compApp_t x) { acb_set (z, x); }
static __inline__ void compApp_set_real_realApp(compApp_t x, const realApp_t re) { arb_set(acb_realref(x), re);}
static __inline__ void compApp_set_imag_realApp(compApp_t x, const realApp_t im) { arb_set(acb_imagref(x), im);}

/* arithmetic */
static __inline__ void compApp_abs   ( realApp_t dest, const compApp_t x, slong prec )                   { acb_abs   (dest, x, prec ); }
static __inline__ void compApp_neg   ( compApp_t dest, const compApp_t x )                               { acb_neg   (dest, x ); }
static __inline__ void compApp_sub   ( compApp_t dest, const compApp_t x, const compApp_t y, slong prec) { acb_sub   (dest, x, y, prec); }
static __inline__ void compApp_add   ( compApp_t dest, const compApp_t x, const compApp_t y, slong prec) { acb_add   (dest, x, y, prec); }
static __inline__ void compApp_mul   ( compApp_t dest, const compApp_t x, const compApp_t y, slong prec) { acb_mul   (dest, x, y, prec); }
static __inline__ void compApp_div   ( compApp_t dest, const compApp_t x, const compApp_t y, slong prec) { acb_div   (dest, x, y, prec); }
static __inline__ void compApp_mul_si( compApp_t dest, const compApp_t x, slong y,           slong prec) { acb_mul_si(dest, x, y, prec); }
static __inline__ void compApp_div_si( compApp_t dest, const compApp_t x, slong y,           slong prec) { acb_div_si(dest, x, y, prec); }
static __inline__ void compApp_addmul( compApp_t dest, const compApp_t x, const compApp_t y, slong prec) { acb_addmul(dest, x, y, prec); }
static __inline__ void compApp_exp_pi_i( compApp_t dest, const compApp_t x, slong prec) { acb_exp_pi_i(dest, x, prec); }
static __inline__ void compApp_pow_si( compApp_t dest, const compApp_t x, slong l, slong prec) { acb_pow_si(dest, x, l, prec); }

/* printing */
static __inline__ void compApp_fprint (FILE * file, const compApp_t x)                           { acb_fprint (file, x               ); }
static __inline__ void compApp_fprintd(FILE * file, const compApp_t x, slong digits)             { acb_fprintd(file, x, digits       ); }
static __inline__ void compApp_fprintn(FILE * file, const compApp_t x, slong digits, ulong flags){ acb_fprintn(file, x, digits, flags); }  

static __inline__ void compApp_print (const compApp_t x)                           { acb_print (x               ); }
static __inline__ void compApp_printd(const compApp_t x, slong digits)             { acb_printd(x, digits       ); }
static __inline__ void compApp_printn(const compApp_t x, slong digits, ulong flags){ acb_printn(x, digits, flags); }

/*interval operations*/
static __inline__ int compApp_contains_zero (const compApp_t ball) {
    return acb_contains_zero(ball);
}
static __inline__ int compApp_intersection(compApp_t z, const compApp_t x, const compApp_t y, slong prec) { 
    if (realApp_intersection( compApp_realref(z), compApp_realref(x), compApp_realref(y), prec) !=0)
        return realApp_intersection( compApp_imagref(z), compApp_imagref(x), compApp_imagref(y), prec);
    else return 0; 
}
#endif