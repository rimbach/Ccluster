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

#ifdef NUMBERS_INLINE_C
#define NUMBERS_INLINE
#else
#define NUMBERS_INLINE static __inline__
#endif

#include <math.h>
#include "arf.h"
#include "arb.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef arb_struct realApp;
typedef realApp realApp_t[1];
typedef realApp * realApp_ptr;
typedef const realApp * realApp_srcptr;

/* memory managment */
NUMBERS_INLINE void realApp_init (realApp_t x) { arb_init (x); }
NUMBERS_INLINE void realApp_clear(realApp_t x) { arb_clear(x); }

/* setting */
NUMBERS_INLINE void realApp_zero      (realApp_t x                                ) { arb_zero     (x); }
NUMBERS_INLINE void realApp_one       (realApp_t x                                ) { arb_one      (x); }
NUMBERS_INLINE void realApp_pi        (realApp_t x, slong prec                    ) { arb_const_pi(x, prec); }
NUMBERS_INLINE void realApp_set       (realApp_t y, const realApp_t x             ) { arb_set      (y, x); }
NUMBERS_INLINE void realApp_set_fmpq  (realApp_t y, const fmpq_t    x, slong prec ) { arb_set_fmpq (y, x, prec); }
NUMBERS_INLINE void realApp_set_fmpz  (realApp_t y, const fmpz_t    x, slong prec ) { arb_set_fmpz (y, x); }
NUMBERS_INLINE void realApp_set_d     (realApp_t y, const double    x )             { arb_set_d (y, x); }

NUMBERS_INLINE void realApp_set_si     (realApp_t y, const slong    x )             { arb_set_si (y, x); }

NUMBERS_INLINE void realApp_set_rad_zero (realApp_t x) { mag_zero(arb_radref(x)); }

/* special cases */
NUMBERS_INLINE void realApp_pos_inf    (realApp_t x                                ) { arb_pos_inf     (x); }
NUMBERS_INLINE void realApp_neg_inf    (realApp_t x                                ) { arb_neg_inf     (x); }

/* comparisons */
NUMBERS_INLINE int realApp_eq(const realApp_t x, const realApp_t y) { return arb_eq(x,y); }
NUMBERS_INLINE int realApp_ne(const realApp_t x, const realApp_t y) { return arb_ne(x,y); }
NUMBERS_INLINE int realApp_lt(const realApp_t x, const realApp_t y) { return arb_lt(x,y); }
NUMBERS_INLINE int realApp_le(const realApp_t x, const realApp_t y) { return arb_le(x,y); }
NUMBERS_INLINE int realApp_gt(const realApp_t x, const realApp_t y) { return arb_gt(x,y); }
NUMBERS_INLINE int realApp_ge(const realApp_t x, const realApp_t y) { return arb_ge(x,y); }
NUMBERS_INLINE int realApp_is_negative(const realApp_t x) { return arb_is_negative(x); }
NUMBERS_INLINE int realApp_is_positive(const realApp_t x) { return arb_is_positive(x); }
NUMBERS_INLINE int realApp_is_zero(const realApp_t x) { return arb_is_zero(x); }

/* ball operations */
NUMBERS_INLINE int  realApp_is_finite(const realApp_t x) { return arb_is_finite( x ); }
NUMBERS_INLINE void realApp_get_rad_realApp(realApp_t z, const realApp_t x) { arb_get_rad_arb( z, x ); }
NUMBERS_INLINE void realApp_get_mid_realApp(realApp_t z, const realApp_t x) { arb_get_mid_arb( z, x ); }

NUMBERS_INLINE int  realApp_get_unique_si(slong * z, const realApp_t x) {
    fmpz_t zf;
    fmpz_init(zf);
    int res = arb_get_unique_fmpz( zf, x);
    res = res && fmpz_fits_si(zf);
    if (res)
        *z = fmpz_get_si(zf);
    fmpz_clear(zf);
    return res;
}

/* performed inplace */
NUMBERS_INLINE void realApp_add_error(realApp_t z, const realApp_t x) { arb_add_error( z, x ); }
/* sets x to a ball centered in zero containing x */
NUMBERS_INLINE void realApp_center_in_zero(realApp_t x) { 
    arb_add_error_arf( x, arb_midref(x) ); 
    arf_zero(arb_midref(x));
}
/* interval operations */
NUMBERS_INLINE int  realApp_intersection(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { 
    return arb_intersection( z, x, y, prec ); 
}
NUMBERS_INLINE void  realApp_union(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { 
    arb_union( z, x, y, prec ); 
}
NUMBERS_INLINE int realApp_contains (const realApp_t x, const realApp_t y) {
    return arb_contains(x,y);
}
NUMBERS_INLINE int realApp_contains_zero (const realApp_t x) {
    return arb_contains_zero(x);
}

/* arithmetic operations */
NUMBERS_INLINE void realApp_abs   ( realApp_t dest, const realApp_t x )                   { arb_abs   (dest, x); }
NUMBERS_INLINE void realApp_inv   ( realApp_t dest, const realApp_t x, slong prec )       { arb_inv   (dest, x, prec); }
NUMBERS_INLINE void realApp_add(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_add(z, x, y, prec); }
NUMBERS_INLINE void realApp_add_si(realApp_t z, const realApp_t x, slong y, slong prec) { arb_add_si(z, x, y, prec); }
NUMBERS_INLINE void realApp_sub(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_sub(z, x, y, prec); }
NUMBERS_INLINE void realApp_sub_si(realApp_t z, const realApp_t x, slong y, slong prec) { arb_sub_si(z, x, y, prec); }
NUMBERS_INLINE void realApp_mul(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_mul(z, x, y, prec); }
NUMBERS_INLINE void realApp_neg(realApp_t z, const realApp_t x) { arb_neg(z, x); }
NUMBERS_INLINE void realApp_div(realApp_t z, const realApp_t x, const realApp_t y, slong prec) { arb_div(z, x, y, prec); }
NUMBERS_INLINE void realApp_div_si(realApp_t z, const realApp_t x, slong y, slong prec) { arb_div_si(z, x, y, prec); }
NUMBERS_INLINE void realApp_div_ui(realApp_t z, const realApp_t x, ulong y, slong prec) { arb_div_ui(z, x, y, prec); }
NUMBERS_INLINE void realApp_mul_si( realApp_t dest, const realApp_t x, slong y,           slong prec) { arb_mul_si(dest, x, y, prec); }
NUMBERS_INLINE void realApp_pow_ui(realApp_t y, const realApp_t x, ulong e, slong prec) { arb_pow_ui(y, x, e, prec); }
NUMBERS_INLINE void realApp_root_ui(realApp_t y, const realApp_t x, ulong e, slong prec) { arb_root_ui(y, x, e, prec); }
NUMBERS_INLINE void realApp_mul_2exp_si(realApp_t y, const realApp_t x, slong e) { arb_mul_2exp_si(y, x, e); }

NUMBERS_INLINE void realApp_sqr (realApp_t z, const realApp_t x, slong prec) { arb_sqr (z, x, prec); }
NUMBERS_INLINE void realApp_sqrt(realApp_t z, const realApp_t x, slong prec) { arb_sqrt(z, x, prec); }

/* logarithm */
NUMBERS_INLINE void realApp_log(realApp_t z, const realApp_t x, slong prec) { arb_log(z, x, prec); }
NUMBERS_INLINE void realApp_log_base_ui(realApp_t z, const realApp_t x, ulong base, slong prec) { arb_log_base_ui(z, x, base, prec); }

/* other */
NUMBERS_INLINE slong realApp_ceil_si(const realApp_t x, slong prec){
    slong res;
    arf_t ubound;
    arf_init(ubound);
    arb_get_ubound_arf(ubound, x, prec);
    res = (slong) ceil(arf_get_d(ubound,  ARF_RND_CEIL));
    arf_clear(ubound);
    return res;
}


NUMBERS_INLINE void realApp_ceil_fmpz(fmpz_t res, const realApp_t x, slong prec){
    arb_t f;
    arb_init(f);
    arb_ceil(f, x, prec);
    arb_get_unique_fmpz(res, f);
    arb_clear(f);
}

NUMBERS_INLINE void realApp_floor_fmpz(fmpz_t res, const realApp_t x, slong prec){
    arb_t f;
    arb_init(f);
    arb_floor(f, x, prec);
    arb_get_unique_fmpz(res, f);
    arb_clear(f);
}

/* Returns 1 if x is strictly positive, -1 if x is strictly negative, and 0 if x is zero or a ball containing zero so that its sign is not determined. */
NUMBERS_INLINE int realApp_sgn_nonzero(const realApp_t z) { return arb_sgn_nonzero(z); }

/* printing */
NUMBERS_INLINE void realApp_fprint (FILE * file, const realApp_t x)                           { arb_fprint (file, x               ); }
NUMBERS_INLINE void realApp_fprintd(FILE * file, const realApp_t x, slong digits)             { arb_fprintd(file, x, digits       ); }
NUMBERS_INLINE void realApp_fprintn(FILE * file, const realApp_t x, slong digits, ulong flags){ arb_fprintn(file, x, digits, flags); }  

NUMBERS_INLINE void realApp_print (const realApp_t x)                           { arb_print (x               ); }
NUMBERS_INLINE void realApp_printd(const realApp_t x, slong digits)             { arb_printd(x, digits       ); }
NUMBERS_INLINE void realApp_printn(const realApp_t x, slong digits, ulong flags){ arb_printn(x, digits, flags); }

NUMBERS_INLINE int realApp_is_exact(const realApp_t x) { return arb_is_exact(x); }

/* for a ball m +/- r, relative accuracy of [max(1,|m|) +/- r]  */
NUMBERS_INLINE int realApp_check_relOne_accuracy( const realApp_t z, slong prec) {
    return ( arb_rel_one_accuracy_bits(z) >= prec );
}

NUMBERS_INLINE slong realApp_get_relOne_accuracy( const realApp_t z) {
    return arb_rel_one_accuracy_bits(z);
}

NUMBERS_INLINE int realApp_check_absolute_accuracy( const realApp_t z, slong prec) {
    return mag_cmp_2exp_si(arb_radref((z)), prec)<=0;
}

NUMBERS_INLINE slong realApp_get_absolute_accuracy( const realApp_t z) {
    double logerror = mag_get_d_log2_approx(arb_radref((z)));
    return - ( (slong) logerror ) -1;
}

#ifdef __cplusplus
}
#endif

#endif
