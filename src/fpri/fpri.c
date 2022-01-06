/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "fpri.h"

void fpri_swap (fpri_t x, fpri_t y) {
    number t1 = x->low, t2 = x->upp;
    x->low = y->low;
    x->upp = y->upp;
    y->low = t1;
    y->upp = t2;
}

void fpri_set_arb   (fpri_t y, const arb_t  x         ){
    number r, ml, mu;
    r  = mag_get_d( arb_radref(x) );
    ml = arf_get_d( arb_midref(x), ARF_RND_FLOOR);
    mu = arf_get_d( arb_midref(x), ARF_RND_CEIL);
    y->low = -ml + r;
    y->upp =  mu + r;
}

int fpri_get_arb   (arb_t  y, const fpri_t x         ){
    int res = 0;
    arf_t l, u;
    arf_init(l);
    arf_init(u);
    arf_set_d(l, x->low);
    arf_neg(l, l);
    arf_set_d(u, x->upp);
//     if (arf_cmp(l, u)<=0){
        res = 1;
        arb_set_interval_arf(y, l, u, FPRI_PREC);
//     }
    arf_clear(l);
    arf_clear(u);
    return res;
}

/* supports aliazing */
void fpri_neg     (fpri_t z, const fpri_t x){
    number neg_low = x->upp;
    z->upp = x->low;
    z->low = neg_low;
}

void fpri_abs     (fpri_t z, const fpri_t x){
    if (x->upp<=0) /* l<=u<=0 */
        fpri_neg(z, x);
    else { /* 0<u */
        if (x->low > 0) { /* l<0<u */
            z->upp = FPRI_MAX( x->low, x->upp );
            z->low = 0.0;
        }
    }
}

/* supports aliazing */
void fpri_inv(fpri_t z, const fpri_t x){
    if (fpri_contains_zero(x)) {
//         z->upp = +INFINITY;
//         z->low = +INFINITY;
        fpri_set_inf(z);
        return;
    }
    number neg_upp = -x->upp;
    number neg_low = -x->low;
    z->upp = 1./neg_low;
    z->low = 1./neg_upp;
}

/* supports aliazing */
void fpri_sqr     (fpri_t z, const fpri_t x){
    number low, upp;
    if (x->low <= 0) { /* case low of x >= 0 */
        upp = (x->upp)*(x->upp);
        low = (x->low)*(-x->low);
    }
    else if (x->upp <= 0) { /* case upp of x <= 0 */
        upp = (-x->low)*(-x->low);
        low = (-x->upp)*(x->upp);
    }
    else { /* case low of x<0 and upp of x > 0 */
        low=0.0;
        number tmp = (x->upp)*(x->upp);
        upp = (x->low)*(x->low);
        upp = FPRI_MAX( upp, tmp );
    }
    z->upp = upp;
    z->low = low;
}

/* supports aliazing */
void fpri_add     (fpri_t res, const fpri_t x, const fpri_t y){
    res->upp = x->upp + y->upp;
    res->low = x->low + y->low;
}
/* supports aliazing */
void fpri_sub     (fpri_t res, const fpri_t x, const fpri_t y){
    number neg_low = x->low + y->upp;
    res->upp = x->upp + y->low;
    res->low = neg_low;
}

/* supports aliazing */
void _fpri_mul   (fpri_t res, const fpri_t x,  const fpri_t y) { 
    
    number upp, low;
    if (x->low <= 0) { /* case low of x >= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(y->upp);
            low = (x->low)*(-y->low);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(y->upp);
            low = (x->upp)*(y->low);
        }
        else { /* case low of y<0 and upp of y > 0 */
            upp = (x->upp)*(y->upp);
            low = (x->upp)*(y->low);
        }
    }
    else if (x->upp <= 0) { /* case upp of x <= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(-y->low);
            low = (x->low)*(y->upp);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(-y->low);
            low = (-x->upp)*(y->upp);
        }
        else { /* case low of y<0 and upp of y > 0 */
            upp = (-x->low)*(-y->low);
            low = (x->low)*(y->upp);
        }
    }
    else { /* case low of x<0 and upp of x > 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(y->upp);
            low = (x->low)*(y->upp);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(-y->low);
            low = (x->upp)*(y->low);
        }
        else { /* case low of y<0 and upp of y > 0 */
            number tmp = (x->upp)*(y->low);
            low = (x->low)*(y->upp);
            low = FPRI_MAX( low, tmp );
            tmp = (x->upp)*(y->upp);
            upp = (x->low)*(y->low);
            upp = FPRI_MAX( upp, tmp );
        }
    }
    res->upp = upp;
    res->low = low;
}

/* fpri_addmul( z, x, y ): sets z to z + x*y */
void fpri_addmul( fpri_t z, const fpri_t x, const fpri_t y ){
    if (x->low <= 0) { /* case low of x >= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            z->upp = fma( (x->upp), (y->upp), z->upp );
            z->low = fma( (x->low), (-y->low), z->low );
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            z->upp = fma( (-x->low), (y->upp), z->upp );
            z->low = fma( (x->upp), (y->low), z->low );
        }
        else { /* case low of y<0 and upp of y > 0 */
            z->upp = fma( (x->upp), (y->upp), z->upp );
            z->low = fma( (x->upp), (y->low), z->low );
        }
    } else if (x->upp <= 0) { /* case upp of x <= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            z->upp = fma( (x->upp), (-y->low), z->upp );
            z->low = fma( (x->low), (y->upp), z->low );
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            z->upp = fma( (-x->low), (-y->low), z->upp );
            z->low = fma( (-x->upp), (y->upp), z->low );
        }
        else { /* case low of y<0 and upp of y > 0 */
            z->upp = fma( (-x->low), (-y->low), z->upp );
            z->low = fma( (x->low), (y->upp), z->low );
        }
    } else { /* case low of x<0 and upp of x > 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            z->upp = fma( (x->upp), (y->upp), z->upp );
            z->low = fma( (x->low), (y->upp), z->low );
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            z->upp = fma( (-x->low), (-y->low), z->upp );
            z->low = fma( (x->upp), (y->low), z->low );
        }
        else { /* case low of y<0 and upp of y > 0 */
            number tmp = fma( (x->upp), (y->low), z->low );
            number low = fma( (x->low), (y->upp), z->low );
            z->low = FPRI_MAX( low, tmp );
                   tmp = fma( (x->upp), (y->upp), z->upp );
            number upp = fma( (x->low), (y->low), z->upp );
            z->upp = FPRI_MAX( upp, tmp );
        }
    }
}
/* does not support aliazing */
void _fpri_div     (fpri_t res, const fpri_t x, const fpri_t y){
                fpri_inv(res, y);
                fpri_mul(res, res, x);
}
/* support aliazing */
void fpri_div     (fpri_t res, const fpri_t x, const fpri_t y){
    if (res!=x) {
        _fpri_div ( res, x, y );
        return;
    }
    fpri_t temp;
    fpri_init(temp);
    _fpri_div ( temp, x, y );
    fpri_set(res, temp);
    fpri_clear(temp);
}

void fpri_pow_ui           (fpri_t res, const fpri_t x, ulong p){
    fpri_one(res);
    if (p==1)
        fpri_set(res, x);
    else if (p>1) {
        slong pp = p;
        fpri_t temp;
        fpri_init(temp);
        fpri_set(temp, x);
        while (pp>0){
            if (pp%2)
                fpri_mul(res, res, temp);
            pp = pp>>1;
            if (pp>0)
                fpri_sqr(temp, temp);
        }
        fpri_clear(temp);
    }
}

// void fpri_pow_ui           (fpri_t res, const fpri_t x, ulong p){
//     if (p==0)
//         fpri_one(res);
//     else if (p==1)
//         fpri_set(res, x);
//     else if (p==2)
//         fpri_sqr(res, x);
//     else if (p==3){
//         fpri_sqr(res, x);
//         fpri_mul(res, res, x);
//     } else if (p==4) {
//         fpri_sqr(res, x);
//         fpri_sqr(res, res);
//     } else {
//         fmpz_t pfmpz;
//         slong bits, i;
//         fmpz_init(pfmpz);
//         fmpz_set_ui(pfmpz, p);
//         fpri_set(res, x);
//         bits = fmpz_bits(pfmpz);
//         for (i = bits - 2; i >= 0; i--) {
//             fpri_sqr(res, res);
//             if (fmpz_tstbit(pfmpz, i))
//                 fpri_mul(res, res, x);
//         }
//         fmpz_clear(pfmpz);
//     }
// }

void fpri_fprint    (FILE * file, const fpri_t x){
    if (fpri_is_zero(x)) {
        fprintf(file, "[ 0. +/- 0. ]");
    }
    else if (fpri_is_exact(x)){
        fprintf(file, "[ %.5f +/- 0. ]", x->upp);
//         fprintf(file, "[ %.5Le +/- 0. ]", x->upp);
    }
    else if (fpri_is_inf(x)) {
        fprintf(file, "[ %.5f, %.5f]", -x->low, x->upp);
    }
    else {
        number rad = (x->upp + x->low)/2;
        number midupp = -((x->low) - rad);
        number midlow = x->low - rad;
        rad = rad + (midupp + midlow);
        fprintf(file, "[ %g +/- %g ] ([-(%g), %g])", midupp, rad, x->low, x->upp);
    }
}

void fpri_fprint_s    (FILE * file, const fpri_t x){
    if (fpri_is_zero(x)) {
        fprintf(file, "[ 0. +/- 0. ]");
    }
    else if (fpri_is_exact(x)){
        fprintf(file, "[ %.5f +/- 0. ]", x->upp);
//         fprintf(file, "[ %.5Le +/- 0. ]", x->upp);
    }
    else if (fpri_is_inf(x)) {
        fprintf(file, "[ %.5f, %.5f]", -x->low, x->upp);
    }
    else {
        number rad = (x->upp + x->low)/2;
        number midupp = -((x->low) - rad);
        number midlow = x->low - rad;
        rad = rad + (midupp + midlow);
        fprintf(file, "[ %g +/- %g ]", midupp, rad);
    }
}
