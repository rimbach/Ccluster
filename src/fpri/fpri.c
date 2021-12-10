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
    double t1 = x->low, t2 = x->upp;
    x->low = y->low;
    x->upp = y->upp;
    y->low = t1;
    y->upp = t2;
}

void fpri_set_arb   (fpri_t y, const arb_t  x         ){
    double r, ml, mu;
    r  = mag_get_d( arb_radref(x) );
    ml = arf_get_d( arb_midref(x), ARF_RND_FLOOR);
    mu = arf_get_d( arb_midref(x), ARF_RND_CEIL);
    y->low = -ml + r;
    y->upp =  mu + r;
}

void fpri_get_arb   (arb_t  y, const fpri_t x         ){
    arf_t l, u;
    arf_init(l);
    arf_init(u);
    arf_set_d(l, x->low);
    arf_neg(l, l);
    arf_set_d(u, x->upp);
    arb_set_interval_arf(y, l, u, FPRI_PREC);
    arf_clear(l);
    arf_clear(u);
}

/* supports aliazing */
void fpri_neg     (fpri_t z, const fpri_t x){
    double neg_low = x->upp;
    z->upp = x->low;
    z->low = neg_low;
}

/* supports aliazing */
void fpri_inv(fpri_t z, const fpri_t x){
    if (fpri_contains_zero(x)) {
        z->upp = +INFINITY;
        z->low = +INFINITY;
        return;
    }
    double neg_upp = -x->upp;
    double neg_low = -x->low;
    z->upp = 1./neg_low;
    z->low = 1./neg_upp;
}

/* supports aliazing */
void fpri_sqr     (fpri_t z, const fpri_t x){
    double low, upp;
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
        double tmp = (x->upp)*(x->upp);
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
    double neg_low = x->low + y->upp;
    res->upp = x->upp + y->low;
    res->low = neg_low;
}

/* supports aliazing */
void _fpri_mul   (fpri_t res, const fpri_t x,  const fpri_t y) { 
    
    double upp, low;
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
            double tmp = (x->upp)*(y->low);
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
        double rad = (x->upp + x->low)/2;
        double midupp = -((x->low) - rad);
        double midlow = x->low - rad;
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
        double rad = (x->upp + x->low)/2;
        double midupp = -((x->low) - rad);
        double midlow = x->low - rad;
        rad = rad + (midupp + midlow);
        fprintf(file, "[ %g +/- %g ]", midupp, rad);
    }
}
