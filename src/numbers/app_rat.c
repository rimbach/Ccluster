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

#include "app_rat.h"

/*converts a disk to a compApp*/
void compApp_set_compDsk( compApp_t res, const compRat_t center, const realRat_t radius, slong prec){
    
    mag_t rad;
    mag_init(rad);
    
    compApp_set_compRat(res, center, prec);
    
    mag_set_fmpz(rad, realRat_numref(radius));
    mag_div_fmpz(rad, rad, realRat_denref(radius));
    mag_add( arb_radref( compApp_realref(res) ), arb_radref( compApp_realref(res) ), rad);
    mag_add( arb_radref( compApp_imagref(res) ), arb_radref( compApp_imagref(res) ), rad);
    
    mag_clear(rad);
}

/* arithmetic */
void realApp_mul_realRat( realApp_t x, const realApp_t y, const realRat_t z, slong prec ){
    arb_mul_fmpz(x, y, fmpq_numref(z), prec);
    arb_div_fmpz(x, x, fmpq_denref(z), prec);
}

void realApp_mul_realRat_in_place( realApp_t x, const realRat_t y, slong prec ){
    arb_mul_fmpz(x, x, fmpq_numref(y), prec);
    arb_div_fmpz(x, x, fmpq_denref(y), prec);
}

void compApp_mul_realRat( compApp_t x, const compApp_t y, const realRat_t z, slong prec ){
    acb_mul_fmpz(x, y, fmpq_numref(z), prec);
    acb_div_fmpz(x, x, fmpq_denref(z), prec);
}

void compApp_mul_realRat_in_place( compApp_t x, const realRat_t y, slong prec ) {
    acb_mul_fmpz(x, x, fmpq_numref(y), prec);
    acb_div_fmpz(x, x, fmpq_denref(y), prec);
}

void _compApp_mul_compRat( compApp_t z, const compApp_t x, const compRat_t y, slong prec ){
#define a compApp_realref(x)
#define b compApp_imagref(x)
#define c compRat_realref(y)
#define d compRat_imagref(y)
#define e compApp_realref(z)
#define f compApp_imagref(z)  
    /* TODO 
    if ( realApp_is_zero(b) ) {
    }
    else if ( realRat_is_zero(d) ) {
    }
    else if ( realApp_is_zero(a) ) {
    }
    else if ( realRat_is_zero(c) ) {
    }
    else {
        if (z == x) {
        }
        else {
    */
            /* Gauss multiplication: e = ac-bd
                                     f = (a+b)(c+d) - ac - bd
            */
            realApp_t t;
            realApp_init(t);
            realRat_t u;
            realRat_init(u);
            
            realApp_mul_realRat(e, a, c, prec);
            realApp_mul_realRat(t, b, d, prec);
            
            realRat_add(u, c, d);
            realApp_add(f, a, b, prec);
            realApp_mul_realRat(f, f, u, prec);
            realApp_sub(f, f, e, prec);
            realApp_sub(f, f, t, prec);
            
            realApp_sub(e, e, t, prec);
            
            realApp_clear(t);
            realRat_clear(u);
    /*
        }
    }
    */
#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
}

void compApp_mul_compRat( compApp_t z, const compApp_t x, const compRat_t y, slong prec ){
    
    if (z==x) {
        compApp_t temp;
        compApp_init(temp);
        _compApp_mul_compRat( temp, x, y, prec );
        compApp_set(z, temp);
        compApp_clear(temp);
    }
    else
        _compApp_mul_compRat( z, x, y, prec );
}

/*getting a realRat lying in the ball defined by a realApp */
void realApp_get_realRat( realRat_t res, realApp_t x){
    arf_srcptr midpoint;
    fmpz_t mantissa, exponent;
    fmpz_init(mantissa);
    fmpz_init(exponent);
    midpoint = arb_mid_ptr(x);
    arf_get_fmpz_2exp( mantissa, exponent, midpoint);
    int sign = fmpz_sgn(exponent);
    ulong exponent_ulong;
    fmpz_abs(exponent, exponent);
    
    realRat_set_fmpz(res, mantissa);
    
    if (fmpz_abs_fits_ui(exponent)==0) { /*TODO: not safe at all...*/
        realRat_clear(res);
    }
    else {
        exponent_ulong = fmpz_get_ui(exponent);
        if (sign==1)
            fmpq_mul_2exp(res, res, exponent_ulong);
        else
            fmpq_div_2exp(res, res, exponent_ulong);
    }
    
    fmpz_clear(mantissa);
    fmpz_clear(exponent);
}
