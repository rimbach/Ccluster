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

#include "fpci.h"

/* does not support aliazing */
void _fpci_sqr   ( fpci_t dest, const fpci_t x ){
#define a fpci_realref(x)
#define b fpci_imagref(x)
#define e fpci_realref(dest)
#define f fpci_imagref(dest)
    fpri_sqr(e, a);
    fpri_sqr(f, b);
    fpri_sub(e, e, f);
    _fpri_mul(f, a, b);
    fpri_mul_si(f, f, 2);
#undef a
#undef b
#undef e
#undef f 
}

/* supports aliazing */
void fpci_sqr   ( fpci_t dest, const fpci_t x ){
    if (dest!=x){
        _fpci_sqr   ( dest, x );
        return;
    }
    fpci_t temp;
    fpci_init(temp);
    _fpci_sqr   ( temp, x );
    fpci_set(dest, temp);
    fpci_clear(temp);
}

/* does not support aliazing */
void _fpci_inv( fpci_t dest, const fpci_t x ) {
#define a fpci_realref(x)
#define b fpci_imagref(x)
#define e fpci_realref(dest)
#define f fpci_imagref(dest)
    if ( fpri_is_zero(b) ) {
        fpri_inv(e, a);
        fpri_zero(f);
    } else if ( fpri_is_zero(a) ) {
        fpri_zero(e);
        fpri_inv(f, b);
        fpri_neg(f, f);
    } else {
        fpri_sqr(e, a);
        fpri_sqr(f, b);
        fpri_add(f, e, f); /* f = a^2 + b^2 */
        fpri_div(e, a, f); /* e = a/(a^2 + b^2) */
        fpri_div(f, b, f); /* f = b/(a^2 + b^2) */
        fpri_neg(f, f);    /* f = -b/(a^2 + b^2) */
    }
#undef a
#undef b
#undef e
#undef f 
}

/* supports aliazing */
void fpci_inv( fpci_t dest, const fpci_t x ) {
    if (dest!=x){
        _fpci_inv   ( dest, x );
        return;
    }
    fpci_t temp;
    fpci_init(temp);
    _fpci_inv( temp, x );
    fpci_set(dest, temp);
    fpci_clear(temp);
}

void fpci_sqrabs   ( fpri_t dest, const fpci_t x ){
    fpri_t temp;
    fpri_sqr( dest, fpci_realref(x) );
    fpri_sqr( temp, fpci_imagref(x) );
    fpri_add( dest, dest, temp );
}

/* does not support aliazing */
void _fpci_mul   ( fpci_t z, const fpci_t x, const fpci_t y ){
#define a fpci_realref(x)
#define b fpci_imagref(x)
#define c fpci_realref(y)
#define d fpci_imagref(y)
#define e fpci_realref(z)
#define f fpci_imagref(z)
    if ( fpri_is_zero(b) ) {
        _fpri_mul(e, c, a);
        _fpri_mul(f, d, a);
    }
    else if ( fpri_is_zero(d) ) {
        _fpri_mul(e, a, c);
        _fpri_mul(f, b, c);
    }
    else if ( fpri_is_zero(a) ) {
        _fpri_mul(e, b, d);
        fpri_neg(e, e);
        _fpri_mul(f, b, c);
    }
    else if ( fpri_is_zero(c) ) {
        _fpri_mul(e, b, d);
        fpri_neg(e, e);
        _fpri_mul(f, a, d);
    }
    else {
        /* Gauss multiplication: e = ac-bd
                                 f = (a+b)(c+d) - ac - bd */
//         fpri_t u, v; /* no need to init */
//         fpri_add(u, a, b);
//         fpri_add(v, c, d);
//         _fpri_mul(f, u, v);
//         _fpri_mul(u, a, c);
//         _fpri_mul(v, b, d);
//         fpri_sub(f, f, u);
//         fpri_sub(f, f, v);
//         fpri_sub(e, u, v);
        /*              classical: e = ac - bd
                                 f = ad + bc */
        fpri_t u, v; /* no need to init */
        _fpri_mul(u, a, c);
        _fpri_mul(v, b, d);
        fpri_sub (e, u, v);
        _fpri_mul(u, a, d);
        _fpri_mul(v, b, c);
        fpri_add (f, u, v);
        /* use fma */
//         fpri_mul(e, b, d);
//         fpri_neg(e, e);
//         fpri_addmul(e, a, c);
//         fpri_mul(f, a, d);
//         fpri_addmul(f, b, c);
    }
#undef a
#undef b
#undef c
#undef d
#undef e
#undef f 
}
/* support aliazing */
void fpci_mul   ( fpci_t dest, const fpci_t x, const fpci_t y ){
    if (x==y) {
        fpci_sqr(dest, x);
    } else {
        if ((dest==x)||(dest==y)) {
            fpci_t temp;
            fpci_init(temp);
            _fpci_mul(temp, x, y);
            fpci_set(dest, temp);
            fpci_clear(temp);
        } else {
            _fpci_mul(dest, x, y);
        }
    }
}

/* does not support aliazing */
void _fpci_div     (fpci_t res, const fpci_t x, const fpci_t y){
     fpci_inv(res, y);
     fpci_mul(res, res, x);
}
/* support aliazing */
void fpci_div     (fpci_t res, const fpci_t x, const fpci_t y){
    if (x==y) {
        fpci_one(res);
        return;
    }
    if (res!=x) {
        _fpci_div ( res, x, y );
        return;
    }
    fpci_t temp;
    fpci_init(temp);
    _fpci_div ( temp, x, y );
    fpci_set(res, temp);
    fpci_clear(temp);
}

void fpci_pow_ui           (fpci_t res, const fpci_t x, ulong p){
    fpci_one(res);
    if (p==1)
        fpci_set(res, x);
    else if (p>1) {
        slong pp = p;
        fpci_t temp;
        fpci_init(temp);
        fpci_set(temp, x);
        while (pp>0){
            if (pp%2)
                fpci_mul(res, res, temp);
            pp = pp>>1;
            if (pp>0)
                fpci_sqr(temp, temp);
        }
        fpci_clear(temp);
    }
}

// void fpci_pow_ui           (fpci_t res, const fpci_t x, ulong p){
//     if (p==0)
//         fpci_one(res);
//     else if (p==1)
//         fpci_set(res, x);
//     else if (p==2)
//         fpci_sqr(res, x);
//     else if (p==3){
//         fpci_sqr(res, x);
//         fpci_mul(res, res, x);
//     } else if (p==4) {
//         fpci_sqr(res, x);
//         fpci_sqr(res, res);
//     } else {
//         fmpz_t pfmpz;
//         slong bits, i;
//         fmpz_init(pfmpz);
//         fmpz_set_ui(pfmpz, p);
//         fpci_set(res, x);
//         bits = fmpz_bits(pfmpz);
//         for (i = bits - 2; i >= 0; i--) {
//             fpci_sqr(res, res);
//             if (fmpz_tstbit(pfmpz, i))
//                 fpci_mul(res, res, x);
//         }
//         fmpz_clear(pfmpz);
//     }
// }

void fpci_fprint (FILE * file, const fpci_t x){
    fprintf(file, "( ");
    fpri_fprint_s (file, fpci_realref(x));
    fprintf(file, " +i ");
    fpri_fprint_s (file, fpci_imagref(x));
    fprintf(file, " )");
}
