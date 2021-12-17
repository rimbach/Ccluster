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

#ifndef FPCI_H
#define FPCI_H

#ifdef FPRI_INLINE_C
#define FPRI_INLINE
#else
#define FPRI_INLINE static __inline__
#endif

#include "acb.h"
#include "fpri.h"


#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct
{
    fpri_struct real;
    fpri_struct imag;
}
fpci_struct;

typedef fpci_struct fpci_t[1];
typedef fpci_struct * fpci_ptr;
typedef const fpci_struct * fpci_srcptr;

#define fpci_realref(x) (&(x)->real)
#define fpci_imagref(x) (&(x)->imag)

/* memory managment */
FPRI_INLINE void fpci_init(fpci_t x) {
    fpri_init(fpci_realref(x));
    fpri_init(fpci_imagref(x));
}

FPRI_INLINE void fpci_clear(fpci_t x){
    fpri_clear(fpci_realref(x));
    fpri_clear(fpci_imagref(x));
}

FPRI_INLINE fpci_ptr _fpci_vec_init( slong n ) {
    return (fpci_ptr) fpri_malloc ( n*sizeof(fpci_struct) );
}
FPRI_INLINE void     _fpci_vec_clear( fpci_ptr v, slong n ) {
    fpri_free(v);
}

FPRI_INLINE void fpci_swap(fpci_t z, fpci_t x) {
    fpri_swap(fpci_realref(z), fpci_realref(x));
    fpri_swap(fpci_imagref(z), fpci_imagref(x));
}

/* setting */
FPRI_INLINE void fpci_zero(fpci_t z) {
    fpri_zero(fpci_realref(z));
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE void fpci_vec_zero ( fpci_ptr v, slong n ) {
    slong i;
    for (i=0; i<n; i++)
        fpci_zero(v + i);
}

FPRI_INLINE void fpci_one(fpci_t z) {
    fpri_one(fpci_realref(z));
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE void fpci_onei(fpci_t z) {
    fpri_zero(fpci_realref(z));
    fpri_one(fpci_imagref(z));
}

FPRI_INLINE void fpci_set(fpci_t z, const fpci_t x) {
    fpri_set(fpci_realref(z), fpci_realref(x));
    fpri_set(fpci_imagref(z), fpci_imagref(x));
}

FPRI_INLINE void fpci_set_d(fpci_t z, const double x) {
    fpri_set_d(fpci_realref(z), x);
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE void fpci_set_d_d(fpci_t z, const double x, const double y) {
    fpri_set_d(fpci_realref(z), x);
    fpri_set_d(fpci_imagref(z), y);
}

FPRI_INLINE void fpci_set_si(fpci_t z, const slong x) {
    fpri_set_si(fpci_realref(z), x);
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE void fpci_set_si_si(fpci_t z, const slong x, const slong y) {
    fpri_set_si(fpci_realref(z), x);
    fpri_set_si(fpci_imagref(z), y);
}

FPRI_INLINE void fpci_set_fpri(fpci_t z, const fpri_t x) {
    fpri_set(fpci_realref(z), x);
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE void fpci_set_acb(fpci_t z, const acb_t x) {
    fpri_set_arb(fpci_realref(z), acb_realref(x));
    fpri_set_arb(fpci_imagref(z), acb_imagref(x));
}

FPRI_INLINE int  fpci_get_acb(acb_t x, const fpci_t z) {
    return ( fpri_get_arb(acb_realref(x), fpci_realref(z) ) &&
             fpri_get_arb(acb_imagref(x), fpci_imagref(z) ) );
}

FPRI_INLINE void fpci_set_arb(fpci_t z, const arb_t x) {
    fpri_set_arb(fpci_realref(z), x);
    fpri_zero(fpci_imagref(z));
}

FPRI_INLINE int fpci_is_zero(const fpci_t z) {
    return ( fpri_is_zero(fpci_realref(z)) && fpri_is_zero(fpci_imagref(z)) );
}

FPRI_INLINE int fpci_contains_zero(const fpci_t z) {
    return ( fpri_contains_zero(fpci_realref(z)) && fpri_contains_zero(fpci_imagref(z)) );
}

/* arithmetic */
FPRI_INLINE void fpci_neg   ( fpci_t dest, const fpci_t x ) { 
    fpri_neg(fpci_realref(dest), fpci_realref(x));
    fpri_neg(fpci_imagref(dest), fpci_imagref(x)); 
}

/* does not support aliazing */
            void _fpci_sqr   ( fpci_t dest, const fpci_t x );
/* support aliazing */
            void fpci_sqr   ( fpci_t dest, const fpci_t x );
/* does not support aliazing */
            void _fpci_inv   ( fpci_t dest, const fpci_t x );
/* support aliazing */
            void fpci_inv   ( fpci_t dest, const fpci_t x );
            /* abs is not implemented because fpri_sqrt is not implemented */
            /* instead compute the square of the abs */
            void fpci_sqrabs   ( fpri_t dest, const fpci_t x );
    

FPRI_INLINE void fpci_sub   ( fpci_t dest, const fpci_t x, const fpci_t y) { 
    fpri_sub(fpci_realref(dest), fpci_realref(x), fpci_realref(y));
    fpri_sub(fpci_imagref(dest), fpci_imagref(x), fpci_imagref(y));
}

FPRI_INLINE void fpci_add   ( fpci_t dest, const fpci_t x, const fpci_t y) { 
    fpri_add(fpci_realref(dest), fpci_realref(x), fpci_realref(y));
    fpri_add(fpci_imagref(dest), fpci_imagref(x), fpci_imagref(y));
}

/* does not support aliazing */
            void _fpci_mul   ( fpci_t dest, const fpci_t x, const fpci_t y );
/* support aliazing */
            void fpci_mul   ( fpci_t dest, const fpci_t x, const fpci_t y );
            
FPRI_INLINE void fpci_mul_fpri ( fpci_t dest, const fpci_t x, const fpri_t y) { 
    fpri_mul(fpci_realref(dest), fpci_realref(x), y);
    fpri_mul(fpci_imagref(dest), fpci_imagref(x), y);
}

FPRI_INLINE void fpci_mul_si   ( fpci_t dest, const fpci_t x, const slong y) { 
    fpri_mul_si(fpci_realref(dest), fpci_realref(x), y);
    fpri_mul_si(fpci_imagref(dest), fpci_imagref(x), y);
}

/* does not support aliazing */
void _fpci_div     (fpci_t res, const fpci_t x, const fpci_t y);
/* support aliazing */
void fpci_div     (fpci_t res, const fpci_t x, const fpci_t y);

FPRI_INLINE void fpci_div_fpri ( fpci_t dest, const fpci_t x, const fpri_t y) { 
    fpri_div(fpci_realref(dest), fpci_realref(x), y);
    fpri_div(fpci_imagref(dest), fpci_imagref(x), y);
}

FPRI_INLINE void fpci_div_si   ( fpci_t dest, const fpci_t x, const slong y) { 
    fpri_div_si(fpci_realref(dest), fpci_realref(x), y);
    fpri_div_si(fpci_imagref(dest), fpci_imagref(x), y);
}

void fpci_pow_ui           (fpci_t res, const fpci_t x, ulong p);

            void fpci_fprint (FILE * file, const fpci_t x);
FPRI_INLINE void fpci_print (const fpci_t x) { fpci_fprint(stdout, x); }

#ifdef __cplusplus
}
#endif

#endif
