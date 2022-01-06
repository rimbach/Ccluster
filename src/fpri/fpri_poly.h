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

#ifndef FPRI_POLY_H
#define FPRI_POLY_H

#ifdef FPRI_INLINE_C
#define FPRI_INLINE
#else
#define FPRI_INLINE static __inline__
#endif

#include "arb_poly.h"
#include "fpri.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    fpri_ptr coeffs;
    slong length;
    slong alloc;
} fpri_poly_struct;

typedef fpri_poly_struct fpri_poly_t[1];

/* Memory management */
            void fpri_poly_init(fpri_poly_t poly);
            void fpri_poly_init2(fpri_poly_t poly, slong len);
            void fpri_poly_clear(fpri_poly_t poly);
            void fpri_poly_fit_length(fpri_poly_t poly, slong len);
            void _fpri_poly_set_length(fpri_poly_t poly, slong len);
            
            
            void _fpri_poly_normalise(fpri_poly_t poly);

FPRI_INLINE void fpri_poly_swap(fpri_poly_t poly1, fpri_poly_t poly2) {
    fpri_poly_struct t= *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

FPRI_INLINE slong fpri_poly_length(const fpri_poly_t poly) {
    return poly->length;
}

FPRI_INLINE slong fpri_poly_degree(const fpri_poly_t poly) {
    return poly->length - 1;
}

/* setting */
FPRI_INLINE void fpri_poly_zero(fpri_poly_t poly){
    poly->length = 0;
}

FPRI_INLINE void fpri_poly_one (fpri_poly_t poly){
    fpri_poly_fit_length(poly, 1);
    fpri_one(poly->coeffs);
    _fpri_poly_set_length(poly, 1);
}

            void fpri_poly_set (fpri_poly_t dest, const fpri_poly_t src);
            void fpri_poly_set_arb_poly (fpri_poly_t dest, const arb_poly_t src);

/* evaluation */
            void fpri_poly_evaluate2_horner(fpri_t y, fpri_t z, const fpri_poly_t p, const fpri_t x); 
/* printing */
            void fpri_poly_fprint (FILE * file, const fpri_poly_t x);
FPRI_INLINE void fpri_poly_print (const fpri_poly_t x) { fpri_poly_fprint(stdout, x); }

#ifdef __cplusplus
}
#endif

#endif
