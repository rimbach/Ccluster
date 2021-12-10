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

#ifndef FPCI_POLY_H
#define FPCI_POLY_H

#ifdef FPRI_INLINE_C
#define FPRI_INLINE
#else
#define FPRI_INLINE static __inline__
#endif

#include "acb_poly.h"
#include "fpci.h"
#include "fpri_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    fpci_ptr coeffs;
    slong length;
    slong alloc;
} fpci_poly_struct;

typedef fpci_poly_struct fpci_poly_t[1];

/* Memory management */
            void fpci_poly_init(fpci_poly_t poly);
            void fpci_poly_init2(fpci_poly_t poly, slong len);
            void fpci_poly_clear(fpci_poly_t poly);
            void fpci_poly_fit_length(fpci_poly_t poly, slong len);
            void _fpci_poly_set_length(fpci_poly_t poly, slong len);
            
            void _fpci_poly_normalise(fpci_poly_t poly);

FPRI_INLINE void fpci_poly_swap(fpci_poly_t poly1, fpci_poly_t poly2) {
    fpci_poly_struct t= *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

FPRI_INLINE slong fpci_poly_length(const fpci_poly_t poly) {
    return poly->length;
}

FPRI_INLINE slong fpci_poly_degree(const fpci_poly_t poly) {
    return poly->length - 1;
}

/* setting */
FPRI_INLINE void fpci_poly_zero(fpci_poly_t poly){
    poly->length = 0;
}

FPRI_INLINE void fpci_poly_one (fpci_poly_t poly){
    fpci_poly_fit_length(poly, 1);
    fpci_one(poly->coeffs);
    _fpci_poly_set_length(poly, 1);
}

FPRI_INLINE void fpci_poly_onei (fpci_poly_t poly){
    fpci_poly_fit_length(poly, 1);
    fpci_onei(poly->coeffs);
    _fpci_poly_set_length(poly, 1);
}

            void fpci_poly_set (fpci_poly_t dest, const fpci_poly_t src);
            void fpci_poly_set_acb_poly (fpci_poly_t dest, const acb_poly_t src);
            void fpci_poly_set_fpri (fpci_poly_t dest, const fpri_poly_t src);
            void fpci_poly_set_arb_poly (fpci_poly_t dest, const arb_poly_t src);

/* evaluation */
            void fpci_poly_evaluate2_horner(fpci_t y, fpci_t z, const fpci_poly_t p, const fpci_t x); 
/* printing */
            void fpci_poly_fprint (FILE * file, const fpci_poly_t x);
FPRI_INLINE void fpci_poly_print (const fpci_poly_t x) { fpci_poly_fprint(stdout, x); }

#ifdef __cplusplus
}
#endif

#endif
