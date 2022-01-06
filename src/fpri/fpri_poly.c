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

#include "fpri_poly.h"

/* Memory management */
void fpri_poly_init(fpri_poly_t poly){
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void fpri_poly_init2(fpri_poly_t poly, slong len){
    fpri_poly_init(poly);
    fpri_poly_fit_length(poly, len);
}

void fpri_poly_clear(fpri_poly_t poly){
    fpri_free(poly->coeffs);
}

void fpri_poly_fit_length(fpri_poly_t poly, slong len){
    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = fpri_realloc(poly->coeffs, len * sizeof(fpri_struct));
        poly->alloc = len;
    }
}

void _fpri_poly_set_length(fpri_poly_t poly, slong len){
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            fpri_zero(poly->coeffs + i);
    }

    poly->length = len;
}

void _fpri_poly_normalise(fpri_poly_t poly){
    slong i;

    for (i = poly->length - 1;
        (i >= 0) && fpri_is_zero(poly->coeffs + i); i--);

    poly->length = i + 1;
}

void fpri_poly_set (fpri_poly_t dest, const fpri_poly_t src){
    slong i;
    slong len = fpri_poly_length(src);

    fpri_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpri_set(dest->coeffs + i, src->coeffs + i);
    
    _fpri_poly_set_length(dest, len);
}

void fpri_poly_set_arb_poly (fpri_poly_t dest, const arb_poly_t src){
    slong i;
    slong len = src->length;

    fpri_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpri_set_arb(dest->coeffs + i, src->coeffs + i);
    
    _fpri_poly_set_length(dest, len);
}

/* evaluation */
void _fpri_poly_evaluate2_horner(fpri_t y, fpri_t z, fpri_srcptr poly, slong len, const fpri_t x) {
    if (len == 0) {
        fpri_zero(y);
        fpri_zero(z);
    }
    else if (len == 1) {
        fpri_set(y, poly + 0);
        fpri_zero(z);
    }
    else if (fpri_is_zero(x)) {
        fpri_set(y, poly + 0);
        fpri_set(z, poly + 1);
    }
    else if (len == 2) {
        fpri_mul(y, x, poly + 1);
        fpri_add(y, y, poly + 0);
        fpri_set(z, poly + 1);
    }
    else {
        fpri_t t, u, v;
        slong i;

        fpri_set(u, poly + len - 1);
        fpri_zero(v);

        for (i = len - 2; i >= 0; i--) {
            fpri_mul(t, v, x);
            fpri_add(v, u, t);
            fpri_mul(t, u, x);
            fpri_add(u, t, poly + i);
        }

        fpri_swap(y, u);
        fpri_swap(z, v);
    }
}
void fpri_poly_evaluate2_horner(fpri_t y, fpri_t z, const fpri_poly_t p, const fpri_t x) {
    _fpri_poly_evaluate2_horner(y, z, p->coeffs, p->length, x);
}
            
/* printing */
void fpri_poly_fprint (FILE * file, const fpri_poly_t x){
    slong i;
    slong len = x->length;
    if (len==0)
        fprintf(file, "degree: 0, coeffs: 0");
    else {
        fprintf(file, "degree: %ld, coeffs: ", len-1);
        for (i = 0; i <len; i++){ 
            fpri_fprint(file, x->coeffs+i);
            if (i<len-1)
                fprintf(file, "\n");
        }
    }
}
