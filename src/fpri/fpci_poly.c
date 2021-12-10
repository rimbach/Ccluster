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

#include "fpci_poly.h"

/* Memory management */
void fpci_poly_init(fpci_poly_t poly){
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void fpci_poly_init2(fpci_poly_t poly, slong len){
    fpci_poly_init(poly);
    fpci_poly_fit_length(poly, len);
}

void fpci_poly_clear(fpci_poly_t poly){
    fpri_free(poly->coeffs);
}

void fpci_poly_fit_length(fpci_poly_t poly, slong len){
    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = fpri_realloc(poly->coeffs, len * sizeof(fpci_struct));
        poly->alloc = len;
    }
}

void _fpci_poly_set_length(fpci_poly_t poly, slong len){
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            fpci_zero(poly->coeffs + i);
    }

    poly->length = len;
}

void _fpci_poly_normalise(fpci_poly_t poly){
    slong i;

    for (i = poly->length - 1;
        (i >= 0) && fpci_is_zero(poly->coeffs + i); i--);

    poly->length = i + 1;
}

void fpci_poly_set (fpci_poly_t dest, const fpci_poly_t src){
    slong i;
    slong len = fpci_poly_length(src);

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_acb_poly (fpci_poly_t dest, const acb_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_acb(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_fpri (fpci_poly_t dest, const fpri_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_fpri(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_arb_poly (fpci_poly_t dest, const arb_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_arb(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

/* evaluation */
void _fpci_poly_evaluate2_horner(fpci_t y, fpci_t z, fpci_srcptr poly, slong len, const fpci_t x) {
    if (len == 0) {
        fpci_zero(y);
        fpci_zero(z);
    }
    else if (len == 1) {
        fpci_set(y, poly + 0);
        fpci_zero(z);
    }
    else if (fpci_is_zero(x)) {
        fpci_set(y, poly + 0);
        fpci_set(z, poly + 1);
    }
    else if (len == 2) {
        fpci_mul(y, x, poly + 1);
        fpci_add(y, y, poly + 0);
        fpci_set(z, poly + 1);
    }
    else {
        fpci_t t, u, v;
        slong i;

        fpci_set(u, poly + len - 1);
        fpci_zero(v);

        for (i = len - 2; i >= 0; i--) {
            fpci_mul(t, v, x);
            fpci_add(v, u, t);
            fpci_mul(t, u, x);
            fpci_add(u, t, poly + i);
        }

        fpci_swap(y, u);
        fpci_swap(z, v);
    }
}
void fpci_poly_evaluate2_horner(fpci_t y, fpci_t z, const fpci_poly_t p, const fpci_t x) {
    _fpci_poly_evaluate2_horner(y, z, p->coeffs, p->length, x);
}
            
/* printing */
void fpci_poly_fprint (FILE * file, const fpci_poly_t x){
    slong i;
    slong len = x->length;
    if (len==0)
        fprintf(file, "degree: 0, coeffs: 0");
    else {
        fprintf(file, "degree: %ld, coeffs: ", len-1);
        for (i = 0; i <len; i++){ 
            fpci_fprint(file, x->coeffs+i);
            if (i<len-1)
                fprintf(file, "\n");
        }
    }
}
