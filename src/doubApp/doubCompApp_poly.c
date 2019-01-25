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

#include "doubCompApp_poly.h"

void doubCompApp_poly_init(doubCompApp_poly_t poly){
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void doubCompApp_poly_init2(doubCompApp_poly_t poly, slong len){
    doubCompApp_poly_init(poly);
    doubCompApp_poly_fit_length(poly, len);
}

void doubCompApp_poly_clear(doubCompApp_poly_t poly){
//     slong i;

//     for (i = 0; i < poly->alloc; i++)
//         doubCompApp_clear(poly->coeffs + i);

    flint_free(poly->coeffs);
}

void doubCompApp_poly_fit_length(doubCompApp_poly_t poly, slong len){
//     slong i;

    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = flint_realloc(poly->coeffs,
            len * sizeof(doubCompApp));
//         for (i = poly->alloc; i < len; i++)
//             doubCompApp_init(poly->coeffs + i);

        poly->alloc = len;
    }
}

void _doubCompApp_poly_set_length(doubCompApp_poly_t poly, slong len){
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            doubCompApp_zero(poly->coeffs + i);
    }

    poly->length = len;
}

void _doubCompApp_poly_normalise(doubCompApp_poly_t poly){
    slong i;

    for (i = poly->length - 1;
        (i >= 0) && doubCompApp_is_zero(poly->coeffs + i); i--);

    poly->length = i + 1;
}

void doubCompApp_poly_set (doubCompApp_poly_t dest, const doubCompApp_poly_t src){
    slong i;
    slong len = doubCompApp_poly_length(src);

    doubCompApp_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        doubCompApp_set(dest->coeffs + i, src->coeffs + i);
    
    _doubCompApp_poly_set_length(dest, len);
}

void doubCompApp_poly_set_compApp_poly (doubCompApp_poly_t dest, const compApp_poly_t src){
    slong i;
    slong len = src->length;

    doubCompApp_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        doubCompApp_set_compApp(dest->coeffs + i, src->coeffs + i);
    
    _doubCompApp_poly_set_length(dest, len);
}

void doubCompApp_poly_fprint (FILE * file, const doubCompApp_poly_t x){
    slong i;
    slong len = doubCompApp_poly_length(x);
    int OK = 0;
    if (len==0) {
        fprintf(file, "deg %5d: 0\n",0);
        return;
    }
    for (i = 0; i <len; i++){
        if ( (!doubCompApp_is_zero(x->coeffs+i))||(i==len-1) ) {
            OK=1;
            fprintf(file, "deg %5ld: ",i);
            doubCompApp_fprint(file, x->coeffs+i);
            fprintf(file, "\n");
        }
    }
    if (!OK){
        fprintf(file, "deg %5d: ",0);
        doubCompApp_fprint(file, x->coeffs);
        fprintf(file, "\n");
    }
}

void doubCompApp_poly_neg (doubCompApp_poly_t y, const doubCompApp_poly_t x){
    slong i;
    slong len = doubCompApp_poly_length(x);

    doubCompApp_poly_fit_length(y, len);
    
    for (i = 0; i < len; i++)
        doubCompApp_neg(y->coeffs + i, x->coeffs + i);
    
    _doubCompApp_poly_set_length(y, len);
}
/*assume z has length >= len
 */
void _doubCompApp_poly_add( doubCompApp_ptr z, doubCompApp_srcptr x, slong lenx, 
                                               doubCompApp_srcptr y, slong leny, slong len){
    slong i;
    slong lenM = CCLUSTER_MAX(lenx, leny);
    slong lenm = CCLUSTER_MIN(lenx, leny);
    
    for (i = 0; i < len && i<lenm; i++)
        doubCompApp_add( z + i, x + i, y + i );
    
    if (lenm==leny)
        for (i = lenm; i < len && i<lenM; i++)
            doubCompApp_set( z + i, x + i);
    else 
        for (i = lenm; i < len && i<lenM; i++)
            doubCompApp_set( z + i, y + i);
}

void doubCompApp_poly_add( doubCompApp_poly_t z, const doubCompApp_poly_t x, const doubCompApp_poly_t y){
    
    slong len = CCLUSTER_MAX(x->length, y->length);
    doubCompApp_poly_fit_length(z,len);
    _doubCompApp_poly_add( z->coeffs, x->coeffs, x->length, y->coeffs, y->length, len);
    _doubCompApp_poly_set_length(z, len);
}

/*assume z has length >= len
 */
void _doubCompApp_poly_sub( doubCompApp_ptr z, doubCompApp_srcptr x, slong lenx, 
                                               doubCompApp_srcptr y, slong leny, slong len){
    slong i;
    slong lenM = CCLUSTER_MAX(lenx, leny);
    slong lenm = CCLUSTER_MIN(lenx, leny);
    
    for (i = 0; i < len && i<lenm; i++)
        doubCompApp_sub( z + i, x + i, y + i );
    
    if (lenm==leny)
        for (i = lenm; i < len && i<lenM; i++)
            doubCompApp_set( z + i, x + i);
    else 
        for (i = lenm; i < len && i<lenM; i++)
            doubCompApp_neg( z + i, y + i);
}

void doubCompApp_poly_sub( doubCompApp_poly_t z, const doubCompApp_poly_t x, const doubCompApp_poly_t y){
    
    slong len = CCLUSTER_MAX(x->length, y->length);
    doubCompApp_poly_fit_length(z,len);
    _doubCompApp_poly_sub( z->coeffs, x->coeffs, x->length, y->coeffs, y->length, len);
    _doubCompApp_poly_set_length(z, len);
}

void _doubCompApp_poly_shift_left(doubCompApp_ptr res, doubCompApp_srcptr poly, slong len, slong n){
    slong i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly) {
        for (i = len; i--; )
            doubCompApp_set(res + n + i, poly + i);
    }
    else {
        for (i = len; i--; )
            doubCompApp_swap(res + n + i, res + i);
    }
    for (i = 0; i < n; i++)
        doubCompApp_zero(res + i);
}

void doubCompApp_poly_shift_left( doubCompApp_poly_t res, const doubCompApp_poly_t poly, slong n){
    if (n == 0) {
        doubCompApp_poly_set(res, poly);
        return;
    }
    if (poly->length == 0) {
        doubCompApp_poly_zero(res);
        return;
    }

    doubCompApp_poly_fit_length(res, poly->length + n);
    _doubCompApp_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n);
    _doubCompApp_poly_set_length(res, poly->length + n);
}

void _doubCompApp_poly_reverse(doubCompApp_ptr res, doubCompApp_srcptr poly, slong len, slong n) {
    if (res == poly) {
        slong i;

        for (i = 0; i < n / 2; i++) {
            doubCompApp t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            doubCompApp_zero(res + i);
    }
    else {
        slong i;

        for (i = 0; i < n - len; i++)
            doubCompApp_zero(res + i);

        for (i = 0; i < len; i++)
            doubCompApp_set(res + (n - len) + i, poly + (len - 1) - i);
    }
}

void doubCompApp_poly_oneGraeffeIteration_in_place( doubCompApp_poly_t f ){
    
    doubCompApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    doubCompApp_poly_t fe, fo;
    doubCompApp_poly_init2(fe, len2);
    doubCompApp_poly_init2(fo, len2);
    doubCompApp_ptr feptr = fe->coeffs;
    doubCompApp_ptr foptr = fo->coeffs;
    doubCompApp_zero( feptr + len2-1);
    doubCompApp_zero( foptr + len2-1);
    
    for (i = 0; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) 
            doubCompApp_set( feptr + quo, fptr+i);
        else
            doubCompApp_set( foptr + quo, fptr+i);
    }
    _doubCompApp_poly_set_length(fe, len2);
    _doubCompApp_poly_set_length(fo, len2);
//     printf("fe \n"); doubCompApp_poly_print(fe); printf("\n\n");
//     printf("fo \n"); doubCompApp_poly_print(fo); printf("\n\n");
    
    doubCompApp_poly_t fes, fos;
    doubCompApp_poly_init2(fes, len1);
    doubCompApp_poly_init2(fos, len1);
//     doubCompApp_poly_sqr_karatsuba( fes, fe);
//     doubCompApp_poly_sqr_karatsuba( fos, fo);
    doubCompApp_poly_mul_classical( fes, fe, fe);
    doubCompApp_poly_mul_classical( fos, fo, fo);
//     printf("fes \n"); doubCompApp_poly_print(fes); printf("\n\n");
//     printf("fos \n"); doubCompApp_poly_print(fos); printf("\n\n");
    doubCompApp_poly_shift_left( fos, fos, 1 );
    doubCompApp_poly_sub(f, fes, fos);
    if ((len1%2)==0)
        doubCompApp_poly_neg(f, f);
    
    doubCompApp_poly_clear(fe);
    doubCompApp_poly_clear(fo);
    doubCompApp_poly_clear(fes);
    doubCompApp_poly_clear(fos);
    
}


