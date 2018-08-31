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
        if (! doubCompApp_is_zero(x->coeffs+i)) {
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

void _doubCompApp_poly_taylor_shift_convolution(doubCompApp_ptr p, const doubCompApp_t c, slong len) {
    slong i, n = len - 1;
    doubCompApp_t d;
    doubRealApp_t f;
    doubCompApp_ptr t, u;

    if (doubCompApp_is_zero(c) || len <= 1)
        return;

    t = _doubCompApp_vec_init(len);
//     u = _doubCompApp_vec_init(len);

//     doubRealApp_init(f);
//     doubCompApp_init(d);

    doubRealApp_one(f);
    for (i = 2; i <= n; i++) {
        doubRealApp_mul_ui(f, f, i);
        doubCompApp_mul_doubRealApp(p + i, p + i, f);
    }

    _doubCompApp_poly_reverse(p, p, len, len);

    doubCompApp_one(t + n);
    for (i = n; i > 0; i--)
        doubCompApp_mul_ui(t + i - 1, t + i, i);

    if (doubCompApp_equal_si(c, -1)) {
        for (i = 1; i <= n; i += 2)
            doubCompApp_neg(t + i, t + i);
    }
    else if (!doubCompApp_is_one(c)){
        doubCompApp_set(d, c);

        for (i = 1; i <= n; i++){
            doubCompApp_mul(t + i, t + i, d);
            doubCompApp_mul(d, d, c);
        }
    }

//     doubCompApp_ptr u2;
//     u2 = _doubCompApp_vec_init(len);
//     for (i = 0; i < len; i++)
//         doubCompApp_zero(u2+i);
//     _doubCompApp_poly_mullow_classical(u2, p, len, t, len, len);
//     printf("classisal: \n"); 
//     for (int i=0; i<len; i++){
//         printf("deg %d, coeff: ",i);
//         doubCompApp_print(u2 + i);
//         printf("\n");
//     }
    
    u = _doubCompApp_vec_init(2*len-1);
    for (i = 0; i < 2*len-1; i++)
        doubCompApp_zero(u+i);
//     _doubCompApp_poly_mullow_karatsuba(u, p, len, t, len);
    _doubCompApp_poly_mullow_classical(u, p, len, t, len, len);
//     printf("\n karatsuba: \n"); 
//     for (int i=0; i<len; i++){
//         printf("deg %d, coeff: ",i);
//         doubCompApp_print(u + i);
//         printf("\n");
//     }
    

    doubRealApp_sqr(f, f);

//     if (doubRealApp_bits(f) > 0.25 * prec)
//     {
//         doubRealApp_inv(f, f, prec);
//     }
//     else
//     {
        for (i = 0; i <= n; i++)
            doubCompApp_div_doubRealApp(u + i, u + i, f);

        doubRealApp_one(f);
//     }

    for (i = n; i >= 0; i--) {
        doubCompApp_mul_doubRealApp(p + i, u + n - i, f);
        doubRealApp_mul_ui(f, f, (i == 0) ? 1 : i);
    }

    _doubCompApp_vec_clear(t, len);
    _doubCompApp_vec_clear(u, len);

    doubRealApp_clear(f);
    doubCompApp_clear(d);
}

// void doubCompApp_poly_taylorShift_in_place( doubCompApp_poly_t f, const compRat_t center, const realRat_t radius){
//     
//     doubCompApp_t c;
//     doubCompApp_init(c);
//     /* first convert center in compApp */
//     compApp_t cc;
//     compApp_init(cc);
//     compApp_set_compRat( cc, center, 53 );
//     doubCompApp_set_compApp(c, cc);
//      
//     doubCompApp_ptr fptr = f->coeffs;
//     const slong len  = f->length;
// //     _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
// 
// //     compApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
//      
//     doubCompApp_clear(c);
//     compApp_clear(cc);
// }


