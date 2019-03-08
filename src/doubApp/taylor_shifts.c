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

/* assume len>0 and p->alloc >= len +1 */
void _doubCompApp_poly_timesRXPC_inplace (doubCompApp_ptr p, const doubCompApp_t c, const doubRealApp_t r, slong len){
    slong i;
    doubCompApp_t t;
    doubCompApp_mul_doubRealApp(p+len, p+len-1, r);
    for (i = len-1; i>0; i--) {
        doubCompApp_mul(p+i, p+i, c);
        doubCompApp_mul_doubRealApp(t, p+i-1, r);
        doubCompApp_add(p+i, p+i, t);
    }
    doubCompApp_mul(p, p, c);
}

void _doubCompApp_poly_taylor_shift_horner( doubCompApp_ptr res, doubCompApp_srcptr src, const doubCompApp_t c, const doubRealApp_t r, slong len){
    
    slong i;
    doubCompApp_set(res, src + len-1);
    for( i=1; i<len; i++) {
        _doubCompApp_poly_timesRXPC_inplace (res, c, r, i);
        doubCompApp_add(res, res, src + (len-1) -i);
    }
    
}

#define TS_CUTOFF 32
void _doubCompApp_poly_taylor_shift_DQ( doubCompApp_ptr res, doubCompApp_poly_ptr pows, 
                                        doubCompApp_srcptr src, const doubCompApp_t c, const doubRealApp_t r, slong len){
    slong i;
    if ( len <= TS_CUTOFF ){
//         _doubCompApp_poly_taylor_shift_horner( res, src, c, r, len);
        /*copy src in res*/
        for (i=0; i<len; i++)
            doubCompApp_set(res + i, src +i);
        _doubCompApp_poly_taylor_shift_convolution( res, c, len);
        return;
    }
    
    slong cut = len/2;
    int pow = 0x1;
    int log = 0;
//     slong i;
    while (pow < cut){
        pow = pow<<1;
        log++;
    }
//     printf("DQ: len: %ld, cut: %ld, log: %d\n", len, cut, log);
    doubCompApp_poly_t t1, t2, t3;
    doubCompApp_poly_init2(t1, cut);
    doubCompApp_poly_init2(t2, cut);
    doubCompApp_poly_init2(t3, len);
    for (i=0;i<len;i++) {
//         if (i<cut){
//             doubCompApp_zero(t1->coeffs+i);
//             doubCompApp_zero(t2->coeffs+i);
//         }
        doubCompApp_zero(t3->coeffs+i);
    }
    
//     printf("src: \n"); 
//     for (i=0; i<len; i++){
//         printf("deg %ld: ",i);
//         doubCompApp_print(src + i);
//         printf("\n");
//     }
//     printf("\n");
    
    _doubCompApp_poly_taylor_shift_DQ( t1->coeffs, pows, src, c, r, cut);
    _doubCompApp_poly_taylor_shift_DQ( t2->coeffs, pows, src+cut, c, r, cut);
    
//     printf("powers: \n"); doubCompApp_poly_print(pows+log); printf("\n\n");
//     _doubCompApp_poly_mullow_classical(t3->coeffs, t2->coeffs, cut, (pows+log)->coeffs, (pows+log)->length, len);
    _doubCompApp_poly_set_length(t3, len);
    _doubCompApp_poly_set_length(t2, cut);
    doubCompApp_poly_mul_karatsuba(t3, t2, pows+log );
    _doubCompApp_poly_add( res, t1->coeffs, cut, t3->coeffs, len, len );
    
//     printf("t1: \n"); 
//     for (i=0; i<cut; i++){
//         printf("deg %ld: ",i);
//         doubCompApp_print(t1->coeffs + i);
//         printf("\n");
//     }
//     printf("\n");
//     
//     printf("t2: \n"); 
//     for (i=0; i<cut; i++){
//         printf("deg %ld: ",i);
//         doubCompApp_print(t2->coeffs + i);
//         printf("\n");
//     }
//     printf("\n");
//     
//     printf("res: \n"); 
//     for (i=0; i<len; i++){
//         printf("deg %ld: ",i);
//         doubCompApp_print(res + i);
//         printf("\n");
//     }
//     printf("\n");
    
    doubCompApp_poly_clear(t1);
    doubCompApp_poly_clear(t2);
    doubCompApp_poly_clear(t3);
}

void doubCompApp_poly_taylor_shift_DQ( doubCompApp_poly_t res, doubCompApp_poly_t f, const doubCompApp_t c, const doubRealApp_t r){
    
    int nlen = 0x1;
    int log = 0;
    slong i;
    doubCompApp_poly_t nf;
    doubCompApp_poly_init(nf);
    doubCompApp_poly_set(nf, f);
    
    /* find the first power of 2 greater than len of f */
    while ( nlen < f->length ) {
        nlen = nlen<<1;
        log++;
    }
//     printf("nlen: %d, log: %d\n", nlen, log);
    
    /* copy f */
    doubCompApp_poly_fit_length(res, nlen);
    doubCompApp_poly_fit_length(nf, nlen);
    _doubCompApp_poly_set_length(nf, nlen);
    
    /* pad nf with 0 */
    for (i=f->length; i<nlen; i++)
        doubCompApp_zero(nf->coeffs + i);
    
    /* compute the powers of (rx+c) */
    doubCompApp_poly_ptr powers = NULL;
    powers = flint_realloc(powers, log * sizeof(doubCompApp_poly));
    doubCompApp_poly_init2(powers,2);
    doubCompApp_set(powers->coeffs, c);
    doubCompApp_set_doubRealApp(powers->coeffs+1, r);
    _doubCompApp_poly_set_length(powers, 2);
    for (i=1; i<log; i++){
        doubCompApp_poly_init(powers + i);
//         doubCompApp_poly_mul_classical(powers + i, powers + (i-1), powers + (i-1));
        doubCompApp_poly_sqr_karatsuba( powers+i, powers + (i-1));   
    }
//     printf("powers: \n");
//     for (i=0; i<log; i++) {
//         printf("--i: %ld \n", i); doubCompApp_poly_print(powers+i); printf("\n");
//     }
//     
//     printf("++  f: "); doubCompApp_poly_print(f); printf("\n");
//     printf("++ nf: "); doubCompApp_poly_print(nf); printf("\n");
    
    /* call the main function */
    _doubCompApp_poly_taylor_shift_DQ( res->coeffs, powers, nf->coeffs, c, r, nlen);
    /* clear powers */
    for (i=0; i<log; i++)
        doubCompApp_poly_clear(powers + i);
    flint_free(powers);
    
    _doubCompApp_poly_set_length(res, f->length);
    _doubCompApp_poly_normalise(res);
    doubCompApp_poly_clear(nf);
}

void doubCompApp_poly_taylor_shift_horner_inplace( doubCompApp_poly_t f, const doubCompApp_t c, const doubRealApp_t r){
    doubCompApp_poly_t t;
    doubCompApp_poly_init2(t, f->length);
    
    _doubCompApp_poly_taylor_shift_horner( t->coeffs, f->coeffs, c, r, f->length);
    
    _doubCompApp_poly_set_length(t, f->length);
    _doubCompApp_poly_normalise(t);
    
    doubCompApp_poly_swap(f,t);
    doubCompApp_poly_clear(t);
    
}

void doubCompApp_poly_taylor_shift_horner( doubCompApp_poly_t res, const doubCompApp_poly_t f, const doubCompApp_t c, const doubRealApp_t r){
    
    doubCompApp_poly_fit_length(res, f->length);
    _doubCompApp_poly_taylor_shift_horner( res->coeffs, f->coeffs, c, r, f->length);
    _doubCompApp_poly_set_length(res, f->length);
    _doubCompApp_poly_normalise(res);
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

void doubCompApp_poly_taylor_shift_convolution( doubCompApp_poly_t res, const doubCompApp_poly_t f, const doubCompApp_t c, const doubRealApp_t r){
    
    doubCompApp_poly_set(res, f);
    _doubCompApp_poly_taylor_shift_convolution( res->coeffs, c, res->length);
    _doubCompApp_poly_set_length(res, f->length);
    _doubCompApp_poly_normalise(res);
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
