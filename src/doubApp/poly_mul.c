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

/* inplace multiplication by an integer */
void _doubCompApp_poly_mul_si(doubCompApp_ptr res, slong lenres, slong s) {
    slong i;
    for (i=0; i<lenres; i++)
        doubCompApp_mul_si(res + i, res + i, s);
}

void doubCompApp_poly_mul_si( doubCompApp_poly_t res, const doubCompApp_poly_t x, slong s){
    
    if (res==x)
        _doubCompApp_poly_mul_si(res->coeffs, res->length, s);
    else {
        doubCompApp_poly_set(res, x);
        _doubCompApp_poly_mul_si(res->coeffs, res->length, s);
    }
        
}

void _doubCompApp_poly_mullow_classical(doubCompApp_ptr res,
    doubCompApp_srcptr x, slong lenx,
    doubCompApp_srcptr y, slong leny, slong len ){
    
//     slong lenres = lenx + leny -1;
    slong i,j;
    doubCompApp_t temp;
    
//     /* initialize the coefficients of res to 0*/
//     for (i = 0; i<lenres; i++)
//         doubCompApp_zero( res + i );

    for (i=0; i<lenx && i<len; i++)
        if (!doubCompApp_is_zero( x + i )){
            for (j=0; j<leny && (i+j)<len; j++) {
                doubCompApp_mul(temp, x+i, y+j);
                doubCompApp_add(res + i + j, res + i + j, temp);
            }
        }
}

void doubCompApp_poly_mul_classical( doubCompApp_poly_t res, const doubCompApp_poly_t x, const doubCompApp_poly_t y){
    
    slong lenres, i;
    
    if (x->length == 0 || y->length == 0){
        doubCompApp_poly_zero(res);
        return;
    }
    
    lenres = x->length + y->length -1;
    
    if (res==x || res==y) {/*aliasing*/
        doubCompApp_poly_t t;
        doubCompApp_poly_init2(t, lenres);
        /* initialize the coefficients of t to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( t->coeffs + i );
    
        _doubCompApp_poly_mullow_classical(t->coeffs, x->coeffs, x->length,
                                                      y->coeffs, y->length, lenres);
        doubCompApp_poly_swap(res, t);
        doubCompApp_poly_clear(t);
    }
    else {
        doubCompApp_poly_fit_length(res, lenres);
        /* initialize the coefficients of res to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( res->coeffs + i );
        _doubCompApp_poly_mullow_classical(res->coeffs, x->coeffs, x->length,
                                                        y->coeffs, y->length, lenres);
    }
    _doubCompApp_poly_set_length(res, lenres);
    _doubCompApp_poly_normalise(res);
    
    
}

/* assume lenx=leny */
#define KARATSUBA_CUTOFF 8
void _doubCompApp_poly_mullow_karatsuba(doubCompApp_ptr res,
    doubCompApp_srcptr x, slong lenx,
    doubCompApp_srcptr y, slong leny){
    
    if (lenx<=KARATSUBA_CUTOFF || leny<=KARATSUBA_CUTOFF){
        _doubCompApp_poly_mullow_classical(res, x, lenx, y, leny, lenx+leny-1);
        return;
    }
    
    slong lenres = lenx + leny -1;
    slong cut = lenx/2;
    slong maxlen = lenx-cut;
    
    _doubCompApp_poly_mullow_karatsuba(res, x, cut, y, cut);
    _doubCompApp_poly_mullow_karatsuba(res + (2*cut), x + cut, maxlen, y+cut, maxlen);
    
    doubCompApp_poly_t t1, t2, t3;
    doubCompApp_poly_init2(t1, maxlen);
    doubCompApp_poly_init2(t2, maxlen);
    doubCompApp_poly_init2(t3, 3*maxlen-1);
    slong i;
    for (i = 0; i<3*maxlen-1; i++)
            doubCompApp_zero( t3->coeffs + i );
    _doubCompApp_poly_add( t1->coeffs, x, cut, x + cut, maxlen, maxlen);
    _doubCompApp_poly_add( t2->coeffs, y, cut, y + cut, maxlen, maxlen);
    _doubCompApp_poly_mullow_karatsuba(t3->coeffs, t1->coeffs, maxlen, t2->coeffs, maxlen);
    _doubCompApp_poly_sub( t3->coeffs, t3->coeffs, 2*maxlen-1, res, 2*cut-1, 2*cut-1);
    _doubCompApp_poly_sub( t3->coeffs, t3->coeffs, 2*maxlen-1, res + (2*cut), 2*maxlen-1, 2*maxlen-1);
    _doubCompApp_poly_shift_left(t3->coeffs, t3->coeffs, 2*maxlen-1, cut);
//     _doubCompApp_poly_set_length(t3, 3*cut-1);
    
    _doubCompApp_poly_add( res, res, lenres, t3->coeffs, 2*maxlen-1 + cut, lenres);
    
    doubCompApp_poly_clear(t1);
    doubCompApp_poly_clear(t2);
    doubCompApp_poly_clear(t3);
}

void doubCompApp_poly_mul_karatsuba( doubCompApp_poly_t res, const doubCompApp_poly_t x, const doubCompApp_poly_t y){
    
    slong lenres, i;
    
    if (x->length == 0 || y->length == 0){
        doubCompApp_poly_zero(res);
        return;
    }
    
    lenres = x->length + y->length -1;
    
    if (res==x || res==y) {/*aliasing*/
        doubCompApp_poly_t t;
        doubCompApp_poly_init2(t, lenres);
        /* initialize the coefficients of t to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( t->coeffs + i );
    
        _doubCompApp_poly_mullow_karatsuba(t->coeffs, x->coeffs, x->length,
                                                      y->coeffs, y->length);
        doubCompApp_poly_swap(res, t);
        doubCompApp_poly_clear(t);
    }
    else {
        doubCompApp_poly_fit_length(res, lenres);
        /* initialize the coefficients of res to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( res->coeffs + i );
        _doubCompApp_poly_mullow_karatsuba(res->coeffs, x->coeffs, x->length,
                                                        y->coeffs, y->length);
    }
    _doubCompApp_poly_set_length(res, lenres);
    _doubCompApp_poly_normalise(res);
    
    
}

void _doubCompApp_poly_square_karatsuba(doubCompApp_ptr res,
    doubCompApp_srcptr x, slong lenx){
    
    if (lenx<=KARATSUBA_CUTOFF){
        _doubCompApp_poly_mullow_classical(res, x, lenx, x, lenx, 2*lenx-1);
        return;
    }
    
    slong lenres = 2*lenx -1;
    slong cut = lenx/2;
    
    _doubCompApp_poly_square_karatsuba(res, x, cut);
    _doubCompApp_poly_square_karatsuba(res + (2*cut), x + cut, lenx-cut);
    
    slong i;
    doubCompApp_poly_t t1,t2;
    doubCompApp_poly_init2(t1, lenx+cut-1);
    doubCompApp_poly_init2(t2, lenx-cut);
    for  (i=0; i<3*cut-1; i++)
        doubCompApp_zero(t1->coeffs + i);
    for  (i=0; i<cut; i++)
        doubCompApp_set(t2->coeffs + i, x + i);
    if (lenx-cut>cut)
        doubCompApp_zero(t2->coeffs + (lenx-cut-1));
    _doubCompApp_poly_mullow_karatsuba(t1->coeffs, t2->coeffs, lenx-cut, x+cut, lenx-cut);
    _doubCompApp_poly_mul_si(t1->coeffs, lenx-1, 2);
    _doubCompApp_poly_shift_left(t1->coeffs, t1->coeffs, lenx-1, cut);
    
    _doubCompApp_poly_add( res, res, lenres, t1->coeffs, lenx+cut -1, lenres);
    
    doubCompApp_poly_clear(t1);
    doubCompApp_poly_clear(t2);
}

void doubCompApp_poly_sqr_karatsuba( doubCompApp_poly_t res, const doubCompApp_poly_t x){
    
    slong lenres, i;
    
    if (x->length == 0){
        doubCompApp_poly_zero(res);
        return;
    }
    
    lenres = 2*(x->length) -1;
    
    if (res==x) {/*aliasing*/
        doubCompApp_poly_t t;
        doubCompApp_poly_init2(t, lenres);
        /* initialize the coefficients of t to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( t->coeffs + i );
    
        _doubCompApp_poly_square_karatsuba(t->coeffs, x->coeffs, x->length);
        doubCompApp_poly_swap(res, t);
        doubCompApp_poly_clear(t);
    }
    else {
        doubCompApp_poly_fit_length(res, lenres);
        /* initialize the coefficients of res to 0*/
        for (i = 0; i<lenres; i++)
            doubCompApp_zero( res->coeffs + i );
        _doubCompApp_poly_square_karatsuba(res->coeffs, x->coeffs, x->length);
    }
    _doubCompApp_poly_set_length(res, lenres);
    _doubCompApp_poly_normalise(res);
    
    
}

/* DEPRECATED */

// void _doubCompApp_poly_mullow_block(doubCompApp_ptr res,
//     doubCompApp_srcptr x, slong lenx,
//     doubCompApp_srcptr y, slong leny){
//     
// //     if (lenx<=KARATSUBA_CUTOFF || leny<=KARATSUBA_CUTOFF){
// //         _doubCompApp_poly_mullow_classical(res, x, lenx, y, leny);
// //         return;
// //     }
//     if (lenx<=1){
//         doubCompApp_mul(res, x, y);
//         return;
//     }
//     
//     /*find the cut*/
//     slong lenres = lenx + leny -1;
//     int cut = 0x1;
//     while ( (cut<<1) < lenx) 
//         cut = cut<<1;
//     
// //     printf("----call, lenx: %d, leny: %d, cut: %d----\n", lenx, leny, cut);
//     
//     _doubCompApp_poly_mullow_block(res, x, cut, y, cut);
//     _doubCompApp_poly_mullow_block(res + (2*cut), x + cut, lenx-cut, y+cut, leny-cut);
//     
// //     printf("res: \n");
// //     for (int i=0; i<3; i++){
// //         printf("deg %d, coeff: ",i);
// //         doubCompApp_print(res + i);
// //         printf("\n");
// //     }
//     
// //     printf("t1: \n");
// //     for (int i=0; i<cut; i++){
// //         printf("deg %d, coeff: ",i);
// //         doubCompApp_print(t1->coeffs + i);
// //         printf("\n");
// //     }
//     
// //     printf("t2: \n");
// //     for (int i=0; i<cut; i++){
// //         printf("deg %d, coeff: ",i);
// //         doubCompApp_print(t2->coeffs + i);
// //         printf("\n");
// //     }
//     
// //     printf("t3: \n");
// //     for (int i=0; i<2*cut-1; i++){
// //         printf("deg %d, coeff: ",i);
// //         doubCompApp_print(t3->coeffs + i);
// //         printf("\n");
// //     }
//     
//     doubCompApp_poly_t t1, t2, t3;
//     doubCompApp_poly_init2(t1, 2*cut-1);
//     doubCompApp_poly_init2(t2, 2*cut-1);
//     doubCompApp_poly_init2(t3, 3*cut-1);
//     slong i;
//     for (i=0; i<2*cut-1; i++) {
//         doubCompApp_zero(t1->coeffs + i);
//         doubCompApp_zero(t2->coeffs + i);
//     }
//     for  (i=0; i<3*cut-1; i++)
//         doubCompApp_zero(t3->coeffs + i);
//     for (i=0; i<cut; i++)
//         if (i< (leny-cut))
//             doubCompApp_set(t3->coeffs + i, y + (cut + i) );
//         else 
//             doubCompApp_zero(t3->coeffs + i);
//         
//     _doubCompApp_poly_mullow_block(t1->coeffs, x, cut, t3->coeffs, cut);
//     for (i=0; i<cut; i++)
//         if (i< (lenx-cut))
//             doubCompApp_set(t3->coeffs + i, x + (cut + i));
//         else 
//             doubCompApp_zero(t3->coeffs + i);
//     _doubCompApp_poly_mullow_block(t2->coeffs, y, cut, t3->coeffs, cut);
//     _doubCompApp_poly_add( t3->coeffs, t1->coeffs, 2*cut-1, t2->coeffs, 2*cut-1, 2*cut-1);
//     _doubCompApp_poly_shift_left(t3->coeffs, t3->coeffs, 2*cut-1, cut);
//     _doubCompApp_poly_set_length(t3, 3*cut-1);
//     
//     _doubCompApp_poly_add( res, res, lenres, t3->coeffs, 3*cut -1, 3*cut -1);
//     
//     doubCompApp_poly_clear(t1);
//     doubCompApp_poly_clear(t2);
//     doubCompApp_poly_clear(t3);
// }
// 
// void doubCompApp_poly_mul_block( doubCompApp_poly_t res, const doubCompApp_poly_t x, const doubCompApp_poly_t y){
//     
//     slong lenres, i;
//     
//     if (x->length == 0 || y->length == 0){
//         doubCompApp_poly_zero(res);
//         return;
//     }
//     
//     lenres = x->length + y->length -1;
//     
//     if (res==x || res==y) {/*aliasing*/
//         doubCompApp_poly_t t;
//         doubCompApp_poly_init2(t, lenres);
//         /* initialize the coefficients of t to 0*/
//         for (i = 0; i<lenres; i++)
//             doubCompApp_zero( t->coeffs + i );
//     
//         _doubCompApp_poly_mullow_block(t->coeffs, x->coeffs, x->length,
//                                                       y->coeffs, y->length);
//         doubCompApp_poly_swap(res, t);
//         doubCompApp_poly_clear(t);
//     }
//     else {
//         doubCompApp_poly_fit_length(res, lenres);
//         /* initialize the coefficients of res to 0*/
//         for (i = 0; i<lenres; i++)
//             doubCompApp_zero( res->coeffs + i );
//         _doubCompApp_poly_mullow_block(res->coeffs, x->coeffs, x->length,
//                                                         y->coeffs, y->length);
//     }
//     _doubCompApp_poly_set_length(res, lenres);
//     _doubCompApp_poly_normalise(res);
//     
//     
// }
