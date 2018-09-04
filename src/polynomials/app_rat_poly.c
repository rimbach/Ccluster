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

#include "app_rat_poly.h"

void compApp_poly_scale_realRat_in_place( compApp_ptr fptr, const realRat_t r, slong len, slong prec ){
    
    slong i;
    realRat_t t;
    
    realRat_init(t);
    realRat_set(t, r);
/*     realRat_canonicalise(t); */
    
    for (i = 1; i < len; i++){
        compApp_mul_realRat_in_place( fptr + i, t, prec);
        if (i + 1 < len)
            realRat_mul(t, t, r);
    }
    realRat_clear(t);
}

/*
void compApp_poly_taylorShift_in_place_new( compApp_poly_t f, const realRat_t creal, const realRat_t cimag, const realRat_t radius, slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, creal, prec);
    compApp_setimag_realRat(c, cimag, prec);
    
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _compApp_poly_taylor_shift_convolution(fptr, c, len, prec);
    compApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
    
}
*/

void compApp_poly_taylor_shift_convolution_without_pre(compApp_poly_t dest, const compApp_poly_t p, 
                                                       realApp_t f, compApp_ptr t, 
                                                       const realRat_t creal, const realRat_t cimag, const realRat_t radius,
                                                       slong prec){
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, creal, prec);
    compApp_setimag_realRat(c, cimag, prec);
    
    compApp_poly_set(dest, p);
    
    _compApp_poly_taylor_shift_convolution(dest->coeffs, f, t, c, dest->length, prec);

    compApp_poly_scale_realRat_in_place( dest->coeffs, radius, dest->length, prec );
    
    compApp_clear(c);
}

void compApp_poly_taylorShift_new( compApp_poly_t res, 
                               const compApp_poly_t f, 
                               const realRat_t creal, const realRat_t cimag, const realRat_t radius, 
                               slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, creal, prec);
    compApp_setimag_realRat(c, cimag, prec);
    
    compApp_poly_taylor_shift_convolution(res, f, c, prec);
    compApp_poly_scale_realRat_in_place( res->coeffs, radius, res->length, prec );
    
    compApp_clear(c);
}

// void compApp_poly_taylorShift_in_place( compApp_poly_t f, const realRat_t creal, const realRat_t cimag, const realRat_t radius, slong prec ) {
//     
//     compApp_t c;
//     compApp_init(c);
//     compApp_setreal_realRat(c, creal, prec);
//     compApp_setimag_realRat(c, cimag, prec);
//     
//     compApp_ptr fptr = f->coeffs;
//     const slong len  = f->length;
//     _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
// /*    _acb_poly_taylor_shift_horner(fptr, c, len, prec);
//     _acb_poly_taylor_shift_divconquer(fptr, c, len, prec);
//     _compApp_poly_taylor_shift_convolution(fptr, c, len, prec);*/
//     compApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
//     
//     compApp_clear(c);
// }
void compApp_poly_taylorShift_in_place( compApp_poly_t f, const compRat_t center, const realRat_t radius, slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, compRat_realref(center), prec);
    compApp_setimag_realRat(c, compRat_imagref(center), prec);
    
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
/*    _acb_poly_taylor_shift_horner(fptr, c, len, prec);
    _acb_poly_taylor_shift_divconquer(fptr, c, len, prec);
    _compApp_poly_taylor_shift_convolution(fptr, c, len, prec);*/
    compApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
    
    compApp_clear(c);
}

void compApp_poly_taylorShift( compApp_poly_t res, 
                               const compApp_poly_t f, 
                               const compRat_t center, const realRat_t radius, 
                               slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, compRat_realref(center), prec);
    compApp_setimag_realRat(c, compRat_imagref(center), prec);
    
    acb_poly_taylor_shift_convolution(res, f, c, prec);
/*     compApp_poly_taylor_shift_convolution(res, f, c, prec);*/
    compApp_poly_scale_realRat_in_place( res->coeffs, radius, res->length, prec );
    
    compApp_clear(c);
}
