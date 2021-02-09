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

void compApp_poly_root_bound_fujiwara(realRat_t bound, compApp_poly_t poly){
    mag_t b;
    mag_init(b);
    acb_poly_root_bound_fujiwara(b, poly);
    mag_get_fmpq(bound, b);
    mag_clear(b);
}

void realApp_poly_scale_realRat_in_place( realApp_ptr fptr, const realRat_t r, slong len, slong prec){
    
    slong i;
    realApp_t temp, factor;
    
    realApp_init(temp);
    realApp_init(factor);
    realApp_set_realRat(temp, r, prec);
    realApp_set(factor, temp);
    
    for (i = 1; i < len; i++){
        realApp_mul( fptr + i, fptr + i, factor, prec);
        if (i + 1 < len)
            realApp_mul(factor, factor, temp, prec);
    }
    realApp_clear(temp);
    realApp_clear(factor);
    
}

void compApp_poly_scale_realRat_in_place( compApp_ptr fptr, const realRat_t r, slong len, slong prec ){
    
    slong i;
    realApp_t temp, factor;
    
    realApp_init(temp);
    realApp_init(factor);
    realApp_set_realRat(temp, r, prec);
    realApp_set(factor, temp);
    
    for (i = 1; i < len; i++){
        compApp_mul_realApp( fptr + i, fptr + i, factor, prec);
        if (i + 1 < len)
            realApp_mul(factor, factor, temp, prec);
    }
    realApp_clear(temp);
    realApp_clear(factor);
}

/* evaluate at order 1 on an interval geven by center + width */
void realApp_poly_evaluate_order_one( realApp_t y, const realApp_poly_t f, const realApp_poly_t fder, const realRat_t c, const realRat_t w, slong prec){
    
    realApp_t center, error, fderval;
    
    realApp_init(center);
    realApp_init(error);
    realApp_init(fderval);
    
    realApp_set_realRat(center, c, prec);
    realApp_set_realRat(error, w, prec);
    realApp_div_ui(     error, error, 2, prec );
    /*evaluate f at center*/
    realApp_poly_evaluate(y, f, center, prec);
    /*evaluate fder at interval*/
    realApp_add_error(center, error);
    realApp_poly_evaluate(fderval, fder, center, prec);
    /*compute interval - center */
    realApp_zero(center);
    realApp_add_error(center, error);
    realApp_mul(fderval, fderval, center, prec);
    realApp_add(y, y, fderval, prec);
    
    
    realApp_clear(center);
    realApp_clear(error);
    realApp_init(fderval);
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

void realApp_poly_taylorShift_in_place( realApp_poly_t f, const realRat_t center, const realRat_t radius, slong prec ) {
    
    realApp_t c;
    realApp_init(c);
    realApp_set_realRat( c, center, prec );
    
    realApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _arb_poly_taylor_shift_convolution(fptr, c, len, prec);
//     arb_poly_taylor_shift_divconquer(f, f, c, prec);

    realApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
    
    realApp_clear(c);
}

void compApp_poly_taylorShift_in_place( compApp_poly_t f, const compRat_t center, const realRat_t radius, slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, compRat_realref(center), prec);
    compApp_setimag_realRat(c, compRat_imagref(center), prec);
    
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
//     acb_poly_taylor_shift_convolution(f, f, c, prec);
/*    _acb_poly_taylor_shift_horner(fptr, c, len, prec);
    _acb_poly_taylor_shift_divconquer(fptr, c, len, prec);
    _compApp_poly_taylor_shift_convolution(fptr, c, len, prec);*/

    compApp_poly_scale_realRat_in_place( fptr, radius, len, prec );
    
    compApp_clear(c);
}


void compApp_poly_taylorShift_in_place_noscale( compApp_poly_t f, 
                                        const compRat_t center, 
                                        slong prec ){
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, compRat_realref(center), prec);
    compApp_setimag_realRat(c, compRat_imagref(center), prec);
    
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
    
    compApp_clear(c);
}

void realApp_poly_taylorShift_in_place_noscale( realApp_poly_t f, 
                                        const realRat_t center, 
                                        slong prec ){
    realApp_t c;
    realApp_init(c);
    realApp_set_realRat(c, center, prec);
    
    realApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _arb_poly_taylor_shift_convolution(fptr, c, len, prec);
    
    realApp_clear(c);
}

void compApp_poly_taylorShift_interval_in_place( compApp_poly_t f, const compRat_t center, const realRat_t radius, slong prec ) {
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, compRat_realref(center), prec);
    compApp_setimag_realRat(c, compRat_imagref(center), prec);
    
    realApp_t error;
    realApp_init(error);
    realApp_set_realRat(error, radius, prec);
    realApp_add_error( compApp_realref(c), error );
    realApp_add_error( compApp_imagref(c), error );
    
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
    
    compApp_clear(c);
    realApp_clear(error);
}

void realApp_poly_taylorShift_interval_in_place( realApp_poly_t f, const realRat_t center, const realRat_t radius, slong prec ){
    realApp_t c;
    realApp_init(c);
    realApp_set_realRat(c, center, prec);
    
    realApp_t error;
    realApp_init(error);
    realApp_set_realRat(error, radius, prec);
    realApp_add_error( c, error );
    
    realApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _arb_poly_taylor_shift_convolution(fptr, c, len, prec);
    
    realApp_clear(c);
    realApp_clear(error);
}

void realApp_poly_taylorShift_in_place_slong( realApp_poly_t f, 
                                              slong centerRe, 
                                              slong prec ){
    realApp_t c;
    realApp_init(c);
    realApp_set_si(c, centerRe);
    realApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _arb_poly_taylor_shift_convolution(fptr, c, len, prec);
    realApp_clear(c);
}

void compApp_poly_taylorShift_in_place_slong( compApp_poly_t f, 
                                              slong centerRe, 
                                              slong centerIm,
                                              slong prec ){
    compApp_t c;
    compApp_init(c);
    compApp_set_sisi(c, centerRe, centerIm);
    compApp_ptr fptr = f->coeffs;
    const slong len  = f->length;
    _acb_poly_taylor_shift_convolution(fptr, c, len, prec);
    compApp_clear(c);

}

void realApp_poly_taylorShift( realApp_poly_t res, 
                               const realApp_poly_t f, 
                               const realRat_t center, const realRat_t radius, 
                               slong prec ) {
    
    realApp_t c;
    realApp_init(c);
    realApp_set_realRat( c, center, prec );
    
    arb_poly_taylor_shift_convolution(res, f, c, prec);
    realApp_poly_scale_realRat_in_place( res->coeffs, radius, res->length, prec );
    
    realApp_clear(c);
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
    compApp_poly_scale_realRat_in_place( res->coeffs, radius, res->length, prec );
    
    compApp_clear(c);
}

///DEPRECATED
// old versions of scale where the rational number is used (instead of a flot. point approx)
// void realApp_poly_scale_realRat_in_place( realApp_ptr fptr, const realRat_t r, slong len, slong prec){
//     
//     slong i;
//     realRat_t t;
//     
//     realRat_init(t);
//     realRat_set(t, r);
// /*     realRat_canonicalise(t); */
//     
//     for (i = 1; i < len; i++){
//         realApp_mul_realRat_in_place( fptr + i, t, prec);
//         if (i + 1 < len)
//             realRat_mul(t, t, r);
//     }
//     realRat_clear(t);
//     
// }
// 
// void compApp_poly_scale_realRat_in_place( compApp_ptr fptr, const realRat_t r, slong len, slong prec ){
//     
//     slong i;
//     realRat_t t;
//     
//     realRat_init(t);
//     realRat_set(t, r);
// /*     realRat_canonicalise(t); */
//     
//     for (i = 1; i < len; i++){
//         compApp_mul_realRat_in_place( fptr + i, t, prec);
//         if (i + 1 < len)
//             realRat_mul(t, t, r);
//     }
//     realRat_clear(t);
// }
