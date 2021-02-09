/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "realApp_poly.h"

/*order one center evaluation */
// void realApp_poly_evaluate_order_one( realApp_t y, const realApp_poly_t f, const realApp_poly_t fder, const realApp_t x, slong prec){
//     
//     realApp_t xmid, fvalder;
//     realApp_init(xmid);
//     realApp_init(fvalder);
//     
//     realApp_get_rad_realApp(xmid, x);
//     realApp_poly_evaluate(y, f, xmid, prec);
//     realApp_poly_evaluate(fvalder, fder, x, prec); 
//     realApp_sub(xmid, x, xmid, prec);
//     realApp_mul(fvalder, fvalder, xmid, prec);
//     realApp_add(y, y, fvalder, prec);
//     
//     realApp_clear(xmid);
//     realApp_clear(fvalder);
// }

int realApp_poly_check_relOne_accuracy( const realApp_poly_t poly, slong prec) {
    int result = 1;
    slong i = 0;
    while ( ( i<=realApp_poly_degree(poly) ) && result ) {
        result = result && realApp_check_relOne_accuracy( realApp_poly_getCoeff( poly, i) , prec );
        i++;
    }
    return result;
}

slong realApp_poly_get_relOne_accuracy_min( const realApp_poly_t poly){
    slong i = 0;
    slong res = realApp_get_relOne_accuracy( realApp_poly_getCoeff( poly, i) );
    while ( i<=realApp_poly_degree(poly) ) {
        res = CCLUSTER_MIN(res, realApp_get_relOne_accuracy( realApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

slong realApp_poly_get_relOne_accuracy_max( const realApp_poly_t poly){
    slong i = 0;
    slong res = realApp_get_relOne_accuracy( realApp_poly_getCoeff( poly, i) );
    while ( i<=realApp_poly_degree(poly) ) {
        res = CCLUSTER_MAX(res, realApp_get_relOne_accuracy( realApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

slong realApp_poly_get_absolute_accuracy_min( const realApp_poly_t poly){
    slong i = 0;
    slong res = realApp_get_absolute_accuracy( realApp_poly_getCoeff( poly, i) );
    while ( i<=realApp_poly_degree(poly) ) {
        res = CCLUSTER_MIN(res, realApp_get_absolute_accuracy( realApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

slong realApp_poly_get_absolute_accuracy_max( const realApp_poly_t poly){
    slong i = 0;
    slong res = realApp_get_absolute_accuracy( realApp_poly_getCoeff( poly, i) );
    while ( i<=realApp_poly_degree(poly) ) {
        res = CCLUSTER_MAX(res, realApp_get_absolute_accuracy( realApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

void realApp_poly_sum_abs_coeffs( realApp_t res, const realApp_poly_t f, slong prec ) {
    
    realApp_srcptr fptr = f->coeffs;
    const slong len = f->length;
    slong i;
    realApp_t temp;
    
    realApp_init(temp);
/*     realApp_zero(res);*/
    realApp_abs(res, fptr );
    for (i = 1; i < len; i++){
        realApp_abs(temp, fptr + i);
        realApp_add(res, res, temp, prec);
    }
    
    realApp_clear(temp);
}

/* computes the coeff of x^(index) of one Graeffe iteration of f, where f has length and degree len-1 */
/* requires: the 2*index first coeffs of f */
void realApp_poly_oneGraeffeIteration_coeff( realApp_ptr coeff, realApp_srcptr f, slong index, slong len, slong prec) {
    
    slong deltaindex, i, deg;
    realApp_t temp;
    
    realApp_init(temp);
    
    realApp_mul(coeff, f+index, f+index, prec ); /* sets coeff to coeff(f,index)*coeff(f,index) */
    if (index%2 == 1)
        realApp_neg(coeff, coeff);               /* sets coeff to (-1)^index * coeff(f,index)*coeff(f,index) */
    
    deg = len-1;
    deltaindex = ( index <= (deg - index) ? index : (deg - index) );
    
    for ( i = 1; i<=deltaindex; i++){
        realApp_mul(temp, f+(index-i), f+(index+i), prec ); /* sets temp to coeff(f,index-i)*coeff(f,index + i) */
        realApp_mul_si(temp, temp, 2, prec);                /* sets temp to 2*coeff(f,index-i)*coeff(f,index + i) */
        if ((index-i)%2 == 1)
            realApp_neg(temp, temp);                        /* sets temp to (-1)^(index-i) * 2*coeff(f,index-i)*coeff(f,index + i) */
        realApp_add( coeff, temp, coeff, prec);
    }
        
    realApp_clear(temp);
}

/* computes the coeffs up to power n of one Graeffe iteration of f, where f has length len and degree len-1 */
/* require: the 2n first coeffs of f */
/*          res has enough space to contain n+1 coeffs */ 
void realApp_poly_oneGraeffeIteration_first_n_coeff( realApp_poly_t res, const realApp_poly_t f, slong n, slong len, slong prec){
    
    slong i;
/*     realApp_poly_fit_length(res, n+1);*/
    realApp_srcptr fptr = f->coeffs;
    realApp_ptr resptr = res->coeffs;
    
    realApp_mul( resptr, fptr, fptr, prec ); /* sets coeff(res,0) to coeff(f,0)*coeff(f,0) */
    for (i=1; i<=n; i++)                     /* computes coeff(res,i) */ 
        realApp_poly_oneGraeffeIteration_coeff( resptr + i, fptr, i, len, prec);
    realApp_poly_set_length(res, n+1);
}

void realApp_poly_oneGraeffeIteration_in_place( realApp_poly_t f, slong prec ){
    
    realApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    realApp_poly_t fe, fo;
    realApp_poly_init2(fe, len2);
    realApp_poly_init2(fo, len2);
    realApp_ptr feptr = fe->coeffs;
    realApp_ptr foptr = fo->coeffs;
    
    for (i = 0; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) 
            realApp_set( feptr + quo, fptr+i);
        else
            realApp_set( foptr + quo, fptr+i);
    }
    realApp_poly_set_length(fe, len2);
    realApp_poly_set_length(fo, len2);
    
    realApp_poly_t fes, fos;
    realApp_poly_init2(fes, len1);
    realApp_poly_init2(fos, len1);
    realApp_poly_mullow( fes, fe, fe, len1, prec);
    realApp_poly_mullow( fos, fo, fo, len1, prec);
    realApp_poly_shift_left( fos, fos, 1 );
    realApp_poly_sub(f, fes, fos, prec);
    if ((len1%2)==0)
        realApp_poly_neg(f, f);
    
    realApp_poly_clear(fe);
    realApp_poly_clear(fo);
    realApp_poly_clear(fes);
    realApp_poly_clear(fos);
    
}
