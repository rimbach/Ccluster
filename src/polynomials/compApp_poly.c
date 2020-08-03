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

#include "compApp_poly.h"

int compApp_poly_check_relOne_accuracy( const compApp_poly_t poly, slong prec) {
    int result = 1;
    slong i = 0;
    while ( ( i<=compApp_poly_degree(poly) ) && result ) {
        result = result && compApp_check_relOne_accuracy( compApp_poly_getCoeff( poly, i) , prec );
        i++;
    }
    return result;
}

slong compApp_poly_get_relOne_accuracy_min( const compApp_poly_t poly){
    slong i = 0;
    slong res = compApp_get_relOne_accuracy( compApp_poly_getCoeff( poly, i) );
    while ( i<=compApp_poly_degree(poly) ) {
        res = CCLUSTER_MIN(res, compApp_get_relOne_accuracy( compApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

slong compApp_poly_get_relOne_accuracy_max( const compApp_poly_t poly){
    slong i = 0;
    slong res = compApp_get_relOne_accuracy( compApp_poly_getCoeff( poly, i) );
    while ( i<=compApp_poly_degree(poly) ) {
        res = CCLUSTER_MAX(res, compApp_get_relOne_accuracy( compApp_poly_getCoeff( poly, i) ) );
        i++;
    }
    return res;
}

void compApp_poly_sum_abs_coeffs( realApp_t res, const compApp_poly_t f, slong prec ) {
    
    compApp_srcptr fptr = f->coeffs;
    const slong len = f->length;
    slong i;
    realApp_t temp;
    
    realApp_init(temp);
/*     realApp_zero(res);*/
    compApp_abs(res, fptr, prec );
    for (i = 1; i < len; i++){
        compApp_abs(temp, fptr + i, prec);
        realApp_add(res, res, temp, prec);
    }
    
    realApp_clear(temp);
}

/* computes the coeff of x^(index) of one Graeffe iteration of f, where f has length and degree len-1 */
/* requires: the 2*index first coeffs of f */
void compApp_poly_oneGraeffeIteration_coeff( compApp_ptr coeff, compApp_srcptr f, slong index, slong len, slong prec) {
    
    slong deltaindex, i, deg;
    compApp_t temp;
    
    compApp_init(temp);
    
    compApp_mul(coeff, f+index, f+index, prec ); /* sets coeff to coeff(f,index)*coeff(f,index) */
    if (index%2 == 1)
        compApp_neg(coeff, coeff);               /* sets coeff to (-1)^index * coeff(f,index)*coeff(f,index) */
    
    deg = len-1;
    deltaindex = ( index <= (deg - index) ? index : (deg - index) );
    
    for ( i = 1; i<=deltaindex; i++){
        compApp_mul(temp, f+(index-i), f+(index+i), prec ); /* sets temp to coeff(f,index-i)*coeff(f,index + i) */
        compApp_mul_si(temp, temp, 2, prec);                /* sets temp to 2*coeff(f,index-i)*coeff(f,index + i) */
        if ((index-i)%2 == 1)
            compApp_neg(temp, temp);                        /* sets temp to (-1)^(index-i) * 2*coeff(f,index-i)*coeff(f,index + i) */
        compApp_add( coeff, temp, coeff, prec);
    }
        
    compApp_clear(temp);
}

/* computes the coeffs up to power n of one Graeffe iteration of f, where f has length len and degree len-1 */
/* require: the 2n first coeffs of f */
/*          res has enough space to contain n+1 coeffs */ 
void compApp_poly_oneGraeffeIteration_first_n_coeff( compApp_poly_t res, const compApp_poly_t f, slong n, slong len, slong prec){
    
    slong i;
/*     compApp_poly_fit_length(res, n+1);*/
    compApp_srcptr fptr = f->coeffs;
    compApp_ptr resptr = res->coeffs;
    
    compApp_mul( resptr, fptr, fptr, prec ); /* sets coeff(res,0) to coeff(f,0)*coeff(f,0) */
    for (i=1; i<=n; i++)                     /* computes coeff(res,i) */ 
        compApp_poly_oneGraeffeIteration_coeff( resptr + i, fptr, i, len, prec);
    compApp_poly_set_length(res, n+1);
}

void compApp_poly_oneGraeffeIteration_in_place( compApp_poly_t f, slong prec ){
    
    compApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    
    compApp_poly_t fe, fo;
    compApp_poly_init2(fe, len2);
    compApp_poly_init2(fo, len2);
    compApp_ptr feptr = fe->coeffs;
    compApp_ptr foptr = fo->coeffs;
    
    for (i = 0; i < len1; i++){
        rem = i%2;
        quo = i>>1;
        if (rem == 0) 
            compApp_set( feptr + quo, fptr+i);
        else
            compApp_set( foptr + quo, fptr+i);
    }
    compApp_poly_set_length(fe, len2);
    compApp_poly_set_length(fo, len2);
    
    compApp_poly_t fes, fos;
    compApp_poly_init2(fes, len1);
    compApp_poly_init2(fos, len1);
    compApp_poly_mullow( fes, fe, fe, len1, prec);
    compApp_poly_mullow( fos, fo, fo, len1, prec);
    compApp_poly_shift_left( fos, fos, 1 );
    compApp_poly_sub(f, fes, fos, prec);
    if ((len1%2)==0)
        compApp_poly_neg(f, f);
    
    compApp_poly_clear(fe);
    compApp_poly_clear(fo);
    compApp_poly_clear(fes);
    compApp_poly_clear(fos);
    
}
