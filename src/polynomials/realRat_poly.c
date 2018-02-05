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

#include "polynomials/realRat_poly.h"

/* here we assume that the coefficients of the pol are interger! */
void realRat_poly_separationBound (realRat_t sep, const realRat_poly_t pol){
    
    realRat_t lcoeff, coefftemp;
    realRat_init(lcoeff);
    realRat_init(coefftemp);
    slong degree = realRat_poly_degree(pol);
    
    /*get the largest coeff*/
    realRat_poly_get_coeff_realRat(lcoeff, pol, 0);
    realRat_abs(lcoeff, lcoeff);
    for (int i=1;i<=degree; i++){
        realRat_poly_get_coeff_realRat(coefftemp, pol, i);
        realRat_abs(coefftemp, coefftemp);
        if (realRat_cmp(coefftemp,lcoeff)>0)
            realRat_set(lcoeff, coefftemp);
    }
    
    realRat_add_si(lcoeff, lcoeff, 1);
    realRat_pow_si(lcoeff, lcoeff, degree);
    
    realRat_set_si(sep, degree, 1);
    realRat_pow_si(sep, sep, ((slong) ((degree+4)/2))+1);
    realRat_mul_si(sep, sep, 2);
    realRat_mul(sep, sep, lcoeff);
    realRat_inv(sep, sep);
    
    realRat_clear(lcoeff);
    realRat_clear(coefftemp);
}

/*bitsize for integer polynomials*/
/* here we assume that the coefficients of the pol are interger! */
slong realRat_poly_bitsize (const realRat_poly_t pol){
    realRat_t lcoeff, coefftemp;
    realRat_init(lcoeff);
    realRat_init(coefftemp);
    slong degree = realRat_poly_degree(pol);
    slong bitsize = 0;
    /*get the largest coeff*/
    realRat_poly_get_coeff_realRat(lcoeff, pol, 0);
    realRat_abs(lcoeff, lcoeff);
    for (int i=1;i<=degree; i++){
        realRat_poly_get_coeff_realRat(coefftemp, pol, i);
        realRat_abs(coefftemp, coefftemp);
        if (realRat_cmp(coefftemp,lcoeff)>0)
            realRat_set(lcoeff, coefftemp);
    }
    bitsize = fmpz_clog_ui( realRat_numref(lcoeff), 2);
    realRat_clear(lcoeff);
    realRat_clear(coefftemp);
    return bitsize;
}

void mignotte_polynomial(realRat_poly_t poly, slong deg, slong bitsize){
    
    realRat_t coeff, two;
    realRat_init(coeff);
    realRat_init(two);
    
    realRat_set_si(two, 2,1);
//     realRat_pow_si(coeff, two, ((slong) bitsize/2)-1 );
    realRat_pow_si(coeff, two, bitsize);
    
    realRat_poly_fit_length(poly,deg+1);
    realRat_poly_set_coeff_realRat(poly, 1, coeff);
    
    realRat_set_si(coeff, -1,1);
    realRat_poly_set_coeff_realRat(poly, 0, coeff);
    realRat_poly_pow(poly, poly, 2);
    realRat_poly_add(poly, poly, poly);
    fmpq_poly_neg(poly, poly);//
    
    realRat_set_si(coeff, 1,1);
    realRat_poly_set_coeff_realRat(poly, deg, coeff);
    
    realRat_clear(coeff);
    realRat_clear(two);
}