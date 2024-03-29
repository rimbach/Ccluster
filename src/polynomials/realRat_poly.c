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

/* here we assume that the coefficients of the pol are integer! */
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
/* here we assume that the coefficients of the pol are integer! */
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

/* multiple of Bernoulli polynomial with integer coefficients*/
void bernoulliInt_polynomial(realRat_poly_t poly, slong deg){
    bernoulli_polynomial( poly, deg);
    realRat_poly_canonicalise(poly);
}

void mignotte_polynomial(realRat_poly_t poly, slong deg, slong bitsize){
    
    realRat_t coeff, two;
    realRat_init(coeff);
    realRat_init(two);
    
    realRat_set_si(two, 2,1);
//     realRat_pow_si(coeff, two, ((slong) bitsize/2)-1 );
    realRat_pow_si(coeff, two, bitsize); /*coeff = 2^(bitsize) */
    
    realRat_poly_fit_length(poly,deg+1);
    realRat_poly_zero(poly);
    realRat_poly_set_coeff_realRat(poly, 1, coeff); /* poly = coeff*x */
    
    realRat_set_si(coeff, -1,1);
    realRat_poly_set_coeff_realRat(poly, 0, coeff); /* poly = coeff*x -1 */
    realRat_poly_pow(poly, poly, 2);                /* poly = (coeff*x -1)^2 */
    realRat_poly_add(poly, poly, poly);             /* poly = 2*(coeff*x -1)^2 */
    fmpq_poly_neg(poly, poly);//                    /* poly = -2*(coeff*x -1)^2 */
    
    realRat_set_si(coeff, 1,1);
    realRat_poly_set_coeff_realRat(poly, deg, coeff); /* poly = x^deg -2*(coeff*x -1)^2 */
    
    realRat_clear(coeff);
    realRat_clear(two);
}

/* assume bitsize has 8 as a factor */
/* assume deg has 4 as a factor */
void nested_mignotte_polynomial(realRat_poly_t poly, slong deg, slong bitsize){
    
    realRat_t coeff;
    realRat_poly_t templ, tempt, temp;
    realRat_init(coeff);
    realRat_poly_init(templ);
    realRat_poly_fit_length(templ,(deg/4)+2);
    realRat_poly_init(tempt);
    realRat_poly_fit_length(tempt,3);
    realRat_poly_init(temp);
    
    realRat_poly_one(poly);
    realRat_poly_one(templ);
    realRat_poly_zero(tempt);
    realRat_poly_shift_left(templ, templ, (deg/4)+1);
    
    realRat_set_si(coeff, 2,1);
    realRat_pow_si(coeff, coeff, bitsize/8); /*coeff = 2^(bitsize/8) */
    realRat_sub_si(coeff, coeff, 1); /*coeff = 2^(bitsize/8)-1 */
    realRat_poly_set_coeff_realRat(temp, 2, coeff); /* tempt = coeff*z^2 */
    realRat_poly_set_coeff_si_ui(temp, 0, -1, 1); /* tempt = coeff*z^2 - 1*/
    
    for (int i=1; i<=4; i++){
        realRat_poly_pow(tempt, temp, 2*i);  /* tempt = (coeff*z^2 - 1)^2*/
        realRat_poly_sub(tempt, templ, tempt);
        realRat_poly_mul(poly, poly, tempt);
    }
    
    realRat_clear(coeff);
    realRat_poly_clear(temp);
    realRat_poly_clear(tempt);
    realRat_poly_clear(templ);
}

void mignotte_generalized(realRat_poly_t poly, slong deg, ulong pow, slong bitsize){
    
    realRat_t coeff, two;
    realRat_init(coeff);
    realRat_init(two);
    realRat_poly_t p1, p2;
    realRat_poly_init(p1);
    realRat_poly_init(p2);
    realRat_poly_fit_length(poly,deg+1);
    realRat_poly_zero(poly);
    
    realRat_set_si(two, 2,1);
    realRat_pow_si(coeff, two, bitsize); /*coeff = 2^(bitsize) */
    realRat_poly_set_coeff_realRat(p1, 1, coeff); /* p1 = coeff*x */
    realRat_poly_set_coeff_realRat(p2, 1, coeff); /* p2 = coeff*x */
    realRat_set_si(coeff, -1,1);
    realRat_poly_set_coeff_realRat(p1, 0, coeff); /* p1 = coeff*x -1 */
    realRat_set_si(coeff, 1,1);
    realRat_poly_set_coeff_realRat(p2, 0, coeff); /* p2 = coeff*x +1 */
    realRat_poly_pow(p1, p1, pow);                /* p1 = (coeff*x -1)^pow */
    realRat_poly_pow(p2, p2, pow);                /* p2 = (coeff*x +1)^pow */
    realRat_poly_mul(poly, p1, p2);               /* poly = (coeff*x -1)^pow*(coeff*x +1)^pow*/
    realRat_poly_add(poly, poly, poly);           /* poly = 2*(coeff*x -1)^pow*(coeff*x +1)^pow*/
    fmpq_poly_neg(poly, poly);                    /* poly = -2*(coeff*x -1)^pow*(coeff*x +1)^pow*/
    
    realRat_set_si(coeff, 1,1);
    realRat_poly_set_coeff_realRat(poly, deg, coeff); /* poly = x^deg -2*(coeff*x -1)^pow*(coeff*x +1)^pow*/
    
    realRat_clear(coeff);
    realRat_clear(two);
    realRat_poly_clear(p1);
    realRat_poly_clear(p2);
}

void wilkinson_polynomial(realRat_poly_t poly, slong degree){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(poly);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(poly, poly, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}
