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

#include <stdio.h>
#include<time.h>
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"

void oneGraeffeIteration_in_place( compApp_poly_t f, slong prec ){
    
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

clock_t oneGraeffeIteration_in_place2( compApp_poly_t f, compApp_poly_t fo, compApp_poly_t fe, compApp_poly_t fos, compApp_poly_t fes, slong prec){
    
    compApp_ptr fptr = f->coeffs;
    const slong len1 = f->length;
    const slong len2 = (len1/2)+1;
    slong i, rem, quo;
    clock_t el;
//     clock_t res;
    
//     compApp_poly_t fe, fo;
//     compApp_poly_init2(fe, len2);
//     compApp_poly_init2(fo, len2);
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
    
//     compApp_poly_t fes, fos;
//     compApp_poly_init2(fes, len1);
//     compApp_poly_init2(fos, len1);
    el = clock();
    compApp_poly_mullow( fes, fe, fe, len1, prec);
    compApp_poly_mullow( fos, fo, fo, len1, prec);
//     acb_poly_pow_ui(fes, fe, 2, prec);
//     acb_poly_pow_ui(fos, fo, 2, prec);
    el = clock() - el;
    compApp_poly_shift_left( fos, fos, 1 );
    compApp_poly_sub(f, fes, fos, prec);
    if ((len1%2)==0)
        compApp_poly_neg(f, f);
    
//     compApp_poly_clear(fe);
//     compApp_poly_clear(fo);
//     compApp_poly_clear(fes);
//     compApp_poly_clear(fos);
    return el;
}

int main() {
    
    slong degree = 362;
    slong prec = 106;
    int N = 9;
    int nbItts = 1000;
    clock_t cumul = clock() - clock();
    clock_t temp;
    clock_t cumul2 = clock() - clock();
    
    realRat_poly_t pbern;
    compRat_poly_t p_global;
    compApp_poly_t polyApprox, polyApprox2;
    
    realRat_poly_init2(pbern, degree+1);
    compRat_poly_init2(p_global, degree+1);
    compApp_poly_init2(polyApprox, degree+1);
    
    bernoulli_polynomial(pbern, degree);
    compRat_poly_set_realRat_poly(p_global,pbern);
    compApp_poly_set_compRat_poly(polyApprox, p_global, prec);
    
    compApp_poly_t fe, fo, fes, fos;
    compApp_poly_init2(fe, ((degree+1)/2)+1);
    compApp_poly_init2(fo, ((degree+1)/2)+1);
    compApp_poly_init2(fes, degree+1);
    compApp_poly_init2(fos, degree+1);
    
    realRat_t creal, cimag, radius;
    realRat_init(creal);
    realRat_init(cimag);
    realRat_init(radius);
    
    for(int i=0; i<nbItts; i++){
        compApp_poly_init2(polyApprox2, degree+1);
        compApp_poly_set(polyApprox2, polyApprox);
        
        realRat_set_si(creal, i, nbItts);
        realRat_set_si(cimag, i, nbItts);
        realRat_set_si(radius, i, nbItts);
        
//         compApp_poly_taylorShift_in_place( polyApprox2, creal, cimag, radius, prec );
        
        for (int j=0; j<N; j++){
//             oneGraeffeIteration_in_place( polyApprox2, prec );
            temp = clock();
            cumul += oneGraeffeIteration_in_place2( polyApprox2, fo, fe, fos, fes, prec );
            temp = clock()-temp;
            cumul2 +=temp;
        }
        
        compApp_poly_clear(polyApprox2);
    }
    
    printf("number of clicks for mullow: %d, time: %f\n", cumul, ((double)cumul)/CLOCKS_PER_SEC);
    printf("number of clicks for graeffe: %d, time: %f\n", cumul2, ((double)cumul2)/CLOCKS_PER_SEC);
    printf("number of clicks for graeffe without mullow: %d, time: %f\n", cumul2-cumul, ((double)(cumul2-cumul))/CLOCKS_PER_SEC);
    
    realRat_poly_clear(pbern);
    compRat_poly_clear(p_global);
    compApp_poly_clear(polyApprox);
    realRat_clear(creal);
    realRat_clear(cimag);
    realRat_clear(radius);
    compApp_poly_clear(fe);
    compApp_poly_clear(fo);
    compApp_poly_clear(fes);
    compApp_poly_clear(fos);
    
    return 0;
}