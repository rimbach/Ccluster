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
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include <time.h>

// void realApp_poly_scale_realRat_in_place( arb_ptr fptr, const realRat_t r, slong len, slong prec){
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

void realApp_poly_taylorShift( arb_poly_t res, 
                               const arb_poly_t f, 
                               const realRat_t center, const realRat_t radius, 
                               slong prec ) {
    
    realApp_t c;
    realApp_init(c);
    realApp_set_realRat( c, center, prec );
    
    arb_poly_taylor_shift_convolution(res, f, c, prec);
/*     compApp_poly_taylor_shift_convolution(res, f, c, prec);*/
    realApp_poly_scale_realRat_in_place( res->coeffs, radius, res->length, prec );
    
    realApp_clear(c);
}

// void realApp_poly_set_realRat_poly(arb_poly_t poly, const realRat_poly_t re, slong prec){
//     arb_poly_set_fmpq_poly(poly, re, prec);
// }

int main() {
    
    
    
    slong prec = 53;
    compApp_poly_t p, pres;
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    bernoulli_polynomial(pbern, 7);
    compApp_poly_init(p);
    compApp_poly_init(pres);
    
    realApp_poly_t rpol, rpolres;
    realApp_poly_init(rpol);
    realApp_poly_init(rpolres);
    
    
    realApp_poly_set_realRat_poly(rpol, pbern, prec);
    
    compApp_poly_set_realRat_poly( p, pbern, prec);
    
//     printf("p: \n"); compApp_poly_printd(p, prec); printf("\n\n");
//     printf("rpol: \n"); arb_poly_printd(rpol, prec); printf("\n\n");
    
    compRat_t center;
    realRat_t radius;
    compRat_init(center);
    realRat_init(radius);
    compRat_set_sisi(center, 1,1,0,1);
    realRat_set_si(radius, 1,16);
    
    compApp_poly_taylorShift(pres, p, center, radius, prec);
//     printf("p shifted in 1 with radius 1/2 (complex version): \n"); compApp_poly_printd(pres, 50); printf("\n\n");
    
    realApp_poly_taylorShift( rpolres, rpol, compRat_real_ptr(center), radius, prec);
//     printf("p shifted in 1 with radius 1/2 (real version): \n"); arb_poly_printd(rpolres, 50); printf("\n\n");
    
//     compApp_poly_taylorShift_new(pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec);
//     printf("p shifted in 1+i with radius 1 (new version): \n"); compApp_poly_printd(pres, 50); printf("\n\n");
    
//     realApp_t f;
//     compApp_ptr t;
//     compApp_poly_t ppre;
//     realApp_init(f);
//     t = _acb_vec_init(p->length);
//     compApp_poly_init(ppre);
//     compApp_poly_taylor_shift_conv_pre(ppre, p, f, t, prec);
//     compApp_poly_taylor_shift_convolution_without_pre(pres, ppre, f, t, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec);
// //     printf("p shifted in 1+i with radius 1 (new version): \n"); compApp_poly_printd(pres, 10); printf("\n\n");
//     _acb_vec_clear(t, p->length);
    
    int degs[7] = {8, 16, 32, 64, 128, 256, 512};
    int nbtests = 100;
    clock_t ti;
    
    for (int i = 0; i<7; i++) {
        bernoulli_polynomial(pbern, degs[i]);
        compApp_poly_set_realRat_poly( p, pbern, prec);
        realApp_poly_set_realRat_poly(rpol, pbern, prec);
//         t = _acb_vec_init(p->length);
//         compApp_poly_taylor_shift_conv_pre(ppre, p, f, t, prec);
        
        ti = clock();
        for (int j = 0; j<nbtests; j++)
            compApp_poly_taylorShift( pres, p, center, radius, prec);
        ti = clock() - ti;
        printf ("time for %d taylor shifts, degree %d: %f seconds.\n",nbtests, degs[i], ((float)ti)/CLOCKS_PER_SEC);
        
        ti = clock();
        for (int j = 0; j<nbtests; j++)
            realApp_poly_taylorShift( rpolres, rpol, compRat_real_ptr(center), radius, prec);
        ti = clock() - ti;
        printf ("time for %d taylor shifts, degree %d: %f seconds.\n",nbtests, degs[i], ((float)ti)/CLOCKS_PER_SEC);
        
//         ti = clock();
//         for (int j = 0; j<nbtests; j++)
//             compApp_poly_taylor_shift_convolution_without_pre(pres, ppre, f, t, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec);
//         ti = clock() - ti;
//         printf ("time for %d taylor shifts, degree %d: %f seconds.\n",nbtests, degs[i], ((float)ti)/CLOCKS_PER_SEC);
        
        printf("\n");
//         _acb_vec_clear(t, p->length);
    }
    
//     compApp_clear(c);
    compRat_clear(center);
    realRat_clear(radius);
    compApp_poly_clear(p);
//     compApp_poly_clear(ppre);
//     realApp_clear(f);
//     _acb_vec_clear(t, p->length);
    realRat_poly_clear(pbern);
    
    realApp_poly_clear(rpol);
    realApp_poly_clear(rpolres);
    return 0;
}
