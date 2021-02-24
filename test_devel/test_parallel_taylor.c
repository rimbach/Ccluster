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
// #include <time.h>
#include <sys/time.h>

int main() {
    
    
    
    slong prec = 53;
    compApp_poly_t p, pres;
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    bernoulli_polynomial(pbern, 12);
    compApp_poly_init(p);
    compApp_poly_init(pres);
    
    
    compApp_poly_set_realRat_poly( p, pbern, prec);
    
    
    compRat_t center;
    realRat_t radius;
    compRat_init(center);
    realRat_init(radius);
    compRat_set_sisi(center, 1,1,1,1);
    realRat_set_si(radius, 1,2);
    
    compApp_poly_taylorShift(pres, p, center, radius, prec);
//     printf("p shifted in 1+i with radius 1: \n"); compApp_poly_printd(pres, 10); printf("\n\n");

//     compApp_poly_parallel_taylor(pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 1);
//     printf("p shifted (1-parallel) in 1+i with radius 1: \n"); compApp_poly_printd(pres, 10); printf("\n\n");
    
//     compApp_poly_parallel_taylor(pres,  p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 2);
//     printf("p shifted (2-parallel) in 1+i with radius 1: \n"); compApp_poly_printd(pres, 10); printf("\n\n");
    
//     compApp_poly_set_realRat_poly( p, pbern, prec);
    compApp_poly_parallel_taylor(pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 4);
//     printf("p shifted (4-parallel) in 1+i with radius 1: \n"); compApp_poly_printd(pres, 10); printf("\n\n");
    
    int degs[8] = {8, 16, 32, 64, 128, 256, 512, 1024};
    int nbtests[8] = {100, 100, 100, 100, 100, 100, 100, 100 };
    struct timeval tbegin, tend;
    
    for (int i = 0; i<8; i++) {
        bernoulli_polynomial(pbern, degs[i]);
        compApp_poly_set_realRat_poly( p, pbern, prec);
        
        gettimeofday(&tbegin, NULL);
        for (int j = 0; j<nbtests[i]; j++)
            compApp_poly_taylorShift( pres, p, center, radius, prec);
        gettimeofday(&tend, NULL);
        
        printf ("time for %d taylor shifts, degree %4d: %10ld useconds.\n",nbtests[i], degs[i], ((tend.tv_sec * 1000000 + tend.tv_usec) - (tbegin.tv_sec * 1000000 + tbegin.tv_usec)));
        
        gettimeofday(&tbegin, NULL);
        for (int j = 0; j<nbtests[i]; j++)
            compApp_poly_parallel_taylor( pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 1);
        gettimeofday(&tend, NULL);
        printf ("time for %d taylor shifts, degree %4d: %10ld useconds.\n",nbtests[i], degs[i], ((tend.tv_sec * 1000000 + tend.tv_usec) - (tbegin.tv_sec * 1000000 + tbegin.tv_usec)));
        
        gettimeofday(&tbegin, NULL);
        for (int j = 0; j<nbtests[i]; j++)
            compApp_poly_parallel_taylor( pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 2);
        gettimeofday(&tend, NULL);
        printf ("time for %d taylor shifts, degree %4d: %10ld useconds.\n",nbtests[i], degs[i], ((tend.tv_sec * 1000000 + tend.tv_usec) - (tbegin.tv_sec * 1000000 + tbegin.tv_usec)));
//         
        gettimeofday(&tbegin, NULL);
        for (int j = 0; j<nbtests[i]; j++)
            compApp_poly_parallel_taylor( pres, p, compRat_real_ptr(center), compRat_imag_ptr(center), radius, prec, 4);
        gettimeofday(&tend, NULL);
        printf ("time for %d taylor shifts, degree %4d: %10ld useconds.\n",nbtests[i], degs[i], ((tend.tv_sec * 1000000 + tend.tv_usec) - (tbegin.tv_sec * 1000000 + tbegin.tv_usec)));
        printf("\n");
    }
    
//     compApp_clear(c);
    compRat_clear(center);
    realRat_clear(radius);
    compApp_poly_clear(p);
    realRat_poly_clear(pbern);
    return 0;
}
