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
#include "polynomials/realRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "caches/cacheApp.h"

void getApp0 (compApp_poly_t dest, slong prec);
void getAppBernoulli (compApp_poly_t dest, slong prec);

int main() {
    
    int max = 15;
    cacheApp cache;
    
    cacheApp_initialize(&cache, getAppBernoulli);
    
    compApp_poly_t p;
    compApp_poly_init(p);
    getAppBernoulli(p,(0x1<<4)*DEFAULT_PREC);
    printf("p bern: "); compApp_poly_printd(p, 30); printf("\n\n");
    compApp_poly_clear(p);
    
    compApp_poly_ptr pApprox;
    for(int i=0; i<max; i++)
        pApprox = cacheApp_getApproximation(&cache, (0x1<<i)*DEFAULT_PREC);
    
    pApprox = cacheApp_getApproximation(&cache, (0x1<<5)*DEFAULT_PREC);
    printf("p bern: "); compApp_poly_printd(pApprox, 30); printf("\n\n");
    
    printf("cache._size: %d, cache._allocsize: %d \n", cache._size, cache._allocsize);
    
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<4)*DEFAULT_PREC);
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<5)*DEFAULT_PREC);
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<6)*DEFAULT_PREC);
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<17)*DEFAULT_PREC);
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<25)*DEFAULT_PREC);
//     
//     printf("cache._size: %d, cache._allocsize: %d \n", cache._size, cache._allocsize);
//     pApprox = cacheApp_getApproximation(&cache, (0x1<<25)*DEFAULT_PREC);
    
    cacheApp_flush(&cache);
    printf("cache._size: %d, cache._allocsize: %d \n", cache._size, cache._allocsize);
//     pApprox = cacheApp_getApproximation(&cache, DEFAULT_PREC);
    
    return 0;
}

// dest is assumed to be initialized
void getApp0 (compApp_poly_t dest, slong prec){
    compApp_poly_zero(dest);
}

void getAppBernoulli (compApp_poly_t dest, slong prec){
    printf("getAppBernoulli, prec: %d\n", prec);
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    bernoulli_polynomial(pbern , 5);
    compApp_poly_set_realRat_poly( dest, pbern, prec);
    realRat_poly_clear(pbern);
//     printf("getAppBernoulli end, prec: %d\n", prec);
}