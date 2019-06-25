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
#include "caches/cacheApp.h"
// #include <time.h>
#include <sys/time.h>

#include <pthread.h>

#define DEGREE 512
#define PREC 424

typedef struct {
//     realRat_poly_t pol;
    cacheApp_ptr cache;
    compApp_t c;
    int nb_TS;
} parallel_test_t;

void * _parallel_test_worker( void * arg_ptr );

void getAppBernoulli (compApp_poly_t dest, slong prec);

int main() {
    
    int nbThreads = 1;
    int nbTS = 10000;
//     slong prec = 424;
//     slong degree = DEGREE;
    
    cacheApp_t cache;
    cacheApp_init(cache, getAppBernoulli);
    
    compRat_t center;
    compRat_init(center);
    compRat_set_sisi(center, 1,1,1,1);
//     realRat_set_si(radius, 1,1);
    
    parallel_test_t * args = (parallel_test_t *) malloc ( sizeof(parallel_test_t) * nbThreads );
    pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nbThreads);
    
    
    for (int i = 0; i< (int) nbThreads; i++) {
        
//         realRat_poly_init(args[i].pol);
//         bernoulli_polynomial(args[i].pol, DEGREE); 
        args[i].cache = (cacheApp_ptr) cache;
        compApp_init(args[i].c);
        compApp_setreal_realRat(args[i].c, compRat_realref(center), PREC);
        compApp_setimag_realRat(args[i].c, compRat_imagref(center), PREC);
        args[i].nb_TS = (int) nbTS/nbThreads;
        /* create the thread */
        pthread_create(&threads[i], NULL, _parallel_test_worker, &args[i]);
    }
    
    for(int i = 0; i< (int) nbThreads; i++) {
        pthread_join(threads[i], NULL);
//         realRat_poly_clear(args[i].pol);
        compApp_clear(args[i].c);
    }
    
    printf("ici!\n");
    
//     compApp_poly_t p;
//     realRat_poly_t pbern;
//     compRat_t center;
//     realRat_t radius;
//     compApp_t c;
//     
//     compApp_poly_init(p);
//     realRat_poly_init(pbern);
//     compRat_init(center);
//     realRat_init(radius);
//     compApp_init(c);
//     
//     compRat_set_sisi(center, 1,1,1,1);
//     realRat_set_si(radius, 1,1);
//     bernoulli_polynomial(pbern, 160);
//     
//     compApp_poly_set_realRat_poly( p, pbern, 53);
//     compApp_setreal_realRat(c, compRat_realref(center), 53);
//     compApp_setimag_realRat(c, compRat_imagref(center), 53);
//     _acb_poly_taylor_shift_convolution(p->coeffs, c, p->length, 53);
    
//     compApp_poly_clear(p);
//     realRat_poly_clear(pbern);
    compRat_clear(center);
//     realRat_clear(radius);
//     compApp_clear(c);
    
    cacheApp_clear(cache);
    
    free(args);
    free(threads);
    
    return 0;
}

void getAppBernoulli (compApp_poly_t dest, slong prec){
//     printf("getAppBernoulli, prec: %d\n", prec);
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    bernoulli_polynomial(pbern , DEGREE);
    compApp_poly_set_realRat_poly( dest, pbern, prec);
    realRat_poly_clear(pbern);
//     printf("getAppBernoulli end, prec: %d\n", prec);
}

void * _parallel_test_worker( void * arg_ptr ){
    
    parallel_test_t * arg = (parallel_test_t *) arg_ptr;
    
    /* arg->status has been set to 1 by caller           */
    /* nb_thread_running has been incremented by caller; */
    
//     compApp_poly_t p;
//     compApp_poly_init2(p, cacheApp_getDegree(arg->cache) +1);
    
//     slong prec = 424;
    
    for (int i=0;i< arg->nb_TS; i++){
        compApp_poly_t p;
        compApp_poly_init2(p, cacheApp_getDegree(arg->cache) +1);
//         compApp_poly_set_realRat_poly( p, arg->pol, PREC);
        compApp_poly_set(p, cacheApp_getApproximation ( arg->cache, PREC));
        _acb_poly_taylor_shift_convolution(p->coeffs, arg->c, p->length, PREC);
        compApp_poly_clear(p);
        
    }
//     compApp_poly_clear(p);
    flint_cleanup();
    
    return NULL;
    
}
