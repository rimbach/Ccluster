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

#include <stdlib.h>
#include "cacheApp.h"

void cacheApp_init ( cacheApp_t cache, void(*getApproximation)(compApp_poly_t, slong) ) {
    cache->_size             = 0;
    cache->_allocsize        = CACHE_DEFAULT_SIZE;
    cache->_cache            = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsize) * sizeof(compApp_poly_t) );
    cache->_getApproximation = getApproximation;
    
    cache->_size_real        = 0;
    cache->_allocsize_real   = CACHE_DEFAULT_SIZE;
    cache->_cache_real       = (realApp_poly_t *) ccluster_malloc ( (cache->_allocsize_real) * sizeof(realApp_poly_t) );
    
    cache->_from_poly = 0;
    
    cache->_degree = -1;
#ifdef CCLUSTER_EXPERIMENTAL    
    cache->_der_size         = (int *) ccluster_malloc ( (cache->_allocsize) * sizeof(int) );
    cache->_cache_der        = (compApp_poly_t **) ccluster_malloc ( (cache->_allocsize) * sizeof(compApp_poly_t *) );
#endif
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(cache->_mutex), NULL);
#endif                   
}

void cacheApp_init_compRat_poly ( cacheApp_t cache, const compRat_poly_t poly){
    cache->_size             = 0;
    cache->_allocsize        = CACHE_DEFAULT_SIZE;
    cache->_cache            = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsize) * sizeof(compApp_poly_t) );
    cache->_getApproximation = NULL;
    
    cache->_size_real        = 0;
    cache->_allocsize_real   = CACHE_DEFAULT_SIZE;
    cache->_cache_real       = (realApp_poly_t *) ccluster_malloc ( (cache->_allocsize_real) * sizeof(realApp_poly_t) );
    
    compRat_poly_init(cache->_poly);
    compRat_poly_set(cache->_poly, poly);
    compRat_poly_canonicalise(cache->_poly);
    cache->_from_poly = 1;
    
    cache->_degree = -1;
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(cache->_mutex), NULL);
#endif
}

void cacheApp_init_realRat_poly ( cacheApp_t cache, const realRat_poly_t poly){
    cache->_size             = 0;
    cache->_allocsize        = CACHE_DEFAULT_SIZE;
    cache->_cache            = (compApp_poly_t *) ccluster_malloc ( (cache->_allocsize) * sizeof(compApp_poly_t) );
    cache->_getApproximation = NULL;
    
    cache->_size_real        = 0;
    cache->_allocsize_real   = CACHE_DEFAULT_SIZE;
    cache->_cache_real       = (realApp_poly_t *) ccluster_malloc ( (cache->_allocsize_real) * sizeof(realApp_poly_t) );
    
    compRat_poly_init(cache->_poly);
    compRat_poly_set_realRat_poly(cache->_poly, poly);
    compRat_poly_canonicalise(cache->_poly);
    
//     printf("\n\n"); compRat_poly_print(cache->_poly); printf("\n\n");
//     fmpz_one( fmpq_poly_denref( compRat_poly_realref(cache->_poly )) );
//     printf("\n\n"); compRat_poly_print(cache->_poly); printf("\n\n");
    
    cache->_from_poly = 1;
    
    cache->_degree = -1;
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(cache->_mutex), NULL);
#endif
}

//requires: prec is 2^n*CCLUSTER_DEFAULT_PREC
compApp_poly_ptr cacheApp_getApproximation ( cacheApp_t cache, slong prec ) {
    //get index in cache
    slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
    int index = 0;
    while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
//     printf("index: %d\n", index); 
    
    if (index < cache->_size)
        return (cache->_cache)[index];

#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_lock(cache);
#endif
    
    if (index < cache->_allocsize) {
        while (index >= cache->_size){
//             printf("initialize %d\n", cache->_size);
            if (cache->_size>=1)
                compApp_poly_init2(cache->_cache[cache->_size], compApp_poly_degree((cache->_cache)[0])+1);
            else
                compApp_poly_init(cache->_cache[cache->_size]);
//             printf("end initialize %d\n", cache->_size);
            slong nprec = (0x1<<(cache->_size))*CCLUSTER_DEFAULT_PREC;
//             printf("call with prec: %d\n", nprec);
            if (cache->_from_poly==0){
                cache->_getApproximation( cache->_cache[cache->_size], nprec);
            }
            else {
                compApp_poly_set_compRat_poly(cache->_cache[cache->_size], cache->_poly, nprec);
            }
//             printf("end call\n");

#ifdef CCLUSTER_EXPERIMENTAL            
            /* initialize cache->_cache_der[cache->_size] */
            slong len = ((cache->_cache)[0])->length;
            cache->_cache_der[cache->_size] = (compApp_poly_t *) ccluster_malloc ( (len-1) * sizeof(compApp_poly_t) );
            for (slong i = 1; i<len; i++)
                compApp_poly_init2((cache->_cache_der)[cache->_size][i-1], len-i);
            (cache->_der_size)[cache->_size] = 0;
            /* end initialize */
#endif            
            cache->_size +=1;
        }
#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_unlock(cache);
#endif
        return (cache->_cache)[index];
    }
    
    while (index >= cache->_allocsize) 
        cache->_allocsize += CACHE_DEFAULT_SIZE;
    
    cache->_cache = (compApp_poly_t *) ccluster_realloc (cache->_cache, (cache->_allocsize) * sizeof(compApp_poly_t) );
#ifdef CCLUSTER_EXPERIMENTAL     
    /* realloc size for cache->_cache_der*/
    cache->_cache_der = (compApp_poly_t **) ccluster_realloc (cache->_cache_der, (cache->_allocsize) * sizeof(compApp_poly_t *) );
    cache->_der_size = (int *) ccluster_realloc (cache->_der_size, (cache->_allocsize) * sizeof(int) );
#endif    
    while (index >= cache->_size){
        compApp_poly_init2(cache->_cache[cache->_size], compApp_poly_degree((cache->_cache)[0])+1);
        slong nprec = (0x1<<(cache->_size))*CCLUSTER_DEFAULT_PREC;
        
        if (cache->_from_poly==0){
            cache->_getApproximation( cache->_cache[cache->_size], nprec);
        }
        else {
            compApp_poly_set_compRat_poly(cache->_cache[cache->_size], cache->_poly, nprec);
        }
        
#ifdef CCLUSTER_EXPERIMENTAL
        /* initialize cache->_cache_der[cache->_size] */
        slong len = ((cache->_cache)[0])->length;
        cache->_cache_der[cache->_size] = (compApp_poly_t *) ccluster_malloc ( (len-1) * sizeof(compApp_poly_t) );
        slong i = 1;
        for (i = 1; i<len; i++)
            compApp_poly_init2(cache->_cache_der[cache->_size][i-1], len-i);
        cache->_der_size[cache->_size] = 0;
        /* end initialize */
#endif            
        cache->_size +=1;
    }
#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_unlock(cache);
#endif
    return (cache->_cache)[index];
    
}

//requires: prec is 2^n*CCLUSTER_DEFAULT_PREC
realApp_poly_ptr cacheApp_getApproximation_real ( cacheApp_t cache, slong prec ) {
    //get index in cache
    slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
    int index = 0;
    while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
    
    if (index < cache->_size_real)
        return (cache->_cache_real)[index];
    
    compApp_poly_t temp;
    if (cache->_from_poly==0) {
//         if (cache->_size_real>=1)
//             compApp_poly_init2( temp, realApp_poly_degree((cache->_cache_real)[0])+1 );
//         else
            compApp_poly_init( temp );
    }

#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_lock(cache);
#endif
    
    if (index < cache->_allocsize_real) {
        while (index >= cache->_size_real){
//             printf("initialize %d\n", cache->_size_real);
            if (cache->_size_real>=1)
                realApp_poly_init2(cache->_cache_real[cache->_size_real], realApp_poly_degree((cache->_cache_real)[0])+1);
            else
                realApp_poly_init(cache->_cache_real[cache->_size_real]);
//             printf("end initialize %d\n", cache->_size);
            slong nprec = (0x1<<(cache->_size_real))*CCLUSTER_DEFAULT_PREC;
//             printf("call with prec: %d\n", nprec);
            if (cache->_from_poly==0){
                cache->_getApproximation( temp, nprec);
                realApp_poly_fit_length(cache->_cache_real[cache->_size_real], temp->length);
                realApp_poly_set_length(cache->_cache_real[cache->_size_real], temp->length);
                for (slong i = 0; i < temp->length; i++)
                    realApp_set( (cache->_cache_real[cache->_size_real])->coeffs + i, compApp_realref(temp->coeffs + i));
            }
            else {
                realApp_poly_set_realRat_poly(cache->_cache_real[cache->_size_real], compRat_poly_realref(cache->_poly), nprec);
            }
//             printf("end call\n");         
            cache->_size_real +=1;
        }
        
#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_unlock(cache);
#endif
    if (cache->_from_poly==0)
        compApp_poly_clear( temp );
    
    return (cache->_cache_real)[index];
    }
    
    while (index >= cache->_allocsize_real) 
        cache->_allocsize_real += CACHE_DEFAULT_SIZE;
    
    cache->_cache_real = (realApp_poly_t *) ccluster_realloc (cache->_cache_real, (cache->_allocsize_real) * sizeof(realApp_poly_t) );
  
    while (index >= cache->_size_real){
        realApp_poly_init2(cache->_cache_real[cache->_size_real], realApp_poly_degree((cache->_cache_real)[0])+1);
        slong nprec = (0x1<<(cache->_size_real))*CCLUSTER_DEFAULT_PREC;
        
        if (cache->_from_poly==0){
            cache->_getApproximation( temp, nprec);
            realApp_poly_fit_length(cache->_cache_real[cache->_size_real], temp->length);
            realApp_poly_set_length(cache->_cache_real[cache->_size_real], temp->length);
            for (slong i = 0; i < temp->length; i++)
                realApp_set( (cache->_cache_real[cache->_size_real])->coeffs + i, compApp_realref(temp->coeffs + i));
        }
        else {
            realApp_poly_set_realRat_poly(cache->_cache_real[cache->_size_real], compRat_poly_realref(cache->_poly), nprec);
        }
        
           
        cache->_size_real +=1;
    }
#ifdef CCLUSTER_HAVE_PTHREAD 
    cacheApp_unlock(cache);
#endif
    
    if (cache->_from_poly==0)
        compApp_poly_clear( temp );
    
    return (cache->_cache_real)[index];
    
}

#ifdef CCLUSTER_EXPERIMENTAL
compApp_poly_ptr cacheApp_getDerivative ( cacheApp_t cache, slong prec, slong order ) {
    
    //get index in cache
    slong log2prec = (slong)(prec/(slong)CCLUSTER_DEFAULT_PREC);
    int index = 0;
    while (log2prec>>=1) index++; //index should contain the log2 of prec/CCLUSTER_DEFAULT_PREC
    
    if (index >= cache->_size)
        cacheApp_getApproximation ( cache, prec );
    
    if ( order == 0 )
        return cacheApp_getApproximation ( cache, prec );
    
    order-=1;
    if ( ((int) order ) < cache->_der_size[index] )
        return cache->_cache_der[index][(int) order];
    
    while ( ((int) order) >= cache->_der_size[index] ) {
//         printf("compute derivative: prec: %d, order: %d\n", (int) prec, (int) order);
        if ( cache->_der_size[index] ==0 )
            compApp_poly_derivative( cache->_cache_der[index][cache->_der_size[index]], cache->_cache[index], prec);
        else 
            compApp_poly_derivative( cache->_cache_der[index][cache->_der_size[index]], cache->_cache_der[index][cache->_der_size[index]-1], prec);
        cache->_der_size[index]+=1;
    }
    
    return cache->_cache_der[index][(int) order];
    
}
#endif

slong cacheApp_getDegree ( cacheApp_t cache ){
    
    if (cache->_degree == -1) {
        
        if (cache->_from_poly) {
            cache->_degree = compRat_poly_degree( cache->_poly );
        }
        
        else {
            if (cache->_size == 0)
                cacheApp_getApproximation (cache, CCLUSTER_DEFAULT_PREC);
            cache->_degree = compApp_poly_degree((cache->_cache)[0]);
        }
        
    }
    
    return cache->_degree;
}

/* returns the number of zero trailing coeffs */
/* it is assumed that _poly has integer coefficients */
slong cacheApp_getMultOfZero ( cacheApp_t cache ){
    realRat_poly_ptr p = compRat_poly_real_ptr( cache->_poly );
    int res = 0;
    while ( (res < realRat_poly_length(p)) && (fmpz_is_zero(p->coeffs + res) ==1 )) 
        res++;
    return res;
}

int cacheApp_is_zero ( cacheApp_t cache ){
    
    if (cache->_from_poly)
        return compRat_poly_is_zero(cache->_poly);
        
    return 0;
}

int cacheApp_is_near_zero ( cacheApp_t cache ){
    
    if (cache->_size == 0)
        cacheApp_getApproximation (cache, CCLUSTER_DEFAULT_PREC);
    
    return compApp_poly_contains_zero((cache->_cache)[0]);
    
}

int cacheApp_is_real ( cacheApp_t cache ){
    
    if (cache->_size == 0)
        cacheApp_getApproximation (cache, CCLUSTER_DEFAULT_PREC);
    return compApp_poly_is_real((cache->_cache)[0]);
}

void cacheApp_root_bound ( realRat_t bound, cacheApp_t cache ){
    
    compApp_poly_ptr p;
    slong deg = cacheApp_getDegree(cache);
    
    /* deal with degree 0 polynomials */
    if (deg<=0) {
        realRat_one(bound);
        return;
    }
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    p = cacheApp_getApproximation (cache, prec);
    
    while ( compApp_contains_zero( compApp_poly_getCoeff(p, deg) ) ) {
        prec = 2*prec;
        p = cacheApp_getApproximation (cache, prec);
    }
    
    compApp_poly_root_bound_fujiwara( bound, p);
}

/* return 1 if no problem */
/* return 0 if lcf vanishes */
int cacheApp_root_bound_unsure ( realRat_t bound, cacheApp_t cache){
    
    int res=1;
    
    compApp_poly_ptr p;
    slong deg = cacheApp_getDegree(cache);
    slong prec = CCLUSTER_DEFAULT_PREC;
    p = cacheApp_getApproximation (cache, prec);
    
    while ( compApp_contains_zero( compApp_poly_getCoeff(p, deg) ) 
            && ( prec <= CCLUSTER_THRESHOLD_FOR_GLOBAL ) ) {
        prec = 2*prec;
        p = cacheApp_getApproximation (cache, prec);
    }
    
    if (compApp_contains_zero( compApp_poly_getCoeff(p, deg))) {
        realRat_set_si( bound, 2, 1 );
        realRat_pow_si( bound, bound, CCLUSTER_THRESHOLD_FOR_GLOBAL );
        res = 0;
//         printf("WARNING FROM CCLUSTER, cacheApp.c, line 232: leading coeff is less than 
    }
    else
        compApp_poly_root_bound_fujiwara( bound, p);
    
    return res;
}

slong cacheApp_getBitsize (cacheApp_t cache){
    
    fmpz_poly_t num;
    fmpz_poly_init(num);
    
    realRat_poly_canonicalise(compRat_poly_realref(cache->_poly));
    fmpq_poly_get_numerator(num, compRat_poly_realref(cache->_poly) );
    /* compute the inverse of the norm of the poly*/
    slong res = fmpz_poly_max_bits(num);
    
    fmpz_poly_clear(num);
    
    if (res<0) res = -res + 1;
    return res;
}

void cacheApp_separation_bound ( realRat_t sepBound, cacheApp_t cache){
    
    slong deg = cacheApp_getDegree ( cache );
    
    /* deal with degree 0 polynomials */
    if (deg==0) {
        realRat_one(sepBound);
        return;
    }
    
    realRat_poly_canonicalise(compRat_poly_realref(cache->_poly));
    /* compute the inverse of the norm of the poly*/
    fmpz_t norm2;
    fmpz_poly_t num;
    fmpz_init(norm2);
    fmpz_poly_init(num);
    fmpq_poly_get_numerator(num, compRat_poly_realref(cache->_poly) );
    fmpz_poly_2norm(norm2, num);
    fmpz_pow_ui(norm2, norm2, deg-1);
    fmpz_set_ui(realRat_numref(sepBound), 1); 
    fmpz_set(realRat_denref(sepBound), norm2);
//     printf("norm: "); fmpz_print(norm2); printf("\n");
    ulong exp = deg;
    if (exp%2==1)
        exp = exp+1;
    exp = (exp+2)/2;
    fmpz_set_si(norm2, deg);
    fmpz_pow_ui(norm2, norm2, exp);
    fmpq_div_fmpz(sepBound, sepBound, norm2);
    
    fmpz_clear(norm2);
    fmpz_poly_clear(num);
}

void cacheApp_clear ( cacheApp_t cache ) {

#ifdef CCLUSTER_EXPERIMENTAL    
    slong len = 0;
    if (cache->_size>=1)
        len = ((cache->_cache)[0])->length;
#endif
    
    for (int i=0; i<cache->_size; i++)
        compApp_poly_clear( (cache->_cache)[i] );
    ccluster_free(cache->_cache);
    
    for (int i=0; i<cache->_size_real; i++)
        realApp_poly_clear( (cache->_cache_real)[i] );
    ccluster_free(cache->_cache_real);

#ifdef CCLUSTER_EXPERIMENTAL
//     printf("ici, len: %d\n", (int) len);
    for (int i=0; i<cache->_size; i++) {
//         printf("i: %d\n", i);
        for (slong j=0; j<len-1; j++) {
//             printf("j: %d\n", (int) j);
            compApp_poly_clear( (cache->_cache_der)[i][(int) j] );
        }
        ccluster_free((cache->_cache_der)[i]);
    }
    ccluster_free(cache->_cache_der);
    
    ccluster_free(cache->_der_size);
#endif    
    cache->_size      = 0;
    cache->_allocsize = 0;
    
    cache->_size_real      = 0;
    cache->_allocsize_real = 0;
    
    if (cache->_from_poly==1) {
        compRat_poly_clear(cache->_poly);
    }
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_destroy( &(cache->_mutex) );
#endif
    
}
