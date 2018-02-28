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
    cache->_cache            = (compApp_poly_t *) malloc ( (cache->_allocsize) * sizeof(compApp_poly_t) );
    cache->_getApproximation = getApproximation;
#ifdef CCLUSTER_EXPERIMENTAL    
    cache->_der_size         = (int *) malloc ( (cache->_allocsize) * sizeof(int) );
    cache->_cache_der        = (compApp_poly_t **) malloc ( (cache->_allocsize) * sizeof(compApp_poly_t *) );
#endif
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(cache->_mutex), NULL);
#endif
/*     compApp_poly_init(cacheApp_workref(cache));    */
/*     cache->_nbIterations = 0;                      */
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
            cache->_getApproximation( cache->_cache[cache->_size], nprec);
//             printf("end call\n");

#ifdef CCLUSTER_EXPERIMENTAL            
            /* initialize cache->_cache_der[cache->_size] */
            slong len = ((cache->_cache)[0])->length;
            cache->_cache_der[cache->_size] = (compApp_poly_t *) malloc ( (len-1) * sizeof(compApp_poly_t) );
            for (slong i = 1; i<len; i++)
                compApp_poly_init2((cache->_cache_der)[cache->_size][i-1], len-i);
            (cache->_der_size)[cache->_size] = 0;
            /* end initialize */
#endif            
            cache->_size +=1;
        }
        return (cache->_cache)[index];
    }
    
    while (index >= cache->_allocsize) 
        cache->_allocsize += CACHE_DEFAULT_SIZE;
    
    cache->_cache = (compApp_poly_t *) realloc (cache->_cache, (cache->_allocsize) * sizeof(compApp_poly_t) );
#ifdef CCLUSTER_EXPERIMENTAL     
    /* realloc size for cache->_cache_der*/
    cache->_cache_der = (compApp_poly_t **) realloc (cache->_cache_der, (cache->_allocsize) * sizeof(compApp_poly_t *) );
    cache->_der_size = (int *) realloc (cache->_der_size, (cache->_allocsize) * sizeof(int) );
#endif    
    while (index >= cache->_size){
        compApp_poly_init2(cache->_cache[cache->_size], compApp_poly_degree((cache->_cache)[0])+1);
        slong nprec = (0x1<<(cache->_size))*CCLUSTER_DEFAULT_PREC;
            cache->_getApproximation( cache->_cache[cache->_size], nprec);
        
#ifdef CCLUSTER_EXPERIMENTAL
        /* initialize cache->_cache_der[cache->_size] */
        slong len = ((cache->_cache)[0])->length;
        cache->_cache_der[cache->_size] = (compApp_poly_t *) malloc ( (len-1) * sizeof(compApp_poly_t) );
        slong i = 1;
        for (i = 1; i<len; i++)
            compApp_poly_init2(cache->_cache_der[cache->_size][i-1], len-i);
        cache->_der_size[cache->_size] = 0;
        /* end initialize */
#endif            
        cache->_size +=1;
    }
    return (cache->_cache)[index];
    
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
    if (cache->_size == 0)
        cacheApp_getApproximation (cache, CCLUSTER_DEFAULT_PREC);
    return compApp_poly_degree((cache->_cache)[0]);
}

void cacheApp_clear ( cacheApp_t cache ) {

#ifdef CCLUSTER_EXPERIMENTAL    
    slong len = 0;
    if (cache->_size>=1)
        len = ((cache->_cache)[0])->length;
#endif
    
    for (int i=0; i<cache->_size; i++)
        compApp_poly_clear( (cache->_cache)[i] );
    free(cache->_cache);

#ifdef CCLUSTER_EXPERIMENTAL
//     printf("ici, len: %d\n", (int) len);
    for (int i=0; i<cache->_size; i++) {
//         printf("i: %d\n", i);
        for (slong j=0; j<len-1; j++) {
//             printf("j: %d\n", (int) j);
            compApp_poly_clear( (cache->_cache_der)[i][(int) j] );
        }
        free((cache->_cache_der)[i]);
    }
    free(cache->_cache_der);
    
    free(cache->_der_size);
#endif    
    cache->_size      = 0;
    cache->_allocsize = 0;
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_destroy( &(cache->_mutex) );
#endif
    
//     compApp_poly_clear(cacheApp_workref(cache));
}