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
    
//     compApp_poly_init(cacheApp_workref(cache));
//     cache->_nbIterations = 0;
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
            cache->_size +=1;
        }
        return (cache->_cache)[index];
    }
    
    while (index >= cache->_allocsize) 
        cache->_allocsize += CACHE_DEFAULT_SIZE;
    
    cache->_cache = (compApp_poly_t *) realloc (cache->_cache, (cache->_allocsize) * sizeof(compApp_poly_t) );
    while (index >= cache->_size){
        compApp_poly_init(cache->_cache[cache->_size]);
        slong nprec = (0x1<<(cache->_size))*CCLUSTER_DEFAULT_PREC;
            cache->_getApproximation( cache->_cache[cache->_size], nprec);
        cache->_size +=1;
    }
    return (cache->_cache)[index];
    
}

slong cacheApp_getDegree ( cacheApp_t cache ){
    if (cache->_size == 0)
        cacheApp_getApproximation (cache, CCLUSTER_DEFAULT_PREC);
    return compApp_poly_degree((cache->_cache)[0]);
}

void cacheApp_clear ( cacheApp_t cache ) {
    for (int i=0; i<cache->_size; i++)
        compApp_poly_clear( (cache->_cache)[i] );
    free(cache->_cache);
    cache->_size      = 0;
    cache->_allocsize = 0;
    
//     compApp_poly_clear(cacheApp_workref(cache));
}