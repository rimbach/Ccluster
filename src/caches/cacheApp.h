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

#ifndef CACHE_APP_H
#define CACHE_APP_H

#include "base/base.h"
#include "numbers/compApp.h"
#include "polynomials/compApp_poly.h"

#define CACHE_DEFAULT_SIZE 10
// #define DEFAULT_PREC 53

typedef struct {
    void(*_getApproximation)(compApp_poly_t, slong);
    compApp_poly_t *_cache;
    int _size;
    int _allocsize;
#ifdef CCLUSTER_EXPERIMENTAL
    compApp_poly_t **_cache_der; /* a table of tables caching derivatives */
    int * _der_size;
#endif
    /*for test: cache the last working polynomial computed and the nb of graeffe iterations*/
/*     compApp_poly _working;
       int _nbIterations; */
} cacheApp;

typedef cacheApp cacheApp_t[1];
typedef cacheApp * cacheApp_ptr;

/* #define cacheApp_workref(X) (&(X)->_working)    */
/* #define cacheApp_nbItref(X) (X->_nbIterations)  */

void cacheApp_init ( cacheApp_t cache, void(*getApproximation)(compApp_poly_t, slong) );

compApp_poly_ptr cacheApp_getApproximation ( cacheApp_t cache, slong prec );
slong cacheApp_getDegree ( cacheApp_t cache );

#ifdef CCLUSTER_EXPERIMENTAL
compApp_poly_ptr cacheApp_getDerivative ( cacheApp_t cache, slong prec, slong order );
#endif

void cacheApp_clear ( cacheApp_t cache );

#endif