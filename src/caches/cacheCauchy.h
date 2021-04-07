/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#ifndef CACHE_CAUCHY_H
#define CACHE_CAUCHY_H

#ifdef CACHE_INLINE_C
#define CACHE_INLINE
#else
#define CACHE_INLINE static __inline__
#endif

#include "base/base.h"
#include "numbers/compApp.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/realApp_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"

#ifdef CCLUSTER_HAVE_PTHREAD
#include <pthread.h>
#endif

#define CACHE_DEFAULT_SIZE 10
// #define DEFAULT_PREC 53

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    
    void(*_evalFast)(compApp_t, compApp_t, const compApp_t, slong);
    slong   _degree;
    
    realRat _isoRatio;
    
    slong   _nbEvalCo;
    slong   _nbEvalCe;
    slong   _quotient; /* _nbEvalCe = _nbEvalCo * _quotient */
    
    realApp _wanErrCo;
    realApp _wanErrCe;
    
    realRat _lBoundUn;
    realRat _uBoundUn;
    
    realRat _curRadiu;
    slong   _precBoun;
    realApp _lBoundAp;
    realApp _uBoundAp;
    
} cacheCauchy;

typedef cacheCauchy cacheCauchy_t[1];
typedef cacheCauchy * cacheCauchy_ptr;

#define cacheCauchy_evalFastref(X)   (X->_evalFast)
#define cacheCauchy_degreeref(X)   (X->_degree)
#define cacheCauchy_isoRatioref(X) (&(X)->_isoRatio)
#define cacheCauchy_nbEvalCoref(X) (X->_nbEvalCo)
#define cacheCauchy_nbEvalCeref(X) (X->_nbEvalCe)
#define cacheCauchy_quotientref(X) (X->_quotient)
#define cacheCauchy_wanErrCoref(X) (&(X)->_wanErrCo)
#define cacheCauchy_wanErrCeref(X) (&(X)->_wanErrCe)
#define cacheCauchy_lBoundUnref(X) (&(X)->_lBoundUn)
#define cacheCauchy_uBoundUnref(X) (&(X)->_uBoundUn)

void cacheCauchy_init ( cacheCauchy_t cache, 
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        slong num,
                        ulong den
                      );

void cacheCauchy_clear ( cacheCauchy_t cache );

// void cacheCauchy_get_lBoundAp( realApp_t lb, 

#ifdef __cplusplus
}
#endif

#endif
