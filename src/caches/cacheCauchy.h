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
#include "geometry/compDsk.h"
#include "metadatas/metadatas.h"

#include "caches/cacheApp.h"

#ifdef CCLUSTER_HAVE_PTHREAD
#include <pthread.h>
#endif

#define CACHE_DEFAULT_SIZE 10

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    cacheApp_ptr _cacheApp;
    
    void(*_evalFast)(compApp_t, compApp_t, const compApp_t, slong);
    slong   _degree;
    
    /* for sparse evaluation */
    slong    _nbNZC;    /* number of non-zero coefficients */
    slong  * _inNZC;  /* table containing the degrees of the non-zero coefficients 
                         when initialized, contains _degree+1 elements*/
    int      _choice; /* -1 per default;  */
    
    realRat _isoRatio; /* assumed isolation ratio */
    realApp _precfdiv; /* si* : 1/4 */
    
    /* uncertified exclusion */
    slong       _nbPwSuEx; /* nb of power sums computed for uncertified exclusion: >=1 */
    slong       _nbEvalEx;
    realApp_ptr _wanErrEx; /* table of length _nbPwSuEx; _wanErrCo[i] is error of s_i */
    
    realRat _lBoundUn;
    realRat _uBoundUn;
    
    realRat _curRadiu;
    slong   _precBoun;
    realApp _lBoundAp;
    realApp _uBoundAp;
    
    /* evaluation points for uncertified exclusion*/
    slong       _precEvalEx; /* precision of evaluation */
    compApp_ptr _pointsEx;
    compApp_ptr _pointsShiftedEx;
    compApp_ptr _fvalsEx;
    compApp_ptr _fdervalsEx;
    compApp_ptr _fdivsEx;
    
} cacheCauchy;

typedef cacheCauchy cacheCauchy_t[1];
typedef cacheCauchy * cacheCauchy_ptr;

#define cacheCauchy_cacheAppref(X)   (X->_cacheApp)
#define cacheCauchy_evalFastref(X)   (X->_evalFast)
#define cacheCauchy_degreeref(X)   (X->_degree)
#define cacheCauchy_nbNZCref(X)   (X->_nbNZC)
#define cacheCauchy_inNZCref(X)   (X->_inNZC)
#define cacheCauchy_choiceref(X)   (X->_choice)
#define cacheCauchy_isoRatioref(X) (&(X)->_isoRatio)
#define cacheCauchy_precfdivref(X) (&(X)->_precfdiv)
#define cacheCauchy_nbPwSuExref(X)   (X->_nbPwSuEx)
#define cacheCauchy_nbEvalExref(X) (X->_nbEvalEx)
#define cacheCauchy_wanErrExref(X) ((X)->_wanErrEx)
#define cacheCauchy_lBoundUnref(X) (&(X)->_lBoundUn)
#define cacheCauchy_uBoundUnref(X) (&(X)->_uBoundUn)
#define cacheCauchy_curRadiuref(X) (&(X)->_curRadiu)
#define cacheCauchy_precBounref(X) (X->_precBoun)
#define cacheCauchy_lBoundApref(X) (&(X)->_lBoundAp)
#define cacheCauchy_uBoundApref(X) (&(X)->_uBoundAp)

#define cacheCauchy_precEvalExref(X)      ((X)->_precEvalEx)
#define cacheCauchy_pointsExref(X)        ((X)->_pointsEx)
#define cacheCauchy_pointsShiftedExref(X) ((X)->_pointsShiftedEx)
#define cacheCauchy_fvalsExref(X)         ((X)->_fvalsEx)
#define cacheCauchy_fdervalsExref(X)      ((X)->_fdervalsEx)
#define cacheCauchy_fdivsExref(X)         ((X)->_fdivsEx)

/* compute q = ceil ( log_isoRatio (4*degree +1) ) + nbPs */
slong cacheCauchy_get_NbOfEvalPoints( slong degree, const realRat_t isoRatio, slong nbPs, slong prec );

/* compute error wP = (d*isoRatio^(-q))/(1-isoRatio^(-q)) */
void cacheCauchy_wantedErrorOnS0 (realApp_t wP, slong degree, slong q, const realRat_t isoRatio, slong prec );

/* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
void cacheCauchy_lowerBoundUnit( realRat_t lowerBoundUnit, slong degree, const realRat_t isoRatio );

/* compute upper bound unit: (d*(isoRatio+1)/(isoRatio-1) */
void cacheCauchy_upperBoundUnit( realRat_t upperBoundUnit, slong degree, const realRat_t isoRatio );

void cacheCauchy_lBoundApp( realApp_t lbApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec);
void cacheCauchy_uBoundApp( realApp_t ubApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec);

void cacheCauchy_sparseEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
                              const compApp_t point, slong prec);
void cacheCauchy_sparseEval2 ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
                              const compApp_t point, slong prec);

void cacheCauchy_rectangularEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache,
                                   const compApp_t point, slong prec);

void cacheCauchy_eval( compApp_ptr fvals, compApp_ptr fdervals, compApp_ptr points, slong nbPoints, 
                       cacheCauchy_t cache, slong prec );

void cacheCauchy_init ( cacheCauchy_t cache,
                        cacheApp_t cacheA,
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        const realRat_t isoRatio,
                        slong nbPows,
                        const metadatas_t meta
                      );

void cacheCauchy_get_lBoundApp( realApp_t lbApp, cacheCauchy_t cache, const realRat_t radius, slong prec);

void cacheCauchy_get_uBoundApp( realApp_t ubApp, cacheCauchy_t cache, const realRat_t radius, slong prec);

void cacheCauchy_set_bounds( cacheCauchy_t cache, const realRat_t radius, slong prec );

void cacheCauchy_clear ( cacheCauchy_t cache );

#ifdef __cplusplus
}
#endif

#endif
