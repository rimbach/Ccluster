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
// #define DEFAULT_PREC 53

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    
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
    
    /* certified exclusion */
    slong   _nbEvalCe;
    slong   _quotient; /* _nbEvalCe = _nbEvalCo * _quotient */
    realApp _wanErrCe;
    
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
    
    /* evaluation points for certification*/
    compApp_ptr _pointsCe;
    compApp_ptr _pointsShiftedCe;
    compApp_ptr _fvalsCe;
    compApp_ptr _fdervalsCe;
    compApp_ptr _fdivsCe;
    
//     /* evaluation points for compression*/
//     slong       _nbEvalComp;
//     compApp_ptr _pointsComp;
//     compApp_ptr _pointsShiftedComp;
//     compApp_ptr _fvalsComp;
//     compApp_ptr _fdervalsComp;
//     compApp_ptr _fdivsComp;
//     realApp_ptr _errorComp;
    
    /* case where the coefficients are known: shifted polynomial */
//     compApp_poly _shiftedPoly;
//     compApp_poly _shiftedPolyDer;
    
//     /* test with p mod (x^q -1) and p' mod (x^q -1) */
//     compRat_poly_t  _RP;
//     compRat_poly_t  _QP;
//     compRat_poly_t  _RPp;
//     compRat_poly_t  _QPp;
//     compApp_poly_t *_cacheRP;
//     compApp_poly_t *_cacheQP;
//     compApp_poly_t *_cacheRPp;
//     compApp_poly_t *_cacheQPp;
//     int _sizeCache;
//     int _allocsizeCache;
    
} cacheCauchy;

typedef cacheCauchy cacheCauchy_t[1];
typedef cacheCauchy * cacheCauchy_ptr;

#define cacheCauchy_evalFastref(X)   (X->_evalFast)
#define cacheCauchy_degreeref(X)   (X->_degree)
#define cacheCauchy_nbNZCref(X)   (X->_nbNZC)
#define cacheCauchy_inNZCref(X)   (X->_inNZC)
#define cacheCauchy_choiceref(X)   (X->_choice)
#define cacheCauchy_isoRatioref(X) (&(X)->_isoRatio)
#define cacheCauchy_precfdivref(X) (&(X)->_precfdiv)
#define cacheCauchy_nbPwSuExref(X)   (X->_nbPwSuEx)
#define cacheCauchy_nbEvalExref(X) (X->_nbEvalEx)
#define cacheCauchy_nbEvalCeref(X) (X->_nbEvalCe)
#define cacheCauchy_quotientref(X) (X->_quotient)
#define cacheCauchy_wanErrExref(X) ((X)->_wanErrEx)
#define cacheCauchy_wanErrCeref(X) (&(X)->_wanErrCe)
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

#define cacheCauchy_pointsCeref(X)        ((X)->_pointsCe)
#define cacheCauchy_pointsShiftedCeref(X) ((X)->_pointsShiftedCe)
#define cacheCauchy_fvalsCeref(X)         ((X)->_fvalsCe)
#define cacheCauchy_fdervalsCeref(X)      ((X)->_fdervalsCe)
#define cacheCauchy_fdivsCeref(X)         ((X)->_fdivsCe)

// #define cacheCauchy_nbEvalCompref(X)         ((X)->_nbEvalComp)
// #define cacheCauchy_pointsCompref(X)         ((X)->_pointsComp)
// #define cacheCauchy_pointsShiftedCompref(X)  ((X)->_pointsShiftedComp)
// #define cacheCauchy_fvalsCompref(X)          ((X)->_fvalsComp)
// #define cacheCauchy_fdervalsCompref(X)       ((X)->_fdervalsComp)
// #define cacheCauchy_fdivsCompref(X)          ((X)->_fdivsComp)
// #define cacheCauchy_errorCompref(X)          ((X)->_errorComp)

// #define cacheCauchy_shiftedPolyref(X)      (&(X)->_shiftedPoly)
// #define cacheCauchy_shiftedPolyDerref(X)   (&(X)->_shiftedPolyDer)

/* compute q = ceil ( log_isoRatio (4*degree +1) ) + nbPs */
slong cacheCauchy_get_NbOfEvalPoints( slong degree, const realRat_t isoRatio, slong nbPs, slong prec );

/* compute q2 = log_isoRatio (2*degree*q1 +1) s.t. q2 multiple of q1 */
slong cacheCauchy_get_NbOfEvalPoints_cert( slong degree, slong q1, const realRat_t isoRatio, slong prec );

/* compute error wP = (d*isoRatio^(-q))/(1-isoRatio^(-q)) */
void cacheCauchy_wantedErrorOnS0 (realApp_t wP, slong degree, slong q, const realRat_t isoRatio, slong prec );

/* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
void cacheCauchy_lowerBoundUnit( realRat_t lowerBoundUnit, slong degree, const realRat_t isoRatio );

/* compute upper bound unit: (d*(isoRatio+1)/(isoRatio-1) */
void cacheCauchy_upperBoundUnit( realRat_t upperBoundUnit, slong degree, const realRat_t isoRatio );

void cacheCauchy_lBoundApp( realApp_t lbApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec);
void cacheCauchy_uBoundApp( realApp_t ubApp, slong degree, const realRat_t isoRatio, const realRat_t radius, slong prec);

void cacheCauchy_init_sparseEval ( cacheCauchy_t cache, 
                                   cacheApp_t cachePol );

void cacheCauchy_sparseEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, cacheApp_t cachePol,
                              const compApp_t point, slong prec);
void cacheCauchy_sparseEval2 ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, cacheApp_t cachePol,
                              const compApp_t point, slong prec);

void cacheCauchy_rectangularEval ( compApp_t fval, compApp_t fderval, cacheCauchy_t cache, cacheApp_t cachePol,
                                   const compApp_t point, slong prec);

void cacheCauchy_init ( cacheCauchy_t cache, 
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        const realRat_t isoRatio,
                        slong nbPows,
                        const metadatas_t meta
                      );

void cacheCauchy_get_lBoundApp( realApp_t lbApp, cacheCauchy_t cache, const realRat_t radius, slong prec);

void cacheCauchy_get_uBoundApp( realApp_t ubApp, cacheCauchy_t cache, const realRat_t radius, slong prec);

void cacheCauchy_set_bounds( cacheCauchy_t cache, const realRat_t radius, slong prec );

// compApp_poly_ptr cacheCauchy_getApproximation_RP      ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec );
// compApp_poly_ptr cacheCauchy_getApproximation_QP      ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec );
// compApp_poly_ptr cacheCauchy_getApproximation_RPp     ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec );
// compApp_poly_ptr cacheCauchy_getApproximation_QPp     ( cacheCauchy_t cacheCau, cacheApp_t cache, slong prec );

// void cacheCauchy_init_comp ( cacheCauchy_t cache, const realRat_t isoRatio, const compDsk_t Delta, slong nbOfRoots, const realRat_t eps );

// void cacheCauchy_clear_comp ( cacheCauchy_t cache );

void cacheCauchy_clear ( cacheCauchy_t cache );

#ifdef __cplusplus
}
#endif

#endif
