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

#ifndef CAUCHYTESTS_H
#define CAUCHYTESTS_H


#include "base/base.h"
#include "numbers/realRat.h"
#include "numbers/realApp.h"
#include "numbers/app_rat.h"
#include "caches/cacheApp.h"
#include "caches/cacheCauchy.h"
#include "metadatas/metadatas.h"

#include "geometry/compDsk.h"

#include "acb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif
    
#define CAUCHYTEST_CERTIFIED 1
#define CAUCHYTEST_UNCERTIFI 0
#define CAUCHYTEST_INCOUNTIN 1
#define CAUCHYTEST_INEXCLUSI 0  
    
typedef struct {
    int nbOfSol;   /* the number of solutions: -1: can not decide, >=0 otherwise */
    slong appPrec; /* the arithmetic precision that has been used to decide      */
} cauchyTest_res;    

void cauchyTest_getEvaluationPoints( const compRat_t center,
                                     const realRat_t radius,
                                     const realRat_t radius2,
                                     slong vangle,
                                     slong vindex,
                                     cacheCauchy_t cacheCau,
                                     int certified,
                                     slong prec );

void cauchyTest_evaluateAtPoints( cacheApp_t cache,
                                  cacheCauchy_t cacheCau,
                                  int certified,
                                  slong prec,
                                  int inCounting,
                                  metadatas_t meta, int depth);

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeFdivs_fromVals(cacheCauchy_t cacheCau,
                                     int certified,
                                     slong prec,
                                     int inCounting,
                                     metadatas_t meta);

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeS0Approx_fromVals(compApp_t ps,
                                        cacheCauchy_t cacheCau,
                                        int certified,
                                        slong rotation,
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta);

/* returns -1: should increase precision
 *         -2: disk is has not expected isolation ratio; should stop
 *          1: OK! */
int cauchyTest_computeSsApprox_fromVals(compApp_ptr ps,
                                        cacheCauchy_t cacheCau,
//                                         int certified,
//                                         slong rotation,
                                        slong prec,
                                        int inCounting,
                                        metadatas_t meta);

cauchyTest_res cauchyTest_computeS0Approx(compApp_t ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,           
                                          slong rotation,
                                          int *alreadyEvaluated,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          int certified,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth );

cauchyTest_res cauchyTest_computeSsApprox(compApp_ptr ps,
                                          const compRat_t center,
                                          const realRat_t radius,
                                          const realRat_t radius2,
                                          slong vangle,           
                                          slong vindex,           
//                                           slong rotation,
//                                           int *alreadyEvaluated,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          int certified,
                                          slong prec,
                                          int inCounting,
                                          metadatas_t meta, int depth );

cauchyTest_res cauchyTest_deterministic_exclusion_test( const compRat_t center,
                                                       const realRat_t radius,
                                                       const realRat_t radius2,
                                                       slong vangle,           
                                                       slong vindex,           
                                                       cacheApp_t cache,
                                                       cacheCauchy_t cacheCau,
                                                       slong prec,
                                                       int inCounting,
                                                       metadatas_t meta, int depth);

cauchyTest_res cauchyTest_probabilistic_exclusion_test( const compRat_t center,
                                                       const realRat_t radius,
                                                       const realRat_t radius2,
                                                       slong vangle,           
                                                       slong vindex,           
                                                       cacheApp_t cache,
                                                       cacheCauchy_t cacheCau,
                                                       slong prec,
                                                       metadatas_t meta, int depth);

cauchyTest_res cauchyTest_probabilistic_counting_test( const compRat_t center,
                                                       const realRat_t radius,           
                                                       cacheApp_t cache,
                                                       cacheCauchy_t cacheCau,
                                                       slong prec,
                                                       metadatas_t meta, int depth);

/* version for separated input disk:
 * D(center, radius) and D(center, 4*radius) 
 * are assumed to contain the same number of roots */
cauchyTest_res cauchyTest_deterministic_counting_test_for_separated_discs( const compRat_t center,
                                                                           const realRat_t radius,           
                                                                           cacheApp_t cache,
                                                                           cacheCauchy_t cacheCau,
                                                                           slong prec,
                                                                           metadatas_t meta, int depth);

/* version for newton: */
cauchyTest_res cauchyTest_deterministic_counting_test_for_newton( const compRat_t center,
                                                                  realRat_t radInf,  
                                                                  const realRat_t radSup,
                                                                  slong nbOfRoots,
                                                                  cacheApp_t cache,
                                                                  cacheCauchy_t cacheCau,
                                                                  slong prec,
                                                                  metadatas_t meta, int depth);

cauchyTest_res cauchyTest_deterministic_counting_test_combinatorial( const compRat_t center,
                                                                     realRat_t radius,
                                                                     slong nbOfRoots,
                                                                     cacheApp_t cache,
                                                                     cacheCauchy_t cacheCau,
                                                                     slong prec,
                                                                     metadatas_t meta, int depth);

    
slong cauchyTest_computeS1compDsk( compDsk_t res,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   metadatas_t meta, int depth);
#ifdef __cplusplus
}
#endif

#endif
