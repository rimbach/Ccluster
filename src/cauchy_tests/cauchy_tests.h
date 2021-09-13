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

#include <acb_dft.h>

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

/* computes approximation sh of h-th powersum */
void cauchyTest_computeShApprox_fromVals(compApp_t sh,
                                         slong h,
                                         cacheCauchy_t cacheCau,
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

cauchyTest_res cauchyTest_deterministic_exclusion_testNEW( const compRat_t center,
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

cauchyTest_res cauchyTest_probabilistic_verification( const compDsk_t Delta,
                                                      slong nbOfRoots,
                                                      const realRat_t a,
                                                      cacheApp_t cache,
                                                      cacheCauchy_t cacheCau,
                                                      slong prec,
                                                      metadatas_t meta, int depth);

cauchyTest_res cauchyTest_probabilistic_counting( const compDsk_t Delta,           
                                                  cacheApp_t cache,
                                                  cacheCauchy_t cacheCau,
                                                  slong prec,
                                                  metadatas_t meta, int depth);

/* this version works for any iso ratio, not necessarily the one
 * defined in CacheCau 
 * the disk is not necessarily isoRatio-isolated
 * can fail, otherwise
 * returns the number of roots in Delta */
cauchyTest_res cauchyTest_probabilistic_counting_withIsoRatio( const realRat_t isoRatio,
                                                               const compDsk_t Delta,
                                                               cacheApp_t cache,
                                                               cacheCauchy_t cacheCau,
                                                               slong prec,
                                                               metadatas_t meta, int depth);

void cauchyTest_fmatheta( realRat_t res, const realRat_t a, const realRat_t theta );
void cauchyTest_fpatheta( realRat_t res, const realRat_t a, const realRat_t theta );
    
cauchyTest_res cauchyTest_deterministic_counting( const compDsk_t Delta, 
                                                  const realRat_t a,
                                                  cacheApp_t cache,
                                                  cacheCauchy_t cacheCau,
                                                  slong prec,
                                                  metadatas_t meta, int depth);

/* this version verifies that there are nbOfRoots roots in Delta */
cauchyTest_res cauchyTest_deterministic_verification( const compDsk_t Delta,
                                                      slong nbOfRoots,
                                                      const realRat_t a,
                                                      cacheApp_t cache,
                                                      cacheCauchy_t cacheCau,
                                                      slong prec,
                                                      metadatas_t meta, int depth);

/* certification that A(center, radInf, radSup) contains no root: */
/* if  0 then A(center, radInf, radSup) contains no root */
/* if -1 then A(center, (radSup+radInf)/2 - (5/4)*isoRatio*(radSup-radInf)/2,  */
/*                      (radSup+radInf)/2 + (5/4)*isoRatio*(radSup-radInf)/2 ) */
/*             contains a root */
cauchyTest_res cauchyTest_rootFreeAnnulus( const compRat_t center,
                                           const realRat_t radInf,  
                                           const realRat_t radSup,
                                           cacheApp_t cache,
                                           cacheCauchy_t cacheCau,
                                           slong prec,
                                           metadatas_t meta, int depth);

// cauchyTest_res cauchyTest_deterministic_counting_with_radInf_radSup( const compRat_t center,
//                                                                      const realRat_t radInf,  
//                                                                      const realRat_t radSup,
//                                                                      cacheApp_t cache,
//                                                                      cacheCauchy_t cacheCau,
//                                                                      slong prec,
//                                                                      metadatas_t meta, int depth);

cauchyTest_res cauchyTest_deterministic_counting_combinatorial( const compRat_t center,
                                                                     realRat_t radius,
                                                                     slong nbOfRoots,
                                                                     cacheApp_t cache,
                                                                     cacheCauchy_t cacheCau,
                                                                     slong prec,
                                                                     metadatas_t meta, int depth);

/* this version works for any iso ratio, not necessarily the one
 * defined in CacheCau 
 * the disk has to be isolated as isoRatio
 * returns the number of roots in Delta */
slong cauchyTest_computeS0compDsk( const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   metadatas_t meta, int depth);

/* this version works for any iso ratio, not necessarily the one
 * defined in CacheCau 
 * the disk IS NOT NECESSARILY isoRatio-isolated;
 * can fail */
cauchyTest_res cauchyTest_computeS0compDsk_probabilistic( const realRat_t isoRatio,
                                                          const compDsk_t Delta,
                                                          cacheApp_t cache,
                                                          cacheCauchy_t cacheCau,
                                                          metadatas_t meta, int depth);

/* Assume Delta is isoRatio-isolated and contains nbOfRoots roots */
/* Computes a disk res with radius less than eps */
/* containing s1(p, Delta) */
slong cauchyTest_computeS1compDsk( compDsk_t res,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   metadatas_t meta, int depth);

/* Assume Delta is isoRatio-isolated and contains nbOfRoots roots */
/* Assume SS is initialized to contain at least nbOfRoots + 1 numbers */
/* Computes nbOfRoots + 1 power sums with error less than eps     */
cauchyTest_res cauchyTest_computeSScompDsk( compApp_ptr SS,
                                   const realRat_t isoRatio,
                                   const compDsk_t Delta,
                                   slong nbOfRoots,
                                   cacheApp_t cache,
                                   cacheCauchy_t cacheCau,
                                   const realRat_t eps,
                                   slong prec,
                                   metadatas_t meta, int depth);
#ifdef __cplusplus
}
#endif

#endif
