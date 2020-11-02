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

#ifndef DEFLATE_H
#define DEFLATE_H

#ifdef DEFLATE_INLINE_C
#define DEFLATE_INLINE
#else
#define DEFLATE_INLINE static __inline__
#endif

#include "metadatas/metadatas.h"
#include "geometry/connCmp.h"
#include "geometry/compDsk.h"
#include "caches/cacheApp.h"

#include "tstar/tstar.h"

#ifdef __cplusplus
extern "C" {
#endif
    
/* deflate */
#define connCmp_isDefref(X) ( (X)->isDef)
#define connCmp_degDeref(X) ( (X)->degDe)
#define connCmp_isDFGref(X) ( (X)->isDFG)
#define connCmp_defPoref(X) (&(X)->defPo)
#define connCmp_defFGref(X) (&(X)->defFG)
    
/* memory managment */
void deflate_connCmp_init  (connCmp_t x);
void deflate_connCmp_clear (connCmp_t x);
/* setting */
void deflate_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta );
void deflate_copy( connCmp_t dest, const connCmp_t src );

tstar_res deflate_tstar_test( connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, int discard, slong prec, metadatas_t meta);

/*DEPRECATED*/
// /* deflate2 */
// #define connCmp_defClref(X) (&(X)->defCl)
// #define connCmp_defiCref(X) ( (X)->defiC)
// #define connCmp_defiDref(X) ( (X)->defiD)
// #define connCmp_defDsref(X) (&(X)->defDs)
// #define connCmp_defPfref(X) (&(X)->defPf)
// #define connCmp_defPrref(X) ( (X)->defPr)
// #define connCmp_defIRref(X) (&(X)->defIR)
//     
// /* memory managment */
// 
// void deflate2_connCmp_init  (connCmp_t x);
// void deflate2_connCmp_clear (connCmp_t x);
// 
// /* setting */
// void deflate2_set_clearance (connCmp_t x, const compDsk_t disk);
// void deflate2_copy_clearance (connCmp_t dest, const connCmp_t src);
// int deflate2_isDef_clearance (const connCmp_t x);
// 
// void deflate2_set (connCmp_t x, const compDsk_t disk, const fmpz_t isoRatio);
// void deflate2_copy (connCmp_t dest, const connCmp_t src);
// 
// slong deflate2_get_Nb_EvalPoints( const fmpz_t isoRatio, slong degree, slong nbPs, slong prec, metadatas_t meta );
// void deflate2_compute_factor( connCmp_t x, cacheApp_t cache, slong prec, metadatas_t meta );
// tstar_res deflate2_tstar_test( const connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, slong prec, metadatas_t meta);
// 
// /* power sums */
// void deflate_evaluateAtPoints( compApp_ptr f_val,
//                                  compApp_ptr fder_val,
//                                  const compApp_ptr points,
//                                  slong nbPoints,
//                                  cacheApp_t cache,
//                                  slong prec,
//                                  metadatas_t meta);
// 
// int deflate_computePsApprox_fromVals(compApp_ptr ps,
//                                         const realRat_t radius,
//                                         compApp_ptr points,
//                                         compApp_ptr fvals,
//                                         compApp_ptr fdervals,
//                                         compApp_ptr fdivs,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
//                                         slong prec,
//                                         metadatas_t meta);
// 
// slong deflate_computePsApprox(compApp_ptr ps,
//                                         const compRat_t center,
//                                         const realRat_t radius,
//                                         compApp_ptr points,
//                                         compApp_ptr pointsShifted,
//                                         compApp_ptr fvals,
//                                         compApp_ptr fdervals,
//                                         compApp_ptr fdivs,
//                                         slong nbPoints,
//                                         slong nbPowerSums,
//                                         cacheApp_t cache,
//                                         slong prec,
//                                         metadatas_t meta);

#ifdef __cplusplus
}
#endif

#endif
