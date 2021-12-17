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

#ifndef CAUCHY_H
#define CAUCHY_H

#include "base/base.h"
#include "geometry/compBox.h"
#include "geometry/box_dsk.h"
#include "geometry/connCmp.h"
#include "geometry/connCmp_dsk.h"
#include "geometry/subdBox.h"
#include "geometry/connCmp_union_find.h"
#include "caches/cacheApp.h"
#include "caches/cacheCauchy.h"
#include "metadatas/chronos.h"
#include "metadatas/metadatas.h"
#include "lists/compBox_list.h"
#include "lists/connCmp_list.h"
#include "tstar/tstar.h"
#include "newton/newton.h"
#include "powerSums/powerSums.h"
#include "ccluster/ccluster.h"

#include "cauchy_rootRadii/cauchy_rootRadii.h"
#include "turan_rootRadii/turan_rootRadii.h"
#include "cauchy_tests/cauchy_tests.h"

#include "fpri/fpri.h"
#include "fpri/fpci.h"
#include "fpri/fpri_poly.h"
#include "fpri/fpci_poly.h"

#ifdef CCLUSTER_HAVE_PTHREAD
#include "ccluster/parallel_discard.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

int cauchy_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                  const realRat_t eps,
                                  const realRat_t isoRatio,
                                  int nbPows,
                                  int certified,
                                  int usefpri,
                                  char * stratstr,
                                  int nbThreads,
                                  int output,
                                  int verb);

int cauchy_global_interface_realRat_poly( const realRat_poly_t poly,
                                          const realRat_t eps,
                                          const realRat_t isoRatio,
                                          int nbPows,
                                          int certified,
                                          int usefpri,
                                          char * stratstr,
                                          int nbThreads,
                                          int output,
                                          int verb);

/* version with function for fast evaluation */
int cauchy_global_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                       void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong), 
                                       const realRat_t eps, 
                                       const realRat_t isoRatio,
                                       int nbPows,
                                       int certified,
                                       int usefpri,
                                       char * stratstr,
                                       int nbThreads,
                                       int output,
                                       int verb);

/* implemented in cauchy.c */
int cauchy_algo_global( connCmp_list_t qResults,
                           compBox_list_t bDiscarded,
                           const realRat_t eps,
                           cacheCauchy_t cacheCau,
                        int certified,
                           metadatas_t meta);

int cauchy_main_loop( connCmp_list_t qResults, 
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheCauchy_t cacheCau,
                      int certified,
                         metadatas_t meta);

void cauchy_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t qPrepLoop, 
                         connCmp_list_t discardedCcs, 
                         cacheCauchy_t cacheCau,
                         metadatas_t meta);

void cauchy_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              compBox_list_t bDiscarded,
                              cacheCauchy_t cacheCau,
                              slong nbMSols, 
                              metadatas_t meta, 
                              slong nbThreads); 

slong cauchy_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
                                     cacheCauchy_t cacheCau,
                                     slong nbMSols, 
                                     slong prec, metadatas_t meta);

/* assume Delta = D(c,r) contains m and has isolation ratio theta >=2 */
/* computes a disk res = D(c',r') such that*/
/* Delta and res contain the same roots */
/* either r' <= eps */
/*     or res is m/((2m-2)*theta) rigid */
slong cauchy_compressionIntoRigidDisk( compDsk_t res, const compDsk_t Delta, slong m, const realRat_t theta, const realRat_t eps,
                                       cacheCauchy_t cacheCau,
                                       slong prec, metadatas_t meta, slong depth);

connCmp_ptr cauchy_actualizeCCafterCompression( connCmp_ptr CC, const compDsk_t Delta, slong appPrec, metadatas_t meta );



int metadatas_cauchy_fprint(FILE * file, metadatas_t meta, const realRat_t eps, cacheCauchy_t cacheCau, int certified);



#ifdef __cplusplus
}
#endif

#endif
