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

#ifndef CCLUSTER_H
#define CCLUSTER_H

#include "base/base.h"
#include "geometry/compBox.h"
#include "geometry/box_dsk.h"
#include "geometry/connCmp.h"
#include "geometry/connCmp_dsk.h"
#include "geometry/subdBox.h"
#include "geometry/connCmp_union_find.h"
#include "caches/cacheApp.h"
#include "metadatas/chronos.h"
#include "metadatas/metadatas.h"
#include "lists/compBox_list.h"
#include "lists/connCmp_list.h"
#include "tstar/tstar.h"
#include "newton/newton.h"
#include "powerSums/powerSums.h"

#include "geometry/compAnn.h"
#include "rootRadii/realIntRootRadii.h"

#ifdef CCLUSTER_HAVE_PTHREAD
#include "ccluster/parallel_discard.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

slong ccluster_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
                                     cacheApp_t cache, 
//                                      int nbSols, 
                                     slong prec, metadatas_t meta);

void ccluster_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              compBox_list_t bDiscarded,
                              cacheApp_t cache, 
                              metadatas_t meta, 
                              slong nbThreads);  

void ccluster_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t qPrepLoop, 
                         connCmp_list_t discardedCcs, 
                         cacheApp_t cache, 
                         metadatas_t meta);

int  ccluster_compDsk_is_separated( const compDsk_t d, 
                                    connCmp_list_t qMainLoop, 
                                    connCmp_list_t discardedCcs );

void ccluster_main_loop( connCmp_list_t qResults, 
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta);

void ccluster_algo( connCmp_list_t qResults, 
                    compBox_list_t bDiscarded, 
                    const compBox_t initialBox, 
                    const realRat_t eps, 
                    cacheApp_t cache, 
                    metadatas_t meta);

void ccluster_algo_global( connCmp_list_t qResults,
                           compBox_list_t bDiscarded,
                           const compBox_t initialBox, 
                           const realRat_t eps, 
                           cacheApp_t cache, 
                           metadatas_t meta);

void ccluster_refine( connCmp_list_t qResults, 
                      connCmp_list_t qMainLoop,
                      const realRat_t eps, 
                      cacheApp_t cache, 
                      metadatas_t meta);

/* rootRadii: implemented in ccluster_rootRadii.c */

/* returns 1 => 2*box contains at least one root */
/*         0 => box contains no root */
/*        -1 => can not decide */
int ccluster_rootRadii_exclusion_test( compBox_t box, slong prec, metadatas_t meta );

void ccluster_algo_global_rootRadii( connCmp_list_t qResults,
                                     compBox_list_t bDiscarded,
                                     compAnn_list_t annulii,
                                     compAnn_list_t annulii1,
                                     compAnn_list_t annulii2,
                                     const compBox_t initialBox, 
                                     const realRat_t eps, 
                                     cacheApp_t cache, 
                                     metadatas_t meta);

void connCmp_print_for_results(FILE * f, 
                               const connCmp_t c, 
                               metadatas_t meta);

void connCmp_list_print_for_results(FILE * f, 
                                    const connCmp_list_t c, 
                                    metadatas_t meta);

void connCmp_print_for_results_withOutput(FILE * f, 
                                          const connCmp_t c, 
                                          int output, 
                                          metadatas_t meta);

void connCmp_list_print_for_results_withOutput(FILE * f, 
                                               const connCmp_list_t l, 
                                               int output, 
                                               metadatas_t meta);

void connCmp_gnuplot(FILE * f, 
                     const connCmp_t c, 
                     metadatas_t meta);

void compBox_gnuplot(FILE * f, 
                     const compBox_t b);

void connCmp_list_gnuplot(FILE * f, 
                          const connCmp_list_t c, 
                          metadatas_t meta,
                          int withInitBox);

void connCmp_list_gnuplot_drawSubdiv(FILE * f, 
                          const connCmp_list_t l, 
                          const compBox_list_t lb,
                          metadatas_t meta);

void connCmp_list_gnuplot_drawSubdiv_rootRadii(FILE * f, 
                          const connCmp_list_t l, 
                          const compBox_list_t lb,
                          const compAnn_list_t la,
                          const compAnn_list_t la1,
                          const compAnn_list_t la2,
                          metadatas_t meta);

/* INTERFACES */

/* default interfaces */

void ccluster_interface_func( void(*func)(compApp_poly_t, slong), 
                              const compBox_t initialBox, 
                              const realRat_t eps, 
                              char * stratstr,
                              int nbThreads,
                              int output,
                              int verb);


void ccluster_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int output,
                                     int verb);

/* name ccluster_global_interface_poly is for singular */
void ccluster_interface_realRat_poly( const realRat_poly_t poly,
                                      const compBox_t initialBox,
                                      const realRat_t eps, 
                                      char * stratstr,
                                      int nbThreads,
                                      int output,
                                      int verb);

void ccluster_global_interface_realRat_poly( const realRat_poly_t poly,
                                             const realRat_t eps, 
                                             char * stratstr,
                                             int nbThreads,
                                             int output,
                                             int verb);

void ccluster_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                                const compBox_t initialBox, 
                                const realRat_t eps, 
//                                 int st,
                                char * stratstr,
                                int nbThreads,
                                int verb);

void ccluster_global_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                   void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong), 
                                   const realRat_t eps, 
                                   char * stratstr,
                                   int nbThreads,
                                   int verb);
/* interfaces for Julia */

void ccluster_forJulia_func( connCmp_list_t qResults, 
                             void(*func)(compApp_poly_t, slong), 
                             const compBox_t initialBox, 
                             const realRat_t eps, 
                             char * stratstr,
                             int nbThreads,
                             int verb);

void ccluster_global_forJulia_func( connCmp_list_t qResults, 
                                    void(*func)(compApp_poly_t, slong),  
                                    compBox_t initialBox,
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb);

void ccluster_forJulia_realRat_poly( connCmp_list_t qResults, 
                                    const realRat_poly_t poly, 
                                    const compBox_t initialBox, 
                                    const realRat_t eps, 
                                    char * stratstr,
                                    int nbThreads,
                                    int verb);

void ccluster_global_forJulia_realRat_poly( connCmp_list_t qResults, 
                                            const realRat_poly_t poly,  
                                            compBox_t initialBox,
                                            const realRat_t eps, 
                                            char * stratstr,
                                            int nbThreads,
                                            int verb);

void ccluster_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                               const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                               const compBox_t initialBox, 
                                               const realRat_t eps, 
                                               char * stratstr,
                                               int nbThreads,
                                               int verb);

void ccluster_global_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                                      const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                                      compBox_t initialBox,
                                                      const realRat_t eps, 
                                                      char * stratstr,
                                                      int nbThreads,
                                                      int verb);

void ccluster_forJulia_refine( connCmp_list_t qResults,
                               connCmp_list_t qMainLoop,
                               void(*func)(compApp_poly_t, slong), 
                               const compBox_t initialBox, 
                               const realRat_t eps, 
                               char * stratstr,
                               int nbThreads,
                               int verb);

/* res = 1: OK */
/* res = -1: lcf vanishes; may miss solutions with huge norm */
int ccluster_global_forJulia_forTcluster_func( connCmp_list_t qResults, 
                                                void(*func)(compApp_poly_t, slong), 
                                                compBox_t initialBox,
                                                const realRat_t eps, 
                                                char * stratstr,
                                                int nbThreads,
                                                int verb);

void ccluster_forJulia_draw( connCmp_list_t qResults, 
                            compBox_list_t qDiscarded, 
                            void(*func)(compApp_poly_t, slong), 
                            const compBox_t initialBox, 
                            const realRat_t eps, 
                            char * stratstr,
                            int nbThreads,
                            int verb);

/* interfaces for Singular */

int ccluster_interface_poly( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                             const compRat_poly_t poly, 
                             const compBox_t initialBox, 
                             const realRat_t eps, 
                             int st, 
                             int verb);

int ccluster_interface_poly_real( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                  const realRat_poly_t poly, 
                                  const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                  const realRat_t eps, 
                                  int st, 
                                  int verb);

int ccluster_interface_poly_real_imag( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                       const realRat_poly_t poly_real, const realRat_poly_t poly_imag,
                                       const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                       const realRat_t eps, 
                                       int st, 
                                       int verb);

///implemented in ccluster_draw.c
/* DEPRECATED */

void ccluster_interface_forJulia_draw( connCmp_list_t qResults, 
                                       compBox_list_t qDiscarded, 
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb);

void ccluster_algo_draw( connCmp_list_t qResults, 
                         compBox_list_t discarded, 
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         cacheApp_t cache, metadatas_t meta);

/* old interfaces for Julia */
/* DEPRECATED */
void ccluster_interface_forJulia( connCmp_list_t qResults, 
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb);

void ccluster_interface_forJulia_realRat_poly( connCmp_list_t qResults, 
                                              const realRat_poly_t poly, 
                                              const compBox_t initialBox, 
                                              const realRat_t eps, 
                                              int st, 
                                              int verb);

void ccluster_interface_forJulia_realRat_poly_real_imag( connCmp_list_t qResults, 
                                                         const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                                         const compBox_t initialBox, 
                                                         const realRat_t eps, 
                                                         int st, 
                                                         int verb);

void ccluster_interface_forJulia_compRat_poly( connCmp_list_t qResults, 
                                              const compRat_poly_t poly,
                                              const compBox_t initialBox, 
                                              const realRat_t eps, 
                                              int st, 
                                              int verb);

void ccluster_refine_forJulia( connCmp_list_t qResults,
                               connCmp_list_t qMainLoop,
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb);

/* utilities */
int  ccluster_compDsk_is_separated( const compDsk_t d, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs );

#ifdef __cplusplus
}
#endif

#endif
