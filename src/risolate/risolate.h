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

#ifndef RISOLATE_H
#define RISOLATE_H

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
#include "deflation/deflate.h"
// #include "powerSums/powerSums.h"
#include "geometry/compAnn.h"
#include "rootRadii/realIntRootRadii.h"

#ifdef CCLUSTER_HAVE_PTHREAD
#include "ccluster/parallel_discard.h"
#endif

#include "ccluster/ccluster.h"

#ifdef __cplusplus
extern "C" {
#endif

void risolate_compBox_get_containing_dsk( compDsk_t d, const compBox_t b);

slong risolate_discard_compBox_list( compBox_list_t boxes,
                                     compBox_list_t bDiscarded,
                                     connCmp_t cc,
                                     cacheApp_t cache, 
                                     slong prec, 
                                     metadatas_t meta);
  
void risolate_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs,
                              compBox_list_t bDiscarded,
                              cacheApp_t cache, 
                              metadatas_t meta, 
                              slong nbThreads);

void risolate_algo( connCmp_list_t qResults, 
                    compBox_list_t bDiscarded,
		            const compBox_t initialBox, 
		            const realRat_t eps, 
		            cacheApp_t cache, 
		            metadatas_t meta);

void risolate_algo_global( connCmp_list_t qResults, 
                           compBox_list_t bDiscarded,
			               const compBox_t initialBox, 
			               const realRat_t eps, 
			               cacheApp_t cache, 
			               metadatas_t meta);

void risolate_main_loop( connCmp_list_t qResults,
                         compBox_list_t bDiscarded,
			             connCmp_list_t qMainLoop, 
			             connCmp_list_t discardedCcs, 
			             const realRat_t eps, 
			             cacheApp_t cache, 
			             metadatas_t meta);

void risolate_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
			             connCmp_list_t qPrepLoop, 
			             connCmp_list_t discardedCcs, 
			             cacheApp_t cache, 
			             metadatas_t meta);

slong risolate_discard_compBox_list_rootRadii( compBox_list_t boxes, 
                                                       compBox_list_t bDiscarded,
                                                       cacheApp_t cache, 
                                                       slong prec, 
                                                       metadatas_t meta);

void risolate_bisect_connCmp_rootRadii( connCmp_list_t dest, 
                                                 connCmp_t cc, 
                                                 connCmp_list_t discardedCcs,
                                                 compBox_list_t bDiscarded, 
                                                 cacheApp_t cache, 
                                                 metadatas_t meta, 
                                                 slong nbThreads);

void risolate_main_loop_rootRadii( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta);

void risolate_algo_global_rootRadii  ( connCmp_list_t qResults, 
                                       compBox_list_t bDiscarded,
                                       compAnn_list_t annulii,
                                       const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta);

/* INTERFACES */

/* default interfaces */

void risolate_global_interface_poly( const realRat_poly_t poly, 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int output,
                                     int verb);

void risolate_interface_poly( const realRat_poly_t poly,
                              const compBox_t initialBox, 
                              const realRat_t eps, 
                              char * stratstr,
                              int nbThreads,
                              int output,
                              int verb);

void connCmp_risolate_print_for_results(FILE * f, const connCmp_t c, metadatas_t meta);

void connCmp_list_risolate_print_for_results(FILE * f, const connCmp_list_t l, metadatas_t meta);

void connCmp_risolate_print_for_results_withOutput(FILE * f, const connCmp_t c, int output, metadatas_t meta);

void connCmp_list_risolate_print_for_results_withOutput(FILE * f, const connCmp_list_t l, int output, metadatas_t meta);

void risolate_connCmp_gnuplot(FILE * f, 
                     const connCmp_t c, 
                     metadatas_t meta);

void risolate_compBox_gnuplot(FILE * f, 
                     const compBox_t b);

void risolate_connCmp_list_gnuplot(FILE * f, 
                          const connCmp_list_t c, 
                          metadatas_t meta,
                          int withInitBox);

void risolate_connCmp_list_gnuplot_drawSubdiv(FILE * f, 
                          const connCmp_list_t l, 
                          const compBox_list_t lb,
                          metadatas_t meta);

void risolate_connCmp_list_gnuplot_drawSubdiv_rootRadii(FILE * f, 
                          const connCmp_list_t l, 
                          const compBox_list_t lb,
                          const compAnn_list_t la,
                          const compAnn_list_t la1,
                          const compAnn_list_t la2,
                          metadatas_t meta);

/* DEPRECATED */

// int risolate_compBox_intersects_only_one( const compBox_t b, int nbList );
// int risolate_compBox_intersects_atLest_one( const compBox_t b, int nbList );

// void risolate_prep_loop_rootRadii( compBox_list_t bDiscarded, 
//                             connCmp_list_t qResult, 
//                            connCmp_list_t qPrepLoop, 
//                            connCmp_list_t discardedCcs, 
//                            cacheApp_t cache, 
//                            metadatas_t meta);
/*
slong risolate_exclusion_rootRadii( connCmp_list_t qCover,
                                   cacheApp_t cache, 
                                   metadatas_t meta);

void risolate_algo_global_rootRadii_old( connCmp_list_t qResults,
                                     compBox_list_t bDiscarded,
                                     compAnn_list_t annulii,
                                     compAnn_list_t annulii1,
                                     compAnn_list_t annulii2,
                                     const compBox_t initialBox, 
                                     const realRat_t eps, 
                                     cacheApp_t cache, 
                                     metadatas_t meta);*/

#ifdef __cplusplus
}
#endif

#endif
