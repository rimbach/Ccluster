/* ************************************************************************** */
/*  Copyright (C) 2019 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#ifndef CCLUSTER_ISSAC20_H
#define CCLUSTER_ISSAC20_H

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
#include "ccluster/ccluster.h"


#ifdef CCLUSTER_HAVE_PTHREAD
#include "ccluster/parallel_discard.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* experimental version */
/* implemented in ccluster_interface.c */
/* output: 0 for OK */
/*         2 for not enough sols*/
/*         3 for incorrect cluster*/
/*         3 for depth of subdivision tree reached */
int ccluster_issac20_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                          const realRat_t eps, 
                                          char * stratstr,
                                          int nbThreads,
                                          int output,
                                          int verb);

int ccluster_issac20_global_interface_func_eval( void(*func)(compApp_poly_t, slong), 
                                              void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                                          const realRat_t eps, 
                                          char * stratstr,
                                          int nbThreads,
                                          int output,
                                          int verb);

/* implemented in ccluster_issac20.c */
int ccluster_issac20_algo_global( connCmp_list_t qResults, 
                                const compBox_t initialBox, 
                                const realRat_t eps, 
                                cacheApp_t cache, 
                                metadatas_t meta);

int ccluster_issac20_main_loop( connCmp_list_t qResults,  
                              connCmp_list_t qMainLoop, 
                              connCmp_list_t discardedCcs, 
                              const realRat_t eps, 
                              cacheApp_t cache, 
                              metadatas_t meta);

void ccluster_issac20_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              cacheApp_t cache, 
                              metadatas_t meta, 
                              slong nbThreads); 

slong ccluster_issac20_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
                                     slong prec, metadatas_t meta);

int metadatas_issac20_fprint(FILE * file, int res, metadatas_t meta, const realRat_t eps);



#ifdef __cplusplus
}
#endif

#endif
