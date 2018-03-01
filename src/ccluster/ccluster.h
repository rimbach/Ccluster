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

#ifdef CCLUSTER_HAVE_PTHREAD
#include <pthread.h>

// typedef struct {
//     int nbsol;
//     slong prec;
//     compBox_ptr box;
//     cacheApp_ptr cache;
//     metadatas_ptr meta;
// } parallel_discard_arg_t;

typedef struct {
    slong prec;
    compBox_list_ptr boxes;
    cacheApp_ptr cache;
    metadatas_ptr meta;
    int status; /* 0: default, 1: is_running, 2: is_finnished */
    pthread_mutex_t mutex;
} parallel_discard_list_arg_t;

int nb_thread_running;
pthread_mutex_t mutex_nb_running;

typedef struct {
    connCmp_list_ptr res;
    connCmp_ptr      cc;
    connCmp_list_ptr dis;
    cacheApp_ptr cache;
    metadatas_ptr meta;
    slong nbThreads;
} parallel_bisect_arg_t;

// void * _parallel_discard_worker( void * arg_ptr );

void * _parallel_discard_list_worker( void * arg_ptr );

slong ccluster_parallel_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
                                        slong prec, metadatas_t meta, slong nbThreads);

void * _parallel_bisect_worker( void * arg_ptr );
void ccluster_parallel_bisect_connCmp_list( connCmp_list_ptr qMainLoop, connCmp_list_ptr discardedCcs,
                                            connCmp_list_ptr toBeBisected, cacheApp_t cache, metadatas_t meta);
#endif

slong ccluster_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
//                                      int nbSols, 
                                     slong prec, metadatas_t meta);

void ccluster_bisect_connCmp( connCmp_list_t dest, connCmp_t cc, connCmp_list_t discardedCcs, cacheApp_t cache, metadatas_t meta, slong nbThreads);  

void ccluster_prep_loop( connCmp_list_t qMainLoop, connCmp_list_t qPrepLoop, connCmp_list_t discardedCcs, cacheApp_t cache, metadatas_t meta);

int  ccluster_compDsk_is_separated( const compDsk_t d, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs );

void ccluster_main_loop( connCmp_list_t qResults,  connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, const realRat_t eps, cacheApp_t cache, metadatas_t meta);

void ccluster_algo( connCmp_list_t qResults, const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta);

void connCmp_print_for_results(FILE * f, const connCmp_t c, metadatas_t meta);
void connCmp_list_print_for_results(FILE * f, const connCmp_list_t c, metadatas_t meta);

void ccluster_interface_func( void(*func)(compApp_poly_t, slong), const compBox_t initialBox, const realRat_t eps, int st, int verb);

void ccluster_interface_forJulia( connCmp_list_t qResults, 
                                  void(*func)(compApp_poly_t, slong), 
                                  const compBox_t initialBox, 
                                  const realRat_t eps, 
                                  int st, 
                                  int verb);

#endif
