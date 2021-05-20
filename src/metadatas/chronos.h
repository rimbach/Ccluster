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

#ifndef CHRONOS_H
#define CHRONOS_H

#ifdef METADATAS_INLINE_C
#define METADATAS_INLINE
#else
#define METADATAS_INLINE static __inline__
#endif

#include <time.h>
#include "base/base.h"
#ifdef CCLUSTER_HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    /* DEPRECATED ?*/
//     clock_t _clicks_Approxi;
//     clock_t _clicks_Graeffe;
//     clock_t _clicks_Taylors;
//     clock_t _clicks_T0Tests;
//     clock_t _clicks_TSTests;
//     clock_t _clicks_Newtons;
//     clock_t _clicks_CclusAl;
//     clock_t _clicks_Derivat;
//     clock_t _clicks_Anticip;
//     clock_t _clicks_Evaluat;
    /* END DEPRECATED ?*/
    double  _clicks_Approxi_cumul;
    double  _clicks_Graeffe_cumul;
    double  _clicks_Taylors_cumul;
    double  _clicks_T0Tests_cumul;
    double  _clicks_TSTests_cumul;
    double  _clicks_Newtons_cumul;
    double  _clicks_CclusAl_cumul;
    double  _clicks_Derivat_cumul;
    double  _clicks_Anticip_cumul;
    /* for powerSums */
    double  _clicks_PSTests_cumul;
    double  _clicks_Evaluat_cumul;
    /* for rootRadii */
    double  _clicks_RRTaylo_cumul;
    double  _clicks_RRGraef_cumul;
    double  _clicks_rootRad_cumul;
    double  _clicks_RRT0Tes_cumul;
    
    double  _clicks_NeTSTes_cumul;
    /* deflation */
    double  _clicks_DefTayl_cumul;
    double  _clicks_DefDeri_cumul;
    double  _clicks_DefEval_cumul;
    double  _clicks_DefScal_cumul;
    double  _clicks_DefGrae_cumul;
    double  _clicks_DefTsta_cumul;
    
    /* Cauchy root finder */
    double  _clicks_CauExTo_cumul; /* time in Cauchy exclusion tests */
    double  _clicks_CauExEP_cumul; /* time in evaluations in probabilistic exclusion tests */
    double  _clicks_CauExED_cumul; /* time in evaluations in certifying    exclusion tests */
    double  _clicks_CauExDS_cumul; /* time in computing divisions in       exclusion tests */
    double  _clicks_CauExCS_cumul; /* time in computing s0 in              exclusion tests */
    
    double  _clicks_CauCoTo_cumul; /* time in Cauchy counting tests */
    double  _clicks_CauCoEP_cumul; /* time in evaluations in probabilistic counting tests */
    double  _clicks_CauCoED_cumul; /* time in evaluations in certifying    counting tests */
    double  _clicks_CauCoDS_cumul; /* time in computing divisions in       counting tests */
    double  _clicks_CauCoCS_cumul; /* time in computing s0 in              counting tests */
    
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_t _mutex;
#endif    
} chronos;

typedef chronos chronos_t[1];
typedef chronos * chronos_ptr;
    
void chronos_init( chronos_t times );
void chronos_clear( chronos_t times );

// void chronos_join ( chronos_t t1, const chronos_t t2 );
METADATAS_INLINE void chronos_lock(chronos_t t){
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_lock (&(t->_mutex));
#endif
}

METADATAS_INLINE void chronos_unlock(chronos_t t){
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_unlock (&(t->_mutex));
#endif
}

METADATAS_INLINE void   chronos_add_time_Approxi( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Approxi_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Approxi_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Approxi_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_Taylors( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Taylors_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Taylors_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Taylors_cumul += d/CLOCKS_PER_SEC;
// #endif
}
METADATAS_INLINE void   chronos_add_time_Graeffe( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Graeffe_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Graeffe_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Graeffe_cumul += d/CLOCKS_PER_SEC;
// #endif
}
METADATAS_INLINE void   chronos_add_time_T0Tests( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_T0Tests_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_T0Tests_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_T0Tests_cumul += d/CLOCKS_PER_SEC;
// #endif
}
METADATAS_INLINE void   chronos_add_time_TSTests( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_TSTests_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_TSTests_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_TSTests_cumul += d/CLOCKS_PER_SEC;
// #endif
}
METADATAS_INLINE void   chronos_add_time_Anticip( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Anticip_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Anticip_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Anticip_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_Newtons( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Newtons_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Newtons_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Newtons_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_PSTests( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_PSTests_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_PSTests_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_PSTests_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_Evaluat( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_Evaluat_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_Evaluat_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_Evaluat_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_CclusAl( chronos_t times, double d, int nbThreads ){
// #ifdef CCLUSTER_HAVE_PTHREAD
//     if (nbThreads>1)
//         times->_clicks_CclusAl_cumul += d/(nbThreads*CLOCKS_PER_SEC);
//     else 
//         times->_clicks_CclusAl_cumul += d/CLOCKS_PER_SEC;
// #else
    times->_clicks_CclusAl_cumul += d/CLOCKS_PER_SEC;
// #endif
}

METADATAS_INLINE void   chronos_add_time_RRTaylo( chronos_t times, double d, int nbThreads ){
    times->_clicks_RRTaylo_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE void   chronos_add_time_RRGraef( chronos_t times, double d, int nbThreads ){
    times->_clicks_RRGraef_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE void   chronos_add_time_rootRad( chronos_t times, double d, int nbThreads ){
    times->_clicks_rootRad_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE void   chronos_add_time_RRT0Tes( chronos_t times, double d, int nbThreads ){
    times->_clicks_RRT0Tes_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE double chronos_get_time_Approxi ( const chronos_t times ) { return times->_clicks_Approxi_cumul; }
METADATAS_INLINE double chronos_get_time_Graeffe ( const chronos_t times ) { return times->_clicks_Graeffe_cumul; }
METADATAS_INLINE double chronos_get_time_Taylors ( const chronos_t times ) { return times->_clicks_Taylors_cumul; }
METADATAS_INLINE double chronos_get_time_T0Tests ( const chronos_t times ) { return times->_clicks_T0Tests_cumul; }
METADATAS_INLINE double chronos_get_time_TSTests ( const chronos_t times ) { return times->_clicks_TSTests_cumul; }
METADATAS_INLINE double chronos_get_time_Newtons ( const chronos_t times ) { return times->_clicks_Newtons_cumul; }
METADATAS_INLINE double chronos_get_time_CclusAl ( const chronos_t times ) { return times->_clicks_CclusAl_cumul; }
METADATAS_INLINE double chronos_get_time_Derivat ( const chronos_t times ) { return times->_clicks_Derivat_cumul; }
METADATAS_INLINE double chronos_get_time_Anticip ( const chronos_t times ) { return times->_clicks_Anticip_cumul; }

METADATAS_INLINE double chronos_get_time_PSTests ( const chronos_t times ) { return times->_clicks_PSTests_cumul; }
METADATAS_INLINE double chronos_get_time_Evaluat ( const chronos_t times ) { return times->_clicks_Evaluat_cumul; }

METADATAS_INLINE void   chronos_add_time_NeTSTes( chronos_t times, double d, int nbThreads ){
    times->_clicks_NeTSTes_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE void   chronos_add_time_DefTayl( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefTayl_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE void   chronos_add_time_DefDeri( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefDeri_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE void   chronos_add_time_DefEval( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefEval_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE void   chronos_add_time_DefScal( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefScal_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE void   chronos_add_time_DefGrae( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefGrae_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE void   chronos_add_time_DefTsta( chronos_t times, double d, int nbThreads ){
    times->_clicks_DefTsta_cumul += d/CLOCKS_PER_SEC;
}

METADATAS_INLINE double chronos_get_time_NeTSTes ( const chronos_t times ) { return times->_clicks_NeTSTes_cumul; }
METADATAS_INLINE double chronos_get_time_DefTayl ( const chronos_t times ) { return times->_clicks_DefTayl_cumul; }
METADATAS_INLINE double chronos_get_time_DefDeri ( const chronos_t times ) { return times->_clicks_DefDeri_cumul; }
METADATAS_INLINE double chronos_get_time_DefEval ( const chronos_t times ) { return times->_clicks_DefEval_cumul; }
METADATAS_INLINE double chronos_get_time_DefScal ( const chronos_t times ) { return times->_clicks_DefScal_cumul; }
METADATAS_INLINE double chronos_get_time_DefGrae ( const chronos_t times ) { return times->_clicks_DefGrae_cumul; }
METADATAS_INLINE double chronos_get_time_DefTsta ( const chronos_t times ) { return times->_clicks_DefTsta_cumul; }

METADATAS_INLINE double chronos_get_time_RRTaylo ( const chronos_t times ) { return times->_clicks_RRTaylo_cumul; }
METADATAS_INLINE double chronos_get_time_RRGraef ( const chronos_t times ) { return times->_clicks_RRGraef_cumul; }
METADATAS_INLINE double chronos_get_time_rootRad ( const chronos_t times ) { return times->_clicks_rootRad_cumul; }
METADATAS_INLINE double chronos_get_time_RRT0Tes ( const chronos_t times ) { return times->_clicks_RRT0Tes_cumul; }

/* Cauchy root finder */
METADATAS_INLINE void   chronos_add_time_CauExTo( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauExTo_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauExTo ( const chronos_t times ) { return times->_clicks_CauExTo_cumul; }

METADATAS_INLINE void   chronos_add_time_CauExEP( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauExEP_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauExEP ( const chronos_t times ) { return times->_clicks_CauExEP_cumul; }

METADATAS_INLINE void   chronos_add_time_CauExED( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauExED_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauExED ( const chronos_t times ) { return times->_clicks_CauExED_cumul; }

METADATAS_INLINE void   chronos_add_time_CauExDS( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauExDS_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauExDS ( const chronos_t times ) { return times->_clicks_CauExDS_cumul; }

METADATAS_INLINE void   chronos_add_time_CauExCS( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauExCS_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauExCS ( const chronos_t times ) { return times->_clicks_CauExCS_cumul; }

METADATAS_INLINE void   chronos_add_time_CauCoTo( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauCoTo_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauCoTo ( const chronos_t times ) { return times->_clicks_CauCoTo_cumul; }

METADATAS_INLINE void   chronos_add_time_CauCoEP( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauCoEP_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauCoEP ( const chronos_t times ) { return times->_clicks_CauCoEP_cumul; }

METADATAS_INLINE void   chronos_add_time_CauCoED( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauCoED_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauCoED ( const chronos_t times ) { return times->_clicks_CauCoED_cumul; }

METADATAS_INLINE void   chronos_add_time_CauCoDS( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauCoDS_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauCoDS ( const chronos_t times ) { return times->_clicks_CauCoDS_cumul; }

METADATAS_INLINE void   chronos_add_time_CauCoCS( chronos_t times, double d, int nbThreads ){
    times->_clicks_CauCoCS_cumul += d/CLOCKS_PER_SEC;
}
METADATAS_INLINE double chronos_get_time_CauCoCS ( const chronos_t times ) { return times->_clicks_CauCoCS_cumul; }

/* DEPRECATED */
// METADATAS_INLINE void   chronos_tic_Approxi      ( chronos_t times ) { times->_clicks_Approxi = clock(); }
// METADATAS_INLINE void   chronos_toc_Approxi      ( chronos_t times ) { times->_clicks_Approxi_cumul += ((double) (clock() - times->_clicks_Approxi))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Graeffe      ( chronos_t times ) { times->_clicks_Graeffe = clock(); }
// METADATAS_INLINE void   chronos_toc_Graeffe      ( chronos_t times ) { times->_clicks_Graeffe_cumul += ((double) (clock() - times->_clicks_Graeffe))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Taylors      ( chronos_t times ) { times->_clicks_Taylors = clock(); }
// METADATAS_INLINE void   chronos_toc_Taylors      ( chronos_t times ) { times->_clicks_Taylors_cumul += ((double) (clock() - times->_clicks_Taylors))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_T0Tests      ( chronos_t times ) { times->_clicks_T0Tests = clock(); }
// METADATAS_INLINE void   chronos_toc_T0Tests      ( chronos_t times ) { times->_clicks_T0Tests_cumul += ((double) (clock() - times->_clicks_T0Tests))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_TSTests      ( chronos_t times ) { times->_clicks_TSTests = clock(); }
// METADATAS_INLINE void   chronos_toc_TSTests      ( chronos_t times ) { times->_clicks_TSTests_cumul += ((double) (clock() - times->_clicks_TSTests))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Newtons      ( chronos_t times ) { times->_clicks_Newtons = clock(); }
// METADATAS_INLINE void   chronos_toc_Newtons      ( chronos_t times ) { times->_clicks_Newtons_cumul += ((double) (clock() - times->_clicks_Newtons))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_CclusAl      ( chronos_t times ) { times->_clicks_CclusAl = clock(); }
// METADATAS_INLINE void   chronos_toc_CclusAl      ( chronos_t times ) { times->_clicks_CclusAl_cumul += ((double) (clock() - times->_clicks_CclusAl))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Evaluat      ( chronos_t times ) { times->_clicks_Evaluat = clock(); }
// METADATAS_INLINE void   chronos_toc_Evaluat      ( chronos_t times ) { times->_clicks_Evaluat_cumul += ((double) (clock() - times->_clicks_Evaluat))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Derivat      ( chronos_t times ) { times->_clicks_Derivat = clock(); }
// METADATAS_INLINE void   chronos_toc_Derivat      ( chronos_t times ) { times->_clicks_Derivat_cumul += ((double) (clock() - times->_clicks_Derivat))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Anticip      ( chronos_t times ) { times->_clicks_Anticip = clock(); }
// METADATAS_INLINE void   chronos_toc_Anticip      ( chronos_t times ) { times->_clicks_Anticip_cumul += ((double) (clock() - times->_clicks_Anticip))/ CLOCKS_PER_SEC; }
// 
// 
// METADATAS_INLINE void   chronos_tic_Test      ( chronos_t times, int discard ) { 
//     if (discard) times->_clicks_T0Tests = clock(); 
//     else         times->_clicks_TSTests = clock();
// }
// 
// METADATAS_INLINE void   chronos_toc_Test      ( chronos_t times, int discard ) { 
//     if (discard) times->_clicks_T0Tests_cumul += ((double) (clock() - times->_clicks_T0Tests))/ CLOCKS_PER_SEC;
//     else         times->_clicks_TSTests_cumul += ((double) (clock() - times->_clicks_TSTests))/ CLOCKS_PER_SEC;
// }


#ifdef __cplusplus
}
#endif

#endif
