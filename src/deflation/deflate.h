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
    
#define connCmp_isDefref(X) ( (X)->isDef)
#define connCmp_degDeref(X) ( (X)->degDe)
#define connCmp_defPoref(X) (&(X)->defPo)
#define connCmp_sumAbref(X) (&(X)->sumAb)
    
/* memory managment */
void deflate_connCmp_init  (connCmp_t x);

void deflate_connCmp_clear (connCmp_t x);

/* setting */
void deflate_set( connCmp_t x, cacheApp_t cache, const compDsk_t disk, int nbSols, slong prec, metadatas_t meta );
void deflate_copy( connCmp_t dest, const connCmp_t src );

tstar_res deflate_tstar_test( const connCmp_t CC, cacheApp_t cache, const compDsk_t d, int max_nb_sols, slong prec, metadatas_t meta);

#ifdef __cplusplus
}
#endif

#endif
