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

#ifndef REALAPPROOTRADII_H
#define REALAPPROOTRADII_H

#ifdef ROOTRAD_INLINE_C
#define ROOTRAD_INLINE
#else
#define ROOTRAD_INLINE static __inline__
#endif

#include "base/base.h"
#include "geometry/compAnn.h"
#include "geometry/subdBox.h"
#include "geometry/connCmp_union_find.h"
#include "lists/compAnn_list.h"
#include "caches/cacheApp.h"
#include "lists/connCmp_list.h"
#include "metadatas/metadatas.h"

#ifdef __cplusplus
extern "C" {
#endif
    
/* returns the precision used to carry out root radii */
slong realApp_rootRadii_fromZero( compAnn_list_t annulii,  /* list of annulii */
                                  cacheApp_t cache,        /* polynomial */
                                  const realRat_t delta,
                                  slong prec,
                                  metadatas_t meta );

void realApp_rootRadii_connectedComponents( compAnn_list_t annulii, slong prec );

#ifdef __cplusplus
}
#endif

#endif
