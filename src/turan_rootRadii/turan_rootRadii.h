/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#ifndef TURANROOTRADII_H
#define TURANROOTRADII_H


#include "base/base.h"
#include "numbers/realRat.h"
#include "numbers/realApp.h"
#include "numbers/app_rat.h"
#include "geometry/compDsk.h"
#include "caches/cacheApp.h"
#include "caches/cacheCauchy.h"
#include "metadatas/metadatas.h"
#include "cauchy_tests/cauchy_tests.h"

#ifdef __cplusplus
extern "C" {
#endif
    
slong turanRootRadii_root_bound( realApp_t rm,
                                 const compDsk_t Delta,
                                 slong           m,
                                 const realRat_t theta,
                                 ulong           N,
                                 const realRat_t eps,
                                 cacheCauchy_t cacheCau,
                                 cacheApp_t cache,
                                 slong prec,
                                 metadatas_t meta, int depth );

#ifdef __cplusplus
}
#endif

#endif
