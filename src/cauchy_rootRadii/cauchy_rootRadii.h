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

#ifndef CAUCHYROOTRADII_H
#define CAUCHYROOTRADII_H


#include "base/base.h"
#include "numbers/realRat.h"
#include "numbers/realApp.h"
#include "numbers/app_rat.h"
#include "caches/cacheApp.h"
#include "caches/cacheCauchy.h"
#include "metadatas/metadatas.h"
#include "cauchy_tests/cauchy_tests.h"

#ifdef __cplusplus
extern "C" {
#endif
    
void cauchyRootRadii_root_bound( realRat_t upperBound,
                                 cacheCauchy_t cacheCau,
                                 cacheApp_t cache,
                                 metadatas_t meta );

// void cauchyRootRadii_root_radius( const compRat_t center,
//                                   realRat_t radInf, /* radInf < r_{d+1-m}(center, p) */
//                                   realRat_t radSup, /* radSup > r_{d+1-m}(center, p) */
//                                   realRat_t relativeError, /* want relativeError*radInf >= radSup */ 
//                                   slong nbOfRoots,
//                                   cacheCauchy_t cacheCau,
//                                  cacheApp_t cache,
//                                  metadatas_t meta );
    
#ifdef __cplusplus
}
#endif

#endif
