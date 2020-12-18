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

#ifndef CAUCHYTESTS_H
#define CAUCHYTESTS_H


#include "base/base.h"
#include "numbers/realRat.h"
#include "numbers/realApp.h"
#include "numbers/app_rat.h"
#include "caches/cacheApp.h"
#include "caches/cacheCauchy.h"
#include "metadatas/metadatas.h"

#include "acb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct {
    int nbOfSol;   /* the number of solutions: -1: can not decide, >=0 otherwise */
    slong appPrec; /* the arithmetic precision that has been used to decide      */
} cauchyTest_res;    

cauchyTest_res cauchyTest_exclusionTest( const compRat_t center,
                                          const realRat_t radius,
                                          cacheApp_t cache,
                                          cacheCauchy_t cacheCau,
                                          slong prec,
                                          metadatas_t meta, int depth);

    
#ifdef __cplusplus
}
#endif

#endif
