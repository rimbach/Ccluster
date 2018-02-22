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

#ifndef NEWTON_H
#define NEWTON_H

// #include <stdio.h>
#include "base/base.h"
#include "numbers/compRat.h"
#include "numbers/compApp.h"
#include "numbers/app_rat.h"
#include "geometry/compDsk.h"
#include "geometry/connCmp.h"
#include "geometry/connCmp_dsk.h"
#include "geometry/subdBox.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "caches/cacheApp.h"
#include "metadatas/metadatas.h"
#include "tstar/tstar.h"

typedef struct {
    int nflag; // 1 if true, 0 if false
    slong appPrec; // the arithmetic precision that has been used to decide
} newton_res;

newton_res newton_first_condition( compApp_t fx, compApp_t fpx, cacheApp_t cache, const compRat_t c, const realRat_t d, slong prec, metadatas_t meta);

newton_res newton_iteration( compApp_t iteration, 
                             cacheApp_t cache, 
                             const connCmp_t CC, 
                             const compRat_t c, 
                             compApp_t fcenter, 
                             compApp_t fpcenter,
                             slong prec, metadatas_t meta);

/*works in place in CC*/
newton_res newton_newton_connCmp( connCmp_t nCC,
                                  connCmp_t CC,
                                  cacheApp_t cache, 
                                  const compRat_t c,
                                  slong prec, 
                                  metadatas_t meta);

/* for julia: DEPRECATED */
/*
newton_res newton_first_condition_forjulia( cacheApp_t cache, const realRat_t cRe, const realRat_t cIm, const realRat_t d, slong prec, metadatas_t meta);

newton_res newton_iteration_forjulia( compApp_t iteration, 
                                      cacheApp_t cache, 
                                      const connCmp_t CC, 
                                      const realRat_t cRe, const realRat_t cIm,
                                      slong prec, metadatas_t meta);

newton_res newton_newton_connCmp_forjulia(  connCmp_t nCC,
                                            connCmp_t CC,
                                            cacheApp_t cache, 
                                            const realRat_t cRe, const realRat_t cIm,
                                            slong prec, 
                                            metadatas_t meta);
*/
#endif