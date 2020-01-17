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

#ifndef PWSUDATAS_H
#define PWSUDATAS_H

#ifdef METADATAS_INLINE_C
#define METADATAS_INLINE
#else
#define METADATAS_INLINE static __inline__
#endif

#include "base/base.h"
#include "numbers/realRat.h"
#include "numbers/compApp.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct{
    
    realRat wantedPrec; /* absulute precision to which are computed power sum approximations ->1/4 */
    slong   nbPntsEval; /* number of evaluation points on the contour of the disk -> set at the beginning */
    /* for experimental version */
    realRat isolaRatio; /* assumed isolation ratio when doing discarding tests */
    slong   nbPwSuComp; /* number of power sums used to confirm discarding tests */
    /* function for fast evaluation */
    void(*evalPoly)(compApp_t, compApp_t, const compApp_t, slong);
    
} pwSuDatas;
    
typedef pwSuDatas pwSuDatas_t[1];
typedef pwSuDatas * pwSuDatas_ptr;

#define pwSuDatas_wantedPrecref(X) (&(X)->wantedPrec)
#define pwSuDatas_isolaRatioref(X) (&(X)->isolaRatio)

void pwSuDatas_init( pwSuDatas_t p );

void pwSuDatas_clear( pwSuDatas_t p );

METADATAS_INLINE slong pwSuDatas_nbPntsEval( const pwSuDatas_t p ) { return p->nbPntsEval; }
METADATAS_INLINE slong pwSuDatas_nbPwSuComp( const pwSuDatas_t p ) { return p->nbPwSuComp; }
METADATAS_INLINE realRat_ptr pwSuDatas_isolaRatio_ptr( pwSuDatas_t p) { return pwSuDatas_isolaRatioref(p); }
METADATAS_INLINE realRat_ptr pwSuDatas_wantedPrec_ptr( pwSuDatas_t p) { return pwSuDatas_wantedPrecref(p); }

METADATAS_INLINE void pwSuDatas_set_nbPntsEval( pwSuDatas_t p, slong nb ) { p->nbPntsEval = nb; }
METADATAS_INLINE void pwSuDatas_set_nbPwSuComp( pwSuDatas_t p, slong nb ) { p->nbPwSuComp = nb; }
METADATAS_INLINE void pwSuDatas_set_isolaRatio_si( pwSuDatas_t p, slong num, ulong den ) { realRat_set_si(pwSuDatas_isolaRatioref(p), num, den); }

#ifdef __cplusplus
}
#endif

#endif
