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

#ifndef CAUCHYDATAS_H
#define CAUCHYDATAS_H

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
    
    realRat isoRatio; /* assumed isolation ratio when doing discarding tests */
    
//     realApp wanPrec1; /* absulute precision to which are computed power sum approximations for counting ->1/4 */
//     slong   nbPntsE1; /* number of evaluation points on the contour of the disk for counting -> set at the beginning */
//     
//     realApp wanPrec2; /* absulute precision to which are computed power sum approximations for certification */
//     slong   nbPntsE2; /* number of evaluation points on the contour of the disk for certification -> set at the beginning */
//     
//     realRat lBoundUn; /* lower bound (isoRatio-1)^d/isoRatio^d */
//     realRat uBoundUn; /* upper bound d*(isoRatio+1)/(isoRatio-1) */
//     
//     /* function for fast evaluation */
//     void(*evalPoly)(compApp_t, compApp_t, const compApp_t, slong);
    
} cauchyDatas;
    
typedef cauchyDatas cauchyDatas_t[1];
typedef cauchyDatas * cauchyDatas_ptr;

#define cauchyDatas_isoRatioref(X) (&(X)->isoRatio)
// #define cauchyDatas_wanPrec1ref(X) (&(X)->wantedPrec1)
// #define cauchyDatas_wanPrec2ref(X) (&(X)->wantedPrec2)
// #define cauchyDatas_wanPrec1ref(X) (&(X)->wantedPrec1)
// #define cauchyDatas_wanPrec2ref(X) (&(X)->wantedPrec2)

void cauchyDatas_init( cauchyDatas_t p );

void cauchyDatas_clear( cauchyDatas_t p );

// METADATAS_INLINE slong cauchyDatas_nbPntsEval1( const pwSuDatas_t p ) { return p->nbPntsEval1; }
// METADATAS_INLINE slong cauchyDatas_nbPntsEval2( const pwSuDatas_t p ) { return p->nbPntsEval2; }
METADATAS_INLINE realRat_ptr cauchyDatas_isoRatio_ptr( cauchyDatas_t p) { return cauchyDatas_isoRatioref(p); }
// METADATAS_INLINE realRat_ptr cauchyDatas_wantedPrec_ptr( pwSuDatas_t p) { return cauchyDatas_wantedPrecref(p); }

// METADATAS_INLINE void cauchyDatas_set_nbPntsEval( pwSuDatas_t p, slong nb ) { p->nbPntsEval = nb; }
// METADATAS_INLINE void cauchyDatas_set_nbPwSuComp( pwSuDatas_t p, slong nb ) { p->nbPwSuComp = nb; }
METADATAS_INLINE void cauchyDatas_set_isolaRatio_si( cauchyDatas_t p, slong num, ulong den ) { realRat_set_si(cauchyDatas_isoRatioref(p), num, den); }

#ifdef __cplusplus
}
#endif

#endif
