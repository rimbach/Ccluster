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

#ifndef REALINTROOTRADII_H
#define REALINTROOTRADII_H

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
    
void realIntRootRadii_getApproximation_real( realApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta);
void realIntRootRadii_getApproximation_comp( compApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta);

void realIntRootRadii_taylor_shift_inplace_real( realApp_poly_t res, slong centerRe, slong prec, metadatas_t meta);
void realIntRootRadii_taylor_shift_inplace_comp( compApp_poly_t res, slong centerRe, slong centerIm, slong prec, metadatas_t meta);

/* assume i<j<k */
/* assume absPi=|Pi|, absPj=|Pj|, absPk=|Pk| are approximations of integers */
/* decide if [j,log|Pj|] lies below of on the line passing trough [i,log|Pi|] and [k,log|Pk|]*/
/* returns 1 if yes */
/*         0 if no  */
/*        -1 if it can not be decided */
int realIntRootRadii_liesBelow( slong i, const realApp_t absPi,
                             slong j, const realApp_t absPj,
                             slong k, const realApp_t absPk,
                             slong prec );

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* returns 0 if needs more precision on the coeffs */
/* otherwise returns the length of the convex hull */
slong realIntRootRadii_convexHull( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec );

/* compute N graeffe iterations and convex hull */
/* if precision is not sufficient, returns a number i<N and an empty convex hull */
/* otherwise returns N and the convex hull */

int   realIntRootRadii_GraeffeAndCH_real ( slong convexHull[], slong * lenCH, slong * nprec, realApp_poly_t absCoeffs,
                                           realApp_poly_t pApprox, int N, slong prec, metadatas_t meta );

int   realIntRootRadii_GraeffeAndCH_comp ( slong convexHull[], slong * lenCH, slong * nprec, realApp_poly_t absCoeffs,
                                           compApp_poly_t pApprox, int N, slong prec, metadatas_t meta );

/* returns the precision used to carry out root radii */
slong realIntRootRadii_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                                  cacheApp_t cache,        /* polynomial */
                                  slong prec,
                                  metadatas_t meta );

/* returns the precision used to carry out root radii */
slong realIntRootRadii_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerIm,
                                  cacheApp_t cache,        /* polynomial */
                                  slong prec,
                                  metadatas_t meta );

void realIntRootRadii_connectedComponents( compAnn_list_t annulii, slong prec );  /* list of annulii */

void realIntRootRadii_containsRealRoot( compAnn_list_t annulii, cacheApp_t cache, slong prec );  /* list of annulii */

/* DEPRECATED version */
int realIntRootRadiiOLD_Ngraeffe_iterations_inplace_real( realApp_poly_t res, int N, slong prec, metadatas_t meta);
int realIntRootRadiiOLD_Ngraeffe_iterations_inplace_comp( compApp_poly_t res, int N, slong prec, metadatas_t meta);

/* assume convexHull is already initialized, and contains enough space for a convex hull of len points */
/* returns 0 if needs more precision on the coeffs */
/* otherwise returns the length of the convex hull */
slong realIntRootRadiiOLD_convexHull( slong * convexHull, const realApp_ptr abscoeffs, slong len, slong prec );

/* returns the precision used to carry out root radii */
slong realIntRootRadiiOLD_rootRadii( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerRe,
                                  cacheApp_t cache,        /* polynomial */
                                  slong prec,
                                  metadatas_t meta );

/* returns the precision used to carry out root radii */
slong realIntRootRadiiOLD_rootRadii_imagCenter( compAnn_list_t annulii,  /* list of annulii */
                                  slong centerIm,
                                  cacheApp_t cache,        /* polynomial */
                                  slong prec,
                                  metadatas_t meta );

#ifdef __cplusplus
}
#endif

#endif
