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

#ifndef TSTAR_H
#define TSTAR_H

#include "base/base.h"
#include "geometry/compDsk.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "caches/cacheApp.h"
#include "tstar/pelletTest.h"
#include "metadatas/metadatas.h"

#include "math.h"

typedef struct {
    int nbOfSol;   /* the number of solutions: -1: can not decide, >=0 otherwise */
    slong appPrec; /* the arithmetic precision that has been used to decide      */
} tstar_res;

static __inline__ int tstar_res_getNbOfSol( tstar_res r) { return r.nbOfSol; }
static __inline__ slong tstar_res_getAppPrec( tstar_res r) { return r.appPrec; }

void tstar_getApproximation( compApp_poly_t res, cacheApp_t cache, slong prec, metadatas_t meta);
void tstar_taylor_shift_inplace( compApp_poly_t res, const compDsk_t d, slong prec, metadatas_t meta);
void tstar_graeffe_iterations_inplace( compApp_poly_t res, int d, slong prec, metadatas_t meta);
void tstar_graeffe_iterations_abs_two_first_coeffs( realApp_t coeff0, realApp_t coeff1, const compApp_poly_t pApprox, int N, slong prec, metadatas_t meta);

int tstar_numerical_test( compApp_poly_t pApprox, const compDsk_t d, slong prec, metadatas_t meta);

tstar_res tstar_interface( cacheApp_t cache,  /*                                               */
                           const compDsk_t d, /*                                               */
                           int max_nb_sols,   /*the maximum number of sols in the disk         */
                           int discard,       /*a flag saying if it is a discarding test or not*/
                           slong prec,        /*the "default" arithmetic precision             */
                           int depth,         /*the depth for counter                          */
                           metadatas_t meta);

tstar_res tstar_asInPaper( cacheApp_t cache,
                           const compDsk_t d,
                           int max_nb_sols,   /*the maximum number of sols in the disk          */
                           int discard,       /*a flag saying if it is a discarding test or not */
                           slong prec,        /*the "default" arithmetic precision              */
                           int depth,         /*the depth for counter                           */
                           metadatas_t meta);

tstar_res tstar_optimized( cacheApp_t cache,
                           const compDsk_t d,
                           int max_nb_sols,   /*the maximum number of sols in the disk           */
                           int discard,       /*a flag saying if it is a discarding test or not  */
                           slong prec,        /*the "default" arithmetic precision               */
                           int depth,         /*the depth for counter                            */
                           metadatas_t meta);

/*EXPERIMENTAL*/
int tstar_inclusion_test( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta);

int tstar_inclusion_radius_th9( cacheApp_t cache,
                          const compDsk_t d,
                          slong prec,
                          int depth, /* just for display*/
                          metadatas_t meta);

int tstar_numerical_test( compApp_poly_t pApprox, const compDsk_t d, slong prec, metadatas_t meta);

/* DEPRECATED */
// tstar_res tstar_count_nb_Sols( cacheApp_t cache,
//                                const compDsk_t d,
//                                int nb_sols,   /*the number of sols in d          */
//                                slong prec,        /*the "default" arithmetic precision              */
//                                int depth,         /*the depth for counter                           */
//                                metadatas_t meta);
/* for julia: DEPRECATED*/
/*int tstar_res_getNbOfSol_forJulia( tstar_res r);   */
/*slong tstar_res_getAppPrec_forJulia( tstar_res r); */

#endif