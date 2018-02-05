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

#ifndef APPRAT_H
#define APPRAT_H

#include "realRat.h"
#include "compRat.h"
#include "realApp.h"
#include "compApp.h"

#include "fmpz.h"

/* converting realRat to realApp */
static __inline__ void realApp_set_realRat( realApp_t x, const realRat_t y, slong prec ) { arb_set_fmpq (x, y, prec); }

/* converting compRat to compApp */
static __inline__ void compApp_set_compRat( compApp_t x, const compRat_t y, slong prec ) { 
    arb_set_fmpq( acb_realref(x), compRat_realref(y), prec);
    arb_set_fmpq( acb_imagref(x), compRat_imagref(y), prec);
}

static __inline__ void compApp_set_realRat( compApp_t x, const realRat_t y, slong prec ) { 
    arb_set_fmpq( acb_realref(x), y, prec);
    arb_zero( acb_imagref(x));
}

static __inline__ void compApp_setreal_realRat( compApp_t x, const realRat_t y, slong prec ) { 
    arb_set_fmpq( acb_realref(x), y, prec);
}

static __inline__ void compApp_setimag_realRat( compApp_t x, const realRat_t y, slong prec ) { 
    arb_set_fmpq( acb_imagref(x), y, prec);
}

/*converts a disk to a compApp*/
void compApp_set_compDsk( compApp_t res, const compRat_t center, const realRat_t radius, slong prec);

/*getting a realRat lying in the ball defined by a realApp */
void realApp_get_realRat( realRat_t res, realApp_t x);
/*getting a compRat lying in the ball defined by a compApp */
static __inline__ void compApp_get_compRat( compRat_t res, compApp_t x){
    realApp_get_realRat( compRat_realref(res), compApp_realref(x));
    realApp_get_realRat( compRat_imagref(res), compApp_imagref(x));
}
void compApp_get_2realRat_forjulia( realRat_t resRe, realRat_t resIm, compApp_t x);

/* arithmetic */
void realApp_mul_realRat( realApp_t x, const realApp_t y, const realRat_t z, slong prec );
void compApp_mul_realRat_in_place( compApp_t x, const realRat_t y, slong prec );

#endif