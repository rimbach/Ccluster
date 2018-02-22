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

#ifndef COMPDSK_H
#define COMPDSK_H

#include "numbers/realRat.h"
#include "numbers/compRat.h"

typedef struct {
    compRat center;
    realRat radius;
} compDsk;

typedef compDsk compDsk_t[1];
typedef compDsk * compDsk_ptr;

#define compDsk_centerref(X) (&(X)->center)
#define compDsk_radiusref(X) (&(X)->radius)

/* memory managment */
static __inline__ void compDsk_init(compDsk_t x) { 
    compRat_init(compDsk_centerref(x)); 
    realRat_init(compDsk_radiusref(x));
}

static __inline__ void compDsk_clear(compDsk_t x) { 
    compRat_clear(compDsk_centerref(x)); 
    realRat_clear(compDsk_radiusref(x));
}

void compDsk_init_forJulia(compDsk_t x);
void compDsk_clear_forJulia(compDsk_t x);

/* members access */
static __inline__ compRat_ptr compDsk_center_ptr(compDsk_t x) {
    return compDsk_centerref(x);
}

static __inline__ realRat_ptr compDsk_radius_ptr(compDsk_t x) {
    return compDsk_radiusref(x);
}

static __inline__ void compDsk_get_center(compRat_t c, const compDsk_t x) {
    compRat_set(c, compDsk_centerref(x));
}

static __inline__ void compDsk_get_radius(realRat_t r, const compDsk_t x) {
    realRat_set(r, compDsk_radiusref(x));
}

void compDsk_get_centerRe_forJulia(realRat_t c, const compDsk_t x);
void compDsk_get_centerIm_forJulia(realRat_t c, const compDsk_t x);
void compDsk_get_radius_forJulia(realRat_t r, const compDsk_t x);

/* printing */
void compDsk_fprint(FILE * file, const compDsk_t x);

static __inline__ void compDsk_print(const compDsk_t x) {
    return compDsk_fprint(stdout, x);
}

/* setting */
static __inline__ void compDsk_set_3realRat(compDsk_t d, const realRat_t cr, const realRat_t ci, const realRat_t r){
    compRat_set_2realRat( compDsk_centerref(d), cr, ci);
    realRat_set( compDsk_radiusref(d), r);
}

void compDsk_set_3realRat_forJulia(compDsk_t d, const realRat_t cr, const realRat_t ci, const realRat_t r);

static __inline__ void compDsk_set_compRat_realRat(compDsk_t d, const compRat_t c, const realRat_t r){
    compRat_set( compDsk_centerref(d), c);
    realRat_set( compDsk_radiusref(d), r);
}

static __inline__ void compDsk_set(compDsk_t d, const compDsk_t e){
    compRat_set( compDsk_centerref(d), compDsk_centerref(e));
    realRat_set( compDsk_radiusref(d), compDsk_radiusref(e));
}

/* Inflate */
//acts in place
static __inline__ void compDsk_inflate_realRat_inplace(compDsk_t d, const realRat_t f){
    realRat_mul( compDsk_radiusref(d), compDsk_radiusref(d), f );
}

void compDsk_inflate_realRat(compDsk_t d, const compDsk_t e, const realRat_t f);

/* geometric predicates */
int compDsk_is_point_in_dsk         ( const compRat_t p, const compDsk_t d);
int compDsk_is_point_strictly_in_dsk( const compRat_t p, const compDsk_t d);
#endif
