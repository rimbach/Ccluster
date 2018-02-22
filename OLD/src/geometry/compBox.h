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

#ifndef COMPBOX_H
#define COMPBOX_H

#include "numbers/realRat.h"
#include "numbers/compRat.h"

typedef struct {
    compRat center;
    realRat bwidth;
    int     nbMSol;
} compBox;

typedef compBox compBox_t[1];
typedef compBox * compBox_ptr;

#define compBox_centerref(X) (&(X)->center)
#define compBox_bwidthref(X) (&(X)->bwidth)
#define compBox_nbMSolref(X) (X->nbMSol)

/* memory managment */
static __inline__ void compBox_init(compBox_t x) { 
    compRat_init(compBox_centerref(x)); 
    realRat_init(compBox_bwidthref(x));
//     x->nbMSol=-1;
}

static __inline__ void compBox_clear(compBox_t x) { 
    compRat_clear(compBox_centerref(x)); 
    realRat_clear(compBox_bwidthref(x));
}

void compBox_init_forJulia(compBox_t x);
void compBox_clear_forJulia(compBox_t x);

/* members access */
static __inline__ compRat_ptr compBox_center_ptr(compBox_t x) {
    return compBox_centerref(x);
}

static __inline__ realRat_ptr compBox_bwidth_ptr(compBox_t x) {
    return compBox_bwidthref(x);
}

static __inline__ void compBox_get_center(compRat_t c, const compBox_t x) {
    compRat_set(c, compBox_centerref(x));
}

static __inline__ void compBox_get_radius(realRat_t r, const compBox_t x) {
    realRat_set(r, compBox_bwidthref(x));
}

static __inline__ int compBox_get_nbMSol(const compBox_t x) {
    return x->nbMSol;
}

void compBox_get_centerRe_forJulia(realRat_t c, const compBox_t x);
void compBox_get_centerIm_forJulia(realRat_t c, const compBox_t x);
void compBox_get_bwidth_forJulia(realRat_t r, const compBox_t x);
int  compBox_get_nbMSol_forJulia(const compBox_t x);

/* printing */
void compBox_fprint(FILE * file, const compBox_t x);

static __inline__ void compBox_print(const compBox_t x) {
    return compBox_fprint(stdout, x);
}

/* setting */
static __inline__ void compBox_set_3realRat_int(compBox_t d, const realRat_t cr, const realRat_t ci, const realRat_t r, int nb){
    compRat_set_2realRat( compBox_centerref(d), cr, ci);
    realRat_set( compBox_bwidthref(d), r);
    d->nbMSol=nb;
}

static __inline__ void compBox_set_3realRat(compBox_t d, const realRat_t cr, const realRat_t ci, const realRat_t r){
    compRat_set_2realRat( compBox_centerref(d), cr, ci);
    realRat_set( compBox_bwidthref(d), r);
    d->nbMSol=-1;
}

static __inline__ void compBox_set_compRat_realRat_int(compBox_t d, const compRat_t c, const realRat_t r, int nb){
    compRat_set( compBox_centerref(d), c);
    realRat_set( compBox_bwidthref(d), r);
    d->nbMSol=nb;
}

static __inline__ void compBox_set_compRat_realRat(compBox_t d, const compRat_t c, const realRat_t r){
    compRat_set( compBox_centerref(d), c);
    realRat_set( compBox_bwidthref(d), r);
    d->nbMSol=-1;
}

static __inline__ void compBox_set_si(compBox_t d, slong cReNum, ulong cReDen, slong cImNum, ulong cImDen, slong widNum, ulong widDen){
    compRat_set_sisi( compBox_centerref(d), cReNum, cReDen, cImNum, cImDen);
    realRat_set_si( compBox_bwidthref(d), widNum, widDen);
    d->nbMSol=-1;
}

static __inline__ int compBox_set_str(compBox_t x, const char * strReN, const char * strReD, 
                                                   const char * strImN, const char * strImD, 
                                                   const char * strWiN, const char * strWiD, 
                                                   int b){
    if (compRat_set_str(compBox_centerref(x), strReN, strReD, strImN, strImD, b) == 0)
        return realRat_set_str(compBox_bwidthref(x), strWiN, strWiD, b);
    else
        return -1;
}

static __inline__ void compBox_set(compBox_t d, const compBox_t b){
    compRat_set( compBox_centerref(d), compBox_centerref(b));
    realRat_set( compBox_bwidthref(d), compBox_bwidthref(b));
    d->nbMSol=b->nbMSol;
}

void compBox_set_3realRat_int_forJulia(compBox_t d, const realRat_t cr, const realRat_t ci, const realRat_t r, int nb);
void compBox_set_3realRat_forJulia    (compBox_t d, const realRat_t cr, const realRat_t ci, const realRat_t r);

/* Inflate */
//acts in place
static __inline__ void compBox_inflate_realRat_inplace(compBox_t d, const realRat_t f){
    realRat_mul( compBox_bwidthref(d), compBox_bwidthref(d), f );
    d->nbMSol=-1;
}

void compBox_inflate_realRat(compBox_t d, const compBox_t e, const realRat_t f);

/* comparison */
static __inline__ int compBox_cmp( const compBox_t x, const compBox_t y ) { 
    return compRat_cmp( compBox_centerref(x), compBox_centerref(y) );
}

/* ordering */
/* Precondition:  b1 and b2 have same width                                                   */
/* Specification: returns 1 if b1<b2,                                                         */
/*                        0 otherwise (b1>b2)                                                 */
/*                b1<b2 iff ( (b1._centerRe<b2._centerRe)                                     */
/*                          or ((b1._centerRe=b2._centerRe) and (b1._centerIm<b2._centerIm)) )*/
static __inline__ int compBox_isless ( const compBox_t b1, const compBox_t b2 ) {
    return (compBox_cmp(b1, b2)<0);
}
int compBox_isless_forJulia ( const compBox_t b1, const compBox_t b2 );

/* geometric predicates */
/* Precondition:  b1 and b2 have disjoint centers and same width                               */
/* Specification: returns 1 if b1 and b2 are connected (8 adjacency),                          */
/*                        0 otherwise                                                          */
int compBox_are_8connected                      ( const compBox_t b1, const compBox_t b2 );

/* Precondition:                                                                               */
/* Specification: returns true if (b1 \cap b2)\neq \emptyset                                   */
/*                        false otherwise                                                      */
int compBox_intersection_is_not_empty           ( const compBox_t b1, const compBox_t b2 );

/* Precondition:                                                                               */
/* Specification: returns true if interior(b1 \cap b2)\neq \emptyset                           */
/*                        false otherwise                                                      */
int compBox_intersection_has_non_empty_interior ( const compBox_t b1, const compBox_t b2 );

/* Precondition:                                                                               */
/* Specification: returns true if b1\subset\interior{b2}                                       */
/*                        false otherwise                                                      */
int compBox_is_strictly_in                      ( const compBox_t b1, const compBox_t b2 );

/* Precondition:                                                                               */
/* Specification: returns true iff p is in b                                                   */
int compBox_is_point_in_box                     ( const compRat_t p,  const compBox_t b  );
int compBox_is_point_in_box_forJulia            ( const realRat_t pre, const realRat_t pim,  const compBox_t b  );

/* geometric operations */
/* precondition                                                                                 */
/* specification returns point res = (resRe, resIm) in b that is the closest to p               */
/*               in particular if p is in b, returns p                                          */
void compBox_closest_point_on_boundary          (compRat_t res, const compRat_t p, const compBox_t b  );
void compBox_closest_point_on_boundary_forJulia          (realRat_t resre, realRat_t resim, const realRat_t pre, const realRat_t pim, const compBox_t b  );

#endif