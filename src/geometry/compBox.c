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

#include <stdio.h>
#include "geometry/compBox.h"

void compBox_fprint(FILE * file, const compBox_t x){
    fprintf(file, "center: ");
    compRat_fprint(file, compBox_centerref(x));
    fprintf(file, " width: ");
    realRat_fprint(file, compBox_bwidthref(x));
    fprintf(file, " maxNbSols: %d", x->nbMSol);
}

void compBox_inflate_realRat(compBox_t d, const compBox_t e, const realRat_t f){
    compBox_set(d, e);
    compBox_inflate_realRat_inplace(d,f);
}

slong compBox_getDepth(const compBox_t b, const compBox_t initialBox) {
    realRat_t depth;
    realRat_init(depth);
    realRat_set(depth, compBox_bwidthref(b));
    realRat_div( depth, compBox_bwidthref(initialBox),  depth);
    slong res = fmpz_clog_ui( realRat_numref(depth), 2);
    
    realRat_clear(depth);
    return res;
}

/* Precondition:  b1 and b2 have disjoint centers and same width                               */
/* Specification: returns 1 if b1 and b2 are connected (8 adjacency),                          */
/*                        0 otherwise                                                          */
int compBox_are_8connected ( const compBox_t b1, const compBox_t b2 ){
    
    int res;
    compRat_t dist;
    compRat_init(dist);
    compRat_comp_distance( dist, compBox_centerref(b1), compBox_centerref(b2) );
    
    res = (realRat_cmp( compRat_realref(dist), compBox_bwidthref(b1) ) <=0) 
       && (realRat_cmp( compRat_imagref(dist), compBox_bwidthref(b1) ) <=0);
       
    compRat_clear(dist);
    
    return res;
}

/* Precondition:                                                                               */
/* Specification: returns true if (b1 \cap b2)\neq \emptyset                                   */
/*                        false otherwise                                                      */
int compBox_intersection_is_not_empty ( const compBox_t b1, const compBox_t b2 ){
    
    int res;
    compRat_t dist;
    realRat_t half, halfwidth;
    compRat_init(dist);
    realRat_init(half);
    realRat_init(halfwidth);
    
    compRat_comp_distance( dist, compBox_centerref(b1), compBox_centerref(b2) );
    realRat_set_si(half, 1,2);
    realRat_add(halfwidth, compBox_bwidthref(b1), compBox_bwidthref(b2) );
    realRat_mul(halfwidth, halfwidth, half);
    
    res = (realRat_cmp( compRat_realref(dist), halfwidth ) <=0) 
       && (realRat_cmp( compRat_imagref(dist), halfwidth ) <=0);
       
    compRat_clear(dist);
    realRat_clear(half);
    realRat_clear(halfwidth);
    
    return res;
}

/* Precondition:                                                                               */
/* Specification: returns true if interior(b1 \cap b2)\neq \emptyset                           */
/*                        false otherwise                                                      */
int compBox_intersection_has_non_empty_interior ( const compBox_t b1, const compBox_t b2 ){
    int res;
    compRat_t dist;
    realRat_t half, halfwidth;
    compRat_init(dist);
    realRat_init(half);
    realRat_init(halfwidth);
    
    compRat_comp_distance( dist, compBox_centerref(b1), compBox_centerref(b2) );
    realRat_set_si(half, 1,2);
    realRat_add(halfwidth, compBox_bwidthref(b1), compBox_bwidthref(b2) );
    realRat_mul(halfwidth, halfwidth, half);
    
    res = (realRat_cmp( compRat_realref(dist), halfwidth ) <0) 
       && (realRat_cmp( compRat_imagref(dist), halfwidth ) <0);
       
    compRat_clear(dist);
    realRat_clear(half);
    realRat_clear(halfwidth);
    
    return res;
}

/* Precondition:                                                                               */
/* Specification: returns true if b1\subset\interior{b2}                                       */
/*                        false otherwise                                                      */
int compBox_is_strictly_in ( const compBox_t b1, const compBox_t b2 ){
    int res;
    compRat_t dist;
    realRat_t halfwidth1, halfwidth2, sum;
    compRat_init(dist);
    realRat_init(halfwidth1);
    realRat_init(halfwidth2);
    realRat_init(sum);
    
    compRat_comp_distance( dist, compBox_centerref(b1), compBox_centerref(b2) );
    realRat_set_si(halfwidth1, 1,2);
    realRat_set(halfwidth2, halfwidth1);
    realRat_mul(halfwidth1, halfwidth1, compBox_bwidthref(b1));
    realRat_mul(halfwidth2, halfwidth2, compBox_bwidthref(b2));
    
    realRat_add(sum, compRat_realref(dist), halfwidth1);
    res = (realRat_cmp( sum, halfwidth2 ) <0);
    if (res) {
        realRat_add(sum, compRat_imagref(dist), halfwidth1);
        res = (realRat_cmp( sum, halfwidth2 ) <0);
    }
    
    compRat_clear(dist);
    realRat_clear(halfwidth1);
    realRat_clear(halfwidth2);
    realRat_clear(sum);
    
    return res;
}

/* Precondition:                                                                               */
/* Specification: returns true iff p is in b                                                   */
int compBox_is_point_in_box                     ( const compRat_t p,  const compBox_t b  ){
    int res;
    compRat_t dist;
    realRat_t halfwidth;
    compRat_init(dist);
    realRat_init(halfwidth);
    
    compRat_comp_distance( dist, p, compBox_centerref(b) );
    realRat_set_si(halfwidth, 1,2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(b));
    
    res = (realRat_cmp( compRat_realref(dist), halfwidth ) <=0);
    if (res) {
        res = (realRat_cmp( compRat_imagref(dist), halfwidth ) <=0);
    }
    
    compRat_clear(dist);
    realRat_clear(halfwidth);
    
    return res;
}

/* precondition                                                                                 */
/* specification returns point res = (resRe, resIm) in b that is the closest to p               */
/*               in particular if p is in b, returns p                                          */
void compBox_closest_point_on_boundary          (compRat_t res, const compRat_t p, const compBox_t b  ){
    
    realRat_t halfwidth, sum;
    realRat_init(halfwidth);
    realRat_init(sum);
    realRat_set_si(halfwidth, 1,2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(b));
    
    /* determine in which area lies the real part of p */
    realRat_add(sum, compRat_realref(compBox_centerref(b)), halfwidth);
    if (realRat_cmp( compRat_realref(p), sum ) <=0) { /* (NW, W, SW) or (N, C, S) */
        realRat_sub(sum, compRat_realref(compBox_centerref(b)), halfwidth);
        if (realRat_cmp( compRat_realref(p), sum ) >=0) { /* (N, C, S) */
            realRat_set( compRat_realref(res), compRat_realref(p) );
        }
        else { /* (NW, W, SW) */
            realRat_set( compRat_realref(res), sum );
        }
    }
    else { /* (NE, E, SE) */
        realRat_set( compRat_realref(res), sum );
    }
    
    /* determine in which area lies the imag part of p */
    realRat_add(sum, compRat_imagref(compBox_centerref(b)), halfwidth);
    if (realRat_cmp( compRat_imagref(p), sum ) <=0) { /* (SW, S, SE) or (W, C, E) */
        realRat_sub(sum, compRat_imagref(compBox_centerref(b)), halfwidth);
        if (realRat_cmp( compRat_imagref(p), sum ) >=0) { /* (W, C, E) */
            realRat_set( compRat_imagref(res), compRat_imagref(p) );
        }
        else { /* (SW, S, SE) */
            realRat_set( compRat_imagref(res), sum );
        }
    }
    else { /* (NW, N, NE) */
        realRat_set( compRat_imagref(res), sum );
    }
    
    realRat_clear(halfwidth);
    realRat_clear(sum);
}

/* RealCoeffs */
/* Precondition:                                                                               */
/* Specification: returns true only if forall p\in b, Im(p)<0                                  */
int compBox_is_imaginary_negative_strict               ( const compBox_t b  ){
    int res;
    realRat_t upper, zero;
    realRat_init(upper);
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    realRat_set_si(upper, 1,2);
    realRat_mul(upper, upper, compBox_bwidthref(b));
    realRat_add(upper, compRat_imagref(compBox_centerref(b)), upper);
    
    res = (realRat_cmp(upper, zero) < 0);
    
    realRat_clear(upper);
    realRat_clear(zero);
    return res;
}
/* Precondition:                                                                               */
/* Specification: returns true only if forall p\in b, Im(p)>0                                  */
int compBox_is_imaginary_positive_strict               ( const compBox_t b  ){
    int res;
    realRat_t lower, zero;
    realRat_init(lower);
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    realRat_set_si(lower, 1,2);
    realRat_mul(lower, lower, compBox_bwidthref(b));
    realRat_sub(lower, compRat_imagref(compBox_centerref(b)), lower);
    
    res = (realRat_cmp(lower, zero) > 0);
    
    realRat_clear(lower);
    realRat_clear(zero);
    return res;
}
/* Precondition: d is initialized                                                                              */
/* Specification: set d to the complex conjugate of b                                          */
void compBox_set_conjugate                      ( compBox_t d, const compBox_t b  ){
    compBox_set(d,b);
    realRat_neg( compRat_imagref(compBox_centerref(d)), compRat_imagref(compBox_centerref(d)) ); 
}
void compBox_conjugate_inplace                  ( compBox_t b  ){
    realRat_neg( compRat_imagref(compBox_centerref(b)), compRat_imagref(compBox_centerref(b)) );
}

/* DEPRECATED

int compBox_is_point_in_box_forJulia            ( const realRat_t pre, const realRat_t pim,  const compBox_t b  ){
    int res;
    compRat_t p;
    compRat_init(p);
    compRat_set_2realRat(p, pre, pim);
    res = compBox_is_point_in_box(p, b);
    compRat_clear(p);
    return res;
}

void compBox_closest_point_on_boundary_forJulia          (realRat_t resre, realRat_t resim, const realRat_t pre, const realRat_t pim, const compBox_t b  ){
    
    compRat_t p, res;
    compRat_init(p);
    compRat_init(res);
    compRat_set_2realRat(p, pre, pim);
    
    compBox_closest_point_on_boundary(res, p, b);
    realRat_set(resre, compRat_realref(res));
    realRat_set(resim, compRat_imagref(res));
    
    compRat_clear(p);
    compRat_clear(res);
}
*/
