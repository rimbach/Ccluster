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
#include "geometry/compDsk.h"

void compDsk_fprint(FILE * file, const compDsk_t x){
    fprintf(file, "center: ");
    compRat_fprint(file, compDsk_centerref(x));
    fprintf(file, " radius: ");
    realRat_fprint(file, compDsk_radiusref(x));
}

void compDsk_inflate_realRat(compDsk_t d, const compDsk_t e, const realRat_t f){
    compDsk_set(d, e);
    compDsk_inflate_realRat_inplace(d,f);
}

int compDsk_is_point_in_dsk( const compRat_t p, const compDsk_t d){
    int res;
    compRat_t dist;
    compRat_init(dist);
    compRat_comp_distance(dist, compDsk_centerref(d), p);
    
    if (  (realRat_cmp(compRat_realref(dist), compDsk_radiusref(d))>0)
        ||(realRat_cmp(compRat_imagref(dist), compDsk_radiusref(d))>0) ){
        res = 0;
    }
    else {
        realRat_t sq;
        realRat_init(sq);
        realRat_mul( compRat_realref(dist), compRat_realref(dist), compRat_realref(dist));
        realRat_mul( compRat_imagref(dist), compRat_imagref(dist), compRat_imagref(dist));
        realRat_add( compRat_realref(dist), compRat_realref(dist), compRat_imagref(dist));
        realRat_mul( sq, compDsk_radiusref(d), compDsk_radiusref(d));
        res = (realRat_cmp(compRat_realref(dist),sq)<=0);
        realRat_clear(sq);
    }
    compRat_clear(dist);
    return res;
}

int compDsk_is_point_strictly_in_dsk( const compRat_t p, const compDsk_t d){
    int res;
    compRat_t dist;
    compRat_init(dist);
    compRat_comp_distance(dist, compDsk_centerref(d), p);
    
    if (  (realRat_cmp(compRat_realref(dist), compDsk_radiusref(d))>=0)
        ||(realRat_cmp(compRat_imagref(dist), compDsk_radiusref(d))>=0) ){
        res = 0;
    }
    else {
        realRat_t sq;
        realRat_init(sq);
        realRat_mul( compRat_realref(dist), compRat_realref(dist), compRat_realref(dist));
        realRat_mul( compRat_imagref(dist), compRat_imagref(dist), compRat_imagref(dist));
        realRat_add( compRat_realref(dist), compRat_realref(dist), compRat_imagref(dist));
        realRat_mul( sq, compDsk_radiusref(d), compDsk_radiusref(d));
        res = (realRat_cmp(compRat_realref(dist),sq)<0);
        realRat_clear(sq);
    }
    compRat_clear(dist);
    return res;
}

int compDsk_intersect_compDsk       ( const compDsk_t a, const compDsk_t b) {
    
    compRat_t cbmca;
    compRat_init(cbmca);
    compRat_sub(cbmca, compDsk_centerref(b), compDsk_centerref(a) );
    realRat_pow_si( compRat_realref(cbmca), compRat_realref(cbmca), 2); 
    realRat_pow_si( compRat_imagref(cbmca), compRat_imagref(cbmca), 2);
    realRat_add( compRat_realref(cbmca), compRat_realref(cbmca), compRat_imagref(cbmca) );
    realRat_add( compRat_imagref(cbmca), compDsk_radiusref(a), compDsk_radiusref(b));
    realRat_pow_si( compRat_imagref(cbmca), compRat_imagref(cbmca), 2);
    realRat_sub( compRat_realref(cbmca), compRat_realref(cbmca), compRat_imagref(cbmca) );
    
    int res = ( realRat_sgn(compRat_realref(cbmca)) <= 0 );
    
    compRat_clear(cbmca);
    
    return res;
}

/* RealCoeffs */
/* Precondition:                                                                               */
/* Specification: returns true only if forall p\in cc, Im(p)>0                                 */
int compDsk_is_imaginary_positive_strict        ( const compDsk_t d  ){
    int res;
    realRat_t lower, zero;
    realRat_init(lower);
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    realRat_set(lower, compDsk_radiusref(d));
    realRat_sub(lower, compRat_imagref(compDsk_centerref(d)), lower);
    
    res = (realRat_cmp(lower, zero) > 0);
    
    realRat_clear(lower);
    realRat_clear(zero);
    return res;
}

/*
int compDsk_is_point_strictly_in_dsk( const compRat_t p, const compDsk_t d){
    int res;
    compRat_t dist;
    realRat_t sq;
    compRat_init(dist);
    realRat_init(sq);
    compRat_comp_distance(dist, compDsk_centerref(d), p);
    realRat_mul( compRat_realref(dist), compRat_realref(dist), compRat_realref(dist));
    realRat_mul( compRat_imagref(dist), compRat_imagref(dist), compRat_imagref(dist));
    realRat_add( compRat_realref(dist), compRat_realref(dist), compRat_imagref(dist));
    realRat_mul( sq, compDsk_radiusref(d), compDsk_radiusref(d));
    
    res = (realRat_cmp(compRat_realref(dist),sq)<0);
    
    compRat_clear(dist);
    realRat_clear(sq);
    return res;
}
*/
