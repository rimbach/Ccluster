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

void compDsk_init_forJulia(compDsk_t x){
    compDsk_init(x);
}
void compDsk_clear_forJulia(compDsk_t x){
    compDsk_clear(x);
}

void compDsk_set_3realRat_forJulia(compDsk_t d, const realRat_t cr, const realRat_t ci, const realRat_t r){
    compDsk_set_3realRat(d, cr, ci, r);
}

void compDsk_get_centerRe_forJulia(realRat_t c, const compDsk_t x) {
    realRat_set(c, compRat_realref(compDsk_centerref(x)));
}

void compDsk_get_centerIm_forJulia(realRat_t c, const compDsk_t x) {
    realRat_set(c, compRat_imagref(compDsk_centerref(x)));
}

void compDsk_get_radius_forJulia(realRat_t r, const compDsk_t x) {
    realRat_set(r, compDsk_radiusref(x));
}

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

// int compDsk_is_point_strictly_in_dsk( const compRat_t p, const compDsk_t d){
//     int res;
//     compRat_t dist;
//     realRat_t sq;
//     compRat_init(dist);
//     realRat_init(sq);
//     compRat_comp_distance(dist, compDsk_centerref(d), p);
//     realRat_mul( compRat_realref(dist), compRat_realref(dist), compRat_realref(dist));
//     realRat_mul( compRat_imagref(dist), compRat_imagref(dist), compRat_imagref(dist));
//     realRat_add( compRat_realref(dist), compRat_realref(dist), compRat_imagref(dist));
//     realRat_mul( sq, compDsk_radiusref(d), compDsk_radiusref(d));
//     
//     res = (realRat_cmp(compRat_realref(dist),sq)<0);
//     
//     compRat_clear(dist);
//     realRat_clear(sq);
//     return res;
// }