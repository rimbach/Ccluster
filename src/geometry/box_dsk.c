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

#include "geometry/box_dsk.h"

slong compDsk_getDepth(const compDsk_t d, const compBox_t initialBox){
    realRat_t depth;
    realRat_init(depth);
    
    realRat_set_si(depth,4,3);
    realRat_mul(depth, depth, compDsk_radiusref(d));
    
    realRat_div( depth, compBox_bwidthref(initialBox),  depth);
    slong res = fmpz_clog_ui( realRat_numref(depth), 2);
    
    realRat_clear(depth);
    return res;
}

void compBox_get_containing_dsk( compDsk_t d, const compBox_t b) {
    realRat_t ratio;
    realRat_init(ratio);
    realRat_set_si(ratio, 3, 4);
    compDsk_set_compRat_realRat(d, compBox_centerref(b), compBox_bwidthref(b) );
    compDsk_inflate_realRat_inplace(d, ratio);
    realRat_clear(ratio);
}

int compBox_intersection_is_not_empty_compDsk ( const compBox_t b, const compDsk_t d){
    
    int res;
    compRat_t p;
    compRat_init(p);
    /*compute the closest point of center of d in b*/
    compBox_closest_point_on_boundary(p, compDsk_centerref(d), b);
    /*check if p==d */
    res = compRat_cmp(p, compDsk_centerref(d));
    if (res==0) 
        res = 1;
    /* the center of d is not in b: check if p is in d */
    else 
        res = compDsk_is_point_in_dsk(p, d);
    
    compRat_clear(p);
    
    return res;
}

int compBox_intersection_has_non_empty_interior_compDsk ( const compBox_t b, const compDsk_t d){
    
    int res;
    compRat_t p;
    compRat_init(p);
    /*compute the closest point of center of d in b*/
    compBox_closest_point_on_boundary(p, compDsk_centerref(d), b);
    /*check if p==d */
    res = compRat_cmp(p, compDsk_centerref(d));
    if (res==0) 
        res = 1;
    /* the center of d is not in b: check if p is in d */
    else 
        res = compDsk_is_point_strictly_in_dsk(p, d);
    compRat_clear(p);
    
    return res;
}

int compBox_is_box_in_dsk ( const compBox_t b, const compDsk_t d){
    
    if (!compDsk_is_point_strictly_in_dsk(compBox_centerref(b), d))
        return 0;
    
    int res=1;
    realRat_t halfwidth;
    compRat_t corner;
    realRat_init(halfwidth);
    compRat_init(corner);
    realRat_set_si(halfwidth, 1,2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(b));
    
    /*set corner to the SW corner of b */
    realRat_sub(compRat_realref(corner), compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_sub(compRat_imagref(corner), compRat_imagref(compBox_centerref(b)), halfwidth);
    res = compDsk_is_point_in_dsk(corner, d);
    
    if (res) {/*set corner to the SE corner of b */
        realRat_add( compRat_realref(corner), compRat_realref(corner), compBox_bwidthref(b) );
        res = compDsk_is_point_in_dsk(corner, d);
    }
    
    if (res) {/*set corner to the NE corner of b */
        realRat_add( compRat_imagref(corner), compRat_realref(corner), compBox_bwidthref(b) );
        res = compDsk_is_point_in_dsk(corner, d);
    }
    
    if (res) {/*set corner to the NW corner of b */
        realRat_sub( compRat_realref(corner), compRat_realref(corner), compBox_bwidthref(b) );
        res = compDsk_is_point_in_dsk(corner, d);
    }
    
    compRat_clear(corner);
    realRat_clear(halfwidth);
    return res;    
}
