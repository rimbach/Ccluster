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

#include "geometry/connCmp_dsk.h"

/*Precondition:                                                   */
/*Specification: returns true if interior(cc \cap d)\neq \emptyset*/ 
/*                       false otherwise                          */
int connCmp_intersection_has_non_empty_interior_compDsk( const connCmp_t cc, const compDsk_t d ){
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_intersection_has_non_empty_interior_compDsk(compBox_list_elmt(it) ,d))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}

/*Precondition:                                            */
/*Specification: returns true if (cc \cap d)\neq \emptyset */
/*                       false otherwise                   */
/*               boxes of cc and b have not necessarily the same width */
int connCmp_intersection_is_not_empty_compDsk( const connCmp_t cc, const compDsk_t d ){
    
    realRat_t temp;
    realRat_init(temp);
    realRat_add(temp, compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    if (realRat_cmp(temp, connCmp_infReref(cc)) < 0) {
        realRat_clear(temp);
        return 0;
    }
    
    realRat_sub(temp, compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    if (realRat_cmp(temp, connCmp_supReref(cc)) > 0) {
        realRat_clear(temp);
        return 0;
    }
    
    realRat_add(temp, compRat_imagref(compDsk_centerref(d)), compDsk_radiusref(d));
    if (realRat_cmp(temp, connCmp_infImref(cc)) < 0) {
        realRat_clear(temp);
        return 0;
    }
    
    realRat_sub(temp, compRat_imagref(compDsk_centerref(d)), compDsk_radiusref(d));
    if (realRat_cmp(temp, connCmp_supImref(cc)) > 0) {
        realRat_clear(temp);
        return 0;
    }
    realRat_clear(temp);
    
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_intersection_is_not_empty_compDsk(compBox_list_elmt(it) ,d))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}
