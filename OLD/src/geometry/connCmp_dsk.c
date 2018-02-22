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
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_intersection_is_not_empty_compDsk(compBox_list_elmt(it) ,d))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}