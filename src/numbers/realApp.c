/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "realApp.h"

/* assume x is positive and compute its e-th root */
void realApp_pos_root_ui   ( realApp_t dest, const realApp_t x, ulong e, slong prec) {
    if (! realApp_contains_zero(x) ) {
        realApp_root_ui(dest,  x, e, prec);
        return; 
    }
    
    realApp_t mid, rad;
    realApp_init(mid);
    realApp_init(rad);
    
    realApp_get_mid_realApp(mid, x);
    realApp_get_rad_realApp(rad, x);
    realApp_add( dest, mid, rad, prec );
    realApp_root_ui(dest,  dest, e, prec);
    realApp_zero(mid);
    realApp_union(dest, dest, mid, prec);
    
    realApp_init(mid);
    realApp_init(rad);
}
