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

#include "pelletTest.h"

int realApp_soft_compare(const realApp_t a, const realApp_t b, slong prec){
    
    if (realApp_gt(a,b)==1) 
        return 1;
    if (realApp_lt(a,b)==1)
        return 0;
    
    realRat_t ratio;
    realApp_t am, bm;
    int res = -2;
    
    realRat_init(ratio);
    realRat_set_si(ratio, 3,2);
    realApp_init(am);
    realApp_init(bm);
    realApp_mul_realRat(bm, b, ratio, prec);
    if (realApp_le(a, bm)==1) {
        realApp_mul_realRat(am, a, ratio, prec);
        if (realApp_le(b, am)==1) {
//             realRat_clear(ratio);
//             realApp_clear(am);
//             realApp_clear(bm);
            res = -1;
        }
    }
            
    realRat_clear(ratio);
    realApp_clear(am);
    realApp_clear(bm);
    return res;
}

int compApp_poly_TkGtilda_with_sum( const compApp_poly_t f, const realApp_t s, const ulong k, slong prec){
    realApp_t abs, diff;
    int res;
    realApp_init(abs);
    realApp_init(diff);
    compApp_abs(abs, (f->coeffs)+k, prec );
    realApp_sub(diff, s, abs, prec);

    res = realApp_soft_compare( abs, diff, prec);
    realApp_clear(abs);
    realApp_clear(diff);
    return res;
}

int realApp_poly_TkGtilda_with_sum( const realApp_poly_t f, const realApp_t s, const ulong k, slong prec){
    realApp_t abs, diff;
    int res;
    realApp_init(abs);
    realApp_init(diff);
    realApp_abs(abs, (f->coeffs)+k);
    realApp_sub(diff, s, abs, prec);
    
    res = realApp_soft_compare( abs, diff, prec);
    realApp_clear(abs);
    realApp_clear(diff);
    return res;
}
