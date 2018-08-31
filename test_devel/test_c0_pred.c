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

#include "numbers/realRat.h"
#include "numbers/compRat.h"
#include "numbers/realApp.h"
#include "numbers/compApp.h"
#include "numbers/app_rat.h"

#include "tstar/tstar.h"

int main() {
    
    slong prec = 53;
    
    realRat_poly_t p1,p2;
    realRat_poly_init(p1);
    realRat_poly_init(p2);
    
    realRat_poly_set_coeff_si_ui(p1, 0,-1,1);
    realRat_poly_set_coeff_si_ui(p1, 1,1,1);
    realRat_poly_set_coeff_si_ui(p2, 0,-2,1);
    realRat_poly_set_coeff_si_ui(p2, 1,1,1);
    realRat_poly_mul(p1,p1,p2);
    
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c,0,1,0,1);
    realRat_t w;
    realRat_init(w);
    realRat_set_si(w,1,2);
    compDsk_t d;
    compDsk_init(d);
    compDsk_set_compRat_realRat(d, c, w);
    
    
    compApp_poly_t p;
    compApp_poly_init(p);
    compApp_poly_set_realRat_poly( p, p1, prec);
    
    printf("C0? %d\n", tstar_predicate_C0( p, d, prec));
    
    compDsk_clear(d);
    realRat_clear(w);
    compRat_clear(c);
    realRat_poly_clear(p1);
    realRat_poly_clear(p2);
    compApp_poly_clear(p);
    
    
}