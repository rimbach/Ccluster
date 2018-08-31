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

int main() {
    
    slong prec = 10;
    realApp_t x,x2;
    realApp_init(x);
    realApp_init(x2);
    realRat_t q;
    realRat_init(q);
    
    realRat_set_si(q,1,3);
    realApp_set_realRat(x, q, prec);
    printf("1/3: "); realApp_printd(x,4); printf("\n");
    realApp_mul_realRat(x2, x, q, prec);
    printf("(1/3)^2: "); realApp_printd(x2,4); printf("\n");
    
    compApp_t y;
    compApp_init(y);
    compRat_t p;
    compRat_init(p);
    compRat_set_sisi(p, 1,3,2,3);
    compApp_set_compRat(y, p, prec);
    printf("1/3 + i2/3: "); compApp_printd(y,4); printf("\n");
    
    realRat_set_si(q,1,9);
    compApp_setreal_realRat(y, q, prec);
    printf("1/9 + i2/3: "); compApp_printd(y,4); printf("\n");
    realRat_set_si(q,1,9);
    compApp_setimag_realRat(y, q, prec);
    printf("1/9 + i1/9: "); compApp_printd(y,4); printf("\n");
    compApp_mul_realRat_in_place(y, q, prec);
    printf("(1/9)^2 + i(1/9)^2: "); compApp_printd(y,4); printf("\n");
    
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 1,2,3,4);
    realRat_t r;
    realRat_init(r);
    realRat_set_si(r,1,5);
    compApp_t ball;
    compApp_init(ball);
    compApp_set_compDsk( ball, c, r, prec);
    printf("ball: "); compApp_printd(ball,4); printf("\n");
    realRat_clear(r);
    compRat_clear(c);
    compApp_clear(ball);
    
    realRat_clear(q);
    realApp_clear(x);
    realApp_clear(x2);
    compApp_clear(y);
    compRat_clear(p);
    
    return 0;
}
