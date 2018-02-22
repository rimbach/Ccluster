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

int main() {
    
    compDsk_t d;
    compDsk_init(d);
    
    printf("d: "); compDsk_print(d); printf("\n");
    
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 1,2,1,3); 
    
    realRat_t r;
    realRat_init(r);
    realRat_set_si(r, 1,5); 
    
    compDsk_set_compRat_realRat(d, c, r);
    printf("d: "); compDsk_print(d); printf("\n");
    
    compDsk_t d2;
    compDsk_init(d2);
    compDsk_set(d2,d);
    printf("d2: "); compDsk_print(d2); printf("\n");
    
    realRat_set_si(r, 1,6);
    compDsk_set_compRat_realRat(d, c, r);
    
    printf("d: "); compDsk_print(d); printf("\n");
    printf("d2: "); compDsk_print(d2); printf("\n");
    
    compDsk_t d3;
    compDsk_init(d3);
    realRat_t r1, r2, r3;
    realRat_init(r1); realRat_init(r3); realRat_init(r3);
    realRat_set_si(r1, 1,10); realRat_set_si(r2, 1,20); realRat_set_si(r3, 1,30);
    compDsk_set_3realRat(d3, r1, r2, r3);
    printf("d3: "); compDsk_print(d3); printf("\n");
    
    /* inflate */
    compDsk_inflate_realRat_inplace(d, r);
    printf("d after inflate by 1/6: "); compDsk_print(d); printf("\n");
    compDsk_inflate_realRat(d2, d, r);
    printf("d2 after inflate d by 1/6: "); compDsk_print(d2); printf("\n");
    printf("d: "); compDsk_print(d); printf("\n");
    
    compRat_clear(c);
    realRat_clear(r); realRat_clear(r1); realRat_clear(r2); realRat_clear(r3);
    compDsk_clear(d);
    compDsk_clear(d2);
    compDsk_clear(d3);
}