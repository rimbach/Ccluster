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
#include "numbers/compRat.h"

int main() {
    
    compRat_t c1, c2;
    compRat_init(c1);
    compRat_init(c2);
    compRat_set_sisi(c1, 1,2,1,3);
    printf("c1: "); compRat_print(c1); printf("\n");
    
    compRat_set(c2, c1);
    printf("c2 set as c1: "); compRat_print(c2); printf("\n");
    printf("\n");
    
    compRat_set_sisi(c1, 1,3,1,2);
    printf("c1: "); compRat_print(c1); printf("\n");
    printf("c2: "); compRat_print(c2); printf("\n");
    printf("cmp of c1 and c2: %d\n", compRat_cmp(c1, c2) );
    compRat_t dist;
    compRat_init(dist);
    compRat_comp_distance(dist, c1, c2);
    printf("complex distance between c1 and c2: "); compRat_print(dist); printf("\n");
    compRat_comp_distance(dist, c2, c1);
    printf("complex distance between c2 and c1: "); compRat_print(dist); printf("\n");
    printf("c1: "); compRat_print(c1); printf("\n");
    printf("c2: "); compRat_print(c2); printf("\n");
    printf("\n");
    
    compRat_set_sisi(c1, 1,1,1,2);
    compRat_set_sisi(c2, 1,2,1,2);
    printf("c1: "); compRat_print(c1); printf("\n");
    printf("c2: "); compRat_print(c2); printf("\n");
    printf("cmp of c1 and c2: %d\n", compRat_cmp(c1, c2) );
    printf("cmp of c2 and c1: %d\n", compRat_cmp(c2, c1) );
    compRat_comp_distance(dist, c1, c2);
    printf("complex distance between c1 and c2: "); compRat_print(dist); printf("\n");
    printf("\n");
    
    compRat_set_sisi(c1, 1,1,1,2);
    compRat_set_sisi(c2, 1,1,1,2);
    printf("c1: "); compRat_print(c1); printf("\n");
    printf("c2: "); compRat_print(c2); printf("\n");
    printf("cmp of c1 and c2: %d\n", compRat_cmp(c1, c2) );
    compRat_comp_distance(dist, c1, c2);
    printf("complex distance between c1 and c2: "); compRat_print(dist); printf("\n");
    printf("\n");
    
    compRat_set_sisi(c1, 1,1,1,2);
    compRat_set_sisi(c2, 1,1,1,3);
    printf("c1: "); compRat_print(c1); printf("\n");
    printf("c2: "); compRat_print(c2); printf("\n");
    printf("cmp of c1 and c2: %d\n", compRat_cmp(c1, c2) );
    compRat_comp_distance(dist, c1, c2);
    printf("complex distance between c1 and c2: "); compRat_print(dist); printf("\n");
    printf("\n");
    
    compRat_clear(c1);
    compRat_clear(c2);
    compRat_clear(dist);
    
    return 0;   
}