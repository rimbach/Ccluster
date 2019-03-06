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
#include "geometry/compBox.h"

int main() {
    
    compBox_t b;
    compBox_init(b);
    printf("b: "); compBox_print(b); printf("\n");
    
    realRat_t cr, ci, w;
    realRat_init(cr); realRat_init(ci); realRat_init(w);
    realRat_set_si(cr, 1,2);
    realRat_set_si(ci, 1,3);
    realRat_set_si(w,  1,4);
    compBox_set_3realRat_int(b, cr, ci, w, 10);
    printf("b: "); compBox_print(b); printf("\n");
    
    realRat_set_si(cr, 1,27);
    realRat_set_si(ci, 1,31);
    realRat_set_si(w,  1,43);
    compBox_set_3realRat(b, cr, ci, w);
    printf("b: "); compBox_print(b); printf("\n");
    
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 1,2,1,3); 
    compBox_set_compRat_realRat_int(b, c, w, 15);
    printf("b: "); compBox_print(b); printf("\n");
    compBox_set_compRat_realRat(b, c, w);
    printf("b: "); compBox_print(b); printf("\n");
    
    compBox_t b2;
    compBox_init(b2);
    compBox_set(b2, b);
    printf("b2 copy of b1: "); compBox_print(b2); printf("\n");
    
    compRat_set_sisi(c, 1,9,1,7);
    realRat_set_si(w,  1,5);
    compBox_set_compRat_realRat_int(b, c, w, 5);
    printf("b: "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    
    compBox_inflate_realRat_inplace(b, w);
    printf("b inflated by 1/5: "); compBox_print(b); printf("\n");
    compBox_inflate_realRat(b2, b, w);
    printf("b2 = b inflated by 1/5: "); compBox_print(b2); printf("\n");
    
    printf("isless(b,b2): %d\n", compBox_isless(b,b2));
    printf("isless(b2,b): %d\n\n", compBox_isless(b2,b));
    
    compRat_set_sisi(c, 1,1,1,1);
    realRat_set_si(w,  1,1);
    compBox_set_compRat_realRat(b, c, w);
    compRat_set_sisi(c, 1,2,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b: "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("isless(b,b2): %d\n\n", compBox_isless(b,b2));
    
    compRat_set_sisi(c, 1,1,1,2);
    compBox_set_compRat_realRat(b, c, w);
    compRat_set_sisi(c, 1,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b: "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("isless(b,b2): %d\n\n", compBox_isless(b,b2));
    
    printf("-----------------------------------------\n");
    printf("Testing predicates:\n");
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 0,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 0,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    printf("----------------\n");
    realRat_set_si(w,  1,2);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 0,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, -1,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 0,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    compRat_set_sisi(c, 1,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_are_8connected b ,b2: %d\n", compBox_are_8connected(b,b2) );
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("----------------\n");
    
    printf("----------------\n");
    printf("----------------\n");
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 3,2,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("compBox_intersection_has_non_empty_interior b ,b2: %d\n", compBox_intersection_has_non_empty_interior(b,b2) );
    realRat_set_si(w,  1,2);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("compBox_intersection_has_non_empty_interior b ,b2: %d\n", compBox_intersection_has_non_empty_interior(b,b2) );
    realRat_set_si(w,  3,2);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_intersection_is_not_empty b ,b2: %d\n", compBox_intersection_is_not_empty(b,b2) );
    printf("compBox_intersection_has_non_empty_interior b ,b2: %d\n", compBox_intersection_has_non_empty_interior(b,b2) );
    
    printf("-------------------------------------------\n");
    printf("testing is strictly in\n");
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b, c, w);
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_strictly_in b ,b2: %d\n", compBox_is_strictly_in(b,b2) );
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b, c, w);
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_strictly_in b ,b2: %d\n", compBox_is_strictly_in(b,b2) );
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 1,2,1,2);
    compBox_set_compRat_realRat(b, c, w);
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_strictly_in b ,b2: %d\n", compBox_is_strictly_in(b,b2) );
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 1,4,0,4);
    compBox_set_compRat_realRat(b, c, w);
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_strictly_in b ,b2: %d\n", compBox_is_strictly_in(b,b2) );
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 3,4,0,4);
    compBox_set_compRat_realRat(b, c, w);
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b : "); compBox_print(b); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_strictly_in b ,b2: %d\n", compBox_is_strictly_in(b,b2) );
    
    printf("-------------------------------------------\n");
    printf("testing is point in box\n");
    compRat_set_sisi(c, 0,1,1,1);
    printf("c : "); compRat_print(c); printf("\n");
    printf("b2: "); compBox_print(b2); printf("\n");
    printf("compBox_is_point_in_box c ,b2: %d\n", compBox_is_point_in_box(c,b2) );
    compRat_set_sisi(c, 0,1,2,1);
    printf("c : "); compRat_print(c); printf("\n");
    printf("compBox_is_point_in_box c ,b2: %d\n", compBox_is_point_in_box(c,b2) );
    compBox_clear(b);
    realRat_clear(cr); realRat_clear(ci); realRat_clear(w);
    compRat_clear(c);
    
    printf("-------------------------------------------\n");
    printf("testing realcoeff box\n");
    
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,0,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf(" is b2 imaginary positive? %d\n", compBox_is_imaginary_positive(b2) );
    printf(" is b2 imaginary negative? %d\n", compBox_is_imaginary_negative(b2) );
    
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf(" is b2 imaginary positive? %d\n", compBox_is_imaginary_positive(b2) );
    printf(" is b2 imaginary negative? %d\n", compBox_is_imaginary_negative(b2) );
    
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf(" is b2 imaginary positive? %d\n", compBox_is_imaginary_positive(b2) );
    printf(" is b2 imaginary negative? %d\n", compBox_is_imaginary_negative(b2) );
    
    realRat_set_si(w,  1,1);
    compRat_set_sisi(c, 0,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf(" is b2 imaginary positive? %d\n", compBox_is_imaginary_positive(b2) );
    printf(" is b2 imaginary negative? %d\n", compBox_is_imaginary_negative(b2) );
    
    realRat_set_si(w,  2,1);
    compRat_set_sisi(c, 0,1,-1,1);
    compBox_set_compRat_realRat(b2, c, w);
    printf("b2: "); compBox_print(b2); printf("\n");
    printf(" is b2 imaginary positive? %d\n", compBox_is_imaginary_positive(b2) );
    printf(" is b2 imaginary negative? %d\n", compBox_is_imaginary_negative(b2) );
    
    compBox_set_conjugate(b,b2);
    printf("conjugate of b2: "); compBox_print(b); printf("\n");
    return 0;
    
}
