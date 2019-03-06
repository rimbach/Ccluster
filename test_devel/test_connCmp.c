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
#include "geometry/connCmp.h"

int main() {
    
    connCmp_t cc;
    connCmp_init(cc);
    printf("cc: "); connCmp_print(cc); printf("\n");
    connCmp_clear_for_tables(cc);
    
    compBox_t b;
    compBox_init(b);
    compRat_t c;
    compRat_init(c);
    realRat_t w;
    realRat_init(w);
    compRat_set_sisi(c, 0,1,0,1); 
    realRat_set_si(w,  2,1);
    compBox_set_compRat_realRat_int(b, c, w, -1);
    
    connCmp_init_compBox(cc, b);
    printf("cc: "); connCmp_print(cc); printf("\n");
    
    compBox_ptr bp = connCmp_pop(cc);
    printf("pop cc: "); compBox_print(bp); printf("\n");
//     connCmp_clear_for_tables(cc);
    
    
//     connCmp_init(cc);
    connCmp_insert_compBox(cc, b);
    printf("cc: "); connCmp_print(cc); printf("\n");
    
    compBox_t b2;
    compBox_init(b2);
    compRat_set_sisi(c, 2,1,2,1); 
    realRat_set_si(w,  2,1);
    compBox_set_compRat_realRat_int(b2, c, w, -1);
//     connCmp_push(cc, b2);
    connCmp_insert_compBox(cc, b2);
    printf("cc: "); connCmp_print(cc); printf("\n");
    
    bp = connCmp_compBox_at_index(cc, 0);
    printf("0-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 1);
    printf("1-th element of cc: "); compBox_print(bp); printf("\n");
    
    compBox_t b3;
    compBox_init(b3);
    compRat_set_sisi(c, 0,1,2,1); 
    realRat_set_si(w,  2,1);
    compBox_set_compRat_realRat_int(b3, c, w, -1);
    connCmp_insert_compBox(cc, b3);
    printf("cc: "); connCmp_print(cc); printf("\n");
    bp = connCmp_compBox_at_index(cc, 0);
    printf("0-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 1);
    printf("1-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 2);
    printf("2-th element of cc: "); compBox_print(bp); printf("\n");
    
    compBox_t b4;
    compBox_init(b4);
    compRat_set_sisi(c, 2,1,4,1); 
    realRat_set_si(w,  2,1);
    compBox_set_compRat_realRat_int(b4, c, w, -1);
    connCmp_insert_compBox(cc, b4);
    printf("cc: "); connCmp_print(cc); printf("\n");
    bp = connCmp_compBox_at_index(cc, 0);
    printf("0-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 1);
    printf("1-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 2);
    printf("2-th element of cc: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc, 3);
    printf("3-th element of cc: "); compBox_print(bp); printf("\n");
    
    compRat_clear(c);
    realRat_clear(w);
//     compBox_clear(b);
//     compBox_clear(b2);
    connCmp_clear_for_tables(cc);
    
    printf("\n--------------test merge 2 cc---------------\n");
    
    connCmp_t cc1, cc2;
    compBox_t bb1, bb2, bb3, bb4;
    compBox_init(bb1);
    compBox_init(bb2);
    compBox_init(bb3);
    compBox_init(bb4);
    compRat_t c1;
    compRat_init(c1);
    realRat_t w1;
    realRat_init(w1);
    
    compRat_set_sisi(c1, -2,1,2,1); 
    realRat_set_si(w1,  2,1);
    compBox_set_compRat_realRat_int(bb1, c1, w1, -1);
    compRat_set_sisi(c1, 0,1,2,1);
    compBox_set_compRat_realRat_int(bb2, c1, w1, -1);
    connCmp_init_compBox(cc1, bb1);
    connCmp_insert_compBox(cc1, bb2);
    
    compRat_set_sisi(c1, 0,1,0,1); 
    compBox_set_compRat_realRat_int(bb3, c1, w1, -1);
    compRat_set_sisi(c1, 2,1,2,1);
    compBox_set_compRat_realRat_int(bb4, c1, w1, -1);
    connCmp_init_compBox(cc2, bb3);
    connCmp_insert_compBox(cc2, bb4);
    
    printf("before merge:\n");
    printf("cc1: "); connCmp_print(cc1); printf("\n");
    printf("cc2: "); connCmp_print(cc2); printf("\n");
    
    printf("after merge:\n");
    connCmp_merge_2_connCmp( cc1, cc2 );
    printf("cc1: "); connCmp_print(cc1); printf("\n");
    printf("cc2: "); connCmp_print(cc2); printf("\n");
    
    compBox_ptr bbp;
    bbp = connCmp_compBox_at_index(cc1, 0);
    printf("0-th element of cc1: "); compBox_print(bbp); printf("\n");
    bbp = connCmp_compBox_at_index(cc1, 1);
    printf("1-th element of cc1: "); compBox_print(bbp); printf("\n");
    bbp = connCmp_compBox_at_index(cc1, 2);
    printf("2-th element of cc1: "); compBox_print(bbp); printf("\n");
    bbp = connCmp_compBox_at_index(cc1, 3);
    printf("3-th element of cc1: "); compBox_print(bbp); printf("\n");
    
    compRat_clear(c1);
    realRat_clear(w1);
    connCmp_clear_for_tables(cc1);
    connCmp_clear_for_tables(cc2);
    
    printf("\n--------------test Real coeff---------------\n");
    
//     connCmp_t cc1, cc2;
//     compBox_t bb1, bb2, bb3, bb4;
    compBox_init(bb1);
    compBox_init(bb2);
    compBox_init(bb3);
    compBox_init(bb4);
//     compRat_t c1;
    compRat_init(c1);
//     realRat_t w1;
    realRat_init(w1);
    
    compRat_set_sisi(c1, -2,1,2,1); 
    realRat_set_si(w1,  2,1);
    compBox_set_compRat_realRat_int(bb1, c1, w1, -1);
    compRat_set_sisi(c1, 0,1,2,1);
    compBox_set_compRat_realRat_int(bb2, c1, w1, -1);
    connCmp_init_compBox(cc1, bb1);
    connCmp_insert_compBox(cc1, bb2);
    
    compRat_set_sisi(c1, 0,1,0,1); 
    compBox_set_compRat_realRat_int(bb3, c1, w1, -1);
    compRat_set_sisi(c1, 2,1,2,1);
    compBox_set_compRat_realRat_int(bb4, c1, w1, -1);
    connCmp_insert_compBox(cc1, bb3);
    connCmp_insert_compBox(cc1, bb4);
    
    connCmp_init(cc2);
    connCmp_set_conjugate(cc2,cc1);
    
    printf("cc1: "); connCmp_print(cc1); printf("\n");
    bp = connCmp_compBox_at_index(cc1, 0);
    printf("0-th element of cc1: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc1, 1);
    printf("1-th element of cc1: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc1, 2);
    printf("2-th element of cc1: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc1, 3);
    printf("3-th element of cc1: "); compBox_print(bp); printf("\n");
    
    printf("conjugate of cc1: "); connCmp_print(cc2); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 0);
    printf("0-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 1);
    printf("1-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 2);
    printf("2-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 3);
    printf("3-th element of cc2: "); compBox_print(bp); printf("\n");
    
    connCmp_clear_for_tables(cc2);
    connCmp_init(cc2);
    connCmp_set_conjugate_closure(cc2, cc1);
    printf("conjugate closure of cc1: "); connCmp_print(cc2); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 0);
    printf("0-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 1);
    printf("1-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 2);
    printf("2-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 3);
    printf("3-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 4);
    printf("4-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 5);
    printf("5-th element of cc2: "); compBox_print(bp); printf("\n");
    bp = connCmp_compBox_at_index(cc2, 6);
    printf("6-th element of cc2: "); compBox_print(bp); printf("\n");
    
    compRat_clear(c1);
    realRat_clear(w1);
    connCmp_clear_for_tables(cc1);
    connCmp_clear_for_tables(cc2);
    
    return 0;
}
