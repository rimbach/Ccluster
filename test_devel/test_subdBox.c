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
#include "lists/gen_list.h"
#include "lists/compBox_list.h"
#include "geometry/subdBox.h"

int main() {
    
    compBox_list_t lbox;
    compBox_list_init(lbox);
    
    compBox_t b;
    compBox_init(b);
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 0,1,0,1); 
    realRat_t w;
    realRat_init(w);
    realRat_set_si(w,  2,1);
    compBox_set_compRat_realRat_int(b, c, w, 15);
    printf("b: "); compBox_print(b); printf("\n");
    
    subdBox_quadrisect(lbox, b);
    printf("lbox: "); compBox_list_fprint(stdout, lbox); printf("\n");
//     printf("lbox->_begin: %d\n", (char *) (&lbox->_begin) );
    
//     compBox_ptr b2;
//     b2 = compBox_list_pop(lbox);
//     printf("lbox: "); compBox_list_fprint(stdout, lbox); printf("\n");
//     printf("lbox->_begin: %d\n", (char *) (&lbox->_begin) );
//     printf("b2: "); compBox_fprint(stdout, b2); printf("\n");
    
//     compBox_clear(b2);
//     free(b2);
    
    compBox_list_clear (lbox);
    
    compBox_list_init(lbox);
    compDsk_t d;
    compDsk_init(d);
    compRat_set_sisi(c, 1,1,1,1); 
    realRat_set_si(w,  1,1);
    compDsk_set_compRat_realRat(d, c, w);
    subdBox_quadrisect_intersect_compDsk(lbox, b, d);
    printf("lbox: "); compBox_list_fprint(stdout, lbox); printf("\n");
    
    compBox_list_clear (lbox);
    realRat_clear(w);
    compRat_clear(c);
    compBox_clear(b);
    compDsk_clear(d);
    
    return 0;
}