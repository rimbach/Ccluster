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

int main() {

//     gen_list_t lbox;
//     gen_list_init(lbox, compBox_clear);
//     printf("lbox: "); gen_list_fprint(stdout, lbox, compBox_fprint ); printf("\n");
//     
//     compBox_t b;
//     compBox_init(b);
//     compRat_t c;
//     compRat_init(c);
//     compRat_set_sisi(c, 1,2,1,3); 
//     realRat_t w;
//     realRat_init(w);
//     realRat_set_si(w,  1,4);
//     
//     compBox_set_compRat_realRat_int(b, c, w, 15);
//     gen_list_push(lbox, b);
//     printf("lbox: "); gen_list_fprint(stdout, lbox, compBox_fprint ); printf("\n");
//     
//     realRat_clear(w);
//     compRat_clear(c);
//     gen_list_clear(lbox);
    
    compBox_list_t lbox;
    compBox_list_init (lbox);
    
    compBox_t b;
    compBox_init(b);
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 1,2,1,3); 
    realRat_t w;
    realRat_init(w);
    realRat_set_si(w,  1,4);
    compBox_set_compRat_realRat_int(b, c, w, 15);
    
    compBox_list_push(lbox, b);
    printf("lbox: "); compBox_list_fprint(stdout, lbox); printf("\n");
    
    compBox_ptr b2;
    b2 = compBox_list_pop(lbox);
    printf("lbox: "); compBox_list_fprint(stdout, lbox); printf("\n");
    
    printf("b2: "); compBox_fprint(stdout, b2); printf("\n");
    
    compBox_list_clear (lbox);
    realRat_clear(w);
    compRat_clear(c);
    compBox_clear(b2);
    
    return 0;
}