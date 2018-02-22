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

void mfree(void * data) { free(data); };

void mprint(FILE * file, const void * x){
    fprintf(file, "%d", *((int *) x));
}

int misless(const void * d1, const void * d2){
    return ( *((int *) d1) <= *((int *) d2) );
}

int main() {
    
    gen_list_t l;
    printf("size of l: %d\n", l[0]._size);
    gen_list_init(l, mfree);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
//     printf("size of l: %d\n", l[0]._size);
    
    int *ele = (int *)malloc(sizeof(int));
    *ele = 4;
    int *ele1;
    int *ele2 = (int *)malloc(sizeof(int));
    *ele2 = 6;
    
    ele1 = (int *) gen_list_pop(l);
    if (ele1==NULL)
        printf("l is empty!!!\n");
    else 
        printf("first el of l: %d\n", *ele1 );
    
    gen_list_push(l, ele);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    ele1 = (int *) gen_list_pop(l);
    if (ele1==NULL)
        printf("l is empty!!!\n");
    else 
        printf("first el of l: %d\n", *ele1 );
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    gen_list_push(l, ele);
    gen_list_push(l, ele2);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    ele1 = (int *) gen_list_data_at_index(l,0);
    printf("0-th el of l: %d\n", *ele1 );
    ele1 = (int *) gen_list_data_at_index(l,1);
    printf("1-th el of l: %d\n", *ele1 );
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
//     printf("2-th el of l: %d\n", *ele1 );
    
    ele1 = (int *)malloc(sizeof(int));
    *ele1 = 1;
    gen_list_insert_sorted(l, ele1, misless);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    ele1 = (int *)malloc(sizeof(int));
    *ele1 = 5;
    gen_list_insert_sorted(l, ele1, misless);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    ele1 = (int *)malloc(sizeof(int));
    *ele1 = 7;
    gen_list_insert_sorted(l, ele1, misless);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    ele1 = (int *) gen_list_pop(l);
    if (ele1==NULL)
        printf("l is empty!!!\n");
    else 
        printf("first el of l: %d\n", *ele1 );
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
     
    printf("\n ------------- test iterator with for --------------\n");
    for ( gen_list_iterator it = gen_list_begin(l); it != gen_list_end(); it = gen_list_next(it) ){
        mprint(stdout, gen_list_elmt( it )); printf(", ");
    }
    printf("\n ------------- test iterator with while --------------\n");
    gen_list_iterator it = gen_list_begin(l);
    while (it != gen_list_end()) {
        mprint(stdout, gen_list_elmt( it )); printf(", ");
        it = gen_list_next(it);
    }
    printf("\n ------------- test iterator with for, empty list --------------\n");
    gen_list_clear(l);
    gen_list_init(l, mfree);
    for ( gen_list_iterator it = gen_list_begin(l); it != gen_list_end(); it = gen_list_next(it) ){
        mprint(stdout, gen_list_elmt( it )); printf(", ");
    }
    printf("\n ------------- test iterator with while, empty list --------------\n");
    while (it != gen_list_end()) {
        mprint(stdout, gen_list_elmt( it )); printf(", ");
        it = gen_list_next(it);
    }
    printf("\n");
    
    gen_list_clear(l);
    gen_list_init(l, mfree);
    
    ele1 = (int *)malloc(sizeof(int));
    *ele1 = 7;
    gen_list_insert_sorted(l, ele1, misless);
    printf("l: ");gen_list_fprint(stdout, l, mprint ); printf("\n");
    
    gen_list_clear(l);
    
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
    
    
    return 0;
}
