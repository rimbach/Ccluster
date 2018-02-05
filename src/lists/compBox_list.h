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

#ifndef COMPBOX_LIST_H
#define COMPBOX_LIST_H

#include "geometry/compBox.h"
#include "lists/gen_list.h"

typedef struct gen_list compBox_list;
typedef struct gen_list compBox_list_t[1];
typedef struct gen_list * compBox_list_ptr;

static __inline__ void compBox_clear_for_list(void * b){
    compBox_clear( (compBox_ptr) b );
}

static __inline__ int compBox_isless_for_list(const void * b1, const void * b2){
    return compBox_isless( (compBox_ptr) b1, (compBox_ptr) b2 );
}

static __inline__ void compBox_fprint_for_list(FILE * file, const void * b){
    compBox_fprint(file, (compBox_ptr) b);
}

static __inline__ void compBox_list_init(compBox_list_t l){
    gen_list_init(l, compBox_clear_for_list);
}

static __inline__ void compBox_list_swap(compBox_list_t l1, compBox_list_t l2) {
    gen_list_swap(l1, l2);
}

static __inline__ void compBox_list_clear(compBox_list_t l) {
    gen_list_clear(l);
}

static __inline__ void compBox_list_clear_for_tables(compBox_list_t l) {
    gen_list_clear_for_tables(l);
}

static __inline__ void compBox_list_push(compBox_list_t l, compBox_ptr b){
    gen_list_push(l, b);
}

static __inline__ compBox_ptr compBox_list_pop(compBox_list_t l){
    return (compBox_ptr) gen_list_pop(l);
}

static __inline__ compBox_ptr compBox_list_first(compBox_list_t l){
    return (compBox_ptr) gen_list_first(l);
}

/* do not check bound!!! */
static __inline__ compBox_ptr compBox_list_compBox_at_index(compBox_list_t l, int index){
    return (compBox_ptr) gen_list_data_at_index(l, index);
}

static __inline__ void compBox_list_insert_sorted(compBox_list_t l, compBox_ptr b){
    gen_list_insert_sorted(l, b, compBox_isless_for_list);
}

static __inline__ void compBox_list_fprint(FILE * file, const compBox_list_t l){
    gen_list_fprint(file, l, compBox_fprint_for_list);
}

static __inline__ void compBox_list_print(const compBox_list_t l){
    compBox_list_fprint(stdout, l);
}

static __inline__ int compBox_list_is_empty(const compBox_list_t l){
    return gen_list_is_empty(l);
}
static __inline__ int compBox_list_get_size(const compBox_list_t l){
    return gen_list_get_size(l);
}

/*iterator */
typedef gen_list_iterator compBox_list_iterator;
static __inline__ compBox_list_iterator compBox_list_begin(const compBox_list_t l){
    return gen_list_begin(l);
}
static __inline__ compBox_list_iterator compBox_list_next(compBox_list_iterator it){
    return gen_list_next(it);
}
static __inline__ compBox_list_iterator compBox_list_end(){
    return gen_list_end();
}
static __inline__ compBox_ptr compBox_list_elmt(compBox_list_iterator it){
    return (compBox_ptr) gen_list_elmt(it);
}

void compBox_list_init_forjulia(compBox_list_t l);
void compBox_list_clear_forjulia(compBox_list_t l);
void compBox_list_push_forjulia(compBox_list_t l, compBox_ptr b);
void compBox_list_pop_forjulia(compBox_t dest, compBox_list_t l);
int compBox_list_is_empty_forjulia(const compBox_list_t l);

#endif