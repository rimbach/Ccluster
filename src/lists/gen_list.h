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

#ifndef GEN_LIST_H
#define GEN_LIST_H

#ifdef LISTS_INLINE_C
#define LISTS_INLINE
#else
#define LISTS_INLINE static __inline__
#endif

#include <stdlib.h>
#include <stdio.h>
#include "base/base.h"
// #include "flint/fmpq.h" /*for definition of FILE*/

#ifdef __cplusplus
extern "C" {
#endif
    
struct gen_elmt {
    void            *_elmt;
    struct gen_elmt *_next;
    /* test for union find in large lists of boxes */
    struct gen_elmt *_prev;
};

typedef void (*clear_func)(void *);

struct gen_list {
    struct gen_elmt *_begin;
    struct gen_elmt *_end;
    int              _size;
    clear_func       _clear;
    
};

typedef struct gen_list gen_list_t[1];

void gen_list_init(gen_list_t l, clear_func clear);

void gen_list_swap(gen_list_t l1, gen_list_t l2);

void gen_list_clear(gen_list_t l);
void gen_list_clear_for_tables(gen_list_t l);

/* empty the list with NO delete of elements */
void gen_list_empty(gen_list_t l);
/* copy the list, not the elements */
void gen_list_copy(gen_list_t ltarget, const gen_list_t lsrc);

void gen_list_push(gen_list_t l, void * data);

void * gen_list_pop(gen_list_t l);

void * gen_list_pop_back(gen_list_t l);

void * gen_list_first(gen_list_t l);

void * gen_list_last(gen_list_t l);

/* do not check bound!!! */
void * gen_list_data_at_index(const gen_list_t l, int index);

void gen_list_insert_sorted(gen_list_t l, void * data, int (isless_func)(const void * d1, const void * d2));
void gen_list_insert_sorted_unique(gen_list_t l, void * data, int (isless_func)(const void * d1, const void * d2),
                                                              int (isequal_func)(const void * d1, const void * d2));

void gen_list_fprint(FILE * file, const gen_list_t l, void (* print_func)(FILE *, const void *) );
void gen_list_fprintd(FILE * file, const gen_list_t l, slong digits, void (* print_func)(FILE *, const void *, slong digits) );

int gen_list_is_empty(const gen_list_t l);
int gen_list_get_size(const gen_list_t l);

/*iterator */
typedef struct gen_elmt *gen_list_iterator;

LISTS_INLINE gen_list_iterator gen_list_begin(const gen_list_t l){
    return l->_begin;
}

/* test for union find in large lists of boxes */
LISTS_INLINE gen_list_iterator gen_list_endEl(const gen_list_t l){
    return l->_end;
}

LISTS_INLINE gen_list_iterator gen_list_next(gen_list_iterator it){
    return it->_next;
}

/* test for union find in large lists of boxes */
LISTS_INLINE gen_list_iterator gen_list_prev(gen_list_iterator it){
    return it->_prev;
}

LISTS_INLINE gen_list_iterator gen_list_end(){
    return NULL;
}

LISTS_INLINE void * gen_list_elmt(gen_list_iterator it){
    return it->_elmt;
}

/* assume itnext is not NULL */
/* assume it is not NULL */
/* remove element just after it, if it is not NULL */
void * gen_list_remove_at( gen_list_t l, gen_list_iterator it);
 
#ifdef __cplusplus
}
#endif


#endif
