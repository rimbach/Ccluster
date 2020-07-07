/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#ifndef COMPANN_LIST_H
#define COMPANN_LIST_H

#ifdef LISTS_INLINE_C
#define LISTS_INLINE
#else
#define LISTS_INLINE static __inline__
#endif

#include "geometry/compAnn.h"
#include "lists/gen_list.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct gen_list compAnn_list;
typedef struct gen_list compAnn_list_t[1];
typedef struct gen_list * compAnn_list_ptr;

LISTS_INLINE void compAnn_clear_for_list(void * b){
    compAnn_clear( (compAnn_ptr) b );
}

LISTS_INLINE int compAnn_isless_for_list(const void * b1, const void * b2){
//     return (compAnn_indMaxref( (compAnn_ptr) b1 ) > compAnn_indMaxref( (compAnn_ptr) b2 ));
    return compAnn_isless( (compAnn_ptr) b1, (compAnn_ptr) b2 );
}

LISTS_INLINE void compAnn_fprint_for_list(FILE * file, const void * b){
    compAnn_fprint(file, (compAnn_ptr) b);
}

LISTS_INLINE void compAnn_fprintd_for_list(FILE * file, const void * b, slong digits){
    compAnn_fprintd(file, (compAnn_ptr) b, digits);
}

LISTS_INLINE void compAnn_list_init(compAnn_list_t l){
    gen_list_init(l, compAnn_clear_for_list);
}

LISTS_INLINE void compAnn_list_swap(compAnn_list_t l1, compAnn_list_t l2) {
    gen_list_swap(l1, l2);
}

LISTS_INLINE void compAnn_list_clear(compAnn_list_t l) {
    gen_list_clear(l);
}

LISTS_INLINE void compAnn_list_clear_for_tables(compAnn_list_t l) {
    gen_list_clear_for_tables(l);
}

/*empty the list with NO delete of elements*/
LISTS_INLINE void compAnn_list_empty(compAnn_list_t l){
    gen_list_empty(l);
}
/* copy the list, not the elements */
LISTS_INLINE void compAnn_list_copy(compAnn_list_t ltarget, const compAnn_list_t lsrc){
    gen_list_copy( ltarget, lsrc );
}

LISTS_INLINE void compAnn_list_push(compAnn_list_t l, compAnn_ptr b){
    gen_list_push(l, b);
}

LISTS_INLINE compAnn_ptr compAnn_list_pop(compAnn_list_t l){
    return (compAnn_ptr) gen_list_pop(l);
}

LISTS_INLINE compAnn_ptr compAnn_list_first(compAnn_list_t l){
    return (compAnn_ptr) gen_list_first(l);
}

LISTS_INLINE compAnn_ptr compAnn_list_last(compAnn_list_t l){
    return (compAnn_ptr) gen_list_last(l);
}

/* do not check bound!!! */
LISTS_INLINE compAnn_ptr compAnn_list_compAnn_at_index(compAnn_list_t l, int index){
    return (compAnn_ptr) gen_list_data_at_index(l, index);
}

LISTS_INLINE void compAnn_list_insert_sorted(compAnn_list_t l, compAnn_ptr b){
    gen_list_insert_sorted(l, b, compAnn_isless_for_list);
}

LISTS_INLINE void compAnn_list_fprint(FILE * file, const compAnn_list_t l){
    gen_list_fprint(file, l, compAnn_fprint_for_list);
}

LISTS_INLINE void compAnn_list_fprintd(FILE * file, const compAnn_list_t l, slong digits){
    gen_list_fprintd(file, l, digits, compAnn_fprintd_for_list);
}

LISTS_INLINE void compAnn_list_print(const compAnn_list_t l){
    compAnn_list_fprint(stdout, l);
}

LISTS_INLINE void compAnn_list_printd(const compAnn_list_t l, slong digits){
    compAnn_list_fprintd(stdout, l, digits);
}

LISTS_INLINE int compAnn_list_is_empty(const compAnn_list_t l){
    return gen_list_is_empty(l);
}
LISTS_INLINE int compAnn_list_get_size(const compAnn_list_t l){
    return gen_list_get_size(l);
}

/*iterator */
typedef gen_list_iterator compAnn_list_iterator;
LISTS_INLINE compAnn_list_iterator compAnn_list_begin(const compAnn_list_t l){
    return gen_list_begin(l);
}
LISTS_INLINE compAnn_list_iterator compAnn_list_next(compAnn_list_iterator it){
    return gen_list_next(it);
}
LISTS_INLINE compAnn_list_iterator compAnn_list_end(){
    return gen_list_end();
}
LISTS_INLINE compAnn_ptr compAnn_list_elmt(compAnn_list_iterator it){
    return (compAnn_ptr) gen_list_elmt(it);
}    

/* assume it is not NULL */
/* remove element just after it, if it is not NULL */
LISTS_INLINE compAnn_ptr compAnn_list_remove_at( compAnn_list_t l, compAnn_list_iterator it){
    return (compAnn_ptr) gen_list_remove_at(l, it);
}

#ifdef __cplusplus
}
#endif

#endif
