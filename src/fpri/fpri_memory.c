/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "fpri.h"


static void  *(*__fpri_allocate_func) (size_t) = _fpri_malloc;
static void  *(*__fpri_callocate_func) (size_t, size_t) = _fpri_calloc;
static void  *(*__fpri_reallocate_func) (void *, size_t) = _fpri_realloc;
static void  (*__fpri_free_func) (void *) = _fpri_free;

void __fpri_set_memory_functions( void *(*alloc_func) (size_t),
                                  void *(*calloc_func) (size_t, size_t), 
                                  void *(*realloc_func) (void *, size_t),
                                  void (*free_func) (void *) ){
    
  __fpri_allocate_func = alloc_func;
  __fpri_callocate_func = calloc_func;
  __fpri_reallocate_func = realloc_func;
  __fpri_free_func = free_func;
}

void * fpri_malloc(size_t size) {
    return (*__fpri_allocate_func)(size);
}
void * fpri_realloc(void * ptr, size_t size) {
    return (*__fpri_reallocate_func)(ptr, size);
}
void * fpri_calloc(size_t num, size_t size) {
    return (*__fpri_callocate_func)(num, size);
}
void fpri_free(void * ptr) {
    (*__fpri_free_func)(ptr);
}
