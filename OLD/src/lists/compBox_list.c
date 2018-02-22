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

#include "lists/compBox_list.h"

void compBox_list_init_forjulia(compBox_list_t l){
    gen_list_init(l, compBox_clear_for_list);
}

void compBox_list_clear_forjulia(compBox_list_t l) {
    gen_list_clear(l);
}

void compBox_list_push_forjulia(compBox_list_t l, compBox_ptr b){
    gen_list_push(l, b);
}

void compBox_list_pop_forjulia(compBox_t dest, compBox_list_t l){
    compBox_set(dest, compBox_list_pop( l ));
}

int compBox_list_is_empty_forjulia(const compBox_list_t l){
    return gen_list_is_empty(l);
}