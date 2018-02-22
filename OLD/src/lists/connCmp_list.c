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

#include "lists/connCmp_list.h"

void connCmp_list_init_forjulia(connCmp_list_t l){
    gen_list_init(l, connCmp_clear_for_list);
}

void connCmp_list_clear_forjulia(connCmp_list_t l) {
    gen_list_clear(l);
}

void connCmp_list_push_forjulia(connCmp_list_t l, connCmp_ptr b){
    gen_list_push(l, b);
}

connCmp_ptr connCmp_list_pop_forjulia(connCmp_list_t l){
    return gen_list_pop(l);
}

int connCmp_list_is_empty_forjulia(const connCmp_list_t l){
    return gen_list_is_empty(l);
}