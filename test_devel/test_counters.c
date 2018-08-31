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
#include "metadatas/counters.h"

int main() {
    
    
    counters_t st;
    counters_init(st);
    
    printf("size of table: %u, size allocated: %u\n", st->size, st->size_allocated);
    counters_add_discarded(st, 0);
    printf("size of table: %u, size allocated: %u, nbDiscarded: %u\n", st->size, st->size_allocated, st->table[st->size-1].nbDiscarded);
    counters_add_discarded(st, 0);
    printf("size of table: %u, size allocated: %u, nbDiscarded: %u\n", st->size, st->size_allocated, st->table[st->size-1].nbDiscarded);
    counters_add_discarded(st, 1000);
    printf("size of table: %u, size allocated: %u, nbDiscarded: %u\n", st->size, st->size_allocated, st->table[st->size-1].nbDiscarded);
    
    counters_clear(st);
    
    return 0;
    
}