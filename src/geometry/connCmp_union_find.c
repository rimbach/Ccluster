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

#include "geometry/connCmp_union_find.h"

void connCmp_union_compBox( connCmp_list_t ccs, compBox_t b){
    
    connCmp_ptr cb;
    cb = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(cb, b);
    
    connCmp_list_t ltemp;
    connCmp_list_init(ltemp);
    
    connCmp_ptr cctemp;
    
    while (!connCmp_list_is_empty(ccs)){
        cctemp = connCmp_list_pop(ccs);
        if (connCmp_are_8connected(cctemp, b)){
            connCmp_merge_2_connCmp(cb, cctemp);
            connCmp_clear(cctemp);
            ccluster_free(cctemp);
        }
        else 
            connCmp_list_push(ltemp, cctemp);
    }
    connCmp_list_push(ltemp, cb);
    connCmp_list_swap(ltemp, ccs);
    connCmp_list_clear(ltemp);
}
