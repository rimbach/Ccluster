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
    
//     clock_t start, start2;
//     double temp = 0;
//     start = clock();
    
//     printf("#connCmp_union_find.c, l 29: begin \n");
//     printf("#connCmp_union_find.c, l 30: nb of ccs: %d\n", connCmp_list_get_size(ccs) );
//     int cnt = 0;
    
    while (!connCmp_list_is_empty(ccs)){
        cctemp = connCmp_list_pop(ccs);
        /* optimization: the box is more probably connected to the last cc in the list */
//         cctemp = connCmp_list_pop_back(ccs);
        if (connCmp_are_8connected(cctemp, b)){
            
//             if (connCmp_list_is_empty(ccs)){
//                 printf("# connCmp_union_find.c, l40: connected to     %d cc in queue\n ", cnt);
//                 printf("                             number of boxes in cb: %d\n ", compBox_list_get_size(connCmp_boxesref(cb)));
//                 printf("                             number of boxes in cctemp: %d\n ", compBox_list_get_size(connCmp_boxesref(cctemp)));
//             } else {
//                 printf("# connCmp_union_find.c, l34: connected to NOT LAST cc in queue\n ");
//                 printf("                             number of boxes in cb: %d\n ", compBox_list_get_size(connCmp_boxesref(cb)));
//                 printf("                             number of boxes in cctemp: %d\n ", compBox_list_get_size(connCmp_boxesref(cctemp)));
//             }
            
//             start2 = clock();
            connCmp_merge_2_connCmp(cb, cctemp);
//             temp += ( (float) clock() - start2 )/CLOCKS_PER_SEC ;
            
            connCmp_clear(cctemp);
            ccluster_free(cctemp);
        }
        else 
            connCmp_list_push(ltemp, cctemp);
//         cnt ++;
    }
    
//     timeIn_merge_2_connCmp += temp;
//     timeIn_are_8connected  += ( ( (float) clock() - start )/CLOCKS_PER_SEC - temp );
    
    connCmp_list_push(ltemp, cb);
    connCmp_list_swap(ltemp, ccs);
    connCmp_list_clear(ltemp);
    
//     printf("#connCmp_union_find.c, l 67: end \n");
}
