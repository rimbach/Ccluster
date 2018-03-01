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

#include "ccluster/ccluster.h"

#ifdef CCLUSTER_HAVE_PTHREAD

#include <pthread.h>

void * _parallel_bisect_worker( void * arg_ptr ){
    
    parallel_bisect_arg_t * arg = (parallel_bisect_arg_t *) arg_ptr;
    ccluster_bisect_connCmp( arg->res, arg->cc, arg->dis, arg->cache, arg->meta, arg->nbThreads);
    flint_cleanup();
    return NULL;
    
}

void ccluster_parallel_bisect_connCmp_list( connCmp_list_ptr qMainLoop, connCmp_list_ptr discardedCcs,
                                            connCmp_list_ptr toBeBisected, cacheApp_t cache, metadatas_t meta){
    
    printf("--ccluster_parallel_bisect_connCmp_list: nb connCmp: %d, nb threads: %d \n", (int) connCmp_list_get_size(toBeBisected), (int) metadatas_useNBThreads(meta) );
    slong nb_threads = metadatas_useNBThreads(meta);
    slong nb_threads_by_task;
    if (connCmp_list_get_size(toBeBisected)>0)
        nb_threads_by_task = (slong) (nb_threads / connCmp_list_get_size(toBeBisected));
    else
        nb_threads_by_task = nb_threads;
    slong nb_args = CCLUSTER_MIN(nb_threads,connCmp_list_get_size(toBeBisected));
    parallel_bisect_arg_t * args = (parallel_bisect_arg_t *) malloc ( sizeof(parallel_bisect_arg_t) * nb_args );
    pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
    
    connCmp_list_ptr ress = (connCmp_list_ptr) malloc (sizeof(connCmp_list) * nb_args);
    connCmp_list_ptr diss = (connCmp_list_ptr) malloc (sizeof(connCmp_list) * nb_args);
    
    slong nb_pop;
    connCmp_ptr ccur;
    while (!connCmp_list_is_empty(toBeBisected)){
        
        nb_pop = CCLUSTER_MIN( connCmp_list_get_size(toBeBisected), nb_args );
        for ( int i=0; i< (int) nb_pop; i++) {
            ccur = connCmp_list_pop(toBeBisected);
            connCmp_list_init(&ress[i]);
            connCmp_list_init(&diss[i]);
            args[i].res = &ress[i];
            args[i].dis = &diss[i];
            args[i].cc  = ccur;
            args[i].meta = (metadatas_ptr) meta;
            args[i].cache = (cacheApp_ptr) cache;
            args[i].nbThreads = nb_threads_by_task;
            pthread_create(&threads[i], NULL, _parallel_bisect_worker, &args[i]);
        }
        for ( int i=0; i< (int) nb_pop; i++) {
            pthread_join(threads[i], NULL);
            while (!connCmp_list_is_empty(&ress[i]))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(&ress[i]));
            while (!connCmp_list_is_empty(&diss[i]))
                connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(&diss[i]));
            connCmp_clear(args[i].cc);
            free(args[i].cc);
            connCmp_list_clear(&ress[i]);
            connCmp_list_clear(&diss[i]);
        }
    }
    free(args);
    free(threads);
    free(ress);
    free(diss);
}

void * _parallel_discard_list_worker( void * arg_ptr ){
    
    parallel_discard_list_arg_t * arg = (parallel_discard_list_arg_t *) arg_ptr;
    
    tstar_res res;
    res.appPrec = arg->prec;
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    compBox_ptr btemp;
    
    while (!compBox_list_is_empty(arg->boxes)){
        btemp = compBox_list_pop(arg->boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref(arg->meta));
        res = tstar_interface( arg->cache, bdisk, compBox_get_nbMSol(btemp), 1, res.appPrec, depth, arg->meta);
        
        if (res.nbOfSol==0) {
            metadatas_lock(arg->meta); 
            metadatas_add_discarded( arg->meta, depth);
            metadatas_unlock(arg->meta);
            compBox_clear(btemp);
            free(btemp);
        }
        
        else{
            if (res.nbOfSol>0)
                    btemp->nbMSol = res.nbOfSol;
                compBox_list_push(ltemp, btemp);
        }
    }
    compBox_list_swap(arg->boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
    arg->prec = res.appPrec;
    
    flint_cleanup();
    return NULL;
    
}

slong ccluster_parallel_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
                                        slong prec, metadatas_t meta, slong nbThreads){
    
//     slong nb_threads = metadatas_useNBThreads(meta);
    slong nb_threads = nbThreads;
    slong precres = prec;
    slong nb_args = CCLUSTER_MIN(nb_threads,compBox_list_get_size(boxes));
    parallel_discard_list_arg_t * args = (parallel_discard_list_arg_t *) malloc ( sizeof(parallel_discard_list_arg_t) * nb_args );
    pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
//     metadatas_ptr metas = meta;
    
    /*splits boxes in nbthreads lists*/
    compBox_list_ptr lists = (compBox_list_ptr) malloc (sizeof(compBox_list)* nb_args);
    int nb_boxes_by_thread = ((int) compBox_list_get_size(boxes))/((int) nb_args);
    
    for (int i = 0; i< (int) nb_args; i++) {
        compBox_list_init(&lists[i]);
        int j=0;
        while ( (!compBox_list_is_empty(boxes))&& ((j<nb_boxes_by_thread)||(i==(nb_args-1)) ) ) {
            compBox_list_push(&lists[i], compBox_list_pop(boxes));
            j++;
        }
        /* create the thread */
        args[i].prec = precres;
        args[i].boxes  = (compBox_list_ptr) &lists[i];
        args[i].cache = (cacheApp_ptr) cache;
        args[i].meta = (metadatas_ptr) meta;
        pthread_create(&threads[i], NULL, _parallel_discard_list_worker, &args[i]);
    }
//     printf("size of boxes: %d\n", (int) compBox_list_get_size(boxes));
    
//     for(int i = 0; i< (int) nb_args; i++) {
//         args[i].prec = precres;
//         args[i].boxes  = (compBox_list_ptr) &lists[i];
//         args[i].cache = (cacheApp_ptr) cache;
//         args[i].meta = (metadatas_ptr) meta;
//         pthread_create(&threads[i], NULL, _parallel_discard_list_worker, &args[i]);
//     }
    
    for(int i = 0; i< (int) nb_args; i++) {
        pthread_join(threads[i], NULL);
        if (args[i].prec > precres)
            precres = args[i].prec;
        /* fill boxes */
        while (!compBox_list_is_empty(&lists[i]))
            compBox_list_push(boxes, compBox_list_pop(&lists[i]));
        compBox_list_clear(&lists[i]);
    }
    
//     /*fill boxes*/
//     for (int i = 0; i< (int) nb_args; i++) {
//         while (!compBox_list_is_empty(&lists[i]))
//             compBox_list_push(boxes, compBox_list_pop(&lists[i]));
//         compBox_list_clear(&lists[i]);
//     }
    
    free(args);
    free(threads);
    free(lists);
    
    return precres;
}

// void * _parallel_discard_worker( void * arg_ptr ){
//     
//     parallel_discard_arg_t * arg = (parallel_discard_arg_t *) arg_ptr;
//     
// //     printf("begin worker, nbsols: %d\n", arg->nbsol);
//     
//     tstar_res res;
//     res.appPrec = arg->prec;
//     
//     slong depth;
//     compDsk_t bdisk;
//     compDsk_init(bdisk);
//     compBox_get_containing_dsk(bdisk, arg->box);
//     depth = compDsk_getDepth(bdisk, metadatas_initBref( arg->meta));
//     
//     res = tstar_interface( arg->cache, bdisk, compBox_get_nbMSol(arg->box), 1, res.appPrec, depth, arg->meta);
//     
//     arg->nbsol = res.nbOfSol;
//     arg->prec  = res.appPrec;
//     compDsk_clear(bdisk);
//     
// //     printf("end   worker, nbsols: %d\n", arg->nbsol);
//     return NULL;
// }

// slong ccluster_parallel_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
//                                         slong prec, metadatas_t meta, slong nb_threads){
//     
//     slong nb_args = CCLUSTER_MIN(nb_threads,compBox_list_get_size(boxes));
//     parallel_discard_arg_t * args = (parallel_discard_arg_t *) malloc ( sizeof(parallel_discard_arg_t) * nb_args );
//     pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
//     metadatas_ptr metas = (metadatas_ptr) malloc (sizeof(metadatas) * nb_args );
//     
//     slong nb_pop;
//     slong precres = prec;
//     
//     slong depth;
//     compBox_list_t ltemp;
//     compBox_list_init(ltemp);
//     compBox_ptr btemp;
//     
//     while (!compBox_list_is_empty(boxes)){
//         nb_pop = CCLUSTER_MIN(nb_threads,compBox_list_get_size(boxes));
//         
//         /* pop box and create threads */
//         for (int i=0; i<(int) nb_pop; i++) {
//             
//             btemp = compBox_list_pop(boxes);
//             
//             metadatas_init( metas + i, metadatas_initBref(meta), metadatas_stratref(meta) , metadatas_getVerbo(meta));
//             args[i].prec = precres;
//             args[i].box  = btemp;
//             args[i].cache = (cacheApp_ptr) cache;
//             args[i].meta = (metadatas_ptr) metas + i;
//             
//             pthread_create(&threads[i], NULL, _parallel_discard_worker, &args[i]);
//         }
//         
//         /* join threads and treat box */
//         for (int i=0; i<(int) nb_pop; i++) {
//             pthread_join(threads[i], NULL);
//             
//             btemp = args[i].box;
//             if (args[i].prec > precres)
//                 precres = args[i].prec;
//             metadatas_join(meta, metas + i);
//             metadatas_clear(metas + i);
//             
//             /*compute depth*/
//             depth = compBox_getDepth(args[i].box, metadatas_initBref(meta));
//             
//             if (args[i].nbsol == 0){
//                 metadatas_add_discarded( meta, depth);
//                 compBox_clear(btemp);
//                 free(btemp);
//             }
//             else {
//                 if (args[i].nbsol>0)
//                     btemp->nbMSol = args[i].nbsol;
//                 compBox_list_push(ltemp, btemp);
//             }
//         }
//         
//     }
//     
//     compBox_list_swap(boxes, ltemp);
//     compBox_list_clear(ltemp);
//     
//     free(metas);
//     free(args);
//     free(threads);
//     
//     return precres;
// }
    
#endif /* CCLUSTER_HAVE_PTHREAD */    