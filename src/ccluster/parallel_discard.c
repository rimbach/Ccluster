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
    
    /* arg->status has been set to 1 by caller           */
    /* nb_thread_running has been incremented by caller; */
    
    parallel_bisect_arg_t * arg = (parallel_bisect_arg_t *) arg_ptr;
//     ccluster_bisect_connCmp( arg->res, arg->cc, arg->dis, arg->cache, arg->meta, arg->nbThreads);
    ccluster_bisect_connCmp_without_quadrisect( arg->res, arg->cc, arg->dis, arg->cache, arg->meta, arg->nbThreads);
    flint_cleanup();
    
    /*actualize datas for the scheduler */
    pthread_mutex_lock (&(arg->mutex));
    arg->status =2; /*is finished*/
    pthread_mutex_unlock (&(arg->mutex));
    pthread_mutex_lock (arg->mutex_nb_running);
//     (*(arg->nb_thread_running))--;
    (*(arg->nb_thread_running)) -= arg->nbThreads;
    pthread_mutex_unlock (arg->mutex_nb_running);
    
    return NULL;
    
}

void ccluster_parallel_bisect_connCmp_list( connCmp_list_ptr qMainLoop, connCmp_list_ptr discardedCcs,
                                            connCmp_list_ptr toBeBisected, cacheApp_t cache, metadatas_t meta){
    
    /* number of threads available for this task */
    slong nb_threads = metadatas_useNBThreads(meta);
    /* maximum number of threads that will be created here */
    slong nb_args = CCLUSTER_MIN(nb_threads,connCmp_list_get_size(toBeBisected));
    /* total number of boxes to be bisected */
    slong nb_boxes = 0;
    
    compBox_ptr btemp;
    compBox_list_t subBoxes;
    compBox_list_init(subBoxes);
    
    connCmp_list_iterator it = connCmp_list_begin(toBeBisected);
    while (it!=connCmp_list_end() ) {
        
        /*quadrisect boxes in it*/
        while (!connCmp_is_empty(connCmp_list_elmt(it))) {
            btemp = connCmp_pop(connCmp_list_elmt(it));
            subdBox_quadrisect( subBoxes, btemp );
            compBox_clear(btemp);
            ccluster_free(btemp);
        }
        compBox_list_swap(subBoxes, connCmp_boxesref(connCmp_list_elmt(it)));
        
        nb_boxes += connCmp_nb_boxes(connCmp_list_elmt(it));
        it = connCmp_list_next(it);
    }
    /*number of boxes that should be given to a thread*/
    slong nb_boxes_per_thread = (slong) ceil( ((double) nb_boxes)/((double) nb_threads ));
    /*sort connected components in toBeBisected by decreasing number of boxes*/
    connCmp_list_t ltemp;
    connCmp_list_init(ltemp);
    while (!connCmp_list_is_empty(toBeBisected))
        connCmp_list_push(ltemp, connCmp_list_pop(toBeBisected));
    connCmp_list_swap(toBeBisected, ltemp);
    connCmp_list_clear(ltemp);
    
    /*create the arguments*/
    parallel_bisect_arg_t * args = (parallel_bisect_arg_t *) malloc ( sizeof(parallel_bisect_arg_t) * nb_args );
    pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
    int nb_thread_running = 0;
    pthread_mutex_t mutex_nb_running;
    pthread_mutex_init ( &mutex_nb_running, NULL);
    
//     printf("{ ");
    
    /*initialize args, lists  */
    for (int i = 0; i< (int) nb_args; i++){
        connCmp_list_init(args[i].res);
        connCmp_list_init(args[i].dis);
        args[i].cache = (cacheApp_ptr) cache;
        args[i].meta = (metadatas_ptr) meta;
        args[i].status=0;
        pthread_mutex_init ( &(args[i].mutex), NULL);
        args[i].nb_thread_running = &nb_thread_running;
        args[i].mutex_nb_running  = &mutex_nb_running;
    }
       
    /*main loop: empty toBeBisected*/
    while (!connCmp_list_is_empty(toBeBisected)) {
        
        if (nb_thread_running<nb_threads) {
            int thread = 0;
            /* find an available thread */
            while( (thread < nb_args) && (args[thread].status==1) ) thread++;
            if (args[thread].status==2) { /* join the thread */
//                 printf("----join: %d\n", thread);
                pthread_join(threads[thread], NULL);
                while (!connCmp_list_is_empty(args[thread].res))
                    connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(args[thread].res));
                while (!connCmp_list_is_empty(args[thread].dis))
                    connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(args[thread].dis));
                
                connCmp_clear(args[thread].cc);
                free(args[thread].cc);
                connCmp_list_clear(args[thread].res);
                connCmp_list_clear(args[thread].dis);
            }
            /* actualize arg and nbrunning */
            args[thread].cc = connCmp_list_pop(toBeBisected);
            args[thread].status=1;
            args[thread].nbThreads = (slong) ceil( ((double) connCmp_nb_boxes(args[thread].cc)) / ((double) nb_boxes_per_thread) ); 
//             printf("nb_boxes: %d\n", connCmp_nb_boxes(args[thread].cc));
//             printf("nb_boxes_per_thread: %ld\n", nb_boxes_per_thread);
//             printf("div: %f\n", ((double) connCmp_nb_boxes(args[thread].cc)) / ((double) nb_boxes_per_thread));
//             printf("args[thread].nbThreads: %ld\n", args[thread].nbThreads);
//             printf("nb dispos: %ld\n", nb_threads - nb_thread_running);
            args[thread].nbThreads = CCLUSTER_MIN(args[thread].nbThreads, nb_threads - nb_thread_running);
//             printf("args[thread].nbThreads: %ld\n", args[thread].nbThreads);
            args[thread].nbThreads = CCLUSTER_MAX(args[thread].nbThreads, 1);
            
            pthread_mutex_lock (&(mutex_nb_running));
            nb_thread_running += args[thread].nbThreads;
            pthread_mutex_unlock (&(mutex_nb_running));
            /* create the thread */
//             printf("----create: %d\n", thread);
//             printf("%d:%ld ", connCmp_nb_boxes(args[thread].cc),args[thread].nbThreads);
            pthread_create(&threads[thread], NULL, _parallel_bisect_worker, &args[thread]);
            
        }
    }
//     printf("} \n");
    
    /* join threads still running */
    for (int i = 0; i< (int) nb_args; i++){
        if (args[i].status>0) { /* join the thread */
//             printf("----join: %d\n", i);
            pthread_join(threads[i], NULL);
            while (!connCmp_list_is_empty(args[i].res))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(args[i].res));
            while (!connCmp_list_is_empty(args[i].dis))
                connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(args[i].dis));
            connCmp_clear(args[i].cc);
            free(args[i].cc);
        }
        pthread_mutex_destroy( &(args[i].mutex) );
        connCmp_list_clear(args[i].res);
        connCmp_list_clear(args[i].dis);
    }
    
    compBox_list_clear(subBoxes);
    
    free(args);
    free(threads);
}

void ccluster_bisect_connCmp_without_quadrisect( connCmp_list_t dest, 
                                                 connCmp_t cc, 
                                                 connCmp_list_t discardedCcs, 
                                                 cacheApp_t cache, 
                                                 metadatas_t meta, 
                                                 slong nbThreads){
    
    slong prec = connCmp_appPr(cc);
    compBox_list_t subBoxes;
    connCmp_list_t ltemp;
    compBox_list_init(subBoxes);
    connCmp_list_init(ltemp);
    
    compBox_ptr btemp;
    connCmp_ptr ctemp;
    
//     while (!connCmp_is_empty(cc)) {
//         btemp = connCmp_pop(cc);
//         subdBox_quadrisect( subBoxes, btemp );
//         compBox_clear(btemp);
//         ccluster_free(btemp);
//     }
    compBox_list_swap(subBoxes, connCmp_boxesref(cc));

// #ifdef CCLUSTER_HAVE_PTHREAD
    if (nbThreads>1) {
//         printf("--ccluster_parallel_bisect_connCmp: nb threads: %d \n", (int) nbThreads );
        prec = ccluster_parallel_discard_compBox_list( subBoxes, cache, prec, meta, nbThreads);
    }
    else
        prec = ccluster_discard_compBox_list( subBoxes, cache, prec, meta);
// #else
//     prec = ccluster_discard_compBox_list( subBoxes, cache, prec, meta);
// #endif
    
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    int specialFlag = 1;
    if (connCmp_list_get_size(ltemp) == 1)
        specialFlag = 0;
    
    slong nprec; 
    if (prec == connCmp_appPrref(cc)) {
        nprec = CCLUSTER_MAX(prec/2,CCLUSTER_DEFAULT_PREC);
//         printf("decrease precision\n");
    }
    else 
        nprec = prec;
    
    
    while (!connCmp_list_is_empty(ltemp)){
        ctemp = connCmp_list_pop(ltemp);
        
        if (connCmp_intersection_is_not_empty(ctemp, metadatas_initBref(meta))){
            connCmp_appPrref(ctemp) = nprec;
            if (specialFlag)
                connCmp_initiali_nwSpd(ctemp);
            else {
                connCmp_initiali_nwSpd_connCmp(ctemp, cc);
                connCmp_decrease_nwSpd(ctemp);
                /* copy the number of sols */
                connCmp_nSolsref(ctemp) = connCmp_nSolsref(cc);
                /* test */
                connCmp_isSep(ctemp) = connCmp_isSep(cc);
                /*end test */
            }
            connCmp_list_push(dest, ctemp);
        }
        else {
            connCmp_appPrref(ctemp) = prec;
            connCmp_list_push(discardedCcs, ctemp);
        }
    }
    
    compBox_list_clear(subBoxes);
    connCmp_list_clear(ltemp);
}



// void ccluster_parallel_bisect_connCmp_list( connCmp_list_ptr qMainLoop, connCmp_list_ptr discardedCcs,
//                                             connCmp_list_ptr toBeBisected, cacheApp_t cache, metadatas_t meta){
//     
// //     printf("--ccluster_parallel_bisect_connCmp_list: nb connCmp: %d, nb threads: %d \n", (int) connCmp_list_get_size(toBeBisected), (int) metadatas_useNBThreads(meta) );
//     slong nb_threads = metadatas_useNBThreads(meta);
//     
// //     /*for test*/
// //     if (metadatas_useNBThreads(meta)==1023) {
// //         nb_threads = connCmp_list_get_size(toBeBisected);
// //     }
// //     /*end test*/
//     slong nb_threads_by_task = nb_threads;
//     if (connCmp_list_get_size(toBeBisected)>0)
//         nb_threads_by_task = (slong) (nb_threads / connCmp_list_get_size(toBeBisected));
//     nb_threads_by_task = CCLUSTER_MAX( nb_threads_by_task, 1 );   
//     slong nb_args = CCLUSTER_MIN(nb_threads,connCmp_list_get_size(toBeBisected));
//     parallel_bisect_arg_t * args = (parallel_bisect_arg_t *) malloc ( sizeof(parallel_bisect_arg_t) * nb_args );
//     pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
//     
//     int nb_thread_running = 0;
//     pthread_mutex_t mutex_nb_running;
//     pthread_mutex_init ( &mutex_nb_running, NULL);
//     
//     printf("{ ");
//     
//     /*initialize args, lists and create nb_args threads */
//     for (int i = 0; i< (int) nb_args; i++){
//         args[i].cc = connCmp_list_pop(toBeBisected);
//         connCmp_list_init(args[i].res);
//         connCmp_list_init(args[i].dis);
//         args[i].cache = (cacheApp_ptr) cache;
//         args[i].meta = (metadatas_ptr) meta;
//         args[i].nbThreads = nb_threads_by_task;
//         args[i].status=1;
//         pthread_mutex_init ( &(args[i].mutex), NULL);
//         args[i].nb_thread_running = &nb_thread_running;
//         args[i].mutex_nb_running  = &mutex_nb_running;
//         
//         pthread_mutex_lock (&(mutex_nb_running));
//         nb_thread_running ++;
//         pthread_mutex_unlock (&(mutex_nb_running));
//         /* create the thread */
// //         printf("----create: %d\n", i);
//         printf("%d:%ld ", 4*connCmp_nb_boxes(args[i].cc),nb_threads_by_task);
//         pthread_create(&threads[i], NULL, _parallel_bisect_worker, &args[i]);
//     }
//     
//     /*main loop: empty toBeBisected*/
//     while (!connCmp_list_is_empty(toBeBisected)) {
//         
//         if (nb_thread_running<nb_args) {
//             int thread = 0;
//             /* find an available thread */
//             while( (thread < nb_args) && (args[thread].status==1) ) thread++;
//             if (args[thread].status==2) { /* join the thread */
// //                 printf("----join: %d\n", thread);
//                 pthread_join(threads[thread], NULL);
//                 while (!connCmp_list_is_empty(args[thread].res))
//                     connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(args[thread].res));
//                 while (!connCmp_list_is_empty(args[thread].dis))
//                     connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(args[thread].dis));
//                 
//                 connCmp_clear(args[thread].cc);
//                 free(args[thread].cc);
//                 connCmp_list_clear(args[thread].res);
//                 connCmp_list_clear(args[thread].dis);
//             }
//             
//             /* actualize arg and nbrunning */
//             args[thread].cc = connCmp_list_pop(toBeBisected);
//             args[thread].status=1;
//             pthread_mutex_lock (&(mutex_nb_running));
//             nb_thread_running ++;
//             pthread_mutex_unlock (&(mutex_nb_running));
//             /* create the thread */
// //             printf("----create: %d\n", thread);
//             printf("%d:%ld ", 4*connCmp_nb_boxes(args[thread].cc),nb_threads_by_task);
//             pthread_create(&threads[thread], NULL, _parallel_bisect_worker, &args[thread]);
//         }
//         
//     }
//     
//     printf("} \n");
//     
//     /* join threads still running */
//     for (int i = 0; i< (int) nb_args; i++){
//         if (args[i].status>0) { /* join the thread */
// //             printf("----join: %d\n", i);
//             pthread_join(threads[i], NULL);
//             while (!connCmp_list_is_empty(args[i].res))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(args[i].res));
//             while (!connCmp_list_is_empty(args[i].dis))
//                 connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(args[i].dis));
//             connCmp_clear(args[i].cc);
//             free(args[i].cc);
//         }
//         pthread_mutex_destroy( &(args[i].mutex) );
//         connCmp_list_clear(args[i].res);
//         connCmp_list_clear(args[i].dis);
//     }
//     
//     free(args);
//     free(threads);
// }

void * _parallel_discard_list_worker( void * arg_ptr ){
    
    parallel_discard_list_arg_t * arg = (parallel_discard_list_arg_t *) arg_ptr;
    
    /* arg->status has been set to 1 by caller           */
    /* nb_thread_running has been incremented by caller; */
    
    arg->prec = ccluster_discard_compBox_list( arg->boxes, arg->cache, arg->prec, arg->meta);
    
    flint_cleanup();
    
    /*actualize datas for the scheduler */
//     pthread_mutex_lock (&(arg->mutex));
//     arg->status =2; /*is finished*/
//     pthread_mutex_unlock (&(arg->mutex));
//     pthread_mutex_lock (arg->mutex_nb_running);
//     (*(arg->nb_thread_running))--;
//     pthread_mutex_unlock (arg->mutex_nb_running);
    
    return NULL;
    
}



slong ccluster_parallel_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
                                        slong prec, metadatas_t meta, slong nbThreads){
    
//     slong nb_threads = metadatas_useNBThreads(meta);
    slong nb_threads = nbThreads;
    
//     /*for test*/
//     if (metadatas_useNBThreads(meta)==1023) {
//         nb_threads = compBox_list_get_size(boxes);
//     }
//     /*end test*/
    
    slong precres = prec;
    slong nb_args = CCLUSTER_MIN(nb_threads,compBox_list_get_size(boxes));
    parallel_discard_list_arg_t * args = (parallel_discard_list_arg_t *) malloc ( sizeof(parallel_discard_list_arg_t) * nb_args );
    pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
    
    int nb_boxes_by_thread = ((int) compBox_list_get_size(boxes))/((int) nb_args);
    
    for (int i = 0; i< (int) nb_args; i++) {
        
        args[i].prec = precres;
        compBox_list_init(args[i].boxes);
        args[i].cache = (cacheApp_ptr) cache;
        args[i].meta = (metadatas_ptr) meta;
        /*splits boxes in nbthreads lists*/
        int j=0;
        while ( (!compBox_list_is_empty(boxes))&& ((j<nb_boxes_by_thread)||(i==(nb_args-1)) ) ) {
            compBox_list_push(args[i].boxes, compBox_list_pop(boxes));
            j++;
        }
        /* create the thread */
        pthread_create(&threads[i], NULL, _parallel_discard_list_worker, &args[i]);
    }
    
    for(int i = 0; i< (int) nb_args; i++) {
        pthread_join(threads[i], NULL);
        if (args[i].prec > precres)
            precres = args[i].prec;
        /* fill boxes */
        while (!compBox_list_is_empty(args[i].boxes))
            compBox_list_push(boxes, compBox_list_pop(args[i].boxes));
        compBox_list_clear(args[i].boxes);
    }
    
    free(args);
    free(threads);
    
    return precres;
}

/* DEPRECATED */

// void ccluster_parallel_bisect_connCmp_list( connCmp_list_ptr qMainLoop, connCmp_list_ptr discardedCcs,
//                                             connCmp_list_ptr toBeBisected, cacheApp_t cache, metadatas_t meta){
//     
// //     printf("--ccluster_parallel_bisect_connCmp_list: nb connCmp: %d, nb threads: %d \n", (int) connCmp_list_get_size(toBeBisected), (int) metadatas_useNBThreads(meta) );
//     slong nb_threads = metadatas_useNBThreads(meta);
//     slong nb_threads_by_task = nb_threads;
//     if (connCmp_list_get_size(toBeBisected)>0)
//         nb_threads_by_task = (slong) (nb_threads / connCmp_list_get_size(toBeBisected));
//     slong nb_args = CCLUSTER_MIN(nb_threads,connCmp_list_get_size(toBeBisected));
//     parallel_bisect_arg_t * args = (parallel_bisect_arg_t *) malloc ( sizeof(parallel_bisect_arg_t) * nb_args );
//     pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
//     
//     
//     slong nb_pop;
//     while (!connCmp_list_is_empty(toBeBisected)){
//         
//         nb_pop = CCLUSTER_MIN( connCmp_list_get_size(toBeBisected), nb_args );
//         for ( int i=0; i< (int) nb_pop; i++) {
//             args[i].cc = connCmp_list_pop(toBeBisected);
//             connCmp_list_init(args[i].res);
//             connCmp_list_init(args[i].dis);
//             args[i].meta = (metadatas_ptr) meta;
//             args[i].cache = (cacheApp_ptr) cache;
//             args[i].nbThreads = nb_threads_by_task;
//             pthread_create(&threads[i], NULL, _parallel_bisect_worker, &args[i]);
//         }
//         for ( int i=0; i< (int) nb_pop; i++) {
//             pthread_join(threads[i], NULL);
//             while (!connCmp_list_is_empty(args[i].res))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(args[i].res));
//             while (!connCmp_list_is_empty(args[i].dis))
//                 connCmp_list_insert_sorted(discardedCcs, connCmp_list_pop(args[i].dis));
//             
//             connCmp_clear(args[i].cc);
//             free(args[i].cc);
//             connCmp_list_clear(args[i].res);
//             connCmp_list_clear(args[i].dis);
//         }
//     }
//     free(args);
//     free(threads);
// }

// slong ccluster_parallel_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
//                                         slong prec, metadatas_t meta, slong nbThreads){
//     
// //     slong nb_threads = metadatas_useNBThreads(meta);
//     slong nb_threads = nbThreads;
//     slong precres = prec;
//     slong nb_args = CCLUSTER_MIN(nb_threads,compBox_list_get_size(boxes));
//     parallel_discard_list_arg_t * args = (parallel_discard_list_arg_t *) malloc ( sizeof(parallel_discard_list_arg_t) * nb_args );
//     pthread_t * threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_args);
//     
//     compBox_list_t ltemp;
//     compBox_list_init(ltemp);
//     
//     int nb_thread_running = 0;
//     pthread_mutex_t mutex_nb_running;
//     pthread_mutex_init ( &mutex_nb_running, NULL);
//     
// //     printf("----ccluster_parallel_discard_compBox_list: nb_boxes: %d, nb_args: %d\n", (int) compBox_list_get_size(boxes), (int) nb_args);
//     
//     /*initialize args, lists and create nb_args threads */
//     for (int i = 0; i< (int) nb_args; i++){
//         args[i].cache = (cacheApp_ptr) cache;
//         args[i].meta = (metadatas_ptr) meta;
// //         args[i].status = 0;
//         pthread_mutex_init ( &(args[i].mutex), NULL);
//         args[i].nb_thread_running = &nb_thread_running;
//         args[i].mutex_nb_running  = &mutex_nb_running;
//         compBox_list_init(args[i].boxes);
//         
//         compBox_list_push(args[i].boxes, compBox_list_pop(boxes));
//         args[i].prec = precres;
//         args[i].status=1;
//         pthread_mutex_lock (&(mutex_nb_running));
//         nb_thread_running ++;
//         pthread_mutex_unlock (&(mutex_nb_running));
//         /* create the thread */
// //         printf("----create: %d\n", i);
//         pthread_create(&threads[i], NULL, _parallel_discard_list_worker, &args[i]);
//     }
//     
//     /*main loop: empty boxes*/
//     while (!compBox_list_is_empty(boxes)) {
//         
//         if (nb_thread_running<nb_args) {
//             int thread = 0;
//             /* find an available thread */
//             while( (thread < nb_args) && (args[thread].status==1) ) thread++;
//             if (args[thread].status==2) { /* join the thread */
// //                 printf("----join: %d\n", thread);
//                 pthread_join(threads[thread], NULL);
//                 if (args[thread].prec > precres)
//                     precres = args[thread].prec;
//                 /* fill boxes */
//                 while (!compBox_list_is_empty(args[thread].boxes))
//                     compBox_list_push(ltemp, compBox_list_pop(args[thread].boxes));
// //                     compBox_list_insert_sorted(ltemp, compBox_list_pop(&lists[thread]));
//             }
//             /* actualize arg and nbrunning */
//             compBox_list_push(args[thread].boxes, compBox_list_pop(boxes));
//             args[thread].prec = precres;
//             args[thread].status=1;
//             pthread_mutex_lock (&(mutex_nb_running));
//             nb_thread_running ++;
//             pthread_mutex_unlock (&(mutex_nb_running));
//             /* create the thread */
// //             printf("----create: %d\n", thread);
//             pthread_create(&threads[thread], NULL, _parallel_discard_list_worker, &args[thread]);
//         }
//         
//     }
//     /* join threads still running */
//     for (int i = 0; i< (int) nb_args; i++){
//         if (args[i].status>0) { /* join the thread */
// //             printf("----join: %d\n", i);
//             pthread_join(threads[i], NULL);
//             if (args[i].prec > precres)
//                 precres = args[i].prec;
//             /* fill boxes */
//             while (!compBox_list_is_empty(args[i].boxes))
//                 compBox_list_push(ltemp, compBox_list_pop(args[i].boxes));
// //                 compBox_list_insert_sorted(ltemp, compBox_list_pop(&lists[i]));
//         }
//         pthread_mutex_destroy( &(args[i].mutex) );
//         compBox_list_clear(args[i].boxes);
//     }
//        
//     compBox_list_swap(boxes, ltemp);
//     
//     free(args);
//     free(threads);
//     compBox_list_clear(ltemp);
//     pthread_mutex_destroy( &mutex_nb_running );
//     
//     return precres;
// }


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
