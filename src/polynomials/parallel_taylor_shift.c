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

#include "app_rat_poly.h"
#include <pthread.h>


#define CCLUSTER_MAX(a,b) (a>b ? a : b)
#define CCLUSTER_MIN(a,b) (a<b ? a : b)

typedef struct {
    slong size;
    slong prec;
    compApp_ptr le_table; /* table of size compApp_t */
    compApp_ptr up_table; /* table of size compApp_t */
    slong le_inf;
    slong le_sup;
    slong up_inf;
    slong up_sup;
    compApp * c;
} taylor_shift_parallel_arg_t;

/*// void * _compApp_poly_parallel_taylor_worker( void * arg_ptr ){
//     
//     taylor_shift_parallel_arg_t arg = *((taylor_shift_parallel_arg_t *) arg_ptr);
//     slong i, j;
//     compApp_t temp;
//     compApp_init(temp);
//     
//     for ( i=arg.le_sup; i>= arg.le_inf; i-- ) {
// //         printf("i: %d ", i); compApp_printd( arg.le_table + i, arg.prec ); printf("\n");
//         j = arg.up_inf;
//         if (j<=i){
// //             j = arg.up_inf;
//             compApp_mul( arg.up_table + j, arg.up_table + j, arg.c, arg.prec);
//             compApp_add( arg.up_table + j, arg.up_table + j, arg.le_table + i, arg.prec);
// //             printf("(%d,%d), ", i,j);
//             for ( j=arg.up_inf+1; j<= CCLUSTER_MIN(arg.up_sup, i); j++ ) {
//                 compApp_mul( arg.up_table + j, arg.up_table + j, arg.c, arg.prec);
//                 compApp_add( arg.up_table + j, arg.up_table + j, arg.up_table + (j-1), arg.prec);
// //                 printf("(%d,%d), ", i,j);
//     //             compApp_set(temp, arg.up_table + (j-1));
//     //             compApp_addmul( temp, arg.up_table + j, arg.c, arg.prec);
//     //             compApp_set(arg.up_table + j, temp);
//     //             printf("j: %d", j); compApp_printd( arg.up_table + j, arg.prec ); printf("\n");
//             }
//     //         compApp_set( arg.le_table + i, arg.up_table + CCLUSTER_MIN(arg.up_sup, i) );
//             compApp_set( arg.le_table + i, arg.up_table + j-1 );
//         }
//     }
// //     printf("\n");
//     compApp_clear(temp);
//     return NULL;
// }
*/
void * _compApp_poly_parallel_taylor_worker( void * arg_ptr ){
    
    taylor_shift_parallel_arg_t arg = *((taylor_shift_parallel_arg_t *) arg_ptr);
    slong i, j;
    
    for ( i=arg.le_inf; i<= arg.le_sup; i++ ) {
        j = arg.up_inf;
        if (j<= ((arg.size-1)-i)){
            compApp_mul( arg.up_table + j, arg.up_table + j, arg.c, arg.prec);
            compApp_add( arg.up_table + j, arg.up_table + j, arg.le_table + i, arg.prec);
/*             printf("(%d,%d), ", i,j);*/
            for ( j=arg.up_inf+1; j<= CCLUSTER_MIN(arg.up_sup, (arg.size-1)-i); j++ ) {
                compApp_mul( arg.up_table + j, arg.up_table + j, arg.c, arg.prec);
                compApp_add( arg.up_table + j, arg.up_table + j, arg.up_table + (j-1), arg.prec);
/*                 printf("(%d,%d), ", i,j);
        //         compApp_set(temp, arg.up_table + (j-1));
        //         compApp_addmul( temp, arg.up_table + j, arg.c, arg.prec);
        //         compApp_set(arg.up_table + j, temp);
        //         printf("j: %d", j); compApp_printd( arg.up_table + j, arg.prec ); printf("\n");
*/
            }
            compApp_set( arg.le_table + i, arg.up_table + j-1 );
        }
    }
/*     printf("\n");*/
    return NULL;
}

void _compApp_poly_parallel_taylor_inplace( compApp_poly_t p, const compApp_t c, const realRat_t radius, slong prec, slong nb_threads){
    
    taylor_shift_parallel_arg_t * args;
    pthread_t * threads;
    compApp_poly_t up;
    compApp_ptr up_table;
    slong size = p->length;
    
    compApp_poly_init2(up, size);
    compApp_poly_set_length(up, size);
    up_table = up->coeffs;
    
    args = (taylor_shift_parallel_arg_t *) malloc ( sizeof(taylor_shift_parallel_arg_t) * nb_threads );
    threads = (pthread_t *) malloc (sizeof(pthread_t) * nb_threads);
    
    /*reverse p*/
    _acb_poly_reverse(p->coeffs, p->coeffs, p->length, p->length);
    /*version with one thread*/
    if (nb_threads==1){
        args[0].size = size;
        args[0].prec = prec;
        args[0].le_table = p->coeffs;
        args[0].up_table = up_table;
        args[0].le_inf = 0;
        args[0].le_sup = size-1;
        args[0].up_inf = 0;
        args[0].up_sup = size-1;
        args[0].c = (compApp *) c;
/*         printf("create thread --- le_inf: %d, le_sup: %d, up_inf: %d, up_sup: %d\n", args[0].le_inf, args[0].le_sup, args[0].up_inf, args[0].up_sup);*/
        _compApp_poly_parallel_taylor_worker( &args[0] );
/*         pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_worker, &args[0]);*/
/*         pthread_join(threads[0], NULL);*/
    }
    /*blocking strategy*/
    else {
        ulong block_size = 64;//test
        ulong nb_gen     = (ulong) p->length/block_size;
/*         printf("nb_gen: %d\n", nb_gen);*/
        /*iterate generations*/
        
        ulong bound_gen = (p->length%block_size==0? nb_gen-1: nb_gen);
        for( ulong gen = 0; gen <= bound_gen; gen++){
            
            args[0].size     = size;
            args[0].prec     = prec;
            args[0].le_table = p->coeffs;
            args[0].up_table = up_table;
            args[0].le_inf   = 0;
            args[0].le_sup   = block_size-1;
            args[0].up_inf   = (gen    )*block_size;
            args[0].up_sup   = (gen + 1)*block_size - 1;
            args[0].c        = (compApp *) c;
            if (gen == nb_gen) //deal with irregular blocks
                args[0].up_sup   = size;
            
/*             printf("gen: %d, run: %d, thread %d: create thread --- le_inf: %d, le_sup: %d, up_inf: %d, up_sup: %d\n", gen, 0, 0, args[0].le_inf, args[0].le_sup, args[0].up_inf, args[0].up_sup);*/
            pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_worker, &args[0]);
/*             _compApp_poly_parallel_taylor_worker( &args[0] );*/
            
            for ( ulong run = 1; run <= gen; run++){
                int thread = run%nb_threads;
                int lastindex = (thread - 1)%nb_threads;
                lastindex = (lastindex>=0 ? lastindex : nb_threads+lastindex);
                
                if (thread == 0) {/*join nb_threads threads*/
/*                     printf("gen: %d, run: %d, thread %d: join %d threads\n", gen, run, thread, nb_threads);*/
                    for (int i = 0; i < nb_threads; i++)
                        pthread_join(threads[i], NULL);
                }
                
                args[thread].size     = size;
                args[thread].prec     = prec;
                args[thread].le_table = p->coeffs;
                args[thread].up_table = up_table;
                args[thread].le_inf   = args[lastindex].le_inf + block_size;
                args[thread].le_sup   = args[lastindex].le_sup + block_size;
                args[thread].up_inf   = args[lastindex].up_inf - block_size;
                args[thread].up_sup   = args[lastindex].up_sup - block_size;
                args[thread].c        = (compApp *) c;
                if (gen == nb_gen){ /*deal with irregular blocks*/
                    args[thread].up_sup   = size; /*the min in worker function deal with this*/
                    if (run==gen) {
                        args[thread].le_sup = size-1;
                    }
                }
/*                 printf("gen: %d, run: %d, thread %d: create thread --- le_inf: %d, le_sup: %d, up_inf: %d, up_sup: %d\n", gen, run, thread, args[thread].le_inf, args[thread].le_sup, args[thread].up_inf, args[thread].up_sup);*/
                pthread_create(&threads[thread], NULL, _compApp_poly_parallel_taylor_worker, &args[thread]);
/*                 _compApp_poly_parallel_taylor_worker( &args[thread] );*/
                
            }
            
            /*join (gen)%nb_threads*/
/*             printf("gen: %d, join %d threads\n", gen, gen%nb_threads +1);*/
            for (int i = 0; i <= gen%nb_threads; i++)
                pthread_join(threads[i], NULL);
            
            
        }
    }
    /*reverse p*/
    _acb_poly_reverse(p->coeffs, p->coeffs, p->length, p->length);
    
    compApp_poly_scale_realRat_in_place( p->coeffs, radius, p->length, prec );
    
    free(args);
    compApp_poly_clear(up);
}

void compApp_poly_parallel_taylor_inplace( compApp_poly_t p, const realRat_t creal, const realRat_t cimag, const realRat_t radius, slong prec, slong num_threads){
    
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, creal, prec);
    compApp_setimag_realRat(c, cimag, prec);
    
    _compApp_poly_parallel_taylor_inplace( p, c, radius, prec, num_threads);
    
    compApp_clear(c);
}

void compApp_poly_parallel_taylor( compApp_poly_t dest, const compApp_poly_t p, const realRat_t creal, const realRat_t cimag, const realRat_t radius, slong prec, slong num_threads){
    compApp_poly_set(dest, p);
    compApp_poly_parallel_taylor_inplace( dest, creal, cimag, radius, prec, num_threads);
}
/*
//     version with two threads
//     if (num_threads==2) {
//         slong block_size = (slong) size/num_threads;
// //         printf("block size: %d\n", block_size);
//         //first block
//         args[0].size = size;
//         args[0].prec = prec;
//         args[0].le_table = p->coeffs;
//         args[0].up_table = up_table;
//         args[0].le_inf = size - block_size;
//         args[0].le_sup = size-1;
//         args[0].up_inf = 0;
//         args[0].up_sup = block_size-1;
//         args[0].c = (compApp *) c;
//         _compApp_poly_parallel_taylor_worker( &args[0] );
//         
//     //     two others blocks
//         args[0].size = size;
//         args[0].prec = prec;
//         args[0].le_table = p->coeffs;
//         args[0].up_table = up_table;
//         args[0].le_inf = size - block_size;
//         args[0].le_sup = size-1;
//         args[0].up_inf = block_size;
//         args[0].up_sup = size-1;
//         args[0].c = (compApp *) c;
// //         _compApp_poly_parallel_taylor_worker( &args[0] );
//         pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_worker, &args[0]);
//     //     
//         args[1].size = size;
//         args[1].prec = prec;
//         args[1].le_table = p->coeffs;
//         args[1].up_table = up_table;
//         args[1].le_inf = 0;
//         args[1].le_sup = size-block_size-1;
//         args[1].up_inf = 0;
//         args[1].up_sup = size-block_size-1;
//         args[1].c = (compApp *) c;
// //         _compApp_poly_parallel_taylor_worker( &args[1] );
//         pthread_create(&threads[1], NULL, _compApp_poly_parallel_taylor_worker, &args[1]);
//         
//         for (int i = 0; i < num_threads; i++)
//             pthread_join(threads[i], NULL);
//     }
// 
// //     version with n=2^n' threads
//     if (num_threads>2){
//         slong block_size = (slong) size/num_threads;
// //         printf("block size: %d\n", block_size);
//         slong num_threads_iter, num_block_iter;
//         for (num_threads_iter =0; num_threads_iter< num_threads; num_threads_iter++){
//             
//             args[0].size     = size;
//             args[0].prec     = prec;
//             args[0].le_table = p->coeffs;
//             args[0].up_table = up_table;
//             args[0].le_inf   = size - block_size;
//             args[0].le_sup   = size-1;
//             args[0].up_inf   = (num_threads_iter    )*block_size;
//             args[0].up_sup   = (num_threads_iter + 1)*block_size - 1;
//             args[0].c        = (compApp *) c;
//             if (num_threads_iter == num_threads-1)
//                 args[0].up_sup   = size;
//             
//             pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_worker, &args[0]);
// //             _compApp_poly_parallel_taylor_worker( &args[0] );
//             
// //             printf("--------------------------------------------\n");
// //             printf("num_threads_iter: %d, num_block_iter: %d, le_inf: %d, le_sup: %d, up_inf: %d, up_sup: %d\n", 
// //                    num_threads_iter, 0, args[0].le_inf, args[0].le_sup, args[0].up_inf, args[0].up_sup);
//                 
//             for (num_block_iter =1; num_block_iter<= num_threads_iter; num_block_iter++){
//                 args[num_block_iter].size     = size;
//                 args[num_block_iter].prec     = prec;
//                 args[num_block_iter].le_table = p->coeffs;
//                 args[num_block_iter].up_table = up_table;
//                 args[num_block_iter].le_inf   = args[num_block_iter -1].le_inf - block_size;
//                 args[num_block_iter].le_sup   = args[num_block_iter -1].le_sup - block_size;
//                 args[num_block_iter].up_inf   = args[num_block_iter -1].up_inf - block_size;
//                 args[num_block_iter].up_sup   = args[num_block_iter -1].up_sup - block_size;
//                 args[num_block_iter].c        = (compApp *) c;
//                 //deal with irregular blocks
//                 if (num_threads_iter == num_threads-1){
//                     args[num_block_iter].up_sup   = size; //the min in worker function deal with this
//                     if (num_block_iter==num_threads_iter) {
//                         args[num_block_iter].le_inf = 0;
//                     }
//                 }
//                 pthread_create(&threads[num_block_iter], NULL, _compApp_poly_parallel_taylor_worker, &args[num_block_iter]);
// //                 _compApp_poly_parallel_taylor_worker( &args[num_block_iter] );
//                 
// //                 printf("num_threads_iter: %d, num_block_iter: %d, le_inf: %d, le_sup: %d, up_inf: %d, up_sup: %d\n", 
// //                    num_threads_iter, num_block_iter, args[num_block_iter].le_inf, args[num_block_iter].le_sup, 
// //                    args[num_block_iter].up_inf, args[num_block_iter].up_sup);
//             }
//             for (int i = 0; i <= num_threads_iter; i++)
//                 pthread_join(threads[i], NULL);
//         }
//     }
*/    