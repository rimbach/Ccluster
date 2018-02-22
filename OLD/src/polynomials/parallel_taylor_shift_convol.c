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

typedef struct {
    compApp_ptr coeffs; /* coefficients */
    slong       length; /* length */
    compApp *   c;      /* center */
    slong       prec;
} taylor_parallel_convol_arg_t;

void * _compApp_poly_parallel_taylor_convol_worker( void * arg_ptr ){
    
    taylor_parallel_convol_arg_t arg = *((taylor_parallel_convol_arg_t *) arg_ptr);
    _acb_poly_taylor_shift_convolution(arg.coeffs, arg.c, arg.length, arg.prec);
    return NULL;
    
}

typedef struct {
    compApp_ptr p1_coeffs; /* coefficients */
    slong       p1_length; /* length */
    compApp_ptr p2_coeffs; /* coefficients */
    slong       p2_length; /* length */
    slong       prec;
} taylor_parallel_add_arg_t;

void * _compApp_poly_parallel_taylor_add_worker( void * arg_ptr ){
    
    taylor_parallel_add_arg_t arg = *((taylor_parallel_add_arg_t *) arg_ptr);
    for (int i=0; i<arg.p2_length;i++) {
        compApp_add((arg.p1_coeffs)+i, (arg.p1_coeffs)+i, (arg.p2_coeffs)+i, arg.prec);
    }
    return NULL;
    
}

void compApp_poly_parallel_taylor_convol(compApp_poly_t dest, const compApp_poly_t p, const compApp_t c, const realRat_t radius, slong prec, slong num_threads){
    
    taylor_parallel_convol_arg_t * args;
    args = (taylor_parallel_convol_arg_t *) malloc ( sizeof(taylor_parallel_convol_arg_t) * num_threads );
    
    taylor_parallel_add_arg_t * args_add;
    args_add = (taylor_parallel_add_arg_t *) malloc ( sizeof(taylor_parallel_add_arg_t) * (num_threads-1) );
    
    pthread_t * threads;
    threads = (pthread_t *) malloc (sizeof(pthread_t) * num_threads);
    
    if (num_threads==1){
        compApp_poly_set(dest, p);
        args[0].coeffs = dest->coeffs;
        args[0].length = dest->length;
        args[0].c = (compApp *) c;
        args[0].prec = prec;
        _compApp_poly_parallel_taylor_convol_worker( &args[0] );
    }
    
    if (num_threads==2){
        slong poly_size = (slong) (p->length)/num_threads;
//         printf("poly_size: %d\n", poly_size);
        compApp_poly_t temp;
        compApp_poly_init2(temp, poly_size);
        
        compApp_poly_set(dest, p);        
        acb_poly_set_trunc(temp, p, poly_size);
        for (int i=0;i<poly_size;i++) 
            compApp_zero((dest->coeffs)+i);
        
//         printf("dest--: \n"); compApp_poly_printd(dest, prec); printf("\n\n");
//         printf("temp--: \n"); compApp_poly_printd(temp, prec); printf("\n\n");
        
        args[0].coeffs = dest->coeffs;
        args[0].length = dest->length;
        args[0].c = (compApp *) c;
        args[0].prec = prec;
//         _compApp_poly_parallel_taylor_convol_worker( &args[0] );
        pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_convol_worker, &args[0]);
        
        args[1].coeffs = temp->coeffs;
        args[1].length = temp->length;
        args[1].c = (compApp *) c;
        args[1].prec = prec;
//         _compApp_poly_parallel_taylor_convol_worker( &args[1] );
        pthread_create(&threads[1], NULL, _compApp_poly_parallel_taylor_convol_worker, &args[1]);
        
//         pthread_join(threads[0], NULL);
        for (int i = 0; i < num_threads; i++)
                pthread_join(threads[i], NULL);
        
        compApp_poly_add(dest, dest, temp, prec);
        compApp_poly_clear(temp);
    }
//     printf("ici: \n");

    if (num_threads==4){
        slong poly_size = (slong) (p->length)/num_threads;
        
//         compApp_poly_set(dest, p);
        compApp_poly_init2(dest, p->length);
        compApp_poly_set_length(dest, p->length);
        for (int j=0;j<(num_threads-1)*poly_size; j++)
            compApp_zero((dest->coeffs)+j);
        for (int j=(num_threads-1)*poly_size;j<p->length; j++)
            compApp_set((dest->coeffs)+j, (p->coeffs)+j);
        
        args[0].coeffs = dest->coeffs;
        args[0].length = dest->length;
        args[0].c = (compApp *) c;
        args[0].prec = prec;
//         _compApp_poly_parallel_taylor_convol_worker( &args[0] );
        pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_convol_worker, &args[0]);
        
        compApp_poly_t * temp;
        temp = (compApp_poly_t *) malloc ((num_threads-1)*sizeof(compApp_poly_t));
        
        for(int i=0; i<num_threads-1; i++) {
            
//             compApp_poly_init(temp[i]);
//             acb_poly_set_trunc(temp[i], p, (i+1)*poly_size);
//             for (int j=0;j<i*poly_size; j++)
//                 compApp_zero((temp[i]->coeffs)+j);
            
            compApp_poly_init2(temp[i], (i+1)*poly_size);
            compApp_poly_set_length(temp[i], (i+1)*poly_size);
            for (int j=0;j<i*poly_size; j++)
                compApp_zero((temp[i]->coeffs)+j);
            for (int j=i*poly_size;j<(i+1)*poly_size; j++)
                compApp_set((temp[i]->coeffs)+j, (p->coeffs)+j);
            
            args[i+1].coeffs = temp[i]->coeffs;
            args[i+1].length = temp[i]->length;
            args[i+1].c = (compApp *) c;
            args[i+1].prec = prec;
//         _compApp_poly_parallel_taylor_convol_worker( &args[0] );
            pthread_create(&threads[i+1], NULL, _compApp_poly_parallel_taylor_convol_worker, &args[i+1]);
        }
        
//         pthread_join(threads[0], NULL);
        
//         for (int i = 0; i < num_threads; i++)
//                 pthread_join(threads[i], NULL);
        
//         for(int i=0; i<num_threads-2; i++) {
//             compApp_poly_add(temp[i+1], temp[i+1], temp[i], prec);
//             compApp_poly_clear(temp[i]);
//         }
//         compApp_poly_add(dest, dest, temp[num_threads-2], prec);
//         compApp_poly_clear(temp[num_threads-2]);
           
           pthread_join(threads[1], NULL);
           pthread_join(threads[2], NULL);
           args_add[0].p1_coeffs = temp[1]->coeffs;
           args_add[0].p1_length = temp[1]->length;
           args_add[0].p2_coeffs = temp[0]->coeffs;
           args_add[0].p2_length = temp[0]->length;
           args_add[0].prec = prec;
//            printf("la: temp[0]->length: %d, temp[1]->length: %d\n", temp[0]->length, temp[1]->length);
//            printf("temp[0]--: \n"); compApp_poly_printd(temp[0], prec); printf("\n\n");
//            printf("temp[1]--: \n"); compApp_poly_printd(temp[1], prec); printf("\n\n");
//            _compApp_poly_parallel_taylor_add_worker(&args_add[0]);
           pthread_create(&threads[2], NULL, _compApp_poly_parallel_taylor_add_worker, &args_add[0]);
//            printf("la\n");
           pthread_join(threads[2], NULL);
           pthread_join(threads[3], NULL);
           
           args_add[1].p1_coeffs = temp[2]->coeffs;
           args_add[1].p1_length = temp[2]->length;
           args_add[1].p2_coeffs = temp[1]->coeffs;
           args_add[1].p2_length = temp[1]->length;
           args_add[1].prec = prec;
//            _compApp_poly_parallel_taylor_add_worker(&args_add[1]);
           pthread_create(&threads[3], NULL, _compApp_poly_parallel_taylor_add_worker, &args_add[1]);
           
           pthread_join(threads[3], NULL);
           
           pthread_join(threads[0], NULL);
//            printf("la\n");
           args_add[2].p1_coeffs = dest->coeffs;
           args_add[2].p1_length = dest->length;
           args_add[2].p2_coeffs = temp[2]->coeffs;
           args_add[2].p2_length = temp[2]->length;
           args_add[2].prec = prec;
           
//            _compApp_poly_parallel_taylor_add_worker(&args_add[2]);
           pthread_create(&threads[0], NULL, _compApp_poly_parallel_taylor_add_worker, &args_add[2]);
           
           pthread_join(threads[0], NULL);
           
           for (int i = 0; i < num_threads-1; i++)
                compApp_poly_clear(temp[i]);
           free(temp);
        
    }
    
    compApp_poly_scale_realRat_in_place( dest->coeffs, radius, dest->length, prec );
    
    free(threads);
    free(args);
    free(args_add);
}

void compApp_poly_parallel_taylor_convol_forJulia(compApp_poly_t dest, const compApp_poly_t p, const realRat_t creal, const realRat_t cimag, const realRat_t radius, slong prec, slong num_threads){
    compApp_t c;
    compApp_init(c);
    compApp_setreal_realRat(c, creal, prec);
    compApp_setimag_realRat(c, cimag, prec);
    
//     compApp_poly_init2(dest, p->length);
    
    compApp_poly_parallel_taylor_convol(dest, p, c,radius, prec, num_threads);
    compApp_clear(c);
}