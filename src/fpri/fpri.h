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

#ifndef FPRI_H
#define FPRI_H

#ifdef FPRI_INLINE_C
#define FPRI_INLINE
#else
#define FPRI_INLINE static __inline__
#endif

#include <stdlib.h>
#include <math.h>
#include <fenv.h>       /* fesetround, FE_* */
#include "arb.h"
#include "flint/fmpq.h"

#define FPRI_MIN(A,B) (A<=B? A : B)
#define FPRI_MAX(A,B) (A>=B? A : B)
#define FPRI_ABS(A)   (A<=0? -A : A)
#define FPRI_PREC 53

#define FPRI_DEFAULT_MALLOC  malloc
#define FPRI_DEFAULT_REALLOC realloc
#define FPRI_DEFAULT_CALLOC  calloc
#define FPRI_DEFAULT_FREE    free

#ifdef __cplusplus
extern "C" {
#endif
    
FPRI_INLINE void fpri_lib_init(){
    fesetround(FE_UPWARD);
}

/* memory functions */
           void * fpri_malloc (size_t size);             
           void * fpri_realloc(void * ptr, size_t size);
           void * fpri_calloc (size_t num, size_t size); 
           void   fpri_free   (void * ptr);                  

/* default memory functions called in functions above */
FPRI_INLINE void * _fpri_malloc(size_t size)              { return FPRI_DEFAULT_MALLOC(size); }
FPRI_INLINE void * _fpri_realloc(void * ptr, size_t size) { return FPRI_DEFAULT_REALLOC(ptr, size); }
FPRI_INLINE void * _fpri_calloc(size_t num, size_t size)  { return FPRI_DEFAULT_CALLOC(num, size); }
FPRI_INLINE void   _fpri_free(void * ptr)                 {        FPRI_DEFAULT_FREE(ptr); }

            void __fpri_set_memory_functions( void *(*alloc_func)   (size_t),
                                              void *(*calloc_func)  (size_t, size_t), 
                                              void *(*realloc_func) (void *, size_t),
                                              void  (*free_func)    (void *) );

typedef double number;
typedef struct {
    /* Assuming rounding upward; 
     *low stores MINUS the lower 
     * bound of the interval;*/
    number low;
    number upp;
} fpri_struct;

typedef fpri_struct fpri_t[1];
typedef fpri_struct * fpri_ptr;
typedef const fpri_struct * fpri_srcptr;

/* memory managment */
// FPRI_INLINE void fpri_init (fpri_t x) { x->low=0.; x->upp=0.; }
FPRI_INLINE void fpri_init (fpri_t x) { }
FPRI_INLINE void fpri_clear(fpri_t x) {  }
            void fpri_swap (fpri_t x, fpri_t y);
            
FPRI_INLINE fpri_ptr _fpri_vec_init( slong n ) {
    return (fpri_ptr) fpri_malloc ( n*sizeof(fpri_struct) );
}
FPRI_INLINE void     _fpri_vec_clear( fpri_ptr v, slong n ) {
    fpri_free(v);
}
/* setting */
FPRI_INLINE void fpri_zero      (fpri_t x                         ) { x->low=0.; x->upp=0.;  }
FPRI_INLINE void fpri_vec_zero ( fpri_ptr v, slong n ) {
    slong i;
    for (i=0; i<n; i++)
        fpri_zero(v + i);
}

FPRI_INLINE void fpri_one       (fpri_t x                         ) { x->low=-1.; x->upp=1.; }
FPRI_INLINE void fpri_set       (fpri_t y, const fpri_t x         ) { y->low = x->low; y->upp=x->upp; }
FPRI_INLINE void fpri_set_d     (fpri_t y, const double x         ) { y->low = -x; y->upp = x; }
FPRI_INLINE void fpri_set_d_d   (fpri_t y, const double l, const double u ) { y->low = -l; y->upp = u; }
            void fpri_set_arb   (fpri_t y, const arb_t  x         );
FPRI_INLINE void fpri_set_si    (fpri_t y, const slong x          ) { y->low = -x; y->upp = x; }
/* special values */
FPRI_INLINE void fpri_set_inf    (fpri_t y) { y->low = +INFINITY; y->upp = +INFINITY; }
FPRI_INLINE void fpri_set_nan    (fpri_t y) { y->low = NAN; y->upp = NAN; }
/* getting */
            int  fpri_get_arb   (arb_t  y, const fpri_t x         );
            
/* comparisons */
FPRI_INLINE int  fpri_is_zero (const fpri_t x) { return (x->low==0.) && (x->upp==0.); }
FPRI_INLINE int  fpri_is_one  (const fpri_t x) { return (x->low==-1) && (x->upp==1); }
FPRI_INLINE int  fpri_is_finite  (const fpri_t x) { return (isfinite(x->low) && isfinite(x->upp)); }
FPRI_INLINE int  fpri_is_inf     (const fpri_t x) { return (isinf   (x->low) || isinf(x->upp)); }
FPRI_INLINE int  fpri_is_nan     (const fpri_t x) { return (isnan   (x->low) || isnan(x->upp)); }
FPRI_INLINE int  fpri_is_exact(const fpri_t x) { return (fpri_is_finite(x) && (-x->low==x->upp)); }
FPRI_INLINE int  fpri_contains_zero(const fpri_t x) { return (x->low>=0.) && (x->upp>=0.); }
FPRI_INLINE int  fpri_lt (const fpri_t x, const fpri_t y) { return (x->upp < (-y->low)); }
FPRI_INLINE int  fpri_ge (const fpri_t x, const fpri_t y) { return ((-x->low) >= y->upp); }
FPRI_INLINE int  fpri_gt (const fpri_t x, const fpri_t y) { return ((-x->low) >  y->upp); }

/* arithmetic operations */
/* support aliazing */
            void fpri_neg     (fpri_t z, const fpri_t x); 
            void fpri_abs     (fpri_t z, const fpri_t x); 
            void fpri_inv     (fpri_t z, const fpri_t x);
            void fpri_sqr     (fpri_t z, const fpri_t x);
            /* sqrt is not implemented because it would require to change the rounding mode */
            /* in order to get a lower bound on sqrt(l) */
//             void fpri_sqrt     (fpri_t z, const fpri_t x);
            void fpri_add     (fpri_t res, const fpri_t x, const fpri_t y);
            void fpri_sub     (fpri_t res, const fpri_t x, const fpri_t y);
            void _fpri_mul    (fpri_t res, const fpri_t x, const fpri_t y);
FPRI_INLINE void fpri_mul     (fpri_t res, const fpri_t x, const fpri_t y){
    if (x==y)
        fpri_sqr(res, x);
    else
        _fpri_mul(res, x, y);
}
FPRI_INLINE void fpri_mul_si  (fpri_t res, const fpri_t x, slong y){ res->low = y*(x->low); res->upp=y*(x->upp); }
/* does not support aliazing */
            void _fpri_div     (fpri_t res, const fpri_t x, const fpri_t y);
/* support aliazing */
            void fpri_div     (fpri_t res, const fpri_t x, const fpri_t y);

            void fpri_pow_ui           (fpri_t res, const fpri_t x, ulong p);
FPRI_INLINE void fpri_div_si  (fpri_t res, const fpri_t x, slong y){ 
    if (y==0) {
        res->upp = +INFINITY;
        res->low = +INFINITY;
    } else {
        res->low = (x->low)/y; 
        res->upp = (x->upp)/y;
    }
}       

/* sets z to z + xy */
void fpri_addmul( fpri_t z, const fpri_t x, const fpri_t y );

/* printing */            
            void fpri_fprint    (FILE * file, const fpri_t x);
FPRI_INLINE void fpri_print     (const fpri_t x)                    { fpri_fprint(stdout, x); }

            void fpri_fprint_s    (FILE * file, const fpri_t x);
FPRI_INLINE void fpri_print_s     (const fpri_t x)                    { fpri_fprint_s(stdout, x); }

#ifdef __cplusplus
}
#endif

#endif
