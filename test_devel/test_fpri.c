#include <stdio.h>
#include "fpri/fpri.h"

int main( ){
    fpri_lib_init();
    
    fpri_t a;
    fpri_init(a);
    printf("un-initialized fpri: "); fpri_print(a); printf("\n");
    fpri_clear(a);
    
    fpri_t b;
    fpri_init(b);
    fpri_zero(b);
    printf("zero: "); fpri_print(b); printf("\n");
    fpri_one(b);
    printf("one: "); fpri_print(b); printf("\n");
    fpri_set_d(b, 0.1);
    printf("0.1 converted in double: "); fpri_print(b); printf("\n");
    
    slong prec = 53;
    arb_t bb;
    arb_init(bb);
    arb_set_si(bb, 10);
    arb_inv(bb, bb, prec);
    printf("arb 1/10: "); arb_printd(bb, 16); printf("\n");
    fpri_set_arb(b, bb);
    printf("arb 1/10 converted in fpri: "); fpri_print(b); printf("\n");
    arb_zero(bb);
    fpri_get_arb(bb, b);
    printf("arb 1/10 converted in fpri, reconverted in arb: "); arb_printd(bb, 16); printf("\n");
    arb_clear(bb);
    
    fpri_zero(b);
    printf(" zero contains 0: %d\n", fpri_contains_zero(b));
    fpri_one(b);
    fpri_inv(b,b);
    printf("inv(one): "); fpri_print(b); printf("\n");
    fpri_zero(b);
    fpri_inv(b,b);
    printf("inv(zero): "); fpri_print(b); printf("\n");
    fpri_set_d(b, 2);
    fpri_inv(b,b);
    printf("inv(2): "); fpri_print(b); printf("\n");
    fpri_set_d(b, 10);
    fpri_inv(b,b);
    printf("inv(10): "); fpri_print(b); printf("\n");
    fpri_neg(b,b);
    printf("neg(inv(10)): "); fpri_print(b); printf("\n");
    fpri_neg(b,b);
    printf("neg(neg(inv(10))): "); fpri_print(b); printf("\n");
    fpri_sqr(b,b);
    printf("sqr(inv(10)): "); fpri_print(b); printf("\n");
    
    fpri_set_d_d(b, 0.4, 0.5);
    fpri_abs(b, b);
    printf("abs[0.4, 0.5]: "); fpri_print(b); printf("\n");
    fpri_set_d_d(b, -0.5, -0.4);
    fpri_abs(b, b);
    printf("abs[-0.5, -0.4]: "); fpri_print(b); printf("\n");
    fpri_set_d_d(b, -0.5, 0.4);
//     printf("[-0.5, 0.4]: "); fpri_print(b); printf("\n");
    fpri_abs(b, b);
    printf("abs[-0.5, 0.4]: "); fpri_print(b); printf("\n");
    fpri_clear(b);
    
    fpri_t c, d, e;
    fpri_init(c);
    fpri_init(d);
    fpri_init(e);
    fpri_set_d(c, 10);
    fpri_inv(c,c);
    fpri_neg(d,c);
    fpri_add(d, c, d);
    printf("inv(10) + neg(inv(10)): "); fpri_print(d); printf("\n");
    fpri_sub(c, c, c);
    printf("inv(10) - inv(10): "); fpri_print(d); printf("\n");
    
    fpri_set_d(e, 10);
    fpri_set_d(c, 3);
    fpri_set_d(d, 7);
    fpri_inv(e, e);
    fpri_mul(c, c, e);
    printf("3/10: "); fpri_print(c); printf("\n");
    fpri_mul(d, d, e);
    printf("7/10: "); fpri_print(d); printf("\n");
    fpri_mul(e, c, d);
    printf("3/10*7/10: "); fpri_print(e); printf("\n");
    
    fpri_set_d(e, 10);
    fpri_set_d(c, 3);
    fpri_set_d(d, 7);
    fpri_div(c, c, e);
    printf("3/10: "); fpri_print(c); printf("\n");
    fpri_div(d, d, e);
    printf("7/10: "); fpri_print(d); printf("\n");
    
    fpri_pow_ui(e, c, 0);
    printf("(3/10)^0: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 1);
    printf("(3/10)^1: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 2);
    printf("(3/10)^2: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 3);
    printf("(3/10)^3: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 4);
    printf("(3/10)^4: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 10);
    printf("(3/10)^10: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 20);
    printf("(3/10)^20: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 100);
    printf("(3/10)^100: "); fpri_print(e); printf("\n");
    fpri_pow_ui(e, c, 1000);
    printf("(3/10)^1000: "); fpri_print(e); printf("\n");
    fpri_clear(c);
    fpri_clear(d);
    fpri_clear(e);
    
//     __fpri_set_memory_functions( flint_malloc, flint_calloc, flint_realloc, flint_free );
    fpri_ptr v = _fpri_vec_init( 10 );
    fpri_vec_zero(v, 10);
//     __fpri_set_memory_functions( malloc, calloc, realloc, free );
    slong i;
    printf("[ ");
    for (i=0;i<10; i++) {
        fpri_print(v+i); printf(", ");
    }
    printf(" ]\n");
        
    _fpri_vec_clear(v, 10);
    
    flint_cleanup();
    
    return 0;
}
