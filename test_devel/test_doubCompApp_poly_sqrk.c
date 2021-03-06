#include <stdio.h>

#include <time.h>

#include "numbers/compApp.h"
#include "doubApp/doubCompApp.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/realRat_poly.h"
#include "doubApp/doubCompApp_poly.h"

#include <fenv.h> 

int main() {
    fesetround(FE_UPWARD);
    
    slong prec = 53;
    
    /* degree 0 */
    doubCompApp_poly_t p, r, q;
    doubCompApp_poly_init(p);
    doubCompApp_poly_init(r);
    doubCompApp_poly_init(q);
    
    doubCompApp_poly_one(p);
    
    doubCompApp_poly_sqr_karatsuba(r, p);
    printf("1 * 1 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    /*degree 1*/
    doubCompApp_poly_shift_left(p, p, 1);
    doubCompApp_mul_si(p->coeffs+1, p->coeffs+1, 2);
    doubCompApp_one(p->coeffs);
    doubCompApp_poly_sqr_karatsuba(r, p);
    printf("2x+1 : \n"); doubCompApp_poly_print(p); printf("\n\n");
    printf("2x+1 * 2x+1 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
    doubCompApp_poly_mul_classical(r, p, p);
    printf("2x+1 * 2x+1 class: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    /* degree 2 */
    doubCompApp_poly_shift_left(p, p, 1);
    doubCompApp_mul_si(p->coeffs+2, p->coeffs+2, 2);
    doubCompApp_mul_si(p->coeffs+1, p->coeffs+1, 3);
    doubCompApp_one(p->coeffs);
    doubCompApp_poly_sqr_karatsuba(r, p);
    printf("4x2+3x+1 : \n"); doubCompApp_poly_print(p); printf("\n\n");
    printf("4x2+3x+1 * 4x2+3x+1 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
 doubCompApp_poly_mul_classical(r, p, p);
    printf("4x2+3x+1 * 4x2+3x+1 class: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    int degree = 40;
    realRat_poly_t pbern, pbern2;
    realRat_poly_init(pbern);
    realRat_poly_init(pbern2);
    bernoulli_polynomial(pbern, degree);
    bernoulli_polynomial(pbern2, degree);
    realRat_poly_add(pbern2, pbern2, pbern);
    compApp_poly_t pp, rr;
    compApp_poly_init(pp);
    compApp_poly_init(rr);
    compApp_poly_set2_fmpq_poly( pp, pbern, pbern2, prec);
//     compApp_poly_set_fmpq_poly( pp, pbern, prec);
    doubCompApp_poly_set_compApp_poly(p,pp);
    doubCompApp_poly_sqr_karatsuba(r, p);
    printf("degree %d, kara: \n", degree); doubCompApp_poly_print(r); printf("\n\n");
    doubCompApp_poly_mul_classical(q, p, p);
    printf("degree %d, class: \n", degree); doubCompApp_poly_print(q); printf("\n\n");
    
    doubCompApp_poly_sub(r, q, r);
    printf("degree %d, diff: \n", degree); doubCompApp_poly_print(r); printf("\n\n");
    
//     acb_poly_mullow_classical( rr, pp, qq, pp->length + qq->length, prec);
//     acb_poly_mul( rr, pp, qq, prec);
//     printf("degree %d, acb_: \n", degree); compApp_poly_printd(rr, 10); printf("\n\n");
    
    doubCompApp_poly_clear(p);
//     doubCompApp_poly_clear(q);
    doubCompApp_poly_clear(r);
    realRat_poly_clear(pbern);
    realRat_poly_clear(pbern2);
    compApp_poly_clear(pp);
//     compApp_poly_clear(qq);
    compApp_poly_clear(rr);
    return 0;
    
}
