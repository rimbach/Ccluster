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
    doubCompApp_poly_t p, q, r;
    doubCompApp_poly_init(p);
    doubCompApp_poly_init(q);
    doubCompApp_poly_init(r);
    
    doubCompApp_poly_one(p);
//     doubCompApp_poly_onei(q);
//     doubCompApp_poly_add(p, p, q);
    doubCompApp_poly_set(q, p);
    
    doubCompApp_mul_si(q->coeffs, q->coeffs, 2);
    
    doubCompApp_poly_mul_karatsuba(r, p, q);
    printf("1 * 2 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    /* degree 1 */
    doubCompApp_poly_shift_left(q, q, 1);
    doubCompApp_poly_add(p, p, q);
    
    doubCompApp_poly_fit_length(q,2);
    doubCompApp_one(q->coeffs);
    doubCompApp_one(q->coeffs+1);
    doubCompApp_mul_si(q->coeffs+1, q->coeffs+1, 3);
    doubCompApp_mul_si(q->coeffs, q->coeffs, 2);
    doubCompApp_poly_mul_karatsuba(r, p, q);
    printf("2x+1 : \n"); doubCompApp_poly_print(p); printf("\n\n");
    printf("3x+2 : \n"); doubCompApp_poly_print(q); printf("\n\n");
    printf("2x+1 * 3x+2 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    doubCompApp_poly_mul_classical(r, p, q);
    printf("2x+1 * 3x+2 class: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    /* degree 2 */
    doubCompApp_poly_shift_left(p, p, 1);
    doubCompApp_mul_si(p->coeffs+2, p->coeffs+2, 2);
    doubCompApp_mul_si(p->coeffs+1, p->coeffs+1, 3);
    doubCompApp_one(p->coeffs);
    doubCompApp_poly_shift_left(q, q, 1);
    doubCompApp_mul_si(q->coeffs+2, q->coeffs+2, 2);
    doubCompApp_one(q->coeffs);
    doubCompApp_mul_si(q->coeffs, q->coeffs, 5);
    doubCompApp_poly_mul_karatsuba(r, p, q);
    printf("4x2+3x+1 : \n"); doubCompApp_poly_print(p); printf("\n\n");
    printf("6x2+2x+5 : \n"); doubCompApp_poly_print(q); printf("\n\n");
    printf("4x2+3x+1 * 6x2+2x+5 kara: \n"); doubCompApp_poly_print(r); printf("\n\n");
    doubCompApp_poly_mul_classical(r, p, q);
    printf("4x2+3x+1 * 6x2+2x+5 class: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    int degree1 = 7;
    int degree2 = 7;
    realRat_poly_t pbern, pbern2;
    realRat_poly_init(pbern);
    realRat_poly_init(pbern2);
    bernoulli_polynomial(pbern, degree1);
    bernoulli_polynomial(pbern2, degree2);
//     realRat_poly_add(pbern2, pbern2, pbern);
    compApp_poly_t pp, qq, rr;
    compApp_poly_init(pp);
    compApp_poly_init(qq);
    compApp_poly_init(rr);
    compApp_poly_set_fmpq_poly( pp, pbern, prec);
    compApp_poly_set_fmpq_poly( qq, pbern2, prec);
    doubCompApp_poly_set_compApp_poly(p,pp);
    doubCompApp_poly_set_compApp_poly(q,qq);
    
    doubCompApp_poly_mul_karatsuba(r, p, q);
    printf("degrees %d,%d, kara: \n", degree1, degree2); doubCompApp_poly_print(r); printf("\n\n");
    doubCompApp_poly_mul_classical(r, p, q);
    printf("degrees %d,%d, class: \n", degree1, degree2); doubCompApp_poly_print(r); printf("\n\n");
    
    
    doubCompApp_poly_clear(p);
    doubCompApp_poly_clear(q);
    doubCompApp_poly_clear(r);
    realRat_poly_clear(pbern);
    realRat_poly_clear(pbern2);
    compApp_poly_clear(pp);
    compApp_poly_clear(qq);
    compApp_poly_clear(rr);
    return 0;
    
}
