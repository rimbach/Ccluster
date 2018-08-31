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
    
    doubCompApp_poly_t p, q, r;
    doubCompApp_poly_init(p);
    doubCompApp_poly_init(q);
    doubCompApp_poly_init(r);
    
    compApp_poly_t pp, qq, rr;
    compApp_poly_init(pp);
    compApp_poly_init(qq);
    compApp_poly_init(rr);
    
    realRat_poly_t pbern, pbern10;
    realRat_poly_init(pbern);
    realRat_poly_init(pbern10);
    bernoulli_polynomial(pbern, 32);
    compApp_poly_set2_fmpq_poly( pp, pbern, pbern, prec);
    doubCompApp_poly_set_compApp_poly(p,pp);
    bernoulli_polynomial(pbern10, 32);
    compApp_poly_set2_fmpq_poly( qq, pbern10, pbern10, prec);
    doubCompApp_poly_set_compApp_poly(q,qq);
    
    acb_poly_mullow_classical( rr, pp, qq, pp->length + qq->length, prec);
    doubCompApp_poly_mul_classical( r, p, q);
    
//     printf("bern * bern10 acb_poly: \n"); compApp_poly_printd(rr, 10); printf("\n\n");
    
//     printf("bern * bern10 doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    clock_t start;
    double time_in_arb, time_in_arb_q, time_in_doub, time_in_doub_q;
    
    int degree;
    int nbops = 1000;
    for (int d = 6; d< 10; d++){
        degree = pow(2,d);
        bernoulli_polynomial(pbern, degree);
        bernoulli_polynomial(pbern10, degree);
        compApp_poly_set2_fmpq_poly( pp, pbern, pbern, prec);
//         compApp_poly_set_fmpq_poly( pp, pbern, prec);
        doubCompApp_poly_set_compApp_poly(p,pp);
        compApp_poly_set2_fmpq_poly( qq, pbern10, pbern10, prec);
//         compApp_poly_set_fmpq_poly( qq, pbern10, prec);
        doubCompApp_poly_set_compApp_poly(q,qq);
        
        start = clock();
        for (int i=0; i<nbops; i++)
            acb_poly_mullow_classical( rr, pp, qq, pp->length + qq->length, prec);
        time_in_arb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("degree: %d, time for arb classical: %lf \n", degree, time_in_arb);
        
        start = clock();
        for (int i=0; i<nbops; i++)
            acb_poly_mul( rr, pp, qq, prec);
        time_in_arb_q = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("degree: %d, time for arb quick: %lf \n", degree, time_in_arb_q);
        
        start = clock();
        for (int i=0; i<nbops; i++)
            doubCompApp_poly_mul_classical( r, p, q);
        time_in_doub = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("degree: %d, time for doub classical: %lf, ratio: %lf, ratio: %lf \n", degree, time_in_doub, 
               time_in_arb/time_in_doub, time_in_arb_q/time_in_doub);
        
        start = clock();
        for (int i=0; i<nbops; i++)
            doubCompApp_poly_mul_karatsuba( r, p, q);
        time_in_doub_q = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("degree: %d, time for doub karatsuba: %lf, ratio: %lf \n", degree, time_in_doub_q, 
               time_in_arb_q/time_in_doub_q);
    }
    
    realRat_poly_clear(pbern);
    realRat_poly_clear(pbern10);
    compApp_poly_clear(pp);
    compApp_poly_clear(qq);
    compApp_poly_clear(rr);
    doubCompApp_poly_clear(p);
    doubCompApp_poly_clear(q);
    doubCompApp_poly_clear(r);
    
    return 0;
    
}
