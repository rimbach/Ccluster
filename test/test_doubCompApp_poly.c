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
    
    doubCompApp_poly_t p,q, r;
    doubCompApp_poly_init(p);
    doubCompApp_poly_init(q);
    doubCompApp_poly_init(r);
    
    printf("zero: \n");
    doubCompApp_poly_zero(p);
    doubCompApp_poly_print(p);
    
    printf("one: \n");
    doubCompApp_poly_one(p);
    doubCompApp_poly_print(p);
    
    compApp_poly_t pp, qq, rr;
    compApp_poly_init(pp);
    compApp_poly_init(qq);
    compApp_poly_init(rr);
    
    realRat_poly_t pbern, pbern10;
    realRat_poly_init(pbern);
    realRat_poly_init(pbern10);
    bernoulli_polynomial(pbern, 5);
//     mignotte_polynomial(pmign, 5, (int)64/2-1);
    
    compApp_poly_set_fmpq_poly( pp, pbern, prec);
    doubCompApp_poly_set_compApp_poly(p,pp);
    printf("bern acb_poly: \n"); compApp_poly_printd(pp, 10); printf("\n\n");
    
    printf("bern doub_poly: \n"); doubCompApp_poly_print(p); printf("\n\n");
    
    doubCompApp_poly_neg(r,p);
    printf("-bern doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    bernoulli_polynomial(pbern10, 10);
    compApp_poly_set_fmpq_poly( qq, pbern10, prec);
    doubCompApp_poly_set_compApp_poly(q,qq);
    printf("bern10 acb_poly: \n"); compApp_poly_printd(qq, 10); printf("\n\n");
    
    printf("bern10 doub_poly: \n"); doubCompApp_poly_print(q); printf("\n\n");
    
    compApp_poly_add( rr, pp, qq, prec);
    doubCompApp_poly_add(r,p,q);
    
    printf("bern + bern10 acb_poly: \n"); compApp_poly_printd(rr, 10); printf("\n\n");
    
    printf("bern + bern10 doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    doubCompApp_poly_add(r,q,p);
    
    printf("bern10 + bern doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    compApp_poly_sub( rr, pp, qq, prec);
    doubCompApp_poly_sub(r,p,q);
    
    printf("bern - bern10 acb_poly: \n"); compApp_poly_printd(rr, 10); printf("\n\n");
    
    printf("bern - bern10 doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    compApp_poly_sub( rr, qq, pp, prec);
    doubCompApp_poly_sub(r,q,p);
    
    printf("bern10 - bern acb_poly: \n"); compApp_poly_printd(rr, 10); printf("\n\n");
    
    printf("bern10 - bern doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
    compApp_poly_shift_left( rr, qq, 5);
    doubCompApp_poly_shift_left(r,q,5);
    
    printf("bern10 shifted acb_poly: \n"); compApp_poly_printd(rr, 10); printf("\n\n");
    
    printf("bern10 shifted doub_poly: \n"); doubCompApp_poly_print(r); printf("\n\n");
    
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
