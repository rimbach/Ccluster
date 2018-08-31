#include <stdio.h>

#include <time.h>

#include "numbers/compApp.h"
#include "doubApp/doubCompApp.h"

#include <fenv.h> 

int main() {
    fesetround(FE_UPWARD);
    
    slong prec = 53;
    
    doubCompApp_t x, y, z;
//     doubCompApp_init(x);
//     doubCompApp_init(y);
    
    printf("--- set: ---\n");
    doubCompApp_zero(x);
    printf("zero: \n"); doubCompApp_print(x); printf("\n");
    doubCompApp_one(x);
    printf("one:  \n"); doubCompApp_print(x); printf("\n");
    doubCompApp_onei(x);
    printf("onei: \n"); doubCompApp_print(x); printf("\n");
    doubCompApp_set(y, x);
    printf("onei: \n"); doubCompApp_print(y); printf("\n");
    
    realApp_t rr;
    realApp_init(rr);
    fmpq_t r;
    fmpq_init(r);
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(rr, r, prec);
    
    compApp_t xx, yy, zz;
    compApp_init(xx);
    compApp_init(yy);
    compApp_init(zz);
    compApp_set_real_realApp(xx, rr);
    compApp_set_imag_realApp(xx, rr);
    doubCompApp_set_compApp   (x, xx );
    
    fmpq_set_si(r,3,10);
    realApp_set_fmpq(rr, r, prec);
    compApp_set_real_realApp(yy, rr);
    compApp_set_imag_realApp(yy, rr);
    doubCompApp_set_compApp   (y, yy );
    
    printf("1/10 + 1/10i arb: \n"); compApp_print(xx); printf("\n");
    printf("1/10 + 1/10i doub: \n"); doubCompApp_print(x); printf("\n");
    
    compApp_mul(zz,xx,yy,prec);
    doubCompApp_mul(z,x,y);
    
    printf("(1/10 + 1/10i)*(3/10 + 3/10i) arb: \n"); compApp_print(yy); printf("\n");
    printf("(1/10 + 1/10i)*(3/10 + 3/10i) doub: \n"); doubCompApp_print(y); printf("\n");
    
    clock_t start;
    
    start = clock();
    int nbops = 100000000;
    for (int i=0; i<nbops; i++)
//         compApp_mul(zz,xx,yy,prec);
        compApp_mul(zz,xx,xx,prec);
    double time_in_compApp = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    
    printf("time in compApp: %f \n", time_in_compApp);
    
    start = clock();
    for (int i=0; i<nbops; i++)
//         doubCompApp_mul(z,x,y);
        doubCompApp_mul(z,x,x);
//         doubCompApp_sqr(z,x);
    double time_in_doubCompApp = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    
    printf("time in doubCompApp: %f \n", time_in_doubCompApp);
    printf("ratio: %f \n", time_in_compApp/time_in_doubCompApp);
    
    printf("4: "); compApp_print(zz); printf("\n");
    printf("4: "); doubCompApp_print(z); printf("\n");
    
    doubCompApp_get_compApp(yy,z);
    printf("5: "); compApp_print(yy); printf("\n");
    printf("OK: %d\n", compApp_contains(yy,zz) );
    
    return 0;
    
}
