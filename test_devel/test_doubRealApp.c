#include <stdio.h>

#include <time.h>

#include "numbers/realApp.h"
#include "doubApp/doubRealApp.h"

#include <fenv.h>       /* fesetround, FE_* */
// #pragma STDC FENV_ACCESS on
// #include <float.h>
// #include <math.h>

int main() {
    fesetround(FE_UPWARD);
    
    slong prec = 53;
    
    realApp_t x, x2, x3;
    realApp_init(x);
    realApp_init(x2);
    realApp_init(x3);
    
    doubRealApp_t y, y2, y3;
    doubRealApp_one(y);
    printf("1 doub: "); doubRealApp_print(y); printf("\n");
    doubRealApp_zero(y);
    printf("0 doub: "); doubRealApp_print(y); printf("\n");
//     doubRealApp_neg(y,y);
    
    fmpq_t r;
    fmpq_init(r);
    
    
    fmpq_set_si(r,1,1);
    realApp_set_fmpq(x, r, prec);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y, x);
    doubRealApp_set_realApp(y2,x);
    printf("--- test with 1: ---\n");
    printf("1 arb: "); realApp_print(x); printf("\n");
    printf("1 doub: "); doubRealApp_print(y); printf("\n");
    realApp_neg(x3, x);
    doubRealApp_neg(y3,y);
    printf("-1 arb: "); realApp_print(x3); printf("\n");
    printf("-1 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_add(x3, x, x2, prec);
    doubRealApp_add(y3,y,y2);
    printf("1+1 arb: "); realApp_print(x3); printf("\n");
    printf("1+1 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_sub(x3, x, x2, prec);
    doubRealApp_sub(y3,y,y2);
    printf("1-1 arb: "); realApp_print(x3); printf("\n");
    printf("1-1 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf("1*1 arb: "); realApp_print(x3); printf("\n");
    printf("1*1 doub: "); doubRealApp_print(y3); printf("\n");
    
//     fmpq_set_si(r,-1,1);
//     realApp_set_fmpq(x, r, prec);
//     realApp_set_fmpq(x2, r, prec);
//     realApp_add(x3, x, x2, prec);
//     doubRealApp_set_realApp(y,x);
//     doubRealApp_set_realApp(y2,x);
//     doubRealApp_add(y3,y,y2);
//     
//     printf("--- test with -1: ---\n");
//     printf("-1 arb: "); realApp_print(x); printf("\n");
//     printf("-1 doub: "); doubRealApp_print(y); printf("\n");
//     printf("-1+-1 arb: "); realApp_print(x3); printf("\n");
//     printf("-1+-1 doub: "); doubRealApp_print(y3); printf("\n");
//     realApp_sub(x3, x, x2, prec);
//     doubRealApp_sub(y3,y,y2);
//     printf("-1--1 arb: "); realApp_print(x3); printf("\n");
//     printf("-1--1 doub: "); doubRealApp_print(y3); printf("\n");
//     realApp_mul(x3, x, x2, prec);
//     doubRealApp_mul(y3,y,y2);
//     printf("-1*-1 arb: "); realApp_print(x3); printf("\n");
//     printf("-1*-1 doub: "); doubRealApp_print(y3); printf("\n");
//     
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x, r, prec);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y,x);
    doubRealApp_set_realApp(y2,x);
    printf("--- test with 1/10: ---\n");
    printf("1/10: "); realApp_print(x); printf("\n");
    printf("1/10: "); doubRealApp_print(y); printf("\n");
    realApp_neg(x3, x);
    doubRealApp_neg(y3,y);
    printf("-1/10 arb: "); realApp_print(x3); printf("\n");
    printf("-1/10 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_add(x3, x, x2, prec);
    doubRealApp_add(y3,y,y2);
    printf("1/10+1/10: "); realApp_print(x3); printf("\n");
    printf("1/10+1/10: "); doubRealApp_print(y3); printf("\n");
    realApp_sub(x3, x, x2, prec);
    doubRealApp_sub(y3,y,y2);
    printf("1/10-1/10 arb: "); realApp_print(x3); printf("\n");
    printf("1/10-1/10 doub: "); doubRealApp_print(y3); printf("\n");
//     realApp_mul(x3, x, x2, prec);
//     doubRealApp_mul(y3,y,y2);
//     printf("1/10*1/10 arb: "); realApp_print(x3); printf("\n");
//     printf("1/10*1/10 doub: "); doubRealApp_print(y3); printf("\n");
    
    fmpq_set_si(r,3,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    printf("3/10: "); realApp_print(x2); printf("\n");
    printf("3/10: "); doubRealApp_print(y2); printf("\n");
    realApp_neg(x3, x);
    doubRealApp_neg(y3,y);
    printf("-3/10 arb: "); realApp_print(x3); printf("\n");
    printf("-3/10 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_sub(x3, x, x2, prec);
    doubRealApp_sub(y3,y,y2);
    printf("1/10-3/10 arb: "); realApp_print(x3); printf("\n");
    printf("1/10-3/10 doub: "); doubRealApp_print(y3); printf("\n");
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf("1/10*3/10 arb: "); realApp_print(x3); printf("\n");
    printf("1/10*3/10 doub: "); doubRealApp_print(y3); printf("\n");
    
    realApp_sub(x, x, x, prec);
    doubRealApp_sub(y,y,y);
    printf("--- test aliasing: ---\n");
    printf("1/10-1/10 arb: "); realApp_print(x); printf("\n");
    printf("1/10-1/10 doub: "); doubRealApp_print(y); printf("\n");
    realApp_neg(x2, x2);
    doubRealApp_neg(y2,y2);
    printf("-3/10 arb: "); realApp_print(x2); printf("\n");
    printf("-3/10 doub: "); doubRealApp_print(y2); printf("\n");
    
    printf("\n--- test mul: ---\n");
    printf("\n--- low of x >=0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x, r, prec);
    doubRealApp_set_realApp(y,x);
    printf("\n--- low of y >=0 : ---\n");
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- upp of y <=0 : ---\n");
    fmpq_set_si(r,-1,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- case low of y<0 and upp of y > 0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x2, r, prec);
    realApp_sub(x2, x2, x2, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- upp of x <=0 : ---\n");
    fmpq_set_si(r,-1,10);
    realApp_set_fmpq(x, r, prec);
    doubRealApp_set_realApp(y,x);
    printf("\n--- low of y >=0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- upp of y <=0 : ---\n");
    fmpq_set_si(r,-1,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- case low of y<0 and upp of y > 0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x2, r, prec);
    realApp_sub(x2, x2, x2, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- case low of x<0 and upp of x > 0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x, r, prec);
    realApp_sub(x, x, x, prec);
    doubRealApp_set_realApp(y,x);
    printf("\n--- low of y >=0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- upp of y <=0 : ---\n");
    fmpq_set_si(r,-1,10);
    realApp_set_fmpq(x2, r, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
    printf("\n--- case low of y<0 and upp of y > 0 : ---\n");
    fmpq_set_si(r,1,10);
    realApp_set_fmpq(x2, r, prec);
    realApp_sub(x2, x2, x2, prec);
    doubRealApp_set_realApp(y2,x2);
    realApp_mul(x3, x, x2, prec);
    doubRealApp_mul(y3,y,y2);
    printf(" arb: "); realApp_print(x3); printf("\n");
    printf("doub: "); doubRealApp_print(y3); printf("\n");
    doubRealApp_get_realApp(x2,y3);
    printf("contains: %d, %d\n", realApp_contains(x2,x3), realApp_contains(x3,x2) );
    
//     realApp_mul(x3, x, x2, prec);
//     doubRealApp_mul(y3,y,y2);
//     printf("1/10*1/10 arb: "); realApp_print(x3); printf("\n");
//     printf("1/10*1/10 doub: "); doubRealApp_print(y3); printf("\n");
//     
//     
//     printf("3: "); realApp_print(x3); printf("\n");
//     printf("3: "); doubRealApp_print(y3); printf("\n");
//     
//     clock_t start;
//     
//     start = clock();
//     int nbops = 100000000;
//     for (int i=0; i<nbops; i++)
//         realApp_mul(x3, x, x2, prec);
//     double time_in_realApp = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//     
//     printf("time in realApp: %f \n", time_in_realApp);
//     
//     start = clock();
//     for (int i=0; i<nbops; i++)
//         doubRealApp_mul(y3,y,y2);
//     double time_in_doubRealApp = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//     
//     printf("time in doubRealApp: %f \n", time_in_doubRealApp);
//     printf("ratio: %f \n", time_in_realApp/time_in_doubRealApp);
//     
//     printf("4: "); realApp_print(x); printf("\n");
//     printf("4: "); doubRealApp_print(y); printf("\n");
//     
//     doubRealApp_get_realApp(x2,y);
//     printf("5: "); realApp_print(x2); printf("\n");
//     printf("OK: %d\n", realApp_contains(x2,x) );
    
    return 0;
}

// int n_b (char *addr, int i){
//   return (((char) 0x1)& ( ( char ) *(addr + i/8) )>>i%8);
// }
// 
// ulong mantissa( double f ){
// //     ulong res=0;
// //     for (int i=0 ;i<DBL_MANT_DIG;i++) res = res + (n_b( (char*)&f, i )<<i);
// //     return res;
//     ulong res = (((ulong)0x1)<<54) -1;
//     res = res & (*(ulong*)&f);
//     return res;
// }
// 
// ulong exponent( double f ){
//     ulong res=0;
//     for (int i=DBL_MANT_DIG ;i<(64-1);i++) res = res + (n_b( (char*)&f, i )<<(i-DBL_MANT_DIG));
//     return res;
// }
// 
// ulong sign( double f ){
//     return f>0;
// }
// 
// // void base2(int n){
// //   if (n==0) return;
// //   else {
// //     base2(n/2);
// //     printf(" %d ", n%2);
// //     return;
// //   }
// // }
// // 
// 
// // 
// // char d2c (int n){
// //  return ( n<0? '?': (n<10? '0'+n : (n<36 ? 'A' + (n-10) : '?') )); 
// // }
// // 
// // #define S_MAN 23
// // #define S_EXP 8
// // 
// // void mantisse (float f, int result[]) {
// //  int i;
// //  for (i=0 ;i<S_MAN;i++) result[S_MAN -1 -i]=n_b( (char*)&f, i );
// //  return;
// // }
// // 
// // void exposant (float f, int result[]) {
// //  int i;
// //  for (i=0;i<S_EXP;i++) result[S_EXP -1 -i]=n_b( (char*)&f, i + S_MAN );
// //  return;
// // }
// // 
// // int signe (float f) {
// //   return n_b( (char*)&f, S_MAN + S_EXP);
// // }
