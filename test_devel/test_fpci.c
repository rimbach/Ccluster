#include <stdio.h>
#include "fpri/fpri.h"
#include "fpri/fpci.h"
#include <time.h>

int main( ){
    fpri_lib_init();
    
    fpci_t a;
    fpci_init(a);
    printf("un-initialized fpci: "); fpci_print(a); printf("\n");
    fpci_zero(a);
    printf("zero fpci: "); fpci_print(a); printf("\n");
    fpci_one(a);
    printf("one fpci: "); fpci_print(a); printf("\n");
    fpci_onei(a);
    printf("onei fpci: "); fpci_print(a); printf("\n");
    
    slong prec = 53;
    acb_t aa;
    acb_init(aa);
    acb_unit_root(aa, 21, prec);
    printf("acb 21-th root of unit: "); acb_printd(aa, 16); printf("\n");
    fpci_set_acb(a, aa);
    printf("acb 21-th root of unit converted: "); fpci_print(a); printf("\n");
    acb_sqr(aa, aa, prec);
    printf("acb  (21-th root of unit)^2: "); acb_printd(aa, 16); printf("\n");
    fpci_sqr(a, a);
    printf("fpci (21-th root of unit)^2: "); fpci_print(a); printf("\n");
    acb_clear(aa);
    fpci_clear(a);
    
    fpci_t b, c, d;
    fpci_init(b);
    fpci_init(c);
    fpci_init(d);
    
    printf("\ntest mult: \n");
    fpci_set_si_si(b, 3, 0);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(d, b, c);
    printf("(3)*(5+i6): should be (15+i18): "); fpci_print(d); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 0);
    fpci_mul(d, b, c);
    printf("(3+i4)*(5): should be (15+i20): "); fpci_print(d); printf("\n");
    fpci_set_si_si(b, 0, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(d, b, c);
    printf("(i4)*(5+i6): should be (-24+i20): "); fpci_print(d); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 0, 6);
    fpci_mul(d, b, c);
    printf("(3+i4)*(i6): should be (-24+i18): "); fpci_print(d); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(d, b, c);
    printf("(3+i4)*(5+i6): should be (-9+i38): "); fpci_print(d); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(d, b, b);
    printf("(3+i4)*(3+i4): should be (-7+i24): "); fpci_print(d); printf("\n");
    printf("\ntest aliazing in mult: \n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(b, b, b);
    printf("(3+i4)*(3+i4): should be (-7+i24): "); fpci_print(b); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(b, b, c);
    printf("(3+i4)*(5+i6): should be (-9+i38): "); fpci_print(b); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_mul(c, b, c);
    printf("(3+i4)*(5+i6): should be (-9+i38): "); fpci_print(c); printf("\n");
    
    
    printf("\ntest inv: \n");
    acb_t bb, cc, dd;
    acb_init(bb);
    acb_init(cc);
    acb_init(dd);
    
    fpci_set_si_si(b, 3, 0);
//     fpci_set_si_si(c, 5, 6);
    fpci_inv(d, b);
    printf("1/3: should be ~(%f) ", 1./3); fpci_print(d); printf("\n");
    acb_set_si_si(bb, 3, 0);
//     acb_set_si_si(cc, 5, 6);
    acb_inv(dd, bb, prec);
    printf("1/3: in acb "); acb_printd(dd, 10); printf("\n");
    
    fpci_set_si_si(b, 0, 4);
//     fpci_set_si_si(c, 5, 6);
    fpci_inv(d, b);
    printf("1/(i4): should be ~(-i%f) ", 1./4); fpci_print(d); printf("\n");
    acb_set_si_si(bb, 0, 4);
//     acb_set_si_si(cc, 5, 6);
    acb_inv(dd, bb, prec);
    printf("1/(i4): in acb "); acb_printd(dd, 10); printf("\n");
    
    fpci_set_si_si(b, 3, 4);
//     fpci_set_si_si(c, 5, 6);
    fpci_inv(d, b);
    printf("1/(3+i4): should be (3-i4)/25: ~(%f -i%f) ", 3./25, 4./25); fpci_print(d); printf("\n");
    acb_set_si_si(bb, 3, 4);
//     acb_set_si_si(cc, 5, 6);
    acb_inv(dd, bb, prec);
    printf("1/(3+i4): in acb "); acb_printd(dd, 10); printf("\n");
    
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_div(d, b, c);
    printf("(3+i4)/(5+i6) in fpci: "); fpci_print(d); printf("\n");
    acb_set_si_si(bb, 3, 4);
    acb_set_si_si(cc, 5, 6);
    acb_div(dd, bb, cc, prec);
    printf("(3+i4)/(5+i6) in  acb: "); acb_printd(dd, 10); printf("\n");
    
    printf("\ntest aliazing in div: \n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_div(b, b, b);
    printf("(3+i4)/(3+i4): should be 1: "); fpci_print(b); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_div(b, b, c);
    printf("(3+i4)/(5+i6): "); fpci_print(b); printf("\n");
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_div(c, b, c);
    printf("(3+i4)/(5+i6): "); fpci_print(c); printf("\n");
    
    fpci_set_si_si(b, 3, 4);
    fpci_set_si_si(c, 5, 6);
    fpci_div(d, b, c);
    fpci_sqrabs(fpci_realref(d), d);
    printf("|(3+i4)/(5+i6)|^2 in fpci: "); fpri_print(fpci_realref(d)); printf("\n");
    acb_set_si_si(bb, 3, 4);
    acb_set_si_si(cc, 5, 6);
    acb_div(dd, bb, cc, prec);
    acb_abs(acb_realref(dd), dd, prec);
    arb_sqr(acb_realref(dd), acb_realref(dd), prec );
    printf("|(3+i4)/(5+i6)|^2 in  acb: "); arb_printd(acb_realref(dd), 10); printf("\n");
    
    
    acb_clear(bb);
    acb_clear(cc);
    acb_clear(dd);
    
    fpci_clear(b);
    fpci_clear(c);
    fpci_clear(d);
    
    fpci_t e, f;
    fpci_init(e);
    fpci_init(f);
    printf("\ntest pow_ui: \n");
    acb_t ee, ff;
    acb_init(ee);
    acb_init(ff);
    
    slong order = 21;
    slong pow = 0;
    acb_unit_root(ee, order, prec);
    fpci_set_acb(e, ee);
    
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = 1;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = 2;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = 3;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = 4;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = 5;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    pow = order;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    slong nbops = 10000;
    clock_t start;
    
    start = clock();
    acb_unit_root(ee, order, prec);
    for (slong i=1; i<=nbops; i++) {
//         acb_unit_root(ee, nbops, prec);
//         acb_pow_si(ff, ee, nbops+1, prec);
        
//         acb_unit_root(ee, order, prec);
        acb_pow_si(ff, ee, 1024, prec);
    }
    double time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in powering with acb : %f \n", time_in_acb);
    printf("acb  (%ld-th root of unit)^%ld: ", nbops, nbops+1); acb_printd(ff, 16); printf("\n");
    
    start = clock();
    acb_unit_root(ee, order, prec);
    fpci_set_acb(e, ee);
    for (slong i=1; i<=nbops; i++) {
//         acb_unit_root(ee, nbops, prec);
//         fpci_set_acb(e, ee);
//         fpci_pow_ui(f, e, nbops+1);
        
//         acb_unit_root(ee, order, prec);
//         fpci_set_acb(e, ee);
        fpci_pow_ui(f, e, 1024);
    }
    double time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in powering with fpci: %f \n", time_in_fpci);
    printf("fpci (%ld-th root of unit)^%ld: ", nbops, nbops+1); fpci_print(f); printf("\n");
//     fpri_print(fpci_realref(f)); printf("\n");
//     fpri_print(fpci_imagref(f)); printf("\n");
    
    fpci_t g;
    fpci_init(g);
    acb_t gg;
    acb_init(gg);
    
    acb_unit_root(ee, nbops, prec);
    acb_pow_si(ff, ee, nbops+1, prec);
    fpci_set_acb(e, ee);
    fpci_set_acb(f, ff);
        
    start = clock();
    for (slong i=1; i<=nbops; i++) {
        acb_zero(gg);
        acb_div(gg, ff, ee, prec);
    }
    time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in division with acb : %f \n", time_in_acb);
    printf("acb  (%ld-th root of unit)/(%ld-th root of unit)^%ld: ", nbops, nbops, nbops+1); acb_printd(ff, 16); printf("\n");
    
    start = clock();
    for (slong i=1; i<=nbops; i++) {
        fpci_zero(g);
        fpci_div(g, f, e);
    }
    time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in division with fpci: %f \n", time_in_fpci);
    printf("fpci (%ld-th root of unit)/(%ld-th root of unit)^%ld: ", nbops, nbops, nbops+1); fpci_print(f); printf("\n");
    
    acb_clear(ee);
    acb_clear(ff);
    acb_clear(gg);
    fpci_clear(e);
    fpci_clear(f);
    fpci_clear(g);
    
    
    return 0;
}
