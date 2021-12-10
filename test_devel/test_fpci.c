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
    
    pow = order;
    acb_pow_si(ff, ee, pow, prec);
    fpci_pow_ui(f, e, pow);
    printf("acb  (%ld-th root of unit)^%ld: ", order, pow); acb_printd(ff, 16); printf("\n");
    printf("fpci (%ld-th root of unit)^%ld: ", order, pow); fpci_print(f); printf("\n");
    
    slong nbops = 10000;
    clock_t start;
    
    start = clock();
    for (slong i=1; i<=nbops; i++) {
        acb_unit_root(ee, nbops, prec);
        acb_pow_si(ff, ee, nbops+1, prec);
    }
    double time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in acb : %f \n", time_in_acb);
    printf("acb  (%ld-th root of unit)^%ld: ", nbops, nbops+1); acb_printd(ff, 16); printf("\n");
    
    start = clock();
    for (slong i=1; i<=nbops; i++) {
        acb_unit_root(ee, nbops, prec);
        fpci_set_acb(e, ee);
        fpci_pow_ui(f, e, nbops+1);
    }
    double time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
    printf("time in fpci: %f \n", time_in_fpci);
    printf("fpci (%ld-th root of unit)^%ld: ", nbops, nbops+1); fpci_print(f); printf("\n");
//     fpri_print(fpci_realref(f)); printf("\n");
//     fpri_print(fpci_imagref(f)); printf("\n");
    
    acb_clear(ee);
    acb_clear(ff);
    fpci_clear(e);
    fpci_clear(f);
    
    
    return 0;
}
