#include <stdio.h>
#include "flint/fmpq_poly.h"
#include "fpri/fpri.h"
#include "fpri/fpri_poly.h"
#include "fpri/fpci.h"
#include "fpri/fpci_poly.h"
#include <time.h>

int main( int argc, char **argv ){
    fpri_lib_init();
    
    fpci_poly_t p;
    fpci_poly_init(p);
    
    printf("zero: \n");
    fpci_poly_zero(p);
    fpci_poly_print(p);
    printf("\n");
    
    printf("one: \n");
    fpci_poly_one(p);
    fpci_poly_print(p);
    printf("\n");
    
    fpci_poly_clear(p);
    
    if (argc>=2) {
        slong prec = 53;
        fmpq_poly_t pfmpq;
        fmpq_poly_init(pfmpq);
        acb_poly_t  pacb;
        acb_poly_init(pacb);
        fpci_poly_t pfpci;
        fpci_poly_init(pfpci);
        
        char * filename = argv[1];
        FILE * curFile  = fopen (filename,"r");
        if (curFile!=NULL) {
            fmpq_poly_fread(curFile, pfmpq);
        }
        fclose (curFile);    
        
        acb_poly_set_fmpq_poly(pacb, pfmpq, prec);
        fpci_poly_set_acb_poly(pfpci, pacb);
        
//         printf("acb  input poly: "); acb_poly_printd(pacb,16); printf("\n");
        
//         printf("fpci input poly: "); fpci_poly_print(pfpci); printf("\n");
        
        acb_t xx, yy, zz;
        acb_init(xx);
        acb_init(yy);
        acb_init(zz);
        acb_one(xx);
        slong logError = -20;
//         arb_add_error_2exp_si(acb_realref(xx), logError);
//         arb_add_error_2exp_si(acb_imagref(xx), logError);
        fpci_t x, y, z;
        fpci_set_acb(x, xx);
//         printf("arb  [1 +/- 2^(%ld)]: ", logError); arb_printd(xx, 16); printf("\n");
//         printf("fpci [1 +/- 2^(%ld)]: ", logError); fpci_print(x); printf("\n");
        slong nbops = 1;
        clock_t start;
        start = clock();
        for (slong i=0; i<nbops; i++){
            acb_poly_evaluate2_horner( yy, zz, pacb, xx, prec);
        }
        double time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in acb : %f \n", time_in_acb);
        start = clock();
        for (slong i=0; i<nbops; i++){
            fpci_poly_evaluate2_horner( y, z, pfpci, x);
        }
        double time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in fpci: %f \n", time_in_fpci);
        printf("acb  values: "); acb_printd(yy, 16); printf(", ");
        acb_printd(zz, 16); printf("\n");
        printf("fpci values: "); fpci_print(y); printf(", ");
        fpci_print(z); printf("\n");
        
        slong order = 32;
        acb_ptr roots = _acb_vec_init(order);
        _acb_vec_unit_roots(roots, order, order, prec);
        for (slong j=0; j<order; j++) {
            acb_div_si(roots+j, roots+j, 10, prec);
//             acb_add(roots+j, roots+j, xx, prec);
        }
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                acb_poly_evaluate2_horner( yy, zz, pacb, roots+j, prec);
                printf("acb  values: "); acb_printd(yy, 16); printf(", ");
                acb_printd(zz, 16); printf("\n");
            }
        }
        time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in acb : %f \n", time_in_acb);
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                fpci_set_acb(x, roots+j);
                fpci_poly_evaluate2_horner( y, z, pfpci, x);
                printf("fpci values: "); fpci_print(y); printf(", ");
                fpci_print(z); printf("\n");
            }
        }
        time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in fpci: %f \n", time_in_fpci);
        
        _acb_vec_clear(roots, order);
        acb_clear(xx);
        acb_clear(yy);
        acb_clear(zz);
        fpci_poly_clear(pfpci);
        acb_poly_clear(pacb);
        fmpq_poly_clear(pfmpq);
    }
    
    flint_cleanup();
    
    return 0;
}
