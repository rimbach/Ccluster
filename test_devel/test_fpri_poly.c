#include <stdio.h>
#include "flint/fmpq_poly.h"
#include "fpri/fpri.h"
#include "fpri/fpri_poly.h"
#include <time.h>

int main( int argc, char **argv ){
    fpri_lib_init();
    
    fpri_poly_t p;
    fpri_poly_init(p);
    
    printf("zero: \n");
    fpri_poly_zero(p);
    fpri_poly_print(p);
    printf("\n");
    
    printf("one: \n");
    fpri_poly_one(p);
    fpri_poly_print(p);
    printf("\n");
    
    fpri_poly_clear(p);
    
    if (argc>=2) {
        slong prec = 53;
        fmpq_poly_t pfmpq;
        fmpq_poly_init(pfmpq);
        arb_poly_t  parb;
        arb_poly_init(parb);
        fpri_poly_t pfpri;
        fpri_poly_init(pfpri);
        
        char * filename = argv[1];
        FILE * curFile  = fopen (filename,"r");
        if (curFile!=NULL) {
            fmpq_poly_fread(curFile, pfmpq);
        }
        fclose (curFile);    
        
        arb_poly_set_fmpq_poly(parb, pfmpq, prec);
        fpri_poly_set_arb_poly(pfpri, parb);
        
//         printf("arb  input poly: "); arb_poly_printd(parb,16); printf("\n");
        
//         printf("fpri input poly: "); fpri_poly_print(pfpri); printf("\n");
        
        arb_t xx, yy, zz;
        arb_init(xx);
        arb_init(yy);
        arb_init(zz);
        arb_one(xx);
        slong logError = -20;
        arb_add_error_2exp_si(xx, logError);
        fpri_t x, y, z;
        fpri_set_arb(x, xx);
//         printf("arb  [1 +/- 2^(%ld)]: ", logError); arb_printd(xx, 16); printf("\n");
//         printf("fpri [1 +/- 2^(%ld)]: ", logError); fpri_print(x); printf("\n");
        slong nbops = 10000;
        clock_t start;
        start = clock();
        for (slong i=0; i<nbops; i++){
            arb_poly_evaluate2_horner( yy, zz, parb, xx, prec);
        }
        double time_in_acb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in acb : %f \n", time_in_acb);
        start = clock();
        for (slong i=0; i<nbops; i++){
            fpri_poly_evaluate2_horner( y, z, pfpri, x);
        }
        double time_in_fpci = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf("time in fpci: %f \n", time_in_fpci);
        printf("arb  values: "); arb_printd(yy, 16); printf(", ");
        arb_printd(zz, 16); printf("\n");
        printf("fpri values: "); fpri_print(y); printf(", ");
        fpri_print(z); printf("\n");
        
        arb_clear(xx);
        arb_clear(yy);
        arb_clear(zz);
        fpri_poly_clear(pfpri);
        arb_poly_clear(parb);
        fmpq_poly_clear(pfmpq);
    }
    
    flint_cleanup();
    
    return 0;
}
