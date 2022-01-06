#include <stdio.h>
#include "flint/fmpq_poly.h"
#include "fpri/fpri.h"
#include "fpri/fpri_poly.h"
#include "fpri/fpci.h"
#include "fpri/fpci_poly.h"
#include "polynomials/app_rat_poly.h"
#include <time.h>

void acb_poly_sparse_eval(acb_t y, acb_t z, const acb_poly_t p, slong inNZC[], slong nbNZC, const acb_t point, slong prec);

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
        fmpq_t factor;
        fmpq_init(factor);
        fmpq_set_si(factor, 1, 16);
        
        compApp_poly_scale_realRat_in_place( pacb->coeffs, factor, pacb->length, prec);
        fmpq_clear(factor);
        
        fpci_poly_set_acb_poly(pfpci, pacb);
        
//         printf("acb  input poly: "); acb_poly_printd(pacb,16); printf("\n");
        
//         printf("fpci input poly: "); fpci_poly_print(pfpci); printf("\n");
        
        acb_t xx, yy, zz;
        acb_init(xx);
        acb_init(yy);
        acb_init(zz);
        acb_one(xx);
//         slong logError = -20;
//         arb_add_error_2exp_si(acb_realref(xx), logError);
//         arb_add_error_2exp_si(acb_imagref(xx), logError);
        fpci_t x, y, z;
        fpci_set_acb(x, xx);
//         printf("arb  [1 +/- 2^(%ld)]: ", logError); arb_printd(xx, 16); printf("\n");
//         printf("fpci [1 +/- 2^(%ld)]: ", logError); fpci_print(x); printf("\n");
        slong nbops = 100;
        clock_t start;
        
        
//         start = clock();
//         for (slong i=0; i<nbops; i++){
//             acb_poly_evaluate2_rectangular( yy, zz, pacb, xx, prec);
//         }
//         double time_in_acb_rectangular = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//         
//         start = clock();
//         for (slong i=0; i<nbops; i++){
//             acb_poly_evaluate2_horner( yy, zz, pacb, xx, prec);
//         }
//         double time_in_acb_horner = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//         
//         
//         start = clock();
//         for (slong i=0; i<nbops; i++){
//             fpci_poly_evaluate2_rectangular( y, z, pfpci, x);
//         }
//         double time_in_fpci_rectangular = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//         
//         
//         start = clock();
//         for (slong i=0; i<nbops; i++){
//             fpci_poly_evaluate2_horner( y, z, pfpci, x);
//         }
//         double time_in_fpci_horner = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//         
//         printf("time in acb  rectangular: %f \n", time_in_acb_rectangular);
//         printf("time in acb  horner     : %f \n", time_in_acb_horner);
//         printf("time in fpci rectangular: %f \n", time_in_fpci_rectangular);
//         printf("time in fpci horner     : %f \n", time_in_fpci_horner);
        
        slong order = 32;
        acb_ptr roots = _acb_vec_init(order);
        _acb_vec_unit_roots(roots, order, order, prec);
        for (slong j=0; j<order; j++) {
//             acb_div_si(roots+j, roots+j, 10, prec);
            acb_mul_si(roots+j, roots+j, 2, prec);
//             acb_add(roots+j, roots+j, xx, prec);
        }
        
        slong *inNZC = NULL;
        slong nbNZC = fpci_init_sparse_eval( &inNZC, pfpci );
//         printf("number of non zero coeffs: %ld\n", nbNZC);
        for (slong j=0; j<order; j++) {
            acb_poly_evaluate2_horner( yy, zz, pacb, roots+j, prec);
            printf("acb  horner      : f(x) : "); acb_printd(yy, 16); printf("\n");
            printf("                   f'(x): "); acb_printd(zz, 16); printf("\n");
            acb_poly_evaluate2_rectangular( yy, zz, pacb, roots+j, prec);
            printf("acb  rectangular : f(x) : "); acb_printd(yy, 16); printf("\n");
            printf("                   f'(x): "); acb_printd(zz, 16); printf("\n");
            acb_poly_sparse_eval(yy, zz, pacb, inNZC, nbNZC, roots+j, prec);
            printf("acb  sparse      : f(x) : "); acb_printd(yy, 16); printf("\n");
            printf("                   f'(x): "); acb_printd(zz, 16); printf("\n");
            fpci_set_acb(x, roots+j);
            fpci_poly_evaluate2_horner( y, z, pfpci, x);
            printf("fpci horner      : f(x) : "); fpci_print(y); printf("\n");
            printf("                   f'(x): "); fpci_print(z); printf("\n");
            fpci_poly_evaluate2_rectangular( y, z, pfpci, x);
            printf("fpci rectangular : f(x) : "); fpci_print(y); printf("\n");
            printf("                   f'(x): "); fpci_print(z); printf("\n");
            fpci_poly_evaluate2_rectangulart( y, z, pfpci, x);
            printf("fpci rectangulart: f(x) : "); fpci_print(y); printf("\n");
            printf("                   f'(x): "); fpci_print(z); printf("\n");
            fpci_poly_sparse_eval(y, z, pfpci, inNZC, nbNZC, x);
            printf("fpci sparse      : f(x) : "); fpci_print(y); printf("\n");
            printf("                   f'(x): "); fpci_print(z); printf("\n");
            printf("\n\n");
        }
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                acb_poly_evaluate2_rectangular( yy, zz, pacb, roots+j, prec);
            }
        }
        double time_in_acb_rectangular = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                acb_poly_evaluate2_horner( yy, zz, pacb, roots+j, prec);
            }
        }
        double time_in_acb_horner = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                acb_poly_sparse_eval(yy, zz, pacb, inNZC, nbNZC, roots+j, prec);
            }
        }
        double time_in_acb_sparse = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                fpci_set_acb(x, roots+j);
                fpci_poly_evaluate2_rectangular( y, z, pfpci, x);
            }
        }
        double time_in_fpci_rectangular = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                fpci_set_acb(x, roots+j);
                fpci_poly_evaluate2_rectangulart( y, z, pfpci, x);
            }
        }
        double time_in_fpci_rectangulart = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                fpci_set_acb(x, roots+j);
                fpci_poly_evaluate2_horner( y, z, pfpci, x);
            }
        }
        double time_in_fpci_horner = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        start = clock();
        for (slong i=0; i<nbops; i++){
            for (slong j=0; j<order; j++) {
                fpci_set_acb(x, roots+j);
                fpci_poly_sparse_eval(y, z, pfpci, inNZC, nbNZC, x);
            }
        }
        double time_in_fpci_sparse = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        
        printf("time in acb  horner      : %f \n", time_in_acb_horner);
        printf("time in acb  rectangular : %f \n", time_in_acb_rectangular);
        printf("time in acb  sparse      : %f \n", time_in_acb_sparse);
        printf("time in fpci horner      : %f \n", time_in_fpci_horner);
        printf("time in fpci rectangular : %f \n", time_in_fpci_rectangular);
        printf("time in fpci rectangulart: %f \n", time_in_fpci_rectangulart);
        printf("time in fpci sparse      : %f \n", time_in_fpci_sparse);
        
        fpci_clear_sparse_eval( &inNZC, nbNZC);
        
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

void acb_poly_sparse_eval(acb_t fval, acb_t fderval, const acb_poly_t app, slong inNZC[], slong nbNZC, const acb_t point, slong prec){
    
    acb_set( fval, (app->coeffs) + 0);
    acb_zero( fderval );
    
    acb_t x, xp, mon;
    acb_init(x);
    acb_init(xp);
    acb_init(mon);
    
    slong ind=0;
    if (inNZC[ind] == 0)
        ind++;
    acb_pow_si( xp, point, inNZC[ind]-1, prec );
    acb_mul(x, xp, point, prec);
    
    while ( ind < nbNZC ){
        acb_mul(mon, (app->coeffs) + inNZC[ind], x, prec);
        acb_add(fval, fval, mon, prec);
//         acb_addmul(fval, (app->coeffs) + inNZC[ind], x, prec);
        
        acb_mul(mon, (app->coeffs) + inNZC[ind], xp, prec);
        acb_mul_si( mon, mon, inNZC[ind], prec);
        acb_add(fderval, fderval, mon, prec);
        
        ind++;
        
        if (ind < nbNZC) {
            acb_pow_si(mon, point, inNZC[ind] - inNZC[ind-1], prec);
            acb_mul( x, x, mon, prec);
            acb_mul( xp, xp, mon, prec);
        }
    }
    
    acb_clear(x);
    acb_clear(xp);
    acb_clear(mon);
    
}
