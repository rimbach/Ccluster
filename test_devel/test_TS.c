/* ************************************************************************** */
/*  Copyright (C) 2018 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include <stdio.h>
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include <time.h>

void realRat_poly_wilkinson( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(res, res, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkRat( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, ((ulong) degree)+1);
        realRat_poly_mul(res, res, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkFac( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_t faci, mfaci;
    realRat_init(faci);
    realRat_init(mfaci);
    realRat_set_si(faci, 1,1);
    
    for (int i=1; i<=degree; i++){
        realRat_mul_si(mfaci, faci, -1);
        realRat_poly_set_coeff_realRat(ptemp, 0, mfaci);
        realRat_poly_mul(res, res, ptemp);
        realRat_mul_si(faci, faci, (slong) i+1 );
    }
    
    realRat_clear(faci);
    realRat_clear(mfaci);
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkFacRat( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_t faci, mfaci, facd;
    realRat_init(faci);
    realRat_init(mfaci);
    realRat_init(facd);
    realRat_set_si(faci, 1,1);
    realRat_set_si(facd, 1,1);
    
    for (int i=1; i<degree; i++)
        realRat_mul_si(facd, facd, (slong) i+1 );
    
    for (int i=1; i<=degree; i++){
        realRat_mul_si(mfaci, faci, -1);
        realRat_div(mfaci, mfaci, facd);
        realRat_poly_set_coeff_realRat(ptemp, 0, mfaci);
        realRat_poly_mul(res, res, ptemp);
        realRat_mul_si(faci, faci, (slong) i+1 );
    }
    
    realRat_clear(faci);
    realRat_clear(mfaci);
    realRat_clear(facd);
    realRat_poly_clear(ptemp);
}


void compApp_poly_div_si ( compApp_poly_t res, const compApp_poly_t p, slong y, slong prec) {
    slong len = p->length;
    while (len>0) {
        compApp_div_si( (res->coeffs) + (len-1), (p->coeffs) + (len-1), y, prec);
        len--;
    }
}

void precompute( compApp_poly_ptr tab, const compApp_poly_t p, slong prec ){
    
    slong d = p->length -1;
    slong i;
    compApp_poly_init2( tab, p->length);
    compApp_poly_set( tab, p);
//     printf("0-th derivative of p: \n"); compApp_poly_printd(tab, prec); printf("\n\n");
    for (i=1; i<=d; i++) {
        compApp_poly_init2( tab+i, p->length -i);
        compApp_poly_derivative( tab + i, tab +(i-1), prec);
        compApp_poly_div_si ( tab + i, tab + i, i, prec);
//         printf("%i-th derivative of p: \n", (int) i); compApp_poly_printd(tab + i, prec); printf("\n\n");
    }
    
}

void compApp_poly_taylorShift_home( compApp_poly_t res, compApp_poly_ptr tab, const compApp_t center, slong prec){
    slong lentable = tab->length;
//     printf("len of table: %d\n", (int) lentable);
    compApp_poly_fit_length(res, lentable);
    slong i;
    for (i=0;i<lentable;i++){
//         printf("i: %d\n", (int) i);
//         printf("i-th coeff: "); compApp_printd( (res->coeffs) + i, prec); printf("\n");
        compApp_poly_evaluate( (res->coeffs) + i, tab+i, center, prec);
//         if (i+1<lentable) {
//             compApp_poly_evaluate2( (res->coeffs) + i,(res->coeffs) + i+1, tab+i, center, prec);
//             compApp_div_si( (res->coeffs) + i+1, (res->coeffs) + i+1, i+1, prec);
//             i++;
//         }
//         else {
//             compApp_poly_evaluate( (res->coeffs) + i, tab+i, center, prec);
//         }
    }
    compApp_poly_set_length(res, lentable);
}

void clear_table( compApp_poly_ptr tab ) {
    slong lentable = tab->length;
    slong i;
    for (i=0; i<lentable; i++) compApp_poly_clear(tab + i);
}

int main() {
    
    slong degree = 512;
    slong prec = 500;
    int nbtests = 1000;
    
    compApp_poly_t p, pshift1, pshift2;
    realRat_poly_t pbern;
    compRat_t center;
    compApp_t c;
    compRat_init(center);
    compApp_init(c);
    compApp_poly_init(p);
    compApp_poly_init(pshift1);
    compApp_poly_init(pshift2);
    realRat_poly_init(pbern);
    
    compApp_poly_ptr ptab = (compApp_poly_ptr) malloc ((degree+1)*sizeof(compApp_poly));
    
//     bernoulli_polynomial(pbern, degree);
    mignotte_polynomial(pbern, degree, 200);
//     realRat_poly_wilkinson(pbern, degree);
//     realRat_poly_wilkRat(pbern, degree);
//     realRat_poly_wilkFac(pbern, degree);
//     realRat_poly_wilkRat(pbern, degree);
    
    compApp_poly_set_realRat_poly( p, pbern, prec);
//     printf("p: \n"); compApp_poly_printd(p, prec); printf("\n\n");
    
//     compApp_poly_div_si(p, p, (slong) 2, prec);
//     printf("p/2: \n"); compApp_poly_printd(p, prec); printf("\n\n");
    
    compRat_set_sisi(center, 1,1,1,1);
    compApp_set_compRat(c, center, prec);
    
    clock_t ti;
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_set(pshift1, p);
//         _acb_poly_taylor_shift_horner(pshift1->coeffs, c, pshift1->length, prec);
        _acb_poly_taylor_shift_convolution(pshift1->coeffs, c, pshift1->length, prec);
//         _acb_poly_taylor_shift_divconquer(pshift1->coeffs, c, pshift1->length, prec);
//         printf("pshift1: \n"); compApp_poly_printd(pshift1, prec); printf("\n\n");
    }
    ti = clock() - ti;
    printf ("time for %d convo taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    precompute( ptab, p, prec );
    ti = clock() - ti;
    printf ("time for precomputing derivatives: %f seconds.\n", ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++)
        compApp_poly_taylorShift_home( pshift2, ptab, c, prec);
    ti = clock() - ti;
    printf ("time for %d home taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
//     printf("pshift2: \n"); compApp_poly_printd(pshift2, prec); printf("\n\n");
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_evaluate( (pshift2->coeffs), ptab, c, prec);
        compApp_poly_evaluate( (pshift2->coeffs)+1, ptab+1, c, prec);
    }
    ti = clock() - ti;
    printf ("time for %d rectangular evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_evaluate( (pshift2->coeffs), ptab+degree/2, c, prec);
        compApp_poly_evaluate( (pshift2->coeffs)+1, ptab+degree/2 +1, c, prec);
    }
    ti = clock() - ti;
    printf ("time for %d rectangular evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_evaluate_horner( (pshift2->coeffs), ptab, c, prec);
        compApp_poly_evaluate_horner( (pshift2->coeffs)+1, ptab+1, c, prec);
    }
    ti = clock() - ti;
    printf ("time for %d horner evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_evaluate_horner( (pshift2->coeffs), ptab+degree/2, c, prec);
        compApp_poly_evaluate_horner( (pshift2->coeffs)+1, ptab+degree/2 +1, c, prec);
    }
    ti = clock() - ti;
    printf ("time for %d horner evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    clear_table( ptab );
    
    free(ptab);
    compRat_clear(center);
    compApp_clear(c);
    compApp_poly_clear(p);
    compApp_poly_clear(pshift1);
    compApp_poly_clear(pshift2);
    realRat_poly_clear(pbern);

    
    return 0;
}
