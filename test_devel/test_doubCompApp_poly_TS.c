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
#include "numbers/compApp.h"
#include "doubApp/doubCompApp.h"
#include "polynomials/realRat_poly.h"
#include "doubApp/doubCompApp_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include <time.h>

#include <fenv.h>
int main() {
    fesetround(FE_UPWARD);
    
    slong degree = 5;
//     slong prec = 53;
    slong prec = 106;
    compApp_poly_t pres, p;
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
//     bernoulli_polynomial(pbern, degree);
    mignotte_polynomial(pbern, degree, 10);
//     printf("pbern: \n");realRat_poly_print_pretty(pbern, "x");printf("\n\n");
    compApp_poly_init(p);
    compApp_poly_init(pres);
    doubCompApp_poly_t pp, ppres;
    doubCompApp_poly_init(pp);
    doubCompApp_poly_init(ppres);
    
    
    compApp_poly_set_realRat_poly( p, pbern, prec);
    doubCompApp_poly_set_compApp_poly(pp,p);
    
//     printf("p: \n"); compApp_poly_printd(p, prec); printf("\n\n");
    
    compRat_t center;
    realRat_t radius;
    compRat_init(center);
    realRat_init(radius);
    compRat_set_sisi(center, 3,2,5,3);
    realRat_set_si(radius, 1,1);
    compApp_t c;
    realApp_t r;
    compApp_init(c);
    realApp_init(r);
    compApp_set_compRat(c, center, prec);
    realApp_set_realRat(r, radius, prec);
    doubCompApp_t cc;
    doubRealApp_t rr;
    doubRealApp_set_realApp(rr, r);
    doubCompApp_set_compApp(cc, c);
//     printf("p : \n"); doubCompApp_poly_print(pp); printf("\n\n");
    
//     _acb_poly_taylor_shift_convolution(p->coeffs, c, p->length, prec);
//     compApp_poly_taylorShift_in_place( p, center, radius, prec );
    compApp_poly_taylorShift( pres, p, center, radius, prec );
    printf("p shifted in 3/2+i5/3 with radius 1 acb: \n"); compApp_poly_printd(pres, 50); printf("\n\n");
    
//     _doubCompApp_poly_taylor_shift_convolution(pp->coeffs, cc, pp->length);
//     doubCompApp_poly_taylor_shift_horner_inplace( pp, cc, rr);
    doubCompApp_poly_taylor_shift_horner( ppres, pp, cc, rr);
    printf("p shifted in 3/2+i5/3 with radius 1 doub horner: \n"); doubCompApp_poly_print(ppres); printf("\n\n");
    
    doubCompApp_poly_taylor_shift_DQ( ppres, pp, cc, rr);
    printf("p shifted in 3/2+i5/3 with radius 1 doub DQ: \n"); doubCompApp_poly_print(ppres); printf("\n\n");
    
//     doubCompApp_poly_taylor_shift_convolution( ppres, pp, cc, rr);
//     printf("p shifted in 3/2+i5/3 with radius 1 doub convol: \n"); doubCompApp_poly_print(ppres); printf("\n\n");
    

    
    int degs[6] = {31, 63, 127, 255, 511, 1023};
//     int degs[6] = {32, 64, 128, 256, 512, 1024};
    int bitsize = 8;
    int nbtests = 1000;
    clock_t start;
    double time_in_arb, time_in_doub, time_in_doub_DQ;
//     , time_in_doub_convol;
    
    for (int i = 0; i<5; i++) {
        bernoulli_polynomial(pbern, degs[i]);
//         mignotte_polynomial(pbern, degs[i], bitsize);
//         printf("pbern: \n");realRat_poly_print_pretty(pbern, "x");printf("\n\n");
        compApp_poly_set_realRat_poly( p, pbern, prec);
        doubCompApp_poly_set_compApp_poly(pp,p);
//         printf("p: \n"); compApp_poly_printd(p,50); printf("\n\n");
        
        start = clock();
        for (int j = 0; j<nbtests; j++)
            compApp_poly_taylorShift( pres, p, center, radius, prec);
        time_in_arb = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf ("time for %d taylor shifts with arb, degree %d: %lf seconds.\n",nbtests, degs[i], time_in_arb);
        
        start = clock();
        for (int j = 0; j<nbtests; j++)
            doubCompApp_poly_taylor_shift_horner( ppres, pp, cc, rr);
        time_in_doub = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf ("time for %d taylor shifts with doub, degree %d: %lf seconds, ratio: %lf\n",nbtests, degs[i], time_in_doub, time_in_arb/time_in_doub);
        
        start = clock();
        for (int j = 0; j<nbtests; j++)
            doubCompApp_poly_taylor_shift_DQ( ppres, pp, cc, rr);
        time_in_doub_DQ = ((double) (clock() - start))/ CLOCKS_PER_SEC;
        printf ("time for %d taylor shifts DQ with doub, degree %d: %lf seconds, ratio: %lf\n",nbtests, degs[i], time_in_doub_DQ, time_in_arb/time_in_doub_DQ);
        
//         start = clock();
//         for (int j = 0; j<nbtests; j++)
//             doubCompApp_poly_taylor_shift_convolution( ppres, pp, cc, rr);
//         time_in_doub_convol = ((double) (clock() - start))/ CLOCKS_PER_SEC;
//         printf ("time for %d taylor shifts convol with doub, degree %d: %lf seconds, ratio: %lf\n",nbtests, degs[i], time_in_doub_convol, time_in_arb/time_in_doub_convol);
    }
    
    compRat_clear(center);
    realRat_clear(radius);
    compApp_poly_clear(p);
    compApp_poly_clear(pres);
    realRat_poly_clear(pbern);
    doubCompApp_poly_clear(ppres);
    doubCompApp_poly_clear(pp);
    compApp_clear(c);
    doubCompApp_clear(cc);
    realApp_clear(r);
    doubRealApp_clear(rr);
    return 0;
}
