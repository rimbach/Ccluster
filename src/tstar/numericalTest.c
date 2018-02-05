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

#include "tstar/tstar.h"

int tstar_numerical_test( compApp_poly_t pApprox, const compDsk_t d, slong prec, metadatas_t meta){
    
    int nbsols = -1;
    
    compApp_t ball, pball, dmc, ppball, pppball;
    compApp_t center, pcenter, ppcenter;
    compApp_poly_t pp;
    realRat_t half, nrad;
    compApp_init(ball);
    compApp_init(pball);
    compApp_init(dmc);
    compApp_init(ppball);
    compApp_init(pppball);
    compApp_init(center);
    compApp_init(pcenter);
    compApp_init(ppcenter);
    compApp_poly_init2(pp, pApprox->length);
    realRat_init(half);
    realRat_init(nrad);
    
    /* order 0 evaluation */
    realRat_set_si(half,3,4);
    realRat_mul(nrad, compDsk_radiusref(d), half);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
//     compApp_set_compDsk( ball, compDsk_centerref(d), compDsk_radiusref(d), prec);
    compApp_poly_evaluate(pball, pApprox, ball, prec);
    if (compApp_contains_zero(pball)==0) {
        nbsols = 0;
    }
//     printf("evaluation at order 0: "); compApp_printd(pball,10); printf("\n");
    
    /* order 2 evaluation */
    if (nbsols==-1) {
        compApp_set_compRat(center, compDsk_centerref(d), prec);
        compApp_sub(dmc, ball, center, prec);
        compApp_poly_evaluate2(pcenter,ppcenter, pApprox, center, prec);
        compApp_poly_derivative(pp,pApprox,prec);
        compApp_poly_evaluate2(ppball, pppball, pp, ball, prec);
        realRat_set_si(half,1,2);
        compApp_mul(pball, dmc, dmc, prec);        /*val =                                    (ball-c)^2*/
        compApp_mul_realRat_in_place(pball, half, prec); /*val =                         (1/2)*(ball-c)^2*/
        compApp_mul(pball, pball, pppball, prec); /*val =                         (1/2)*(ball-c)^2*p''(ball)*/
        compApp_mul(ppcenter, ppcenter, dmc, prec); /*valt=        p'(c)*(ball-c) */
        compApp_add(pball, pball,ppcenter, prec); /*val =        p'(c)*(ball-c) + (1/2)*(ball-c)^2*p''(ball)*/
        compApp_add(pball, pball, pcenter, prec); /*val = p(c) + p'(c)*(ball-c) + (1/2)*(ball-c)^2*p''(ball)*/
        if (compApp_contains_zero(pball)==0) {
            nbsols = 0;
        }
        
    }
    
    /* try to show that there is 1 root in the box */
    if ((nbsols==-1)&&CCLUSTER_EXP_NUM_T1(meta)) {
        realApp_t coeff0, coeff1, othercoeffs;
        realApp_init(coeff0);
        realApp_init(coeff1);
        realApp_init(othercoeffs);
        
        compApp_abs( coeff0, pcenter, prec);
        compApp_abs( coeff1, ppcenter, prec);
        compApp_mul(pball, dmc, dmc, prec);        /*val =                                    (ball-c)^2*/
        compApp_mul_realRat_in_place(pball, half, prec); /*val =                         (1/2)*(ball-c)^2*/
        compApp_mul(pball, pball, pppball, prec); /*val =                         (1/2)*(ball-c)^2*p''(ball)*/
        compApp_abs( othercoeffs, pball, prec);
        realApp_add( othercoeffs, coeff0, othercoeffs, prec);
        nbsols = realApp_soft_compare( coeff1, othercoeffs, prec );
        if (nbsols==-2) {
//             printf("NOT ENOUGH PRECISION IN NUMERICAL T1TEST\n");
        }
        if (nbsols==1) nbsols = 1;
        else nbsols = -1;
        
        realApp_clear(coeff0);
        realApp_clear(coeff1);
        realApp_clear(othercoeffs);
    }
    
    compApp_clear(ball);
    compApp_clear(pball);
    compApp_clear(dmc);
    compApp_clear(ppball);
    compApp_clear(pppball);
    compApp_clear(center);
    compApp_clear(pcenter);
    compApp_clear(ppcenter);
    compApp_poly_clear(pp);
    realRat_clear(half);
    realRat_clear(nrad);
    
    return nbsols;
}