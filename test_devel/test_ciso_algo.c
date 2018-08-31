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
#include "polynomials/compRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "cisolate/cisolate.h"

compRat_poly_t p_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

int main(int argc, char **argv){
    
    if (argc<3){
        printf("usage: %s, degree, bitsize\n", argv[0]);
        return -1;
    }
    slong degree;
    slong bitsize;
    sscanf(argv[1], "%d", &degree);
    sscanf(argv[2], "%d", &bitsize);
//     printf("degree: %d, bitsize: %d\n", degree, bitsize);
    
    realRat_poly_t pmign;
    realRat_poly_init(pmign);
    compRat_poly_init(p_global);
    
    mignotte_polynomial(pmign, degree, bitsize);
    compRat_poly_set_realRat_poly(p_global,pmign);
    
    realRat_t eps;
    realRat_init(eps);
    realRat_set_si(eps, 1,100);
    
    compBox_t bInit;
    compBox_init(bInit);
    compBox_set_si(bInit, 0,1,0,1,100,1);
    
    int st = 0;
    st+=(0x1<<0);//newton
    st+=(0x1<<1);//tstar optim
    st+=(0x1<<2);//predict prec
//     st+=(0x1<<3);//stop when compact
//     st+=(0x1<<4);//anticipate
//     st+=(0x1<<5);//count sols
    
    cisolate_interface_func( getApprox, bInit, eps, st, 2);
    
    realRat_poly_clear(pmign);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    return 0;
}
