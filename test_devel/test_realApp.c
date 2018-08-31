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

#include "numbers/realApp.h"

int main() {
    
    realApp_t x, y, z;
    realApp_init(x);
    realApp_init(y);
    realApp_init(z);
    
    printf("--- set: ---\n");
    realApp_zero(x);
    printf("zero: "); realApp_print(x); printf("\n");
    printf("zero: "); realApp_printd(x,6); printf("\n");
    
    realApp_one(x);
    printf("one: "); realApp_print(x); printf("\n");
    printf("one: "); realApp_printd(x,6); printf("\n");
    realApp_set(y, x);
    printf("one: "); realApp_print(y); printf("\n");
    printf("one: "); realApp_printd(y,6); printf("\n");
    
    
    fmpq_t r1,r2;
    fmpq_init(r1);
    fmpq_set_si(r1,2,3);
    fmpq_init(r2);
    fmpq_set_si(r2,1,3);
    slong prec = 10;
    
    printf("--- set_fmpq: ---\n");
    realApp_set_fmpq(x, r1, prec);
    printf("2/3: "); realApp_print(x); printf("\n");
    printf("2/3: "); realApp_printd(x,6); printf("\n");
    realApp_set_fmpq(y, r2, prec);
    printf("1/3: "); realApp_print(y); printf("\n");
    printf("1/3: "); realApp_printd(y,6); printf("\n");
    
    printf("--- comparisons: ---\n");
    printf("eq 2/3 1/3? %d\n", realApp_eq(x,y));
    printf("ne 2/3 1/3? %d\n", realApp_ne(x,y));
    printf("lt 2/3 1/3? %d\n", realApp_lt(x,y));
    printf("le 2/3 1/3? %d\n", realApp_le(x,y));
    printf("gt 2/3 1/3? %d\n", realApp_gt(x,y));
    printf("ge 2/3 1/3? %d\n", realApp_ge(x,y));
    
    printf("--- arithmetic operations: ---\n");
    realApp_add(z, x, y, prec);
    printf(" 2/3 + 1/3: "); realApp_printd(z, 6); printf("\n");
    realApp_sub(z, x, y, prec);
    printf(" 2/3 - 1/3: "); realApp_printd(z, 6); printf("\n");
    
    realApp_clear(x);
    realApp_clear(y);
    realApp_clear(z);
    fmpq_clear(r1);
    fmpq_clear(r2);
    
    return 0;
}