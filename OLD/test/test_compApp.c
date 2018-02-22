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
#include "numbers/compApp.h"

int main() {
    compApp_t x, y;
    compApp_init(x);
    compApp_init(y);
    
    printf("--- set: ---\n");
    compApp_zero(x);
    printf("zero: \n"); compApp_print(x); printf("\n");
    compApp_one(x);
    printf("one:  \n"); compApp_print(x); printf("\n");
    compApp_onei(x);
    printf("onei: \n"); compApp_print(x); printf("\n");
    compApp_set(y, x);
    printf("onei: \n"); compApp_print(y); printf("\n");
    
    realApp_t re;
    realApp_init(re);
    compApp_get_imag(re, x);
    compApp_set_real_realApp(y, re);
    printf("one + onei: \n"); compApp_print(y); printf("\n");
    
    compApp_clear(x);
    compApp_clear(y);
    realApp_clear(re);
    return 0;
}