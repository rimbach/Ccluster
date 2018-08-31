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

#include "numbers/realRat.h"


int main() {
    
    realRat_t r, s;
    realRat_init(r);
    realRat_init(s);
    fmpq_one(r);
    realRat_set(s,r);
    realRat_print(r);
    printf("\n");
    realRat_print(s);
    printf("\n");
    realRat_set_si(r,1,3);
    realRat_print(r);
    printf("\n");
    realRat_print(s);
    printf("\n");
    
    realRat_set_si(r,2,4);
    printf("2/4: "); realRat_print(r); printf("\n");
    realRat_canonicalise(r);
    printf("1/2: "); realRat_print(r); printf("\n");
    realRat_set_si(s,2,1);
    realRat_mul(r,r,s);
    printf("(1/2)*2: "); realRat_print(r); printf("\n");
    
    realRat_clear(r);
    realRat_clear(s);
    return 0;   
}


