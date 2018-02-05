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
#include <unistd.h>
#include "metadatas/chronos.h"

int main() {
    
    chronos_t times;
    chronos_init( times );
    
    
    printf("time for graeffe: %f\n", chronos_get_time_Graeffe(times));
    
    chronos_tic_Graeffe(times);
    chronos_toc_Graeffe(times);
    printf("time for graeffe: %f\n", chronos_get_time_Graeffe(times));
    
    double a = 0;
    clock_t t;
    t = clock();
    chronos_tic_Graeffe(times);
    for(int i=1; i<10000; i++)
        for(int j=1; j<10000; j++) {
           a+= (double) i/j;
        }
    printf("a: %f\n", a);
    chronos_toc_Graeffe(times);
    t = clock() - t;
    printf("time for graeffe: %f\n", chronos_get_time_Graeffe(times));
    printf("time for graeffe: %f\n", ( (double) t)/ CLOCKS_PER_SEC );
    
    chronos_clear(times);
    
    return 0;
    
}