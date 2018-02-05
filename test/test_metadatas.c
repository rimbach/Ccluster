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
// #include <stdlib.h>
#include "metadatas/metadatas.h"

int main() {
    
    strategies_t strat;
    strategies_init(strat);
    strategies_set_int( strat, 1, 1, 1, 1, 0, 1);
    
    compBox_t initBox;
    compBox_init(initBox);
    compRat_t c;
    compRat_init(c);
    compRat_set_sisi(c, 1,2,1,3);
    realRat_t w;
    realRat_init(w);
    realRat_set_si(w,  1,4);
    compBox_set_compRat_realRat(initBox, c, w);
    
    metadatas_t mt;
    metadatas_init(mt, initBox, strat, 1);
    
    printf("verbosity         : %d\n", metadatas_getVerbo  ( mt ) );
    printf("useNewton         : %d\n", metadatas_useNewton ( mt ) );
    
    
    metadatas_clear(mt);
    realRat_clear(w);
    compRat_clear(c);
    compBox_clear(initBox);
    strategies_clear(strat);
    return 0;
    
}