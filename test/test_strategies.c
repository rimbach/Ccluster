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
#include <stdlib.h>
#include "metadatas/strategies.h"

int main() {
    
    strategies strat;
    
    strategies_set_int( &strat, 1, 1, 1, 1, 0, 1);
    
    printf("useNewton         : %d\n", strategies_useNewton         ( &strat ) );
    printf("useTstarOptim     : %d\n", strategies_useTstarOptim     ( &strat ) );
    printf("usePredictPrec    : %d\n", strategies_usePredictPrec    ( &strat ) );
    printf("useStopWhenCompact: %d\n", strategies_useStopWhenCompact( &strat ) );
    printf("useAnticipate     : %d\n", strategies_useAnticipate     ( &strat ) );
    printf("useCountSols      : %d\n", strategies_useCountSols      ( &strat ) );
    
    
    return 0;
    
}