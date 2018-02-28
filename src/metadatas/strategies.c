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

#include "metadatas/strategies.h"

void strategies_init( strategies_t strat ){}
void strategies_clear(strategies_t strat) {}

void strategies_set_int( strategies_t strat, int useNewton, 
                                             int useTstarOptim, 
                                             int usePredictPrec, 
                                             int useStopWhenCompact, 
                                             int useAnticipate, 
//                                              int useCountSols,
                                             int useNBThreads,
                                             int additionalFlags    
                                             ){
    strat->_useNewton             = useNewton;
    strat->_useTstarOptim         = useTstarOptim;
    strat->_usePredictPrec        = usePredictPrec;
    strat->_useStopWhenCompact    = useStopWhenCompact;
    strat->_useAnticipate         = useAnticipate;
//     strat->_useCountSols          = useCountSols;
    strat->_useNBThreads          = useNBThreads;
    strat->_additionalFlags       = additionalFlags;
}

void strategies_set( strategies_t strat, const strategies_t strat2) {
    strat->_useNewton             = strat2->_useNewton         ;
    strat->_useTstarOptim         = strat2->_useTstarOptim     ;
    strat->_usePredictPrec        = strat2->_usePredictPrec    ;
    strat->_useStopWhenCompact    = strat2->_useStopWhenCompact;
    strat->_useAnticipate         = strat2->_useAnticipate     ;
//     strat->_useCountSols          = strat2->_useCountSols      ;
    strat->_useNBThreads          = strat2->_useNBThreads      ;
    strat->_additionalFlags       = strat2->_additionalFlags   ;
}