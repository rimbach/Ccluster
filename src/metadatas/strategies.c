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

void strategies_init( strategies_t strat ){
    realRat_init(strategies_isoRatioRef(strat));
}
void strategies_clear(strategies_t strat) {
    realRat_clear(strategies_isoRatioRef(strat));
}

void strategies_set_int( strategies_t strat, int useNewton, 
                                             int useTstarOptim, 
                                             int usePredictPrec, 
                                             int useStopWhenCompact, 
                                             int useAnticipate, 
                                             int realCoeffs,
                                             int pwSuTest,
//                                              int useCountSols,
                                             int useNBThreads,
                                             int additionalFlags
                                             ){
    strat->_useNewton             = 0;
    strat->_useTstarOptim         = 0;
    strat->_usePredictPrec        = 0;
    strat->_useStopWhenCompact    = 0;
    strat->_useAnticipate         = 0;
    strat->_useNBThreads          = 0;
    strat->_additionalFlags       = 0;
    strat->_realCoeffs            = 0;
    strat->_forTests              = 0;
    strat->_pwSuTest              = 0;
    strat->_pwSuNbPs              = 0;
    
    strat->_useNewton             = useNewton;
    strat->_useTstarOptim         = useTstarOptim;
    strat->_usePredictPrec        = usePredictPrec;
    strat->_useStopWhenCompact    = useStopWhenCompact;
    strat->_useAnticipate         = useAnticipate;
//     strat->_useCountSols          = useCountSols;
    strat->_useNBThreads          = useNBThreads;
    strat->_additionalFlags       = additionalFlags;
    strat->_realCoeffs            = realCoeffs;
    strat->_pwSuTest              = pwSuTest;
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
    strat->_realCoeffs            = strat2->_realCoeffs   ;
    strat->_forTests              = strat2->_forTests   ;
    strat->_pwSuTest              = strat2->_pwSuTest   ;
    strat->_pwSuNbPs           = strat2->_pwSuNbPs;
    realRat_set(strategies_isoRatioRef(strat), strategies_isoRatioRef(strat2));
}

void strategies_set_str ( strategies_t strat, char * stratName, int nbThreads){
    
    strat->_useNewton             = 0;
    strat->_useTstarOptim         = 0;
    strat->_usePredictPrec        = 0;
    strat->_useStopWhenCompact    = 0;
    strat->_useAnticipate         = 0;
    strat->_useNBThreads          = 0;
    strat->_additionalFlags       = 0;
    strat->_realCoeffs            = 0;
    strat->_pwSuTest             = 0;
    strat->_forTests             = 0;
    strat->_pwSuNbPs          = 0;
    
    strat->_useNBThreads          = nbThreads;
    
    if (strcmp( stratName, STRAT_STR_DEFAULT ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        strat->_pwSuTest              = 0;
        strat->_forTests              = 0;
//         printf("default strategy\n");
        return;
    }
    if (strcmp( stratName, STRAT_STR_V1 ) == 0) {
        strat->_useNewton             = 1;
        strat->_usePredictPrec        = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_V3 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_V4 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_V5 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_PWSUTESTV4 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 0;
        strat->_pwSuTest              = 1;
        strat->_forTests              = 0;
        return;
    }
    if (strcmp( stratName, STRAT_STR_V6 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        strat->_pwSuTest              = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_FORTESTS ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        strat->_pwSuTest              = 1;
        strat->_forTests              = 1;
        strat->_pwSuNbPs              = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_FORTESTS1 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        strat->_pwSuTest              = 1;
        strat->_forTests              = 1;
        strat->_pwSuNbPs              = 2;
        return;
    }
    if (strcmp( stratName, STRAT_STR_FORTESTS2 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_realCoeffs            = 1;
        strat->_pwSuTest              = 1;
        strat->_forTests              = 1;
        strat->_pwSuNbPs              = 3;
        return;
    }
    
    /* otherwise set strategy to default */
    strategies_set_str ( strat, STRAT_STR_DEFAULT, nbThreads);
}
