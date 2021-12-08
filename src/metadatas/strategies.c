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

void strategies_init( strategies_t strat ){ }
void strategies_clear(strategies_t strat) { }

void strategies_set_int( strategies_t strat, int useNewton, 
                                             int useTstarOptim, 
                                             int usePredictPrec, 
                                             int useAnticipate, 
                                             int useRealCoeffs,
                                             int usePowerSums,
//                                              int useCountSols,
                                             int useNBThreads,
                                             int additionalFlags
                                             ){
    strat->_useNewton             = 0;
    strat->_useTstarOptim         = 0;
    strat->_usePredictPrec        = 0;
    strat->_useAnticipate         = 0;
    strat->_useNBThreads          = 0;
    strat->_additionalFlags       = 0;
    strat->_useRealCoeffs         = 0;
    strat->_usePowerSums          = 0;
    strat->_forTests              = 0;
//     strat->_pwSuTest              = 0;
    strat->_pwSuNbPs              = 0;
    strat->_useRootRadii          = 0;
    strat->_useScaAndRou          =0;
    
    strat->_useNewton             = useNewton;
    strat->_useTstarOptim         = useTstarOptim;
    strat->_usePredictPrec        = usePredictPrec;
    strat->_useAnticipate         = useAnticipate;
//     strat->_useCountSols          = useCountSols;
    strat->_useNBThreads          = useNBThreads;
    strat->_additionalFlags       = additionalFlags;
    strat->_useRealCoeffs            = useRealCoeffs;
    strat->_usePowerSums             = usePowerSums;
//     strat->_pwSuTest              = pwSuTest;
}

void strategies_set( strategies_t strat, const strategies_t strat2) {
    strat->_useNewton             = strat2->_useNewton         ;
    strat->_useTstarOptim         = strat2->_useTstarOptim     ;
    strat->_usePredictPrec        = strat2->_usePredictPrec    ;
    strat->_useAnticipate         = strat2->_useAnticipate     ;
//     strat->_useCountSols          = strat2->_useCountSols      ;
    strat->_useNBThreads          = strat2->_useNBThreads      ;
    strat->_additionalFlags       = strat2->_additionalFlags   ;
    strat->_useRealCoeffs            = strat2->_useRealCoeffs   ;
    strat->_useDeflation          = strat2->_useDeflation       ;
    strat->_usePowerSums             = strat2->_usePowerSums   ;
    strat->_forTests              = strat2->_forTests   ;
//     strat->_pwSuTest              = strat2->_pwSuTest   ;
    strat->_pwSuNbPs           = strat2->_pwSuNbPs;
    strat->_useRootRadii             = strat2->_useRootRadii  ;
    strat->_useScaAndRou             = strat2->_useScaAndRou  ;
    strat->_useCompression             = strat2->_useCompression  ;
}

void strategies_set_str ( strategies_t strat, char * stratName, int nbThreads){
    
    strat->_useNewton             = 0;
    strat->_useTstarOptim         = 0;
    strat->_usePredictPrec        = 0;
    strat->_useAnticipate         = 0;
    strat->_useNBThreads          = 0;
    strat->_additionalFlags       = 0;
    strat->_useRealCoeffs         = 0;
    strat->_useDeflation          = 0;
    strat->_usePowerSums         = 0;
//     strat->_pwSuTest             = 0;
    strat->_forTests             = 0;
    strat->_pwSuNbPs          = 0;
    strat->_useRootRadii         = 0;
    strat->_useScaAndRou         = 0;
    
    strat->_useCompression         = 0;
    
    strat->_useNBThreads          = nbThreads;
    
    if (strcmp( stratName, STRAT_STR_DEFAULT ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useDeflation          = 1;
        strat->_usePowerSums          = 0;
        strat->_forTests              = 0;
        strat->_useRootRadii          = 1;
        strat->_useScaAndRou          = 1;
        return;
    }
    
    if (strcmp( stratName, STRAT_STR_ONLSUBD ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useDeflation          = 1;
        strat->_usePowerSums          = 0;
        strat->_forTests              = 0;
        strat->_useRootRadii          = 0;
        strat->_useScaAndRou          = 1;
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
        strat->_useRealCoeffs         = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_PWSUTESTV4 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 0;
        strat->_usePowerSums          = 1;
        strat->_forTests              = 0;
        return;
    }
    if (strcmp( stratName, STRAT_STR_V6 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_usePowerSums          = 1;
        return;
    }
    if (strcmp( stratName, STRAT_STR_FORTESTS ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useDeflation          = 1;
        strat->_usePowerSums          = 0;
//         strat->_forTests              = 1;
        strat->_pwSuNbPs              = 0;
        return;
    }

    if (strcmp( stratName, STRAT_STR_V7 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useDeflation          = 1;
        strat->_usePowerSums          = 0;
        strat->_forTests              = 0;
        strat->_useRootRadii          = 0;
        strat->_useScaAndRou          = 0;
        return;
    }
    
    if (strcmp( stratName, STRAT_STR_V8 ) == 0) {
        strat->_useNewton             = 1;
        strat->_useTstarOptim         = 1;
        strat->_usePredictPrec        = 1;
        strat->_useAnticipate         = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useDeflation          = 1;
        strat->_usePowerSums          = 0;
        strat->_forTests              = 0;
        strat->_useRootRadii          = 1;
        strat->_useScaAndRou          = 0;
        return;
    }
    
    if (strcmp( stratName, STRAT_STR_C1 ) == 0) {
        strat->_useNewton             = 1;
        strat->_usePredictPrec        = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useCompression        = 0;
        return;
    }
    
    if (strcmp( stratName, STRAT_STR_C2 ) == 0) {
        strat->_useNewton             = 1;
        strat->_usePredictPrec        = 1;
        strat->_useRealCoeffs         = 1;
        strat->_useCompression        = 1;
        return;
    }
    
    /* otherwise set strategy to default */
    strategies_set_str ( strat, STRAT_STR_DEFAULT, nbThreads);
}
