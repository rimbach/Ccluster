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

#ifndef STRATEGIES_H
#define STRATEGIES_H

#ifdef METADATAS_INLINE_C
#define METADATAS_INLINE
#else
#define METADATAS_INLINE static __inline__
#endif

#include <stdlib.h>
#include <string.h>

/*for isoRatio*/
#include "numbers/realRat.h"

#ifdef __cplusplus
extern "C" {
#endif
    
/* strategies names */
#define STRAT_NEWTON   1
#define STRAT_TSTAROPT 2
#define STRAT_PREDPREC 4
// #define STRAT_STOPWHCO 8
#define STRAT_ANTICIPA 16
#define STRAT_REALCOEF 32
#define STRAT_PWSUTEST 64
#define STRAT_ROOTRADI 128
#define STRAT_SCAANDRO 256
#define STRAT_FORTESTS 512

#define STRAT_STR_ONLSUBD "onlySubd"
#define STRAT_STR_DEFAULT "default"

#define STRAT_STR_V1 "V1"
#define STRAT_STR_V3 "V3"
#define STRAT_STR_V4 "V4"
#define STRAT_STR_V5 "V5"
#define STRAT_STR_PWSUTESTV4 "psV4"
#define STRAT_STR_V6 "V6"
/* CASC 2021 only subd version */
#define STRAT_STR_V7 "V7"
/* CASC 2021 root radii + subd version */
#define STRAT_STR_V8 "V8"

#define STRAT_INT_FORTESTS 247
#define STRAT_STR_FORTESTS "test"

#define STRAT_STR_FORTESTS1 "test1"
#define STRAT_STR_FORTESTS2 "test2"

/* for Cauchy */
#define STRAT_STR_C1 "C1"
#define STRAT_STR_C2 "C2"
#define STRAT_STR_C3 "C3"

typedef struct {
    int _useNewton;
    int _useTstarOptim;
    int _usePredictPrec;
    int _useAnticipate;
//     int _useCountSols;
    int _useNBThreads;
    int _additionalFlags;
    int _useRealCoeffs;
    int _useDeflation;
    int _usePowerSums;
    int _forTests;
//     int _pwSuTest;
    int _pwSuNbPs;
//     realRat _pwSuIsoRatio;
    int _useRootRadii;
    int _useScaAndRou;
    int _useCompression;
    int _usefpri;
} strategies;

typedef strategies strategies_t[1];
typedef strategies * strategies_ptr;

void strategies_init( strategies_t strat );
/* try to get rid of this interface... */
void strategies_set_int ( strategies_t strat, int useNewton, 
                                              int useTstarOptim, 
                                              int usePredictPrec, 
                                              int useAnticipate, 
                                              int useRealCoeffs,
                                              int usePowerSums,
//                                               int useCountSols,
                                              int useNBThreads,
                                              int additionalFlags
                        );

void strategies_set( strategies_t strat, const strategies_t strat2);

void strategies_clear(strategies_t strat);

void strategies_set_str ( strategies_t strat, char * stratName, int nbThreads);
    
METADATAS_INLINE int strategies_useNewton         ( const strategies_t strat ) { return strat->_useNewton         ; } 
METADATAS_INLINE void strategies_set_useNewton         ( strategies_t strat, int flag ) { strat->_useNewton = flag; } 

METADATAS_INLINE int strategies_useTstarOptim     ( const strategies_t strat ) { return strat->_useTstarOptim     ; }
METADATAS_INLINE int strategies_usePredictPrec    ( const strategies_t strat ) { return strat->_usePredictPrec    ; }
METADATAS_INLINE int strategies_useAnticipate     ( const strategies_t strat ) { return strat->_useAnticipate     ; }
// METADATAS_INLINE int strategies_useCountSols      ( const strategies_t strat ) { return strat->_useCountSols      ; }
METADATAS_INLINE int strategies_useNBThreads      ( const strategies_t strat ) { return strat->_useNBThreads      ; }

METADATAS_INLINE int strategies_useRealCoeffs        ( const strategies_t strat ) { return strat->_useRealCoeffs      ; }
METADATAS_INLINE void strategies_set_realCoeffs   ( strategies_t strat, int flag ) { strat->_useRealCoeffs=flag      ; }

METADATAS_INLINE int strategies_useDeflation        ( const strategies_t strat ) { return strat->_useDeflation      ; }
METADATAS_INLINE void strategies_set_Deflation   ( strategies_t strat, int flag ) { strat->_useDeflation=flag      ; }

METADATAS_INLINE int strategies_usePowerSums        ( const strategies_t strat ) { return strat->_usePowerSums      ; }
METADATAS_INLINE void strategies_set_powerSums   ( strategies_t strat, int flag ) { strat->_usePowerSums=flag      ; }

METADATAS_INLINE int strategies_forTests        ( const strategies_t strat ) { return strat->_forTests      ; }
METADATAS_INLINE void strategies_set_forTests   ( strategies_t strat, int flag ) { strat->_forTests=flag      ; }

// METADATAS_INLINE int strategies_pwSuTest        ( const strategies_t strat ) { return strat->_pwSuTest      ; }
// METADATAS_INLINE void strategies_set_pwSuTest   ( strategies_t strat, int flag ) { strat->_pwSuTest=flag      ; }

METADATAS_INLINE int strategies_pwSuNbPs       ( const strategies_t strat ) { return strat->_pwSuNbPs      ; }
METADATAS_INLINE void strategies_set_pwSuNbPs   ( strategies_t strat, int nb ) { strat->_pwSuNbPs=nb      ; }

METADATAS_INLINE int strategies_useRootRadii        ( const strategies_t strat ) { return strat->_useRootRadii      ; }
METADATAS_INLINE void strategies_set_RootRadii   ( strategies_t strat, int flag ) { strat->_useRootRadii=flag      ; }

METADATAS_INLINE int strategies_useScaAndRou        ( const strategies_t strat ) { return strat->_useScaAndRou      ; }
METADATAS_INLINE void strategies_set_ScaAndRou   ( strategies_t strat, int flag ) { strat->_useScaAndRou=flag      ; }

METADATAS_INLINE int strategies_useCompression        ( const strategies_t strat ) { return strat->_useCompression      ; }
METADATAS_INLINE void strategies_set_Compression   ( strategies_t strat, int flag ) { strat->_useCompression=flag      ; }

METADATAS_INLINE int strategies_usefpri        ( const strategies_t strat ) { return strat->_usefpri      ; }
METADATAS_INLINE void strategies_set_fpri   ( strategies_t strat, int flag ) { strat->_usefpri=flag      ; }

#ifdef __cplusplus
}
#endif

#endif
