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

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int _useNewton;
    int _useTstarOptim;
    int _usePredictPrec;
    int _useStopWhenCompact;
    int _useAnticipate;
//     int _useCountSols;
    int _useNBThreads;
    int _additionalFlags;
    int _realCoeffs;
} strategies;

typedef strategies strategies_t[1];
typedef strategies * strategies_ptr;

void strategies_init( strategies_t strat );
void strategies_set_int ( strategies_t strat, int useNewton, 
                                              int useTstarOptim, 
                                              int usePredictPrec, 
                                              int useStopWhenCompact, 
                                              int useAnticipate, 
                                              int realCoeffs,
//                                               int useCountSols,
                                              int useNBThreads,
                                              int additionalFlags
                        );
void strategies_set( strategies_t strat, const strategies_t strat2);
void strategies_clear(strategies_t strat);
    
METADATAS_INLINE int strategies_useNewton         ( const strategies_t strat ) { return strat->_useNewton         ; } 
METADATAS_INLINE int strategies_useTstarOptim     ( const strategies_t strat ) { return strat->_useTstarOptim     ; }
METADATAS_INLINE int strategies_usePredictPrec    ( const strategies_t strat ) { return strat->_usePredictPrec    ; }
METADATAS_INLINE int strategies_useStopWhenCompact( const strategies_t strat ) { return strat->_useStopWhenCompact; }
METADATAS_INLINE int strategies_useAnticipate     ( const strategies_t strat ) { return strat->_useAnticipate     ; }
// METADATAS_INLINE int strategies_useCountSols      ( const strategies_t strat ) { return strat->_useCountSols      ; }
METADATAS_INLINE int strategies_useNBThreads      ( const strategies_t strat ) { return strat->_useNBThreads      ; }

METADATAS_INLINE int strategies_realCoeffs        ( const strategies_t strat ) { return strat->_realCoeffs      ; }
METADATAS_INLINE void strategies_set_realCoeffs   ( strategies_t strat, int flag ) { strat->_realCoeffs=flag      ; }

#ifdef __cplusplus
}
#endif

#endif
