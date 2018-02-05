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

#include <stdlib.h>

typedef struct {
    int _useNewton;
    int _useTstarOptim;
    int _usePredictPrec;
    int _useStopWhenCompact;
    int _useAnticipate;
    int _useCountSols;
    int _additionalFlags;
} strategies;

typedef strategies strategies_t[1];
typedef strategies * strategies_ptr;

void strategies_init( strategies_t strat );
void strategies_set_int ( strategies_t strat, int useNewton, 
                                              int useTstarOptim, 
                                              int usePredictPrec, 
                                              int useStopWhenCompact, 
                                              int useAnticipate, 
                                              int useCountSols,
                                              int additionalFlags
                        );
void strategies_set( strategies_t strat, const strategies_t strat2);
void strategies_clear(strategies_t strat);
    
static __inline__ int strategies_useNewton         ( const strategies_t strat ) { return strat->_useNewton         ; } 
static __inline__ int strategies_useTstarOptim     ( const strategies_t strat ) { return strat->_useTstarOptim     ; }
static __inline__ int strategies_usePredictPrec    ( const strategies_t strat ) { return strat->_usePredictPrec    ; }
static __inline__ int strategies_useStopWhenCompact( const strategies_t strat ) { return strat->_useStopWhenCompact; }
static __inline__ int strategies_useAnticipate     ( const strategies_t strat ) { return strat->_useAnticipate     ; }
static __inline__ int strategies_useCountSols      ( const strategies_t strat ) { return strat->_useCountSols      ; }

#endif