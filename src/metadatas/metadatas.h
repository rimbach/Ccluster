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

#ifndef METADATAS_H
#define METADATAS_H

#ifdef METADATAS_INLINE_C
#define METADATAS_INLINE
#else
#define METADATAS_INLINE static __inline__
#endif

#include "base/base.h"
#include "geometry/compBox.h"
#include "metadatas/strategies.h"
#include "metadatas/counters.h"
#include "metadatas/chronos.h"

#include <string.h>

typedef struct{
    compBox    initB;
    int        verbo;
    strategies strat;
    counters   count;
    chronos    chron;
} metadatas;

typedef metadatas metadatas_t[1];
typedef metadatas * metadatas_ptr;

#define metadatas_initBref(X) (&(X)->initB)
// #define metadatas_verboref(X) (&(X)->verbo)
#define metadatas_stratref(X) (&(X)->strat)
#define metadatas_countref(X) (&(X)->count)
#define metadatas_chronref(X) (&(X)->chron)

void metadatas_init(metadatas_t m, const compBox_t initialBox, const strategies_t strategy, int verbosity);
void metadatas_clear(metadatas_t m);

METADATAS_INLINE int  metadatas_getVerbo(const metadatas_t m) { return m->verbo; }

// /* strategies */
METADATAS_INLINE int metadatas_useNewton         ( const metadatas_t m ) { return strategies_useNewton         (metadatas_stratref(m)); } 
METADATAS_INLINE int metadatas_useTstarOptim     ( const metadatas_t m ) { return strategies_useTstarOptim     (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_usePredictPrec    ( const metadatas_t m ) { return strategies_usePredictPrec    (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useStopWhenCompact( const metadatas_t m ) { return strategies_useStopWhenCompact(metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useAnticipate     ( const metadatas_t m ) { return strategies_useAnticipate     (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useCountSols      ( const metadatas_t m ) { return strategies_useCountSols      (metadatas_stratref(m)); }
// /* counters */
METADATAS_INLINE void metadatas_add_discarded( metadatas_t m, int depth ) { counters_add_discarded(metadatas_countref(m), depth); }
METADATAS_INLINE void metadatas_add_validated( metadatas_t m, int depth, int nbSols ) { counters_add_validated(metadatas_countref(m), depth, nbSols);}
METADATAS_INLINE void metadatas_add_Test     ( metadatas_t m, int depth, int res, int discard, int nbTaylors, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted) {
    counters_add_Test( metadatas_countref(m), depth, res, discard, nbTaylors, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
}
// METADATAS_INLINE void metadatas_add_discarding_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_discarding_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
// METADATAS_INLINE void metadatas_add_validating_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_validating_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
METADATAS_INLINE void metadatas_add_Newton   ( metadatas_t m, int depth, int res ) { counters_add_Newton( metadatas_countref(m), depth, res);}
METADATAS_INLINE void metadatas_count ( metadatas_t m ) { counters_count(metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getDepth( const metadatas_t m) {return counters_getDepth (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbDiscarded                 ( const metadatas_t m ){ return counters_getNbDiscarded                 (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbValidated                 ( const metadatas_t m ){ return counters_getNbValidated                 (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbSolutions                 ( const metadatas_t m ){ return counters_getNbSolutions                 (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbT0Tests                   ( const metadatas_t m ){ return counters_getNbT0Tests                   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbFailingT0Tests            ( const metadatas_t m ){ return counters_getNbFailingT0Tests            (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbGraeffeInT0Tests          ( const metadatas_t m ){ return counters_getNbGraeffeInT0Tests          (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbGraeffeRepetedInT0Tests   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInT0Tests   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTaylorsInT0Tests          ( const metadatas_t m ){ return counters_getNbTaylorsInT0Tests          (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTaylorsRepetedInT0Tests   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInT0Tests   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTSTests                   ( const metadatas_t m ){ return counters_getNbTSTests                   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbFailingTSTests            ( const metadatas_t m ){ return counters_getNbFailingTSTests            (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbGraeffeInTSTests          ( const metadatas_t m ){ return counters_getNbGraeffeInTSTests          (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbGraeffeRepetedInTSTests   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInTSTests   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTaylorsInTSTests          ( const metadatas_t m ){ return counters_getNbTaylorsInTSTests          (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTaylorsRepetedInTSTests   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInTSTests   (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbNewton                    ( const metadatas_t m ){ return counters_getNbNewton                    (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbFailingNewton             ( const metadatas_t m ){ return counters_getNbFailingNewton             (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbEval             ( const metadatas_t m ){ return counters_getNbEval             (metadatas_countref(m));}

// 
// /* chronos */
METADATAS_INLINE double metadatas_get_time_Approxi ( const metadatas_t m ) { return chronos_get_time_Approxi (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Graeffe ( const metadatas_t m ) { return chronos_get_time_Graeffe (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Taylors ( const metadatas_t m ) { return chronos_get_time_Taylors (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_T0Tests ( const metadatas_t m ) { return chronos_get_time_T0Tests (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_TSTests ( const metadatas_t m ) { return chronos_get_time_TSTests (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Newtons ( const metadatas_t m ) { return chronos_get_time_Newtons (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_CclusAl ( const metadatas_t m ) { return chronos_get_time_CclusAl (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Evaluat ( const metadatas_t m ) { return chronos_get_time_Evaluat (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Derivat ( const metadatas_t m ) { return chronos_get_time_Derivat (metadatas_chronref(m)); }

/* printing */
int metadatas_fprint(FILE * file, const metadatas_t meta, const realRat_t eps);
int metadatas_print(const metadatas_t meta, const realRat_t eps);

#endif