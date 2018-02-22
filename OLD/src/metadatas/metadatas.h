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

static __inline__ int  metadatas_getVerbo(const metadatas_t m) { return m->verbo; }

// /* strategies */
static __inline__ int metadatas_useNewton         ( const metadatas_t m ) { return strategies_useNewton         (metadatas_stratref(m)); } 
static __inline__ int metadatas_useTstarOptim     ( const metadatas_t m ) { return strategies_useTstarOptim     (metadatas_stratref(m)); }
static __inline__ int metadatas_usePredictPrec    ( const metadatas_t m ) { return strategies_usePredictPrec    (metadatas_stratref(m)); }
static __inline__ int metadatas_useStopWhenCompact( const metadatas_t m ) { return strategies_useStopWhenCompact(metadatas_stratref(m)); }
static __inline__ int metadatas_useAnticipate     ( const metadatas_t m ) { return strategies_useAnticipate     (metadatas_stratref(m)); }
static __inline__ int metadatas_useCountSols      ( const metadatas_t m ) { return strategies_useCountSols      (metadatas_stratref(m)); }
// /* counters */
static __inline__ void metadatas_add_discarded( metadatas_t m, int depth ) { return counters_add_discarded(metadatas_countref(m), depth);}
static __inline__ void metadatas_add_validated( metadatas_t m, int depth, int nbSols ) { return counters_add_validated(metadatas_countref(m), depth, nbSols);}
static __inline__ void metadatas_add_Test     ( metadatas_t m, int depth, int res, int discard, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted) {
    return counters_add_Test( metadatas_countref(m), depth, res, discard, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
}
// static __inline__ void metadatas_add_discarding_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_discarding_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
// static __inline__ void metadatas_add_validating_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_validating_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
static __inline__ void metadatas_add_Newton   ( metadatas_t m, int depth, int res ) { return counters_add_Newton( metadatas_countref(m), depth, res);}
static __inline__ void metadatas_count ( metadatas_t m ) { counters_count(metadatas_countref(m));}
static __inline__ int  metadatas_getDepth( const metadatas_t m) {return counters_getDepth (metadatas_countref(m));}
static __inline__ int  metadatas_getNbDiscarded                 ( const metadatas_t m ){ return counters_getNbDiscarded                 (metadatas_countref(m));}
static __inline__ int  metadatas_getNbValidated                 ( const metadatas_t m ){ return counters_getNbValidated                 (metadatas_countref(m));}
static __inline__ int  metadatas_getNbSolutions                 ( const metadatas_t m ){ return counters_getNbSolutions                 (metadatas_countref(m));}
static __inline__ int  metadatas_getNbT0Tests                   ( const metadatas_t m ){ return counters_getNbT0Tests                   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbFailingT0Tests            ( const metadatas_t m ){ return counters_getNbFailingT0Tests            (metadatas_countref(m));}
static __inline__ int  metadatas_getNbGraeffeInT0Tests          ( const metadatas_t m ){ return counters_getNbGraeffeInT0Tests          (metadatas_countref(m));}
static __inline__ int  metadatas_getNbGraeffeRepetedInT0Tests   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInT0Tests   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbTaylorsRepetedInT0Tests   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInT0Tests   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbTSTests                   ( const metadatas_t m ){ return counters_getNbTSTests                   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbFailingTSTests            ( const metadatas_t m ){ return counters_getNbFailingTSTests            (metadatas_countref(m));}
static __inline__ int  metadatas_getNbGraeffeInTSTests          ( const metadatas_t m ){ return counters_getNbGraeffeInTSTests          (metadatas_countref(m));}
static __inline__ int  metadatas_getNbGraeffeRepetedInTSTests   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInTSTests   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbTaylorsRepetedInTSTests   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInTSTests   (metadatas_countref(m));}
static __inline__ int  metadatas_getNbNewton                    ( const metadatas_t m ){ return counters_getNbNewton                    (metadatas_countref(m));}
static __inline__ int  metadatas_getNbFailingNewton             ( const metadatas_t m ){ return counters_getNbFailingNewton             (metadatas_countref(m));}

// 
// /* chronos */
static __inline__ double metadatas_get_time_Approxi ( const metadatas_t m ) { return chronos_get_time_Approxi (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_Graeffe ( const metadatas_t m ) { return chronos_get_time_Graeffe (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_Taylors ( const metadatas_t m ) { return chronos_get_time_Taylors (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_T0Tests ( const metadatas_t m ) { return chronos_get_time_T0Tests (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_TSTests ( const metadatas_t m ) { return chronos_get_time_TSTests (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_Newtons ( const metadatas_t m ) { return chronos_get_time_Newtons (metadatas_chronref(m)); }
static __inline__ double metadatas_get_time_CclusAl ( const metadatas_t m ) { return chronos_get_time_CclusAl (metadatas_chronref(m)); }

/* printing */
int metadatas_fprint(FILE * file, const metadatas_t meta, const realRat_t eps);
int metadatas_print(const metadatas_t meta, const realRat_t eps);

//for julia...
void metadatas_getInitB_forJulia(compBox_t b, metadatas_t m);
int  metadatas_getVerbo_forJulia(metadatas_t m);

int metadatas_useNewton_forJulia         ( const metadatas_t m ); 
int metadatas_useTstarOptim_forJulia     ( const metadatas_t m );
int metadatas_usePredictPrec_forJulia    ( const metadatas_t m );
int metadatas_useStopWhenCompact_forJulia( const metadatas_t m );
int metadatas_useAnticipate_forJulia     ( const metadatas_t m );
int metadatas_useCountSols_forJulia      ( const metadatas_t m );

void metadatas_add_discarded_forJulia( metadatas_t m, int depth );
void metadatas_add_validated_forJulia( metadatas_t m, int depth, int nbSols );
void metadatas_add_Test_forJulia     ( metadatas_t m, int depth, int res, int discard, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
// void metadatas_add_discarding_test_forJulia( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
// void metadatas_add_validating_test_forJulia( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
void metadatas_add_Newton_forJulia   ( metadatas_t m, int depth, int res );
void metadatas_count_forJulia ( metadatas_t m );
int  metadatas_getDepth_forJulia( const metadatas_t m);
int  metadatas_getNbDiscarded_forJulia                 ( const metadatas_t m );
int  metadatas_getNbValidated_forJulia                 ( const metadatas_t m );
int  metadatas_getNbSolutions_forJulia                 ( const metadatas_t m );
int  metadatas_getNbT0Tests_forJulia                   ( const metadatas_t m );
int  metadatas_getNbFailingT0Tests_forJulia            ( const metadatas_t m );
int  metadatas_getNbGraeffeInT0Tests_forJulia          ( const metadatas_t m );
int  metadatas_getNbGraeffeRepetedInT0Tests_forJulia   ( const metadatas_t m );
int  metadatas_getNbTaylorsRepetedInT0Tests_forJulia   ( const metadatas_t m );
int  metadatas_getNbTSTests_forJulia                   ( const metadatas_t m );
int  metadatas_getNbFailingTSTests_forJulia            ( const metadatas_t m );
int  metadatas_getNbGraeffeInTSTests_forJulia          ( const metadatas_t m );
int  metadatas_getNbGraeffeRepetedInTSTests_forJulia   ( const metadatas_t m );
int  metadatas_getNbTaylorsRepetedInTSTests_forJulia   ( const metadatas_t m );
int  metadatas_getNbNewton_forJulia                    ( const metadatas_t m );
int  metadatas_getNbFailingNewton_forJulia             ( const metadatas_t m );

double metadatas_get_time_Approxi_for_julia ( const metadatas_t m );
double metadatas_get_time_Graeffe_for_julia ( const metadatas_t m );
double metadatas_get_time_Taylors_for_julia ( const metadatas_t m );
double metadatas_get_time_T0Tests_for_julia ( const metadatas_t m );
double metadatas_get_time_TSTests_for_julia ( const metadatas_t m );
double metadatas_get_time_Newtons_for_julia ( const metadatas_t m );
double metadatas_get_time_CclusAl_for_julia ( const metadatas_t m );

void metadatas_add_time_Approxi_for_julia ( metadatas_t m, double ellapsedTime ); 
void metadatas_add_time_Graeffe_for_julia ( metadatas_t m, double ellapsedTime ); 
void metadatas_add_time_Taylors_for_julia ( metadatas_t m, double ellapsedTime ); 
void metadatas_add_time_T0Tests_for_julia ( metadatas_t m, double ellapsedTime ); 
void metadatas_add_time_TSTests_for_julia ( metadatas_t m, double ellapsedTime ); 
void metadatas_add_time_Newtons_for_julia ( metadatas_t m, double ellapsedTime );
void metadatas_add_time_CclusAl_for_julia ( metadatas_t m, double ellapsedTime );
#endif