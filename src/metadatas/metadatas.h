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
#include "numbers/compApp.h"
#include "geometry/compBox.h"
#include "metadatas/strategies.h"
#include "metadatas/counters.h"
#include "metadatas/chronos.h"

#include "metadatas/pwSuDatas.h"

#include "metadatas/cauchyDatas.h"

#include <string.h>

// #ifdef CCLUSTER_HAVE_PTHREAD
// #include <pthread.h>
// #endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
    compBox    initB;
    int        verbo;
    strategies strat;
    counters   count;
    chronos    chron;
// #ifdef CCLUSTER_HAVE_PTHREAD
//     pthread_mutex_t _mutex;
// #endif
    /* for power sums */
    pwSuDatas  pwSum;
//     slong      nbEvalPoints;
//     slong      nbPowerSums; /* >=1, how many power sums, 
//                                including 0-th, are computed in the discarding test*/
//     void(*evalPoly)(compApp_t, compApp_t, const compApp_t, slong);
//     slong      appPrec;
    int        drSub;
    /* for Risolate */
    realRat    spBnd;
} metadatas;

typedef metadatas metadatas_t[1];
typedef metadatas * metadatas_ptr;

#define metadatas_initBref(X) (&(X)->initB)
// #define metadatas_verboref(X) (&(X)->verbo)
#define metadatas_stratref(X) (&(X)->strat)
#define metadatas_countref(X) (&(X)->count)
#define metadatas_chronref(X) (&(X)->chron)
#define metadatas_pwSumref(X) (&(X)->pwSum)
#define metadatas_spBndref(X) (&(X)->spBnd)

void metadatas_init(metadatas_t m, const compBox_t initialBox, const strategies_t strategy, int verbosity);
void metadatas_clear(metadatas_t m);

// void metadatas_join(metadatas_t m1, const metadatas_t m2);

METADATAS_INLINE void metadatas_lock(metadatas_t m){
//     pthread_mutex_lock (&(m->_mutex));
    counters_lock(metadatas_countref(m));
    chronos_lock(metadatas_chronref(m));
}

METADATAS_INLINE void metadatas_unlock(metadatas_t m){
//     pthread_mutex_unlock (&(m->_mutex));
    counters_unlock(metadatas_countref(m));
    chronos_unlock(metadatas_chronref(m));
}  

METADATAS_INLINE int  metadatas_getVerbo(const metadatas_t m) { return m->verbo; }
METADATAS_INLINE void  metadatas_setVerbo(metadatas_t m, int v) { m->verbo=v; }
METADATAS_INLINE int  metadatas_haveToCount(const metadatas_t m) { return (m->verbo > 1); }

METADATAS_INLINE int  metadatas_getDrSub(const metadatas_t m) { return m->drSub; }
METADATAS_INLINE void  metadatas_setDrSub(metadatas_t m, int v) { m->drSub=v; }

METADATAS_INLINE slong  metadatas_getNbEvalPoints   (const metadatas_t m) { 
    return pwSuDatas_nbPntsEval(metadatas_pwSumref(m)); 
}
METADATAS_INLINE void   metadatas_setNbEvalPoints   (metadatas_t m, slong nbEvalPoints) { 
    pwSuDatas_set_nbPntsEval(metadatas_pwSumref(m), nbEvalPoints); 
}
METADATAS_INLINE slong  metadatas_getNbPowerSums    (const metadatas_t m) { 
    return pwSuDatas_nbPwSuComp(metadatas_pwSumref(m)); 
}
METADATAS_INLINE void   metadatas_setNbPowerSums    (metadatas_t m, slong nbPowerSums) { 
    pwSuDatas_set_nbPwSuComp(metadatas_pwSumref(m), nbPowerSums); 
}
METADATAS_INLINE void   metadatas_setIsoRatio_si    (metadatas_t m, slong num, ulong den) { 
    pwSuDatas_set_isolaRatio_si(metadatas_pwSumref(m), num, den); 
}
METADATAS_INLINE realRat_ptr metadatas_getIsoRatio  (metadatas_t m) { 
    return pwSuDatas_isolaRatio_ptr( metadatas_pwSumref(m));
}
METADATAS_INLINE realRat_ptr metadatas_getWantedPrec  (metadatas_t m) { 
    return pwSuDatas_wantedPrec_ptr( metadatas_pwSumref(m));
}


METADATAS_INLINE realRat_ptr metadatas_getSepBound  (metadatas_t m) { 
    return metadatas_spBndref(m);
}

METADATAS_INLINE void metadatas_setSepBound  (metadatas_t m, const realRat_t sepBound) { 
    realRat_set(metadatas_spBndref(m), sepBound);
}


// METADATAS_INLINE slong  metadatas_getAppPrec(const metadatas_t m) { return m->appPrec; }
// METADATAS_INLINE void   metadatas_setAppPrec(metadatas_t m, slong appPrec) { m->appPrec = appPrec; }

// /* strategies */
METADATAS_INLINE int metadatas_useNewton         ( const metadatas_t m ) { return strategies_useNewton         (metadatas_stratref(m)); } 
METADATAS_INLINE int metadatas_useTstarOptim     ( const metadatas_t m ) { return strategies_useTstarOptim     (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_usePredictPrec    ( const metadatas_t m ) { return strategies_usePredictPrec    (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useAnticipate     ( const metadatas_t m ) { return strategies_useAnticipate     (metadatas_stratref(m)); }
// METADATAS_INLINE int metadatas_useCountSols      ( const metadatas_t m ) { return strategies_useCountSols      (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useNBThreads      ( const metadatas_t m ) { return strategies_useNBThreads      (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useDeflation      ( const metadatas_t m ) { return strategies_useDeflation     (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_useRealCoeffs     ( const metadatas_t m ) { return strategies_useRealCoeffs     (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_usePowerSums      ( const metadatas_t m ) { return strategies_usePowerSums      (metadatas_stratref(m)); }
// METADATAS_INLINE int metadatas_pwSuTest         ( const metadatas_t m ) { return strategies_pwSuTest           (metadatas_stratref(m)); }
METADATAS_INLINE int metadatas_forTests          ( const metadatas_t m ) { return strategies_forTests          (metadatas_stratref(m)); }
// /* counters */
METADATAS_INLINE void metadatas_add_discarded( metadatas_t m, int depth ) { 
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_discarded(metadatas_countref(m), depth);
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_validated( metadatas_t m, int depth, int nbSols ) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_validated(metadatas_countref(m), depth, nbSols);
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_explored ( metadatas_t m, int depth ) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_explored(metadatas_countref(m), depth);
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_PsCountingTest ( metadatas_t m, int depth ) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_PsCountingTest(metadatas_countref(m), depth );
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_Test     ( metadatas_t m, int depth, int res, int discard, int inNewton, 
                                               int nbTaylors, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted, int prec, 
                                               double d) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_Test( metadatas_countref(m), depth, res, discard, inNewton, nbTaylors, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted, prec);
    if (discard)
        chronos_add_time_T0Tests( metadatas_chronref(m), d, metadatas_useNBThreads(m));
    else
        chronos_add_time_TSTests( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}
// METADATAS_INLINE void metadatas_add_discarding_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_discarding_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
// METADATAS_INLINE void metadatas_add_validating_test( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_validating_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
METADATAS_INLINE void metadatas_add_Newton   ( metadatas_t m, int depth, int res, double d ) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_Newton( metadatas_countref(m), depth, res);
    chronos_add_time_Newtons( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_Evals( metadatas_t m, int depth, int nbEvals, double d ) {
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    counters_add_Eval( metadatas_countref(m), nbEvals, depth);
    chronos_add_time_Evaluat( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif    
}

METADATAS_INLINE void metadatas_count ( metadatas_t m ) { counters_count(metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getDepth( const metadatas_t m) {return counters_getDepth (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbDiscarded                 ( const metadatas_t m ){ return counters_getNbDiscarded                 (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbValidated                 ( const metadatas_t m ){ return counters_getNbValidated                 (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbExplored                  ( const metadatas_t m ){ return counters_getNbExplored                 (metadatas_countref(m));}
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
METADATAS_INLINE int  metadatas_getNbTSTestsInNewton           ( const metadatas_t m ){ return counters_getNbTSTestsInNewton             (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbTaylorsInNewton           ( const metadatas_t m ){ return counters_getNbTaylorsInNewton             (metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbGraeffeInNewton           ( const metadatas_t m ){ return counters_getNbGraeffeInNewton             (metadatas_countref(m));}

METADATAS_INLINE int  metadatas_getNbPsCountingTest            ( const metadatas_t m ){ return counters_getNbPsCountingTest(metadatas_countref(m));}
METADATAS_INLINE int  metadatas_getNbEval             ( const metadatas_t m ){ return counters_getNbEval             (metadatas_countref(m));}

METADATAS_INLINE int metadatas_boxes_by_prec_fprint ( FILE * file, const metadatas_t m ) {
    return counters_boxes_by_prec_fprint ( file, metadatas_countref(m) );
}

// 
// /* chronos */
METADATAS_INLINE double metadatas_get_time_Approxi ( const metadatas_t m ) { return chronos_get_time_Approxi (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Graeffe ( const metadatas_t m ) { return chronos_get_time_Graeffe (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Taylors ( const metadatas_t m ) { return chronos_get_time_Taylors (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_T0Tests ( const metadatas_t m ) { return chronos_get_time_T0Tests (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_TSTests ( const metadatas_t m ) { return chronos_get_time_TSTests (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Newtons ( const metadatas_t m ) { return chronos_get_time_Newtons (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_CclusAl ( const metadatas_t m ) { return chronos_get_time_CclusAl (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Derivat ( const metadatas_t m ) { return chronos_get_time_Derivat (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Anticip ( const metadatas_t m ) { return chronos_get_time_Anticip (metadatas_chronref(m)); }

METADATAS_INLINE double metadatas_get_time_PSTests ( const metadatas_t m ) { return chronos_get_time_PSTests (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_Evaluat ( const metadatas_t m ) { return chronos_get_time_Evaluat (metadatas_chronref(m)); }

METADATAS_INLINE void metadatas_add_time_Approxi(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
    chronos_add_time_Approxi( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

// METADATAS_INLINE void metadatas_add_time_Newtons(metadatas_t m, double d){
//     chronos_add_time_Newtons( metadatas_chronref(m), d, metadatas_useNBThreads(m));
// }

METADATAS_INLINE void metadatas_add_time_Taylors(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
    chronos_add_time_Taylors( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_Graeffe(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_Graeffe( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

// METADATAS_INLINE void metadatas_add_time_T0Tests(metadatas_t m, double d){
//     chronos_add_time_T0Tests( metadatas_chronref(m), d, metadatas_useNBThreads(m));
// }
// 
// METADATAS_INLINE void metadatas_add_time_TSTests(metadatas_t m, double d){
//     chronos_add_time_TSTests( metadatas_chronref(m), d, metadatas_useNBThreads(m));
// }

METADATAS_INLINE void metadatas_add_time_Anticip(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_Anticip( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_PSTests(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_PSTests( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

// #ifdef CCLUSTER_STATS_PS
// METADATAS_INLINE void metadatas_add_time_PSTestV(metadatas_t m, double d){
// #ifdef CCLUSTER_HAVE_PTHREAD
//                 if (metadatas_useNBThreads(m) >1)
//                     metadatas_lock(m);
// #endif
//     chronos_add_time_PSTestV( metadatas_chronref(m), d, metadatas_useNBThreads(m));
// #ifdef CCLUSTER_HAVE_PTHREAD
//                 if (metadatas_useNBThreads(m) >1)
//                     metadatas_unlock(m);
// #endif
// }
// #endif

METADATAS_INLINE void metadatas_add_time_CclusAl(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_CclusAl( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

// METADATAS_INLINE void metadatas_add_time_Evaluat(metadatas_t m, double d){
// #ifdef CCLUSTER_HAVE_PTHREAD
//                 if (metadatas_useNBThreads(m) >1)
//                     metadatas_unlock(m);
// #endif
//     chronos_add_time_Evaluat( metadatas_chronref(m), d, metadatas_useNBThreads(m));
// #ifdef CCLUSTER_HAVE_PTHREAD
//                 if (metadatas_useNBThreads(m) >1)
//                     metadatas_unlock(m);
// #endif
// }

METADATAS_INLINE void metadatas_add_time_NeTSTes(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_NeTSTes( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefTayl(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefTayl( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefDeri(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefDeri( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefEval(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefEval( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefScal(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefScal( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefGrae(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefGrae( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE void metadatas_add_time_DefTsta(metadatas_t m, double d){
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_lock(m);
#endif
    chronos_add_time_DefTsta( metadatas_chronref(m), d, metadatas_useNBThreads(m));
#ifdef CCLUSTER_HAVE_PTHREAD
                if (metadatas_useNBThreads(m) >1)
                    metadatas_unlock(m);
#endif
}

METADATAS_INLINE double metadatas_get_time_NeTSTes ( const metadatas_t m ) { return chronos_get_time_NeTSTes (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefTayl ( const metadatas_t m ) { return chronos_get_time_DefTayl (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefDeri ( const metadatas_t m ) { return chronos_get_time_DefDeri (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefEval ( const metadatas_t m ) { return chronos_get_time_DefEval (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefScal ( const metadatas_t m ) { return chronos_get_time_DefScal (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefGrae ( const metadatas_t m ) { return chronos_get_time_DefGrae (metadatas_chronref(m)); }
METADATAS_INLINE double metadatas_get_time_DefTsta ( const metadatas_t m ) { return chronos_get_time_DefTsta (metadatas_chronref(m)); }

/* printing */
char * compBox_sprint_for_stat(char * out, const compBox_t x);
char * realRat_sprint_for_stat(char * out, const realRat_t x);

int metadatas_fprint(FILE * file, metadatas_t meta, const realRat_t eps);
int metadatas_print(metadatas_t meta, const realRat_t eps);

int metadatas_risolate_fprint(FILE * file, metadatas_t meta, const realRat_t eps);

#ifdef __cplusplus
}
#endif

#endif
