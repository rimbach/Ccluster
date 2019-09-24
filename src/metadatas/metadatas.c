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

#include "metadatas/metadatas.h"

void metadatas_init(metadatas_t m, const compBox_t initialBox, const strategies_t strategy, int verbosity) {
    compBox_init ( metadatas_initBref(m) );
    compBox_set(metadatas_initBref(m), initialBox);
    m->verbo = verbosity;
    strategies_init( metadatas_stratref(m) );
    strategies_set( metadatas_stratref(m), strategy );
    counters_init( metadatas_countref(m) );
    chronos_init( metadatas_chronref(m) );
}

void metadatas_clear(metadatas_t m) {
    compBox_clear( metadatas_initBref(m) );
    strategies_clear( metadatas_stratref(m) );
    counters_clear( metadatas_countref(m) );
    chronos_clear( metadatas_chronref(m) );
}

// void metadatas_join(metadatas_t m1, const metadatas_t m2){
//     chronos_join( metadatas_chronref(m1), metadatas_chronref(m2) );
//     counters_join( metadatas_countref(m1), metadatas_countref(m2) );
// }

/* printing */

char * compBox_sprint_for_stat(char * out, const compBox_t x){
    char cRe[100];
    char cIm[100];
    char *wid = NULL;
    fmpq_get_str( cRe, 10, compRat_realref(compBox_centerref(x)));
    fmpq_get_str( cIm, 10, compRat_imagref(compBox_centerref(x)));
    wid = fmpq_get_str( wid, 10, compBox_bwidthref(x));
    if (strlen(wid)>=15) {
        slong lnum = fmpz_clog_ui( realRat_numref(compBox_bwidthref(x)), (ulong) 2);
        slong lden = fmpz_clog_ui( realRat_denref(compBox_bwidthref(x)), (ulong) 2);
        sprintf(wid, "2^(%d)/2^(%d)", (int) lnum, (int) lden);
    }
    sprintf(out, " cRe: %-16s cIm: %-16s wid: %-15s|", cRe, cIm, wid);
    ccluster_free(wid);
    return out;
}

char * realRat_sprint_for_stat(char * out, const realRat_t x){
    char *temp = NULL;
    temp = fmpq_get_str( temp, 10, x);
    if (strlen(temp)<=10) {
        sprintf(out, "%s", temp);
    }
    else {
        slong l = fmpz_clog_ui( realRat_denref(x), (ulong) 2);
        fmpz_get_str( temp, 10, realRat_numref(x));
        sprintf(out, "%s/2^(%d)", temp, (int) l);
    }
    ccluster_free(temp);
    return out;
}

int metadatas_fprint(FILE * file, metadatas_t meta, const realRat_t eps){
    int r=1;
    int nbTaylorShifts  = metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsInTSTests(meta);
    int nbTaylorShiftsR = metadatas_getNbTaylorsRepetedInT0Tests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta);
    int nbGraeffe       = metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeInTSTests(meta);
    int nbGraeffeR      = metadatas_getNbGraeffeRepetedInT0Tests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta);
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, " -------------------Ccluster: ----------------------------------------\n");
    r = fprintf(file, " -------------------Input:    ----------------------------------------\n");
    char temp[1000];
    compBox_sprint_for_stat( temp, metadatas_initBref(meta) );
    r = fprintf(file, "|box:%-65s\n", temp);
    realRat_sprint_for_stat( temp, eps );
    r = fprintf(file, "|eps: %-64s|\n", temp);
    int len = 0;
    if (metadatas_useNewton(meta)) len += sprintf( temp + len, " newton");
    if (metadatas_useTstarOptim(meta)) len += sprintf( temp + len, " tstarOpt");
    if (metadatas_usePredictPrec(meta)) len += sprintf( temp + len, " predPrec");
    if (metadatas_useStopWhenCompact(meta)) len += sprintf( temp + len, " stopWhenCompact");
    if (metadatas_useAnticipate(meta)) len += sprintf( temp + len, " anticip");
    if (metadatas_realCoeffs(meta)) len += sprintf( temp + len, " realCoeffs");
    if (metadatas_pwSuTest(meta)) len += sprintf( temp + len, " pwSuTest");
    if (metadatas_forTests(meta)) len += sprintf( temp + len, " test");
#ifdef CCLUSTER_HAVE_PTHREAD
    if (metadatas_useNBThreads(meta)>1) len += sprintf( temp + len, " %d threads", metadatas_useNBThreads(meta));
#endif
    if (metadatas_stratref(meta)->_additionalFlags !=0) 
        len += sprintf(temp +len, " %d", metadatas_stratref(meta)->_additionalFlags);
    r = fprintf(file, "|strat:%-63s|\n", temp);
    
    if (metadatas_getVerbo(meta)>=2) {
//         metadatas_count(meta);
    r = fprintf(file, " -------------------TSTest used to discard boxes----------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number:",                       metadatas_getNbT0Tests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingT0Tests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests:",          metadatas_get_time_T0Tests(meta),    " " );
    r = fprintf(file, " -------------------TSTest used to validate clusters------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number:",                       metadatas_getNbTSTests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingTSTests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests:",          metadatas_get_time_TSTests(meta),    " " );
    r = fprintf(file, " -------------------Taylor shifts-------------------------------------\n");
    r = fprintf(file, "|%-39s %14d |%13d|\n", "total number TS:",                    nbTaylorShifts + nbTaylorShiftsR, nbTaylorShiftsR );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in discarding TSTests TS:",    metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsRepetedInT0Tests(meta), metadatas_getNbTaylorsRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in validating TSTests TS:",    metadatas_getNbTaylorsInTSTests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta), metadatas_getNbTaylorsRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Taylor shifts:",  metadatas_get_time_Taylors(meta),    " " );
    r = fprintf(file, " -------------------Graeffe Iterations--------------------------------\n");
    r = fprintf(file, "|%-39s %14d |%13d|\n", "total number GR:",                       nbGraeffe + nbGraeffeR, nbGraeffeR );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in discarding TSTests GR:",       metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeRepetedInT0Tests(meta), metadatas_getNbGraeffeRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in validating TSTests GR:",       metadatas_getNbGraeffeInTSTests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta), metadatas_getNbGraeffeRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Graeffe Iterations:", metadatas_get_time_Graeffe(meta),    " " );
    if (metadatas_useNewton(meta)){
    r = fprintf(file, " -------------------Newton Iterations---------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number NE:",                       metadatas_getNbNewton(meta),         " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of fails:",                    metadatas_getNbFailingNewton(meta),  " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in newton:",         metadatas_get_time_Newtons(meta),    " " );
    }
    r = fprintf(file, " -------------------Other---------------------------------------------\n");
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in getApproximation:",           metadatas_get_time_Approxi(meta),    " " );
    if (metadatas_useAnticipate(meta)){
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Anticipate:",                 metadatas_get_time_Anticip(meta),    " " );
    }
    if (metadatas_pwSuTest(meta)){
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of Ps counting tests:",  metadatas_getNbPsCountingTest(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests:",          metadatas_get_time_PSTests(meta),    " " );
#ifdef CCLUSTER_STATS_PS
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests V:",        metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests D:",        metadatas_get_time_PSTests(meta)-metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Evaluation:",                 metadatas_get_time_Evaluat(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of evaluations:",        metadatas_getNbEval(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -2:",                 metadatas_getNbM2(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -1:",                 metadatas_getNbM1(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of errors:",             metadatas_getNbEr(meta),    " " );
#endif 
    }
    r = fprintf(file, " -------------------Precision-----------------------------------------\n");
    r = metadatas_boxes_by_prec_fprint ( file, meta );
    
#ifdef CCLUSTER_EXPERIMENTAL    
    if (CCLUSTER_EXP_NUM_T0(meta)||CCLUSTER_EXP_NUM_T1(meta)||CCLUSTER_INC_TEST(meta)) {
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in getDerivative:",              metadatas_get_time_Derivat(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in evaluate:",                   metadatas_get_time_Evaluat(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of evaluations:",              metadatas_getNbEval(meta),    " " );
    }
#endif 
    }
   
    r = fprintf(file, " -------------------Output:   ----------------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of clusters:",                 metadatas_getNbValidated(meta),      " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of solutions:",                metadatas_getNbSolutions(meta),      " " );
    r = fprintf(file, " -------------------Stats:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree depth:",                         metadatas_getDepth(meta),            " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree size:",                          metadatas_getNbExplored(meta),       " " );
    }
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time:",                         metadatas_get_time_CclusAl(meta),    " " );
    r = fprintf(file, " ---------------------------------------------------------------------\n");
    }
    return r;
}

/* DEPRECATED
int metadatas_print(const metadatas_t meta, const realRat_t eps) {
    return metadatas_fprint(stdout, meta, eps);
}
*/
