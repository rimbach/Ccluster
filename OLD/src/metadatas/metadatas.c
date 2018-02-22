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
//     printf("strategy: %d, %d, %d, %d, %d, %d \n", 
//            metadatas_useNewton(m), 
//            metadatas_useTstarOptim(m), 
//            metadatas_usePredictPrec(m), 
//            metadatas_useStopWhenCompact(m),
//            metadatas_useAnticipate(m),
//            metadatas_useCountSols(m)
//           );
}

void metadatas_clear(metadatas_t m) {
    compBox_clear( metadatas_initBref(m) );
    strategies_clear( metadatas_stratref(m) );
    counters_clear( metadatas_countref(m) );
    chronos_clear( metadatas_chronref(m) );
}

/* printing */

char * compBox_sprint_for_stat(char * out, const compBox_t x){
    char cRe[100];
    char cIm[100];
    char wid[100];
    fmpq_get_str( cRe, 10, compRat_realref(compBox_centerref(x)));
    fmpq_get_str( cIm, 10, compRat_imagref(compBox_centerref(x)));
    fmpq_get_str( wid, 10, compBox_bwidthref(x));
    sprintf(out, " cRe: %-16s cIm: %-16s wid: %-15s|", cRe, cIm, wid);
    return out;
}

char * realRat_sprint_for_stat(char * out, const realRat_t x){
    char * temp = NULL;
    temp = fmpq_get_str( temp, 10, x);
    if (strlen(temp)<=63) {
        sprintf(out, "%s", temp);
    }
    else {
        sprintf(out, "too long representation...");
    }
//     free(temp);
    return out;
}

int metadatas_fprint(FILE * file, const metadatas_t meta, const realRat_t eps){
    int r=1;
    int nbTaylorShifts  = metadatas_getNbT0Tests(meta) + metadatas_getNbTSTests(meta);
    int nbTaylorShiftsR = metadatas_getNbTaylorsRepetedInT0Tests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta);
    int nbGraeffe       = metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeInTSTests(meta);
    int nbGraeffeR      = metadatas_getNbGraeffeRepetedInT0Tests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta);
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, " -------------------Ccluster: ----------------------------------------\n");
    r = fprintf(file, " -------------------Input:    ----------------------------------------\n");
    char temp[100];
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
    if (metadatas_useCountSols(meta)) len += sprintf( temp + len, " count");
    if (metadatas_stratref(meta)->_additionalFlags !=0) 
        len += sprintf(temp +len, " %d", metadatas_stratref(meta)->_additionalFlags);
    r = fprintf(file, "|strat:%-63s|\n", temp);
    
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, " -------------------TSTest used to discard boxes----------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number:",                       metadatas_getNbT0Tests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingT0Tests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests:",          metadatas_get_time_T0Tests(meta),    " " );
    r = fprintf(file, " -------------------TSTest used to validate clusters------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number:",                       metadatas_getNbTSTests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingTSTests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests:",          metadatas_get_time_TSTests(meta),    " " );
    r = fprintf(file, " -------------------Taylor shifts-------------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14d|\n", "total number:",                       nbTaylorShifts + nbTaylorShiftsR, nbTaylorShiftsR );
    r = fprintf(file, "|%-39s %14d %14d|\n", "number in discarding TSTests:",       metadatas_getNbT0Tests(meta) + metadatas_getNbTaylorsRepetedInT0Tests(meta), metadatas_getNbTaylorsRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d %14d|\n", "number in validating TSTests:",       metadatas_getNbTSTests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta), metadatas_getNbTaylorsRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Taylor shifts:",  metadatas_get_time_Taylors(meta),    " " );
    r = fprintf(file, " -------------------Graeffe Iterations--------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14d|\n", "total number:",                       nbGraeffe + nbGraeffeR, nbGraeffeR );
    r = fprintf(file, "|%-39s %14d %14d|\n", "number in discarding TSTests:",       metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeRepetedInT0Tests(meta), metadatas_getNbGraeffeRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d %14d|\n", "number in validating TSTests:",       metadatas_getNbGraeffeInTSTests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta), metadatas_getNbGraeffeRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Graeffe Iterations:", metadatas_get_time_Graeffe(meta),    " " );
    if (metadatas_useNewton(meta)){
    r = fprintf(file, " -------------------Newton Iterations---------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number:",                       metadatas_getNbNewton(meta),         " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of fails:",                    metadatas_getNbFailingNewton(meta),  " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in newton:",         metadatas_get_time_Newtons(meta),    " " );
    }
    r = fprintf(file, " -------------------Other---------------------------------------------\n");
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in getApproximation:",           metadatas_get_time_Approxi(meta),    " " );
    }
    
    r = fprintf(file, " -------------------Output:   ----------------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of clusters:",                 metadatas_getNbValidated(meta),      " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of solutions:",                metadatas_getNbSolutions(meta),      " " );
    r = fprintf(file, " -------------------Stats:    ----------------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree depth:",                         metadatas_getDepth(meta),            " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree size:",                          metadatas_getNbT0Tests(meta),        " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time:",                         metadatas_get_time_CclusAl(meta),    " " );
    r = fprintf(file, " ---------------------------------------------------------------------\n");
    }
    return r;
}

int metadatas_print(const metadatas_t meta, const realRat_t eps) {
    return metadatas_fprint(stdout, meta, eps);
}

void metadatas_getInitB_forJulia(compBox_t b, metadatas_t m) {
    compBox_set(b, metadatas_initBref(m));
}

int  metadatas_getVerbo_forJulia(metadatas_t m){
    return m->verbo;
}

int metadatas_useNewton_forJulia         ( const metadatas_t m ) { return strategies_useNewton         (metadatas_stratref(m)); } 
int metadatas_useTstarOptim_forJulia     ( const metadatas_t m ) { return strategies_useTstarOptim     (metadatas_stratref(m)); }
int metadatas_usePredictPrec_forJulia    ( const metadatas_t m ) { return strategies_usePredictPrec    (metadatas_stratref(m)); }
int metadatas_useStopWhenCompact_forJulia( const metadatas_t m ) { return strategies_useStopWhenCompact(metadatas_stratref(m)); }
int metadatas_useAnticipate_forJulia     ( const metadatas_t m ) { return strategies_useAnticipate     (metadatas_stratref(m)); }
int metadatas_useCountSols_forJulia      ( const metadatas_t m ) { return strategies_useCountSols      (metadatas_stratref(m)); }

void metadatas_add_discarded_forJulia( metadatas_t m, int depth ) { return counters_add_discarded(metadatas_countref(m), depth);}
void metadatas_add_validated_forJulia( metadatas_t m, int depth, int nbSols ) { return counters_add_validated(metadatas_countref(m), depth, nbSols);}
void metadatas_add_Test_forJulia     ( metadatas_t m, int depth, int res, int discard, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted) {
    return counters_add_Test( metadatas_countref(m), depth, res, discard, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
}
// void metadatas_add_discarding_test_forJulia( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_discarding_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
// void metadatas_add_validating_test_forJulia( metadatas_t m, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted){
//     return counters_add_validating_test( metadatas_countref(m), depth, res, nbTaylorsRepeted, nbGraeffe, nbGraeffeRepeted);
// }
void metadatas_add_Newton_forJulia   ( metadatas_t m, int depth, int res ) { return counters_add_Newton( metadatas_countref(m), depth, res);}
void metadatas_count_forJulia ( metadatas_t m ) { counters_count(metadatas_countref(m));}
int  metadatas_getDepth_forJulia( const metadatas_t m) {return counters_getDepth (metadatas_countref(m));}
int  metadatas_getNbDiscarded_forJulia                 ( const metadatas_t m ){ return counters_getNbDiscarded                 (metadatas_countref(m));}
int  metadatas_getNbValidated_forJulia                 ( const metadatas_t m ){ return counters_getNbValidated                 (metadatas_countref(m));}
int  metadatas_getNbSolutions_forJulia                 ( const metadatas_t m ){ return counters_getNbSolutions                 (metadatas_countref(m));}
int  metadatas_getNbT0Tests_forJulia                   ( const metadatas_t m ){ return counters_getNbT0Tests                   (metadatas_countref(m));}
int  metadatas_getNbFailingT0Tests_forJulia            ( const metadatas_t m ){ return counters_getNbFailingT0Tests            (metadatas_countref(m));}
int  metadatas_getNbGraeffeInT0Tests_forJulia          ( const metadatas_t m ){ return counters_getNbGraeffeInT0Tests          (metadatas_countref(m));}
int  metadatas_getNbGraeffeRepetedInT0Tests_forJulia   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInT0Tests   (metadatas_countref(m));}
int  metadatas_getNbTaylorsRepetedInT0Tests_forJulia   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInT0Tests   (metadatas_countref(m));}
int  metadatas_getNbTSTests_forJulia                   ( const metadatas_t m ){ return counters_getNbTSTests                   (metadatas_countref(m));}
int  metadatas_getNbFailingTSTests_forJulia            ( const metadatas_t m ){ return counters_getNbFailingTSTests            (metadatas_countref(m));}
int  metadatas_getNbGraeffeInTSTests_forJulia          ( const metadatas_t m ){ return counters_getNbGraeffeInTSTests          (metadatas_countref(m));}
int  metadatas_getNbGraeffeRepetedInTSTests_forJulia   ( const metadatas_t m ){ return counters_getNbGraeffeRepetedInTSTests   (metadatas_countref(m));}
int  metadatas_getNbTaylorsRepetedInTSTests_forJulia   ( const metadatas_t m ){ return counters_getNbTaylorsRepetedInTSTests   (metadatas_countref(m));}
int  metadatas_getNbNewton_forJulia                    ( const metadatas_t m ){ return counters_getNbNewton                    (metadatas_countref(m));}
int  metadatas_getNbFailingNewton_forJulia             ( const metadatas_t m ){ return counters_getNbFailingNewton             (metadatas_countref(m));}

double metadatas_get_time_Approxi_for_julia ( const metadatas_t m ) { return chronos_get_time_Approxi_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_Graeffe_for_julia ( const metadatas_t m ) { return chronos_get_time_Graeffe_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_Taylors_for_julia ( const metadatas_t m ) { return chronos_get_time_Taylors_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_T0Tests_for_julia ( const metadatas_t m ) { return chronos_get_time_T0Tests_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_TSTests_for_julia ( const metadatas_t m ) { return chronos_get_time_TSTests_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_Newtons_for_julia ( const metadatas_t m ) { return chronos_get_time_Newtons_for_julia (metadatas_chronref(m)); }
double metadatas_get_time_CclusAl_for_julia ( const metadatas_t m ) { return chronos_get_time_CclusAl_for_julia (metadatas_chronref(m)); }

void metadatas_add_time_Approxi_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_Approxi_for_julia (metadatas_chronref(m),ellapsedTime);} 
void metadatas_add_time_Graeffe_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_Graeffe_for_julia (metadatas_chronref(m),ellapsedTime);} 
void metadatas_add_time_Taylors_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_Taylors_for_julia (metadatas_chronref(m),ellapsedTime);} 
void metadatas_add_time_T0Tests_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_T0Tests_for_julia (metadatas_chronref(m),ellapsedTime);} 
void metadatas_add_time_TSTests_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_TSTests_for_julia (metadatas_chronref(m),ellapsedTime);}
void metadatas_add_time_Newtons_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_Newtons_for_julia (metadatas_chronref(m),ellapsedTime);}
void metadatas_add_time_CclusAl_for_julia ( metadatas_t m, double ellapsedTime ) { return chronos_add_time_CclusAl_for_julia (metadatas_chronref(m),ellapsedTime);}