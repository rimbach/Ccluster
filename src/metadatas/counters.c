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

#include "metadatas/counters.h"

void counters_by_depth_init( counters_by_depth_t st) {
    st->nbDiscarded               = 0;
    st->nbValidated               = 0;
    st->nbSolutions               = 0;
    st->nbT0Tests                 = 0;
    st->nbFailingT0Tests          = 0;
    st->nbGraeffeInT0Tests        = 0;
    st->nbGraeffeRepetedInT0Tests = 0;
    st->nbTaylorsInT0Tests        = 0;
    st->nbTaylorsRepetedInT0Tests = 0;
    st->nbTSTests                 = 0;
    st->nbFailingTSTests          = 0;
    st->nbGraeffeInTSTests        = 0;
    st->nbGraeffeRepetedInTSTests = 0;
    st->nbTaylorsInTSTests        = 0;
    st->nbTaylorsRepetedInTSTests = 0;
    st->nbNewton                  = 0;
    st->nbFailingNewton           = 0;
    st->nbEval           = 0;
}

void counters_init( counters * st) {
    st->size = 0;
    st->size_allocated = INIT_SIZE_STATS;
    st->table = (counters_by_depth_ptr) malloc (INIT_SIZE_STATS*sizeof(counters_by_depth));
    counters_by_depth_init( st->total );
}

void counters_clear( counters * st) {
    free(st->table);
    counters_by_depth_clear( st->total );
}

void counters_adjust_table( counters_t st, int depth ){
    /* realocate size if needed */
    while (depth+1 > st->size_allocated) {
        st->size_allocated += INIT_SIZE_STATS;
        st->table = (counters_by_depth_ptr) realloc ( (void *) st->table, (st->size_allocated)*sizeof(counters_by_depth) );
    }
    /* set to 0 intermediate cases if needed */
    while (depth+1 > st->size ){
        counters_by_depth_init( st->table + st->size );
        st->size += 1;
    }
}

void counters_add_discarded( counters_t st, int depth ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbDiscarded +=1;
}

void counters_add_validated( counters_t st, int depth, int nbSols ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbValidated +=1;
    (st->table[depth]).nbSolutions +=nbSols;
}

void counters_add_Test     ( counters_t st, int depth, int res, int discard, 
                             int nbTaylors, int nbTaylorsRepeted, 
                             int nbGraeffe, int nbGraeffeRepeted){
    counters_adjust_table(st, depth);
    if (discard) {
        (st->table[depth]).nbT0Tests                           += 1;
        if (!res) (st->table[depth]).nbFailingT0Tests          += 1;
        (st->table[depth]).nbGraeffeInT0Tests                  += nbGraeffe;
        (st->table[depth]).nbGraeffeRepetedInT0Tests           += nbGraeffeRepeted;
        (st->table[depth]).nbTaylorsInT0Tests                           += nbTaylors;
        (st->table[depth]).nbTaylorsRepetedInT0Tests           += nbTaylorsRepeted;
    }
    else {
        (st->table[depth]).nbTSTests                           += 1;
        if (!res) (st->table[depth]).nbFailingTSTests          += 1;
        (st->table[depth]).nbGraeffeInTSTests                  += nbGraeffe;
        (st->table[depth]).nbGraeffeRepetedInTSTests           += nbGraeffeRepeted;
        (st->table[depth]).nbTaylorsInTSTests                           += nbTaylors;
        (st->table[depth]).nbTaylorsRepetedInTSTests           += nbTaylorsRepeted;
    }
        
}

void counters_add_Newton   ( counters_t st, int depth, int res ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbNewton                  +=1;
    if (!res) (st->table[depth]).nbFailingNewton +=1;
}

void counters_add_Eval( counters_t st, int nbEvals, int depth ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbEval                  +=1;
}

void counters_count ( counters_t st ) {
    for (int i = 0; i< st->size; i++) {
       st->total->nbDiscarded               += (st->table)[i].nbDiscarded               ; 
       st->total->nbValidated               += (st->table)[i].nbValidated               ; 
       st->total->nbSolutions               += (st->table)[i].nbSolutions               ;
       st->total->nbT0Tests                 += (st->table)[i].nbT0Tests                 ; 
       st->total->nbFailingT0Tests          += (st->table)[i].nbFailingT0Tests          ; 
       st->total->nbGraeffeInT0Tests        += (st->table)[i].nbGraeffeInT0Tests        ; 
       st->total->nbGraeffeRepetedInT0Tests += (st->table)[i].nbGraeffeRepetedInT0Tests ;
       st->total->nbTaylorsInT0Tests        += (st->table)[i].nbTaylorsInT0Tests        ;
       st->total->nbTaylorsRepetedInT0Tests += (st->table)[i].nbTaylorsRepetedInT0Tests ; 
       st->total->nbTSTests                 += (st->table)[i].nbTSTests                 ; 
       st->total->nbFailingTSTests          += (st->table)[i].nbFailingTSTests          ; 
       st->total->nbGraeffeInTSTests        += (st->table)[i].nbGraeffeInTSTests        ; 
       st->total->nbGraeffeRepetedInTSTests += (st->table)[i].nbGraeffeRepetedInTSTests ;
       st->total->nbTaylorsInTSTests        += (st->table)[i].nbTaylorsInTSTests        ;
       st->total->nbTaylorsRepetedInTSTests += (st->table)[i].nbTaylorsRepetedInTSTests ; 
       st->total->nbNewton                  += (st->table)[i].nbNewton                  ; 
       st->total->nbFailingNewton           += (st->table)[i].nbFailingNewton           ; 
       st->total->nbEval                    += (st->table)[i].nbEval           ;
    }

}

int counters_getDepth( const counters_t st) { return st->size;}

int counters_getNbDiscarded                 ( const counters_t st ){ return st->total->nbDiscarded               ;}
int counters_getNbValidated                 ( const counters_t st ){ return st->total->nbValidated               ;}
int counters_getNbSolutions                 ( const counters_t st ){ return st->total->nbSolutions               ;}
int counters_getNbT0Tests                   ( const counters_t st ){ return st->total->nbT0Tests                 ;}
int counters_getNbFailingT0Tests            ( const counters_t st ){ return st->total->nbFailingT0Tests          ;}
int counters_getNbGraeffeInT0Tests          ( const counters_t st ){ return st->total->nbGraeffeInT0Tests        ;}
int counters_getNbGraeffeRepetedInT0Tests   ( const counters_t st ){ return st->total->nbGraeffeRepetedInT0Tests ;}
int counters_getNbTaylorsInT0Tests          ( const counters_t st ){ return st->total->nbTaylorsInT0Tests        ;}
int counters_getNbTaylorsRepetedInT0Tests   ( const counters_t st ){ return st->total->nbTaylorsRepetedInT0Tests ;}
int counters_getNbTSTests                   ( const counters_t st ){ return st->total->nbTSTests                 ;}
int counters_getNbFailingTSTests            ( const counters_t st ){ return st->total->nbFailingTSTests          ;}
int counters_getNbGraeffeInTSTests          ( const counters_t st ){ return st->total->nbGraeffeInTSTests        ;}
int counters_getNbGraeffeRepetedInTSTests   ( const counters_t st ){ return st->total->nbGraeffeRepetedInTSTests ;}
int counters_getNbTaylorsInTSTests          ( const counters_t st ){ return st->total->nbTaylorsInTSTests        ;}
int counters_getNbTaylorsRepetedInTSTests   ( const counters_t st ){ return st->total->nbTaylorsRepetedInTSTests ;}
int counters_getNbNewton                    ( const counters_t st ){ return st->total->nbNewton                  ;}
int counters_getNbFailingNewton             ( const counters_t st ){ return st->total->nbFailingNewton           ;}
int counters_getNbEval                      ( const counters_t st ){ return st->total->nbEval                    ;}

/* DEPRECATED
void counters_by_depth_get_lenghts_of_str( counters_by_depth_t res, counters_by_depth_t st){
    
    st->nbDiscarded                = 1 + (int) floor(log10( st->nbDiscarded               +1 ) );
    st->nbValidated                = 1 + (int) floor(log10( st->nbValidated               +1 ) );
    st->nbT0Tests                  = 1 + (int) floor(log10( st->nbT0Tests                 +1 ) );
    st->nbFailingT0Tests           = 1 + (int) floor(log10( st->nbFailingT0Tests          +1 ) );
    st->nbGraeffeInT0Tests         = 1 + (int) floor(log10( st->nbGraeffeInT0Tests        +1 ) );
    st->nbGraeffeRepetedInT0Tests  = 1 + (int) floor(log10( st->nbGraeffeRepetedInT0Tests +1 ) );
    st->nbTaylorsRepetedInT0Tests  = 1 + (int) floor(log10( st->nbTaylorsRepetedInT0Tests +1 ) );
    st->nbTSTests                  = 1 + (int) floor(log10( st->nbTSTests                 +1 ) );
    st->nbFailingTSTests           = 1 + (int) floor(log10( st->nbFailingTSTests          +1 ) );
    st->nbGraeffeInTSTests         = 1 + (int) floor(log10( st->nbGraeffeInTSTests        +1 ) );
    st->nbGraeffeRepetedInTSTests  = 1 + (int) floor(log10( st->nbGraeffeRepetedInTSTests +1 ) );
    st->nbTaylorsRepetedInTSTests  = 1 + (int) floor(log10( st->nbTaylorsRepetedInTSTests +1 ) );
    st->nbNewton                   = 1 + (int) floor(log10( st->nbNewton                  +1 ) );
    st->nbFailingNewton            = 1 + (int) floor(log10( st->nbFailingNewton           +1 ) );
    
}

void counters_get_lenghts_of_str( counters_by_depth_t res, counters_by_depth_t st){
    
    counters_by_depth_t temp;
    
    if (st->size == 0) {
        counters_by_depth_init(temp);
        counters_by_depth_get_lenghts_of_str(res, temp);
        return;
    }
    
    counters_by_depth_get_lenghts_of_str(res, (st->table)[0]);
    for (int i = 1; i< size; i++) {
        counters_by_depth_get_lenghts_of_str(temp, (st->table)[i]);
        
        if (temp->nbDiscarded               > res-> nbDiscarded              )  res-> nbDiscarded               = temp->nbDiscarded              ;
        if (temp->nbValidated               > res-> nbValidated              )  res-> nbValidated               = temp->nbValidated              ;
        if (temp->nbT0Tests                 > res-> nbT0Tests                )  res-> nbT0Tests                 = temp->nbT0Tests                ;
        if (temp->nbFailingT0Tests          > res-> nbFailingT0Tests         )  res-> nbFailingT0Tests          = temp->nbFailingT0Tests         ;
        if (temp->nbGraeffeInT0Tests        > res-> nbGraeffeInT0Tests       )  res-> nbGraeffeInT0Tests        = temp->nbGraeffeInT0Tests       ;
        if (temp->nbGraeffeRepetedInT0Tests > res-> nbGraeffeRepetedInT0Tests)  res-> nbGraeffeRepetedInT0Tests = temp->nbGraeffeRepetedInT0Tests;
        if (temp->nbTaylorsRepetedInT0Tests > res-> nbTaylorsRepetedInT0Tests)  res-> nbTaylorsRepetedInT0Tests = temp->nbTaylorsRepetedInT0Tests;
        if (temp->nbTSTests                 > res-> nbTSTests                )  res-> nbTSTests                 = temp->nbTSTests                ;
        if (temp->nbFailingTSTests          > res-> nbFailingTSTests         )  res-> nbFailingTSTests          = temp->nbFailingTSTests         ;
        if (temp->nbGraeffeInTSTests        > res-> nbGraeffeInTSTests       )  res-> nbGraeffeInTSTests        = temp->nbGraeffeInTSTests       ;
        if (temp->nbGraeffeRepetedInTSTests > res-> nbGraeffeRepetedInTSTests)  res-> nbGraeffeRepetedInTSTests = temp->nbGraeffeRepetedInTSTests;
        if (temp->nbTaylorsRepetedInTSTests > res-> nbTaylorsRepetedInTSTests)  res-> nbTaylorsRepetedInTSTests = temp->nbTaylorsRepetedInTSTests;
        if (temp->nbNewton                  > res-> nbNewton                 )  res-> nbNewton                  = temp->nbNewton                 ;
        if (temp->nbFailingNewton           > res-> nbFailingNewton          )  res-> nbFailingNewton           = temp->nbFailingNewton          ;
       
    }
    
}
*/