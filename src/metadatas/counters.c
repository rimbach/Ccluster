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

void boxes_by_prec_init( boxes_by_prec_t bt ){
    bt->size = 0;
    bt->size_allocated = INIT_SIZE_STATS;
    bt->table = (int *) ccluster_malloc (INIT_SIZE_STATS*sizeof(int));
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(bt->_mutex), NULL);
#endif 
}

void boxes_by_prec_clear( boxes_by_prec_t bt ){
    ccluster_free(bt->table);
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_destroy( &(bt->_mutex) );
#endif
}

void boxes_by_prec_adjust_table( boxes_by_prec_t bt, int index ){
    
    /* realocate size if needed */
    while (index+1 > bt->size_allocated) {
        bt->size_allocated += INIT_SIZE_STATS;
        bt->table = (int *) ccluster_realloc ( (void *) bt->table, (bt->size_allocated)*sizeof(int) );
    }
    /* set to 0 intermediate cases if needed */
    while (index+1 > bt->size ){
        bt->table[bt->size] = 0;
        bt->size += 1;
    }
}

void boxes_by_prec_add_int( boxes_by_prec_t bt, slong prec, int nbBoxes ){
    int indexpow2 = ((int) prec/CCLUSTER_DEFAULT_PREC);
    int index=0;
    while (0x1<<index < indexpow2) index+=1;
    boxes_by_prec_adjust_table(bt, index);
    bt->table[index] +=nbBoxes;
}

void boxes_by_prec_add_boxes_by_prec( boxes_by_prec_t bt, boxes_by_prec_t t ){
    
    for (int index = 0; index < t->size; index ++ ){
        boxes_by_prec_adjust_table(bt, index);
        bt->table[index] += t->table[index];
    }
}

int  boxes_by_prec_fprint( FILE * file, const boxes_by_prec_t bt ){
    int r = 0;
    for (int index = 0; index < bt->size; index ++ ){
        char buffer[50];
        r = sprintf (buffer, "boxes with %d:", ( 0x1<<(index) )*CCLUSTER_DEFAULT_PREC);
        r = fprintf(file, "|%-39s %14d %14s|\n", buffer,           bt->table[index],    " " );
    }
    return r;
}



void counters_by_depth_init( counters_by_depth_t st) {
    st->nbDiscarded               = 0;
    st->nbValidated               = 0;
    st->nbSolutions               = 0;
    st->nbExplored                = 0;
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
    st->nbEval                    = 0;
    st->nbPsCountingTest          = 0;

    boxes_by_prec_init( st->bpc );
}

// void counters_by_depth_join( counters_by_depth_t c1, const counters_by_depth_t c2){
//     c1->nbDiscarded                += c2->nbDiscarded               ;
//     c1->nbValidated                += c2->nbValidated               ;
//     c1->nbSolutions                += c2->nbSolutions               ;
//     c1->nbT0Tests                  += c2->nbT0Tests                 ;
//     c1->nbFailingT0Tests           += c2->nbFailingT0Tests          ;
//     c1->nbGraeffeInT0Tests         += c2->nbGraeffeInT0Tests        ;
//     c1->nbGraeffeRepetedInT0Tests  += c2->nbGraeffeRepetedInT0Tests ;
//     c1->nbTaylorsInT0Tests         += c2->nbTaylorsInT0Tests        ;
//     c1->nbTaylorsRepetedInT0Tests  += c2->nbTaylorsRepetedInT0Tests ;
//     c1->nbTSTests                  += c2->nbTSTests                 ;
//     c1->nbFailingTSTests           += c2->nbFailingTSTests          ;
//     c1->nbGraeffeInTSTests         += c2->nbGraeffeInTSTests        ;
//     c1->nbGraeffeRepetedInTSTests  += c2->nbGraeffeRepetedInTSTests ;
//     c1->nbTaylorsInTSTests         += c2->nbTaylorsInTSTests        ;
//     c1->nbTaylorsRepetedInTSTests  += c2->nbTaylorsRepetedInTSTests ;
//     c1->nbNewton                   += c2->nbNewton                  ;
//     c1->nbFailingNewton            += c2->nbFailingNewton           ;
//     c1->nbEval                     += c2->nbEval                    ;
// }

void counters_init( counters * st) {
    st->size = 0;
    st->size_allocated = INIT_SIZE_STATS;
    st->table = (counters_by_depth_ptr) ccluster_malloc (INIT_SIZE_STATS*sizeof(counters_by_depth));
    counters_by_depth_init( st->total );
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(st->_mutex), NULL);
#endif    
}

void counters_clear( counters * st) {
    ccluster_free(st->table);
    counters_by_depth_clear( st->total );
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_destroy( &(st->_mutex) );
#endif
}

// void counters_join_depth( counters_t c1, const counters_by_depth_t c2, int depth){
//     counters_adjust_table(c1, depth);
//     counters_by_depth_join( &(c1->table[depth]), c2);
// }
// 
// void counters_join( counters_t c1, const counters_t c2){
//     for (int i=0; i<c2->size; i++)
//         counters_join_depth(c1, &(c2->table[i]), i);
// }

void counters_adjust_table( counters_t st, int depth ){
    /* realocate size if needed */
    while (depth+1 > st->size_allocated) {
        st->size_allocated += INIT_SIZE_STATS;
        st->table = (counters_by_depth_ptr) ccluster_realloc ( (void *) st->table, (st->size_allocated)*sizeof(counters_by_depth) );
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

void counters_add_explored ( counters_t st, int depth ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbExplored +=1;
}

void counters_add_Test     ( counters_t st, int depth, int res, int discard, 
                             int nbTaylors, int nbTaylorsRepeted, 
                             int nbGraeffe, int nbGraeffeRepeted,
                             slong prec
                           ){
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
    boxes_by_prec_add_int( (st->table[depth]).bpc, prec, 1);    
}

void counters_add_Newton   ( counters_t st, int depth, int res ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbNewton                  +=1;
    if (!res) (st->table[depth]).nbFailingNewton +=1;
}

void counters_add_Eval( counters_t st, int nbEvals, int depth ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbEval                  +=nbEvals;
}

void counters_add_PsCountingTest( counters_t st, int depth ){
    counters_adjust_table(st, depth);
    (st->table[depth]).nbPsCountingTest                  +=1;
}

void counters_count ( counters_t st ) {
    for (int i = 0; i< st->size; i++) {
       st->total->nbDiscarded               += (st->table)[i].nbDiscarded               ; 
       st->total->nbValidated               += (st->table)[i].nbValidated               ; 
       st->total->nbSolutions               += (st->table)[i].nbSolutions               ;
       st->total->nbExplored                += (st->table)[i].nbExplored                ;
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
       st->total->nbPsCountingTest          += (st->table)[i].nbPsCountingTest           ;
       boxes_by_prec_add_boxes_by_prec( st->total->bpc, (st->table)[i].bpc ); 
    }

}

int counters_getDepth( const counters_t st) { return st->size;}

int counters_getNbDiscarded                 ( const counters_t st ){ return st->total->nbDiscarded               ;}
int counters_getNbValidated                 ( const counters_t st ){ return st->total->nbValidated               ;}
int counters_getNbSolutions                 ( const counters_t st ){ return st->total->nbSolutions               ;}
int counters_getNbExplored                  ( const counters_t st ){ return st->total->nbExplored                ;}
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
int counters_getNbPsCountingTest            ( const counters_t st ){ return st->total->nbPsCountingTest                    ;}

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
