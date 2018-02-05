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

#ifndef COUNTERS_H
#define COUNTERS_H

#include <stdlib.h>
// #include <stdio.h>
// #include <math.h>

typedef struct {
    int nbDiscarded;
    int nbValidated;
    int nbSolutions;
    /* T0Tests */
    int nbT0Tests;
    int nbFailingT0Tests;
    int nbGraeffeInT0Tests;
    int nbGraeffeRepetedInT0Tests;
    int nbTaylorsRepetedInT0Tests;
    /* TSTests */
    int nbTSTests;
    int nbFailingTSTests;
    int nbGraeffeInTSTests;
    int nbGraeffeRepetedInTSTests;
    int nbTaylorsRepetedInTSTests;
    /* Newton steps */
    int nbNewton;
    int nbFailingNewton;
} counters_by_depth;

typedef counters_by_depth counters_by_depth_t[1];
typedef counters_by_depth * counters_by_depth_ptr;

void counters_by_depth_init( counters_by_depth_t st);
static __inline__ void counters_by_depth_clear( counters_by_depth_t st) {}

// void counters_by_depth_get_lenghts_of_str( counters_by_depth_t res, counters_by_depth_t st);

typedef struct {
        int size;
        int size_allocated;
        counters_by_depth_ptr table;
        counters_by_depth_t total;
} counters; /* a stat is a table of counters by depth */

typedef counters counters_t[1];

#define INIT_SIZE_STATS 1000

void counters_init( counters_t st);
void counters_clear( counters_t st);

void counters_add_discarded( counters_t st, int depth );
void counters_add_validated( counters_t st, int depth, int nbSols );
void counters_add_Test     ( counters_t st, int depth, int res, int discard, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
// void counters_add_discarding_test( counters_t st, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
// void counters_add_validating_test( counters_t st, int depth, int res, int nbTaylorsRepeted, int nbGraeffe, int nbGraeffeRepeted);
void counters_add_Newton   ( counters_t st, int depth, int res );

// void counters_get_lenghts_of_str( counters_by_depth_t res, counters_by_depth_t st);

void counters_count ( counters_t st );
int counters_getDepth( const counters_t st);

int counters_getNbDiscarded                 ( const counters_t st );
int counters_getNbValidated                 ( const counters_t st );
int counters_getNbSolutions                 ( const counters_t st );
int counters_getNbT0Tests                   ( const counters_t st );
int counters_getNbFailingT0Tests            ( const counters_t st );
int counters_getNbGraeffeInT0Tests          ( const counters_t st );
int counters_getNbGraeffeRepetedInT0Tests   ( const counters_t st );
int counters_getNbTaylorsRepetedInT0Tests   ( const counters_t st );
int counters_getNbTSTests                   ( const counters_t st );
int counters_getNbFailingTSTests            ( const counters_t st );
int counters_getNbGraeffeInTSTests          ( const counters_t st );
int counters_getNbGraeffeRepetedInTSTests   ( const counters_t st );
int counters_getNbTaylorsRepetedInTSTests   ( const counters_t st );
int counters_getNbNewton                    ( const counters_t st );
int counters_getNbFailingNewton             ( const counters_t st );
#endif