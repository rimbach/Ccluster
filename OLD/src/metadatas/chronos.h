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

#ifndef CHRONOS_H
#define CHRONOS_H

#include<time.h>

typedef struct {
    clock_t _clicks_Approxi;
    clock_t _clicks_Graeffe;
    clock_t _clicks_Taylors;
    clock_t _clicks_T0Tests;
    clock_t _clicks_TSTests;
    clock_t _clicks_Newtons;
    clock_t _clicks_CclusAl;
    double  _clicks_Approxi_cumul;
    double  _clicks_Graeffe_cumul;
    double  _clicks_Taylors_cumul;
    double  _clicks_T0Tests_cumul;
    double  _clicks_TSTests_cumul;
    double  _clicks_Newtons_cumul;
    double  _clicks_CclusAl_cumul;
} chronos;

typedef chronos chronos_t[1];
typedef chronos * chronos_ptr;
    
void chronos_init( chronos_t times );
void chronos_clear( chronos_t times );

static __inline__ void   chronos_tic_Approxi      ( chronos_t times ) { times->_clicks_Approxi = clock(); }
static __inline__ void   chronos_toc_Approxi      ( chronos_t times ) { times->_clicks_Approxi_cumul += ((double) (clock() - times->_clicks_Approxi))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_Approxi ( const chronos_t times ) { return times->_clicks_Approxi_cumul; }

static __inline__ void   chronos_tic_Graeffe      ( chronos_t times ) { times->_clicks_Graeffe = clock(); }
static __inline__ void   chronos_toc_Graeffe      ( chronos_t times ) { times->_clicks_Graeffe_cumul += ((double) (clock() - times->_clicks_Graeffe))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_Graeffe ( const chronos_t times ) { return times->_clicks_Graeffe_cumul; }

static __inline__ void   chronos_tic_Taylors      ( chronos_t times ) { times->_clicks_Taylors = clock(); }
static __inline__ void   chronos_toc_Taylors      ( chronos_t times ) { times->_clicks_Taylors_cumul += ((double) (clock() - times->_clicks_Taylors))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_Taylors ( const chronos_t times ) { return times->_clicks_Taylors_cumul; }

static __inline__ void   chronos_tic_T0Tests      ( chronos_t times ) { times->_clicks_T0Tests = clock(); }
static __inline__ void   chronos_toc_T0Tests      ( chronos_t times ) { times->_clicks_T0Tests_cumul += ((double) (clock() - times->_clicks_T0Tests))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_T0Tests ( const chronos_t times ) { return times->_clicks_T0Tests_cumul; }

static __inline__ void   chronos_tic_TSTests      ( chronos_t times ) { times->_clicks_TSTests = clock(); }
static __inline__ void   chronos_toc_TSTests      ( chronos_t times ) { times->_clicks_TSTests_cumul += ((double) (clock() - times->_clicks_TSTests))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_TSTests ( const chronos_t times ) { return times->_clicks_TSTests_cumul; }

static __inline__ void   chronos_tic_Newtons      ( chronos_t times ) { times->_clicks_Newtons = clock(); }
static __inline__ void   chronos_toc_Newtons      ( chronos_t times ) { times->_clicks_Newtons_cumul += ((double) (clock() - times->_clicks_Newtons))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_Newtons ( const chronos_t times ) { return times->_clicks_Newtons_cumul; }

static __inline__ void   chronos_tic_CclusAl      ( chronos_t times ) { times->_clicks_CclusAl = clock(); }
static __inline__ void   chronos_toc_CclusAl      ( chronos_t times ) { times->_clicks_CclusAl_cumul += ((double) (clock() - times->_clicks_CclusAl))/ CLOCKS_PER_SEC; }
static __inline__ double chronos_get_time_CclusAl ( const chronos_t times ) { return times->_clicks_CclusAl_cumul; }

static __inline__ void   chronos_tic_Test      ( chronos_t times, int discard ) { 
    if (discard) times->_clicks_T0Tests = clock(); 
    else         times->_clicks_TSTests = clock();
}

static __inline__ void   chronos_toc_Test      ( chronos_t times, int discard ) { 
    if (discard) times->_clicks_T0Tests_cumul += ((double) (clock() - times->_clicks_T0Tests))/ CLOCKS_PER_SEC;
    else         times->_clicks_TSTests_cumul += ((double) (clock() - times->_clicks_TSTests))/ CLOCKS_PER_SEC;
}

double chronos_get_time_Approxi_for_julia ( const chronos_t times ) ;
double chronos_get_time_Graeffe_for_julia ( const chronos_t times ) ;
double chronos_get_time_Taylors_for_julia ( const chronos_t times ) ;
double chronos_get_time_T0Tests_for_julia ( const chronos_t times ) ;
double chronos_get_time_TSTests_for_julia ( const chronos_t times ) ;
double chronos_get_time_Newtons_for_julia ( const chronos_t times ) ;
double chronos_get_time_CclusAl_for_julia ( const chronos_t times ) ;

void chronos_add_time_Approxi_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_Graeffe_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_Taylors_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_T0Tests_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_TSTests_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_Newtons_for_julia ( chronos_t times, double ellapsedTime ) ;
void chronos_add_time_CclusAl_for_julia ( chronos_t times, double ellapsedTime ) ;

#endif