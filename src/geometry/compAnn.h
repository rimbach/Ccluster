/* ************************************************************************** */
/*  Copyright (C) 2020 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#ifndef COMPANN_H
#define COMPANN_H

#ifdef GEOMETRY_INLINE_C
#define GEOMETRY_INLINE
#else
#define GEOMETRY_INLINE static __inline__
#endif

#include "numbers/realApp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    slong   indMax;
    slong   indMin;
    realApp radInf;
    realApp radSup;
    int     rrInPo; /* 1 if it contains a unique real root in its intersection with R+ */
    int     rrInNe; /* 1 if it contains a unique real root in its intersection with R- */
} compAnn;

typedef compAnn compAnn_t[1];
typedef compAnn * compAnn_ptr;

#define compAnn_indMaxref(X) (X->indMax)
#define compAnn_indMinref(X) (X->indMin)
#define compAnn_radInfref(X) (&(X)->radInf)
#define compAnn_radSupref(X) (&(X)->radSup)
#define compAnn_rrInPoref(X) (X->rrInPo)
#define compAnn_rrInNeref(X) (X->rrInNe)

GEOMETRY_INLINE void compAnn_init( compAnn_t x ){
    realApp_init( compAnn_radInfref(x) );
    realApp_init( compAnn_radSupref(x) );
    compAnn_rrInPoref(x) = -1;
    compAnn_rrInNeref(x) = -1;
}

GEOMETRY_INLINE void compAnn_clear( compAnn_t x ){
    realApp_clear( compAnn_radInfref(x) );
    realApp_clear( compAnn_radSupref(x) );
}

GEOMETRY_INLINE void compAnn_set( compAnn_t x, slong indMax, slong indMin, const realApp_t radInf, const realApp_t radSup){
    compAnn_indMaxref(x) = indMax;
    compAnn_indMinref(x) = indMin;
    realApp_set( compAnn_radInfref(x), radInf);
    realApp_set( compAnn_radSupref(x), radSup);
}

void compAnn_fprintd( FILE * file, const compAnn_t x, slong digits );
GEOMETRY_INLINE void compAnn_printd( const compAnn_t x, slong digits ){ compAnn_fprintd( stdout, x, digits ); }

void compAnn_fprint( FILE * file, const compAnn_t x);
GEOMETRY_INLINE void compAnn_print( const compAnn_t x ){ compAnn_fprint( stdout, x); }

#ifdef __cplusplus
}
#endif

#endif    
