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

#ifndef DEFLATE_H
#define DEFLATE_H

#ifdef DEFLATE_INLINE_C
#define DEFLATE_INLINE
#else
#define DEFLATE_INLINE static __inline__
#endif

#include "geometry/connCmp.h"

#ifdef __cplusplus
extern "C" {
#endif
    
#define connCmp_isDefref(X) ( (X)->isDef)
#define connCmp_degDeref(X) ( (X)->degDe)
#define connCmp_defPoref(X) (&(X)->defPo)
#define connCmp_sumAbref(X) (&(X)->sumAb)


#ifdef __cplusplus
}
#endif

#endif
