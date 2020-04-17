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

#include <stdio.h>
#include "geometry/compAnn.h"

void compAnn_fprintd( FILE * file, const compAnn_t x, slong digits ){
    fprintf(file, "indMax: %ld, indMin: %ld, rrInPo: %d, rrInNe: %d \n", 
            compAnn_indMaxref(x), compAnn_indMinref(x), compAnn_rrInPoref(x), compAnn_rrInNeref(x) );
    fprintf(file, "radInf: ");
    realApp_fprintd(file, compAnn_radInfref(x), digits);
    fprintf(file, "  radSup: ");
    realApp_fprintd(file, compAnn_radSupref(x), digits);
    fprintf(file, "\n");
}

void compAnn_fprint( FILE * file, const compAnn_t x){
    fprintf(file, "indMax: %ld, indMin: %ld, rrInPo: %d, rrInNe: %d \n", 
            compAnn_indMaxref(x), compAnn_indMinref(x), compAnn_rrInPoref(x), compAnn_rrInNeref(x) );
    fprintf(file, "radInf: ");
    realApp_fprint(file, compAnn_radInfref(x));
    fprintf(file, "  radSup: ");
    realApp_fprint(file, compAnn_radSupref(x));
    fprintf(file, "\n");
}
