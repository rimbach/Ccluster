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

#include "compRat.h"

// void compRat_init_forJulia(compRat_t x){
//     compRat_init(x);
// }
// void compRat_clear_forJulia(compRat_t x){
//     compRat_clear(x);
// }

void compRat_fprint(FILE * file, const compRat_t x){
    
    if ( realRat_is_zero(compRat_imagref(x)) ) {
        realRat_fprint( file, compRat_realref(x) );
    }
    else if ( realRat_is_zero(compRat_realref(x)) ) {
        fputc('j', file);
        realRat_fprint( file, compRat_imagref(x) );
    }
    else {
        fputc('(', file);
        realRat_fprint( file, compRat_realref(x) );
        fputc('+', file);
        fputc('j', file);
        realRat_fprint( file, compRat_imagref(x) );
        fputc(')', file);
    }
}

/* sets dest to abs(x.real-y.real) + i*abs(x.imag-y.imag) */
void compRat_comp_distance( compRat_t dest, const compRat_t x, const compRat_t y ){
    realRat_sub( compRat_realref(dest), compRat_realref(x), compRat_realref(y) );
    realRat_sub( compRat_imagref(dest), compRat_imagref(x), compRat_imagref(y) );
    realRat_abs( compRat_realref(dest), compRat_realref(dest));
    realRat_abs( compRat_imagref(dest), compRat_imagref(dest));
}

/* returns negative if x.real<y.real or (x.real=y.real and x.imag<y.imag) */
/*                0 if x.real=y.real and x.imag=y.imag                    */
/*         positive if x.real>y.real or (x.real=y.real and x.imag>y.imag) */
int  compRat_cmp    (const compRat_t x, const compRat_t y){
    int res = realRat_cmp( compRat_realref(x), compRat_realref(y) );
    if (res==0)
        return realRat_cmp( compRat_imagref(x), compRat_imagref(y) );
    return res;
}