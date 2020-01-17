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

#include "metadatas/pwSuDatas.h"

void pwSuDatas_init( pwSuDatas_t p ){
    realRat_init(pwSuDatas_wantedPrecref(p));
    realRat_init(pwSuDatas_isolaRatioref(p));
    realRat_set_si(pwSuDatas_wantedPrecref(p), 1, 4);
    p->nbPntsEval = 0;
    p->nbPwSuComp = 0;
    p->evalPoly = NULL;
}

void pwSuDatas_clear( pwSuDatas_t p ){
    realRat_clear(pwSuDatas_wantedPrecref(p));
    realRat_clear(pwSuDatas_isolaRatioref(p));
}
