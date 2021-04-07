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

#include "metadatas/cauchyDatas.h"

void cauchyDatas_init( cauchyDatas_t p ){
    
    realRat_init(cauchyDatas_isoRatioref(p));
    
}

void cauchyDatas_clear( cauchyDatas_t p ){
    realRat_clear(cauchyDatas_isoRatioref(p));
}
