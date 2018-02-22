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

#include "metadatas/chronos.h"

void chronos_init( chronos_t times ){
    times->_clicks_Approxi_cumul = 0.0;
    times->_clicks_Graeffe_cumul = 0.0;
    times->_clicks_Taylors_cumul = 0.0;
    times->_clicks_T0Tests_cumul = 0.0;
    times->_clicks_TSTests_cumul = 0.0;
    times->_clicks_Newtons_cumul = 0.0;
    times->_clicks_CclusAl_cumul = 0.0;
    times->_clicks_Evaluat_cumul = 0.0;
    times->_clicks_Derivat_cumul = 0.0;
}

void chronos_clear( chronos_t times ) {}