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
    times->_clicks_Derivat_cumul = 0.0;
    times->_clicks_Anticip_cumul = 0.0;
    
    times->_clicks_PSTests_cumul = 0.0;
#ifdef CCLUSTER_STATS_PS
    times->_clicks_PSTestV_cumul = 0.0;
    times->_clicks_Evaluat_cumul = 0.0;
#endif
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_init ( &(times->_mutex), NULL);
#endif
}

void chronos_clear( chronos_t times ) {
#ifdef CCLUSTER_HAVE_PTHREAD
    pthread_mutex_destroy( &(times->_mutex) );
#endif    
}

// void chronos_join ( chronos_t t1, const chronos_t t2 ){
//     t1->_clicks_Approxi_cumul += t2->_clicks_Approxi_cumul;
//     t1->_clicks_Graeffe_cumul += t2->_clicks_Graeffe_cumul;
//     t1->_clicks_Taylors_cumul += t2->_clicks_Taylors_cumul;
//     t1->_clicks_T0Tests_cumul += t2->_clicks_T0Tests_cumul;
//     t1->_clicks_TSTests_cumul += t2->_clicks_TSTests_cumul;
//     t1->_clicks_Newtons_cumul += t2->_clicks_Newtons_cumul;
//     t1->_clicks_CclusAl_cumul += t2->_clicks_CclusAl_cumul;
//     t1->_clicks_Evaluat_cumul += t2->_clicks_Evaluat_cumul;
//     t1->_clicks_Derivat_cumul += t2->_clicks_Derivat_cumul;
//     t1->_clicks_Anticip_cumul += t2->_clicks_Anticip_cumul;
// }
