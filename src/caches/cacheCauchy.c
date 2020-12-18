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

#include <stdlib.h>
#include "cacheCauchy.h"

void cacheCauchy_init ( cacheCauchy_t cache, 
                        void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                        slong degree,
                        slong num,
                        ulong den
                      ){
    cacheCauchy_evalFastref(cache) = evalFast;
    cacheCauchy_degreeref(cache) = degree;
    
    realRat_init(cacheCauchy_isoRatioref(cache));
    realRat_set_si(cacheCauchy_isoRatioref(cache), num, den);
    
    realApp_init(cacheCauchy_wanErrCoref(cache));
    realApp_init(cacheCauchy_wanErrCeref(cache));
    
    realRat_init(cacheCauchy_lBoundUnref(cache));
    realRat_init(cacheCauchy_uBoundUnref(cache));
    
    /* compute nbEvalCo = log_isoRatio (2*degree +1) */
    realApp_t liR;
    realApp_init(liR);
    realApp_set_realRat(liR, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_log(liR, liR, CCLUSTER_DEFAULT_PREC );
    realApp_t q1App;
    realApp_init(q1App);
    realApp_set_si(q1App, 4*degree+1);
    realApp_log(q1App, q1App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q1App, q1App, liR, CCLUSTER_DEFAULT_PREC );
    slong q1 = realApp_ceil_si( q1App, CCLUSTER_DEFAULT_PREC ) +1;
    cacheCauchy_nbEvalCoref(cache) = q1;
    
    /* compute nbEvalCe = max ( log_isoRatio (2*degree*nbEvalCo +1), degree +1 ) s.t. nbEvalCe multiple of nbEvalCo */
    realApp_t q2App;
    realApp_init(q2App);
    realApp_set_si(q2App, 4*degree*q1+1);
    realApp_log(q2App, q2App, CCLUSTER_DEFAULT_PREC );
    realApp_div(q2App, q2App, liR, CCLUSTER_DEFAULT_PREC );
    slong q2 = realApp_ceil_si( q2App, CCLUSTER_DEFAULT_PREC ) +1;
    q2 = CCLUSTER_MAX(q2, degree +1);
    slong quo = ((slong) q2/q1) +1;
    q2 = q1*quo;
    cacheCauchy_quotientref(cache) = quo;
    cacheCauchy_nbEvalCeref(cache) = q2;
    
    /* compute error cacheCauchy_wanErrCoref(cache) = (d*isoRatio^(-q1))/(1-isoRatio^(-q1)) */
    realApp_ptr wP = cacheCauchy_wanErrCoref(cache);
    realApp_t tempApp;
    realApp_init(tempApp);
    realApp_set_realRat(wP, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP, wP, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP, wP, q1, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP, wP, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP, wP, degree, CCLUSTER_DEFAULT_PREC);
    
    /* compute error cacheCauchy_wanErrCeref(cache) = (d*isoRatio^(-q2))/(1-isoRatio^(-q2)) */
    realApp_ptr wP2 = cacheCauchy_wanErrCeref(cache);
    realApp_set_realRat(wP2, cacheCauchy_isoRatioref(cache), CCLUSTER_DEFAULT_PREC);
    realApp_inv(wP2, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_pow_ui(wP2, wP2, q2, CCLUSTER_DEFAULT_PREC);
    realApp_set_si(tempApp, 1);
    realApp_sub(tempApp, tempApp, wP2, CCLUSTER_DEFAULT_PREC);
    realApp_div(wP2, wP2, tempApp, CCLUSTER_DEFAULT_PREC);
    realApp_mul_si(wP2, wP2, degree, CCLUSTER_DEFAULT_PREC);
    
    /* compute lower bound unit: (isoRatio-1)^d/isoRatio^d */
    realRat_ptr lb = cacheCauchy_lBoundUnref(cache);
    realRat_add_si(lb, cacheCauchy_isoRatioref(cache), -1);
    realRat_div(lb, lb, cacheCauchy_isoRatioref(cache));
    realRat_pow_si(lb, lb, degree);
    
    /* compute upper bound: (d*(isoRatio+1)/(isoRatio-1) */
    realRat_ptr ub = cacheCauchy_uBoundUnref(cache);
    realRat_t temp;
    realRat_init(temp);
    realRat_add_si(ub, cacheCauchy_isoRatioref(cache), +1);
    realRat_mul_si(ub, ub, degree);
    realRat_add_si(temp, cacheCauchy_isoRatioref(cache), -1);
    realRat_div(ub, ub, temp);
    
    realApp_clear(liR);
    realApp_clear(q1App);
    realApp_clear(q2App);
    realApp_clear(tempApp);
    
    realRat_clear(temp);
}

void cacheCauchy_clear ( cacheCauchy_t cache ){
    
    realRat_clear(cacheCauchy_isoRatioref(cache));
    
    realApp_clear(cacheCauchy_wanErrCoref(cache));
    realApp_clear(cacheCauchy_wanErrCeref(cache));
    
    realRat_clear(cacheCauchy_lBoundUnref(cache));
    realRat_clear(cacheCauchy_uBoundUnref(cache));
}
