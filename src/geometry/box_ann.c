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
#include "geometry/compBox.h"
#include "geometry/compAnn.h"

void compBox_init_annulii(compBox_t x) { 
    for (int i=0; i<GEOMETRY_NB_ANN_PER_BOX; i++)
        compAnn_list_init(compBox_annuliref(x,i));
}

void compBox_copy_annulii(compBox_t x, int ind, const compAnn_list_t anulii){
     compAnn_list_copy( compBox_annuliref(x,ind), anulii);
}

void compBox_clear_annulii(compBox_t x){
    for (int i=0; i<GEOMETRY_NB_ANN_PER_BOX; i++){
        compAnn_list_empty(compBox_annuliref(x,i));
        compAnn_list_clear(compBox_annuliref(x,i));
    }
}

/* assume compBox_annuliref(x) contains anulii in increasing order of their radii */
void compBox_actualize_anulii_compAnn_list_risolate( compBox_t x, int ind, const compAnn_list_t lmother ){
    
    realApp_t center;
    realApp_t left, right, rad;
    realApp_init(center);
    realApp_init(left);
    realApp_init(right);
    realApp_init(rad);
    
    compAnn_list_iterator it = compAnn_list_begin( lmother );
    
    realApp_set_realRat( center, compRat_realref(compBox_centerref(x)), CCLUSTER_DEFAULT_PREC );
    realApp_add_si(center, center, -compAnn_getCenterRe(compAnn_list_elmt( it )), CCLUSTER_DEFAULT_PREC );
    realApp_abs        ( left,   center );
    realApp_set        ( right,  left );
    realApp_set_realRat( rad,    compBox_bwidthref(x), CCLUSTER_DEFAULT_PREC );
    realApp_div_si     ( rad,    rad,                  2, CCLUSTER_DEFAULT_PREC );
    realApp_sub        ( left,   left,                 rad, CCLUSTER_DEFAULT_PREC );
    realApp_add        ( right,  right,                rad, CCLUSTER_DEFAULT_PREC );
    
    /* go to the first annulus with rsup >= left */
    while ( (it!=compAnn_list_end()) && (realApp_lt( compAnn_radSupref(compAnn_list_elmt( it )), left )==1)  )
        it = compAnn_list_next(it);
    
    /* push annulii satisfying rinf <= right (and rsup >= left) */
    while ( (it!=compAnn_list_end()) && (!(realApp_gt( compAnn_radInfref(compAnn_list_elmt( it )), right )==1))  ) {
            
        compAnn_list_push( compBox_annuliref(x,ind), compAnn_list_elmt( it ) );
        
        it = compAnn_list_next(it);
    }
    
    realApp_clear(center);
    realApp_clear(left);
    realApp_clear(right);
    realApp_clear(rad);
}

void compBox_actualize_anulii_risolate( compBox_t x, const compBox_t bmother ){
        for (int ind=0; ind<GEOMETRY_NB_ANN_PER_BOX; ind++)
        if (compAnn_list_get_size(compBox_annuliref(bmother,ind))>=1)
            compBox_actualize_anulii_compAnn_list_risolate( x, ind, compBox_annuliref(bmother, ind) );
}

/* assume lmother contains at least one element */
void compBox_actualize_anulii_compAnn_list_real( compBox_t x, int ind, const compAnn_list_t lmother ){
    
//     if (compAnn_list_get_size(lmother) >= 1) {
        compRat_t shiftedCenter;
        compApp_t temp;
        realApp_t left, right;
        compRat_t closest, furthest;
        realRat_t halfwidth;
        compApp_init(temp);
        realApp_init(left);
        realApp_init(right);
        compRat_init(shiftedCenter);
        compRat_init(closest);
        compRat_init(furthest);
        realRat_init(halfwidth);
        
        
        compAnn_list_iterator it;
        it = compAnn_list_begin( lmother );
        
        compRat_set(shiftedCenter, compBox_centerref(x));
        realRat_add_si(compRat_realref(shiftedCenter), compRat_realref(shiftedCenter), -compAnn_getCenterRe(compAnn_list_elmt( it )));
        realRat_set(halfwidth, compBox_bwidthref(x));
        realRat_div_ui(halfwidth, halfwidth, 2 );
        
        realRat_sub( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
        realRat_add( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_realref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_realref(closest) )<0) {
            if (realRat_sgn( compRat_realref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
                realRat_sub( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_realref(closest), compRat_realref(closest) );
                if (realRat_cmp( compRat_realref(closest), compRat_realref(furthest) ) >0)
                    realRat_neg( compRat_realref(furthest), compRat_realref(closest) );
                realRat_set_si( compRat_realref(closest), 0, 1);
            }
        }
        
        realRat_sub( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
        realRat_add( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_imagref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_imagref(closest) )<0) {
            if (realRat_sgn( compRat_imagref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
                realRat_sub( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_imagref(closest), compRat_imagref(closest) );
                if (realRat_cmp( compRat_imagref(closest), compRat_imagref(furthest) ) >0)
                    realRat_neg( compRat_imagref(furthest), compRat_imagref(closest) );
                realRat_set_si( compRat_imagref(closest), 0, 1);
            }
        }
        
        compApp_set_compRat(temp, closest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (left, temp,    CCLUSTER_DEFAULT_PREC );
        compApp_set_compRat(temp, furthest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (right, temp,    CCLUSTER_DEFAULT_PREC );
        
        /* go to the first annulus with rsup >= left */
        while ( (it!=compAnn_list_end()) && (realApp_lt( compAnn_radSupref(compAnn_list_elmt( it )), left )==1)  ){
            it = compAnn_list_next(it);
        }
        
        /* push annulii satisfying rinf <= right (and rsup >= left) */
        while ( (it!=compAnn_list_end()) 
            && (!(realApp_gt( compAnn_radInfref(compAnn_list_elmt( it )), right )==1))
            ) {
            compAnn_list_push( compBox_annuliref(x,ind), compAnn_list_elmt( it ) );
            it = compAnn_list_next(it);
        }
        
        
        compApp_clear(temp);
        realApp_clear(left);
        realApp_clear(right);
        compRat_clear(shiftedCenter);
        compRat_clear(closest);
        compRat_clear(furthest);
        realRat_clear(halfwidth);
        
//         printf("#actualize anulii: length of list: %d\n", compAnn_list_get_size(compBox_annuliref(x)) );
//     }
}

void compBox_actualize_anulii_compAnn_list_imag( compBox_t x, int ind, const compAnn_list_t lmother ){
    
//     if (compAnn_list_get_size(lmother) >= 1) {
        compRat_t shiftedCenter;
        compApp_t temp;
        realApp_t left, right;
        compRat_t closest, furthest;
        realRat_t halfwidth;
        compApp_init(temp);
        realApp_init(left);
        realApp_init(right);
        compRat_init(shiftedCenter);
        compRat_init(closest);
        compRat_init(furthest);
        realRat_init(halfwidth);
        
        
        compAnn_list_iterator it;
        it = compAnn_list_begin( lmother );
        
        compRat_set(shiftedCenter, compBox_centerref(x));
        realRat_add_si(compRat_imagref(shiftedCenter), compRat_imagref(shiftedCenter), -compAnn_getCenterIm(compAnn_list_elmt( it )));
        realRat_set(halfwidth, compBox_bwidthref(x));
        realRat_div_ui(halfwidth, halfwidth, 2 );
        
        realRat_sub( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
        realRat_add( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_realref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_imagref(closest) )<0) {
            if (realRat_sgn( compRat_imagref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
                realRat_sub( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_imagref(closest), compRat_imagref(closest) );
                if (realRat_cmp( compRat_imagref(closest), compRat_imagref(furthest) ) >0)
                    realRat_neg( compRat_imagref(furthest), compRat_imagref(closest) );
                realRat_set_si( compRat_imagref(closest), 0, 1);
            }
        }
        
        realRat_sub( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
        realRat_add( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_imagref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_realref(closest) )<0) {
            if (realRat_sgn( compRat_realref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
                realRat_sub( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_realref(closest), compRat_realref(closest) );
                if (realRat_cmp( compRat_realref(closest), compRat_realref(furthest) ) >0)
                    realRat_neg( compRat_realref(furthest), compRat_realref(closest) );
                realRat_set_si( compRat_realref(closest), 0, 1);
            }
        }
        
        compApp_set_compRat(temp, closest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (left, temp,    CCLUSTER_DEFAULT_PREC );
        compApp_set_compRat(temp, furthest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (right, temp,    CCLUSTER_DEFAULT_PREC );
        
        /* go to the first annulus with rsup >= left */
        while ( (it!=compAnn_list_end()) && (realApp_lt( compAnn_radSupref(compAnn_list_elmt( it )), left )==1)  ){
            it = compAnn_list_next(it);
        }
        
        /* push annulii satisfying rinf <= right (and rsup >= left) */
        while ( (it!=compAnn_list_end()) 
            && (!(realApp_gt( compAnn_radInfref(compAnn_list_elmt( it )), right )==1))
            ) {
            compAnn_list_push( compBox_annuliref(x,ind), compAnn_list_elmt( it ) );
            it = compAnn_list_next(it);
        }
        
        
        compApp_clear(temp);
        realApp_clear(left);
        realApp_clear(right);
        compRat_clear(shiftedCenter);
        compRat_clear(closest);
        compRat_clear(furthest);
        realRat_clear(halfwidth);
        
//         printf("#actualize anulii: length of list: %d\n", compAnn_list_get_size(compBox_annuliref(x)) );
//     }
}

void compBox_actualize_anulii_compAnn_list_imag_conj( compBox_t x, int ind, const compAnn_list_t lmother ){
    
//     if (compAnn_list_get_size(lmother) >= 1) {
        compRat_t shiftedCenter;
        compApp_t temp;
        realApp_t left, right;
        compRat_t closest, furthest;
        realRat_t halfwidth;
        compApp_init(temp);
        realApp_init(left);
        realApp_init(right);
        compRat_init(shiftedCenter);
        compRat_init(closest);
        compRat_init(furthest);
        realRat_init(halfwidth);
        
        
        compAnn_list_iterator it;
        it = compAnn_list_begin( lmother );
        
        compRat_set(shiftedCenter, compBox_centerref(x));
        realRat_neg( compRat_imagref(shiftedCenter), compRat_imagref(shiftedCenter) );
        realRat_add_si(compRat_imagref(shiftedCenter), compRat_imagref(shiftedCenter), -compAnn_getCenterIm(compAnn_list_elmt( it )));
        realRat_set(halfwidth, compBox_bwidthref(x));
        realRat_div_ui(halfwidth, halfwidth, 2 );
        
        realRat_sub( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
        realRat_add( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_realref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_imagref(closest) )<0) {
            if (realRat_sgn( compRat_imagref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_imagref(closest), compRat_imagref(shiftedCenter), halfwidth );
                realRat_sub( compRat_imagref(furthest), compRat_imagref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_imagref(closest), compRat_imagref(closest) );
                if (realRat_cmp( compRat_imagref(closest), compRat_imagref(furthest) ) >0)
                    realRat_neg( compRat_imagref(furthest), compRat_imagref(closest) );
                realRat_set_si( compRat_imagref(closest), 0, 1);
            }
        }
        
        realRat_sub( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
        realRat_add( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
        
        /* if (realRat_sgn( compRat_imagref(closest) )>=0) do nothing */
        if (realRat_sgn( compRat_realref(closest) )<0) {
            if (realRat_sgn( compRat_realref(furthest) )<=0) {
                /* swap closest and furthest*/
                realRat_add( compRat_realref(closest), compRat_realref(shiftedCenter), halfwidth );
                realRat_sub( compRat_realref(furthest), compRat_realref(shiftedCenter), halfwidth );
            } else {
                /* check which one is farthest from zero */
                realRat_abs( compRat_realref(closest), compRat_realref(closest) );
                if (realRat_cmp( compRat_realref(closest), compRat_realref(furthest) ) >0)
                    realRat_neg( compRat_realref(furthest), compRat_realref(closest) );
                realRat_set_si( compRat_realref(closest), 0, 1);
            }
        }
        
        compApp_set_compRat(temp, closest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (left, temp,    CCLUSTER_DEFAULT_PREC );
        compApp_set_compRat(temp, furthest, CCLUSTER_DEFAULT_PREC );
        compApp_abs        (right, temp,    CCLUSTER_DEFAULT_PREC );
        
        /* go to the first annulus with rsup >= left */
        while ( (it!=compAnn_list_end()) && (realApp_lt( compAnn_radSupref(compAnn_list_elmt( it )), left )==1)  ){
            it = compAnn_list_next(it);
        }
        
        /* push annulii satisfying rinf <= right (and rsup >= left) */
        while ( (it!=compAnn_list_end()) 
            && (!(realApp_gt( compAnn_radInfref(compAnn_list_elmt( it )), right )==1))
            ) {
            compAnn_list_push( compBox_annuliref(x,ind), compAnn_list_elmt( it ) );
            it = compAnn_list_next(it);
        }
        
        
        compApp_clear(temp);
        realApp_clear(left);
        realApp_clear(right);
        compRat_clear(shiftedCenter);
        compRat_clear(closest);
        compRat_clear(furthest);
        realRat_clear(halfwidth);
        
//         printf("#actualize anulii: length of list: %d\n", compAnn_list_get_size(compBox_annuliref(x)) );
//     }
}
    
/* assume compBox_annuliref(x) contains anulii in increasing order of their radii */
void compBox_actualize_anulii( compBox_t x, const compBox_t bmother ){
    
    for (int ind=0; ind<GEOMETRY_NB_ANN_PER_BOX; ind++)
        if (compAnn_list_get_size(compBox_annuliref(bmother,ind))>=1){
            compAnn_list_iterator it = compAnn_list_begin(compBox_annuliref(bmother,ind));
            if ( compAnn_centerImref( compAnn_list_elmt(it) )==0 ) {
                compBox_actualize_anulii_compAnn_list_real( x, ind, compBox_annuliref(bmother, ind) );
            } else {
                if (ind==3)
                    compBox_actualize_anulii_compAnn_list_imag_conj( x, ind, compBox_annuliref(bmother, ind) );
                else
                    compBox_actualize_anulii_compAnn_list_imag( x, ind, compBox_annuliref(bmother, ind) );
            }
        }
        
//     if (compAnn_list_get_size(compBox_annuliref(bmother,0))>=1)
//         compBox_actualize_anulii_compAnn_list_real( x, 0, compBox_annuliref(bmother, 0) );
//     if (compAnn_list_get_size(compBox_annuliref(bmother,1))>=1)
//         compBox_actualize_anulii_compAnn_list_real( x, 1, compBox_annuliref(bmother, 1) );
//     if (compAnn_list_get_size(compBox_annuliref(bmother,2))>=1)
//         compBox_actualize_anulii_compAnn_list_imag( x, 2, compBox_annuliref(bmother, 2) );
//     if (compAnn_list_get_size(compBox_annuliref(bmother,3))>=1)
//         compBox_actualize_anulii_compAnn_list_imag_conj( x, 3, compBox_annuliref(bmother, 3) );
}
