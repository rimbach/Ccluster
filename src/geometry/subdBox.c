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

#include "geometry/subdBox.h"


void subdBox_quadrisect( compBox_list_t res, const compBox_t b ){
    
    clock_t start = clock();
    
    realRat_t shift, width;
    realRat_init(shift);
    realRat_init(width);
    
    compBox_ptr bNE, bSE, bSW, bNW;
    bNE = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    bSE = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    bSW = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    bNW = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    compBox_init(bNE);
    compBox_init(bSE);
    compBox_init(bSW);
    compBox_init(bNW);
    
    realRat_set_si(shift, 1, 4);
    realRat_set_si(width, 1, 2);
    realRat_mul(shift, shift, compBox_bwidthref(b));
    realRat_mul(width, width, compBox_bwidthref(b));
    
    compBox_set_compRat_realRat_int(bNE, compBox_centerref(b), width, b->nbMSol);
    realRat_add( compRat_realref(compBox_centerref(bNE)), compRat_realref(compBox_centerref(bNE)), shift);
    realRat_add( compRat_imagref(compBox_centerref(bNE)), compRat_imagref(compBox_centerref(bNE)), shift);
    
    compBox_set_compRat_realRat_int(bSE, compBox_centerref(b), width, b->nbMSol);
    realRat_add( compRat_realref(compBox_centerref(bSE)), compRat_realref(compBox_centerref(bSE)), shift);
    realRat_sub( compRat_imagref(compBox_centerref(bSE)), compRat_imagref(compBox_centerref(bSE)), shift);
    
    compBox_set_compRat_realRat_int(bSW, compBox_centerref(b), width, b->nbMSol);
    realRat_sub( compRat_realref(compBox_centerref(bSW)), compRat_realref(compBox_centerref(bSW)), shift);
    realRat_sub( compRat_imagref(compBox_centerref(bSW)), compRat_imagref(compBox_centerref(bSW)), shift);
    
    compBox_set_compRat_realRat_int(bNW, compBox_centerref(b), width, b->nbMSol);
    realRat_sub( compRat_realref(compBox_centerref(bNW)), compRat_realref(compBox_centerref(bNW)), shift);
    realRat_add( compRat_imagref(compBox_centerref(bNW)), compRat_imagref(compBox_centerref(bNW)), shift);
    
    /* root radii */
    clock_t start2 = clock();
    
    compBox_actualize_anulii( bNE, b);
    compBox_actualize_anulii( bSE, b);
    compBox_actualize_anulii( bSW, b);
    compBox_actualize_anulii( bNW, b);
    
    timeIn_actualize_anulii += ( (float) clock() - start2 )/CLOCKS_PER_SEC;
    
    /*printf("bNE: "); compBox_print(bNE); printf("\n");*/
    /*printf("bSE: "); compBox_print(bSE); printf("\n");*/
    /*printf("bSW: "); compBox_print(bSW); printf("\n");*/
    /*printf("bNW: "); compBox_print(bNW); printf("\n");*/
    
    compBox_list_push(res, bSW);
    compBox_list_push(res, bNW);
    compBox_list_push(res, bSE);
    compBox_list_push(res, bNE);
    
    
    realRat_clear(shift);
    realRat_clear(width);
    
    timeIn_quadrisect += ( (float) clock() - start )/CLOCKS_PER_SEC ;
}

void subdBox_quadrisect_with_compDsk( compBox_list_t res, const compBox_t b, const compDsk_t d, const realRat_t nwidth){
    
    realRat_t quo;
    realRat_t half, halfwidth, halfnwidth, temp;
    compRat_t downLeftCorner, upRightCorner;
    compRat_t distInf, distSup, infDsk, supDsk, infCenter, supCenter;
    compRat_t boxCenter;
    
    realRat_init(quo);
    realRat_init(half);
    realRat_init(temp);
    realRat_init(halfwidth);
    realRat_init(halfnwidth);
    compRat_init(downLeftCorner);
    compRat_init(upRightCorner);
    compRat_init(distInf);
    compRat_init(distSup);
    compRat_init(infDsk);
    compRat_init(supDsk);
    compRat_init(infCenter);
    compRat_init(supCenter);
    compRat_init(boxCenter);
    
    realRat_set_si(quo, 1,1);
    realRat_set_si(half, 1,2);
    realRat_mul(halfwidth, half, compBox_bwidthref(b));
    realRat_mul(halfnwidth, half, nwidth);
    
    realRat_sub( compRat_realref(downLeftCorner), compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_sub( compRat_imagref(downLeftCorner), compRat_imagref(compBox_centerref(b)), halfwidth);
    realRat_add( compRat_realref(upRightCorner), compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_add( compRat_imagref(upRightCorner), compRat_imagref(compBox_centerref(b)), halfwidth);
    
    realRat_sub( compRat_realref(infDsk), compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    realRat_sub( compRat_imagref(infDsk), compRat_imagref(compDsk_centerref(d)), compDsk_radiusref(d));
    realRat_add( compRat_realref(supDsk), compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    realRat_add( compRat_imagref(supDsk), compRat_imagref(compDsk_centerref(d)), compDsk_radiusref(d));
    
    realRat_sub( compRat_realref(distInf), compRat_realref(infDsk), compRat_realref(downLeftCorner));
    realRat_sub( compRat_imagref(distInf), compRat_imagref(infDsk), compRat_imagref(downLeftCorner));
    realRat_sub( compRat_realref(distSup), compRat_realref(upRightCorner), compRat_realref(supDsk));
    realRat_sub( compRat_imagref(distSup), compRat_imagref(upRightCorner), compRat_imagref(supDsk));
    
    if (fmpq_sgn( compRat_realref(distInf)) <=0 ){
        realRat_add( compRat_realref(infCenter), compRat_realref(downLeftCorner), halfnwidth);
    }
    else {
        realRat_div(temp, compRat_realref(distInf), nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_add(compRat_realref(infCenter), compRat_realref(downLeftCorner), temp);
    }
    
    if (fmpq_sgn( compRat_imagref(distInf)) <=0 ){
        realRat_add( compRat_imagref(infCenter), compRat_imagref(downLeftCorner), halfnwidth);
    }
    else {
        realRat_div(temp, compRat_imagref(distInf), nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_add(compRat_imagref(infCenter), compRat_imagref(downLeftCorner), temp);
    }
    
    if (fmpq_sgn( compRat_realref(distSup)) <=0 ){
        realRat_sub( compRat_realref(supCenter), compRat_realref(upRightCorner), halfnwidth);
    }
    else {
        realRat_div(temp, compRat_realref(distSup), nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_sub(compRat_realref(supCenter), compRat_realref(upRightCorner), temp);
    }
    
    if (fmpq_sgn( compRat_imagref(distSup)) <=0 ){
        realRat_sub( compRat_imagref(supCenter), compRat_imagref(upRightCorner), halfnwidth);
    }
    else {
        realRat_div(temp, compRat_imagref(distSup), nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_sub(compRat_imagref(supCenter), compRat_imagref(upRightCorner), temp);
    }
    
    compRat_set(boxCenter, infCenter);
    compBox_ptr bnew;
    
    while (realRat_cmp(compRat_realref(boxCenter), compRat_realref(supCenter)) <=0){
        
        while (realRat_cmp(compRat_imagref(boxCenter), compRat_imagref(supCenter)) <=0){
        
            bnew = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
            compBox_init(bnew);
            compBox_set_compRat_realRat_int(bnew, boxCenter, nwidth, b->nbMSol);
            
            if (compBox_intersection_has_non_empty_interior_compDsk (bnew, d)){
                compBox_actualize_anulii( bnew, b);
                compBox_list_push(res, bnew);
            }
            else {
                compBox_clear(bnew);
                ccluster_free(bnew);
            }
            /* compBox_list_push(res, bnew); */
            
            realRat_add( compRat_imagref(boxCenter), compRat_imagref(boxCenter), nwidth);
            
        }
        
        realRat_set(compRat_imagref(boxCenter), compRat_imagref(infCenter)); 
        realRat_add(compRat_realref(boxCenter), compRat_realref(boxCenter), nwidth);
        
    }
    
    realRat_clear(quo);
    realRat_clear(half);
    realRat_clear(temp);
    realRat_clear(halfwidth);
    realRat_clear(halfnwidth);
    compRat_clear(downLeftCorner);
    compRat_clear(upRightCorner);
    compRat_clear(distInf);
    compRat_clear(distSup);
    compRat_clear(infDsk);
    compRat_clear(supDsk);
    compRat_clear(infCenter);
    compRat_clear(supCenter);
    compRat_clear(boxCenter);
    
}

void subdBox_risolate_bisect( compBox_list_t res, const compBox_t b ){
    realRat_t shift, width;
    realRat_init(shift);
    realRat_init(width);
    
    compBox_ptr bE, bW;
    bE = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    bW = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
    compBox_init(bE);
    compBox_init(bW);
    
    realRat_set_si(shift, 1, 4);
    realRat_set_si(width, 1, 2);
    realRat_mul(shift, shift, compBox_bwidthref(b));
    realRat_mul(width, width, compBox_bwidthref(b));
    
    compBox_set_compRat_realRat_int(bE, compBox_centerref(b), width, b->nbMSol);
    realRat_add( compRat_realref(compBox_centerref(bE)), compRat_realref(compBox_centerref(bE)), shift);
    
    compBox_set_compRat_realRat_int(bW, compBox_centerref(b), width, b->nbMSol);
    realRat_sub( compRat_realref(compBox_centerref(bW)), compRat_realref(compBox_centerref(bW)), shift);
    
//     printf("b : "); compBox_print(b ); printf("\n");
//     printf("bE: "); compBox_print(bE); printf("\n");
//     printf("bW: "); compBox_print(bW); printf("\n");
    
    /* root radii */
    compBox_actualize_anulii_risolate( bW, b);
    compBox_actualize_anulii_risolate( bE, b);
    
    compBox_list_push(res, bW);
    compBox_list_push(res, bE);
    
    
    realRat_clear(shift);
    realRat_clear(width);
    
}

void subdBox_risolate_bisect_with_compDsk( compBox_list_t res, const compBox_t b, const compDsk_t d, const realRat_t nwidth){
    
    realRat_t quo;
    realRat_t half, halfwidth, halfnwidth, temp;
    realRat_t leftCorner, rightCorner;
    realRat_t distInf, distSup, infDsk, supDsk, infCenter, supCenter;
    compRat_t boxCenter;
    
    realRat_init(quo);
    realRat_init(half);
    realRat_init(temp);
    realRat_init(halfwidth);
    realRat_init(halfnwidth);
    realRat_init(leftCorner);
    realRat_init(rightCorner);
    realRat_init(distInf);
    realRat_init(distSup);
    realRat_init(infDsk);
    realRat_init(supDsk);
    realRat_init(infCenter);
    realRat_init(supCenter);
    compRat_init(boxCenter);
    
    realRat_set_si(quo, 1,1);
    realRat_set_si(half, 1,2);
    realRat_mul(halfwidth, half, compBox_bwidthref(b));
    realRat_mul(halfnwidth, half, nwidth);
    
    realRat_sub( leftCorner, compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_add( rightCorner, compRat_realref(compBox_centerref(b)), halfwidth);
    
    realRat_sub( infDsk, compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    realRat_add( supDsk, compRat_realref(compDsk_centerref(d)), compDsk_radiusref(d));
    
    realRat_sub( distInf, infDsk, leftCorner);
    realRat_sub( distSup, rightCorner, supDsk);
    
    if (fmpq_sgn( distInf) <=0 ){ /*infDsk is on the left of leftCorner*/
        realRat_add( infCenter, leftCorner, halfnwidth);
    }
    else {
        realRat_div(temp, distInf, nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_add(infCenter, leftCorner, temp);
    }
    
    if (fmpq_sgn( distSup) <=0 ){ /*supDsk is on the right of rightCorner*/
        realRat_sub( supCenter, rightCorner, halfnwidth);
    }
    else {
        realRat_div(temp, distSup, nwidth);
        fmpz_fdiv_q( realRat_numref(quo), realRat_numref(temp), realRat_denref(temp));
        realRat_mul(temp, nwidth, quo);
        realRat_add(temp, temp, halfnwidth);
        realRat_sub(supCenter, rightCorner, temp);
    }
    
    realRat_set(compRat_realref(boxCenter), infCenter);
    realRat_set_si(compRat_imagref(boxCenter), 0,1);
    compBox_ptr bnew;
    
    while (realRat_cmp(compRat_realref(boxCenter), supCenter) <=0){
        
        bnew = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
        compBox_init(bnew);
        compBox_set_compRat_realRat_int(bnew, boxCenter, nwidth, b->nbMSol);
        
        if (compBox_intersection_has_non_empty_interior_compDsk (bnew, d)){
            compBox_actualize_anulii_risolate( bnew, b);
            compBox_list_push(res, bnew);
        }
        else {
            compBox_clear(bnew);
            ccluster_free(bnew);
        } 
        realRat_add(compRat_realref(boxCenter), compRat_realref(boxCenter), nwidth);
        
    }
    
    realRat_clear(quo);
    realRat_clear(half);
    realRat_clear(temp);
    realRat_clear(halfwidth);
    realRat_clear(halfnwidth);
    realRat_clear(leftCorner);
    realRat_clear(rightCorner);
    realRat_clear(distInf);
    realRat_clear(distSup);
    realRat_clear(infDsk);
    realRat_clear(supDsk);
    realRat_clear(infCenter);
    realRat_clear(supCenter);
    compRat_clear(boxCenter);
    
}

/* DEPRECATED 

void subdBox_quadrisect_intersect_compDsk( compBox_list_t res, const compBox_t b, const compDsk_t d){
    
    compBox_list_t temp;
    compBox_list_init(temp);
    compBox_ptr btemp;
    
    subdBox_quadrisect(temp, b);
    while (!compBox_list_is_empty(temp)) {
        btemp = compBox_list_pop(temp);
        if (compBox_intersection_has_non_empty_interior_compDsk (btemp, d)){
            compBox_list_push(res, btemp);
        }
        else {
            compBox_clear(btemp);
            ccluster_free(btemp);
        }
    }
    
    compBox_list_clear(temp);
    
}
*/
