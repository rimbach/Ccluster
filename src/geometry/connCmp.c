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

#include "geometry/connCmp.h"

void connCmp_init(connCmp_t x){
    compBox_list_init(connCmp_boxesref(x));
    realRat_init(connCmp_widthref(x));
    realRat_init(connCmp_infReref(x));
    realRat_init(connCmp_supReref(x));
    realRat_init(connCmp_infImref(x));
    realRat_init(connCmp_supImref(x));
    fmpz_init(connCmp_nwSpdref(x));
    connCmp_nSolsref(x) = -1;
    fmpz_set_si(connCmp_nwSpdref(x), 4);
    connCmp_appPrref(x) = CCLUSTER_DEFAULT_PREC;
    connCmp_newSuref(x) = 0;
    connCmp_isSepref(x) = 0;
}

void connCmp_init_compBox(connCmp_t x, compBox_t b){
    connCmp_init(x);
    compBox_list_push(connCmp_boxesref(x), b);
    realRat_set(connCmp_widthref(x), compBox_bwidthref(b));
    
    realRat_t halfwidth;
    realRat_init(halfwidth);
    realRat_set_si(halfwidth, 1,2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(b));
    
    realRat_sub(connCmp_infReref(x), compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_add(connCmp_supReref(x), compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_sub(connCmp_infImref(x), compRat_imagref(compBox_centerref(b)), halfwidth);
    realRat_add(connCmp_supImref(x), compRat_imagref(compBox_centerref(b)), halfwidth);
    
    realRat_clear(halfwidth);
}

void connCmp_clear(connCmp_t x){
    compBox_list_clear(connCmp_boxesref(x));
    realRat_clear(connCmp_widthref(x));
    realRat_clear(connCmp_infReref(x));
    realRat_clear(connCmp_supReref(x));
    realRat_clear(connCmp_infImref(x));
    realRat_clear(connCmp_supImref(x));
    fmpz_clear(connCmp_nwSpdref(x));
}

void connCmp_clear_for_tables(connCmp_t x){
    compBox_list_clear_for_tables(connCmp_boxesref(x));
    realRat_clear(connCmp_widthref(x));
    realRat_clear(connCmp_infReref(x));
    realRat_clear(connCmp_supReref(x));
    realRat_clear(connCmp_infImref(x));
    realRat_clear(connCmp_supImref(x));
    fmpz_clear(connCmp_nwSpdref(x));
}

void connCmp_set(connCmp_t dest, const connCmp_t src){
    realRat_set( connCmp_widthref(dest), connCmp_widthref(src) );
    realRat_set( connCmp_infReref(dest), connCmp_infReref(src) );
    realRat_set( connCmp_supReref(dest), connCmp_supReref(src) );
    realRat_set( connCmp_infImref(dest), connCmp_infImref(src) );
    realRat_set( connCmp_supImref(dest), connCmp_supImref(src) );
    connCmp_nSolsref(dest) = connCmp_nSolsref(src);
    fmpz_set( connCmp_nwSpdref(dest), connCmp_nwSpdref(src) );
    connCmp_appPrref(dest) = connCmp_appPrref(src);
    connCmp_newSuref(dest) = connCmp_newSuref(src);
    connCmp_isSepref(dest) = connCmp_isSepref(src);
    
    /* copy the boxes */
    compBox_ptr nBox;
    
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (src) );
    while (it!=compBox_list_end()) {
        
        nBox = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
        compBox_init(nBox);
        compBox_set(nBox, compBox_list_elmt(it));
        compBox_list_push(connCmp_boxesref (dest), nBox);
        
        it = compBox_list_next(it);
    }
}

connCmp_ptr connCmp_copy(connCmp_t src){
    connCmp_ptr nCmp;
    nCmp = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
    connCmp_init(nCmp);
    connCmp_set(nCmp, src);
    return nCmp;
}

void connCmp_initiali_nwSpd(connCmp_t x){
    fmpz_set_si(connCmp_nwSpdref(x), 4);
}

void connCmp_initiali_nwSpd_connCmp(connCmp_t x, const connCmp_t src){
    fmpz_set(connCmp_nwSpdref(x), connCmp_nwSpdref(src));
}

void connCmp_increase_nwSpd(connCmp_t x){
    fmpz_mul( connCmp_nwSpdref(x), connCmp_nwSpdref(x), connCmp_nwSpdref(x));
}

void connCmp_decrease_nwSpd(connCmp_t x){
    
    if (fmpz_cmp_si(connCmp_nwSpdref(x), 4)>0){ /* means that connCmp_nwSpdref(x) > 4*/
        fmpz_sqrt(connCmp_nwSpdref(x), connCmp_nwSpdref(x));
    }
    else
        fmpz_set_si(connCmp_nwSpdref(x), 4);
}

slong connCmp_getDepth(const connCmp_t c, const compBox_t initialBox){
    realRat_t depth;
    realRat_init(depth);
    
    realRat_div( depth, compBox_bwidthref(initialBox),  connCmp_widthref(c));
    slong res = fmpz_clog_ui( realRat_numref(depth), 2);
    
    realRat_clear(depth);
    return res;
}

int connCmp_is_confined(const connCmp_t c, const compBox_t initialBox){
    int res;
    compBox_t inflatedBox;
    realRat_t factor;
    compBox_init(inflatedBox);
    realRat_init(factor);
    
    realRat_set_si(factor, 5,4);
    compBox_inflate_realRat(inflatedBox, initialBox, factor);
    res = connCmp_is_strictly_in_compBox( c, inflatedBox);
    
    compBox_clear(inflatedBox);
    realRat_clear(factor);
    return res;
}


void connCmp_fprint(FILE * file, const connCmp_t x){
    fprintf(file, "Conn Comp with %d boxes of width ", connCmp_nb_boxes(x));
    realRat_fprint(file, connCmp_widthref(x));
    fprintf(file, "; hull: [ ");
    realRat_fprint(file, connCmp_infReref(x));
    fprintf(file, ", ");
    realRat_fprint(file, connCmp_supReref(x));
    fprintf(file, " ] + i[ ");
    realRat_fprint(file, connCmp_infImref(x));
    fprintf(file, ", ");
    realRat_fprint(file, connCmp_supImref(x));
    fprintf(file, " ]");
    fprintf(file, "; nSols: %d", connCmp_nSols(x));
    fprintf(file, "; nwSpd: ");
    fmpz_fprint(file, connCmp_nwSpdref(x));
    fprintf(file, "; appPr: %ld", connCmp_appPr(x));
}

void connCmp_insert_compBox(connCmp_t x, compBox_t b){
    if (connCmp_is_empty(x)) {
        connCmp_clear(x);
        connCmp_init_compBox(x, b);
        return;
    }
    
    compBox_list_insert_sorted( connCmp_boxesref(x), b);
    
    realRat_t halfwidth;
    realRat_init(halfwidth);
    realRat_set_si(halfwidth, 1,2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(b));
    
    realRat_t bound;
    realRat_init(bound);
    
    realRat_sub(bound, compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_min_2_realRat(connCmp_infReref(x), bound);
    
    realRat_add(bound, compRat_realref(compBox_centerref(b)), halfwidth);
    realRat_max_2_realRat(connCmp_supReref(x), bound);
    
    realRat_sub(bound, compRat_imagref(compBox_centerref(b)), halfwidth);
    realRat_min_2_realRat(connCmp_infImref(x), bound);
    
    realRat_add(bound, compRat_imagref(compBox_centerref(b)), halfwidth);
    realRat_max_2_realRat(connCmp_supImref(x), bound);
    
    realRat_clear(halfwidth);
    realRat_clear(bound);
}

/* Preconditions: cc1, cc2 are not empty, have the same width and are connected */
/* Specifs      : set cc1 to the union of cc1 and cc2 */
/*              : let cc2 empty but do not clear it   */
void connCmp_merge_2_connCmp( connCmp_t cc1, connCmp_t cc2 ){
//     printf("merge, nb boxes in cc1: %d, nb boxes in cc2: %d\n", connCmp_nb_boxes(cc1),connCmp_nb_boxes(cc2) );
//     if( (connCmp_nb_boxes(cc1)==0 )||(connCmp_nb_boxes(cc2)==0 ) )
//         printf("preconditions not satisfied!\n");
    /* sets hull */
    realRat_min_2_realRat(connCmp_infReref(cc1), connCmp_infReref(cc2));
    realRat_max_2_realRat(connCmp_supReref(cc1), connCmp_supReref(cc2));
    realRat_min_2_realRat(connCmp_infImref(cc1), connCmp_infImref(cc2));
    realRat_max_2_realRat(connCmp_supImref(cc1), connCmp_supImref(cc2));
    
    /* save cc1.boxes in temp and sets cc1.boxes to an empty list */
    compBox_list_t temp;
    compBox_list_init(temp);
    compBox_list_swap(temp, connCmp_boxesref(cc1));
    
    /* fill cc1.boxes while preserving order */
    compBox_ptr b1, b2;
    
    while ( (!compBox_list_is_empty(temp)) && (!compBox_list_is_empty( connCmp_boxesref(cc2) )) ) {
        b1 = compBox_list_first( temp );
        b2 = compBox_list_first( connCmp_boxesref(cc2) );
        if (compBox_isless(b1, b2))
            compBox_list_push( connCmp_boxesref(cc1), compBox_list_pop( temp ));
        else
            compBox_list_push( connCmp_boxesref(cc1), compBox_list_pop( connCmp_boxesref(cc2) ));
    }
    /* now either temp or cc2.boxes are empty: fill cc1.boxes with the remaining boxes */
    while (!compBox_list_is_empty(temp))
        compBox_list_push( connCmp_boxesref(cc1), compBox_list_pop( temp ));
    
    while (!compBox_list_is_empty(connCmp_boxesref(cc2)))
        compBox_list_push( connCmp_boxesref(cc1), compBox_list_pop( connCmp_boxesref(cc2) ));
    
    compBox_list_clear(temp);
/*     printf("merge, nb boxes in cc1: %d, nb boxes in cc2: %d\n", connCmp_nb_boxes(cc1),connCmp_nb_boxes(cc2) );*/
}

/* specifs      : define diam(cc) as the max of the width of the real part and the one of the imag part */
void connCmp_diameter( realRat_t diam, const connCmp_t cc){
    realRat_t width;
    realRat_init(width);
    
    realRat_sub( width, connCmp_supReref(cc), connCmp_infReref(cc) );
    realRat_abs( width, width );
    realRat_sub( diam,  connCmp_supImref(cc), connCmp_infImref(cc) );
    realRat_abs( diam,  diam );
    realRat_max_2_realRat(diam, width);
    
    realRat_clear(width);
}

/* ordering */
/* specifs      : cc1<=cc2 <=> diam(cc1)<=diam(cc2) */
int connCmp_isless ( const connCmp_t cc1, const connCmp_t cc2 ){
    realRat_t diam1, diam2;
    int res;
    realRat_init(diam1);
    realRat_init(diam2);
    
    connCmp_diameter( diam1, cc1);
    connCmp_diameter( diam2, cc2);
    res = (realRat_cmp(diam1, diam2) <=0);
    
    realRat_clear(diam1);
    realRat_clear(diam2);
    
    return res;
}

int connCmp_isgreater ( const connCmp_t cc1, const connCmp_t cc2 ){
    realRat_t diam1, diam2;
    int res;
    realRat_init(diam1);
    realRat_init(diam2);
    
    connCmp_diameter( diam1, cc1);
    connCmp_diameter( diam2, cc2);
    res = (realRat_cmp(diam1, diam2) >=0);
    
    realRat_clear(diam1);
    realRat_clear(diam2);
    
    return res;
}

void connCmp_componentBox( compBox_t res, const connCmp_t cc, const compBox_t initialBox){
    
    realRat_t width, halfwidthenlarged, frac;
    compRat_t RightBottomBordEn; /*rightbotton corner of enlarged box*/
    compRat_t RightBottomBordCb; /*rightbotton corner of component box*/
    
    realRat_init(frac);
    realRat_init(width);
    realRat_init(halfwidthenlarged);
    compRat_init(RightBottomBordEn);
    compRat_init(RightBottomBordCb);
    
    realRat_set_si(frac, 5, 4);
    compBox_inflate_realRat(res, initialBox, frac);
    
    connCmp_diameter(width, cc);
    realRat_set_si(frac, 1, 2);
    
    realRat_mul(halfwidthenlarged, frac, compBox_bwidthref(res));
    
    realRat_add( compRat_realref(RightBottomBordEn), compRat_realref(compBox_centerref(res)), halfwidthenlarged);
    realRat_sub( compRat_imagref(RightBottomBordEn), compRat_imagref(compBox_centerref(res)), halfwidthenlarged);
    realRat_add( compRat_realref(RightBottomBordCb), connCmp_infReref(cc), width);
    realRat_sub( compRat_imagref(RightBottomBordCb), connCmp_supImref(cc), width);
    
    realRat_min_2_realRat( compRat_realref(RightBottomBordCb), compRat_realref(RightBottomBordEn));
    realRat_max_2_realRat( compRat_imagref(RightBottomBordCb), compRat_imagref(RightBottomBordEn));
    
    realRat_set( compBox_bwidthref(res), width );
    realRat_mul(width, frac, width);
    realRat_sub( compRat_realref(compBox_centerref(res)), compRat_realref(RightBottomBordCb), width);
    realRat_add( compRat_imagref(compBox_centerref(res)), compRat_imagref(RightBottomBordCb), width);
        
    realRat_clear(frac);
    realRat_clear(width);
    realRat_clear(halfwidthenlarged);
    compRat_clear(RightBottomBordEn);
    compRat_clear(RightBottomBordCb);
}

void connCmp_risolate_componentBox( compBox_t res, const connCmp_t cc, const compBox_t initialBox){
    
    realRat_ptr centerRe = compRat_realref(compBox_centerref(res));
    realRat_ptr centerIm = compRat_imagref(compBox_centerref(res));
    realRat_ptr width    = compBox_bwidthref(res);
    
    realRat_set_si(centerIm,0,1);
    realRat_set_si(centerRe,1,2);
    realRat_sub( width, connCmp_supReref(cc), connCmp_infReref(cc));
    realRat_mul( centerRe, centerRe, width );
    realRat_add( centerRe, connCmp_infReref(cc), centerRe);
    
    
}

int connCmp_is_strictly_in_compBox( const connCmp_t cc, const compBox_t b ){
    
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (!compBox_is_strictly_in(compBox_list_elmt(it) ,b))
            return 0;
        it = compBox_list_next(it);
    }
    return 1;
}

/*Precondition:                                                   */
/*Specification: returns true if interior(cc \cap b)\neq \emptyset*/ 
/*                       false otherwise                          */
int connCmp_intersection_has_non_empty_interior( const connCmp_t cc, const compBox_t b ){
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_intersection_has_non_empty_interior(compBox_list_elmt(it) ,b))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}

/*Precondition:                                            */
/*Specification: returns true if (cc \cap b)\neq \emptyset */
/*                       false otherwise                   */
int connCmp_intersection_is_not_empty( const connCmp_t cc, const compBox_t b ){
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_intersection_is_not_empty(compBox_list_elmt(it) ,b))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}

/*Precondition:  boxes of cc and b HAVE the same width  */                               
/*Specification: returns true if (cc \cap b)\neq \emptyset */
/*                       false otherwise                   */
int connCmp_are_8connected( const connCmp_t cc, const compBox_t b ){
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        if (compBox_are_8connected(compBox_list_elmt(it) ,b))
            return 1;
        it = compBox_list_next(it);
    }
    return 0;
}

void connCmp_find_point_outside_connCmp( compRat_t res, const connCmp_t cc, const compBox_t initialBox ){
    
    int resOK = 0;
    realRat_t halfwidth, halfwidthenlarged, frac;
    compBox_t componentBox;
    
    realRat_init(frac);
    realRat_init(halfwidth);
    realRat_init(halfwidthenlarged);
    compBox_init(componentBox);
    
    realRat_set_si(frac, 5, 4);
    realRat_mul(halfwidthenlarged, frac, compBox_bwidthref(initialBox) );
    realRat_set_si(frac, 1, 2);
    realRat_mul(halfwidthenlarged, frac, halfwidthenlarged);
    connCmp_componentBox(componentBox, cc, initialBox);
    realRat_mul(halfwidth, frac, connCmp_widthref(cc));
    
    /* try a point in the west of the center of componentBox */
    realRat_add( compRat_realref(res), connCmp_supReref(cc), halfwidth );
    realRat_set( compRat_imagref(res), compRat_imagref( compBox_centerref(componentBox) ) );
    /* compute the distance from the west boundary of component box -> store it in frac*/ 
    realRat_add( frac, compRat_realref( compBox_centerref(initialBox) ), halfwidthenlarged);
    realRat_sub( frac, frac, compRat_realref(res));
    resOK = realRat_cmp( frac, halfwidth);
    
    if (resOK<=0) { /* try a point in the east of the center of componentBox */
        realRat_sub( compRat_realref(res), connCmp_infReref(cc), halfwidth );
        /* compute the distance from the east boundary of component box -> store it in frac*/
        realRat_sub( frac, compRat_realref( compBox_centerref(initialBox) ), halfwidthenlarged);
        realRat_sub( frac, compRat_realref(res), frac);
        resOK = realRat_cmp( frac, halfwidth);
    }
    
    if (resOK<=0) { /* try a point in the north of the center of componentBox */
        realRat_set( compRat_realref(res), compRat_realref( compBox_centerref(componentBox) ) );
        realRat_add( compRat_imagref(res), connCmp_supImref(cc), halfwidth );
        /* compute the distance from the north boundary of component box -> store it in frac*/
        realRat_add( frac, compRat_imagref( compBox_centerref(initialBox) ), halfwidthenlarged);
        realRat_sub( frac, frac, compRat_imagref(res));
        resOK = realRat_cmp( frac, halfwidth);
    }
    
    if (resOK<=0) { /* use a point in the south of the center of componentBox */
        realRat_sub( compRat_imagref(res), connCmp_infImref(cc), halfwidth );
         /* compute the distance from the north boundary of component box -> store it in frac*/
/*         realRat_sub( frac, compRat_imagref( compBox_centerref(initialBox) ), halfwidthenlarged);
         realRat_sub( frac, compRat_imagref(res), frac);
         resOK = realRat_cmp( frac, halfwidth);*/
    }
    
    realRat_clear(frac);
    realRat_clear(halfwidth);
    realRat_clear(halfwidthenlarged);
    compBox_clear(componentBox);
}

void connCmp_risolate_find_point_outside_connCmp( compRat_t res, const connCmp_t cc, const compBox_t initialBox ){
    
    realRat_ptr Rres = compRat_realref(res);
    realRat_ptr Ires = compRat_imagref(res);
    realRat_set_si(Ires, 0,0);
    realRat_sub(Rres, connCmp_supReref(cc), connCmp_infReref(cc));
    realRat_sub(Rres, connCmp_infReref(cc), Rres);
    
}

/* RealCoeffs */
/* Precondition:                                                                               */
/* Specification: returns true only if forall p\in b, Im(p)<0                                  */
int connCmp_is_imaginary_negative_strict               ( const connCmp_t cc  ){
    int res;
    realRat_t zero;
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    
    res = (realRat_cmp(connCmp_supImref(cc), zero) < 0);
    
    realRat_clear(zero);
    return res;
}
/* Precondition:                                                                               */
/* Specification: returns true only if forall p\in b, Im(p)>0                                  */
int connCmp_is_imaginary_positive_strict               ( const connCmp_t cc  ){
    int res;
    realRat_t zero;
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    
    res = (realRat_cmp(connCmp_infImref(cc), zero) > 0);
    
    realRat_clear(zero);
    return res;
}
/* Precondition:                                                                                */
/* Specification: returns true only if forall p\in cc, Im(p)>=0                                 */
int connCmp_is_imaginary_positive               ( const connCmp_t cc  ){
    int res;
    realRat_t zero;
    realRat_init(zero);
    realRat_set_si(zero, 0, 1);
    
    res = (realRat_cmp(connCmp_infImref(cc), zero) >= 0);
    
    realRat_clear(zero);
    return res;
}
/* Precondition: res is initialized                                                            */
/* Specification: set res to the complex conjugate of cc                                       */
void connCmp_set_conjugate                      ( connCmp_t res, const connCmp_t cc  ){
    /* copy the boxes */
    compBox_ptr nBox;
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        
        nBox = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
        compBox_init(nBox);
        compBox_set_conjugate(nBox, compBox_list_elmt(it));
        connCmp_insert_compBox(res, nBox);
        
        it = compBox_list_next(it);
    }
    /* the other properties */
    realRat_set( connCmp_widthref(res), connCmp_widthref(cc) );
    connCmp_nSolsref(res) = connCmp_nSolsref(cc);
}

void connCmp_set_conjugate_closure              ( connCmp_t res, const connCmp_t cc, const compBox_t initBox ){
    connCmp_set( res, cc);
    compBox_ptr nBox;
    
    /* test if initBox is symetric with respect to the real line */
    int is_sym = realRat_is_zero(compRat_imagref(compBox_centerref(initBox)));
    realRat_t shift;
    realRat_init(shift);
    realRat_set(shift, connCmp_infImref(cc) );
    realRat_abs(shift, shift);
    realRat_mul_si(shift,shift,2);
    
    compBox_list_iterator it = compBox_list_begin( connCmp_boxesref (cc) );
    while (it!=compBox_list_end()) {
        
        /* get the conjugate */
        nBox = ( compBox_ptr ) ccluster_malloc (sizeof(compBox));
        compBox_init(nBox);
        
        compBox_set_conjugate(nBox, compBox_list_elmt(it));
        
        if (!is_sym)
            realRat_sub( compRat_imagref(compBox_centerref(nBox)), compRat_imagref(compBox_centerref(nBox)), shift );
//         printf(" box: "); compBox_print(compBox_list_elmt(it)); printf("\n"); 
//         printf(" conjugate: "); compBox_print(nBox); printf("\n");
        /* check if the conjugate is in the list */
        compBox_list_iterator it2 = compBox_list_begin( connCmp_boxesref (res) );
        int isIn = 1;
        while ( (it2!=compBox_list_end()) && (isIn>0) ) {
//             printf(" box in res: "); compBox_print(compBox_list_elmt(it2)); printf("\n");
            isIn = compBox_cmp( nBox, compBox_list_elmt(it2) );
            it2 = compBox_list_next(it2);
        }
        if (isIn==0) { /* the conjugate is already in the CC */
            compBox_clear(nBox);
            ccluster_free(nBox);
        }
        else { /* the conjugate is not already in the CC */
        
            connCmp_insert_compBox(res, nBox);
        }
        
        it = compBox_list_next(it);
        
    }
    
    realRat_clear(shift);
}

/* DEPRECATED
void connCmp_find_point_outside_connCmp_for_julia( realRat_t resRe, realRat_t resIm, const connCmp_t cc, const compBox_t initialBox ){
    compRat_t temp;
    compRat_init(temp);
    connCmp_find_point_outside_connCmp( temp, cc, initialBox );
    realRat_set(resRe, compRat_realref(temp));
    realRat_set(resIm, compRat_imagref(temp));
    compRat_clear(temp);
}
*/
