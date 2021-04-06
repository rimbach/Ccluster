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

#include "risolate/risolate.h"

/* assume the cc contains 1 boxes */
slong risolate_connCmp_getZoomRatio_1b( connCmp_t cc, slong prec ) {
    
    compBox_ptr b;
    compAnn_ptr a;
    int sgn;
    slong resleft = 0;
    slong resright = 0;
    
    realRat_t center, width, temp;
    realRat_init(center);
    realRat_init(width);
    realRat_init(temp);
    realApp_t tempapp, boundann;
    realApp_init(tempapp);
    realApp_init(boundann);
    
    b = compBox_list_first(connCmp_boxesref(cc)); /* unique box */
    sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    /* leftmost anulus */
    if (sgn < 0)
        a = compAnn_list_last( compBox_annuli0ref(b) );
    else 
        a = compAnn_list_first( compBox_annuli0ref(b) );
    
    /* if zero annulus return zero */
    if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
        realRat_clear(center);
        realRat_clear(width);
        realRat_clear(temp);
        realApp_clear(tempapp);
        realApp_clear(boundann);
        return 0;
    }
    
    realRat_set(width, compBox_bwidthref(b) );
    realRat_div_ui( width, width, 2 );
    realRat_add(center, compRat_realref(compBox_centerref(b)), width );
    realRat_sub(temp, center, width);
    realApp_set_realRat(tempapp, temp, prec);
    if (sgn<0) realApp_neg(boundann, compAnn_radSupref(a));
    while (  ( (sgn<0) && ( realApp_lt( tempapp, boundann ) == 1 ) )
          || ( (sgn>0) && ( realApp_lt( tempapp, compAnn_radInfref(a) ) == 1 ) ) ) {
        resleft += 1;
        realRat_div_ui( width, width, 2 );
        realRat_sub(temp, center, width);
        realApp_set_realRat(tempapp, temp, prec);
    }
    
    sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    /* rightmost anulus */
    if (sgn < 0)
        a = compAnn_list_first( compBox_annuli0ref(b) );
    else 
        a = compAnn_list_last( compBox_annuli0ref(b) );
    realRat_set(width, compBox_bwidthref(b) );
    realRat_div_ui( width, width, 2 );
    realRat_sub(center, compRat_realref(compBox_centerref(b)), width );
    realRat_add(temp, center, width);
    realApp_set_realRat(tempapp, temp, prec);
    if (sgn<0) realApp_neg(boundann, compAnn_radInfref(a));
    while (  ( (sgn<0) && ( realApp_gt( tempapp, boundann ) == 1 ) )
          || ( (sgn>0) && ( realApp_gt( tempapp, compAnn_radSupref(a) ) == 1 ) ) ) {
        resright += 1;
        realRat_div_ui( width, width, 2 );
        realRat_add(temp, center, width);
        realApp_set_realRat(tempapp, temp, prec);
    }
    
    
    realRat_clear(center);
    realRat_clear(width);
    realRat_clear(temp);
    realApp_clear(tempapp);
    realApp_clear(boundann);
    
    if (resleft>resright)
        return -resleft;
    else 
        return resright;
}

/* assume the cc contains 2 boxes */
slong risolate_connCmp_getZoomRatio_2b( connCmp_t cc, slong prec ) {
    
    compBox_ptr b;
    compAnn_ptr a;
    int sgn;
    slong resleft = 0;
    slong resright = 0;
    
    realRat_t center, width, temp;
    realRat_init(center);
    realRat_init(width);
    realRat_init(temp);
    realApp_t tempapp, boundann;
    realApp_init(tempapp);
    realApp_init(boundann);
    
    b = compBox_list_first(connCmp_boxesref(cc)); /* leftmost box */
    sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    if (sgn < 0)
        a = compAnn_list_last( compBox_annuli0ref(b) );
    else 
        a = compAnn_list_first( compBox_annuli0ref(b) );
    
    /* if zero annulus return zero */
    if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
        realRat_clear(center);
        realRat_clear(width);
        realRat_clear(temp);
        realApp_clear(tempapp);
        realApp_clear(boundann);
        return 0;
    }
        
    realRat_set(width, compBox_bwidthref(b) );
    realRat_div_ui( width, width, 2 );
    realRat_add(center, compRat_realref(compBox_centerref(b)), width );
    realRat_sub(temp, center, width);
    realApp_set_realRat(tempapp, temp, prec);
    if (sgn<0) realApp_neg(boundann, compAnn_radSupref(a));
    while (  ( (sgn<0) && ( realApp_lt( tempapp, boundann ) == 1 ) )
          || ( (sgn>0) && ( realApp_lt( tempapp, compAnn_radInfref(a) ) == 1 ) ) ) {
        resleft += 1;
        realRat_div_ui( width, width, 2 );
        realRat_sub(temp, center, width);
        realApp_set_realRat(tempapp, temp, prec);
    }
    
    b = compBox_list_last(connCmp_boxesref(cc)); /* rightmost box */
    sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    if (sgn < 0)
        a = compAnn_list_first( compBox_annuli0ref(b) );
    else 
        a = compAnn_list_last( compBox_annuli0ref(b) );
    realRat_set(width, compBox_bwidthref(b) );
    realRat_div_ui( width, width, 2 );
    realRat_sub(center, compRat_realref(compBox_centerref(b)), width );
    realRat_add(temp, center, width);
    realApp_set_realRat(tempapp, temp, prec);
    if (sgn<0) realApp_neg(boundann, compAnn_radInfref(a));
    while (  ( (sgn<0) && ( realApp_gt( tempapp, boundann ) == 1 ) )
          || ( (sgn>0) && ( realApp_gt( tempapp, compAnn_radSupref(a) ) == 1 ) ) ) {
        resright += 1;
        realRat_div_ui( width, width, 2 );
        realRat_add(temp, center, width);
        realApp_set_realRat(tempapp, temp, prec);
    }
    
    
    realRat_clear(center);
    realRat_clear(width);
    realRat_clear(temp);
    realApp_clear(tempapp);
    realApp_clear(boundann);
    
    return CCLUSTER_MIN( resleft, resright );
}

void risolate_compBox_getNumbers( int * nbA,  /* nb of annulii intersecting the box */
                                    int * nbA0, /* nb of segments intersecting the box where the nb of real roots is 0 */
                                    int * nbA1, /* nb of segments contained in  b where the nb of real roots is >=1 */
                                    int * nbA2, /* nb of segments contained in 2b where the nb of real roots is >=1 */
                                    const compBox_t b,
                                    slong prec ) {
    *nbA = 0;
    *nbA0 = 0;
    *nbA1 = 0;
    *nbA2 = 0;
    
    realApp_t center, rad, left, right, left2, right2;
    realApp_init( center );
    realApp_init( rad );
    realApp_init( left );
    realApp_init( right );
    realApp_init( left2 );
    realApp_init( right2 );
    /* compute segment b=[binf,bsup]*/
    realApp_set_realRat( center, compRat_realref(compBox_centerref(b)), prec );
    realApp_abs        ( center,   center );
    realApp_set_realRat( rad,    compBox_bwidthref(b), prec );
    realApp_div_si     ( rad,    rad,                  2, prec );
    realApp_sub        ( left,   center,                 rad, prec );
    realApp_add        ( right,  center,                rad, prec );
    /* compute segment 2b=[2binf,2bsup]*/
    realApp_set_realRat( rad,    compBox_bwidthref(b), prec );
    realApp_sub        ( left2,   center,                 rad, prec );
    realApp_add        ( right2,  center,                rad, prec );
    
    int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
    while ( ita != compAnn_list_end() ) {
        compAnn_ptr a = compAnn_list_elmt(ita);
        *nbA+=1;
        if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
                /* it is the zero annulus */
                *nbA1+=1;
        } else {
            if ( ( sgn <= 0 ) && (compAnn_rrInNeref( a ) == 0) ) *nbA0+=1;
            if ( ( sgn >= 0 ) && (compAnn_rrInPoref( a ) == 0) ) *nbA0+=1;
            if ( ( ( sgn <= 0 ) && (compAnn_rrInNeref( a ) >= 1) )
               ||( ( sgn >= 0 ) && (compAnn_rrInPoref( a ) >= 1) ) ) {
               /* check if b contains the segment */
               if (  ( realApp_lt( left, compAnn_radInfref( a ) ) ==1 )
                   &&( realApp_gt( right, compAnn_radSupref( a ) ) ==1 ) ) {
                   *nbA1+=1;
               } else if (  ( realApp_lt( left2, compAnn_radInfref( a ) ) ==1 )
                   &&( realApp_gt( right2, compAnn_radSupref( a ) ) ==1 ) ) {
                   *nbA2+=1;
               }
                   
            }
        }
        ita = compAnn_list_next(ita);
    }
    
    realApp_clear( center );
    realApp_clear( rad );
    realApp_clear( left );
    realApp_clear( right );
    realApp_clear( left2 );
    realApp_clear( right2 );
    
    return;
}

void risolate_compBox_getNumbersForCounting( int * nbA,  /* nb of annulii intersecting the box */
                                    int * nbA0, /* nb of segments intersecting the box where the nb of real roots is 0 */
                                    int * nbA1, /* nb of segments contained in  b where the nb of real roots is =1 */
                                    int * nbA2, /* nb of segments contained in  b where the nb of real roots is >=1 */
                                    const compBox_t b,
                                    slong prec ) {
    *nbA = 0;
    *nbA0 = 0;
    *nbA1 = 0;
    *nbA2 = 0;
    
    realApp_t center, rad, left, right, neginf, negsup;
    realApp_init( center );
    realApp_init( rad );
    realApp_init( left );
    realApp_init( right );
    realApp_init( neginf );
    realApp_init( negsup );
    /* compute segment b=[binf,bsup]*/
    realApp_set_realRat( center, compRat_realref(compBox_centerref(b)), prec );
//     realApp_abs        ( center,   center );
    realApp_set_realRat( rad,    compBox_bwidthref(b), prec );
    realApp_div_si     ( rad,    rad,                  2, prec );
    realApp_sub        ( left,   center,                 rad, prec );
    realApp_add        ( right,  center,                rad, prec );
    
    int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
    compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
    while ( ita != compAnn_list_end() ) {
        
        
        compAnn_ptr a = compAnn_list_elmt(ita);
        /*check if the neg and pos segments are included in the cc*/
        int posContainedInB = ( realApp_lt( left, compAnn_radInfref( a ) ) ==1 )
                                 &&( realApp_gt( right, compAnn_radSupref( a ) ) ==1 );
        realApp_neg( neginf, compAnn_radSupref( a ) );
        realApp_neg( negsup, compAnn_radInfref( a ) );
        int negContainedInB = ( realApp_lt( left, neginf ) ==1 )
                                &&( realApp_gt( right, negsup ) ==1 );
        
        if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
                /* it is the zero annulus */
                *nbA+=1;
                *nbA1+=(compAnn_indMaxref(a) - compAnn_indMinref(a) +1 );
        } else {
            
            *nbA+=1;
            if (sgn==0) *nbA+=1; /* special case where it is the initial box */
            if ( sgn <= 0 ) {/* it is on the left of 0 */
                if (compAnn_rrInNeref( a ) == 0) *nbA0+=1;
                if ( (compAnn_rrInNeref( a ) == 1) && negContainedInB ) *nbA1+=1;
                if ( (compAnn_rrInNeref( a ) == 2) && negContainedInB ) *nbA2+=1;
            }
            if ( sgn >= 0 ) {/* it is on the right of 0 */
                if (compAnn_rrInPoref( a ) == 0) *nbA0+=1;
                if ( (compAnn_rrInPoref( a ) == 1) && posContainedInB ) *nbA1+=1;
                if ( (compAnn_rrInPoref( a ) == 2) && posContainedInB ) *nbA2+=1;
            }    
            
        }
        ita = compAnn_list_next(ita);
    }
    
    realApp_clear( center );
    realApp_clear( rad );
    realApp_clear( left );
    realApp_clear( right );
    realApp_clear( neginf );
    realApp_clear( negsup );
    
    return;
}

// void risolate_connCmp_getNumbers( int * nbA,    /* nb of segments intersecting the cc */
//                                     int * nbA0, /* nb of segments intersecting the cc where the nb of real roots is 0 */
//                                     int * nbA1, /* nb of segments contained in the cc where the nb of real roots is 1 */
//                                     int * nbA2, /* nb of segments contained in the cc where the nb of real roots is >=1 */
//                                     const connCmp_t c,
//                                     cacheApp_t cache, /* need the degree */
//                                     slong prec ) {
//     *nbA = 0;
//     *nbA0 = 0;
//     *nbA1 = 0;
//     *nbA2 = 0;
//     
//     compAnn_ptr acommon = NULL;
//     
//     realApp_t left, right, neginf,negsup;
//     realApp_init( left );
//     realApp_init( right );
//     realApp_init( neginf );
//     realApp_init( negsup );
//     realApp_set_realRat( left, connCmp_infReref(c), prec );
//     realApp_set_realRat( right, connCmp_supReref(c), prec );
//     
//     compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( c ) );
//     while ( itb != compBox_list_end() ) {
//         printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers \n");
//         compBox_ptr b = compBox_list_elmt(itb);
//         int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
//         printf("# --- new box: "); compBox_print(b); printf("\n");
//         /* check if it is a common annulus */
//         compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
//         
//         if (acommon !=NULL) {
//             printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers, common annulus: ");
//             compAnn_printd(acommon, 10); printf("\n");
//             printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers, cur annulus: ");
//             compAnn_printd(compAnn_list_elmt(ita), 10); printf("\n");
//         }
//         if ( (ita != NULL) && ( compAnn_list_elmt(ita) == acommon )
//                            && ( compAnn_indMaxref(compAnn_list_elmt(ita)) != cacheApp_getDegree(cache) ) 
//            ){
//            
//             printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers \n");
//             printf("# common annulus\n");
//             ita = compAnn_list_next(ita);
//         }
//             
//         while ( ita != compAnn_list_end() ) {
//             compAnn_ptr a = compAnn_list_elmt(ita);
//             printf("# ------ new annulus"); compAnn_printd(a, 10); printf("\n");
//             /*check if the neg and pos segments are included in the cc*/
//             int posContainedInCC = ( realApp_lt( left, compAnn_radInfref( a ) ) ==1 )
//                                  &&( realApp_gt( right, compAnn_radSupref( a ) ) ==1 );
//             realApp_neg( neginf, compAnn_radSupref( a ) );
//             realApp_neg( negsup, compAnn_radInfref( a ) );
//             int negContainedInCC = ( realApp_lt( left, neginf ) ==1 )
//                                  &&( realApp_gt( right, negsup ) ==1 );
//                 
//             if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
//                 printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers zero annulus \n");
//                 /* it is the zero annulus */
//                 if (sgn<=0) { /*add it only once*/
// //                     *nbA+=1;
//                     *nbA+= (compAnn_indMaxref(a) - compAnn_indMinref(a) +1 );
//                     *nbA1+= (compAnn_indMaxref(a) - compAnn_indMinref(a) +1 );
// //                     printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers, zero annulus \n");
//                     printf("# nbA: %d, nbA0: %d, nbA1: %d, nbA2: %d\n", *nbA, *nbA0, *nbA1, *nbA2);
//                 }
//             } else {
//                 printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers not zero annulus \n");
//                 *nbA+=1;
//                 if (sgn==0) *nbA+=1; /* special case where it is the initial box */
//                 if ( sgn <= 0 ) {
//                     if (compAnn_rrInNeref( a ) == 0) *nbA0+=1;
//                     if ( (compAnn_rrInNeref( a ) == 1) && negContainedInCC ) *nbA1+=1;
//                     if ( (compAnn_rrInNeref( a ) == 2) && negContainedInCC ) *nbA2+=1;
//                 }
//                 if ( sgn >= 0 ) {
//                     if (compAnn_rrInPoref( a ) == 0) *nbA0+=1;
//                     if ( (compAnn_rrInPoref( a ) == 1) && posContainedInCC ) *nbA1+=1;
//                     if ( (compAnn_rrInPoref( a ) == 2) && posContainedInCC ) *nbA2+=1;
//                 }
//                 printf("# nbA: %d, nbA0: %d, nbA1: %d, nbA2: %d\n", *nbA, *nbA0, *nbA1, *nbA2);
//             }
//             
//             ita = compAnn_list_next(ita);
//         }
//         /* actualize common annulus */
//         if (sgn<0) {
//             acommon = compAnn_list_first(compBox_annuliref( b, 0 ));
//         } else if (sgn>0) {
//             acommon = compAnn_list_last(compBox_annuliref( b, 0 ));
//         } else {
//             acommon = NULL; //it is the first box of the subdivision tree
//         }
//         itb = compBox_list_next(itb);
//     }
//     
//     realApp_clear(left);
//     realApp_clear(right);
//     realApp_clear( neginf );
//     realApp_clear( negsup );
//     return;
// }

void risolate_connCmp_getNumbers( int * nbA,    /* nb of segments intersecting the cc */
                                    int * nbA0, /* nb of segments intersecting the cc where the nb of real roots is 0 */
                                    int * nbA1, /* nb of segments contained in the cc where the nb of real roots is 1 */
                                    int * nbA2, /* nb of segments contained in the cc where the nb of real roots is >=1 */
                                    const connCmp_t c,
                                    cacheApp_t cache, /* need the degree */
                                    slong prec ) {
    *nbA = 0;
    *nbA0 = 0;
    *nbA1 = 0;
    *nbA2 = 0;
    
    compAnn_ptr acommon = NULL;
    
    realApp_t left, right, neginf,negsup;
    realApp_init( left );
    realApp_init( right );
    realApp_init( neginf );
    realApp_init( negsup );
    realApp_set_realRat( left, connCmp_infReref(c), prec );
    realApp_set_realRat( right, connCmp_supReref(c), prec );
    
    compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( c ) );
    while ( itb != compBox_list_end() ) {
        
//         printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers \n");
        
        compBox_ptr b = compBox_list_elmt(itb);
        int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
//         printf("# --- new box: "); compBox_print(b); printf("\n");
        /* check if it is a common annulus */
        compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
            
        while ( ita != compAnn_list_end() ) {
            compAnn_ptr a = compAnn_list_elmt(ita);
//             printf("# ------ new annulus"); compAnn_printd(a, 10); printf("\n");
            /* check if it is a common annulus */
            if ( (acommon != NULL) && ( compAnn_list_elmt(ita) == acommon )
                           && ( compAnn_indMaxref(compAnn_list_elmt(ita)) != cacheApp_getDegree(cache) ) 
               ){
           
//                     printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers \n");
//                     printf("# it is the common annulus\n");
                    ita = compAnn_list_next(ita);
                }
            
            if (ita == compAnn_list_end())
                continue;
            /*check if the neg and pos segments are included in the cc*/
            int posContainedInCC = ( realApp_lt( left, compAnn_radInfref( a ) ) ==1 )
                                 &&( realApp_gt( right, compAnn_radSupref( a ) ) ==1 );
            realApp_neg( neginf, compAnn_radSupref( a ) );
            realApp_neg( negsup, compAnn_radInfref( a ) );
            int negContainedInCC = ( realApp_lt( left, neginf ) ==1 )
                                 &&( realApp_gt( right, negsup ) ==1 );
                
            if ( realApp_is_zero(compAnn_radInfref(a)) && realApp_is_zero(compAnn_radSupref(a)) ) {
//                 printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers zero annulus \n");
                /* it is the zero annulus */
                if (sgn<=0) { /*add it only once*/
//                     *nbA+=1;
                    *nbA+= (compAnn_indMaxref(a) - compAnn_indMinref(a) +1 );
                    *nbA1+= (compAnn_indMaxref(a) - compAnn_indMinref(a) +1 );
//                     printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers, zero annulus \n");
//                     printf("# nbA: %d, nbA0: %d, nbA1: %d, nbA2: %d\n", *nbA, *nbA0, *nbA1, *nbA2);
                }
            } else {
//                 printf("# risolate_rootRadii.c: risolate_connCmp_getNumbers not zero annulus \n");
                *nbA+=1;
                if (sgn==0) *nbA+=1; /* special case where it is the initial box */
                if ( sgn <= 0 ) {
                    if (compAnn_rrInNeref( a ) == 0) *nbA0+=1;
                    if ( (compAnn_rrInNeref( a ) == 1) && negContainedInCC ) *nbA1+=1;
                    if ( (compAnn_rrInNeref( a ) == 2) && negContainedInCC ) *nbA2+=1;
                }
                if ( sgn >= 0 ) {
                    if (compAnn_rrInPoref( a ) == 0) *nbA0+=1;
                    if ( (compAnn_rrInPoref( a ) == 1) && posContainedInCC ) *nbA1+=1;
                    if ( (compAnn_rrInPoref( a ) == 2) && posContainedInCC ) *nbA2+=1;
                }
//                 printf("# nbA: %d, nbA0: %d, nbA1: %d, nbA2: %d\n", *nbA, *nbA0, *nbA1, *nbA2);
            }
            
            ita = compAnn_list_next(ita);
        }
        /* actualize common annulus */
        if (sgn<0) {
            acommon = compAnn_list_first(compBox_annuliref( b, 0 ));
        } else if (sgn>0) {
            acommon = compAnn_list_last(compBox_annuliref( b, 0 ));
        } else {
            acommon = NULL; //it is the first box of the subdivision tree
        }
        itb = compBox_list_next(itb);
    }
    
    realApp_clear(left);
    realApp_clear(right);
    realApp_clear( neginf );
    realApp_clear( negsup );
    return;
}

slong risolate_discard_compBox_list_rootRadii( compBox_list_t boxes, 
                                                       compBox_list_t bDiscarded,
                                                       connCmp_t cc,
                                                       cacheApp_t cache, 
                                                       slong prec, 
                                                       metadatas_t meta){
    
    tstar_res res;
    res.appPrec = prec;
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
    compBox_ptr btemp;
    
//     if (metadatas_getVerbo(meta)>3) {
//             printf("#risolate_discard_compBox_list_rootRadii: begin \n");
//     }
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        risolate_compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
//         if (metadatas_getVerbo(meta)>3) {
//             printf("#---depth: %d\n", (int) depth);
//             printf("#---Box:       "); compBox_print(btemp); printf("\n");
//         }
        
        int nbA;  /* nb of annulii intersecting the box */
        int nbA0; /* nb of segments intersecting the box where the nb of real roots is 0 */
        int nbA1; /* nb of segments contained in  b where the nb of real roots is >=1 */
        int nbA2; /* nb of segments contained in 2b where the nb of real roots is >=1 */
        risolate_compBox_getNumbers( &nbA, &nbA0, &nbA1, &nbA2, btemp, CCLUSTER_DEFAULT_PREC);
        if (nbA==nbA0) {
            if (metadatas_haveToCount(meta)){
                metadatas_add_discarded( meta, depth);
            }
//             if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
//             } else {
//                 compBox_list_push(bDiscarded, btemp);
//             }
            continue;
        } else if ( (nbA==1) && ( (nbA1==1)||(nbA2==1) ) ) {
            if (nbA1==1){
                compAnn_ptr atemp = compAnn_list_first(compBox_annuli0ref(btemp));
                btemp->nbMSol = compAnn_indMaxref(atemp) - compAnn_indMinref(atemp) +1;
            }
//             if (metadatas_getVerbo(meta)>3) {
//                 printf("#---bisect compBox list: continue because\n");
//                 printf("#---bisect compBox list: nbA : %d\n", nbA);
//                 printf("#---bisect compBox list: nbA0: %d\n", nbA0);
//                 printf("#---bisect compBox list: nbA1: %d\n", nbA1);
//                 printf("#---bisect compBox list: nbA2: %d\n", nbA2);
//             }
            compBox_list_push(ltemp, btemp);
            continue;
        } else if ( (nbA1>=1)||(nbA2>=1) ) {
//             if (metadatas_getVerbo(meta)>3) {
//                 printf("#---bisect compBox list: continue because\n");
//                 printf("#---bisect compBox list: nbA : %d\n", nbA);
//                 printf("#---bisect compBox list: nbA0: %d\n", nbA0);
//                 printf("#---bisect compBox list: nbA1: %d\n", nbA1);
//                 printf("#---bisect compBox list: nbA2: %d\n", nbA2);
//             }
            compBox_list_push(ltemp, btemp);
            continue;
        }
        
        
        res.nbOfSol = -2;
        /* deflation */
        if (metadatas_useDeflation(meta)) {
                if (connCmp_isDefref(cc)==1) {
                    slong precsave = res.appPrec;
                    
//                     if (metadatas_getVerbo(meta)>=3)
//                         printf("---bisect compBox list: deflation is defined\n");
                    
                    res = deflate_tstar_test( cc, cache, bdisk, connCmp_nSolsref(cc), 1, res.appPrec, meta);
                    if (metadatas_getVerbo(meta)>=3)
                        printf("---tstar with deflation        : nbSols: %d, prec: %ld \n", res.nbOfSol, res.appPrec);
                    if (res.nbOfSol == -2)
                        res.appPrec = precsave;
                }
        } 
        if (res.nbOfSol == -2) {
            res = tstar_real_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta); 
        }
        
        
        if (res.nbOfSol==0) {
            if (metadatas_haveToCount(meta)){
                metadatas_add_discarded( meta, depth);
            }
            if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
            } else {
                compBox_list_push(bDiscarded, btemp);
            }
        }
        
        else{
            
                if (res.nbOfSol>0)
                    btemp->nbMSol = res.nbOfSol;
                compBox_list_push(ltemp, btemp);
        }
    }   
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
//     if (metadatas_getVerbo(meta)>3) {
//         printf("risolate_discard_compBox_list_rootRadii: end \n");
//     }
    
    return res.appPrec;
    
    
}

void risolate_bisect_connCmp_rootRadii( connCmp_list_t dest, 
                                                 connCmp_t cc, 
                                                 connCmp_list_t discardedCcs,
                                                 compBox_list_t bDiscarded, 
                                                 cacheApp_t cache, 
                                                 metadatas_t meta, 
                                                 slong nbThreads){
    
    slong prec = connCmp_appPr(cc);
    compBox_list_t subBoxes;
    connCmp_list_t ltemp;
    compBox_list_init(subBoxes);
    connCmp_list_init(ltemp);
    
    compBox_ptr btemp;
    connCmp_ptr ctemp;
    
    while (!connCmp_is_empty(cc)) {
        btemp = connCmp_pop(cc);
        subdBox_risolate_bisect( subBoxes, btemp );
        compBox_clear(btemp);
        ccluster_free(btemp);
    }

    prec = risolate_discard_compBox_list_rootRadii( subBoxes, bDiscarded, cc, cache, prec, meta);
    
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    int specialFlag = 1;
    if (connCmp_list_get_size(ltemp) == 1)
        specialFlag = 0;
    
    slong nprec; 
    if (prec == connCmp_appPrref(cc)) {
        nprec = CCLUSTER_MAX(prec/2,CCLUSTER_DEFAULT_PREC);
//         printf("decrease precision\n");
    }
    else 
        nprec = prec;
    
    
    while (!connCmp_list_is_empty(ltemp)){
        ctemp = connCmp_list_pop(ltemp);
        
        if (connCmp_intersection_is_not_empty(ctemp, metadatas_initBref(meta))){
            connCmp_appPrref(ctemp) = nprec;
            if (specialFlag)
                connCmp_initiali_nwSpd(ctemp);
            else {
                connCmp_initiali_nwSpd_connCmp(ctemp, cc);
                connCmp_decrease_nwSpd(ctemp);
                /* copy the number of sols */
                connCmp_nSolsref(ctemp) = connCmp_nSolsref(cc);
                /* test */
                connCmp_isSep(ctemp) = connCmp_isSep(cc);
                /*end test */
                
                if (metadatas_useDeflation(meta)) {
                    if (connCmp_isDefref(cc)==1) {
//                         printf("---bisect conn comp: deflation is defined; copy deflation data\n");
                        deflate_connCmp_init(ctemp);
                        deflate_copy(ctemp, cc);
                    }
                }
            }
            connCmp_list_push(dest, ctemp);
        }
        else {
            connCmp_appPrref(ctemp) = prec;
            connCmp_list_push(discardedCcs, ctemp);
        }
    }
    
    compBox_list_clear(subBoxes);
    connCmp_list_clear(ltemp);
    
    
}

void risolate_bisect_connCmp_with_ratio( connCmp_list_t dest, 
                                         connCmp_t cc, 
                                         slong ratio){
    
    compBox_ptr btemp;
    connCmp_ptr ctemp;
    compBox_list_t subBoxes;
    connCmp_list_t ltemp;
    compBox_list_init(subBoxes);
    connCmp_list_init(ltemp);
    
    if (compBox_list_get_size(connCmp_boxesref(cc))==1) {
        btemp = connCmp_pop(cc);
        int side = 1;
        if (ratio<0) {
            side = -1;
            ratio = -ratio;
        }
        subdBox_risolate_bisect_with_ratio( subBoxes, btemp, ratio, side );
        compBox_clear(btemp);
        ccluster_free(btemp);
    } else {
    /* assume ccur has two boxes */
        btemp = connCmp_pop(cc);
        subdBox_risolate_bisect_with_ratio( subBoxes, btemp, ratio, -1 );
        compBox_clear(btemp);
        ccluster_free(btemp);
        
        btemp = connCmp_pop(cc);
        subdBox_risolate_bisect_with_ratio( subBoxes, btemp, ratio, 1 );
        compBox_clear(btemp);
        ccluster_free(btemp);
    }
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    
    /* ltemp contains only one cc */
    ctemp = connCmp_list_pop(ltemp);
    connCmp_appPrref(ctemp) = connCmp_appPr(cc);
    connCmp_initiali_nwSpd_connCmp(ctemp, cc);
    connCmp_decrease_nwSpd(ctemp);
    /* copy the number of sols */
    connCmp_nSolsref(ctemp) = connCmp_nSolsref(cc);
    /* test */
    connCmp_isSep(ctemp) = connCmp_isSep(cc);
    /*end test */
    connCmp_list_push(dest, ctemp);
    
    compBox_list_clear(subBoxes);
    connCmp_list_clear(ltemp);
    
    
}

/* return values: -2: no conclusion => run Pellet test */
/*                -1: 2B\(1/2)B contains at least one real root => skip Pellet test*/
/*                 0: ??? */
/*               >=1: B contains that numb. of sol => skip Pellet test */ 
int risolate_rootRadii_countingTest( connCmp_t ccur,
                                     slong prec,
                                     cacheApp_t cache, 
                                     metadatas_t meta ) {
    
    int nbA, nbA0, nbA1, nbA2;
    risolate_connCmp_getNumbers( &nbA, &nbA0, &nbA1, &nbA2, ccur, cache, prec );
    /* nbA:  nb of segments intersecting the cc */
    /* nbA0: nb of segments intersecting the cc where the nb of real roots is 0 */
    /* nbA1: nb of segments contained in the cc where the nb of real roots is 1 */
    /* nbA2: nb of segments contained in the cc where the nb of real roots is >=1 */
//     if (metadatas_getVerbo(meta)>3) {
//         printf("# risolate_rootRadii.c: risolate_rootRadii_countingTest \n");
//         printf("# nbA: %d, nbA0: %d, nbA1: %d, nbA2: %d\n", nbA, nbA0, nbA1, nbA2);
//     }
    
    int determined_by_rootRadii = (nbA==(nbA0+nbA1)); /* in this case nb of Sols = nbA1 */
    
    if (determined_by_rootRadii)
        return nbA1;
    
    /* check if there are roots in 2B\(1/2)B where B is the component box */
    int nbAL, nbAL0, nbAL1, nbAL2;
    int nbAR, nbAR0, nbAR1, nbAR2;
    /* BL-> left  part of 2B\(1/2)B */
    /* BR-> right part of 2B\(1/2)B */
    compBox_t componentBox, B;
    compBox_init(componentBox);
    compBox_init(B);
    connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
    realRat_t nrad;
    realRat_init(nrad);
    /* radii (i.e. halfwidth) of left and right part : (3/8)*w(B) */
    realRat_mul_si(nrad, compBox_bwidthref( componentBox ), 3);
    realRat_div_ui(nrad, nrad, 8);
    compBox_set(B, componentBox);
    realRat_mul_si(compBox_bwidthref( B ), nrad, 2);
    compBox_init_annulii(B);
    /* copy list of annulii */
    compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( ccur ) );
    while ( itb != compBox_list_end() ) {
        compBox_ptr b = compBox_list_elmt(itb);
        compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
        while ( ita != compAnn_list_end() ) {
            compAnn_list_insert_sorted_unique(compBox_annuliref(B,0), compAnn_list_elmt(ita));
             ita = compAnn_list_next(ita);
        }
        itb = compBox_list_next(itb);
    }
    /* compute the left part of 2B\(1/2)B */
    realRat_sub( compRat_realref( compBox_centerref(B) ), compRat_realref( compBox_centerref(componentBox) ), compBox_bwidthref( componentBox ) );
    realRat_add( compRat_realref( compBox_centerref(B) ), compRat_realref( compBox_centerref(B) ), nrad );
    risolate_compBox_getNumbersForCounting( &nbAL, &nbAL0, &nbAL1, &nbAL2, B, prec );
    /* compute the right part of 2B\B */
    realRat_add( compRat_realref( compBox_centerref(B) ), compRat_realref( compBox_centerref(componentBox) ), compBox_bwidthref( componentBox ) );
    realRat_sub( compRat_realref( compBox_centerref(B) ), compRat_realref( compBox_centerref(B) ), nrad );
//     realRat_mul_si(compBox_bwidthref( B ), nrad, 2);
    risolate_compBox_getNumbersForCounting( &nbAR, &nbAR0, &nbAR1, &nbAR2, B, prec );
    
    int skipTstar = ( (nbA1 + nbAL2 + nbAR2) >1 );
    
    compBox_clear_annulii(B);
    compBox_clear(B);
    compBox_clear(componentBox);
    realRat_clear(nrad);
    
    if (skipTstar)
        return -1;
    else 
        return -2;
    
}

int  risolate_rootRadii_connCmp_same_annulii( connCmp_ptr ccur1, connCmp_ptr ccur2 ) {
    
    compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( ccur1 ) );
    if (compAnn_list_get_size( compBox_annuliref( compBox_list_elmt( itb), 0 ) ) > 1)
        return 0;
    
    /* get first annulus of first box of first component */
    compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( compBox_list_elmt( itb), 0 ) );
    compAnn_ptr a = compAnn_list_elmt( ita );
    
    itb = compBox_list_next(itb);
    while ( itb!=compBox_list_end() ) {
        if (compAnn_list_get_size( compBox_annuliref( compBox_list_elmt( itb), 0 ) ) > 1)
            return 0;
        ita = compAnn_list_begin( compBox_annuliref( compBox_list_elmt( itb), 0 ) );
        if ( !(a==compAnn_list_elmt( ita )) )
            return 0;
        itb = compBox_list_next(itb);
    }
    /* here ccur1 intersects only one annulus */
    itb = compBox_list_begin( connCmp_boxesref( ccur2 ) );
    while ( itb!=compBox_list_end() ) {
        if (compAnn_list_get_size( compBox_annuliref( compBox_list_elmt( itb), 0 ) ) > 1)
            return 0;
        ita = compAnn_list_begin( compBox_annuliref( compBox_list_elmt( itb), 0 ) );
        if ( !(a==compAnn_list_elmt( ita )) )
            return 0;
        itb = compBox_list_next(itb);
    }
    
    return 1;
    
}

int  risolate_rootRadii_connCmp_is_separated( connCmp_ptr ccur, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, metadatas_t meta ){
    
    compBox_t componentBox;
    compDsk_t fourCCDisk;
    realRat_t four;
    compBox_init(componentBox);
    compDsk_init(fourCCDisk);
    realRat_init(four);
    
    realRat_set_si(four, 4, 1);
    connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
    risolate_compBox_get_containing_dsk(fourCCDisk, componentBox);
    compDsk_inflate_realRat(fourCCDisk, fourCCDisk, four);
    
    int res = 1;
    connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
    while ( res && (it!=connCmp_list_end()) ) {
        int restemp = connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , fourCCDisk);
        if ( restemp ) {
            int sameAnnulii = risolate_rootRadii_connCmp_same_annulii( ccur, connCmp_list_elmt(it) );
            if (sameAnnulii)
                restemp = 0;
        }
        res = res && (! restemp);
        it = connCmp_list_next(it);
    }
    it = connCmp_list_begin(discardedCcs);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , fourCCDisk));
        it = connCmp_list_next(it);
    }
    
    compBox_clear(componentBox);
    compDsk_clear(fourCCDisk);
    realRat_clear(four);
    return res;
}

void risolate_main_loop_rootRadii( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta){
    
    int separationFlag;
    int widthFlag;
    int compactFlag;
    int sepBoundFlag;
    
    slong prec, depth;
    tstar_res resTstar;
    newton_res resNewton;
    
    compBox_t componentBox;
    compDsk_t ccDisk, fourCCDisk;
    realRat_t three, four, threeWidth, two;
    compRat_t initPoint;
    connCmp_list_t ltemp;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    compDsk_init(fourCCDisk);
    realRat_init(three);
    realRat_init(two);
    realRat_init(four);
    realRat_init(threeWidth);
    compRat_init(initPoint);
    connCmp_list_init(ltemp);
    
    connCmp_ptr ccur;
    
    clock_t start=clock();
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    realRat_set_si(two, 2, 1);
    
    while (!connCmp_list_is_empty(qMainLoop)) {
        
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }

        resNewton.nflag = 0;
        ccur = connCmp_list_pop(qMainLoop);
        
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//             connCmp_print(ccur);
//             printf("\n");
//         }
        
        connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
        risolate_compBox_get_containing_dsk(ccDisk, componentBox);
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
        /* for test */
//         compDsk_inflate_realRat(twoCCDisk, ccDisk, two);
        /* end for test*/
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        /* do not need natural clusters, but waiting for separation before pellet's test speed'up the process */
        separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
        if (separationFlag)
            compDsk_inflate_realRat(ccDisk, ccDisk, two);
        else
            separationFlag = risolate_rootRadii_connCmp_is_separated( ccur, qMainLoop, discardedCcs, meta );
//         separationFlag = ccluster_compDsk_is_separated(ccDisk, qMainLoop, discardedCcs);
//         separationFlag = 1;
      
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        widthFlag      = widthFlag || (fmpz_is_zero(realRat_denref(eps)));
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        sepBoundFlag   = (realRat_cmp( compBox_bwidthref(componentBox), metadatas_getSepBound(meta))<=0);
        
/* return values: -2: no conclusion => run Pellet test */
/*                -1: 2B\(1/2)B contains at least one real root => skip Pellet test*/
/*                 0: ??? */
/*               >=1: B contains that numb. of sol => skip Pellet test */ 
        
        int rootRadiiCountingTestRes = risolate_rootRadii_countingTest( ccur, CCLUSTER_DEFAULT_PREC, cache, meta );
        int skipTstarFlag = (rootRadiiCountingTestRes>=-1);
        if (rootRadiiCountingTestRes >=0) {
            connCmp_nSolsref(ccur) = rootRadiiCountingTestRes;
        }
        
        //         
        if (metadatas_getVerbo(meta)>3) {
            printf("#---depth: %d\n", (int) depth);
            compApp_t centerApp;
            realApp_t widthApp;
            compApp_init(centerApp);
            realApp_init(widthApp);
            compApp_set_compRat(centerApp, compBox_centerref(componentBox), CCLUSTER_DEFAULT_PREC );
            realApp_set_realRat(widthApp,  compBox_bwidthref(componentBox), CCLUSTER_DEFAULT_PREC );
            printf("------component Box:       center: ");
            compApp_printd(centerApp, 10); printf(" width: ");
            realApp_printd(widthApp, 10); printf("\n");
            compApp_clear(centerApp);
            realApp_clear(widthApp);
            printf("#------nb of boxes:         %d\n", compBox_list_get_size(connCmp_boxesref(ccur)));
            printf("#------nb of roots:         %d\n", connCmp_nSolsref(ccur));
            printf("#------last newton success: %d\n", connCmp_newSuref(ccur));
            printf("------newton speed       : "); fmpz_print(connCmp_nwSpdref(ccur)); printf("\n");
            printf("#------last prec used     : %ld\n", prec);
            printf("#------nb of CC in m queue: %d\n", connCmp_list_get_size(qMainLoop));
            printf("#------nb of CC in r queue: %d\n", connCmp_list_get_size(qResults));
            printf("#------separation Flag:     %d\n", separationFlag);
            printf("#------widthFlag:           %d, (width: )", widthFlag); realRat_print(eps); printf(")\n");
            printf("#------compactFlag:         %d\n", compactFlag);
            printf("#------sepBoundFlag:        %d\n", sepBoundFlag);
            printf("#------rootRadiiCountingTestRes:    %d\n", rootRadiiCountingTestRes);
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)&&(skipTstarFlag==0)) {
            
//             if ( (connCmp_nSolsref(ccur)!=1) && (determined_by_rootRadii == 0) ) {
            if (connCmp_nSolsref(ccur)!=1) {
                
                resTstar = tstar_real_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
//                 resTstar = tstar_real_interface( cache, twoCCDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                 if (metadatas_getVerbo(meta)>3)
//                     printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
//                 ???
                prec = resTstar.appPrec;
                
                if (metadatas_getVerbo(meta)>3) {
                    printf("#--- Tstar test, res: %d, prec: %ld\n", resTstar.nbOfSol, prec);
                }
            }
//             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
        }
        
        /* special case where zero is a root with mult>1 */
        /* and current cc is separated and contains zero */
        if ( (separationFlag) && 
             (connCmp_nSolsref(ccur) > 1 ) &&
             (realRat_sgn(connCmp_infReref(ccur)) < 0) &&
             (realRat_sgn(connCmp_supReref(ccur)) > 0) ) {
            if ( connCmp_nSolsref(ccur) == cacheApp_getMultOfZero( cache ) ) {
//                 printf("ici\n");
                sepBoundFlag = 1;
            }
            
        }
             
        
        if ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
             ( (!widthFlag) 
               || ( (!sepBoundFlag) &&
                                    ( (connCmp_nSols(ccur)>1) || (connCmp_nSols(ccur) == cacheApp_getDegree(cache)) )
                  )
             ) ) {
            
            if (metadatas_getVerbo(meta)>3) {
                    printf("#--- Newton iteration\n");
            }
                
            if (metadatas_haveToCount(meta)){
                start = clock();
            }
        
            if (connCmp_nSols(ccur)==1) 
                compRat_set(initPoint, compBox_centerref(componentBox));
            else
                connCmp_risolate_find_point_outside_connCmp( initPoint, ccur, metadatas_initBref(meta) );
        
            connCmp_ptr nCC;
            nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
            connCmp_init(nCC);
            resNewton = newton_risolate_newton_connCmp( nCC, ccur, cache, initPoint, prec, meta);

//             printf("+++depth: %ld, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
            if (resNewton.nflag) {
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = nCC;
                connCmp_increase_nwSpd(ccur);
                connCmp_newSuref(ccur) = 1;
                connCmp_appPrref(nCC) = resNewton.appPrec;
                
                connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
//                 compBox_get_containing_dsk(ccDisk, componentBox);
                risolate_compBox_get_containing_dsk( ccDisk, componentBox);
                
                if (metadatas_useDeflation(meta)) {
                    if ( (connCmp_nSolsref(ccur) > 1 ) 
//                         && (connCmp_nSolsref(ccur) < (cacheApp_getDegree(cache)/10)) 
                       ) {
//                         if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)==0) {
                        if (connCmp_isDefref(ccur)==0) {
                            if (realRat_cmp_ui(compDsk_radiusref(ccDisk), 1 ) < 0) {
//                                 if (metadatas_getVerbo(meta)>=3){
//                                     printf("\n\n\n ------Success of Newton Iteration for this Component with a cluster of %d roots------\n", connCmp_nSolsref(ccur) );
//                                 }
                                deflate_connCmp_init(ccur);
                                deflate_set( ccur, cache, ccDisk, connCmp_nSolsref(ccur), connCmp_appPrref(ccur), meta );
                                
                            }
                        }
                    }
                }
                
                connCmp_increase_nwSpd(ccur);
    
            }
            else {
                connCmp_newSuref(ccur) = 2;
                connCmp_clear(nCC);
                ccluster_free(nCC);
            }
            if (metadatas_haveToCount(meta)){
                metadatas_add_Newton   ( meta, depth, resNewton.nflag, (double) (clock() - start) );
            }
        }
        
        if ( (connCmp_nSols(ccur)==1) && widthFlag && compactFlag ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
            if (metadatas_getVerbo(meta)>3) {
                printf("#------validated with %d roots\n", connCmp_nSols(ccur));
            }
        } else
        if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag && sepBoundFlag ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
            if (metadatas_getVerbo(meta)>3) {
                printf("#------validated with %d roots\n", connCmp_nSols(ccur));
            }
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
            
            if (metadatas_getVerbo(meta)>3) {
                    printf("#--- insert in queue\n");
            }
            
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            
            if (metadatas_getVerbo(meta)>3) {
                    printf("#--- insert in queue and decrease Newton speed\n");
            }
            
            connCmp_decrease_nwSpd(ccur);
//             if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0)
//                 connCmp_decrease_nwSpd(ccur);
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        else {
//             if (connCmp_nSols(ccur)==0) 
//                 printf("la, %d\n", compBox_list_get_size(connCmp_boxesref(ccur)));
// #ifdef CCLUSTER_HAVE_PTHREAD
//             ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ccur);
//             ccluster_free(ccur);
// #else
//             risolate_bisect_connCmp_rootRadii( ltemp, ccur, discardedCcs, bDiscarded, cache, meta,1);
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ccur);
//             ccluster_free(ccur);
            
            slong ratio = 0;
            if (compBox_list_get_size(connCmp_boxesref(ccur))==1)
                ratio = risolate_connCmp_getZoomRatio_1b( ccur, CCLUSTER_DEFAULT_PREC );
            if (compBox_list_get_size(connCmp_boxesref(ccur))==2)
                ratio = risolate_connCmp_getZoomRatio_2b( ccur, CCLUSTER_DEFAULT_PREC );
            if (metadatas_getVerbo(meta)>3) {
                    printf("#--- bisect; ratio: %ld\n", ratio);
            }
            if ( (ratio > 1)||(ratio<-1) )
                risolate_bisect_connCmp_with_ratio( ltemp, ccur, ratio);
            else 
                risolate_bisect_connCmp_rootRadii( ltemp, ccur, discardedCcs, bDiscarded, cache, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
// #endif
        }
    }
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    realRat_clear(two);
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(threeWidth);
    compRat_clear(initPoint);
    connCmp_list_clear(ltemp);
}

int risolate_nbGIt_rootRadii( slong degree, const realRat_t delta ){
    
    realRat_t oneplusdelta, oneplusdeltainv;
    realRat_init(oneplusdelta);
    realRat_init(oneplusdeltainv);
    
    realRat_add_si(oneplusdelta, delta, 1);
    realRat_inv( oneplusdeltainv, oneplusdelta );
    
    double log2_1pdelta = fmpz_dlog( realRat_numref(oneplusdelta) ) - fmpz_dlog( realRat_denref(oneplusdelta) );
    log2_1pdelta = log2_1pdelta / log(2);
    int N = (int) ceil( log2( log2(2*degree)/log2_1pdelta ) );
    return N;
}

void risolate_algo_global_rootRadii  ( connCmp_list_t qResults, 
                                       compBox_list_t bDiscarded,
                                       compAnn_list_t annulii,
                                       const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    clock_t start = clock();
    clock_t start2 = clock();
    
    slong degree = cacheApp_getDegree( cache );
    slong bitsize = cacheApp_getBitsize (cache);
    
    realRat_t upperBound;
    realRat_init(upperBound);
    
    /* set relative precision for root radii */
    metadatas_setRelPr_si( meta, 1, degree*degree);
    /* compute number of Graeffe Iterations for root radii */
    int N = risolate_nbGIt_rootRadii( degree, metadatas_getRelPr(meta) );
    metadatas_setNbGIt ( meta, N);
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#degree  of input polynomial: %ld\n", degree);
        printf("#bitsize of input polynomial: %ld\n", bitsize);
        printf("#number of Graeffe iterations for root radii: %d\n", N);
    }
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    /* heuristic to predict the precision: at least the degree */
    while (prec<degree/2)
        prec = 2*prec;
//     slong precpred = prec;
    /* */
    start2 = clock();
    prec = realIntRootRadii_rootRadii( annulii, 0, cache, prec, meta );
    
    /* use annulii to get a sharp upper bound on the roots */
//     slong upperBound = realApp_ceil_si(compAnn_radSupref(compAnn_list_last(annulii)), prec);
    realApp_ceil_realRat(upperBound, compAnn_radSupref(compAnn_list_last(annulii)), prec);
    
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing 1st RootRadii       : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed                      : %ld\n", prec);
//         printf("#precision predicted                   : %ld\n", precpred);
//         if (metadatas_getVerbo(meta)>=2) {
//             printf("#Annulii: ");
//             compAnn_list_printd(annulii, 10);
//             printf("\n\n");
//         }
// //         printf("#upperBound on the norm of the roots:");
// //         realRat_print(upperBound);
// //         printf("\n");
//     }
    
    start2 = clock();
    realIntRootRadii_connectedComponents( annulii, prec );
    
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#Annulii: ");
//         compAnn_list_printd(annulii, 10);
//         printf("\n\n");
//     }
    
//     /* approximated version */
//     start2 = clock();
//     slong prec2 = CCLUSTER_DEFAULT_PREC;
// //     slong prec2 = 424;
//     compAnn_list_t annuliiApp;
//     compAnn_list_init(annuliiApp);
// //     realRat_set_si(delta, 2*degree, 1);
//     prec2 = realApp_rootRadii_fromZero( annuliiApp, cache, delta, prec2, meta );
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing 1st RootRadii       : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed                      : %ld\n", prec2);
// //         if (metadatas_getVerbo(meta)>=3) {
//             printf("#Annulii: ");
//             compAnn_list_printd(annuliiApp, 10);
//             printf("\n\n");
// //         }
// //         printf("#upperBound on the norm of the roots:");
// //         realRat_print(upperBound);
// //         printf("\n");
//     }
//     start2 = clock();
//     realApp_rootRadii_connectedComponents( annuliiApp, prec );
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#Annulii: ");
//         compAnn_list_printd(annuliiApp, 10);
//         printf("\n\n");
//     }
//     compAnn_list_clear(annuliiApp);
    
    start2 = clock();
    
    realIntRootRadii_containsRealRoot( annulii, cache, prec );
    
//     if (metadatas_getVerbo(meta)>=2) {
// //         printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        if (metadatas_getVerbo(meta)>=3) {
            printf("#Annulii: ");
            compAnn_list_printd(annulii, 10);
            printf("\n\n");
        }
//     }
    start2 = clock();
    if (metadatas_haveToCount(meta))
        metadatas_add_time_rootRad(meta, (double) (start2 - start) );
//     
    compBox_ptr box;
    box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(box);
    compBox_set(box, initialBox);
    compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
    /* if upperBound is zero, set it to one */
    if ( realRat_is_zero(upperBound) )
        realRat_set_si(upperBound, 1, 1);
    /* set width to 2*upperBound */
    realRat_mul_si( upperBound, upperBound, 2);
    realRat_set( compBox_bwidthref(box), upperBound );
    
    compBox_set(metadatas_initBref(meta), box);
//     
    compBox_copy_annulii(box, 0, annulii);
//     
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, box);
//     
//     connCmp_list_t qPrepLoop, qMainLoop, discardedCcs;
    connCmp_list_t qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
//     connCmp_list_init(qPrepLoop);
    connCmp_list_init(discardedCcs);
//     
    connCmp_list_push(qMainLoop, initialCC);
//     
    risolate_main_loop_rootRadii( qResults,  bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
//     
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    realRat_clear(upperBound);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

/* DEPRECATED */

// int risolate_connCmp_intersects_only_one( const connCmp_t cc, int nbList ){
//     
//     int res = 1;
//     slong indMins[GEOMETRY_NB_ANN_PER_BOX];
//     compBox_ptr bcur;
//     compBox_list_iterator itb;
//     
//     itb = compBox_list_begin( connCmp_boxesref( cc ) );
//     bcur= compBox_list_elmt(itb);
//     /* each list of annulii contains at least one annulus */
//     for (int ind = 0; ind < nbList; ind++) {
//         /* check if only one annulus */
//         if ( compAnn_list_get_size(compBox_annuliref(bcur, ind))>1 )
//             res = 0;
//         indMins[ind] = compAnn_indMinref(compAnn_list_first(compBox_annuliref(bcur, ind)));
//     }
//     itb = compBox_list_next(itb);
//     while ( (res==1) && ( itb != compBox_list_end() ) ){
//         bcur= compBox_list_elmt(itb);
//         for (int ind = 0; ind < nbList; ind++) {
//             if ( compAnn_list_get_size(compBox_annuliref(bcur, ind))>1 )
//                 res = 0;
//             if ( compAnn_indMinref(compAnn_list_first(compBox_annuliref(bcur, ind))) != indMins[ind] )
//                 res = 0;
//         }
//         itb = compBox_list_next(itb);
//     }
//     
//     return res;
// }
// 
// int risolate_compBox_nbMsols( const compBox_t b, int nbList ){
//     int res = compBox_get_nbMSol(b);
// //     printf("nb of sols in the box begin: %d\n", res);
//     int ind=0;
//     while ( (ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )==1 ) ){
//         /* the maximum number of sol is the min of the number of sols in the annulii */
//         compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b,ind ) );
//         compAnn_ptr ann = compAnn_list_elmt(ita);
//         int temp = compAnn_indMaxref( ann ) - compAnn_indMinref( ann ) +1;           
//         if (temp < res)
//             res = temp;
//         ind ++;
//     }
// //     printf("nb of sols in the box end  : %d\n", res);
//     return res;    
// }
// 
// /* assume that annulii in the first list are centered in 0  */
// /* as a consequence, if a cc intersects a unique annulus of the first list that contains one root, */
// /* this root is real, and it is in the intersection of that annulus with the real axis */
// /* in the cc */
// /* NOT USED */
// int risolate_connCmp_ContainsOneRealRoot( const connCmp_t cc ) {
//     
//     
//     if (! risolate_connCmp_intersects_only_one( cc, 1 ) )
//         return -1;
//     
//     int infReSgn = realRat_sgn( connCmp_infReref( cc ) );
//     int supReSgn = realRat_sgn( connCmp_supReref( cc ) );
//     
//     /* get the annulus */
//     compBox_list_iterator itb = compBox_list_begin(connCmp_boxesref(cc));
//     compBox_ptr btemp = compBox_list_elmt(itb);
//     compAnn_ptr atemp = compAnn_list_first(compBox_annuliref(btemp, 0));
//     
//     if ( (supReSgn < 0)&&( compAnn_rrInNeref( atemp ) == 1 ) ) {
//         return 1;
//     }
//     if ( (infReSgn > 0)&&( compAnn_rrInPoref( atemp ) == 1 ) ) {
//         return 1;
//     }
//     if ( (supReSgn > 0)&&(infReSgn < 0)
//          && realApp_is_zero(compAnn_radInfref(atemp)) && realApp_is_zero(compAnn_radSupref(atemp)) ) {
//         /* contains 0 with multiplicity given in compAnn_indMaxref(atemp) - compAnn_indMinref(atemp) +1 */
//         return compAnn_indMaxref(atemp) - compAnn_indMinref(atemp) + 1;
//     }
//     
//     return -1;
// }
// 
// void connCmp_mergeAnnulii( compAnn_list_t dest, int ind, const connCmp_t cc ){
//     
//     compAnn_list_iterator itadest, itab;
//     compBox_list_iterator itb = compBox_list_begin(connCmp_boxesref(cc));
//     
// //     printf("---Merge: \n");
//     
//     while( itb != compBox_list_end() ) {
//         
// //         printf("nbox; list: \n");
// //         compAnn_list_printd(compBox_annuliref(compBox_list_elmt(itb), ind), 10);
// //         printf("\n\n");
//     
//         itab = compAnn_list_begin( compBox_annuliref(compBox_list_elmt(itb), ind) );
//         while (itab!= compAnn_list_end() ){
//             
//             int isIn = 0;
//             itadest = compAnn_list_begin(dest);
//             while ( (isIn==0) && (itadest!= compAnn_list_end()) ){
//                 if ( compAnn_list_elmt(itadest) == compAnn_list_elmt(itab) )
//                     isIn = 1;
//                 itadest = compAnn_list_next(itadest);
//             }
//             if (isIn==0)
//                 compAnn_list_insert_sorted(dest, compAnn_list_elmt(itab));
//             itab = compAnn_list_next(itab);
//         }
//         
//         itb = compBox_list_next(itb);
//     }
// //     printf("dest at the end: \n");
// //     compAnn_list_printd(dest, 10);
// //     printf("\n\n");
// }
// 
// int risolate_compBox_intersects_only_one( const compBox_t b, int nbList ){
//     int ind = 0;
//     while ( ( ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )==1 ) )
//         ind ++;
//     return (ind==nbList);
// }
// 
// int risolate_compBox_intersects_atLest_one( const compBox_t b, int nbList ){
//     int ind = 0;
//     while ( ( ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )>0 ) )
//         ind ++;
//     return (ind==nbList);
// }
// 
// int risolate_compBox_NbNonZero( const compBox_t b ){
//     int res = 0;
//     int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
//     compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
//     while ( ita != compAnn_list_end() ) {
//         compAnn_ptr a = compAnn_list_elmt(ita);
//         if ( ( sgn < 0 ) && (compAnn_rrInNeref( a )!=0) )
//             res ++;
//         if ( ( sgn > 0 ) && (compAnn_rrInPoref( a )!=0) )
//             res ++;
//         ita = compAnn_list_next(ita);
//     }
//     return res;
// }
// 
// int risolate_compBox_NbDetermined( const compBox_t b ){
//     int res = 0;
//     int sgn = realRat_sgn( compRat_realref( compBox_centerref( b ) ) );
//     compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, 0 ) );
//     while ( ita != compAnn_list_end() ) {
//         compAnn_ptr a = compAnn_list_elmt(ita);
//         if ( ( sgn < 0 ) && (compAnn_rrInNeref( a ) > -1) )
//             res ++;
//         if ( ( sgn > 0 ) && (compAnn_rrInPoref( a ) > -1) )
//             res ++;
//         ita = compAnn_list_next(ita);
//     }
//     return res;
// }
// 
// void risolate_discard_compBox_list_prepLoop_rootRadii( compBox_list_t boxes, 
//                                                        compBox_list_t bDiscarded,
//                                                        cacheApp_t cache, 
//                                                        slong prec, 
//                                                        metadatas_t meta){
//     
//     slong depth;
//     
//     compBox_list_t ltemp;
//     compDsk_t bdisk;
//     compBox_list_init(ltemp);
//     compDsk_init(bdisk);
//     
//     compBox_ptr btemp;
//     
//     while (!compBox_list_is_empty(boxes)){
//         
//         btemp = compBox_list_pop(boxes);
//         
//         risolate_compBox_get_containing_dsk(bdisk, btemp);
//         depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
//         metadatas_add_explored( meta, depth);
//             
//         /* check if btemp intersects at least an annulus of each list */
//         if ( risolate_compBox_intersects_atLest_one( btemp, 1 ) == 0 ){
//             if (metadatas_haveToCount(meta)){
//                 metadatas_add_discarded( meta, depth);
//             }
// //             if (metadatas_getDrSub(meta)==0){
//                 compBox_clear(btemp);
//                 ccluster_free(btemp);
// //             } else {
// //                 compBox_list_push(bDiscarded, btemp);
// //             }
//             continue;
//         }
//         
//         btemp->nbMSol = risolate_compBox_nbMsols( btemp, 1 );
//         compBox_list_push(ltemp, btemp);
//     }
//         
//  
//     compBox_list_swap(boxes, ltemp);
//     compBox_list_clear(ltemp);
//     compDsk_clear(bdisk);
//     
// }
// 
// void risolate_bisect_connCmp_prepLoop_rootRadii( connCmp_list_t dest, 
//                                                  connCmp_t cc, 
//                                                  connCmp_list_t discardedCcs,
//                                                  compBox_list_t bDiscarded, 
//                                                  cacheApp_t cache, 
//                                                  metadatas_t meta, 
//                                                  slong nbThreads){
//     
//     slong prec = connCmp_appPr(cc);
//     compBox_list_t subBoxes;
//     connCmp_list_t ltemp;
//     compBox_list_init(subBoxes);
//     connCmp_list_init(ltemp);
//     
//     compBox_ptr btemp;
//     connCmp_ptr ctemp;
//     
//     while (!connCmp_is_empty(cc)) {
//         btemp = connCmp_pop(cc);
//         subdBox_risolate_bisect( subBoxes, btemp );
//         compBox_clear(btemp);
//         ccluster_free(btemp);
//     }
//     
//     risolate_discard_compBox_list_prepLoop_rootRadii(subBoxes, bDiscarded, cache, prec, meta);
//     
//     while (!compBox_list_is_empty(subBoxes)) {
//         btemp = compBox_list_pop(subBoxes);
//         connCmp_union_compBox( ltemp, btemp);
//     }
//     
//     while (!connCmp_list_is_empty(ltemp)){
//         ctemp = connCmp_list_pop(ltemp);
//         connCmp_list_push(dest, ctemp);
//     }
//     
//     compBox_list_clear(subBoxes);
//     connCmp_list_clear(ltemp);
// }
// 
// 
// 
// slong risolate_exclusion_rootRadii( connCmp_list_t qCover,
//                                    cacheApp_t cache, 
//                                    metadatas_t meta){
//     
//     connCmp_list_t ltemp;
//     connCmp_ptr ccur;
//     compBox_ptr box;
//     
//     compDsk_t ccDisk;
//     compDsk_init(ccDisk);
//     
//     connCmp_list_init(ltemp);
//     
//     compBox_t componentBox;
//     compBox_init(componentBox);
//     
//     slong depth = 0;
//     tstar_res res;
//     res.appPrec = CCLUSTER_DEFAULT_PREC;
//     
//     while (!connCmp_list_is_empty(qCover)) {
//           ccur = connCmp_list_pop(qCover);
//           
//           if (connCmp_nSolsref(ccur) == -1) {
//               
//             connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));  
//               
//             box  = compBox_list_first( connCmp_boxesref(ccur) );
//             /*get containing disk */
//             risolate_compBox_get_containing_dsk(ccDisk, componentBox);
//             /* do an exclusion test */
//             res = tstar_real_interface( cache, ccDisk, compBox_get_nbMSol(box), 1, 0, res.appPrec , depth, meta);
// //             res.appPrec = res.appPrec/2;
//             
//             if (res.nbOfSol==0) { /* clear ccur */
//                 if (metadatas_haveToCount(meta)){
//                         metadatas_add_discarded( meta, depth);
//                 }
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//             } else {
//                 if (res.nbOfSol>0) {
//                     compBox_nbMSolref(box) = res.nbOfSol;
//                     connCmp_nSolsref(ccur) = res.nbOfSol;
//                 }
//                 connCmp_appPrref(ccur) = res.appPrec;
//                 connCmp_list_insert_sorted(ltemp, ccur);
//             }
//           } else {
//               connCmp_list_insert_sorted(ltemp, ccur);
//           }
//     }
//     
//     
//     connCmp_list_swap(ltemp, qCover);
//     
//     connCmp_list_clear(ltemp);
//     compDsk_clear(ccDisk);
//     
//     compBox_clear(componentBox);
//     
//     return res.appPrec;
// }
// 
// void risolate_prep_loop_rootRadii( compBox_list_t bDiscarded, 
//                                    connCmp_list_t qResult, 
//                                    connCmp_list_t qPrepLoop, 
//                                    connCmp_list_t discardedCcs, 
//                                    cacheApp_t cache, 
//                                    metadatas_t meta) {
//     
//     connCmp_ptr ctemp;
//     connCmp_list_t ltemp;
//     connCmp_list_init(ltemp);
//     
//     int widthFlag;
//     int intersectOnlyOneFlag;
//     int containsOneRealRootFlag;
//     int separationFlag;
//     int compactFlag;
//     int sgn;
//     
//     compBox_t componentBox;
//     compDsk_t ccDisk, fourCCDisk;
//     realRat_t four;
//     realApp_t widthCtemp;
//     realApp_t widthAnn;
//     compBox_init(componentBox);
//     compDsk_init(ccDisk);
//     compDsk_init(fourCCDisk);
//     realRat_init(four);
//     realApp_init(widthCtemp);
//     realApp_init(widthAnn);
//     
//     compBox_ptr bcur, nbox;
//     compAnn_ptr acur;
//     
//     realRat_set_si(four, 4, 1);
//     
//     
//     while (!connCmp_list_is_empty(qPrepLoop)) {
//         ctemp = connCmp_list_pop(qPrepLoop);
//         
//         intersectOnlyOneFlag = risolate_connCmp_intersects_only_one( ctemp, 1 );
//         containsOneRealRootFlag = -1;
//         widthFlag = 0;
//         separationFlag = 0;
//         compactFlag    = compBox_list_get_size(connCmp_boxesref(ctemp)) <= 3;
//         
//         
//         if (intersectOnlyOneFlag == 1) {
//             bcur = compBox_list_first(connCmp_boxesref(ctemp));
//             acur = compAnn_list_first(compBox_annuliref(bcur, 0));
//             
//             /* check if ctemp contains the Zero annulus */
//             if ( realApp_is_zero(compAnn_radInfref(acur)) && realApp_is_zero(compAnn_radSupref(acur)) ) {
// //                 containsOneRealRootFlag = 1;
//                 /* in this case, let containsOneRealRootFlag = multiplicity of 0 */
//                 containsOneRealRootFlag = compAnn_indMaxref(acur) - compAnn_indMinref(acur) +1;
//             }
//             
//             sgn = realRat_sgn(connCmp_supReref(ctemp));
//             if ( (sgn<0) && (compAnn_rrInNeref(acur) == 1) )
//                 containsOneRealRootFlag = 1;
//             if ( (sgn<0) && (compAnn_rrInNeref(acur) == 0) ) {
//                 connCmp_clear(ctemp);
//                 ccluster_free(ctemp);
//                 continue;
//             }
//             
//             sgn = realRat_sgn(connCmp_infReref(ctemp));
//             if ( (sgn>0) && (compAnn_rrInPoref(acur) == 1) )
//                 containsOneRealRootFlag = 1;
//             if ( (sgn>0) && (compAnn_rrInPoref(acur) == 0) ) {
//                 connCmp_clear(ctemp);
//                 ccluster_free(ctemp);
//                 continue;
//             }
//             
//             realApp_set_realRat( widthCtemp, compBox_bwidthref(bcur), CCLUSTER_DEFAULT_PREC );
//             realApp_set(widthAnn, compAnn_radSupref(acur));
//             realApp_sub(widthAnn, widthAnn, compAnn_radInfref(acur), CCLUSTER_DEFAULT_PREC );
//             widthFlag = realApp_lt(widthCtemp, widthAnn);
//             if ((widthFlag != 1)&&( realApp_ge(widthCtemp, widthAnn) != 1 ))
//                 widthFlag = 1;
//         }
//         
// //         if ((intersectOnlyOneFlag == 1) && (containsOneRealRootFlag==1)){
//         if ((intersectOnlyOneFlag == 1) && (containsOneRealRootFlag >= 1)){
//             connCmp_risolate_componentBox(componentBox, ctemp, metadatas_initBref(meta));
//             compBox_get_containing_dsk(ccDisk, componentBox);
//             compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
// //             separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qPrepLoop, discardedCcs);
//             separationFlag = ccluster_compDsk_is_separated(ccDisk, qPrepLoop, discardedCcs);
//         }
//         
//         if (intersectOnlyOneFlag == 0) {
//             risolate_bisect_connCmp_prepLoop_rootRadii( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ctemp);
//             ccluster_free(ctemp);
//         } 
//         else {
//             if ( ( widthFlag == 0 )
//                &&( (containsOneRealRootFlag==-1)||(separationFlag==0) )
//             ){
//                 risolate_bisect_connCmp_prepLoop_rootRadii( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
//                 while (!connCmp_list_is_empty(ltemp))
//                     connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
//                 connCmp_clear(ctemp);
//                 ccluster_free(ctemp);
//             }
//             else {
//                 
//                 if (compactFlag==0){ /* replace the CC by a cc containing only one box */
// //                     
// //                     printf("# Non compact cc: nb of boxes: %d\n", compBox_list_get_size(connCmp_boxesref(ctemp)) );
//                     bcur = compBox_list_first(connCmp_boxesref(ctemp));
// //                     
//                     nbox = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//                     compBox_init(nbox);
//                     connCmp_risolate_componentBox( nbox, ctemp, metadatas_initBref(meta));
//                     compBox_nbMSolref( nbox ) = risolate_compBox_nbMsols( bcur, 1 );
//                     /* copy annulii 0 */
//                     compBox_copy_annulii(nbox, 0, compBox_annuliref(bcur, 0));
//                     /*delete ctemp*/
//                     connCmp_clear(ctemp);
//                     ccluster_free(ctemp);
//                     /* create the new connected component */
//                     ctemp = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//                     connCmp_init_compBox(ctemp, nbox);
//                     
//                 }
//                 
//                 connCmp_nSolsref(ctemp) = containsOneRealRootFlag;
//                 connCmp_list_insert_sorted(qResult, ctemp);
//             }
//         }
//             
//     }
//     connCmp_list_clear(ltemp);
//     
//     compBox_clear(componentBox);
//     compDsk_clear(ccDisk);
//     compDsk_clear(fourCCDisk);
//     realRat_clear(four);
//     realApp_clear(widthCtemp);
//     realApp_clear(widthAnn);
//     
// }
// 
// void risolate_algo_global_rootRadii  ( connCmp_list_t qResults, 
//                                        compBox_list_t bDiscarded,
//                                        compAnn_list_t annulii,
//                                        const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
//     
//     clock_t start = clock();
//     clock_t start2 = clock();
//     
//     realRat_t delta;
//     realRat_init(delta);
//     
//     slong degree = cacheApp_getDegree( cache );
// //     realRat_set_si(delta, 1, 1);
//     realRat_set_si(delta, 1, degree);
// //     realRat_set_si(delta, 1, degree*degree);
// //     realRat_set_si(delta, 1, degree*degree*degree);
//     
//     slong prec = CCLUSTER_DEFAULT_PREC;
//     /* heuristic to predict the precision: at least the degree */
//     while (prec<degree)
//         prec = 2*prec;
//     prec = realIntRootRadii_rootRadii( annulii, 0, cache, delta, prec, meta );
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing 1st RootRadii       : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed                      : %ld\n", prec);
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#Annulii: ");
//             compAnn_list_printd(annulii, 10);
//             printf("\n\n");
//         }
//     }
//     
//     /* use annulii to get a sharp upper bound on the roots */
//     slong upperBound = realApp_ceil_si(compAnn_radSupref(compAnn_list_last(annulii)), prec);
//     start2 = clock();
//     
//     realIntRootRadii_connectedComponents( annulii, prec );
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//     }
//     start2 = clock();
//     
//     realIntRootRadii_containsRealRoot( annulii, cache, prec );
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#Annulii: ");
//             compAnn_list_printd(annulii, 10);
//             printf("\n\n");
//         }
//     }
//     start2 = clock();
//     if (metadatas_haveToCount(meta))
//         metadatas_add_time_rootRad(meta, (double) (start2 - start) );
//     
//     compBox_ptr box;
//     box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//     compBox_init(box);
//     compBox_set(box, initialBox);
//     compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
//     /* set width to 2*upperBound */
//     realRat_set_si( compBox_bwidthref(box), 2*upperBound, 1 );
//     
//     compBox_copy_annulii(box, 0, annulii);
//     
//     connCmp_ptr initialCC;
//     initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//     connCmp_init_compBox(initialCC, box);
//     
//     connCmp_list_t qPrepLoop, qMainLoop, discardedCcs;
//     connCmp_list_init(qMainLoop);
//     connCmp_list_init(qPrepLoop);
//     connCmp_list_init(discardedCcs);
//     
//     connCmp_list_push(qPrepLoop, initialCC);
//     risolate_prep_loop_rootRadii(bDiscarded, qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in prep loop: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#Nb of CC: %d\n", connCmp_list_get_size(qMainLoop));
//         if (metadatas_getVerbo(meta)>=3) {
//         connCmp_list_iterator itc = connCmp_list_begin(qMainLoop);
//         while ( itc!= connCmp_list_end() ){
//             printf("#--- CC with %2d sols:\n", connCmp_nSolsref( connCmp_list_elmt(itc) ) );
// // //               compBox_ptr b= compBox_list_first(connCmp_boxesref( connCmp_list_elmt(itc) ));
// // //               for (int ind = 0; ind < GEOMETRY_NB_ANN_PER_BOX; ind++){
// // //                   if (compAnn_list_get_size(compBox_annuliref(b, ind))>0){
// // //                       compAnn_ptr ann = compAnn_list_first( compBox_annuliref(b, ind) );
// // //                       compAnn_printd(ann, 10);
// // //                       printf("#\n");
// // //                   }
// // //               }
//                 itc = connCmp_list_next(itc);
//             }
//         }
//     }
//     
//     start2 = clock();
//     
//     prec = risolate_exclusion_rootRadii( qMainLoop, cache, meta);
//     
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in exclusion tests               : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed: %ld\n", prec);
//         printf("#Nb of CC: %d\n", connCmp_list_get_size(qMainLoop));
//         if (metadatas_getVerbo(meta)>=3) {
//         connCmp_list_iterator itc = connCmp_list_begin(qMainLoop);
//         while ( itc!= connCmp_list_end() ){
//             printf("#--- CC with %2d sols:\n", connCmp_nSolsref( connCmp_list_elmt(itc) ) );
// // //               compBox_ptr b= compBox_list_first(connCmp_boxesref( connCmp_list_elmt(itc) ));
// // //               for (int ind = 0; ind < GEOMETRY_NB_ANN_PER_BOX; ind++){
// // //                   if (compAnn_list_get_size(compBox_annuliref(b, ind))>0){
// // //                       compAnn_ptr ann = compAnn_list_first( compBox_annuliref(b, ind) );
// // //                       compAnn_printd(ann, 10);
// // //                       printf("#\n");
// // //                   }
// // //               }
//                 itc = connCmp_list_next(itc);
//             }
//         }
//     }
//     start2 = clock();
//     
//     risolate_main_loop( qResults,  bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
// //     connCmp_list_swap(qResults, qMainLoop);
//     
// //     realRat_clear(factor);
//     connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
//     connCmp_list_clear(discardedCcs);
//     realRat_clear(delta);
//     
// //     chronos_toc_CclusAl(metadatas_chronref(meta));
//     metadatas_add_time_CclusAl(meta, (double) (clock() - start));
// }
// 
// // OLD RISOLATE ROOTRADII
// 
// // void risolate_prep_loop_rootRadii2( connCmp_list_t qCover, 
// //                                    const compBox_t initialBox,
// //                                    cacheApp_t cache, 
// //                                    metadatas_t meta){
// //     
// //     /* qCover contains a list of cc intersecting one annulus */
// //     connCmp_ptr ccur;
// //     connCmp_list_t ltemp;
// //     connCmp_list_init(ltemp);
// //     
// //     compBox_ptr box;
// //     compAnn_ptr ann;
// //     compBox_list_iterator itb;
// //     compAnn_list_ptr annl;
// //     
// //     while (!connCmp_list_is_empty(qCover)) {
// //         ccur = connCmp_list_pop(qCover);
// //         itb  = compBox_list_begin( connCmp_boxesref( ccur ) );
// //         box  = compBox_list_elmt( itb );
// //         annl = compBox_annuliref( box );
// //         ann  = compAnn_list_first(annl);
// //         /* transform ccur in a cc with a unique box */
// //         box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
// //         compBox_init(box);
// //         connCmp_risolate_componentBox( box, ccur, initialBox);
// //         compBox_init_annuli(box);
// //         compBox_copy_annuli(box, annl);
// //         compBox_nbMSolref( box ) = compAnn_indMaxref(ann) - compAnn_indMinref(ann) + 1;
// //         /*delete Ccur*/
// //         while ( itb!=compBox_list_end() ){
// //           compBox_clear_annuli(compBox_list_elmt(itb));
// //           itb = compBox_list_next(itb);
// //         }
// //         connCmp_clear(ccur);
// //         ccluster_free(ccur);
// //         /* create the new connected component */
// //         ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
// //         connCmp_init_compBox(ccur, box);
// //         /*check if ccur contains a unique real root */
// //         if (compAnn_rrInPoref(ann)>-1)
// //           connCmp_nSolsref(ccur) = 1;
// //           
// //         /* push ccur in ltemp */
// //         connCmp_list_insert_sorted(ltemp, ccur);
// //         
// //     }
// //     
// //     connCmp_list_swap(ltemp, qCover);
// //     connCmp_list_clear(ltemp);
// //     
// // }
// 
// void risolate_prep_loop_rootRadii_old( compBox_list_t bDiscarded,
//                                    connCmp_list_t qCover, 
//                                    const compBox_t initialBox,
//                                    const compAnn_list_t annulii,
//                                    cacheApp_t cache, 
//                                    metadatas_t meta) {
//     
//     connCmp_list_t qLoop, ltemp;
//     connCmp_ptr ccur;
//     compBox_ptr box;
//     compAnn_ptr ann;
//     
//     connCmp_list_init(qLoop);
//     connCmp_list_init(ltemp);
//     
//     box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//     compBox_init(box);
//     compBox_set(box, initialBox);
//     compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
// //     compBox_init_annuli(box);
//     compBox_copy_annulii(box, 0, annulii);
//     
//     ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//     connCmp_init_compBox(ccur, box);
//     connCmp_list_push(qLoop, ccur);
//     
//     compBox_list_iterator itb;
//     compAnn_list_ptr annl;
//     int intersectUniqueAnnulus = 0;
//     
//     while (!connCmp_list_is_empty(qLoop)) {
//           ccur = connCmp_list_pop(qLoop);
//           /* check if ccur intersects a unique annulus */  
//           itb = compBox_list_begin( connCmp_boxesref( ccur ) );
//           box = compBox_list_elmt( itb ) ;                     /* at least one box in ccur */
//           intersectUniqueAnnulus = (compAnn_list_get_size( compBox_annuli0ref( box ) ) == 1);
//           ann = compAnn_list_first(compBox_annuli0ref( box ));  /* box intersects at least one annulus */
//           itb = compBox_list_next( itb );
//           while ( ( intersectUniqueAnnulus == 1 ) && 
//                     (itb != compBox_list_end()   ) ) {
//               
//               box = compBox_list_elmt( itb ) ;
//               annl = compBox_annuli0ref( box ); 
//               intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_get_size( annl ) == 1) ;
//               intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_first(annl) == ann ) ;
//               itb = compBox_list_next( itb );
//           }
//           
//           if (intersectUniqueAnnulus==1) {
//               
//               itb  = compBox_list_begin( connCmp_boxesref( ccur ) );
//               box  = compBox_list_elmt( itb );
//               annl = compBox_annuli0ref( box );
//               ann  = compAnn_list_first(annl);
//               
//               /* transform ccur in a cc with a unique box */
//               box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//               compBox_init(box);
//               connCmp_risolate_componentBox( box, ccur, initialBox);
//               compBox_copy_annulii(box, 0, annl);
//               compBox_nbMSolref( box ) = compAnn_indMaxref(ann) - compAnn_indMinref(ann) + 1;
//               /*delete Ccur*/
//               connCmp_clear(ccur);
//               ccluster_free(ccur);
//               /* create the new connected component */
//               ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//               connCmp_init_compBox(ccur, box);
//               /*check if ccur contains a unique real root */
//               if ( (compAnn_indMaxref(ann)==compAnn_indMinref(ann)) && (compAnn_rrInPoref(ann)>-1) )
//                 connCmp_nSolsref(ccur) = 1;
//                 
//               /* push ccur in qCover */
//                 connCmp_list_insert_sorted(qCover, ccur);
//               
//           } else { /* bisect ccur */
//                 
//                 realIntRootRadii_bisect_connCmp( ltemp, ccur);
//                 while (!connCmp_list_is_empty(ltemp))
//                     connCmp_list_insert_sorted(qLoop, connCmp_list_pop(ltemp));
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//                 
//           }
//     }
//     
//     connCmp_list_clear(ltemp);
//     connCmp_list_clear(qLoop);
// }
// 
// void risolate_algo_global_rootRadii_old( connCmp_list_t qResults, 
//                                      compBox_list_t bDiscarded,
//                                      compAnn_list_t annulii,
//                                      compAnn_list_t annulii1,
//                                      compAnn_list_t annulii2,
//                                      const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
//     
//     clock_t start = clock();
//     clock_t start2 = clock();
//     
// //     compAnn_list_t annulii;
//     connCmp_list_t qCover;
//     realRat_t delta;
//     
// //     compAnn_list_init(annulii);
//     connCmp_list_init(qCover);
//     realRat_init(delta);
//     
//     slong degree = cacheApp_getDegree( cache );
// //     realRat_set_si(delta, 1, degree);
//     realRat_set_si(delta, 1, degree*degree);
// //     realRat_set_si(delta, 1, degree*degree*degree);
//     
//     slong prec = CCLUSTER_DEFAULT_PREC;
//     prec = realIntRootRadii_rootRadii( annulii, 0, cache, delta, prec, meta );
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#time in computing RootRadii           : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed: %ld\n", prec);
// //     }
//     start2 = clock();
//     
//     realIntRootRadii_connectedComponents( annulii, prec );
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
// //     }
//     start2 = clock();
//     
// //     printf("Annulii: ");
// //     compAnn_list_printd(annulii, 10);
// //     printf("\n\n");
//     
//     realIntRootRadii_containsRealRoot( annulii, cache, prec );
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
// //     }
//     start2 = clock();
//     
// //     printf("Annulii: ");
// //     compAnn_list_printd(annulii, 10);
// //     printf("\n\n");
//     
//     risolate_prep_loop_rootRadii_old( bDiscarded, qCover, initialBox, annulii, cache, meta);
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#time in covering intervals with boxes : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
// //     }
//     start2 = clock();
//     
// //     connCmp_list_iterator itc;
// // //     display qCover
// //     itc = connCmp_list_begin(qCover);
// // //     if (metadatas_getVerbo(meta)>3) {
// //         printf("#Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
// // //     }
// //     while( itc!= connCmp_list_end() ){
// //         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
// //         printf("\n");
// //         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
// //         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))),0)) );
// // //            
// //         itc = connCmp_list_next( itc );
// //     }
//     printf("\n\n");
//         
//     prec = risolate_exclusion_rootRadii( qCover, cache, meta);
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#time in exclusion tests               : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         printf("#precision needed: %ld\n", prec);
// //     }
//     start2 = clock();
//     
//     /* display qCover */
// //     itc = connCmp_list_begin(qCover);
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
// //     }
// //     while( itc!= connCmp_list_end() ){
// //         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
// //         printf("\n");
// //         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
// //         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))))) ); 
// //         itc = connCmp_list_next( itc );
// //     }
// //     printf("\n\n");
//     
//     
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#total time in root radii              : %f \n", ((double) (clock() - start))/CLOCKS_PER_SEC );
// //     }
//     
//     connCmp_list_t discardedCcs;
//     connCmp_list_init(discardedCcs);
//     
//     /* main loop */
//     risolate_main_loop( qResults, bDiscarded,  qCover, discardedCcs, eps, cache, meta);
//     
//     connCmp_list_clear(qCover);
//     realRat_clear(delta);
// //     compAnn_list_clear(annulii);
//     connCmp_list_clear(discardedCcs);
//     metadatas_add_time_CclusAl(meta, (double) (clock() - start));
//     
// //     compBox_ptr box;
// //     box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
// //     compBox_init(box);
// //     compBox_set(box, initialBox);
// //     compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
//     
// //     connCmp_ptr initialCC;
// //     initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
// //     connCmp_init_compBox(initialCC, box);
// /*    
//     connCmp_list_t qMainLoop, discardedCcs;
//     connCmp_list_init(qMainLoop);
//     connCmp_list_init(discardedCcs);
//     
//     connCmp_list_push(qMainLoop, initialCC);
//     risolate_main_loop( qResults,  qMainLoop, discardedCcs, eps, cache, meta);
//     
//     
//     connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(discardedCcs);
//     
//     metadatas_add_time_CclusAl(meta, (double) (clock() - start));*/
// }
