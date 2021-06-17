/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "ccluster/ccluster.h"

/* assume b does not intersect real axis:
 * can get rid of all annulii with only one root */
void ccluster_actualize_annulii_real( compBox_t b ) {
    
    compAnn_list_t ltemp;
    compAnn_list_init(ltemp);
    
    for (int ind = 0; ind < GEOMETRY_NB_ANN_PER_BOX; ind++){
        if (compAnn_list_get_size( compBox_annuliref( b, ind ) ) >=1 ){
            compAnn_list_iterator it = compAnn_list_begin( compBox_annuliref( b, ind ) );
            if ( compAnn_centerImref(compAnn_list_elmt( it )) == 0 ) {
                while (it!=compAnn_list_end() ) {
                    if ( compAnn_indMaxref(compAnn_list_elmt( it )) > compAnn_indMinref(compAnn_list_elmt( it )) )
                        compAnn_list_push( ltemp, compAnn_list_elmt( it ) );
                    it = compAnn_list_next(it);
                }
                compAnn_list_swap(  compBox_annuliref( b, ind ), ltemp );
                compAnn_list_empty(  ltemp );
            }
        }
    }
                
    compAnn_list_clear(ltemp);

}

void ccluster_compBox_getNumbers_real( int * nbA,  /* nb of annulii intersecting the box */
                                       int * nbA0, /* nb of segments intersecting the box where the nb of real roots is 0 */
                                       int * nbA1, /* nb of segments contained in  b where the nb of real roots is >=1 */
                                       int * nbA2, /* nb of segments contained in 2b where the nb of real roots is >=1 */
                                       const compBox_t b,
                                       int ind, /* index of the group of annulii that is used: should be 0 */
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
    compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b, ind ) );
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

int ccluster_compBox_intersects_atLest_one( const compBox_t b, int nbList ){
    
    int ind = 0;
    while ( ( ind < nbList)
        &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )>0 ) )
        ind ++;
    return (ind==nbList);
}

/* Precondition:                                                                  */
/* Specification: returns false only if p is not in b                             */
int ccluster_is_compApp_in_box                     ( const compApp_t p,  const compBox_t b, slong prec  ){
    
    int res = 1;
    compApp_t dist;
    realApp_t halfwidth;
    compApp_init(dist);
    realApp_init(halfwidth);
    
    realApp_set_realRat( halfwidth, compBox_bwidthref(b), prec );
    realApp_div_ui( halfwidth, halfwidth, 2, prec );
    compApp_set_compRat( dist, compBox_centerref(b), prec );
    compApp_sub( dist, dist, p, prec);
    realApp_abs( compApp_realref( dist ), compApp_realref( dist ) );
    realApp_abs( compApp_imagref( dist ), compApp_imagref( dist ) );
    
    if ( ( realApp_gt( compApp_realref( dist ), halfwidth ) == 1 )
      || ( realApp_gt( compApp_imagref( dist ), halfwidth ) == 1 ) )
        res = 0;
    
    compApp_clear(dist);
    realApp_clear(halfwidth);
    
    return res;
}

int ccluster_is_compApp_in_compAnn (const compApp_t p, const compAnn_t ann, slong prec ){
    
    int res = 1;
    compApp_t dist;
    realApp_t mod;
    
    compApp_init(dist);
    realApp_init(mod);
    
    compApp_set( dist, p );
    realApp_add_si( compApp_realref(dist), compApp_realref(dist), -compAnn_centerReref(ann), prec );
    realApp_add_si( compApp_imagref(dist), compApp_imagref(dist), -compAnn_centerImref(ann), prec );
    compApp_abs( mod, dist, prec );
    
    if ( ( realApp_gt( mod, compAnn_radSupref( ann ) ) == 1 )
      || ( realApp_lt( mod, compAnn_radInfref( ann ) ) == 1 ) )
        res = 0;
        
    compApp_clear(dist);
    realApp_clear(mod);
    
    return res;
}

/* returns 1 => 2*btemp contains at least one root */
/*         0 => btemp contains no root */
/*        -1 => can not decide */
int ccluster_rootRadii_exclusion_test( compBox_t box, slong prec, metadatas_t meta ) {
    clock_t start = clock();
    int res = -1; 
    
    if ( compBox_is_imaginary_positive_strict(box)
              || compBox_is_imaginary_negative_strict(box) ) {
                ccluster_actualize_annulii_real( box );
//                 printf("# do not intersect real axis \n");
    } else {
        
        int nbA;  /* nb of annulii intersecting the box */
        int nbA0; /* nb of segments intersecting the box where the nb of real roots is 0 */
        int nbA1; /* nb of segments contained in  b where the nb of real roots is >=1 */
        int nbA2; /* nb of segments contained in 2b where the nb of real roots is >=1 */
        ccluster_compBox_getNumbers_real( &nbA, &nbA0, &nbA1, &nbA2, box, 0, prec);
        
        if ( (nbA==1) && ( (nbA1==1)||(nbA2==1) ) ) {
            res = 1;
        } else if ( (nbA1>=1)||(nbA2>=1) ) {
            res = 1;
        }
    }
    
    if ( (res==-1) && (ccluster_compBox_intersects_atLest_one( box, 3 ) == 0 ))
        res = 0;
    
    if (res==-1) {
        /* count nb of intersection */
        int nbInt = 0;
        compApp_t inter;
        compApp_init(inter);
        
        compAnn_list_iterator it = compAnn_list_begin( compBox_annuliref( box, 0 ) );
        while ( (it!=compAnn_list_end() ) && (nbInt==0) ) {
            compAnn_list_iterator it1 = compAnn_list_begin( compBox_annuliref( box, 1 ) );
            while ( (it1!=compAnn_list_end() ) && (nbInt==0) ) {
//                 printf("###\n");
                int intersect = compAnn_intersect_realCenter( inter, compAnn_list_elmt( it ), compAnn_list_elmt( it1 ),
                                                               CCLUSTER_DEFAULT_PREC);
//                 printf("## intersect: %d\n", intersect );
                /* check if the intersection is in the box */
                if (intersect)
                    intersect = ccluster_is_compApp_in_box (inter,  box, CCLUSTER_DEFAULT_PREC);
//                 printf("## intersection is in box: %d\n", intersect);
                
                /* check if the intersection is in at least one annulus of the third group */
                if (intersect) {
                    intersect = 0;
                    compAnn_list_iterator it2 = compAnn_list_begin( compBox_annuliref( box, 2 ) );
                    while ((it2!=compAnn_list_end() )&&(intersect==0)) {
                        intersect = ccluster_is_compApp_in_compAnn( inter, compAnn_list_elmt( it2 ), CCLUSTER_DEFAULT_PREC);
                        it2 = compAnn_list_next(it2);
                    }
                }
                /* check if the conjugate of the intersection is in at least one annulus of the fourth group */
                if (intersect) {
                    intersect = 0;
                    compApp_t interconj;
                    compApp_init(interconj);
                    compApp_set(interconj, inter);
                    realApp_neg( compApp_imagref(interconj), compApp_imagref(interconj) );
                    compAnn_list_iterator it3 = compAnn_list_begin( compBox_annuliref( box, 3 ) );
                    while ((it3!=compAnn_list_end() )&&(intersect==0)) {
                        intersect = ccluster_is_compApp_in_compAnn( interconj, compAnn_list_elmt( it3 ), CCLUSTER_DEFAULT_PREC);
                        it3 = compAnn_list_next(it3);
                    }
                    
                    compApp_clear(interconj);
                }
                
                if (intersect)
                    nbInt ++;
                it1 = compAnn_list_next(it1);
            }
            it = compAnn_list_next(it);
        }
        
        compApp_clear(inter);
        
        if (nbInt==0) {
            res = 0;
        }
        
    }
               
    if (metadatas_haveToCount(meta))
        metadatas_add_time_RRT0Tes(meta, (float) clock() - start );
    
    return res;
}

int ccluster_nbGIt_rootRadii( slong degree, const realRat_t delta ){
    
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

void ccluster_algo_global_rootRadii( connCmp_list_t qResults,
                                     compBox_list_t bDiscarded,
                                     compAnn_list_t annulii,
                                     compAnn_list_t annulii1,
                                     compAnn_list_t annulii2,
                                     const compBox_t initialBox, 
                                     const realRat_t eps, 
                                     cacheApp_t cache, 
                                     metadatas_t meta){
    clock_t start = clock();
    clock_t start2 = clock();
    
    slong degree = cacheApp_getDegree( cache );
    slong bitsize = cacheApp_getBitsize (cache);
    
    /* upper bound for the norm of the roots */
    realRat_t upperBound;
    realRat_init(upperBound);
    
    /* set relative precision for root radii */
    metadatas_setRelPr_si( meta, 1, degree*degree);
//     metadatas_setRelPr_si( meta, 1, degree*bitsize);
    /* compute number of Graeffe Iterations for root radii */
    int N = ccluster_nbGIt_rootRadii( degree, metadatas_getRelPr(meta) );
    metadatas_setNbGIt ( meta, N);
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#degree  of input polynomial: %ld\n", degree);
        printf("#bitsize of input polynomial: %ld\n", bitsize);
        printf("#number of Graeffe iterations for root radii: %d\n", N);
    }
    
    /* assume 0 is not a root of the poly */
    slong prec = CCLUSTER_DEFAULT_PREC;
    /* heuristic to predict the precision: at least the degree */
    while (prec<degree/2)
        prec = 2*prec;
    /* */
    slong prec1 = realIntRootRadii_rootRadii( annulii, 0, cache, prec, meta );
//     slong prec1 = realIntRootRadii_rootRadii_2( annulii, 0, cache, prec, meta );
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("time in first root radii: %f\n", (double) (clock() - start)/CLOCKS_PER_SEC );
    }
    
    /* derive upper bound for the norm of the roots */
//     slong centerRe1 = realApp_ceil_si(compAnn_radSupref(compAnn_list_last(annulii)), prec);
    realApp_ceil_realRat(upperBound, compAnn_radSupref(compAnn_list_last(annulii)), prec);
    prec = 2*prec1;
    slong prec2 = realIntRootRadii_rootRadii( annulii1, 1, cache, prec, meta );
    slong prec3 = realIntRootRadii_rootRadii_imagCenter( annulii2, 1, cache, prec, meta );
    
    prec3 = CCLUSTER_MAX( prec2, prec3);
    prec3 = CCLUSTER_MAX( prec1, prec3);
    
    realIntRootRadii_connectedComponents( annulii, prec );
    realIntRootRadii_connectedComponents( annulii1, prec );
    realIntRootRadii_connectedComponents( annulii2, prec );
    
    realIntRootRadii_containsRealRoot( annulii, cache, prec );
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        if (metadatas_getVerbo(meta)>=3) {
            printf("#Annulii: ");
            compAnn_list_printd(annulii, 10);
            printf("\n\n");
        }
//     }
    
    start2 = clock();
    if (metadatas_haveToCount(meta))
        metadatas_add_time_rootRad(meta, (double) (start2 - start) );
    
//     printf("#Annulii 0: ");
//     compAnn_list_printd(annulii, 10);
//     printf("\n\n");
// //     
//     printf("#Annulii 1: ");
//     compAnn_list_printd(annulii1, 10);
//     printf("\n\n");
//     
//     printf("#Annulii 2: ");
//     compAnn_list_printd(annulii2, 10);
//     printf("\n\n");
    
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
    
    compBox_copy_annulii(box, 0, annulii);
    compBox_copy_annulii(box, 1, annulii1);
    compBox_copy_annulii(box, 2, annulii2);
    /* fourth list of annulii contains annulii from i having intersection with the */
    /* complex conjugate of the box */
    compBox_copy_annulii(box, 3, annulii2);
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, box);
    
//     connCmp_list_t qPrepLoop;
//     connCmp_list_init(qPrepLoop);
    
    connCmp_list_t qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
    connCmp_list_init(discardedCcs);
    
//     connCmp_list_push(qPrepLoop, initialCC);
    
    connCmp_list_push(qMainLoop, initialCC);
    ccluster_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cache, meta);
    
    
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(discardedCcs);
    realRat_clear(upperBound);
    
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

/* DEPRECATED */
    /* test: compute all the intersections */
//     int nbInt = 0;
//     compApp_t inter, interConj;
//     compApp_init(inter);
//     compApp_init(interConj);
//     compAnn_list_iterator it = compAnn_list_begin( annulii );
//     while ( it!=compAnn_list_end() ) {
//         
//         int intersect = 0;
//         int nbInterA0 = 0;
//         /* find first annulus of annulii1 intersecting that annulus */
//         compAnn_list_iterator it1 = compAnn_list_begin( annulii1 );
//         while ( it1!=compAnn_list_end() ) {
//             
//             intersect = compAnn_intersect_realCenter( inter, compAnn_list_elmt( it ), compAnn_list_elmt( it1 ),
//                                                                        CCLUSTER_DEFAULT_PREC);
//             if (intersect) {
//                 int containsRealLine = realApp_contains_zero ( compApp_imagref( inter ) );
//                 int nbSolsInIt = compAnn_indMaxref(compAnn_list_elmt( it )) - compAnn_indMinref(compAnn_list_elmt( it )) + 1;
//                 int nbSolsInIt1 = compAnn_indMaxref(compAnn_list_elmt( it1 )) - compAnn_indMinref(compAnn_list_elmt( it1 )) + 1;
//                 if ( (containsRealLine==0) && ( (nbSolsInIt == 1) || (nbSolsInIt1 == 1) ) ) {
//                     it1 = compAnn_list_next(it1);
//                     continue;
//                 }
//                 
//                 /* check if the intersection intersects at least one annulus of the third group */
//                 /*   and if the conjugate of the intersection intersects at least one annulus of the third group */
//                 
//                 /* find first annulus of annulii2 intersecting that annulus of annulii1 */
//                 compAnn_list_iterator it2 = compAnn_list_begin( annulii2 );
//                 
//                 intersect = 0;
//                 int intersectConj = 1;
//                 compApp_set(interConj, inter);
//                 realApp_neg(compApp_imagref(interConj), compApp_imagref(interConj));
//                 
//                 while ( (it2!=compAnn_list_end() )&& ( (intersect==0) || (intersectConj==0) ) ) {
//                     if (intersect == 0)
//                         intersect = ccluster_is_compApp_in_compAnn( inter, compAnn_list_elmt( it2 ), CCLUSTER_DEFAULT_PREC);
//                     if (intersectConj == 0)
//                         intersectConj = ccluster_is_compApp_in_compAnn( interConj, compAnn_list_elmt( it2 ), CCLUSTER_DEFAULT_PREC);
//                     it2 = compAnn_list_next(it2);
//                 }
//                 if (intersect && intersectConj) {
//                     nbInt ++;
//                     nbInterA0++;
//                 }
//             }
//             it1 = compAnn_list_next(it1);
//         }
//         printf("### number of intersections: %d\n", nbInterA0);
//         it = compAnn_list_next(it);
//     }
//     
//     compApp_clear(inter);
//     compApp_clear(interConj);
//     
//     printf("# number of intersections: %d\n", nbInt);
//     printf("#time in computing intersections: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    /* fin test */
    
// int ccluster_compBox_intersects_atMost_one( const compBox_t b, int nbList ){
//     int ind = 0;
//     while ( ( ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )<=1 ) )
//         ind ++;
//     return (ind==nbList);
// }
// 
// int ccluster_compBox_intersects_atMost( const compBox_t b, int nb, int nbList ){
//     int ind = 0;
//     while ( ( ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )<=nb ) )
//         ind ++;
//     return (ind==nbList);
// }
// 
// int ccluster_compBox_intersects_atMost_spe( const compBox_t b, slong nb1, int nb2, int nbList ){
//     if (compAnn_list_get_size( compBox_annuliref( b,0 ) )>nb1 )
//         return 0;
//     int ind = 1;
//     while ( ( ind < nbList)
//         &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )<=nb2 ) )
//         ind ++;
//     return (ind==nbList);
// }

// int ccluster_connCmp_intersects_only_one( const connCmp_t cc, int nbList ){
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

/* prep loop for root radii */

// void ccluster_prep_loop_rootRadii( compBox_list_t bDiscarded,
//                          connCmp_list_t qMainLoop, 
//                          connCmp_list_t qPrepLoop, 
//                          connCmp_list_t discardedCcs, 
//                          cacheApp_t cache, 
//                          metadatas_t meta) {
//     
//     connCmp_list_t ltemp;
//     connCmp_list_init(ltemp);
//     
//     connCmp_ptr ctemp;
//     
// //     slong depth = 0;
//     
//     while (!connCmp_list_is_empty(qPrepLoop)) {
// //         while ((!connCmp_list_is_empty(qPrepLoop)) && (depth < 4) ) {
//         
//         ctemp = connCmp_list_pop(qPrepLoop);
//         
// //         depth = connCmp_getDepth(ctemp, metadatas_initBref(meta));
//         
//         int intersectOnlyOneFlag = ccluster_connCmp_intersects_only_one( ctemp, 3);
//         
// //         printf(" depth: %ld, intersectOnlyOne: %d\n", depth, intersectOnlyOneFlag);
// //         printf(" nb of cc in Preploop: %d, nb in mainloop %d\n", connCmp_list_get_size(qPrepLoop), connCmp_list_get_size(qMainLoop));
//         
//         if ( intersectOnlyOneFlag )
//             connCmp_list_insert_sorted(qMainLoop, ctemp);
//         else {
//             ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
//             
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ctemp);
//             ccluster_free(ctemp);
//         }        
//     }
//     
//     connCmp_list_clear(ltemp);
// }
