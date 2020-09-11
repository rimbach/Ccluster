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

int risolate_connCmp_intersects_only_one( const connCmp_t cc, int nbList ){
    
    int res = 1;
    slong indMins[GEOMETRY_NB_ANN_PER_BOX];
    compBox_ptr bcur;
    compBox_list_iterator itb;
    
    itb = compBox_list_begin( connCmp_boxesref( cc ) );
    bcur= compBox_list_elmt(itb);
    /* each list of annulii contains at least one annulus */
    for (int ind = 0; ind < nbList; ind++) {
        /* check if only one annulus */
        if ( compAnn_list_get_size(compBox_annuliref(bcur, ind))>1 )
            res = 0;
        indMins[ind] = compAnn_indMinref(compAnn_list_first(compBox_annuliref(bcur, ind)));
    }
    itb = compBox_list_next(itb);
    while ( (res==1) && ( itb != compBox_list_end() ) ){
        bcur= compBox_list_elmt(itb);
        for (int ind = 0; ind < nbList; ind++) {
            if ( compAnn_list_get_size(compBox_annuliref(bcur, ind))>1 )
                res = 0;
            if ( compAnn_indMinref(compAnn_list_first(compBox_annuliref(bcur, ind))) != indMins[ind] )
                res = 0;
        }
        itb = compBox_list_next(itb);
    }
    
    return res;
}

int risolate_compBox_intersects_only_one( const compBox_t b, int nbList ){
    int ind = 0;
    while ( ( ind < nbList)
        &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )==1 ) )
        ind ++;
    return (ind==nbList);
}

int risolate_compBox_intersects_atLest_one( const compBox_t b, int nbList ){
    int ind = 0;
    while ( ( ind < nbList)
        &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )>0 ) )
        ind ++;
    return (ind==nbList);
}

int risolate_compBox_nbMsols( const compBox_t b, int nbList ){
    int res = compBox_get_nbMSol(b);
//     printf("nb of sols in the box begin: %d\n", res);
    int ind=0;
    while ( (ind < nbList)
        &&  (compAnn_list_get_size( compBox_annuliref( b,ind ) )==1 ) ){
        /* the maximum number of sol is the min of the number of sols in the annulii */
        compAnn_list_iterator ita = compAnn_list_begin( compBox_annuliref( b,ind ) );
        compAnn_ptr ann = compAnn_list_elmt(ita);
        int temp = compAnn_indMaxref( ann ) - compAnn_indMinref( ann ) +1;           
        if (temp < res)
            res = temp;
        ind ++;
    }
//     printf("nb of sols in the box end  : %d\n", res);
    return res;    
}

/* assume that annulii in the first list are centered in a+0i where a> any root */
/* as a consequence, if a cc intersects a unique annulus of the first list that contains one root, */
/* this root is real, and it is in the intersection of that annulus with the real axis */
/* in the cc */
/* NOT USED */
int risolate_connCmp_ContainsOneRealRoot( const connCmp_t cc ) {
    
    
    if (! risolate_connCmp_intersects_only_one( cc, 1 ) )
        return -1;
    
    compBox_list_iterator itb = compBox_list_begin(connCmp_boxesref(cc));
    compBox_ptr btemp = compBox_list_elmt(itb);
    compAnn_ptr atemp = compAnn_list_first(compBox_annuliref(btemp, 0));
    if ( compAnn_indMaxref(atemp) == compAnn_indMinref(atemp) )
        return 1;
    else 
        return -1;
}

void connCmp_mergeAnnulii( compAnn_list_t dest, int ind, const connCmp_t cc ){
    
    compAnn_list_iterator itadest, itab;
    compBox_list_iterator itb = compBox_list_begin(connCmp_boxesref(cc));
    
//     printf("---Merge: \n");
    
    while( itb != compBox_list_end() ) {
        
//         printf("nbox; list: \n");
//         compAnn_list_printd(compBox_annuliref(compBox_list_elmt(itb), ind), 10);
//         printf("\n\n");
    
        itab = compAnn_list_begin( compBox_annuliref(compBox_list_elmt(itb), ind) );
        while (itab!= compAnn_list_end() ){
            
            int isIn = 0;
            itadest = compAnn_list_begin(dest);
            while ( (isIn==0) && (itadest!= compAnn_list_end()) ){
                if ( compAnn_list_elmt(itadest) == compAnn_list_elmt(itab) )
                    isIn = 1;
                itadest = compAnn_list_next(itadest);
            }
            if (isIn==0)
                compAnn_list_insert_sorted(dest, compAnn_list_elmt(itab));
            itab = compAnn_list_next(itab);
        }
        
        itb = compBox_list_next(itb);
    }
//     printf("dest at the end: \n");
//     compAnn_list_printd(dest, 10);
//     printf("\n\n");
}

void risolate_discard_compBox_list_prepLoop_rootRadii( compBox_list_t boxes, 
                                                       compBox_list_t bDiscarded,
                                                       cacheApp_t cache, 
                                                       slong prec, 
                                                       metadatas_t meta){
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
    compBox_ptr btemp;
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        
        risolate_compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
            
        /* check if btemp intersects at least an annulus of each list */
        if ( risolate_compBox_intersects_atLest_one( btemp, 1 ) == 0 ){
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
        }
        
        btemp->nbMSol = risolate_compBox_nbMsols( btemp, 1 );
        compBox_list_push(ltemp, btemp);
    }
        
 
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
}

void risolate_bisect_connCmp_prepLoop_rootRadii( connCmp_list_t dest, 
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
    
    risolate_discard_compBox_list_prepLoop_rootRadii(subBoxes, bDiscarded, cache, prec, meta);
    
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    
    while (!connCmp_list_is_empty(ltemp)){
        ctemp = connCmp_list_pop(ltemp);
        connCmp_list_push(dest, ctemp);
    }
    
    compBox_list_clear(subBoxes);
    connCmp_list_clear(ltemp);
}

void risolate_prep_loop_rootRadii( compBox_list_t bDiscarded, 
                                   connCmp_list_t qResult, 
                                   connCmp_list_t qPrepLoop, 
                                   connCmp_list_t discardedCcs, 
                                   cacheApp_t cache, 
                                   metadatas_t meta) {
    
    connCmp_ptr ctemp;
    connCmp_list_t ltemp;
    connCmp_list_init(ltemp);
    
    int widthFlag;
    int intersectOnlyOneFlag;
    int containsOneRealRootFlag;
    int separationFlag;
    int compactFlag;
    int sgn;
    
    compBox_t componentBox;
    compDsk_t ccDisk, fourCCDisk;
    realRat_t four;
    realApp_t widthCtemp;
    realApp_t widthAnn;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    compDsk_init(fourCCDisk);
    realRat_init(four);
    realApp_init(widthCtemp);
    realApp_init(widthAnn);
    
    compBox_ptr bcur, nbox;
    compAnn_ptr acur;
    
    realRat_set_si(four, 4, 1);
    
    
    while (!connCmp_list_is_empty(qPrepLoop)) {
        ctemp = connCmp_list_pop(qPrepLoop);
        
        intersectOnlyOneFlag = risolate_connCmp_intersects_only_one( ctemp, 1 );
        containsOneRealRootFlag = -1;
        widthFlag = 0;
        separationFlag = 0;
        compactFlag    = compBox_list_get_size(connCmp_boxesref(ctemp)) <= 3;
        
        
        if (intersectOnlyOneFlag == 1) {
            bcur = compBox_list_first(connCmp_boxesref(ctemp));
            acur = compAnn_list_first(compBox_annuliref(bcur, 0));
            
            /* check if ctemp contains the Zero annulus */
            if ( realApp_is_zero(compAnn_radInfref(acur)) && realApp_is_zero(compAnn_radSupref(acur)) ) {
//                 containsOneRealRootFlag = 1;
                /* in this case, let containsOneRealRootFlag = multiplicity of 0 */
                containsOneRealRootFlag = compAnn_indMaxref(acur) - compAnn_indMinref(acur) +1;
            }
            
            sgn = realRat_sgn(connCmp_supReref(ctemp));
            if ( (sgn<0) && (compAnn_rrInNeref(acur) == 1) )
                containsOneRealRootFlag = 1;
            if ( (sgn<0) && (compAnn_rrInNeref(acur) == 0) ) {
                connCmp_clear(ctemp);
                ccluster_free(ctemp);
                continue;
            }
            
            sgn = realRat_sgn(connCmp_infReref(ctemp));
            if ( (sgn>0) && (compAnn_rrInPoref(acur) == 1) )
                containsOneRealRootFlag = 1;
            if ( (sgn>0) && (compAnn_rrInPoref(acur) == 0) ) {
                connCmp_clear(ctemp);
                ccluster_free(ctemp);
                continue;
            }
            
            realApp_set_realRat( widthCtemp, compBox_bwidthref(bcur), CCLUSTER_DEFAULT_PREC );
            realApp_set(widthAnn, compAnn_radSupref(acur));
            realApp_sub(widthAnn, widthAnn, compAnn_radInfref(acur), CCLUSTER_DEFAULT_PREC );
            widthFlag = realApp_lt(widthCtemp, widthAnn);
            if ((widthFlag != 1)&&( realApp_ge(widthCtemp, widthAnn) != 1 ))
                widthFlag = 1;
        }
        
//         if ((intersectOnlyOneFlag == 1) && (containsOneRealRootFlag==1)){
        if ((intersectOnlyOneFlag == 1) && (containsOneRealRootFlag >= 1)){
            connCmp_risolate_componentBox(componentBox, ctemp, metadatas_initBref(meta));
            compBox_get_containing_dsk(ccDisk, componentBox);
            compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
            separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qPrepLoop, discardedCcs);
        }
        
        if (intersectOnlyOneFlag == 0) {
            risolate_bisect_connCmp_prepLoop_rootRadii( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ctemp);
            ccluster_free(ctemp);
        } 
        else {
            if ( ( widthFlag == 0 )
               &&( (containsOneRealRootFlag==-1)||(separationFlag==0) )
            ){
                risolate_bisect_connCmp_prepLoop_rootRadii( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
                while (!connCmp_list_is_empty(ltemp))
                    connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
                connCmp_clear(ctemp);
                ccluster_free(ctemp);
            }
            else {
                
                if (compactFlag==0){ /* replace the CC by a cc containing only one box */
//                     
//                     printf("# Non compact cc: nb of boxes: %d\n", compBox_list_get_size(connCmp_boxesref(ctemp)) );
                    bcur = compBox_list_first(connCmp_boxesref(ctemp));
//                     
                    nbox = (compBox_ptr) ccluster_malloc (sizeof(compBox));
                    compBox_init(nbox);
                    connCmp_risolate_componentBox( nbox, ctemp, metadatas_initBref(meta));
                    compBox_nbMSolref( nbox ) = risolate_compBox_nbMsols( bcur, 1 );
                    /* copy annulii 0 */
                    compBox_copy_annulii(nbox, 0, compBox_annuliref(bcur, 0));
                    /*delete ctemp*/
                    connCmp_clear(ctemp);
                    ccluster_free(ctemp);
                    /* create the new connected component */
                    ctemp = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
                    connCmp_init_compBox(ctemp, nbox);
                    
                }
                
                connCmp_nSolsref(ctemp) = containsOneRealRootFlag;
                connCmp_list_insert_sorted(qResult, ctemp);
            }
        }
            
    }
    connCmp_list_clear(ltemp);
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    realRat_clear(four);
    realApp_clear(widthCtemp);
    realApp_clear(widthAnn);
    
}

slong risolate_exclusion_rootRadii( connCmp_list_t qCover,
                                   cacheApp_t cache, 
                                   metadatas_t meta){
    
    connCmp_list_t ltemp;
    connCmp_ptr ccur;
    compBox_ptr box;
    
    compDsk_t ccDisk;
    compDsk_init(ccDisk);
    
    connCmp_list_init(ltemp);
    
    compBox_t componentBox;
    compBox_init(componentBox);
    
    slong depth = 0;
    tstar_res res;
    res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    while (!connCmp_list_is_empty(qCover)) {
          ccur = connCmp_list_pop(qCover);
          
          if (connCmp_nSolsref(ccur) == -1) {
              
            connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));  
              
            box  = compBox_list_first( connCmp_boxesref(ccur) );
            /*get containing disk */
            risolate_compBox_get_containing_dsk(ccDisk, componentBox);
            /* do an exclusion test */
            res = tstar_real_interface( cache, ccDisk, compBox_get_nbMSol(box), 1, 0, res.appPrec , depth, meta);
//             res.appPrec = res.appPrec/2;
            
            if (res.nbOfSol==0) { /* clear ccur */
                if (metadatas_haveToCount(meta)){
                        metadatas_add_discarded( meta, depth);
                }
                connCmp_clear(ccur);
                ccluster_free(ccur);
            } else {
                if (res.nbOfSol>0) {
                    compBox_nbMSolref(box) = res.nbOfSol;
                    connCmp_nSolsref(ccur) = res.nbOfSol;
                }
                connCmp_appPrref(ccur) = res.appPrec;
                connCmp_list_insert_sorted(ltemp, ccur);
            }
          } else {
              connCmp_list_insert_sorted(ltemp, ccur);
          }
    }
    
    
    connCmp_list_swap(ltemp, qCover);
    
    connCmp_list_clear(ltemp);
    compDsk_clear(ccDisk);
    
    compBox_clear(componentBox);
    
    return res.appPrec;
}

void risolate_algo_global_rootRadii  ( connCmp_list_t qResults, 
                                       compBox_list_t bDiscarded,
                                       compAnn_list_t annulii,
                                       const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    clock_t start = clock();
    clock_t start2 = clock();
    
    realRat_t delta;
    realRat_init(delta);
    
    slong degree = cacheApp_getDegree( cache );
//     realRat_set_si(delta, 1, degree);
    realRat_set_si(delta, 1, degree*degree);
//     realRat_set_si(delta, 1, degree*degree*degree);
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    /* heuristic to predict the precision: at least the degree */
    while (prec<degree)
        prec = 2*prec;
    prec = realIntRootRadii_rootRadii( annulii, 0, cache, delta, prec );
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#time in computing 1st RootRadii       : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        printf("#precision needed                      : %ld\n", prec);
        if (metadatas_getVerbo(meta)>=3) {
            printf("#Annulii: ");
            compAnn_list_printd(annulii, 10);
            printf("\n\n");
        }
    }
    
    /* use annulii to get a sharp upper bound on the roots */
    slong upperBound = realApp_ceil_si(compAnn_radSupref(compAnn_list_last(annulii)), prec);
    start2 = clock();
    
    realIntRootRadii_connectedComponents( annulii, prec );
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    }
    start2 = clock();
    
    realIntRootRadii_containsRealRoot( annulii, cache, prec );
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        if (metadatas_getVerbo(meta)>=3) {
            printf("#Annulii: ");
            compAnn_list_printd(annulii, 10);
            printf("\n\n");
        }
    }
    start2 = clock();
    
    compBox_ptr box;
    box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(box);
    compBox_set(box, initialBox);
    compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
    /* set width to 2*upperBound */
    realRat_set_si( compBox_bwidthref(box), 2*upperBound, 1 );
    
    compBox_copy_annulii(box, 0, annulii);
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, box);
    
    connCmp_list_t qPrepLoop, qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
    connCmp_list_init(qPrepLoop);
    connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qPrepLoop, initialCC);
    risolate_prep_loop_rootRadii(bDiscarded, qMainLoop, qPrepLoop, discardedCcs, cache, meta);
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#time in prep loop: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        printf("#Nb of CC: %d\n", connCmp_list_get_size(qMainLoop));
        if (metadatas_getVerbo(meta)>=3) {
        connCmp_list_iterator itc = connCmp_list_begin(qMainLoop);
        while ( itc!= connCmp_list_end() ){
            printf("#--- CC with %2d sols:\n", connCmp_nSolsref( connCmp_list_elmt(itc) ) );
// //               compBox_ptr b= compBox_list_first(connCmp_boxesref( connCmp_list_elmt(itc) ));
// //               for (int ind = 0; ind < GEOMETRY_NB_ANN_PER_BOX; ind++){
// //                   if (compAnn_list_get_size(compBox_annuliref(b, ind))>0){
// //                       compAnn_ptr ann = compAnn_list_first( compBox_annuliref(b, ind) );
// //                       compAnn_printd(ann, 10);
// //                       printf("#\n");
// //                   }
// //               }
                itc = connCmp_list_next(itc);
            }
        }
    }
    
    start2 = clock();
    
    prec = risolate_exclusion_rootRadii( qMainLoop, cache, meta);
    
    if (metadatas_getVerbo(meta)>=2) {
        printf("#time in exclusion tests               : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        printf("#precision needed: %ld\n", prec);
        printf("#Nb of CC: %d\n", connCmp_list_get_size(qMainLoop));
        if (metadatas_getVerbo(meta)>=3) {
        connCmp_list_iterator itc = connCmp_list_begin(qMainLoop);
        while ( itc!= connCmp_list_end() ){
            printf("#--- CC with %2d sols:\n", connCmp_nSolsref( connCmp_list_elmt(itc) ) );
// //               compBox_ptr b= compBox_list_first(connCmp_boxesref( connCmp_list_elmt(itc) ));
// //               for (int ind = 0; ind < GEOMETRY_NB_ANN_PER_BOX; ind++){
// //                   if (compAnn_list_get_size(compBox_annuliref(b, ind))>0){
// //                       compAnn_ptr ann = compAnn_list_first( compBox_annuliref(b, ind) );
// //                       compAnn_printd(ann, 10);
// //                       printf("#\n");
// //                   }
// //               }
                itc = connCmp_list_next(itc);
            }
        }
    }
    start2 = clock();
    
    risolate_main_loop( qResults,  bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
//     connCmp_list_swap(qResults, qMainLoop);
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    realRat_clear(delta);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

// OLD RISOLATE ROOTRADII

// void risolate_prep_loop_rootRadii2( connCmp_list_t qCover, 
//                                    const compBox_t initialBox,
//                                    cacheApp_t cache, 
//                                    metadatas_t meta){
//     
//     /* qCover contains a list of cc intersecting one annulus */
//     connCmp_ptr ccur;
//     connCmp_list_t ltemp;
//     connCmp_list_init(ltemp);
//     
//     compBox_ptr box;
//     compAnn_ptr ann;
//     compBox_list_iterator itb;
//     compAnn_list_ptr annl;
//     
//     while (!connCmp_list_is_empty(qCover)) {
//         ccur = connCmp_list_pop(qCover);
//         itb  = compBox_list_begin( connCmp_boxesref( ccur ) );
//         box  = compBox_list_elmt( itb );
//         annl = compBox_annuliref( box );
//         ann  = compAnn_list_first(annl);
//         /* transform ccur in a cc with a unique box */
//         box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//         compBox_init(box);
//         connCmp_risolate_componentBox( box, ccur, initialBox);
//         compBox_init_annuli(box);
//         compBox_copy_annuli(box, annl);
//         compBox_nbMSolref( box ) = compAnn_indMaxref(ann) - compAnn_indMinref(ann) + 1;
//         /*delete Ccur*/
//         while ( itb!=compBox_list_end() ){
//           compBox_clear_annuli(compBox_list_elmt(itb));
//           itb = compBox_list_next(itb);
//         }
//         connCmp_clear(ccur);
//         ccluster_free(ccur);
//         /* create the new connected component */
//         ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//         connCmp_init_compBox(ccur, box);
//         /*check if ccur contains a unique real root */
//         if (compAnn_rrInPoref(ann)>-1)
//           connCmp_nSolsref(ccur) = 1;
//           
//         /* push ccur in ltemp */
//         connCmp_list_insert_sorted(ltemp, ccur);
//         
//     }
//     
//     connCmp_list_swap(ltemp, qCover);
//     connCmp_list_clear(ltemp);
//     
// }

void risolate_prep_loop_rootRadii_old( compBox_list_t bDiscarded,
                                   connCmp_list_t qCover, 
                                   const compBox_t initialBox,
                                   const compAnn_list_t annulii,
                                   cacheApp_t cache, 
                                   metadatas_t meta) {
    
    connCmp_list_t qLoop, ltemp;
    connCmp_ptr ccur;
    compBox_ptr box;
    compAnn_ptr ann;
    
    connCmp_list_init(qLoop);
    connCmp_list_init(ltemp);
    
    box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(box);
    compBox_set(box, initialBox);
    compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
//     compBox_init_annuli(box);
    compBox_copy_annulii(box, 0, annulii);
    
    ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(ccur, box);
    connCmp_list_push(qLoop, ccur);
    
    compBox_list_iterator itb;
    compAnn_list_ptr annl;
    int intersectUniqueAnnulus = 0;
    
    while (!connCmp_list_is_empty(qLoop)) {
          ccur = connCmp_list_pop(qLoop);
          /* check if ccur intersects a unique annulus */  
          itb = compBox_list_begin( connCmp_boxesref( ccur ) );
          box = compBox_list_elmt( itb ) ;                     /* at least one box in ccur */
          intersectUniqueAnnulus = (compAnn_list_get_size( compBox_annuli0ref( box ) ) == 1);
          ann = compAnn_list_first(compBox_annuli0ref( box ));  /* box intersects at least one annulus */
          itb = compBox_list_next( itb );
          while ( ( intersectUniqueAnnulus == 1 ) && 
                    (itb != compBox_list_end()   ) ) {
              
              box = compBox_list_elmt( itb ) ;
              annl = compBox_annuli0ref( box ); 
              intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_get_size( annl ) == 1) ;
              intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_first(annl) == ann ) ;
              itb = compBox_list_next( itb );
          }
          
          if (intersectUniqueAnnulus==1) {
              
              itb  = compBox_list_begin( connCmp_boxesref( ccur ) );
              box  = compBox_list_elmt( itb );
              annl = compBox_annuli0ref( box );
              ann  = compAnn_list_first(annl);
              
              /* transform ccur in a cc with a unique box */
              box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
              compBox_init(box);
              connCmp_risolate_componentBox( box, ccur, initialBox);
              compBox_copy_annulii(box, 0, annl);
              compBox_nbMSolref( box ) = compAnn_indMaxref(ann) - compAnn_indMinref(ann) + 1;
              /*delete Ccur*/
              connCmp_clear(ccur);
              ccluster_free(ccur);
              /* create the new connected component */
              ccur = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
              connCmp_init_compBox(ccur, box);
              /*check if ccur contains a unique real root */
              if ( (compAnn_indMaxref(ann)==compAnn_indMinref(ann)) && (compAnn_rrInPoref(ann)>-1) )
                connCmp_nSolsref(ccur) = 1;
                
              /* push ccur in qCover */
                connCmp_list_insert_sorted(qCover, ccur);
              
          } else { /* bisect ccur */
                
                realIntRootRadii_bisect_connCmp( ltemp, ccur);
                while (!connCmp_list_is_empty(ltemp))
                    connCmp_list_insert_sorted(qLoop, connCmp_list_pop(ltemp));
                connCmp_clear(ccur);
                ccluster_free(ccur);
                
          }
    }
    
    connCmp_list_clear(ltemp);
    connCmp_list_clear(qLoop);
}

void risolate_algo_global_rootRadii_old( connCmp_list_t qResults, 
                                     compBox_list_t bDiscarded,
                                     compAnn_list_t annulii,
                                     compAnn_list_t annulii1,
                                     compAnn_list_t annulii2,
                                     const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    clock_t start = clock();
    clock_t start2 = clock();
    
//     compAnn_list_t annulii;
    connCmp_list_t qCover;
    realRat_t delta;
    
//     compAnn_list_init(annulii);
    connCmp_list_init(qCover);
    realRat_init(delta);
    
    slong degree = cacheApp_getDegree( cache );
//     realRat_set_si(delta, 1, degree);
    realRat_set_si(delta, 1, degree*degree);
//     realRat_set_si(delta, 1, degree*degree*degree);
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    prec = realIntRootRadii_rootRadii( annulii, 0, cache, delta, prec );
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#time in computing RootRadii           : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        printf("#precision needed: %ld\n", prec);
//     }
    start2 = clock();
    
    realIntRootRadii_connectedComponents( annulii, prec );
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//     }
    start2 = clock();
    
//     printf("Annulii: ");
//     compAnn_list_printd(annulii, 10);
//     printf("\n\n");
    
    realIntRootRadii_containsRealRoot( annulii, cache, prec );
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//     }
    start2 = clock();
    
//     printf("Annulii: ");
//     compAnn_list_printd(annulii, 10);
//     printf("\n\n");
    
    risolate_prep_loop_rootRadii_old( bDiscarded, qCover, initialBox, annulii, cache, meta);
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#time in covering intervals with boxes : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//     }
    start2 = clock();
    
//     connCmp_list_iterator itc;
// //     display qCover
//     itc = connCmp_list_begin(qCover);
// //     if (metadatas_getVerbo(meta)>3) {
//         printf("#Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
// //     }
//     while( itc!= connCmp_list_end() ){
//         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
//         printf("\n");
//         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
//         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))),0)) );
// //            
//         itc = connCmp_list_next( itc );
//     }
    printf("\n\n");
        
    prec = risolate_exclusion_rootRadii( qCover, cache, meta);
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#time in exclusion tests               : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
        printf("#precision needed: %ld\n", prec);
//     }
    start2 = clock();
    
    /* display qCover */
//     itc = connCmp_list_begin(qCover);
//     if (metadatas_getVerbo(meta)>3) {
        printf("#Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
//     }
//     while( itc!= connCmp_list_end() ){
//         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
//         printf("\n");
//         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
//         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))))) ); 
//         itc = connCmp_list_next( itc );
//     }
//     printf("\n\n");
    
    
//     if (metadatas_getVerbo(meta)>3) {
        printf("#total time in root radii              : %f \n", ((double) (clock() - start))/CLOCKS_PER_SEC );
//     }
    
    connCmp_list_t discardedCcs;
    connCmp_list_init(discardedCcs);
    
    /* main loop */
    risolate_main_loop( qResults, bDiscarded,  qCover, discardedCcs, eps, cache, meta);
    
    connCmp_list_clear(qCover);
    realRat_clear(delta);
//     compAnn_list_clear(annulii);
    connCmp_list_clear(discardedCcs);
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
//     compBox_ptr box;
//     box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//     compBox_init(box);
//     compBox_set(box, initialBox);
//     compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
    
//     connCmp_ptr initialCC;
//     initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//     connCmp_init_compBox(initialCC, box);
/*    
    connCmp_list_t qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
    connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qMainLoop, initialCC);
    risolate_main_loop( qResults,  qMainLoop, discardedCcs, eps, cache, meta);
    
    
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(discardedCcs);
    
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));*/
}
