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

void risolate_compBox_get_containing_dsk( compDsk_t d, const compBox_t b) {
    compDsk_set_compRat_realRat(d, compBox_centerref(b), compBox_bwidthref(b) );
    realRat_div_ui( compDsk_radiusref(d), compDsk_radiusref(d), 2 );
}

slong risolate_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
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
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        
        risolate_compBox_get_containing_dsk(bdisk, btemp);
        /* compBox_get_containing_dsk(bdisk, btemp); */
        
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
//         printf("nbMSol: %d\n", (int) compBox_get_nbMSol(btemp) );
        
        res = tstar_real_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta);  
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
    
    return res.appPrec;
}

void risolate_bisect_connCmp( connCmp_list_t dest, 
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

    prec = risolate_discard_compBox_list( subBoxes, bDiscarded, cache, prec, meta);
    
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

void risolate_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t qPrepLoop, 
                         connCmp_list_t discardedCcs, 
                         cacheApp_t cache, 
                         metadatas_t meta) {
    
    connCmp_list_t ltemp;
    realRat_t halfwidth, diam;
    connCmp_list_init(ltemp);
    realRat_init(halfwidth);
    realRat_init(diam);
    
    realRat_set_si(halfwidth, 1, 2);
    realRat_mul(halfwidth, halfwidth, compBox_bwidthref(metadatas_initBref(meta)) );
    
    connCmp_ptr ctemp;
    
    while (!connCmp_list_is_empty(qPrepLoop)) {
        
        ctemp = connCmp_list_pop(qPrepLoop);
        connCmp_diameter(diam, ctemp);
        
        if ( connCmp_is_confined(ctemp, metadatas_initBref(meta)) && (realRat_cmp(diam, halfwidth)<0) )
            connCmp_list_insert_sorted(qMainLoop, ctemp);
        else {
            risolate_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
            
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ctemp);
            ccluster_free(ctemp);
        }        
    }
    
    connCmp_list_clear(ltemp);
    realRat_clear(halfwidth);
    realRat_clear(diam);
}

void risolate_prep_loop_rootRadii( connCmp_list_t qCover, 
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
    compBox_init_annuli(box);
    compBox_copy_annuli(box, annulii);
    
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
          intersectUniqueAnnulus = (compAnn_list_get_size( compBox_annuliref( box ) ) == 1);
          ann = compAnn_list_first(compBox_annuliref( box ));  /* box intersects at least one annulus */
          itb = compBox_list_next( itb );
          while ( ( intersectUniqueAnnulus == 1 ) && 
                    (itb != compBox_list_end()   ) ) {
              
              box = compBox_list_elmt( itb ) ;
              annl = compBox_annuliref( box ); 
              intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_get_size( annl ) == 1) ;
              intersectUniqueAnnulus = intersectUniqueAnnulus & (compAnn_list_first(annl) == ann ) ;
              itb = compBox_list_next( itb );
          }
          
          if (intersectUniqueAnnulus==1) {
              
              itb  = compBox_list_begin( connCmp_boxesref( ccur ) );
              box  = compBox_list_elmt( itb );
              annl = compBox_annuliref( box );
              ann  = compAnn_list_first(annl);
              
              /* transform ccur in a cc with a unique box */
              box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
              compBox_init(box);
              connCmp_risolate_componentBox( box, ccur, initialBox);
              compBox_init_annuli(box);
              compBox_copy_annuli(box, annl);
              compBox_nbMSolref( box ) = compAnn_indMaxref(ann) - compAnn_indMinref(ann) + 1;
              /*delete Ccur*/
              while ( itb!=compBox_list_end() ){
                compBox_clear_annuli(compBox_list_elmt(itb));
                itb = compBox_list_next(itb);
              }
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

slong risolate_exclusion_rootRadii( connCmp_list_t qCover,
                                   cacheApp_t cache, 
                                   metadatas_t meta){
    
    connCmp_list_t ltemp;
    connCmp_ptr ccur;
    compBox_ptr box;
    
    compDsk_t ccDisk;
    compDsk_init(ccDisk);
    
    connCmp_list_init(ltemp);
    
    slong depth = 0;
    tstar_res res;
    res.appPrec = CCLUSTER_DEFAULT_PREC;
    
    while (!connCmp_list_is_empty(qCover)) {
          ccur = connCmp_list_pop(qCover);
          
          if (connCmp_nSolsref(ccur) == -1) {
              
            box  = compBox_list_first( connCmp_boxesref(ccur) );
            /*get containing disk */
            risolate_compBox_get_containing_dsk(ccDisk, box);
            /* do an exclusion test */
            res = tstar_real_interface( cache, ccDisk, compBox_get_nbMSol(box), 1, 0, res.appPrec , depth, meta);
//             res.appPrec = res.appPrec/2;
            
            if (res.nbOfSol==0) { /* clear ccur */
                if (metadatas_haveToCount(meta)){
                        metadatas_add_discarded( meta, depth);
                }
                compBox_clear_annuli(box);
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
    
    return res.appPrec;
}

void risolate_main_loop( connCmp_list_t qResults,  
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
    realRat_t three, four, threeWidth;
    compRat_t initPoint;
    connCmp_list_t ltemp;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    compDsk_init(fourCCDisk);
    realRat_init(three);
    realRat_init(four);
    realRat_init(threeWidth);
    compRat_init(initPoint);
    connCmp_list_init(ltemp);
    
    connCmp_ptr ccur;
    
    clock_t start=clock();
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    
    while (!connCmp_list_is_empty(qMainLoop)) {
        
        //         if (metadatas_getVerbo(meta)>0) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }

        resNewton.nflag = 0;
        
        ccur = connCmp_list_pop(qMainLoop);
        
        connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
        compBox_get_containing_dsk(ccDisk, componentBox);
        
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
        separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
        
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
      
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        
        sepBoundFlag   = (realRat_cmp( compBox_bwidthref(componentBox), metadatas_getSepBound(meta))<=0);
        
        if (metadatas_getVerbo(meta)>3) {
            printf("---depth: %d\n", (int) depth);
            printf("------component Box:       "); compBox_print(componentBox); printf("\n");
            printf("------nb of roots:         %d\n", connCmp_nSolsref(ccur));
            printf("------last newton success: %d\n", connCmp_newSuref(ccur));
            printf("------separation Flag:     %d\n", separationFlag);
            printf("------widthFlag:           %d\n", widthFlag); 
            printf("------compactFlag:         %d\n", compactFlag);
            printf("------sepBoundFlag:        %d\n", sepBoundFlag);
            compBox_list_print( connCmp_boxesref(ccur) );
            printf("\n");
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
//         if ((separationFlag)) {
//             printf("depth: %d, connCmp_nSolsref(ccur): %d, prec: %d\n", (int) depth, (int) connCmp_nSolsref(ccur), (int) prec);
//             
//             if (connCmp_nSolsref(ccur)==-1){ /* do not do that because one can loose complex roots */
               if (connCmp_nSolsref(ccur)!=1){
                
//                 resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                resTstar = tstar_real_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                 if (metadatas_getVerbo(meta)>3)
//                     printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
//                 ???
                prec = resTstar.appPrec;
            }
                
//             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
//             ???
//             prec = resTstar.appPrec;
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
             ( (!widthFlag) 
               ||( connCmp_nSols(ccur)== cacheApp_getDegree(cache) )
               ||( !sepBoundFlag && (connCmp_nSols(ccur)>1))   
            )  )
//             &&!( metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) ) 
        ) {
        
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
    
            }
            else {
                connCmp_newSuref(ccur) = 0;
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
                printf("------validated with %d roots\n", connCmp_nSols(ccur));
            }
        } else
        if ( (connCmp_nSols(ccur)>0) && (connCmp_nSols(ccur)<cacheApp_getDegree(cache)) 
             && separationFlag && widthFlag && compactFlag
             && ( sepBoundFlag || (connCmp_nSols(ccur)==1)) ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
            if (metadatas_getVerbo(meta)>3) {
                printf("------validated with %d roots\n", connCmp_nSols(ccur));
            }
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            connCmp_decrease_nwSpd(ccur);
//             if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0)
//                 connCmp_decrease_nwSpd(ccur);
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        else {
//             if (connCmp_nSols(ccur)==0) 
//                 printf("ici\n");
// #ifdef CCLUSTER_HAVE_PTHREAD
//             ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ccur);
//             ccluster_free(ccur);
// #else
            risolate_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cache, meta,1);
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
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(threeWidth);
    compRat_clear(initPoint);
    connCmp_list_clear(ltemp);
    
}

void risolate_algo( connCmp_list_t qResults, 
                    compBox_list_t bDiscarded,
                    const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
//     chronos_tic_CclusAl(metadatas_chronref(meta));
    clock_t start = clock();
    
    realRat_t factor;
    realRat_init(factor);
    realRat_set_si(factor, 5, 4);
    
    compBox_ptr bEnlarged;
    bEnlarged = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(bEnlarged);
    compBox_inflate_realRat(bEnlarged, initialBox, factor);
    compBox_nbMSolref(bEnlarged) = cacheApp_getDegree ( cache );
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, bEnlarged);
    
    connCmp_list_t qMainLoop, qPrepLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
    connCmp_list_init(qPrepLoop);
    connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qPrepLoop, initialCC);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster preploop: \n");
    risolate_prep_loop( bDiscarded, qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster mainloop: \n");
    risolate_main_loop( qResults, bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
    
    
    realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

void risolate_algo_global( connCmp_list_t qResults, 
                           compBox_list_t bDiscarded,
                           const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    clock_t start = clock();
    
//     realRat_t factor;
//     realRat_init(factor);
//     realRat_set_si(factor, 5, 4);
    
    compBox_ptr box;
    box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
    compBox_init(box);
    compBox_set(box, initialBox);
    compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, box);
    
    connCmp_list_t qMainLoop, discardedCcs;
    connCmp_list_init(qMainLoop);
//     connCmp_list_init(qPrepLoop);
    connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qMainLoop, initialCC);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster preploop: \n");
//     ccluster_prep_loop( qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster mainloop: \n");
    risolate_main_loop( qResults,  bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
    
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

void risolate_algo_global_rootRadii( connCmp_list_t qResults, const compBox_t initialBox, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    clock_t start = clock();
    clock_t start2 = clock();
    
    compAnn_list_t annulii;
    connCmp_list_t qCover;
    realRat_t delta;
    
    compAnn_list_init(annulii);
    connCmp_list_init(qCover);
    realRat_init(delta);
    
    slong degree = cacheApp_getDegree( cache );
//     realRat_set_si(delta, 1, degree);
    realRat_set_si(delta, 1, degree*degree);
//     realRat_set_si(delta, 1, degree*degree*degree);
    
    slong prec = realIntRootRadii_rootRadii( annulii, cache, delta );
    
    printf("time in computing RootRadii           : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    printf("precision needed: %ld\n", prec);
    start2 = clock();
    
    realIntRootRadii_connectedComponents( annulii, prec );
    
    printf("time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    start2 = clock();
    
//     printf("Annulii: ");
//     compAnn_list_printd(annulii, 10);
//     printf("\n\n");
    
    realIntRootRadii_containsRealRoot( annulii, cache, prec );
    
    printf("time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    start2 = clock();
    
//     printf("Annulii: ");
//     compAnn_list_printd(annulii, 10);
//     printf("\n\n");
    
    risolate_prep_loop_rootRadii( qCover, initialBox, annulii, cache, meta);
    
    printf("time in covering intervals with boxes : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    start2 = clock();
    
    connCmp_list_iterator itc;
    /* display qCover */
//     itc = connCmp_list_begin(qCover);
    printf("Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
//     while( itc!= connCmp_list_end() ){
//         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
//         printf("\n");
//         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
//         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))))) );
//            
//         itc = connCmp_list_next( itc );
//     }
//     printf("\n\n");
        
    prec = risolate_exclusion_rootRadii( qCover, cache, meta);
    
    printf("time in exclusion tests               : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
    printf("precision needed: %ld\n", prec);
    start2 = clock();
    
    /* display qCover */
//     itc = connCmp_list_begin(qCover);
    printf("Number of CC in qCover: %d \n", connCmp_list_get_size(qCover));
//     while( itc!= connCmp_list_end() ){
//         printf("--- Box: "); compBox_print( compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))) );
//         printf("\n");
//         printf("--- nb of sols in CC: %d \n", connCmp_nSolsref(connCmp_list_elmt(itc)) );
//         printf("--- nb of inter annulus: %d \n", compAnn_list_get_size(compBox_annuliref(compBox_list_first(connCmp_boxesref(connCmp_list_elmt(itc))))) ); 
//         itc = connCmp_list_next( itc );
//     }
//     printf("\n\n");
    
    /* clear annulii information */
    itc = connCmp_list_begin(qCover);
    while( itc!= connCmp_list_end() ){
        compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( connCmp_list_elmt(itc) ) );
        while ( itb!=compBox_list_end() ){
            compBox_clear_annuli(compBox_list_elmt(itb));
            itb = compBox_list_next(itb);
        }
        itc = connCmp_list_next( itc );
    }
    
    printf("total time in root radii              : %f \n", ((double) (clock() - start))/CLOCKS_PER_SEC );
    
    connCmp_list_t discardedCcs;
    connCmp_list_init(discardedCcs);
    
    /* main loop */
    risolate_main_loop( qResults,  qCover, discardedCcs, eps, cache, meta);
    
    connCmp_list_clear(qCover);
    realRat_clear(delta);
    compAnn_list_clear(annulii);
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

void connCmp_risolate_print_for_results(FILE * f, const connCmp_t c, metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    
    fprintf(f, "--solution with mult. %5d: ", connCmp_nSols(c));
    
    connCmp_risolate_componentBox( containingBox, c, metadatas_initBref(meta));
    
    risolate_compBox_get_containing_dsk( containingDisk, containingBox);
    
    slong d = fmpz_clog_ui(realRat_denref(compDsk_radiusref(containingDisk)),10) - fmpz_clog_ui(realRat_numref(compDsk_radiusref(containingDisk)),10); 
    slong p = fmpz_clog_ui(realRat_denref(compDsk_radiusref(containingDisk)),2) - fmpz_clog_ui(realRat_numref(compDsk_radiusref(containingDisk)),2)+50; 
    
    realApp_t cRe, rad;
    realApp_init(cRe);
    realApp_init(rad);
    
    realApp_set_realRat(cRe, compRat_realref(compDsk_centerref(containingDisk)), p);
    realApp_set_realRat(rad, compDsk_radiusref(containingDisk), p);
    
    fprintf(f, "center: ");
    realApp_fprintn(f, cRe, d, ARB_STR_NO_RADIUS);
    fprintf(f, " radius: ");
    realApp_fprintn(f, rad, 5, ARB_STR_NO_RADIUS);
    
    compBox_init(containingBox);
    compDsk_init(containingDisk);
    realApp_clear(cRe);
//     realApp_clear(cIm);
    realApp_clear(rad);
}

void connCmp_list_risolate_print_for_results(FILE * f, const connCmp_list_t l, metadatas_t meta){
    connCmp_list_iterator it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_risolate_print_for_results(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
}

void connCmp_risolate_print_for_results_withOutput(FILE * f, const connCmp_t c, int output, metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    
    if (connCmp_nSols(c) <= (10^6)-1)
        fprintf(f, "#--solution with mult. %5d: ", connCmp_nSols(c));
    else
        fprintf(f, "#--solution with mult. %5d: ", connCmp_nSols(c));
    
    connCmp_componentBox( containingBox, c, metadatas_initBref(meta));
    
    risolate_compBox_get_containing_dsk( containingDisk, containingBox);
    
    if (output == -1) { /* rational output */
        fprintf(f, "center: ");
        realRat_print(compRat_realref(compDsk_centerref(containingDisk)));
//         fprintf(f, " + ");
//         realRat_print(compRat_imagref(compDsk_centerref(containingDisk)));
        fprintf(f, "\n#%28s radius: ", " ");
        realRat_print(compDsk_radiusref(containingDisk));
    } else if (output>0) {
    
        int coeff = 4; /* = ceil(log(10)/log(2)) */
        slong prec = coeff*output; 
        
        realApp_t cRe, rad;
        realApp_init(cRe);
        realApp_init(rad);
        
        realApp_set_realRat(cRe, compRat_realref(compDsk_centerref(containingDisk)), prec);
//         realApp_set_realRat(cIm, compRat_imagref(compDsk_centerref(containingDisk)), prec);
        realApp_set_realRat(rad, compDsk_radiusref(containingDisk), prec);
        
        
    //     printf("d: %d, prec: %d\n", (int) d, (int) p);
        fprintf(f, "center: ");
        realApp_fprintn(f, cRe, output, ARB_STR_MORE);
//         fprintf(f, " + ");
//         realApp_fprintn(f, cIm, output, ARB_STR_MORE);
        fprintf(f, "\n#%28s radius: ", " ");
        realApp_fprintn(f, rad, 5, ARB_STR_MORE);
        
        realApp_clear(cRe);
        realApp_clear(rad);
    
    }
    
    compBox_clear(containingBox);
    compDsk_clear(containingDisk);
}

void connCmp_list_risolate_print_for_results_withOutput(FILE * f, const connCmp_list_t l, int output, metadatas_t meta){
    connCmp_list_iterator it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_risolate_print_for_results_withOutput(f, connCmp_list_elmt(it), output, meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
}
