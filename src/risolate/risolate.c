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
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
//         /* Real Coeffs */
//         if (( metadatas_realCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
//             compBox_clear(btemp);
//             ccluster_free(btemp);
//             continue;
//         }
//         printf("nbMSol: %d\n", (int) compBox_get_nbMSol(btemp) );
           
//             res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta);  
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
    
//     /* RealCoeffs */
//     int cc_contains_real_line = 0;
//     /* Check if cc contains the real line */
//     if ( (metadatas_realCoeffs(meta)) && (!connCmp_is_imaginary_positive(cc)) )
//         cc_contains_real_line = 1;
//     /* end RealCoeffs */
    
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
    
//     /* RealCoeffs */
//     if ( (metadatas_realCoeffs(meta)) && (connCmp_list_get_size(ltemp) == 1) && (cc_contains_real_line == 1) ){
//         ctemp = connCmp_list_first(ltemp);
//         /* test if cc has been separated from real case;
//          in which case reset everything*/
//         if ( connCmp_is_imaginary_positive(ctemp) )
//             specialFlag = 1;
//     }
//     /* end RealCoeffs */
    
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
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
      
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        
        sepBoundFlag   = (realRat_cmp( compBox_bwidthref(componentBox), metadatas_getSepBound(meta))<=0);
        
        if (metadatas_getVerbo(meta)>3) {
            printf("---depth: %d\n", (int) depth);
//             printf("------component Box:       "); compBox_print(componentBox); printf("\n");
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
            printf("------nb of boxes:         %d\n", compBox_list_get_size(connCmp_boxesref(ccur)));
            printf("------nb of roots:         %d\n", connCmp_nSolsref(ccur));
            printf("------last newton success: %d\n", connCmp_newSuref(ccur));
            printf("------newton speed       : "); fmpz_print(connCmp_nwSpdref(ccur)); printf("\n");
            printf("------separation Flag:     %d\n", separationFlag);
            printf("------widthFlag:           %d\n", widthFlag); 
            printf("------compactFlag:         %d\n", compactFlag);
            printf("------sepBoundFlag:        %d\n", sepBoundFlag);
            printf("------current prec:        %ld\n",connCmp_appPrref(ccur));
//             compBox_list_print( connCmp_boxesref(ccur) );
//             printf("\n");
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
//         if ((separationFlag)) {
//             printf("depth: %d, connCmp_nSolsref(ccur): %d, prec: %d\n", (int) depth, (int) connCmp_nSolsref(ccur), (int) prec);
//             
//             if (connCmp_nSolsref(ccur)==-1){
                
//                 resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                resTstar = tstar_real_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
                if (metadatas_getVerbo(meta)>3)
                    printf("------run tstar: nbSols: %d, required precision: %ld\n", (int) connCmp_nSolsref(ccur), resTstar.appPrec);
//                 ???
                prec = resTstar.appPrec;
//             }
                
//             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
//             ???
//             prec = resTstar.appPrec;
                
//                 if ( ( (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && (!widthFlag) )
//                      ||( connCmp_nSols(ccur)== cacheApp_getDegree(cache) )
//                      ||( !sepBoundFlag && (connCmp_nSols(ccur)>1))   )
//                 {
//                     
//                         compDsk_t nDisk;
//                         compDsk_init(nDisk);
//                         realRat_t factor, nwidth;
//                         realRat_init(factor);
//                         realRat_init(nwidth);
//                         printf("Try to DOWNSIZE!!!\n");
//                         int resDownsize = tstar_real_downsize( nDisk, factor, cache, ccDisk, (int) connCmp_nSolsref(ccur), prec, depth, meta );
//                         
//                         if (resDownsize==1) {
//                             
//                             printf("res: %d, old disk: ", resDownsize);
//                             compDsk_print(ccDisk);
//                             printf(" new disk: ");
//                             compDsk_print(nDisk);
//                             printf("\n");
//                             printf("ratio: "); realRat_print(factor); printf("\n"); 
//                             
//                             connCmp_ptr nCC;
//                             nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//                             connCmp_init(nCC);
//             
//                             realRat_mul(nwidth, connCmp_widthref(ccur), factor);
//                             
//                             compBox_list_ptr ltemp;
//                             compBox_ptr btemp;
//                             ltemp = connCmp_boxesref(ccur);
//                             
//                             compBox_list_t ltemp2;
//                             compBox_list_init(ltemp2);
//                             
//                             while (compBox_list_get_size(ltemp)>0){
//                                 btemp = compBox_list_pop(ltemp);
//                                 subdBox_risolate_bisect_with_compDsk( ltemp2, btemp, nDisk, nwidth);
//                                 compBox_clear(btemp);
//                                 ccluster_free(btemp); /*comment it for julia...*/
//                             }
//                             compBox_list_swap(ltemp, ltemp2);
//                             compBox_list_clear(ltemp2);
//                             
// //                             printf("length of list: %d\n", compBox_list_get_size(ltemp));
// //                             compBox_list_print(ltemp); printf("\n"); 
//                             
//                             btemp = compBox_list_pop(ltemp);
//                             realRat_set(connCmp_widthref(nCC), compBox_bwidthref(btemp));
//                             connCmp_insert_compBox(nCC, btemp);
//                             while (!compBox_list_is_empty(ltemp))
//                                 connCmp_insert_compBox(nCC, compBox_list_pop(ltemp));
//                             connCmp_nSols(nCC) = connCmp_nSols(ccur);
//                             fmpz_set(connCmp_nwSpdref(nCC), connCmp_nwSpdref(ccur));
//                             /* test */
//                             connCmp_isSep(nCC) = connCmp_isSep(ccur);
//                             
//                             connCmp_clear(ccur);
//                             ccluster_free(ccur);
//                             ccur = nCC;
//                             connCmp_newSuref(ccur) = 0;
//                             connCmp_appPrref(ccur) = prec;
//                             
//                             connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
//                             compBox_get_containing_dsk(ccDisk, componentBox);
//                             realRat_mul(threeWidth, three, connCmp_widthref(ccur));
//                             depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
//         
//       
//                             widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
//                             compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
//         
//                             sepBoundFlag   = (realRat_cmp( compBox_bwidthref(componentBox), metadatas_getSepBound(meta))<=0);
//                         }
//     
//                         compDsk_clear(nDisk);
//                         realRat_clear(factor);
//                         realRat_clear(nwidth);
//                 }
        }
        
        /* special case where zero is a root with mult>1 */
        /* and current cc is separated and contains zero */
        if ( (separationFlag) && 
             (connCmp_nSolsref(ccur) > 1 ) &&
             (realRat_sgn(connCmp_infReref(ccur)) < 0) &&
             (realRat_sgn(connCmp_supReref(ccur)) > 0) ) {
            if ( connCmp_nSolsref(ccur) == cacheApp_getMultOfZero( cache ) )
                sepBoundFlag = 1;
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

            if (metadatas_getVerbo(meta)>3)
                    printf("------run Newton: res: %d, required precision: %ld\n", resNewton.nflag, resNewton.appPrec);
            
//             printf("+++depth: %ld, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
            if (resNewton.nflag) {
                
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = nCC;
                connCmp_newSuref(ccur) = 1;
                connCmp_appPrref(ccur) = resNewton.appPrec;
                
                connCmp_risolate_componentBox(componentBox, ccur, metadatas_initBref(meta));
                compBox_get_containing_dsk(ccDisk, componentBox);
                
                if (metadatas_useDeflation(meta))
                if (metadatas_getVerbo(meta)>1) {
                    if (connCmp_nSolsref(ccur) > 1 ) {
//                         if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)==0) {
                        if (connCmp_isDefref(ccur)==0) {
                            if (realRat_cmp_ui(compDsk_radiusref(ccDisk), 1 ) < 0) {
                                
                                printf("\n\n\n ------Success of Newton Iteration for this Component with a cluster of %d roots------\n", connCmp_nSolsref(ccur) );
//                                 printf("------Compute new deflation, newton speed: ");
//                                 fmpz_print(connCmp_nwSpdref(ccur));
//                                 printf("\n");
                                deflate_connCmp_init(ccur);
                                deflate_set( ccur, cache, ccDisk, connCmp_nSolsref(ccur), connCmp_appPrref(ccur), meta );
                                printf("precision: %ld\n", connCmp_appPrref(ccur) );
                                printf("------\n\n\n");
                            }
                        }
                    }
                }
                
                connCmp_increase_nwSpd(ccur);
    
            }
            else {
                /* newton has been run but wasn't successful */
                connCmp_newSuref(ccur) = 2;
                connCmp_clear(nCC);
                ccluster_free(nCC);
            }
            if (metadatas_haveToCount(meta)){
                metadatas_add_Newton   ( meta, depth, resNewton.nflag, (double) (clock() - start) );
            }
        }
        
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
            if (metadatas_getVerbo(meta)>3) {
                printf("------bisect\n" );
            }
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

void connCmp_risolate_print_for_results(FILE * f, const connCmp_t c, metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    
    fprintf(f, "--solution with mult. %5d: ", connCmp_nSols(c));
    
    connCmp_risolate_componentBox( containingBox, c, metadatas_initBref(meta));
    compBox_get_containing_dsk( containingDisk, containingBox);
    
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
    compBox_get_containing_dsk( containingDisk, containingBox);
    
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
