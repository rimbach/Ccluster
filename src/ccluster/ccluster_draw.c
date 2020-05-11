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

#include "ccluster/ccluster.h"

#include <pthread.h>

slong ccluster_discard_compBox_list_draw( compBox_list_t boxes, compBox_list_t discarded, 
                                          cacheApp_t cache, 
//                                      int nbSols, 
                                     slong prec, metadatas_t meta){
    
//     printf("ccluster_draw.c: ccluster_discard_compBox_list_draw: \n");
    
    tstar_res res;
    res.appPrec = prec;
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
    compBox_ptr btemp;
    
    /* For test */
//     int nbSolsAlreadyfound = 0;
//     compBox_list_t ltempDetermined;
//     compBox_list_init(ltempDetermined);
    /* End For test */
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
        /* Real Coeffs */
        if (( metadatas_useRealCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
            compBox_list_push(discarded, btemp);
//             compBox_clear(btemp);
//             ccluster_free(btemp);
            continue;
        }
        
        res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta); 
        if (res.nbOfSol==0) {
            if (metadatas_haveToCount(meta)){
                metadatas_add_discarded( meta, depth);
            }
            compBox_list_push(discarded, btemp);
//             compBox_clear(btemp);
//             ccluster_free(btemp);
        }
        
        else{
            
//             if (metadatas_useCountSols(meta)&&(nbSols>0)) {                                                                                    
//                 if (res.nbOfSol==-1) {                                                                                                         
//                     compBox_list_push(ltemp, btemp);                                                                                           
//                 }                                                                                                                              
//                 else { /*res.nbOfSol>0*/                                                                                                       
//                     btemp->nbMSol = res.nbOfSol;                                                                                               
//                     /* try to validate that btemp contains btemp->nbMSol roots */                                                              
//                     res = tstar_count_nb_Sols(cache, bdisk, res.nbOfSol, res.appPrec, depth, meta);                                            
//                     if (res.nbOfSol==1){ /* btemp contains btemp->nbMSol roots */                                                              
//                         nbSolsAlreadyfound+=btemp->nbMSol;                                                                                     
//                         compBox_list_push(ltempDetermined, btemp);                                                                             
//                         if (nbSolsAlreadyfound==nbSols) {                                                                                      
//                              printf("&&&&&&&& la &&&&&&&&&& depth: %d nb of boxes still in list: %d\n", depth, compBox_list_get_size(boxes));  
//                             /*empty boxes*/                                                                                                    
//                             while (!compBox_list_is_empty(boxes)) {                                                                            
//                                 btemp = compBox_list_pop(boxes);                                                                               
//                                 metadatas_add_discarded( meta, depth);                                                                         
//                                 compBox_clear(btemp);                                                                                          
//                                 free(btemp);                                                                                                   
//                             }                                                                                                                  
//                             /*empty ltemp*/                                                                                                    
//                             while (!compBox_list_is_empty(ltemp)) {                                                                            
//                                 btemp = compBox_list_pop(ltemp);                                                                               
//                                 metadatas_add_discarded( meta, depth);                                                                         
//                                 compBox_clear(btemp);                                                                                          
//                                 free(btemp);                                                                                                   
//                             }                                                                                                                  
//                         }                                                                                                                      
//                     }                                                                                                                          
//                     else {                                                                                                                     
//                         compBox_list_push(ltemp, btemp);                                                                                       
//                     }                                                                                                                          
//                 }                                                                                                                              
//             }                                                                                                                                  
            
//             else {
                if (res.nbOfSol>0)
                    btemp->nbMSol = res.nbOfSol;
                compBox_list_push(ltemp, btemp);
//                 compBox_list_insert_sorted(ltemp, btemp);
//             }
        }
    }
/*    
    if (metadatas_useCountSols(meta)&&(nbSols>0)){
        while (!compBox_list_is_empty(ltempDetermined)) 
            compBox_list_push(ltemp, compBox_list_pop(ltempDetermined));
    }
*/    
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
    /* For test */
//     compBox_list_clear(ltempDetermined);
    /* End For test */
    
    return res.appPrec;
}

void ccluster_bisect_connCmp_draw( connCmp_list_t dest, connCmp_t cc, connCmp_list_t discardedCcs, compBox_list_t discarded, cacheApp_t cache, metadatas_t meta, slong nbThreads){
    
//     printf("ccluster_draw.c: ccluster_bisect_connCmp_draw: \n");
    
    slong prec = connCmp_appPr(cc);
    compBox_list_t subBoxes;
    connCmp_list_t ltemp;
    compBox_list_init(subBoxes);
    connCmp_list_init(ltemp);
    
    compBox_ptr btemp;
    connCmp_ptr ctemp;
    
    /* RealCoeffs */
    int cc_contains_real_line = 0;
    /* Check if cc contains the real line */
    if ( (metadatas_useRealCoeffs(meta)) && (!connCmp_is_imaginary_positive(cc)) )
        cc_contains_real_line = 1;
    /* end RealCoeffs */
    
    while (!connCmp_is_empty(cc)) {
        btemp = connCmp_pop(cc);
        subdBox_quadrisect( subBoxes, btemp, 0 );
        compBox_clear(btemp);
        ccluster_free(btemp);
    }
    
    prec = ccluster_discard_compBox_list_draw( subBoxes, discarded, cache, prec, meta);
    
    while (!compBox_list_is_empty(subBoxes)) {
        btemp = compBox_list_pop(subBoxes);
        connCmp_union_compBox( ltemp, btemp);
    }
    int specialFlag = 1;
    if (connCmp_list_get_size(ltemp) == 1)
        specialFlag = 0;
    
    /* RealCoeffs */
    if ( (metadatas_useRealCoeffs(meta)) && (connCmp_list_get_size(ltemp) == 1) && (cc_contains_real_line == 1) ){
        ctemp = connCmp_list_first(ltemp);
        /* test if cc has been separated from real case;
         in which case reset everything*/
        if ( connCmp_is_imaginary_positive(ctemp) )
            specialFlag = 1;
    }
    /* end RealCoeffs */
    
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

void ccluster_prep_loop_draw( connCmp_list_t qMainLoop, connCmp_list_t qPrepLoop, connCmp_list_t discardedCcs, compBox_list_t discarded, cacheApp_t cache, metadatas_t meta) {
    
//     printf("ccluster_draw.c: ccluster_prep_loop_draw: \n");
    
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
            ccluster_bisect_connCmp_draw( ltemp, ctemp, discardedCcs, discarded, cache, meta, metadatas_useNBThreads(meta));
//             ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, cache, meta,1);
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

// int  ccluster_compDsk_is_separated( const compDsk_t d, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs ){
//     int res = 1;
//     connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
//     while ( res && (it!=connCmp_list_end()) ) {
//         res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
//         it = connCmp_list_next(it);
//     }
//     it = connCmp_list_begin(discardedCcs);
//     while ( res && (it!=connCmp_list_end()) ) {
//         res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
//         it = connCmp_list_next(it);
//     }
//     return res;
// }
    
void ccluster_main_loop_draw( connCmp_list_t qResults,  connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, compBox_list_t discarded, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
    int separationFlag;
    int widthFlag;
    int compactFlag;
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
    
    /* Real Coeff */
    connCmp_ptr ccurConjClo, ccurConj;
    ccurConjClo = NULL;
    ccurConj = NULL;
    int pushConjugFlag = 0;
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    
    while (!connCmp_list_is_empty(qMainLoop)) {
        
        resNewton.nflag = 0;
        
        ccur = connCmp_list_pop(qMainLoop);
        
        /* Real Coeff */
        pushConjugFlag = 0;
        if (metadatas_useRealCoeffs(meta)){
            /* test if the component contains the real line in its interior */
//             printf("number of boxes before conjugate closure: %d\n", connCmp_nb_boxes(ccur));
            if (!connCmp_is_imaginary_positive(ccur)) {
//                 printf("number of boxes before conjugate closure: %d\n", connCmp_nb_boxes(ccur));
                ccurConjClo = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
                connCmp_init( ccurConjClo );
                connCmp_set_conjugate_closure(ccurConjClo, ccur, metadatas_initBref(meta));
                
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = ccurConjClo;
//                 printf("number of boxes after  conjugate closure: %d\n", connCmp_nb_boxes(ccur));
            }
        }
        
        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
        compBox_get_containing_dsk(ccDisk, componentBox);
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
        
        /* Real Coeff */
        if ( (separationFlag)&&(metadatas_useRealCoeffs(meta)) ) {
            if (connCmp_is_imaginary_positive(ccur)) {
                /* check if ccur is separated from its complex conjugate */
                realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
                separationFlag = separationFlag&&(!compBox_intersection_is_not_empty_compDsk ( componentBox, fourCCDisk));
                realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
            }
        }
        
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
//         if ((separationFlag)) {
//             printf("depth: %d, connCmp_nSolsref(ccur): %d\n", depth, connCmp_nSolsref(ccur));
            if (connCmp_nSolsref(ccur)==-1){
                resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
            }
//             printf("validate: prec avant: %d prec apres: %d\n", prec, resTstar.appPrec);
            prec = resTstar.appPrec;
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
               ( (!widthFlag)||( connCmp_nSols(ccur)== cacheApp_getDegree(cache) ) )  )
//             &&!( metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) ) //this is DEPRECATED: pass eps = 1/0 instead
        ) {
            
            if (metadatas_haveToCount(meta)){
                start = clock();
            }
            
            if (connCmp_nSols(ccur)==1) 
                compRat_set(initPoint, compBox_centerref(componentBox));
            else
                connCmp_find_point_outside_connCmp( initPoint, ccur, metadatas_initBref(meta) );
            
            connCmp_ptr nCC;
            nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
            connCmp_init(nCC);
            resNewton = newton_newton_connCmp( nCC, ccur, cache, initPoint, prec, meta);
//             metadatas_add_Newton   ( meta, depth, resNewton.nflag );
//             printf("+++depth: %d, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
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
        
        /* Real Coeff */
        if (metadatas_useRealCoeffs(meta)
            && ( 
//             (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag)
//                ||//this is DEPRECATED: pass eps = 1/0 instead 
               ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) ) ) {
            
            if (connCmp_is_imaginary_positive(ccur)) {
                pushConjugFlag = 1;
                /*compute the complex conjugate*/
                ccurConj = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
                connCmp_init( ccurConj );
                connCmp_set_conjugate(ccurConj, ccur);
                
                /* test if initial box is symetric relatively to real axe */
                if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
                    /* test if the cc intersects initial box */
                    if ( connCmp_intersection_is_not_empty(ccurConj, metadatas_initBref(meta)) ) {
                        /* test if the cc is confined */
//                         if (connCmp_is_confined(ccurConj, metadatas_initBref(meta))) {
//                             pushConjugFlag = 1;
//                         }
//                         else {
                        if (!connCmp_is_confined(ccurConj, metadatas_initBref(meta))) {
                            pushConjugFlag = 0;
                            separationFlag = 0; 
                          /* delete ccurConj*/
                            connCmp_clear(ccurConj);
                            ccluster_free(ccurConj);
                        }
                    }
                    else {
                        pushConjugFlag = 0;
                        /* delete ccurConj*/
                        connCmp_clear(ccurConj);
                        ccluster_free(ccurConj);
                    }
                } 
            }
            else {
                /* test if initial box is symetric relatively to real axe */
                if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
                    /* test if the cc is confined and intersects initial box */
                    if (! ( connCmp_is_confined(ccur, metadatas_initBref(meta)) 
                         && connCmp_intersection_is_not_empty(ccur, metadatas_initBref(meta)) ) ){
//                         /* bisect ccur until this hold */
                        separationFlag = 0;
                    }
                }
            }
        }
        
//         if (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag){
//             metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
//             connCmp_list_push(qResults, ccur);
// //             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
//             /* Real Coeff */
//             if ((metadatas_useRealCoeffs(meta))&&(pushConjugFlag)){
//                 /*compute the complex conjugate*/
//                 metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
//                 connCmp_list_push(qResults, ccurConj);
//             }
//         }
//         else
        if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag && (connCmp_nSols(ccur)<cacheApp_getDegree(cache)) ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
//             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
            /* Real Coeff */
            if ((metadatas_useRealCoeffs(meta))&&(pushConjugFlag)){
                /*compute the complex conjugate*/
                metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
                connCmp_list_push(qResults, ccurConj);
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
#ifdef CCLUSTER_HAVE_PTHREAD
                ccluster_bisect_connCmp_draw( ltemp, ccur, discardedCcs, discarded, cache, meta,metadatas_useNBThreads(meta));
                while (!connCmp_list_is_empty(ltemp))
                    connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
                connCmp_clear(ccur);
                ccluster_free(ccur);
#else
            ccluster_bisect_connCmp_draw( ltemp, ccur, discardedCcs, discarded, cache, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
#endif
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

void ccluster_algo_draw( connCmp_list_t qResults, 
                         compBox_list_t discarded, 
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         cacheApp_t cache, metadatas_t meta){
    
//     chronos_tic_CclusAl(metadatas_chronref(meta));
    clock_t start = clock();
    
//     printf("ccluster_draw.c: ccluster_algo_draw \n");
    
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
//     printf("ccluster_draw.c: ccluster_algo_draw, preploop: \n");
    ccluster_prep_loop_draw( qMainLoop, qPrepLoop, discardedCcs, discarded, cache, meta);
//     printf("ccluster_draw.c: ccluster_algo_draw, mainloop: \n");
    ccluster_main_loop_draw( qResults,  qMainLoop, discardedCcs, discarded, eps, cache, meta);
    
    
    realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}




