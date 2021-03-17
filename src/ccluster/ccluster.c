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

// #include <pthread.h>

slong ccluster_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
                                     cacheApp_t cache, 
//                                      int nbSols, 
                                     slong prec, metadatas_t meta){
    
//     clock_t start = clock();
    
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
    
    slong degree = cacheApp_getDegree(cache);
    slong dtemp = 1;
    int log2d = 0;
    while (degree>=dtemp){
        dtemp = 2*dtemp;
        log2d++;
    }
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
        /* Real Coeffs */
        if (( metadatas_useRealCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
//             if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
//             } else {
//                 compBox_list_push(bDiscarded, btemp);
//             }
            continue;
        }
        
        /* Root Radii */
        if ( metadatas_useRootRadii(meta) ){
/* returns 1 => 2*box contains at least one root */
/*         0 => box contains no root */
/*        -1 => can not decide */
            int RRExcl = ccluster_rootRadii_exclusion_test( btemp, CCLUSTER_DEFAULT_PREC, meta );
            if (RRExcl==0) {
                if (metadatas_haveToCount(meta)){
                    metadatas_add_discarded( meta, depth);
                }
//                 if (metadatas_getDrSub(meta)==0){
                    compBox_clear(btemp);
                    ccluster_free(btemp);
//                 } else {
//                     compBox_list_push(bDiscarded, btemp);
//                 }
                
//                 printf("# do not intersect at least one annulii of each group\n");
                
                continue;
            } else if (RRExcl==1) {
                compBox_list_push(ltemp, btemp);
                continue;
            }
        }
        
        if ( metadatas_usePowerSums(meta) ){
        
            powerSums_res resp;
            resp.appPrec = CCLUSTER_DEFAULT_PREC;
            
            /* resp = powerSums_discardingTest( compDsk_centerref(bdisk), compDsk_radiusref(bdisk),
                                                        cache,
                                                        metadatas_getNbEvalPoints(meta),
                                                        metadatas_getNbPowerSums(meta),
                                                        resp.appPrec, meta, depth );
            
            res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1,0, res.appPrec, depth, meta);
            metadatas_add_PsCountingTest (meta, depth, resp.nbOfSol, res.nbOfSol); */
            
            resp = powerSums_countingTest( compDsk_centerref(bdisk), compDsk_radiusref(bdisk),
                                                        cache,
                                                        metadatas_getNbEvalPoints(meta),
                                                        0,
                                                        resp.appPrec, meta, depth );
            metadatas_add_PsCountingTest (meta, depth);
            if ((resp.nbOfSol==0)||(resp.nbOfSol==-2)) {
                res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1,0, res.appPrec, depth, meta); 
            }
            else
                res.nbOfSol = -1;
        }
        else {
  
            res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta);  
//             printf("## tstar result: %d \n", res.nbOfSol);
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
//                                 ccluster_free(btemp);                                                                                                   
//                             }                                                                                                                  
//                             /*empty ltemp*/                                                                                                    
//                             while (!compBox_list_is_empty(ltemp)) {                                                                            
//                                 btemp = compBox_list_pop(ltemp);                                                                               
//                                 metadatas_add_discarded( meta, depth);                                                                         
//                                 compBox_clear(btemp);                                                                                          
//                                 ccluster_free(btemp);                                                                                                   
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

void ccluster_bisect_connCmp( connCmp_list_t dest, 
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
    
    /* RealCoeffs */
    int cc_contains_real_line = 0;
    /* Check if cc contains the real line */
    if ( (metadatas_useRealCoeffs(meta)) && (!connCmp_is_imaginary_positive(cc)) )
        cc_contains_real_line = 1;
    /* end RealCoeffs */
    
    while (!connCmp_is_empty(cc)) {
        btemp = connCmp_pop(cc);
        subdBox_quadrisect( subBoxes, btemp );
        compBox_clear(btemp);
        ccluster_free(btemp);
    }
    
#ifdef CCLUSTER_HAVE_PTHREAD
    if (nbThreads>1) {
        prec = ccluster_parallel_discard_compBox_list( subBoxes, cache, prec, meta, nbThreads);
    }
    else {
        prec = ccluster_discard_compBox_list( subBoxes, bDiscarded, cache, prec, meta);
    }
#else
    prec = ccluster_discard_compBox_list( subBoxes, bDiscarded, cache, prec, meta);
#endif
    
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

void ccluster_prep_loop( compBox_list_t bDiscarded,
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
            ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
            
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

int  ccluster_compDsk_is_separated( const compDsk_t d, connCmp_list_t qMainLoop, connCmp_list_t discardedCcs ){
    int res = 1;
    connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    it = connCmp_list_begin(discardedCcs);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    return res;
}
                         
void ccluster_main_loop( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta){
    
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
    
    realRat_t epsSave;
    realRat_init( epsSave );
    realRat_set_si( epsSave, 2,1  );
    realRat_pow_si( epsSave, epsSave, -53);
    int widthFlagSave;
    
    clock_t start=clock();
    
    /* Real Coeff */
    connCmp_ptr ccurConjClo, ccurConj;
    ccurConjClo = NULL;
    ccurConj = NULL;
    int pushConjugFlag = 0;
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    
    while (!connCmp_list_is_empty(qMainLoop)) {
        
//         if (metadatas_getVerbo(meta)>=2) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }

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
        
        widthFlagSave      = (realRat_cmp( compBox_bwidthref(componentBox), epsSave)<=0);
        
        if (metadatas_getVerbo(meta)>3) {
            printf("---depth: %d\n", (int) depth);
            printf("------component Box:"); compBox_print(componentBox); printf("\n");
            printf("------separation Flag: %d\n", separationFlag);
            printf("------widthFlag: %d\n", widthFlag); 
            printf("------compactFlag: %d\n", compactFlag); 
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
//         if ((separationFlag)) {
//             printf("depth: %d, connCmp_nSolsref(ccur): %d, prec: %d\n", (int) depth, (int) connCmp_nSolsref(ccur), (int) prec);
            if (connCmp_nSolsref(ccur)==-1){
                
                if (metadatas_usePowerSums(meta)) {
                    realRat_t temp;
                    realRat_init(temp);
                    realRat_set_si(temp, 2, 1);
                    realRat_mul(temp, compDsk_radiusref(ccDisk), temp);
                    
                    powerSums_res resp;
                    
                    resp = powerSums_countingTest( compDsk_centerref(ccDisk), temp,
                                                        cache,
                                                        metadatas_getNbEvalPoints(meta), 
                                                        1,
                                                        prec, meta, depth );
                    
                    connCmp_nSolsref(ccur) = resp.nbOfSol;
                    prec = resp.appPrec;     
                    realRat_clear(temp);
                }    
                else {
                    resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                    connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                     if (metadatas_getVerbo(meta)>3)
//                         printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
//                 ???
                    prec = resTstar.appPrec;
                }
                
            }
//             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
//             ???
//             prec = resTstar.appPrec;
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
               ( (!widthFlag) || 
                 ( (widthFlag)&&
                   (connCmp_nSols(ccur)>1)&&
                   (connCmp_nSols(ccur)== cacheApp_getDegree(cache) )&&
                   (!widtFlagSave) )  
               ) )
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
//             (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag) //this is DEPRECATED: pass eps = 1/0 instead 
//                ||
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
        //this is DEPRECATED: pass eps = 1/0 instead 
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
        if ( (connCmp_nSols(ccur)>0) && 
             separationFlag && 
             widthFlag && 
             compactFlag && 
             ( ( connCmp_nSols(ccur)==1 ) || 
               (connCmp_nSols(ccur)<cacheApp_getDegree(cache)) ||
               ( (connCmp_nSols(ccur)==cacheApp_getDegree(cache)) && widthFlagSave)  
            ) ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
//             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
//             printf("metadatas_useRealCoeffs(meta): %d, pushConjugFlag: %d\n", metadatas_useRealCoeffs(meta), pushConjugFlag);
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
            if (connCmp_nSols(ccur)==0) { /* can occur only when using root radii */
                connCmp_clear(ccur);
                ccluster_free(ccur);
                continue;
            }
#ifdef CCLUSTER_HAVE_PTHREAD
            ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cache, meta, metadatas_useNBThreads(meta));
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
#else
            ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cache, meta,1);
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
    
    realRat_init( epsSave );
}

void ccluster_algo( connCmp_list_t qResults, 
                    compBox_list_t bDiscarded, 
                    const compBox_t initialBox, 
                    const realRat_t eps, 
                    cacheApp_t cache, 
                    metadatas_t meta){
    
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
    ccluster_prep_loop( bDiscarded, qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     if (metadatas_getVerbo(meta)>3) printf("Ccluster mainloop: \n");
    ccluster_main_loop( qResults, bDiscarded, qMainLoop, discardedCcs, eps, cache, meta);
    
    
    realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

void ccluster_algo_global( connCmp_list_t qResults, 
                           compBox_list_t bDiscarded,
                           const compBox_t initialBox, 
                           const realRat_t eps, 
                           cacheApp_t cache, 
                           metadatas_t meta){
    
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
    
    ccluster_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cache, meta);
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

// void ccluster_algo_global_rootRadii( connCmp_list_t qResults,
//                                      compBox_list_t bDiscarded,
//                                      compAnn_list_t annulii,
//                                      compAnn_list_t annulii1,
//                                      compAnn_list_t annulii2,
//                                      const compBox_t initialBox, 
//                                      const realRat_t eps, 
//                                      cacheApp_t cache, 
//                                      metadatas_t meta){
//     clock_t start = clock();
//     clock_t start2 = clock();
//     
//     realRat_t delta;
//     realRat_init(delta);
//     
//     slong degree = cacheApp_getDegree( cache );
// //     realRat_set_si(delta, 1, degree);
//     realRat_set_si(delta, 1, degree*degree);
// //     realRat_set_si(delta, 1, degree*degree*degree);
//     
//     slong bitsize = cacheApp_getBitsize (cache);
//     if (metadatas_getVerbo(meta)>=2) {
//         printf("#bitsize of input polynomial: %ld\n", bitsize);
//     }
//     
//     /* assume 0 is not a root of the poly */
//     slong prec = CCLUSTER_DEFAULT_PREC;
//     /* heuristic to predict the precision: at least the degree */
//     while (prec<degree/2)
//         prec = 2*prec;
//     /* */
//     slong prec1 = realIntRootRadii_rootRadii( annulii, 0, cache, delta, prec, meta );
//     
//     /* derive upper bound for the norm of the roots */
//     slong centerRe1 = realApp_ceil_si(compAnn_radSupref(compAnn_list_last(annulii)), prec);
//     prec = 2*prec1;
//     slong prec2 = realIntRootRadii_rootRadii( annulii1, 1, cache, delta, prec, meta );
//     slong prec3 = realIntRootRadii_rootRadii_imagCenter( annulii2, 1, cache, delta, prec, meta );
//     
//     prec3 = CCLUSTER_MAX( prec2, prec3);
//     prec3 = CCLUSTER_MAX( prec1, prec3);
//     
// //     if (metadatas_getVerbo(meta)>3) {
// //         printf("#time in computing RootRadii           : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
// //         printf("#precision needed: %ld\n", prec3 );
// //     }
// //     start2 = clock();
//     
//     realIntRootRadii_connectedComponents( annulii, prec );
//     realIntRootRadii_connectedComponents( annulii1, prec );
//     realIntRootRadii_connectedComponents( annulii2, prec );
//     
// //     compAnn_list_empty(annulii2);
//     
// //     if (metadatas_getVerbo(meta)>3) {
// //         printf("#time in computing connected components: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
// //     }
//         
// //     start2 = clock();
//     
//     realIntRootRadii_containsRealRoot( annulii, cache, prec );
// //     if (metadatas_getVerbo(meta)>=2) {
// //         printf("#time in discarding real intervals     : %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#Annulii: ");
//             compAnn_list_printd(annulii, 10);
//             printf("\n\n");
//         }
// //     }
//     
//     start2 = clock();
//     if (metadatas_haveToCount(meta))
//         metadatas_add_time_rootRad(meta, (double) (start2 - start) );
//     
// //     printf("#Annulii 0: ");
// //     compAnn_list_printd(annulii, 10);
// //     printf("\n\n");
// // //     
// //     printf("#Annulii 1: ");
// //     compAnn_list_printd(annulii1, 10);
// //     printf("\n\n");
// //     
// //     printf("#Annulii 2: ");
// //     compAnn_list_printd(annulii2, 10);
// //     printf("\n\n");
//     
//     /* test: compute all the intersections */
// //     int nbInt = 0;
// //     compApp_t inter, interConj;
// //     compApp_init(inter);
// //     compApp_init(interConj);
// //     compAnn_list_iterator it = compAnn_list_begin( annulii );
// //     while ( it!=compAnn_list_end() ) {
// //         
// //         int intersect = 0;
// //         int nbInterA0 = 0;
// //         /* find first annulus of annulii1 intersecting that annulus */
// //         compAnn_list_iterator it1 = compAnn_list_begin( annulii1 );
// //         while ( it1!=compAnn_list_end() ) {
// //             
// //             intersect = compAnn_intersect_realCenter( inter, compAnn_list_elmt( it ), compAnn_list_elmt( it1 ),
// //                                                                        CCLUSTER_DEFAULT_PREC);
// //             if (intersect) {
// //                 int containsRealLine = realApp_contains_zero ( compApp_imagref( inter ) );
// //                 int nbSolsInIt = compAnn_indMaxref(compAnn_list_elmt( it )) - compAnn_indMinref(compAnn_list_elmt( it )) + 1;
// //                 int nbSolsInIt1 = compAnn_indMaxref(compAnn_list_elmt( it1 )) - compAnn_indMinref(compAnn_list_elmt( it1 )) + 1;
// //                 if ( (containsRealLine==0) && ( (nbSolsInIt == 1) || (nbSolsInIt1 == 1) ) ) {
// //                     it1 = compAnn_list_next(it1);
// //                     continue;
// //                 }
// //                 
// //                 /* check if the intersection intersects at least one annulus of the third group */
// //                 /*   and if the conjugate of the intersection intersects at least one annulus of the third group */
// //                 
// //                 /* find first annulus of annulii2 intersecting that annulus of annulii1 */
// //                 compAnn_list_iterator it2 = compAnn_list_begin( annulii2 );
// //                 
// //                 intersect = 0;
// //                 int intersectConj = 1;
// //                 compApp_set(interConj, inter);
// //                 realApp_neg(compApp_imagref(interConj), compApp_imagref(interConj));
// //                 
// //                 while ( (it2!=compAnn_list_end() )&& ( (intersect==0) || (intersectConj==0) ) ) {
// //                     if (intersect == 0)
// //                         intersect = ccluster_is_compApp_in_compAnn( inter, compAnn_list_elmt( it2 ), CCLUSTER_DEFAULT_PREC);
// //                     if (intersectConj == 0)
// //                         intersectConj = ccluster_is_compApp_in_compAnn( interConj, compAnn_list_elmt( it2 ), CCLUSTER_DEFAULT_PREC);
// //                     it2 = compAnn_list_next(it2);
// //                 }
// //                 if (intersect && intersectConj) {
// //                     nbInt ++;
// //                     nbInterA0++;
// //                 }
// //             }
// //             it1 = compAnn_list_next(it1);
// //         }
// //         printf("### number of intersections: %d\n", nbInterA0);
// //         it = compAnn_list_next(it);
// //     }
// //     
// //     compApp_clear(inter);
// //     compApp_clear(interConj);
// //     
// //     printf("# number of intersections: %d\n", nbInt);
// //     printf("#time in computing intersections: %f \n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
//     /* fin test */
//     
//     compBox_ptr box;
//     box = (compBox_ptr) ccluster_malloc (sizeof(compBox));
//     compBox_init(box);
//     compBox_set(box, initialBox);
//     compBox_nbMSolref(box) = cacheApp_getDegree ( cache );
//     /* set width to 2*centerRe1 */
//     realRat_set_si( compBox_bwidthref(box), 2*centerRe1, 1 );
//     compBox_set(metadatas_initBref(meta), box);
//     
//     compBox_copy_annulii(box, 0, annulii);
//     compBox_copy_annulii(box, 1, annulii1);
//     compBox_copy_annulii(box, 2, annulii2);
//     /* fourth list of annulii contains annulii from i having intersection with the */
//     /* complex conjugate of the box */
//     compBox_copy_annulii(box, 3, annulii2);
//     
//     connCmp_ptr initialCC;
//     initialCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//     connCmp_init_compBox(initialCC, box);
//     
//     connCmp_list_t qPrepLoop;
//     connCmp_list_init(qPrepLoop);
//     
//     connCmp_list_t qMainLoop, discardedCcs;
//     connCmp_list_init(qMainLoop);
//     connCmp_list_init(discardedCcs);
//     
//     connCmp_list_push(qPrepLoop, initialCC);
//     
//     connCmp_list_push(qMainLoop, initialCC);
//     ccluster_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cache, meta);
//     
//     
//     connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(discardedCcs);
//     
//     realRat_clear(delta);
//     
//     metadatas_add_time_CclusAl(meta, (double) (clock() - start));
// }

void ccluster_refine( connCmp_list_t qResults, 
                      connCmp_list_t qMainLoop,
//                       const compBox_t initialBox, 
                      const realRat_t eps, 
                      cacheApp_t cache, 
                      metadatas_t meta){
    
//     printf("ccluster_refine: begin, eps = "); realRat_print(eps); printf("\n");
    
//     chronos_tic_CclusAl(metadatas_chronref(meta));
    clock_t start = clock();
    
    
    connCmp_list_t discardedCcs;
    connCmp_list_init(discardedCcs);

    ccluster_main_loop( qResults,  NULL, qMainLoop, discardedCcs, eps, cache, meta);
    
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
//     printf("ccluster_refine: end \n");
}

void connCmp_print_for_results(FILE * f, const connCmp_t c, metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    
//     int lensols = (int) log10(metadatas_getNbSolutions(meta)) +1;
//     int lens = (int) log10(connCmp_nSols(c)) +1;
//     char temp[100], temp2[100];
//     sprintf(temp, "%d", connCmp_nSols(c));
//     for (int i = lens; i<=lensols; i++) sprintf(temp2, " ");
//     fprintf(f, "--cluster with %s sols: ", temp);
    fprintf(f, "#--cluster with %5d sols: ", connCmp_nSols(c));
    
    connCmp_componentBox( containingBox, c, metadatas_initBref(meta));
    compBox_get_containing_dsk( containingDisk, containingBox);
    
    slong d = fmpz_clog_ui(realRat_denref(compDsk_radiusref(containingDisk)),10) - fmpz_clog_ui(realRat_numref(compDsk_radiusref(containingDisk)),10); 
    slong p = fmpz_clog_ui(realRat_denref(compDsk_radiusref(containingDisk)),2) - fmpz_clog_ui(realRat_numref(compDsk_radiusref(containingDisk)),2)+50; 
    
    realApp_t cRe, cIm, rad;
    realApp_init(cRe);
    realApp_init(cIm);
    realApp_init(rad);
    
    realApp_set_realRat(cRe, compRat_realref(compDsk_centerref(containingDisk)), p);
    realApp_set_realRat(cIm, compRat_imagref(compDsk_centerref(containingDisk)), p);
    realApp_set_realRat(rad, compDsk_radiusref(containingDisk), p);
    
    
//     printf("d: %d, prec: %d\n", (int) d, (int) p);
    fprintf(f, "center: ");
    realApp_fprintn(f, cRe, d, ARB_STR_NO_RADIUS);
    fprintf(f, " + ");
    realApp_fprintn(f, cIm, d, ARB_STR_NO_RADIUS);
    fprintf(f, "j, radius: ");
    realApp_fprintn(f, rad, 5, ARB_STR_NO_RADIUS);
    
    compBox_init(containingBox);
    compDsk_init(containingDisk);
    realApp_clear(cRe);
    realApp_clear(cIm);
    realApp_clear(rad);
}

void connCmp_list_print_for_results(FILE * f, const connCmp_list_t l, metadatas_t meta){
    connCmp_list_iterator it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_print_for_results(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
}

void connCmp_print_for_results_withOutput(FILE * f, const connCmp_t c, int output, metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    
//     int lensols = (int) log10(metadatas_getNbSolutions(meta)) +1;
//     int lens = (int) log10(connCmp_nSols(c)) +1;
//     char temp[100], temp2[100];
//     sprintf(temp, "%d", connCmp_nSols(c));
//     for (int i = lens; i<=lensols; i++) sprintf(temp2, " ");
//     fprintf(f, "--cluster with %s sols: ", temp);
    if (connCmp_nSols(c) <= (10^6)-1)
        fprintf(f, "#--cluster with %5d sols: ", connCmp_nSols(c));
    else
        fprintf(f, "#--cluster with %d sols: ", connCmp_nSols(c));
    
    connCmp_componentBox( containingBox, c, metadatas_initBref(meta));
    compBox_get_containing_dsk( containingDisk, containingBox);
    
    if (output == -1) { /* rational output */
        fprintf(f, "center: ");
        realRat_print(compRat_realref(compDsk_centerref(containingDisk)));
        fprintf(f, " + ");
        realRat_print(compRat_imagref(compDsk_centerref(containingDisk)));
        fprintf(f, "j\n#%26s radius: ", " ");
        realRat_print(compDsk_radiusref(containingDisk));
    } else if (output>0) {
    
        int coeff = 4; /* = ceil(log(10)/log(2)) */
        slong prec = coeff*output; 
        
        realApp_t cRe, cIm, rad;
        realApp_init(cRe);
        realApp_init(cIm);
        realApp_init(rad);
        
        realApp_set_realRat(cRe, compRat_realref(compDsk_centerref(containingDisk)), prec);
        realApp_set_realRat(cIm, compRat_imagref(compDsk_centerref(containingDisk)), prec);
        realApp_set_realRat(rad, compDsk_radiusref(containingDisk), prec);
        
        
    //     printf("d: %d, prec: %d\n", (int) d, (int) p);
        fprintf(f, "center: ");
        realApp_fprintn(f, cRe, output, ARB_STR_MORE);
        fprintf(f, " + ");
        realApp_fprintn(f, cIm, output, ARB_STR_MORE);
        fprintf(f, "j\n#%26s radius: ", " ");
        realApp_fprintn(f, rad, 5, ARB_STR_MORE);
        
        realApp_clear(cRe);
        realApp_clear(cIm);
        realApp_clear(rad);
    
    }
    
    compBox_clear(containingBox);
    compDsk_clear(containingDisk);
}

void connCmp_list_print_for_results_withOutput(FILE * f, const connCmp_list_t l, int output, metadatas_t meta){
    connCmp_list_iterator it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_print_for_results_withOutput(f, connCmp_list_elmt(it), output, meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
}

/* old version with complicated multithreading */
// void ccluster_main_loop( connCmp_list_t qResults,  connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
//     
//     int separationFlag;
//     int widthFlag;
//     int compactFlag;
//     slong prec, depth;
//     tstar_res resTstar;
//     newton_res resNewton;
//     
//     compBox_t componentBox;
//     compDsk_t ccDisk, fourCCDisk;
//     realRat_t three, four, threeWidth;
//     compRat_t initPoint;
//     connCmp_list_t ltemp;
//     compBox_init(componentBox);
//     compDsk_init(ccDisk);
//     compDsk_init(fourCCDisk);
//     realRat_init(three);
//     realRat_init(four);
//     realRat_init(threeWidth);
//     compRat_init(initPoint);
//     connCmp_list_init(ltemp);
//     
//     connCmp_ptr ccur;
//     
//     clock_t start=clock();
//     
//     /* Real Coeff */
//     connCmp_ptr ccurConjClo, ccurConj;
//     ccurConjClo = NULL;
//     ccurConj = NULL;
//     int pushConjugFlag = 0;
//     
//     realRat_set_si(four, 4, 1);
//     realRat_set_si(three, 3, 1);
//     
// // #ifdef CCLUSTER_HAVE_PTHREAD
// //     connCmp_list_t toBeBisected;
// //     connCmp_list_init(toBeBisected);
// //     connCmp_list_t dummy;
// //     connCmp_list_init(dummy);
// // #endif
//     
//     while (!connCmp_list_is_empty(qMainLoop)) {
//         
//         //         if (metadatas_getVerbo(meta)>0) {
// //             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
// //         }
// 
//         resNewton.nflag = 0;
//         
//         ccur = connCmp_list_pop(qMainLoop);
//         
//         /* Real Coeff */
//         pushConjugFlag = 0;
//         if (metadatas_useRealCoeffs(meta)){
//             /* test if the component contains the real line in its interior */
// //             printf("number of boxes before conjugate closure: %d\n", connCmp_nb_boxes(ccur));
//             if (!connCmp_is_imaginary_positive(ccur)) {
// //                 printf("number of boxes before conjugate closure: %d\n", connCmp_nb_boxes(ccur));
//                 ccurConjClo = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
//                 connCmp_init( ccurConjClo );
//                 connCmp_set_conjugate_closure(ccurConjClo, ccur, metadatas_initBref(meta));
//                 
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//                 ccur = ccurConjClo;
// //                 printf("number of boxes after  conjugate closure: %d\n", connCmp_nb_boxes(ccur));
//             }
//         }
//         
//         connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
//         compBox_get_containing_dsk(ccDisk, componentBox);
//         compDsk_inflate_realRat(fourCCDisk, ccDisk, four);
//         realRat_mul(threeWidth, three, connCmp_widthref(ccur));
//         prec = connCmp_appPr(ccur);
//         depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
//         
//         separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
//         
//         /* Real Coeff */
//         if ( (separationFlag)&&(metadatas_useRealCoeffs(meta)) ) {
//             if (connCmp_is_imaginary_positive(ccur)) {
//                 /* check if ccur is separated from its complex conjugate */
//                 realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
//                 separationFlag = separationFlag&&(!compBox_intersection_is_not_empty_compDsk ( componentBox, fourCCDisk));
//                 realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
//             }
//         }
//         
// /*#ifdef CCLUSTER_HAVE_PTHREAD
//         if ((metadatas_useNBThreads(meta) >1)&&(!connCmp_list_is_empty(toBeBisected)))
//             separationFlag = separationFlag&&ccluster_compDsk_is_separated(fourCCDisk, toBeBisected, dummy);
// #endif */       
//         widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
//         compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
//         
//         if (metadatas_getVerbo(meta)>3) {
//             printf("---depth: %d\n", (int) depth);
//             printf("------component Box:"); compBox_print(componentBox); printf("\n");
//             printf("------separation Flag: %d\n", separationFlag);
//             printf("------widthFlag: %d\n", widthFlag); 
//             printf("------compactFlag: %d\n", compactFlag); 
//         }
//         
//         if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
// //         if ((separationFlag)) {
// //             printf("depth: %d, connCmp_nSolsref(ccur): %d, prec: %d\n", (int) depth, (int) connCmp_nSolsref(ccur), (int) prec);
//             if (connCmp_nSolsref(ccur)==-1){
//                 resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, prec, depth, meta);
//                 connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                 if (metadatas_getVerbo(meta)>3)
//                     printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
// //                 ???
//                 prec = resTstar.appPrec;
//             }
// //             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
// //             ???
// //             prec = resTstar.appPrec;
//         }
//         
//         if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && !widthFlag )
//             &&!( metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) ) ) {
//         
//             if (metadatas_haveToCount(meta)){
//                 start = clock();
//             }
//         
//             if (connCmp_nSols(ccur)==1) 
//                 compRat_set(initPoint, compBox_centerref(componentBox));
//             else
//                 connCmp_find_point_outside_connCmp( initPoint, ccur, metadatas_initBref(meta) );
//         
//             connCmp_ptr nCC;
//             nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//             connCmp_init(nCC);
//             resNewton = newton_newton_connCmp( nCC, ccur, cache, initPoint, prec, meta);
// 
// //             printf("+++depth: %d, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
//             if (resNewton.nflag) {
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//                 ccur = nCC;
//                 connCmp_increase_nwSpd(ccur);
//                 connCmp_newSuref(ccur) = 1;
//                 connCmp_appPrref(nCC) = resNewton.appPrec;
//     
//             }
//             else {
//                 connCmp_newSuref(ccur) = 0;
//                 connCmp_clear(nCC);
//                 ccluster_free(nCC);
//             }
//             if (metadatas_haveToCount(meta)){
//                 metadatas_add_Newton   ( meta, depth, resNewton.nflag, (double) (clock() - start) );
//             }
//         }
//         
//         /* Real Coeff */
//         if (metadatas_useRealCoeffs(meta)
//             && ( (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag)
//                ||( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) ) ) {
//             
//             if (connCmp_is_imaginary_positive(ccur)) {
//                 pushConjugFlag = 1;
//                 /*compute the complex conjugate*/
//                 ccurConj = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
//                 connCmp_init( ccurConj );
//                 connCmp_set_conjugate(ccurConj, ccur);
//                 
//                 /* test if initial box is symetric relatively to real axe */
//                 if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
//                     /* test if the cc intersects initial box */
//                     if ( connCmp_intersection_is_not_empty(ccurConj, metadatas_initBref(meta)) ) {
//                         /* test if the cc is confined */
// //                         if (connCmp_is_confined(ccurConj, metadatas_initBref(meta))) {
// //                             pushConjugFlag = 1;
// //                         }
// //                         else {
//                         if (!connCmp_is_confined(ccurConj, metadatas_initBref(meta))) {
//                             pushConjugFlag = 0;
//                             separationFlag = 0; 
//                           /* delete ccurConj*/
//                             connCmp_clear(ccurConj);
//                             ccluster_free(ccurConj);
//                         }
//                     }
//                     else {
//                         pushConjugFlag = 0;
//                         /* delete ccurConj*/
//                         connCmp_clear(ccurConj);
//                         ccluster_free(ccurConj);
//                     }
//                 } 
//             }
//             else {
//                 /* test if initial box is symetric relatively to real axe */
//                 if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
//                     /* test if the cc is confined and intersects initial box */
//                     if (! ( connCmp_is_confined(ccur, metadatas_initBref(meta)) 
//                          && connCmp_intersection_is_not_empty(ccur, metadatas_initBref(meta)) ) ){
// //                         /* bisect ccur until this hold */
//                         separationFlag = 0;
//                     }
//                 }
//             }
//         }
//         
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
//         else if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) {
//             metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
//             connCmp_list_push(qResults, ccur);
// //             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
// //             printf("metadatas_useRealCoeffs(meta): %d, pushConjugFlag: %d\n", metadatas_useRealCoeffs(meta), pushConjugFlag);
//             /* Real Coeff */
//             if ((metadatas_useRealCoeffs(meta))&&(pushConjugFlag)){
//                 /*compute the complex conjugate*/
//                 metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
//                 connCmp_list_push(qResults, ccurConj);
//             }
//         }
//         else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
//             connCmp_list_insert_sorted(qMainLoop, ccur);
//         }
//         
//         else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
//             connCmp_decrease_nwSpd(ccur);
// //             if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0)
// //                 connCmp_decrease_nwSpd(ccur);
//             connCmp_list_insert_sorted(qMainLoop, ccur);
//         }
//         else {
// //             if (connCmp_nSols(ccur)==0) 
// //                 printf("ici\n");
// #ifdef CCLUSTER_HAVE_PTHREAD
// //             if (metadatas_useNBThreads(meta) >1)
// //                 connCmp_list_push(toBeBisected, ccur);
// //             else {
// //                 ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
// //                 while (!connCmp_list_is_empty(ltemp))
// //                     connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
// //                 connCmp_clear(ccur);
// //                 ccluster_free(ccur);
// //             }
//             ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ccur);
//             ccluster_free(ccur);
// #else
//             ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
//             while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//             connCmp_clear(ccur);
//             ccluster_free(ccur);
// #endif
//         }
// // #ifdef CCLUSTER_HAVE_PTHREAD
// //         if ( (!connCmp_list_is_empty(toBeBisected))&&
// //                  ( (connCmp_list_is_empty(qMainLoop))||
// //                    ( realRat_cmp(connCmp_widthref(connCmp_list_first(qMainLoop)), connCmp_widthref(connCmp_list_first(toBeBisected))) !=0 ) ) ) {
// //             
// // //             printf("parallel bisect connCmp: size qMainLoop: %d, size toBeBisected: %d\n", connCmp_list_get_size(qMainLoop), connCmp_list_get_size(toBeBisected));
// //         
// //             if (connCmp_list_get_size(toBeBisected)>1) {
// //                 ccluster_parallel_bisect_connCmp_list(qMainLoop, discardedCcs, toBeBisected, cache, meta);
// //             }
// //             else { /* toBeBisected has only one element */
// //                 ccur = connCmp_list_pop(toBeBisected);
// //                 ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
// //                 while (!connCmp_list_is_empty(ltemp))
// //                         connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
// //                 connCmp_clear(ccur);
// //                 ccluster_free(ccur);
// //             }
// //         }
// // #endif
//     }
//     
//     compBox_clear(componentBox);
//     compDsk_clear(ccDisk);
//     compDsk_clear(fourCCDisk);
//     realRat_clear(three);
//     realRat_clear(four);
//     realRat_clear(threeWidth);
//     compRat_clear(initPoint);
//     connCmp_list_clear(ltemp);
// // #ifdef CCLUSTER_HAVE_PTHREAD
// //     connCmp_list_clear(toBeBisected);
// // #endif
// }


