/* ************************************************************************** */
/*  Copyright (C) 2019 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "ISSAC20/ccluster_issac20.h"

slong ccluster_issac20_discard_compBox_list( compBox_list_t boxes, cacheApp_t cache, 
//                                      int nbSols, 
                                     slong prec, metadatas_t meta){
    
    tstar_res res;
    res.appPrec = prec;
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
    compBox_ptr btemp;
    
//     /* for powerSums */
    powerSums_res resp;
    resp.appPrec = CCLUSTER_DEFAULT_PREC;
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
        /* Real Coeffs */
        if (( metadatas_useRealCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
            compBox_clear(btemp);
            ccluster_free(btemp);
            continue;
        }

        resp = powerSums_discardingTest( compDsk_centerref(bdisk), compDsk_radiusref(bdisk),
                                                        cache,
                                                        metadatas_getNbEvalPoints(meta),
                                                        metadatas_getNbPowerSums(meta),
                                                        resp.appPrec, meta, depth );
        
//         res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1,
// #ifdef CCLUSTER_STATS_PS
//                                                                                0,
// #endif
//                                                                                res.appPrec, depth, meta);
//         if (((res.nbOfSol >0) || (res.nbOfSol ==-1))&& (resp.nbOfSol==0)){
//             
//             printf("ICI!!!\n");
//             printf("------ test for disk centered in "); compRat_print(compDsk_centerref(bdisk)); printf("\n");
//             printf("------ with radius "); realRat_print( compDsk_radiusref(bdisk) ); printf("\n");
//             printf("--- tstar counting test: nbSols: %d, prec: %d \n", (int) res.nbOfSol, (int) res.appPrec );
//             
//             int verbsave = metadatas_getVerbo(meta);
//             metadatas_setVerbo(meta, 4);
//             
//             powerSums_discardingTest( compDsk_centerref(bdisk), compDsk_radiusref(bdisk),
//                                                         cache,
//                                                         metadatas_getNbEvalPoints(meta),
//                                                         metadatas_getNbPowerSums(meta),
//                                                         resp.appPrec, meta, depth );
//             
//             metadatas_setVerbo(meta, verbsave);
//             
//         }
        
        metadatas_add_PsCountingTest (meta, depth
#ifdef CCLUSTER_STATS_PS
                                                , resp.nbOfSol, -1
#endif
                                     );
        
        
        res.nbOfSol = resp.nbOfSol;
        res.appPrec = resp.appPrec;
            
        if (res.nbOfSol==0) {
            if (metadatas_haveToCount(meta)){
                metadatas_add_discarded( meta, depth);
            }
            compBox_clear(btemp);
            ccluster_free(btemp);
        }
        
        else {
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

void ccluster_issac20_bisect_connCmp( connCmp_list_t dest, connCmp_t cc, connCmp_list_t discardedCcs, cacheApp_t cache, metadatas_t meta, slong nbThreads){
    
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

    prec = ccluster_issac20_discard_compBox_list( subBoxes, cache, prec, meta);
    
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

int ccluster_issac20_main_loop( connCmp_list_t qResults,  
                              connCmp_list_t qMainLoop, 
                              connCmp_list_t discardedCcs, 
                              const realRat_t eps, 
                              cacheApp_t cache, 
                              metadatas_t meta){
    
    int res=0;
    
    int separationFlag;
    int widthFlag;
    int compactFlag;
    slong prec, depth;
//     tstar_res resTstar;
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
    
    while ((res==0)&&(!connCmp_list_is_empty(qMainLoop))) {
        
        //         if (metadatas_getVerbo(meta)>0) {
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

#ifdef WITHTSTAR
                    tstar_res resTstar;
                    resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
                    connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                     if (metadatas_getVerbo(meta)>3)
//                         printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
//                 ???
                    prec = resTstar.appPrec;
#else
                    powerSums_res resp;
                    realRat_t temp;
                    realRat_init(temp);
                    realRat_set_si(temp, 2, 1);
                    realRat_mul(temp, compDsk_radiusref(ccDisk), temp);
//                     
                    resp = powerSums_countingTest( compDsk_centerref(ccDisk), temp,
                                                        cache,
                                                        metadatas_getNbEvalPoints(meta), 
                                                        1,
                                                        prec, meta, depth );
//                     
                    connCmp_nSolsref(ccur) = resp.nbOfSol;
                    prec = resp.appPrec;
                    realRat_clear(temp);
                    
                    if (resp.nbOfSol == -1) {
                        if (metadatas_getVerbo(meta)>=2) {
                            printf("FAILURE OF A COUNTING TEST\n");
                        }
                        res = 1;
                    }   
#endif                
            }
//             printf("validate: prec avant: %d prec apres: %d\n", (int) prec, (int) resTstar.appPrec);
//             ???
//             prec = resTstar.appPrec;
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
            && ( //this is DEPRECATED: pass eps = 1/0 instead
//             (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag)
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
        if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag && (connCmp_nSols(ccur)<cacheApp_getDegree(cache)) ) {
            
            /*experimental version: check result with a Pellet test*/
            compBox_get_containing_dsk(ccDisk, componentBox);
            tstar_res resTstar;
            resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
            if (resTstar.nbOfSol != connCmp_nSols(ccur)) {
                if (metadatas_getVerbo(meta)>=2) {
                    printf("FAILURE: INCORRECT NUMBER OF SOLS IN A CLUSTER\n");
                }
                        
                res = 1;
            }
            
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
//             if (connCmp_nSols(ccur)==0) 
//                 printf("ici\n");
#ifdef CCLUSTER_HAVE_PTHREAD
            ccluster_issac20_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
#else
            ccluster_issac20_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
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
    
    return res;
}

int ccluster_issac20_algo_global( connCmp_list_t qResults, 
                                const compBox_t initialBox, 
                                const realRat_t eps, 
                                cacheApp_t cache, 
                                metadatas_t meta){
    
    clock_t start = clock();
    
    int res = 0;
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
    res = ccluster_issac20_main_loop( qResults,  qMainLoop, discardedCcs, eps, cache, meta);
    
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
    return res;
}

int metadatas_issac20_fprint(FILE * file, int res, metadatas_t meta, const realRat_t eps){
    int r=1;
    int nbTaylorShifts  = metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsInTSTests(meta);
    int nbTaylorShiftsR = metadatas_getNbTaylorsRepetedInT0Tests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta);
    int nbGraeffe       = metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeInTSTests(meta);
    int nbGraeffeR      = metadatas_getNbGraeffeRepetedInT0Tests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta);
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, " -------------------Ccluster Expe: -----------------------------------\n");
    r = fprintf(file, " -------------------Input:    ----------------------------------------\n");
    char temp[1000];
    compBox_sprint_for_stat( temp, metadatas_initBref(meta) );
    r = fprintf(file, "|box:%-65s\n", temp);
    if (realRat_is_den_zero( eps ))
        r = fprintf(file, "|eps: %-64s|\n", "+inf");
    else {
        realRat_sprint_for_stat( temp, eps );
        r = fprintf(file, "|eps: %-64s|\n", temp);
    }
    int len = 0;
    //TODO find a better way for this...
    if ( metadatas_useNewton(meta) &&
         metadatas_useTstarOptim(meta) &&
         metadatas_usePredictPrec(meta) &&
         metadatas_useAnticipate(meta) &&
         metadatas_useRealCoeffs(meta) ) len += sprintf( temp + len, " default");
    else {    
        if (metadatas_useNewton(meta)) len += sprintf( temp + len, " newton");
        if (metadatas_useTstarOptim(meta)) len += sprintf( temp + len, " tstarOpt");
        if (metadatas_usePredictPrec(meta)) len += sprintf( temp + len, " predPrec");
        if (metadatas_useAnticipate(meta)) len += sprintf( temp + len, " anticip");
        if (metadatas_useRealCoeffs(meta)) len += sprintf( temp + len, " realCoeffs");
    }
    if (metadatas_usePowerSums(meta)) len += sprintf( temp + len, " + powerSums");
    if (metadatas_forTests(meta)) len += sprintf( temp + len, " + test");
#ifdef CCLUSTER_HAVE_PTHREAD
    if (metadatas_useNBThreads(meta)>1) len += sprintf( temp + len, " %d threads", metadatas_useNBThreads(meta));
#endif
    if (metadatas_stratref(meta)->_additionalFlags !=0) 
        len += sprintf(temp +len, " %d", metadatas_stratref(meta)->_additionalFlags);
    r = fprintf(file, "|strat:%-63s|\n", temp);
    
    if (metadatas_getVerbo(meta)>=2) {
//         metadatas_count(meta);
    r = fprintf(file, " -------------------TSTest used to discard boxes----------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number DT:",                    metadatas_getNbT0Tests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingT0Tests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests DT:",       metadatas_get_time_T0Tests(meta),    " " );
    r = fprintf(file, " -------------------TSTest used to validate clusters------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number VT:",                    metadatas_getNbTSTests(meta),        " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number in Newton iterations:",        metadatas_getNbTSTestsInNewton(meta), " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of tests without conclusion:", metadatas_getNbFailingTSTests(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in tests VT:",       metadatas_get_time_TSTests(meta),    " " );
    r = fprintf(file, " -------------------Taylor shifts-------------------------------------\n");
    r = fprintf(file, "|%-39s %14d |%13d|\n", "total number TS:",                    nbTaylorShifts + nbTaylorShiftsR, nbTaylorShiftsR );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in discarding TSTests TS:",    metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsRepetedInT0Tests(meta), metadatas_getNbTaylorsRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in validating TSTests TS:",    metadatas_getNbTaylorsInTSTests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta), metadatas_getNbTaylorsRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number in Newton iterations:",        metadatas_getNbTaylorsInNewton(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Taylor shifts:",  metadatas_get_time_Taylors(meta),    " " );
    r = fprintf(file, " -------------------Graeffe Iterations--------------------------------\n");
    r = fprintf(file, "|%-39s %14d |%13d|\n", "total number GR:",                       nbGraeffe + nbGraeffeR, nbGraeffeR );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in discarding TSTests GR:",       metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeRepetedInT0Tests(meta), metadatas_getNbGraeffeRepetedInT0Tests(meta) );
    r = fprintf(file, "|%-39s %14d |%13d|\n", "number in validating TSTests GR:",       metadatas_getNbGraeffeInTSTests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta), metadatas_getNbGraeffeRepetedInTSTests(meta) );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number in Newton iterations:",        metadatas_getNbGraeffeInNewton(meta), " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in Graeffe Iterations:", metadatas_get_time_Graeffe(meta),    " " );
    if (metadatas_useNewton(meta)){
    r = fprintf(file, " -------------------Newton Iterations---------------------------------\n");
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number NE:",                       metadatas_getNbNewton(meta),         " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of fails:",                    metadatas_getNbFailingNewton(meta),  " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time spent in newton:",         metadatas_get_time_Newtons(meta),    " " );
    }
    r = fprintf(file, " -------------------Other---------------------------------------------\n");
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in getApproximation:",           metadatas_get_time_Approxi(meta),    " " );
    if (metadatas_useAnticipate(meta)){
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Anticipate:",                 metadatas_get_time_Anticip(meta),    " " );
    }
    if (metadatas_usePowerSums(meta)){
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of Ps counting tests:",  metadatas_getNbPsCountingTest(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests:",          metadatas_get_time_PSTests(meta),    " " );
#ifdef CCLUSTER_STATS_PS_MACIS
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests V:",        metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests D:",        metadatas_get_time_PSTests(meta)-metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Evaluation:",                 metadatas_get_time_Evaluat(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of evaluations:",        metadatas_getNbEval(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -2:",                 metadatas_getNbM2(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -1:",                 metadatas_getNbM1(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of errors:",             metadatas_getNbEr(meta),    " " );
#endif 
#ifdef CCLUSTER_STATS_PS
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests V:",        metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests D:",        metadatas_get_time_PSTests(meta)-metadatas_get_time_PSTestV(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in Evaluation:",                 metadatas_get_time_Evaluat(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of evaluations:",        metadatas_getNbEval(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative:",      metadatas_getNbTN(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive:",     metadatas_getNbFP(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative 1:",      metadatas_getNbTN1(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive 1:",     metadatas_getNbFP1(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative 2:",      metadatas_getNbTN2(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive 2:",     metadatas_getNbFP2(meta),    " " );
#endif 
    }
    r = fprintf(file, " -------------------Precision-----------------------------------------\n");
    r = metadatas_boxes_by_prec_fprint ( file, meta );
    
#ifdef CCLUSTER_EXPERIMENTAL    
    if (CCLUSTER_EXP_NUM_T0(meta)||CCLUSTER_EXP_NUM_T1(meta)||CCLUSTER_INC_TEST(meta)) {
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in getDerivative:",              metadatas_get_time_Derivat(meta),    " " );
    r = fprintf(file, "|%-39s %14f %14s|\n", "time in evaluate:",                   metadatas_get_time_Evaluat(meta),    " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of evaluations:",              metadatas_getNbEval(meta),    " " );
    }
#endif 
    }
   
    r = fprintf(file, " -------------------Output:   ----------------------------------------\n");
    if (res==0)
    r = fprintf(file, "|%-39s %14s %14s|\n", "failure:",                            "0",      " " );
    else if (res==1)
    r = fprintf(file, "|%-39s %14s %14s|\n", "failure:",                            "1 1",      " " );
    else if (res==2)
    r = fprintf(file, "|%-39s %14s %14s|\n", "failure:",                            "1 2",      " " ); 
    else if (res==3)
    r = fprintf(file, "|%-39s %14s %14s|\n", "failure:",                            "1 3",      " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of clusters:",                 metadatas_getNbValidated(meta),      " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "number of solutions:",                metadatas_getNbSolutions(meta),      " " );
    r = fprintf(file, " -------------------Stats:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree depth:",                         metadatas_getDepth(meta),            " " );
    r = fprintf(file, "|%-39s %14d %14s|\n", "tree size:",                          metadatas_getNbExplored(meta),       " " );
    }
    r = fprintf(file, "|%-39s %14f %14s|\n", "total time:",                         metadatas_get_time_CclusAl(meta),    " " );
    r = fprintf(file, " ---------------------------------------------------------------------\n");
    }
    return r;
}
