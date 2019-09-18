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
/* test */
#include "ccluster/cclusterDAC.h"

int  ccluster_compDsk_is_separated_DAC( const compDsk_t d, 
                                        connCmp_list_t qMainLoop, 
                                        connCmp_list_t qResults,
                                        connCmp_list_t qAllResults, 
                                        connCmp_list_t discardedCcs ){
    int res = 1;
    connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    it = connCmp_list_begin(qResults);
    while ( res && (it!=connCmp_list_end()) ) {
        res = res && (! connCmp_intersection_is_not_empty_compDsk( connCmp_list_elmt(it) , d));
        it = connCmp_list_next(it);
    }
    it = connCmp_list_begin(qAllResults);
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

void ccluster_prep_loop_DAC( connCmp_list_t qMainLoop, connCmp_list_t qPrepLoop, connCmp_list_t discardedCcs, cacheApp_t cache, metadatas_t meta) {
    
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
            connCmp_list_insert_sorted_inv(qMainLoop, ctemp);
        else {
            ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, cache, meta, 1);
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

void ccluster_main_loop_DAC( connCmp_list_t qResults, 
                             connCmp_list_t qAllResults,
                             connCmp_list_t qMainLoop, 
                             connCmp_list_t discardedCcs, 
                             int nbSols,
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
    
    clock_t start=clock();
    
    /* Real Coeff */
    connCmp_ptr ccurConjClo, ccurConj;
    ccurConjClo = NULL;
    ccurConj = NULL;
    int pushConjugFlag = 0;
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    
// #ifdef CCLUSTER_HAVE_PTHREAD
//     connCmp_list_t toBeBisected;
//     connCmp_list_init(toBeBisected);
//     connCmp_list_t dummy;
//     connCmp_list_init(dummy);
// #endif
    
    int nbSolsAlreadyFound = 0;
    while ((!connCmp_list_is_empty(qMainLoop))&&(nbSolsAlreadyFound<nbSols)) {
        
        //         if (metadatas_getVerbo(meta)>0) {
//             printf("cclusterDAC.c, ccluster_main_loop_DAC, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }
            
        resNewton.nflag = 0;
        
        ccur = connCmp_list_pop(qMainLoop); //pop the CC with smallest boxes
        
         /* Real Coeff */
        pushConjugFlag = 0;
        if (metadatas_realCoeffs(meta)){
            /* test if the component contains the real line in its interior */
            if (!connCmp_is_imaginary_positive(ccur)) {
                ccurConjClo = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
                connCmp_init( ccurConjClo );
                connCmp_set_conjugate_closure(ccurConjClo, ccur,metadatas_initBref(meta));
                
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = ccurConjClo;
            }
        }
        
        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
        compBox_get_containing_dsk(ccDisk, componentBox);
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four); //??? why four?
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        if (connCmp_isSep(ccur)==0) {
            separationFlag = ccluster_compDsk_is_separated_DAC(fourCCDisk, qMainLoop, qResults, qAllResults, discardedCcs);
            /* Real Coeff */
            if ( (separationFlag)&&(metadatas_realCoeffs(meta)) ) {
                if (connCmp_is_imaginary_positive(ccur)) {
                    /* check if ccur is separated from its complex conjugate */
                    realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
                    separationFlag = separationFlag&&(!compBox_intersection_is_not_empty_compDsk ( componentBox, fourCCDisk));
                    realRat_neg( compRat_imagref(compDsk_centerref(fourCCDisk)), compRat_imagref(compDsk_centerref(fourCCDisk)) );
                }
            }
            if (separationFlag)
                connCmp_isSep(ccur)=1;
        }
        else
            separationFlag = 1;
        
/*#ifdef CCLUSTER_HAVE_PTHREAD
        if ((metadatas_useNBThreads(meta) >1)&&(!connCmp_list_is_empty(toBeBisected)))
            separationFlag = separationFlag&&ccluster_compDsk_is_separated(fourCCDisk, toBeBisected, dummy);
#endif*/        
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
            
            if (connCmp_nSolsref(ccur)==-1){
                resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0,0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
                prec = resTstar.appPrec;
            }
            
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && !widthFlag )
            &&!( metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) ) ) {
            
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
        if (metadatas_realCoeffs(meta)){
            if (connCmp_is_imaginary_positive(ccur)) {
                /*compute the complex conjugate*/
                pushConjugFlag = 1;
                ccurConj = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
                connCmp_init( ccurConj );
                connCmp_set_conjugate(ccurConj, ccur);
                
                /* test if initial box is symetric relatively to real axe */
                if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
                    /* test if the cc intersects initial box */
                    if ( connCmp_intersection_is_not_empty(ccurConj, metadatas_initBref(meta)) ) {
                        /* test if the cc is confined */
                        if (connCmp_is_confined(ccurConj, metadatas_initBref(meta))) {
                            pushConjugFlag = 1;
                        }
                        else {
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
        
        if (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag){
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
            nbSolsAlreadyFound += connCmp_nSols(ccur);

            if ((metadatas_realCoeffs(meta))&&(pushConjugFlag)){
                /*compute the complex conjugate*/
                metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
                connCmp_list_push(qResults, ccurConj);
                nbSolsAlreadyFound += connCmp_nSols(ccurConj);
            }
            
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
            nbSolsAlreadyFound += connCmp_nSols(ccur);
            
            /* Real Coeff */
            if ((metadatas_realCoeffs(meta))&&(pushConjugFlag)){
                /*compute the complex conjugate*/
                metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
                connCmp_list_push(qResults, ccurConj);
                nbSolsAlreadyFound += connCmp_nSols(ccurConj);
            }
            
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
            connCmp_list_insert_sorted_inv(qMainLoop, ccur);
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            connCmp_decrease_nwSpd(ccur);
            connCmp_list_insert_sorted_inv(qMainLoop, ccur);
        }
        else {
//             if (connCmp_nSols(ccur)==0) 
//                 printf("ici\n");
// #ifdef CCLUSTER_HAVE_PTHREAD
//             if (metadatas_useNBThreads(meta) >1)
//                 connCmp_list_push(toBeBisected, ccur);
//             else {
//                 ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
//                 while (!connCmp_list_is_empty(ltemp))
// //                     connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//                     connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//             }
// #else
            ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
// #endif
        }
// #ifdef CCLUSTER_HAVE_PTHREAD
//         if ( (!connCmp_list_is_empty(toBeBisected))&&
//                  ( (connCmp_list_is_empty(qMainLoop))||
//                    ( realRat_cmp(connCmp_widthref(connCmp_list_first(qMainLoop)), connCmp_widthref(connCmp_list_first(toBeBisected))) !=0 ) ) ) {
//             
//             if (connCmp_list_get_size(toBeBisected)>1) {
//                 ccluster_parallel_bisect_connCmp_list(qMainLoop, discardedCcs, toBeBisected, cache, meta);
//             }
//             else { /* toBeBisected has only one element */
//                 ccur = connCmp_list_pop(toBeBisected);
//                 ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
//                 while (!connCmp_list_is_empty(ltemp))
// //                         connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
//                         connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
//                 connCmp_clear(ccur);
//                 ccluster_free(ccur);
//             }
//         }
// #endif
        
    }
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(threeWidth);
    compRat_clear(initPoint);
    connCmp_list_clear(ltemp);
// #ifdef CCLUSTER_HAVE_PTHREAD
//     connCmp_list_clear(toBeBisected);
// #endif
}

void ccluster_DAC_first( connCmp_list_t qResults, 
                         connCmp_list_t qAllResults,
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         int nbSols,
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
    
    connCmp_list_t qPrepLoop;
//     connCmp_list_init(qMainLoop);
    connCmp_list_init(qPrepLoop);
//     connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qPrepLoop, initialCC);
//     printf("preploop: \n");
    ccluster_prep_loop_DAC( qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     printf("mainloop: \n");
    ccluster_main_loop_DAC( qResults, qAllResults, qMainLoop, discardedCcs, nbSols, eps, cache, meta);
    
    
    realRat_clear(factor);
//     connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
//     connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
}

void ccluster_DAC_next( connCmp_list_t qResults, 
                        connCmp_list_t qAllResults,
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         int nbSols,
//                          const compBox_t initialBox, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta){
    
//     chronos_tic_CclusAl(metadatas_chronref(meta)); 
    clock_t start = clock();
    ccluster_main_loop_DAC( qResults, qAllResults, qMainLoop, discardedCcs, nbSols, eps, cache, meta);
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
}

void ccluster_DAC_first_interface_forJulia( connCmp_list_t qResults,
                                            connCmp_list_t qAllResults,
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         void(*func)(compApp_poly_t, slong), 
                         int nbSols,
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         int st, 
                         int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
    ccluster_DAC_first( qResults, qAllResults, qMainLoop, discardedCcs, nbSols, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
}

void ccluster_DAC_next_interface_forJulia( connCmp_list_t qResults,
                                           connCmp_list_t qAllResults,
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         void(*func)(compApp_poly_t, slong), 
                         int nbSols,
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         int st, 
                         int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), st>>7);
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
    ccluster_DAC_next( qResults, qAllResults, qMainLoop, discardedCcs,nbSols, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qResults, meta);
//         connCmp_list_print_for_results(stdout, qResults, 500, 40, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
}

// void ccluster_DAC_first_interface_func( void(*func)(compApp_poly_t, slong), const compBox_t initialBox, const realRat_t eps, int st, int verb){
// 
//     cacheApp_t cache;
//     strategies_t strat;
//     metadatas_t meta;
//     connCmp_list_t qRes;
//     connCmp_list_t qMainLoop, discardedCcs;
//     
//     cacheApp_init(cache, func);
//     strategies_init(strat);
// //     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     
//     metadatas_init(meta, initialBox, strat, verb);
//     connCmp_list_init(qRes);
//     connCmp_list_init(qMainLoop);
//     connCmp_list_init(discardedCcs);
//     
//     ccluster_DAC_first( qRes, qMainLoop, discardedCcs, initialBox, eps, cache, meta);
//     metadatas_count(meta);
//     metadatas_fprint(stdout, meta, eps);
//     
//     if (verb>=3) {
//         connCmp_list_print_for_results(stdout, qRes, meta);
// //         connCmp_list_print_for_results(stdout, qRes, 500, 40, meta);
//     }
//     
//     cacheApp_clear(cache);
//     strategies_clear(strat);
//     metadatas_clear(meta);
//     connCmp_list_clear(qRes);
//     
//     connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(discardedCcs);
// }
