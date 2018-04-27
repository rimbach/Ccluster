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

#include "ccluster/cclusterDAC.h"

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
            ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
//             ccluster_bisect_connCmp( ltemp, ctemp, discardedCcs, cache, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_push(qPrepLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ctemp);
            free(ctemp);
        }       
    }
    
    connCmp_list_clear(ltemp);
    realRat_clear(halfwidth);
    realRat_clear(diam);
}

void ccluster_main_loop_DAC( connCmp_list_t qResults,  connCmp_list_t qMainLoop, connCmp_list_t discardedCcs, const realRat_t eps, cacheApp_t cache, metadatas_t meta){
    
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
    
    realRat_set_si(four, 4, 1);
    realRat_set_si(three, 3, 1);
    
#ifdef CCLUSTER_HAVE_PTHREAD
    connCmp_list_t toBeBisected;
    connCmp_list_init(toBeBisected);
    connCmp_list_t dummy;
    connCmp_list_init(dummy);
#endif
    
//     while (!connCmp_list_is_empty(qMainLoop)) {
    while (connCmp_list_is_empty(qResults)) {    
        
        resNewton.nflag = 0;
        
        ccur = connCmp_list_pop(qMainLoop); //pop the CC with smallest boxes
        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
        compBox_get_containing_dsk(ccDisk, componentBox);
        compDsk_inflate_realRat(fourCCDisk, ccDisk, four); //??? why four?
        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
        prec = connCmp_appPr(ccur);
        depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
        
        separationFlag = ccluster_compDsk_is_separated(fourCCDisk, qMainLoop, discardedCcs);
#ifdef CCLUSTER_HAVE_PTHREAD
        if (!connCmp_list_is_empty(toBeBisected))
            separationFlag = separationFlag&&ccluster_compDsk_is_separated(fourCCDisk, toBeBisected, dummy);
#endif        
        widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
        compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
//         if ((separationFlag)) {
//             printf("depth: %d, connCmp_nSolsref(ccur): %d\n", depth, connCmp_nSolsref(ccur));
            if (connCmp_nSolsref(ccur)==-1){
                resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, prec, depth, meta);
                connCmp_nSolsref(ccur) = resTstar.nbOfSol;
            }
//             printf("validate: prec avant: %d prec apres: %d\n", prec, resTstar.appPrec);
            prec = resTstar.appPrec;
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && !widthFlag )
            &&!( metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) ) ) {
            
            if (connCmp_nSols(ccur)==1) 
                compRat_set(initPoint, compBox_centerref(componentBox));
            else
                connCmp_find_point_outside_connCmp( initPoint, ccur, metadatas_initBref(meta) );
            
            connCmp_ptr nCC;
            nCC = (connCmp_ptr) malloc (sizeof(connCmp));
            connCmp_init(nCC);
            resNewton = newton_newton_connCmp( nCC, ccur, cache, initPoint, prec, meta);
            metadatas_add_Newton   ( meta, depth, resNewton.nflag );
//             printf("+++depth: %d, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
            if (resNewton.nflag) {
                connCmp_clear(ccur);
                free(ccur);
                ccur = nCC;
                connCmp_increase_nwSpd(ccur);
                connCmp_newSuref(ccur) = 1;
//                 printf("newton: prec avant: %d prec apres: %d\n", prec, resNewton.appPrec);
                connCmp_appPrref(nCC) = resNewton.appPrec;
                /* adjust precision: try to make it lower */
//                 if (prec == resNewton.appPrec) {
// //                     printf("ICI\n");
//                     connCmp_appPrref(ccur) = CCLUSTER_MAX(prec/2,CCLUSTER_DEFAULT_PREC);
//                 }
//                 else 
//                     connCmp_appPrref(ccur) = resNewton.appPrec;
    
            }
            else {
                connCmp_newSuref(ccur) = 0;
                connCmp_clear(nCC);
                free(nCC);
            }
                
        }
        
        if (metadatas_useStopWhenCompact(meta) && compactFlag && (connCmp_nSols(ccur)==1) && separationFlag){
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
//             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            connCmp_list_push(qResults, ccur);
//             printf("+++depth: %d, validated with %d roots\n", (int) depth, connCmp_nSols(ccur));
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
//             connCmp_list_insert_sorted(qMainLoop, ccur);
            connCmp_list_insert_sorted_inv(qMainLoop, ccur);
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            connCmp_decrease_nwSpd(ccur);
//             if (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0)
//                 connCmp_decrease_nwSpd(ccur);
//             connCmp_list_insert_sorted(qMainLoop, ccur);
            connCmp_list_insert_sorted_inv(qMainLoop, ccur);
        }
        else {
//             if (connCmp_nSols(ccur)==0) 
//                 printf("ici\n");
#ifdef CCLUSTER_HAVE_PTHREAD
            if (metadatas_useNBThreads(meta) >1)
                connCmp_list_push(toBeBisected, ccur);
            else {
                ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
                while (!connCmp_list_is_empty(ltemp))
//                     connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
                    connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
                connCmp_clear(ccur);
                free(ccur);
            }
#else
            ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta,1);
            while (!connCmp_list_is_empty(ltemp))
//                 connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
                connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            free(ccur);
#endif
        }
#ifdef CCLUSTER_HAVE_PTHREAD
        if ( (!connCmp_list_is_empty(toBeBisected))&&
                 ( (connCmp_list_is_empty(qMainLoop))||
                   ( realRat_cmp(connCmp_widthref(connCmp_list_first(qMainLoop)), connCmp_widthref(connCmp_list_first(toBeBisected))) !=0 ) ) ) {
            
            if (connCmp_list_get_size(toBeBisected)>1) {
                ccluster_parallel_bisect_connCmp_list(qMainLoop, discardedCcs, toBeBisected, cache, meta);
            }
            else { /* toBeBisected has only one element */
                ccur = connCmp_list_pop(toBeBisected);
                ccluster_bisect_connCmp( ltemp, ccur, discardedCcs, cache, meta, metadatas_useNBThreads(meta));
                while (!connCmp_list_is_empty(ltemp))
//                         connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
                        connCmp_list_insert_sorted_inv(qMainLoop, connCmp_list_pop(ltemp));
                connCmp_clear(ccur);
                free(ccur);
            }
        }
#endif
        
    }
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(threeWidth);
    compRat_clear(initPoint);
    connCmp_list_clear(ltemp);
#ifdef CCLUSTER_HAVE_PTHREAD
    connCmp_list_clear(toBeBisected);
#endif
}

void ccluster_DAC_first( connCmp_list_t qResults, 
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta){
    
    chronos_tic_CclusAl(metadatas_chronref(meta));
    
    realRat_t factor;
    realRat_init(factor);
    realRat_set_si(factor, 5, 4);
    
    compBox_ptr bEnlarged;
    bEnlarged = (compBox_ptr) malloc (sizeof(compBox));
    compBox_init(bEnlarged);
    compBox_inflate_realRat(bEnlarged, initialBox, factor);
    compBox_nbMSolref(bEnlarged) = cacheApp_getDegree ( cache );
    
    connCmp_ptr initialCC;
    initialCC = (connCmp_ptr) malloc (sizeof(connCmp));
    connCmp_init_compBox(initialCC, bEnlarged);
    
    connCmp_list_t qPrepLoop;
//     connCmp_list_init(qMainLoop);
    connCmp_list_init(qPrepLoop);
//     connCmp_list_init(discardedCcs);
    
    connCmp_list_push(qPrepLoop, initialCC);
//     printf("preploop: \n");
    ccluster_prep_loop_DAC( qMainLoop, qPrepLoop, discardedCcs, cache, meta);
//     printf("mainloop: \n");
    ccluster_main_loop_DAC( qResults,  qMainLoop, discardedCcs, eps, cache, meta);
    
    
    realRat_clear(factor);
//     connCmp_list_clear(qMainLoop);
    connCmp_list_clear(qPrepLoop);
//     connCmp_list_clear(discardedCcs);
    
    chronos_toc_CclusAl(metadatas_chronref(meta));
    
}

void ccluster_DAC_next( connCmp_list_t qResults, 
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
//                          const compBox_t initialBox, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         metadatas_t meta){
    
    chronos_tic_CclusAl(metadatas_chronref(meta)); 
    ccluster_main_loop_DAC( qResults,  qMainLoop, discardedCcs, eps, cache, meta);
    chronos_toc_CclusAl(metadatas_chronref(meta));
    
}

void ccluster_DAC_first_interface_forJulia( connCmp_list_t qResults, 
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         void(*func)(compApp_poly_t, slong), 
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         int st, 
                         int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
    ccluster_DAC_first( qResults, qMainLoop, discardedCcs, initialBox, eps, cache, meta);
    
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
                         connCmp_list_t qMainLoop,
                         connCmp_list_t discardedCcs,
                         void(*func)(compApp_poly_t, slong), 
                         const compBox_t initialBox, 
                         const realRat_t eps, 
                         int st, 
                         int verb){
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
    metadatas_init(meta, initialBox, strat, verb);
    
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
    ccluster_DAC_next( qResults, qMainLoop, discardedCcs, eps, cache, meta);
    
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