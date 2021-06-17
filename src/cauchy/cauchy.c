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

#include "cauchy/cauchy.h"

slong cauchy_discard_compBox_list( compBox_list_t boxes, 
                                     compBox_list_t bDiscarded,
                                     cacheApp_t cache, 
                                     cacheCauchy_t cacheCau,
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
    
    while (!compBox_list_is_empty(boxes)){
        
        btemp = compBox_list_pop(boxes);
        compBox_get_containing_dsk(bdisk, btemp);
        depth = compDsk_getDepth(bdisk, metadatas_initBref( meta));
        metadatas_add_explored( meta, depth);
        
        /* Real Coeffs */
        if (( metadatas_useRealCoeffs(meta) ) && ( compBox_is_imaginary_negative_strict(btemp) ) ) {
//             printf("ici\n");
//             if (metadatas_getDrSub(meta)==0){
                compBox_clear(btemp);
                ccluster_free(btemp);
//             } else {
//                 compBox_list_push(bDiscarded, btemp);
//             }
            continue;
        }
        
        /* test */
//         if ( compBox_skipref(btemp)==1 ) {
// //             compBox_skipref(btemp) = 0;
//             compBox_list_push(ltemp, btemp);
//             continue;
//         }
        
        /* deterministic exclusion test */
        cauchyTest_res resCauchy;
        resCauchy.appPrec = CCLUSTER_DEFAULT_PREC;
        resCauchy = cauchyTest_deterministic_exclusion_test( compDsk_centerref(bdisk), compDsk_radiusref(bdisk), 
                                                            NULL, 0,0,
                                                            cache, cacheCau, resCauchy.appPrec, CAUCHYTEST_INEXCLUSI, meta, depth);
        
//         resCauchy = cauchyTest_probabilistic_exclusion_test( compDsk_centerref(bdisk), compDsk_radiusref(bdisk), 
//                                                             NULL, 0,0,
//                                                             cache, cacheCau, resCauchy.appPrec, meta, depth);
        
        res.appPrec = resCauchy.appPrec;
        res.nbOfSol = resCauchy.nbOfSol;
        /* fin test */
//         res = tstar_interface( cache, bdisk, compBox_get_nbMSol(btemp), 1, 0, res.appPrec, depth, meta);  
//         
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("------ result of tstar test: %d\n", res.nbOfSol);
//         }
        
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
                if (res.nbOfSol>0) {
                    btemp->nbMSol = res.nbOfSol;
                    /* mark adjacent boxes in boxes for skip */
//                     compBox_list_iterator it = compBox_list_begin(boxes);
//                     compRat_t center;
//                     compRat_init(center);
//                     compRat_set(center, compBox_centerref(btemp));
//                     realRat_add(compRat_realref(center), compRat_realref(center), compBox_bwidthref(btemp));
//                     realRat_add(compRat_imagref(center), compRat_imagref(center), compBox_bwidthref(btemp));
//                     while(it!= compBox_list_end()){
//                         compBox_ptr b = compBox_list_elmt(it);
//                         
//                         if (realRat_cmp( compRat_realref( compBox_centerref( b ) ), compRat_realref( center ) ) > 0)
//                             it = compBox_list_end();
//                         else {
//                             if (   (realRat_cmp( compRat_realref( compBox_centerref( b ) ), compRat_realref( center ) ) <= 0)
//                                 && (realRat_cmp( compRat_imagref( compBox_centerref( b ) ), compRat_imagref( center ) ) <= 0) )
//                                 if (compBox_skipref(b) == 0 )
//                                     compBox_skipref(b) = 1;
//                             it = compBox_list_next(it);
//                         }
//                     }
//                     compRat_clear(center);
                }
                compBox_list_push(ltemp, btemp);
        }
    }
    
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
    return res.appPrec;
}

void cauchy_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              compBox_list_t bDiscarded, 
                              cacheApp_t cache, 
                              cacheCauchy_t cacheCau,
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
    
    prec = cauchy_discard_compBox_list( subBoxes, bDiscarded, cache, cacheCau, prec, meta);
    
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

void cauchy_prep_loop( compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t qPrepLoop, 
                         connCmp_list_t discardedCcs, 
                         cacheApp_t cache,
                         cacheCauchy_t cacheCau,
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
            cauchy_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cache, cacheCau, meta, metadatas_useNBThreads(meta));
            
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

/* let res = D(c',r'), Delta=D(c,r) */
/* set r' st |c-c'|+r/theta <= r' <= 5/4 r */
void cauchy_setRadSup( compDsk_t res, const compDsk_t Delta, const realRat_t theta ){
    
    compRat_t temp;
    compRat_init(temp);
    compRat_sub(temp, compDsk_centerref(Delta), compDsk_centerref(res) );
    
    compApp_t tempApp;
    compApp_init(tempApp);
    
    slong prec = CCLUSTER_DEFAULT_PREC;
    compApp_set_compRat(tempApp, temp, prec);
    realApp_t distApp;
    realApp_init(distApp);
    compApp_abs( distApp, tempApp, prec );
    
    realRat_t dist, rOtheta, ubound;
    realRat_init(dist);
    realRat_init(rOtheta);
    realRat_init(ubound);
    
    realApp_get_realRat(dist, distApp);
    realRat_div( rOtheta, compDsk_radiusref(Delta), theta );
    realRat_set_si(ubound, 5,4);
    realRat_mul(ubound, ubound, compDsk_radiusref(Delta) );
    
    realRat_add(compDsk_radiusref(res), dist, rOtheta );
    
    while( realRat_cmp( compDsk_radiusref(res), ubound ) >0 ) {
        prec = 2*prec;
        compApp_set_compRat(tempApp, temp, prec);
        compApp_abs( distApp, tempApp, prec );
        realApp_get_realRat(dist, distApp);
        realRat_add(compDsk_radiusref(res), dist, rOtheta );
    }
    
    compRat_clear(temp);
    compApp_clear(tempApp);
    realApp_clear(distApp);
    realRat_clear(dist);
    realRat_init(rOtheta);
    realRat_init(ubound);
}

/* assume Delta = D(c,r) contains m and has isolation ratio theta >=2 */
/* computes a disk res = D(c',r') such that*/
/* Delta and res contain the same roots */
/* either r' <= eps */
/*     or res is m/((2m-2)*theta) rigid */
slong cauchy_compressionIntoRigidDisk( compDsk_t res, const compDsk_t Delta, slong m, const realRat_t theta, const realRat_t eps,
                                       cacheApp_t cache,
                                       cacheCauchy_t cacheCau,
                                       slong prec, metadatas_t meta, slong depth) {
    
    clock_t start=clock();
    
    realRat_t epsp;
    realRat_init(epsp);
    realRat_mul_si(epsp, eps, 2);
    /* if 2*eps \geq r, return the disk D(c,r/2) */
    if (metadatas_getVerbo(meta)>=3) {
            printf("# --- cauchy.c: cauchy_compressionIntoRigidDisk; eps:");
            realRat_print(eps);
            printf("\n# --- cauchy.c: cauchy_compressionIntoRigidDisk; 2*eps:");
            realRat_print(epsp);
            printf("\n# --- cauchy.c: cauchy_compressionIntoRigidDisk; r:");
            realRat_print(compDsk_radiusref(Delta));
            printf("\n");
    }
        
    if ( (realRat_cmp( epsp, compDsk_radiusref(Delta) )>=0) ){
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("# --- cauchy.c: cauchy_compressionIntoRigidDisk;  2*eps >= radius of input disc; \n");
        }
        
        compDsk_set(res, Delta );
        realRat_div_ui( compDsk_radiusref(res), compDsk_radiusref(res), 2);
        
        realRat_clear(epsp);
        return prec;
    }
    /* here eps < r/2 */ 
    /* set epsp = m*eps/theta */
    realRat_div(epsp, eps, theta);
    realRat_mul_si(epsp, epsp, m);
    
    /* compute a disk of radius less than m*eps/theta containing s1(p,Delta) */
    slong appPrec = cauchyTest_computeS1compDsk( res, theta, Delta, m, cache, cacheCau, epsp, meta, depth );
    
    /* if m=1 then res is equivalent to Delta and has radius less than eps/2 */
    /* return the disk res */
    if (m==1) {
        
        if (metadatas_getVerbo(meta)>=3) {
            printf("# --- cauchy.c: cauchy_compressionIntoRigidDisk;  m=1; use s1 \n");
        }
        
        realRat_clear(epsp);
        return appPrec;
    }
    /* here m>1 */
    
    realRat_t relativeError;
    realRat_init(relativeError);
    realRat_set_si(relativeError, 19, 10);
    /* set c' = center(s1)/m */
    compRat_div_ui( compDsk_centerref(res), compDsk_centerref(res), (ulong) m );
    /* set u st |c-c'|+r/theta <= u <= 5/4 r */
    cauchy_setRadSup( res, Delta, theta );
    /* compute r such that 1 <= r / r_{d+1-m} (s1/m, p) <= relativeError */ 
    realRat_t radInf;
    realRat_init(radInf);
    realRat_zero(radInf);
    
    if (metadatas_getVerbo(meta)>=3) {
            printf("# --- cauchy.c: cauchy_compressionIntoRigidDisk; call root radii \n");
    }
    cauchyRootRadii_root_radius( compDsk_centerref(res),
                                 radInf,        /* radInf = 0 < r_{d+1-m}(center, p) */
                                 compDsk_radiusref(res),        /* radSup > r_{d+1-m}(center, p) */
                                 relativeError, /* want relativeError*radInf >= radSup */ 
                                 epsp,
                                 theta,   /*isolation ratio of the disk in which is computed rr */ 
                                 m,
                                 cacheCau, cache, meta );
    
    realRat_clear(radInf);
    realRat_clear(epsp);
    realRat_clear(relativeError);
    
    metadatas_addComp_nb( meta, 1);
    metadatas_add_time_CompTot(meta, (double) (clock() - start));
    
    return appPrec;
}

void cauchy_main_loop( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         cacheCauchy_t cacheCau,
                         metadatas_t meta){
    
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
    
    while (!connCmp_list_is_empty(qMainLoop)) {
        
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
            if (connCmp_nSolsref(ccur)==-1){
                    
//                     cauchyTest_res resCauchy = cauchyTest_deterministic_counting_test( compDsk_centerref(ccDisk), compDsk_radiusref(ccDisk),
//                                                                                   cache, cacheCau, prec, meta, depth);
                    realRat_mul_si(compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 2);
                    
                    cauchyTest_res resCauchy = cauchyTest_probabilistic_counting( ccDisk, cache, cacheCau, prec, meta, depth);
                    realRat_div_ui(compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 2);
                    connCmp_nSolsref(ccur) = resCauchy.nbOfSol;
                    prec = resCauchy.appPrec;
                    if (metadatas_getVerbo(meta)>=3)
                        printf("#------nb sols after Cauchy: %d\n", (int) connCmp_nSolsref(ccur));
                     
//                     resTstar = tstar_interface( cache, ccDisk, cacheApp_getDegree(cache), 0, 0, prec, depth, meta);
//                     connCmp_nSolsref(ccur) = resTstar.nbOfSol;
//                     prec = resTstar.appPrec;
//                     if (metadatas_getVerbo(meta)>=3)
//                         printf("------nb sols after tstar: %d\n", (int) connCmp_nSolsref(ccur));
                    
                    realRat_t epsComp;
                    realRat_init(epsComp);
                    realRat_set_si( epsComp, (slong) connCmp_nSolsref(ccur), (slong) connCmp_nSolsref(ccur) + 3*cacheApp_getDegree(cache) );
                    realRat_mul( epsComp, epsComp, epsComp );
                    realRat_mul_si(compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 2);
                    
                    if ( realRat_cmp (epsComp, compDsk_radiusref(ccDisk) ) < 0  ) {
                        
                        if (metadatas_getVerbo(meta)>=3) {
                            printf("\n# Test compression into rigid disc for a CC with %d roots \n", connCmp_nSolsref(ccur));
                            printf("# Delta: "); compDsk_print( ccDisk ); printf("\n");
                            printf("# Required radius: "); realRat_print(epsComp); printf("\n");
                        }
                        
                        compDsk_t res;
                        compDsk_init(res);
                        
                        
                        /* ccDisk is at least 2 isolated */
                        realRat_t theta;
                        realRat_init(theta);
                        realRat_set_si(theta, 2, 1);
                        
                        slong precres = cauchy_compressionIntoRigidDisk( res, ccDisk, connCmp_nSolsref(ccur), theta, epsComp,
                                                                         cache, cacheCau, prec, meta, depth);
                        
                        if (metadatas_getVerbo(meta)>=3) {
                            printf("# Precision after compression: %ld\n", precres);
                            printf("# res: "); compDsk_print( res ); printf("\n\n");
                        }
                        realRat_clear(theta);
                        compDsk_clear(res);
                    }
                    realRat_div_ui(compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), 2);
                    realRat_init(epsComp);
                    
//                  if ( (connCmp_nSolsref(ccur) == 1) && (!widthFlag) ) {
//                      
//                      if (metadatas_getVerbo(meta)>=3) {
//                          printf("# Refine CC with 1st power sum: \n");
//                      }
//                          
//                      realRat_t halfeps;
//                      realRat_init(halfeps);
//                      realRat_div_ui(halfeps, eps, 4);
//                      
//                      realRat_t isoRatio;
//                      realRat_init(isoRatio);
//                      realRat_set_si(isoRatio, 2, 1);
//                      
//                      compDsk_t res;
//                      compDsk_init(res);
//                      
//                      
//                      slong appPrec = cauchyTest_computeS1compDsk( res, isoRatio, ccDisk,
//                                                                   connCmp_nSolsref(ccur), cache, cacheCau,
//                                                                   halfeps, meta, depth );
//                      if (metadatas_getVerbo(meta)>=3)
//                         printf("#------appPrec: %ld\n", appPrec);
//                      
//                      connCmp_ptr nCC;
//                      nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
//                      connCmp_init(nCC);
//                      
//                      compBox_list_ptr ltemp;
//                      ltemp = connCmp_boxesref(ccur);
//                      compBox_list_t ltemp2;
//                      compBox_list_init(ltemp2);
//                      compBox_ptr btemp;
//                      
//                      while (compBox_list_get_size(ltemp)>0){
//                         btemp = compBox_list_pop(ltemp);
//                         subdBox_quadrisect_with_compDsk( ltemp2, btemp, res, halfeps);
//                         compBox_clear(btemp);
//                         ccluster_free(btemp);
//                      }
//                      compBox_list_swap(ltemp, ltemp2);
//                      compBox_list_clear(ltemp2);
//                      btemp = compBox_list_pop(ltemp);
//                      realRat_set(connCmp_widthref(nCC), compBox_bwidthref(btemp));
//                      connCmp_insert_compBox(nCC, btemp);
//                      while (!compBox_list_is_empty(ltemp))
//                      connCmp_insert_compBox(nCC, compBox_list_pop(ltemp));
//                      connCmp_nSols(nCC) = connCmp_nSols(ccur);
//                      connCmp_isSep(nCC) = connCmp_isSep(ccur);
//                      fmpz_set(connCmp_nwSpdref(nCC), connCmp_nwSpdref(ccur));
//                      
//                      connCmp_clear(ccur);
//                      ccluster_free(ccur);
//                      ccur = nCC;
//                      connCmp_appPrref(ccur) = appPrec;
//                      
//                      connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
//                      realRat_mul(threeWidth, three, connCmp_widthref(ccur));
//                      widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), eps)<=0);
//                      compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
//                      
//                      if (metadatas_getVerbo(meta)>=3) {
//                         printf("#------widthFlag: %d\n", widthFlag);
//                         printf("#------compactFlag: %d\n", compactFlag);
//                      }
//                      
//                      realRat_clear(halfeps);
//                      compDsk_clear(res);
//                      realRat_clear(isoRatio);
//                      compBox_list_clear(ltemp2);
//                  }
//                  if ( (metadatas_getVerbo(meta)>=3) && (connCmp_nSolsref(ccur) == 1) ) {
//                      compDsk_t res;
//                      compDsk_init(res);
//                      realRat_t isoRatio;
//                      realRat_init(isoRatio);
//                      realRat_set_si(isoRatio, 2, 1);
//                      printf("#Test cauchyTest_computeS1compDsk: \n");
//                      slong appPrec = cauchyTest_computeS1compDsk( res, isoRatio, ccDisk,
//                                                                   connCmp_nSolsref(ccur), cache, cacheCau,
//                                                                   eps, meta, depth );
//                      
//                      tstar_res resTstar = tstar_interface( cache, res, cacheApp_getDegree(cache), 0, 0, CCLUSTER_DEFAULT_PREC, depth, meta);
//                      
//                      compApp_t app;
//                      compApp_init(app);
//                      compApp_set_compDsk(app, compDsk_centerref(res), compDsk_radiusref(res), CCLUSTER_DEFAULT_PREC);
//                      printf("#------appPrec: %ld\n", appPrec);
// //                      printf("#------res    : "); compDsk_print( res ); printf("\n");
//                      printf("#------res    : "); compApp_printd( app, 10 ); printf("\n");
//                      printf("#------result of tstar : prec: %ld, nbofSols: %d \n", resTstar.appPrec, resTstar.nbOfSol);
//                      compDsk_clear(res);
//                      realRat_clear(isoRatio);
//                      compApp_clear(app);
//                  }
                
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
//             resNewton = newton_newton_connCmp( nCC, ccur, cache, initPoint, prec, meta);
            resNewton = newton_cauchy_newton_connCmp( nCC, ccur, cache, cacheCau, initPoint, prec, meta);

//             printf("+++depth: %d, connCmp_nSolsref(ccur): %d, res_newton: %d \n", depth, connCmp_nSols(ccur), resNewton.nflag);
            if (resNewton.nflag) {
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = nCC;
                connCmp_increase_nwSpd(ccur);
                connCmp_newSuref(ccur) = 1;
                connCmp_appPrref(ccur) = resNewton.appPrec;
    
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
        if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag && (connCmp_nSols(ccur)<cacheApp_getDegree(cache)) ) {
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
            cauchy_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cache, cacheCau, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
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

void cauchy_algo_global( connCmp_list_t qResults, 
                           compBox_list_t bDiscarded,
                           const compBox_t initialBox, 
                           const realRat_t eps, 
                           cacheApp_t cache, 
                           cacheCauchy_t cacheCau,
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
    cauchy_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cache, cacheCau, meta);
    
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
}

int metadatas_cauchy_fprint(FILE * file, metadatas_t meta, const realRat_t eps){
    int r=1;
//     int nbTaylorShifts  = metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsInTSTests(meta);
//     int nbTaylorShiftsR = metadatas_getNbTaylorsRepetedInT0Tests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta);
//     int nbGraeffe       = metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeInTSTests(meta);
//     int nbGraeffeR      = metadatas_getNbGraeffeRepetedInT0Tests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta);
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, "# -------------------Ccluster Expe: -----------------------------------\n");
    r = fprintf(file, "# -------------------Input:    ----------------------------------------\n");
    char temp[1000];
    compBox_sprint_for_stat( temp, metadatas_initBref(meta) );
    r = fprintf(file, "#|box:%-65s\n", temp);
    if (realRat_is_den_zero( eps ))
        r = fprintf(file, "#|eps: %-64s|\n", "+inf");
    else {
        realRat_sprint_for_stat( temp, eps );
        r = fprintf(file, "#|eps: %-64s|\n", temp);
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
    r = fprintf(file, "#|strat:%-63s|\n", temp);
    
    if (metadatas_getVerbo(meta)>=2) {
//         metadatas_count(meta);
    r = fprintf(file, "# -------------------Cauchy exclusion tests----------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number ET:",                       metadatas_getNbCauchyExTests(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals for proba. ET:",         metadatas_getNbCauchyExEvalsP(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals for certi. ET:",         metadatas_getNbCauchyExEvalsD(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in tests  ET:",          metadatas_get_time_CauExTo(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals for proba. ET:",         metadatas_get_time_CauExEP(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals for certi. ET:",         metadatas_get_time_CauExED(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing divs     ET:",         metadatas_get_time_CauExDS(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing s0s      ET:",         metadatas_get_time_CauExCS(meta),    " " );
    r = fprintf(file, "# -------------------Cauchy counting tests----------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number CT:",                       metadatas_getNbCauchyCoTests(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals for proba. CT:",         metadatas_getNbCauchyCoEvalsP(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of evals for certi. CT:",         metadatas_getNbCauchyCoEvalsD(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in tests  CT:",          metadatas_get_time_CauCoTo(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals for proba. CT:",         metadatas_get_time_CauCoEP(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time   in evals for certi. CT:",         metadatas_get_time_CauCoED(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing divs     CT:",         metadatas_get_time_CauCoDS(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in computing s0s      ET:",         metadatas_get_time_CauCoCS(meta),    " " );
    r = fprintf(file, "# -------------------Compression into rigid discs---------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number of compression:",           metadatas_getComp_nb(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in compression:",       metadatas_get_time_CompTot(meta),    " " );
    if (metadatas_useNewton(meta)){
    r = fprintf(file, "# -------------------Newton Iterations---------------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number NE:",                       metadatas_getNbNewton(meta),         " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of fails:",                    metadatas_getNbFailingNewton(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in newton:",         metadatas_get_time_Newtons(meta),    " " );
    }
    r = fprintf(file, "# -------------------Other---------------------------------------------\n");
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in getApproximation:",           metadatas_get_time_Approxi(meta),    " " );
//     if (metadatas_useAnticipate(meta)){
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time in Anticipate:",                 metadatas_get_time_Anticip(meta),    " " );
//     }
//     if (metadatas_usePowerSums(meta)){
// //     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of Ps counting tests:",  metadatas_getNbPsCountingTest(meta),    " " );
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests:",          metadatas_get_time_PSTests(meta),    " " );
// #ifdef CCLUSTER_STATS_PS_MACIS
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests V:",        metadatas_get_time_PSTestV(meta),    " " );
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests D:",        metadatas_get_time_PSTests(meta)-metadatas_get_time_PSTestV(meta),    " " );
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time in Cauchy exclusion tests:",     metadatas_get_time_PSTests(meta),    " " );
    double timeInEval = metadatas_get_time_CauExEP(meta) + metadatas_get_time_CauExED(meta)
                      + metadatas_get_time_CauCoEP(meta) + metadatas_get_time_CauCoED(meta);
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time in Evaluation:",                 timeInEval,    " " );
    int    nbOfEvals  = metadatas_getNbCauchyExEvalsP(meta) + metadatas_getNbCauchyExEvalsD(meta)
                      + metadatas_getNbCauchyCoEvalsP(meta) + metadatas_getNbCauchyCoEvalsD(meta);
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number of evaluations:",        nbOfEvals,    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -2:",                 metadatas_getNbM2(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of -1:",                 metadatas_getNbM1(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of errors:",             metadatas_getNbEr(meta),    " " );
// #endif 
// #ifdef CCLUSTER_STATS_PS
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests V:",        metadatas_get_time_PSTestV(meta),    " " );
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Ps counting tests D:",        metadatas_get_time_PSTests(meta)-metadatas_get_time_PSTestV(meta),    " " );
//     r = fprintf(file, "|%-39s %14f %14s|\n", "time in Evaluation:",                 metadatas_get_time_Evaluat(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of evaluations:",        metadatas_getNbEval(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative:",      metadatas_getNbTN(meta),    " " );
//     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive:",     metadatas_getNbFP(meta),    " " );
// //     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative 1:",      metadatas_getNbTN1(meta),    " " );
// //     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive 1:",     metadatas_getNbFP1(meta),    " " );
// //     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of True Negative 2:",      metadatas_getNbTN2(meta),    " " );
// //     r = fprintf(file, "|%-39s %14d %14s|\n", "total number of False Positive 2:",     metadatas_getNbFP2(meta),    " " );
// #endif 
//     }
    r = fprintf(file, "# -------------------Precision-----------------------------------------\n");
    r = metadatas_boxes_by_prec_fprint ( file, meta );
     
    }
   
    r = fprintf(file, "# -------------------Output:   ----------------------------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of clusters:",                 metadatas_getNbValidated(meta),      " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "number of solutions:",                metadatas_getNbSolutions(meta),      " " );
    r = fprintf(file, "# -------------------Stats:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    r = fprintf(file, "#|%-39s %14d %14s|\n", "tree depth:",                         metadatas_getDepth(meta),            " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "tree size:",                          metadatas_getNbExplored(meta),       " " );
    }
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time:",                         metadatas_get_time_CclusAl(meta),    " " );
    r = fprintf(file, "# ---------------------------------------------------------------------\n");
    }
    return r;
}
