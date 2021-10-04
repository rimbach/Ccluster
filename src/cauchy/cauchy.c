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
                                     slong nbMSols, 
                                     slong prec, metadatas_t meta){
    
    int level = 4;
    
    tstar_res res;
    res.appPrec = prec;
    
    slong depth;
    
    compBox_list_t ltemp;
    compDsk_t bdisk;
    compBox_list_init(ltemp);
    compDsk_init(bdisk);
    
//     compRat_t compDist;
//     compRat_init(compDist);
    
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
        
        /* exclusion test */
        cauchyTest_res resCauchy;
        resCauchy.nbOfSol = 1;
        resCauchy.appPrec = CCLUSTER_DEFAULT_PREC;

#ifdef CERTIFIED        
        resCauchy = cauchyTest_deterministic_exclusion_test( compDsk_centerref(bdisk), compDsk_radiusref(bdisk), 
                                                            NULL, 0,0,
                                                            cache, cacheCau, resCauchy.appPrec, CAUCHYTEST_INEXCLUSI, meta, depth);
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---cauchy deterministic exclusion test: res: %d\n", resCauchy.nbOfSol );
        }
#else
        resCauchy = cauchyTest_probabilistic_exclusion_test( compDsk_centerref(bdisk), compDsk_radiusref(bdisk), 
                                                            NULL, 0,0,
                                                            cache, cacheCau, resCauchy.appPrec, meta, depth);
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---cauchy probabilistic exclusion test: res: %d\n", resCauchy.nbOfSol );
        }
#endif
        
        res.appPrec = resCauchy.appPrec;
        res.nbOfSol = resCauchy.nbOfSol;
        
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
                }
                
                /* there is at least a root in 2bdisk */
//                 realRat_mul_si( compDsk_radiusref(bdisk), compDsk_radiusref(bdisk), 2);
//                 compBox_ptr btemp2;
//                 compBox_list_t ltemp2;
//                 compBox_list_init(ltemp2);
//                 while (!compBox_list_is_empty(boxes)){
//                     btemp2 = compBox_list_pop(boxes);
//                     /* check if btemp2 is 4-adjacent to btemp */
//                     compRat_comp_distance( compDist, compBox_centerref(btemp), compBox_centerref(btemp2) );
//                     if ( ( ( realRat_cmp( compRat_realref(compDist), compBox_bwidthref(btemp) ) ==0 ) &&
//                            ( realRat_is_zero( compRat_imagref(compDist) ) ) ) ||
//                          ( ( realRat_cmp( compRat_imagref(compDist), compBox_bwidthref(btemp) ) ==0 ) &&
//                            ( realRat_is_zero( compRat_realref(compDist) ) ) )  ) {
//                         compBox_list_push(ltemp, btemp2);
//                     } else {
//                         compBox_list_push(ltemp2, btemp2);
//                     }
//                 }
//                 compBox_list_swap(boxes, ltemp2);
//                 compBox_list_clear(ltemp2);
                
                if (nbMSols==1) {
                    /* there is at least a root in 2bdisk */
                    realRat_mul_si( compDsk_radiusref(bdisk), compDsk_radiusref(bdisk), 2);
                    compBox_ptr btemp2;
                    compBox_list_t ltemp2;
                    compBox_list_init(ltemp2);
                    while (!compBox_list_is_empty(boxes)){
                        btemp2 = compBox_list_pop(boxes);
                        /* test if btemp2 intersects 2bdisk */
                        if (compBox_intersection_is_not_empty_compDsk ( btemp2, bdisk)){
                            /* keep btemp2 */
                            compBox_list_push(ltemp2, btemp2);
                        } else {
//                             printf("##############################ICI##########################\n");
                            /* discard btemp2 */
                            if (metadatas_haveToCount(meta)){
                                metadatas_add_discarded( meta, depth);
                            }
                            if (metadatas_getDrSub(meta)==0){
                                compBox_clear(btemp2);
                                ccluster_free(btemp2);
                            } else {
                                compBox_list_push(bDiscarded, btemp2);
                            }
                        }
                    }
                    compBox_list_swap(boxes, ltemp2);
                    compBox_list_clear(ltemp2);
                }
                compBox_list_push(ltemp, btemp);
        }
    }
    
    compBox_list_swap(boxes, ltemp);
    compBox_list_clear(ltemp);
    compDsk_clear(bdisk);
    
//     compRat_clear(compDist);
    
    return res.appPrec;
}

void cauchy_bisect_connCmp( connCmp_list_t dest, 
                              connCmp_t cc, 
                              connCmp_list_t discardedCcs, 
                              compBox_list_t bDiscarded, 
                              cacheApp_t cache, 
                              cacheCauchy_t cacheCau,
                              slong nbMsols, 
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
    
    prec = cauchy_discard_compBox_list( subBoxes, bDiscarded, cache, cacheCau, nbMsols, prec, meta);
    
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
                connCmp_isRigref(ctemp) = connCmp_isRigref(cc);
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
            cauchy_bisect_connCmp( ltemp, ctemp, discardedCcs, bDiscarded, cache, cacheCau, -1, meta, metadatas_useNBThreads(meta));
            
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

/* try to upper bound the number of roots in ccur */
/* at least one sol per connected comp in qMainLoop, */
/* 2 if use real coeffs and the CC does not contain the real line */
slong cauchy_MaxNumberOfRootsInCC( const connCmp_t ccur,
                                   connCmp_list_t qMainLoop,
                                   slong nbSolsInQResults,
                                   cacheApp_t cache, 
                                   metadatas_t meta ){
    slong nbMaxSol=connCmp_nSolsref(ccur);
    if (nbMaxSol==-1) {
        nbMaxSol=0;
        connCmp_list_iterator it = connCmp_list_begin(qMainLoop);
        while (it!=connCmp_list_end()){
            
            if (connCmp_nSolsref( connCmp_list_elmt( it ) ) > -1) {
                nbMaxSol += connCmp_nSolsref( connCmp_list_elmt( it ) ) ;
                if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(connCmp_list_elmt( it ))) )
                    nbMaxSol += connCmp_nSolsref( connCmp_list_elmt( it ) ) ; 
            } else {
                nbMaxSol += 1 ;
                if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(connCmp_list_elmt( it ))) )
                    nbMaxSol += 1 ; 
            }
                
            it = connCmp_list_next(it);
        }
        nbMaxSol += nbSolsInQResults;
//         printf("nbMaxSol: %ld\n", nbMaxSol);
        nbMaxSol = cacheApp_getDegree(cache) - nbMaxSol;
        if ( (metadatas_useRealCoeffs(meta)) && (connCmp_is_imaginary_positive(ccur)) )
            nbMaxSol = nbMaxSol/2;
    }
    return nbMaxSol;
}

/* let ccDisk = D(c,r) be at least 2-isolated and contain m roots and ccur*/
/* let epsp = max( eps, r/( (1 + 3d/m)^2 ) ) */
/* compute a CC with containing disk D(c',r') s.t   */
/* either D(c',r') is rigid with radius r' > epsp */
/* or     D(c',r')           has radius r' <= epsp */
/* in the latter case D(c', (1 + 3d/m)r') is (1 + 3d/m)-isolated */
connCmp_ptr cauchy_compression( connCmp_ptr ccur,
                                int *rigidFlag,
                                int *widthFlag,
                                int *isolaFlag,
                                const compDsk_t ccDisk,
                                slong m,
                                const realRat_t eps,
                                cacheApp_t cache, 
                                cacheCauchy_t cacheCau,
                                slong prec,
                                metadatas_t meta,
                                slong depth ) {
    
    int level = 3;
    
    /* let epsp = max( eps, r/( (1 + 3d/m)^2 ) ) */
    realRat_t epsp;
    realRat_init(epsp);
    realRat_set_si( epsp, (slong) connCmp_nSolsref(ccur), (slong) connCmp_nSolsref(ccur) + 3*cacheApp_getDegree(cache) );
    realRat_mul( epsp, epsp, epsp );
    realRat_mul( epsp, epsp, compDsk_radiusref(ccDisk) );
    /* case where eps = +inf; eps = 1/0 */
    if (realRat_is_den_zero(eps)==0)
        realRat_max_2_realRat(epsp, eps);
    
    /* let epspp = (2/9) epsp */
    /* compute a disk D(c'',r'') s.t. */
    /* either D(c'',r'') is rigid with radius r'' > epspp = (2/9) epsp */
    /*     or D(c'',r'')           has radius r'' <= epspp = (2/9) epsp */
    realRat_t epspp;
    realRat_init(epspp);
    realRat_mul_si(epspp, epsp, 2);
    realRat_div_ui(epspp, epspp, 9);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("\n#------Test compression into rigid disc for a CC with depth %ld with %d roots \n", depth, connCmp_nSolsref(ccur));
        printf("#---------Delta: "); compDsk_print( ccDisk ); printf("\n");
        printf("#---------Required radius of containing disk of CC: "); realRat_print(epsp); printf("\n");
        printf("#---------Required radius of disk:                  "); realRat_print(epspp); printf("\n");
        printf("#---------Epsilon                :                  "); realRat_print(eps); printf("\n");
    }
    
    compDsk_t res;
    compDsk_init(res);
    /* ccDisk is at least 2 isolated */
    realRat_t theta;
    realRat_init(theta);
    realRat_set_si(theta, 2, 1);
    
    slong precres = cauchy_compressionIntoRigidDisk( res, ccDisk, m, theta, epspp, cache, cacheCau, prec, meta, depth);
    
    if (metadatas_getVerbo(meta)>=level) {
        printf("#---------Precision after compression: %ld\n", precres);
        printf("#---------res: "); compDsk_print( res ); printf("\n");
        compApp_t c;
        compApp_init(c);
        compApp_set_compRat(c, compDsk_centerref(res), precres);
        printf("#---------app: center: "); compApp_printd( c, 10 );
        realApp_set_realRat(compApp_realref(c), compDsk_radiusref(res), precres);
        printf("# radius: "); realApp_printd( compApp_realref(c), 10 ); printf("\n");
        compApp_clear(c);
    }
    
    /* actualize connected component: */
    if ( realRat_cmp( compDsk_radiusref(res), epspp) > 0) { /* if r'' > epspp = (2/9)epsp */
        /* cover D(c'',r'') with at most 9 sub-boxes of the subdivision tree of width >r'' and <= 2r'' */
        /* resulting CC has width at most 6*r'' */
        /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*6*r'' = (9/2)r'' */
        ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
        *rigidFlag = 1;
        connCmp_isRigref(ccur) = 1;
        if (metadatas_getVerbo(meta)>=level) {
            printf("#---------rigid disk\n");
        }
    } else {
        if ( realRat_cmp( epsp, eps) == 0) { /* r'' <= (2/9)eps */
            /* cover D(c'',r'') with at most 9 sub-boxes of the subdivision tree of width >r'' and <= 2r'' */
            /* resulting CC has width at most 6*r'' <= 6*(2/9)eps */
            /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*6*(2/9)eps = eps */
            ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
            *widthFlag = 1;
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------disk with size smaller than epsilon\n");
            }
        } else {
            /* the Disk D(c'', (1 + 3d/m)r'') is (1 + 3d/m) isolated */
            realRat_mul_si( compDsk_radiusref(res), compDsk_radiusref(res), 
                            connCmp_nSolsref(ccur) + 3*cacheApp_getDegree(cache) );
            realRat_div_ui( compDsk_radiusref(res), compDsk_radiusref(res), 
                            connCmp_nSolsref(ccur) );
            /* cover D(c'',(1 + 3d/m)r'') with at most 9 sub-boxes of the subdivision tree */
            /* of size (1 + 3d/m)r''=(1/3)(1 + 3d/m)epsp */
            /* resulting CC has width at most 3*(1/3)(1 + 3d/m)*epsp = (1 + 3d/m)*epsp */
            /* the containing disk D(c',r') of the resulting CC has radius at most (3/4)*(1 + 3d/m)*epsp */
            /* and is (1 + 3d/m) isolated */
            ccur = cauchy_actualizeCCafterCompression( ccur, res, precres, meta );
            *isolaFlag = 1;
            if (metadatas_getVerbo(meta)>=level) {
                printf("#---------isolated disk\n");
            }
        }
        
    }
    
    realRat_clear(theta);
    compDsk_clear(res);
    realRat_clear(epsp);
    realRat_clear(epspp);
    
    return ccur;
}

newton_res cauchy_newton ( connCmp_ptr * ccur,
                           const compBox_t componentBox,
                           cacheApp_t cache, 
                           cacheCauchy_t cacheCau,
                           slong prec,
                           metadatas_t meta){
    
    newton_res resNewton;
    
    compRat_t initPoint;
    compRat_init (initPoint);
//     if (connCmp_nSols(*ccur)==1) 
//         compRat_set(initPoint, compBox_centerref(componentBox));
//     else
        connCmp_find_point_outside_connCmp( initPoint, *ccur, metadatas_initBref(meta) );
    
    connCmp_ptr nCC;
    nCC = (connCmp_ptr) ccluster_malloc (sizeof(connCmp));
    connCmp_init(nCC);
    resNewton = newton_cauchy_newton_connCmp( nCC, *ccur, cache, cacheCau, initPoint, prec, meta);
            
    compRat_clear(initPoint);
    
    if (resNewton.nflag) {
        connCmp_clear(*ccur);
        ccluster_free(*ccur);
        *ccur = nCC;
        connCmp_increase_nwSpd(*ccur);
        connCmp_newSuref(*ccur) = 1;
        connCmp_appPrref(*ccur) = resNewton.appPrec;
    
    }
    else {
        connCmp_newSuref(*ccur) = 0;
        connCmp_clear(nCC);
        ccluster_free(nCC);
    }
            
    return resNewton;
}

void cauchy_conjugate ( connCmp_ptr * ccurConj, int * pushConjugFlag, int * separationFlag, 
                        connCmp_ptr ccur, metadatas_t meta ){
    
    if (connCmp_is_imaginary_positive(ccur)) {
        * pushConjugFlag = 1;
        /*compute the complex conjugate*/
        * ccurConj = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
        connCmp_init( *ccurConj );
        connCmp_set_conjugate(*ccurConj, ccur);
        
        /* test if initial box is symetric relatively to real axe */
        if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
            /* test if the cc intersects initial box */
            if ( connCmp_intersection_is_not_empty(*ccurConj, metadatas_initBref(meta)) ) {
                /* test if the cc is confined */
                if (!connCmp_is_confined(*ccurConj, metadatas_initBref(meta))) {
                    *pushConjugFlag = 0;
                    *separationFlag = 0; 
                  /* delete ccurConj*/
                    connCmp_clear(*ccurConj);
                    ccluster_free(*ccurConj);
                }
            }
            else {
                *pushConjugFlag = 0;
                /* delete ccurConj*/
                connCmp_clear(*ccurConj);
                ccluster_free(*ccurConj);
            }
        } 
    } else {
        /* test if initial box is symetric relatively to real axe */
        if ( !realRat_is_zero(compRat_imagref(compBox_centerref(metadatas_initBref(meta)))) ) {
            /* test if the cc is confined and intersects initial box */
            if (! ( connCmp_is_confined(ccur, metadatas_initBref(meta)) 
                 && connCmp_intersection_is_not_empty(ccur, metadatas_initBref(meta)) ) ){
                /* bisect ccur until this hold */
                *separationFlag = 0;
            }
        }
    }
}

int cauchy_main_loop( connCmp_list_t qResults,  
                         compBox_list_t bDiscarded,
                         connCmp_list_t qMainLoop, 
                         connCmp_list_t discardedCcs, 
                         const realRat_t eps, 
                         cacheApp_t cache, 
                         cacheCauchy_t cacheCau,
                         metadatas_t meta){
    
    /* for prints */
    int level = 3;
    
    int separationFlag;
    int widthFlag;
    int compactFlag;
    int rigidFlag; /* the disk is m/( (2m-2)theta )-rigid */
    int isolaFlag; /* the disk is (1 + 3d/m)-isolated */
    
    slong prec, depth;
//     tstar_res resTstar;
    newton_res resNewton;
    slong nbSolsInQResults = 0;
    
    compBox_t componentBox;
    compDsk_t ccDisk, twoCCDisk, fourCCDisk;
    realRat_t two, three, four, mp2, threeWidth, neps;
    connCmp_list_t ltemp;
    compBox_init(componentBox);
    compDsk_init(ccDisk);
    compDsk_init(twoCCDisk);
    compDsk_init(fourCCDisk);
    realRat_init(two);
    realRat_init(three);
    realRat_init(four);
    realRat_init(mp2);
    realRat_init(neps);
    realRat_init(threeWidth);
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
    realRat_set_si(two, 2, 1);
    realRat_set(neps, eps);
    
    int failure = 0;
    
    while ( (failure==0)&&(!connCmp_list_is_empty(qMainLoop)) ) {
        
        //         if (metadatas_getVerbo(meta)>0) {
//             printf("ccluster.c, ccluster_main_loop, size of queue: %d \n", connCmp_list_get_size(qMainLoop) );
//         }
        widthFlag      = 0;
        compactFlag    = 0;
        rigidFlag      = 0;
        isolaFlag      = 0;
        resNewton.nflag = 0;
        
        ccur = connCmp_list_pop(qMainLoop);
        
        /* try to upper bound the number of roots in ccur */
        /* at least one sol per connected comp in qMainLoop, */
        /* 2 if use real coeffs and the CC does not contain the real line */
        slong nbMaxSol= cauchy_MaxNumberOfRootsInCC( ccur, qMainLoop, nbSolsInQResults, cache,  meta );
        
        /* Real Coeff */
        pushConjugFlag = 0;
        if (metadatas_useRealCoeffs(meta)){
            /* test if the component contains the real line in its interior */
            if (!connCmp_is_imaginary_positive(ccur)) {
                ccurConjClo = ( connCmp_ptr ) ccluster_malloc (sizeof(connCmp));
                connCmp_init( ccurConjClo );
                connCmp_set_conjugate_closure(ccurConjClo, ccur, metadatas_initBref(meta));
                connCmp_clear(ccur);
                ccluster_free(ccur);
                ccur = ccurConjClo;
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
        rigidFlag      = 0;
        isolaFlag      = 0;
        
        if (metadatas_getVerbo(meta)>=level) {
            printf("\n#---depth: %d\n", (int) depth);
//             printf("------component Box:"); compBox_print(componentBox); printf("\n");
//             printf("#------length of working queue :         %d\n", connCmp_list_get_size(qMainLoop));
//             printf("#------length of results queue :         %d\n", connCmp_list_get_size(qResults));
            printf("#------number of boxes in ccur:          %d\n", connCmp_nb_boxes(ccur));
            printf("#------connCmp_nSolsref(ccur):           %d\n", connCmp_nSolsref(ccur)); 
            if (connCmp_nSolsref(ccur) == -1)
                printf("#------max number of roots in ccur:      %ld\n", nbMaxSol); 
            printf("#------separation Flag: %d\n", separationFlag);
            printf("#------widthFlag: %d\n", widthFlag); 
            printf("#------compactFlag: %d\n", compactFlag);
            printf("#------rigidFlag: %d\n", rigidFlag); 
            printf("#------isolaFlag: %d\n", isolaFlag); 
        }
        
        if ((separationFlag)&&(connCmp_newSu(ccur)==0)) {
            if (connCmp_nSolsref(ccur)==-1){

                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------run Cauchy probabilistic counter:\n");
                    }
                    
                    compDsk_inflate_realRat(twoCCDisk, ccDisk, two);
                    cauchyTest_res resCauchy = cauchyTest_probabilistic_counting( twoCCDisk, cache, cacheCau, prec, meta, depth);
                    connCmp_nSolsref(ccur) = resCauchy.nbOfSol;
                    prec = resCauchy.appPrec;
                    slong m = connCmp_nSols(ccur);
                    
#ifndef CERTIFIED
                    if (resCauchy.nbOfSol == -1) {
                        failure = 1;
                        continue;
                        printf("#FAILURE: A DISC IS NOT 2-ISOLATED \n");
                    }
#endif                    
                    if (metadatas_getVerbo(meta)>=level) {
                        printf("#------nb sols: %d\n", (int) connCmp_nSolsref(ccur));
                    }
                    
                    if (metadatas_useCompression(meta)){
                        clock_t start2 = clock();
                        
                        if (metadatas_getVerbo(meta)>=level) {
                            printf("\n#---Compression into rigid disc for a CC with depth %ld with %d roots \n", depth, connCmp_nSolsref(ccur));
                            printf("#------Disk: "); compDsk_print(twoCCDisk); printf("\n");
                        }

#ifdef CAUCHY_CERTIFIED
                        realRat_div_ui(neps, eps, m+2);
#endif
                        ccur = cauchy_compression( ccur, &rigidFlag, &widthFlag, &isolaFlag, 
                                                   twoCCDisk, m, neps, cache, cacheCau, prec, meta, depth );
                        
                        connCmp_componentBox(componentBox, ccur, metadatas_initBref(meta));
                        compBox_get_containing_dsk(ccDisk, componentBox);
                        realRat_mul(threeWidth, three, connCmp_widthref(ccur));
                        
                        if (metadatas_getVerbo(meta)>=level) {
                            printf("#------time spent in compression: %f\n", ((double) (clock() - start2))/CLOCKS_PER_SEC );
                            printf("\n");
                        }
                        
                    }
                    
                
            }
        }
        
        if (connCmp_nSolsref(ccur)!=-1){
#ifdef CAUCHY_CERTIFIED
             slong m = connCmp_nSolsref(ccur);
             realRat_div_ui(neps, eps, m+2);
#endif             
             realRat_mul(threeWidth, three, connCmp_widthref(ccur));
             widthFlag      = (realRat_cmp( compBox_bwidthref(componentBox), neps)<=0);
             compactFlag    = (realRat_cmp( compBox_bwidthref(componentBox), threeWidth)<=0);
             depth = connCmp_getDepth(ccur, metadatas_initBref(meta));
#ifdef CAUCHY_CERTIFIED
             realRat_mul_si( compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), m+2);
             separationFlag = ccluster_compDsk_is_separated(ccDisk, qMainLoop, discardedCcs);
             /* Real Coeff */
             if ( (separationFlag)&&(metadatas_useRealCoeffs(meta)) ) {
                 if (connCmp_is_imaginary_positive(ccur)) {
                     /* check if ccur is separated from its complex conjugate */
                     realRat_neg( compRat_imagref(compDsk_centerref(ccDisk)), compRat_imagref(compDsk_centerref(ccDisk)) );
                     separationFlag = separationFlag&&(!compBox_intersection_is_not_empty_compDsk ( componentBox, ccDisk));
                     realRat_neg( compRat_imagref(compDsk_centerref(ccDisk)), compRat_imagref(compDsk_centerref(ccDisk)) );
                 }
             }
             realRat_div_ui( compDsk_radiusref(ccDisk), compDsk_radiusref(ccDisk), m+2);
#endif          
             if (metadatas_getVerbo(meta)>=level) {
                 printf("#------New Disk: "); compDsk_print( ccDisk ); printf("\n");
                 printf("#------separation Flag: %d\n", separationFlag);
                 printf("#------widthFlag      : %d\n", widthFlag); 
                 printf("#------compactFlag    : %d\n", compactFlag);
                 printf("#------depth          : %ld\n", depth);
                 printf("#------rigidFlag: %d\n", rigidFlag); 
                 printf("#------isolaFlag: %d\n", isolaFlag); 
                 printf("\n");
             }
        }
        
        if ( ( separationFlag && (connCmp_nSols(ccur) >0) && metadatas_useNewton(meta) && 
               ( (!widthFlag)||( connCmp_nSols(ccur)== cacheApp_getDegree(cache) ) )
               && ( (!metadatas_useCompression(meta))||(connCmp_isRigref(ccur)==0)) ) ) {
            
            if (metadatas_getVerbo(meta)>=level)
                printf("#------run Newton:\n");
        
            start = clock();
        
            resNewton = cauchy_newton ( &ccur, componentBox, cache, cacheCau, prec, meta);
            
            if (metadatas_getVerbo(meta)>=level)
                printf("#---------res_newton: %d \n", resNewton.nflag);
            
            if (metadatas_haveToCount(meta)){
                metadatas_add_Newton   ( meta, depth, resNewton.nflag, (double) (clock() - start) );
            }
        }
        
        /* Real Coeff */
        if ( metadatas_useRealCoeffs(meta)
            && ( ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag ) ) ) {
            
            cauchy_conjugate ( &ccurConj, &pushConjugFlag, &separationFlag, ccur, meta );

        }

        if ( (connCmp_nSols(ccur)>0) && separationFlag && widthFlag && compactFlag && (connCmp_nSols(ccur)<cacheApp_getDegree(cache))) {
            metadatas_add_validated( meta, depth, connCmp_nSols(ccur) );
            
#ifdef CAUCHY_CERTIFIED 
        slong m = connCmp_nSols(ccur);
        
//         compDsk_inflate_realRat(twoCCDisk, ccDisk, two);
//         cauchyTest_res resCert = cauchyTest_deterministic_counting_combinatorial_with_isoRatio( compDsk_centerref(twoCCDisk),
//                                                                                                 compDsk_radiusref(twoCCDisk),
//                                                                                                 two, m, cache, cacheCau, prec, 
//                                                                                                 meta, depth);
        realRat_mul_si(compDsk_radiusref(twoCCDisk), compDsk_radiusref(ccDisk), m+2);
        cauchyTest_res resCert; 
        resCert.nbOfSol = 1;
        if (m>1) {
            resCert = cauchyTest_deterministic_counting_combinatorial_with_rinfrsup( compDsk_centerref(ccDisk),
                                                                                                compDsk_radiusref(ccDisk),
                                                                                                compDsk_radiusref(twoCCDisk),
                                                                                                m, cache, cacheCau, prec, 
                                                                                                meta, depth);
        }
        if (resCert.nbOfSol == -1) {
                        failure = 1;
                        printf("#FAILURE: CERTIFICATION FAILED: A DISK is NOT 2-ISOLATED \n");
                        continue;
                    }
                    
        /* else inflate the cc */
        realRat_set_si(mp2, m+2, 1);
        connCmp_infate_realRat_inPlace(ccur,  mp2, metadatas_initBref(meta) );
#endif
            
            connCmp_list_push(qResults, ccur);
            nbSolsInQResults += connCmp_nSols(ccur);
            
            if (metadatas_getVerbo(meta)>=level)
                printf("------validated with %d roots\n", connCmp_nSols(ccur));
            
//             printf("metadatas_useRealCoeffs(meta): %d, pushConjugFlag: %d\n", metadatas_useRealCoeffs(meta), pushConjugFlag);
            /* Real Coeff */
            if ((metadatas_useRealCoeffs(meta))&&(pushConjugFlag)){
                /*compute the complex conjugate*/
                metadatas_add_validated( meta, depth, connCmp_nSols(ccurConj) );
                connCmp_list_push(qResults, ccurConj);
                nbSolsInQResults += connCmp_nSols(ccur);
            }
        }
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && resNewton.nflag ) {
            
            if (metadatas_getVerbo(meta)>=level)
                printf("------push in working queue after newton success %d roots\n", connCmp_nSols(ccur));
            
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        
        else if ( (connCmp_nSols(ccur)>0) && separationFlag && (resNewton.nflag==0) && (fmpz_cmp_si(connCmp_nwSpdref(ccur),4)>0) ){
            connCmp_decrease_nwSpd(ccur);

            if (metadatas_getVerbo(meta)>=level)
                printf("------push in working queue after newton fail %d roots\n", connCmp_nSols(ccur));
            
            connCmp_list_insert_sorted(qMainLoop, ccur);
        }
        else {

            if (metadatas_getVerbo(meta)>=level)
                printf("------bisect and push in working queue\n");
            
            cauchy_bisect_connCmp( ltemp, ccur, discardedCcs, bDiscarded, cache, cacheCau, nbMaxSol, meta,1);
            while (!connCmp_list_is_empty(ltemp))
                connCmp_list_insert_sorted(qMainLoop, connCmp_list_pop(ltemp));
            connCmp_clear(ccur);
            ccluster_free(ccur);
        }
    }
    
    compBox_clear(componentBox);
    compDsk_clear(ccDisk);
    compDsk_clear(fourCCDisk);
    compDsk_clear(twoCCDisk);
    realRat_clear(two);
    realRat_clear(three);
    realRat_clear(four);
    realRat_clear(mp2);
    realRat_clear(neps);
    realRat_clear(threeWidth);
    connCmp_list_clear(ltemp);
    
    return failure;
}

int cauchy_algo_global( connCmp_list_t qResults, 
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
    int failure = cauchy_main_loop( qResults, bDiscarded,  qMainLoop, discardedCcs, eps, cache, cacheCau, meta);
    
    
//     realRat_clear(factor);
    connCmp_list_clear(qMainLoop);
//     connCmp_list_clear(qPrepLoop);
    connCmp_list_clear(discardedCcs);
    
//     chronos_toc_CclusAl(metadatas_chronref(meta));
    metadatas_add_time_CclusAl(meta, (double) (clock() - start));
    
    return failure;
}

int metadatas_cauchy_fprint(FILE * file, metadatas_t meta, const realRat_t eps, cacheApp_t cache, cacheCauchy_t cacheCau){
    int r=1;
//     int nbTaylorShifts  = metadatas_getNbTaylorsInT0Tests(meta) + metadatas_getNbTaylorsInTSTests(meta);
//     int nbTaylorShiftsR = metadatas_getNbTaylorsRepetedInT0Tests(meta) + metadatas_getNbTaylorsRepetedInTSTests(meta);
//     int nbGraeffe       = metadatas_getNbGraeffeInT0Tests(meta) + metadatas_getNbGraeffeInTSTests(meta);
//     int nbGraeffeR      = metadatas_getNbGraeffeRepetedInT0Tests(meta) + metadatas_getNbGraeffeRepetedInTSTests(meta);
    char temp[1000];
    
    if (metadatas_getVerbo(meta)>=1) {
    r = fprintf(file, "# -------------------Cauchy Root Finder: ------------------------------\n");
    r = fprintf(file, "# -------------------Input:    ----------------------------------------\n");
    if (metadatas_getVerbo(meta)>=2) {
    slong degree = cacheApp_getDegree( cache );
    if (cacheCauchy_evalFastref(cacheCau) == NULL) {
        if (cache->_from_poly) {
            slong bitsize = cacheApp_getBitsize (cache);
            r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10ld|\n", "exact coefficients", degree, bitsize);
        } else
            r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10c|\n", "approximated coefficients", degree, '?');
    }
    else
        r = fprintf(file, "#|pol: %-25s degree: %-10ld bitsize: %10c|\n", "evaluation procedure", degree, '?');
    }
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
//     r = fprintf(file, "# ---------------------------------------------------------------------\n");
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time in shift for FFT        :",         metadatas_get_time_Taylors(meta),    " " );
    r = fprintf(file, "# -------------------Compression into rigid discs---------------------\n");
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number for clus of  1 root:",       metadatas_getComp_nb_1(meta),  " " );
    r = fprintf(file, "#|%-39s %14d %14s|\n", "total number for clus of >1 root:",       metadatas_getComp_nb_p(meta),  " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in comp s1/m:",                metadatas_get_time_CompCen(meta),    " " );
#ifdef DEFLATION_TURAN
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR with Turan's theorem:", metadatas_get_time_CompRR2(meta),    " " );
#else
    r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR with doub. exp. sieve:", metadatas_get_time_CompRR1(meta),    " " );
#endif
//     r = fprintf(file, "#|%-39s %14f %14s|\n", "time spent in RR algo 3:",                metadatas_get_time_CompRR3(meta),    " " );
    r = fprintf(file, "#|%-39s %14f %14s|\n", "total time spent in compression:",        metadatas_get_time_CompTot(meta),    " " );
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
