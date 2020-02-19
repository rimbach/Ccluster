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

#include "ccluster/ccluster.h"
    
void ccluster_interface_func( void(*func)(compApp_poly_t, slong), 
                              const compBox_t initialBox, 
                              const realRat_t eps, 
                              char * stratstr,
                              int nbThreads,
                              int output,
                              int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    connCmp_list_init(qRes);
    
    ccluster_algo( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 1);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
}

void ccluster_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                     const realRat_t eps, 
                                     char * stratstr,
                                     int nbThreads,
                                     int output,
                                     int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automaticly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    connCmp_list_init(qRes);
    
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    ccluster_algo_global( qRes, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (output==-2) {
//         printf("gnuplot output: not yet implemented\n");
        connCmp_list_gnuplot(stdout, qRes, meta, 0);
    } else if (output!=0) {
//         printf("cluster output: not yet implemented\n");
        connCmp_list_print_for_results_withOutput(stdout, qRes, output, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_clear(initialBox);
}

void connCmp_gnuplot(FILE * f, 
                     const connCmp_t c, 
                     metadatas_t meta){
    
    compBox_t containingBox;
    compBox_init(containingBox);
    compDsk_t containingDisk;
    compDsk_init(containingDisk);
    realApp_t cRe, cIm, rad;
    realApp_init(cRe);
    realApp_init(cIm);
    realApp_init(rad);
    
    connCmp_componentBox( containingBox, c, metadatas_initBref(meta));
    compBox_get_containing_dsk( containingDisk, containingBox);
    
    slong l = fmpz_clog_ui( realRat_denref(compDsk_radiusref(containingDisk)), (ulong) 2);
//     printf("l: %ld\n", l);
    int prec = ( 53 > l ? 53 : l );
    int nbdigits = (int) ceil( prec/4 ) ;
    
    realApp_set_realRat(cRe, compRat_realref(compDsk_centerref(containingDisk)), prec);
    realApp_set_realRat(cIm, compRat_imagref(compDsk_centerref(containingDisk)), prec);
    realApp_set_realRat(rad, compDsk_radiusref(containingDisk), prec);
    
    realApp_fprintn(f, cRe, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, "   ");
    realApp_fprintn(f, cIm, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, "   ");
    realApp_fprintn(f, rad, nbdigits, ARB_STR_NO_RADIUS);
    
    realApp_clear(cRe);
    realApp_clear(cIm);
    realApp_clear(rad);
    compBox_clear(containingBox);
    compDsk_clear(containingDisk);
}

void connCmp_list_gnuplot(FILE * f, 
                          const connCmp_list_t l, 
                          metadatas_t meta,
                          int withInitBox){
    
    char preamble[100] = "# Ccluster output for GNUPLOT\n#Pipe it to gnuplot!\n";
//     char Bcommand[100] = "set pointsize 0.3\nplot '-' title 'Computed clusters' with xyerrorbars\n";
    char Bcommand1[100] = "set pointsize 1\n";
    char Bcommand2[1000] = "plot '-' title 'Computed clusters' with circles lc rgb \"#008080\" fs transparent solid 0.15 noborder,\\\n";
    char Bcommand3[1000] = "     '-' u 1:2 title 'centers of clusters' pt 2 lc rgb \"#008080\"";
    char Ecommand[100] = "e\npause mouse close\n";
    
    /*set X,Y ranges to (5/4) initbox*/
    realRat_t xinf, xsup, yinf, ysup;
    realRat_t factor;
    realApp_t xinfa, xsupa, yinfa, ysupa;
    int nbdigits = 12;
    int prec = 53;
    
    realRat_init(factor);
    realRat_init(xinf);
    realRat_init(xsup);
    realRat_init(yinf);
    realRat_init(ysup);
    realApp_init(xinfa);
    realApp_init(xsupa);
    realApp_init(yinfa);
    realApp_init(ysupa);
    
    realRat_set_si(factor, 5, 8);
    realRat_mul(factor, factor, compBox_bwidthref(metadatas_initBref(meta)));
    realRat_sub(xinf, compRat_realref(compBox_centerref(metadatas_initBref(meta))), factor);
    realRat_add(xsup, compRat_realref(compBox_centerref(metadatas_initBref(meta))), factor);
    realRat_sub(yinf, compRat_imagref(compBox_centerref(metadatas_initBref(meta))), factor);
    realRat_add(ysup, compRat_imagref(compBox_centerref(metadatas_initBref(meta))), factor);
    realApp_set_realRat(xinfa, xinf, prec);
    realApp_set_realRat(xsupa, xsup, prec);
    realApp_set_realRat(yinfa, yinf, prec);
    realApp_set_realRat(ysupa, ysup, prec);
    

    fprintf(f, "%s", preamble);
    if (withInitBox) {
        fprintf(f, "set xrange["); realApp_fprintn(f, xinfa, nbdigits, ARB_STR_NO_RADIUS);
        fprintf(f, ":");realApp_fprintn(f, xsupa, nbdigits, ARB_STR_NO_RADIUS);
        fprintf(f, "]\n");
        fprintf(f, "set yrange["); realApp_fprintn(f, yinfa, nbdigits, ARB_STR_NO_RADIUS);
        fprintf(f, ":");realApp_fprintn(f, ysupa, nbdigits, ARB_STR_NO_RADIUS);
        fprintf(f, "]\n");
    }
//     fprintf(f, "%s", Bcommand);
    fprintf(f, "%s", Bcommand1);
    fprintf(f, "%s", Bcommand2);
    fprintf(f, "%s", Bcommand3);
    if (withInitBox) {
        fprintf(f, ",\\\n     '-' title 'initial box' with lines lw 2 lc rgb \"black\"");
    }
    fprintf(f, "\n");
    
    connCmp_list_iterator it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
    fprintf(f, "e\n");
//     fprintf(f, "%s", Bcommand3);
    
    it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
    if (withInitBox) {
        realRat_set_si(factor, 1, 2);
        realRat_mul(factor, factor, compBox_bwidthref(metadatas_initBref(meta)));
        realRat_sub(xinf, compRat_realref(compBox_centerref(metadatas_initBref(meta))), factor);
        realRat_add(xsup, compRat_realref(compBox_centerref(metadatas_initBref(meta))), factor);
        realRat_sub(yinf, compRat_imagref(compBox_centerref(metadatas_initBref(meta))), factor);
        realRat_add(ysup, compRat_imagref(compBox_centerref(metadatas_initBref(meta))), factor);
        realApp_set_realRat(xinfa, xinf, prec);
        realApp_set_realRat(xsupa, xsup, prec);
        realApp_set_realRat(yinfa, yinf, prec);
        realApp_set_realRat(ysupa, ysup, prec);
        fprintf(f, "e\n");
        realApp_fprintn(f, xinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "   ");
        realApp_fprintn(f, yinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "\n");
        realApp_fprintn(f, xsupa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "   ");
        realApp_fprintn(f, yinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "\n");
        realApp_fprintn(f, xsupa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "   ");
        realApp_fprintn(f, ysupa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "\n");
        realApp_fprintn(f, xinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "   ");
        realApp_fprintn(f, ysupa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "\n");
        realApp_fprintn(f, xinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "   ");
        realApp_fprintn(f, yinfa, nbdigits, ARB_STR_NO_RADIUS); fprintf(f, "\n");
    }
    
    
    fprintf(f, "%s", Ecommand);
    
    realRat_clear(factor);
    realRat_clear(xinf);
    realRat_clear(xsup);
    realRat_clear(yinf);
    realRat_clear(ysup);
    realApp_clear(xinfa);
    realApp_clear(xsupa);
    realApp_clear(yinfa);
    realApp_clear(ysupa);
}

/* version with function for fast evaluation */
void ccluster_interface_func_eval( void(*func)(compApp_poly_t, slong),
                                   void(*evalFast)(compApp_t, compApp_t, const compApp_t, slong),
                                   const compBox_t initialBox, 
                                   const realRat_t eps, 
                                   char * stratstr,
                                   int nbThreads,
                                   int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
 
    connCmp_list_init(qRes);
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, evalFast, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    ccluster_algo( qRes, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
}

/* experimental version */
void ccluster_expe_global_interface_func( void(*func)(compApp_poly_t, slong), 
                                          const realRat_t eps, 
                                          char * stratstr,
                                          int nbThreads,
                                          int verb){

    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init(cache, func);
    strategies_init(strat);
    
    /* automaticly set initialBox */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_si(initialBox, 0,1,0,1,0,1);
    cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
    if (verb>=3) {
        printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
    }
    realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
    
    strategies_set_str( strat, stratstr, nbThreads );
    /* automaticly set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    
    connCmp_list_init(qRes);
    strategies_set_powerSums( strat, 1 );
    metadatas_init(meta, initialBox, strat, verb);
    
    /* initialize power sums */
//     metadatas_setNbPowerSums(meta, 2);
//     metadatas_setIsoRatio_si(meta, 2, 1);
// //     metadatas_setIsoRatio_si(meta, 4, 3);
// //     if ( metadatas_pwSuTest(meta) )
//     ccluster_initialize_pwSuTest(NULL, meta, cache, verb);
//     if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 2, verb );
    
    ccluster_expe_algo_global( qRes, initialBox, eps, cache, meta);
    
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
    }
    
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    compBox_clear(initialBox);
}

// void ccluster_forJulia_func( connCmp_list_t qResults, 
//                              void(*func)(compApp_poly_t, slong), 
//                              const compBox_t initialBox, 
//                              const realRat_t eps, 
//                              char * stratstr,
//                              int nbThreads,
//                              int verb){
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
//     
//     cacheApp_t cache;
//     strategies_t strat;
//     metadatas_t meta;
//     
//     cacheApp_init(cache, func);
//     strategies_init(strat);
//     
//     strategies_set_str( strat, stratstr, nbThreads );
//     /* automatically set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
//     
//     metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//     
//     ccluster_algo( qResults, initialBox, eps, cache, meta);
//     
//     metadatas_count(meta);
//     metadatas_fprint(stdout, meta, eps);
//     if (verb>=3) {
//         connCmp_list_print_for_results(stdout, qResults, meta);
//     }
//     
//     cacheApp_clear(cache);
//     strategies_clear(strat);
//     metadatas_clear(meta);
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
// }

// void ccluster_global_forJulia_func( connCmp_list_t qResults, 
//                                     void(*func)(compApp_poly_t, slong),  
//                                     const realRat_t eps, 
//                                     char * stratstr,
//                                     int nbThreads,
//                                     int verb){
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: begin\n");
//     
//     cacheApp_t cache;
//     strategies_t strat;
//     metadatas_t meta;
//     
//     cacheApp_init(cache, func);
//     strategies_init(strat);
//     
//     /* automaticly set initialBox */
// //     compBox_t initialBox;
// //     compBox_init(initialBox);
//     compBox_set_si(initialBox, 0,1,0,1,0,1);
//     cacheApp_root_bound ( compBox_bwidthref(initialBox), cache );
//     if (verb>=3) {
//         printf("root bound: "); realRat_print(compBox_bwidthref(initialBox)); printf("\n");
//     }
//     realRat_mul_si(compBox_bwidthref(initialBox), compBox_bwidthref(initialBox), 2);
//     
//     strategies_set_str( strat, stratstr, nbThreads );
//     /* automatically set realCoeffs */
//     if (cacheApp_is_real(cache)==0
//         || compBox_contains_real_line_in_interior(initialBox)==0 )
//         strategies_set_realCoeffs(strat, 0);
//     
//     metadatas_init(meta, initialBox, strat, verb);
//     /* initialize power sums */
//     if (metadatas_usePowerSums(meta))
//         metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
//     
//     ccluster_algo_global( qResults, initialBox, eps, cache, meta);
//     
//     metadatas_count(meta);
//     metadatas_fprint(stdout, meta, eps);
//     if (verb>=3) {
//         connCmp_list_print_for_results(stdout, qResults, meta);
//     }
//     
//     cacheApp_clear(cache);
//     strategies_clear(strat);
//     metadatas_clear(meta);
// //     compBox_clear(initialBox);
//     
// //     printf("ccluster.c: ccluster_interface_forJulia_func: end\n");
// }

int ccluster_interface_poly( realRat_t * centerRe, realRat_t * centerIm, int * mults, 
                             const compRat_poly_t poly, 
                             const compBox_t initialBox, 
                             const realRat_t eps, 
                             int st, 
                             int verb){
    
    cacheApp_t cache;
    strategies_t strat;
    metadatas_t meta;
    connCmp_list_t qRes;
    
    cacheApp_init_compRat_poly(cache, poly);
    strategies_init(strat);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st>>6);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4),0, (st&( ((0x1<<10)-1)<<5 ))>>5, st>>16);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), (st&( ((0x1<<10)-1)<<6 ))>>6, st>>17);
//     strategies_set_int ( strat, st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<3), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), (st&( ((0x1<<10)-1)<<7 ))>>7, st>>18);
    strategies_set_int ( strat, 
                         st&(0x1), st&(0x1<<1), st&(0x1<<2), st&(0x1<<4), st&(0x1<<5), st&(0x1<<6), (st&( ((0x1<<10)-1)<<7 ))>>7, st>>18);
    
    
    /* automatically set realCoeffs */
    if (cacheApp_is_real(cache)==0
        || compBox_contains_real_line_in_interior(initialBox)==0 )
        strategies_set_realCoeffs(strat, 0);
    metadatas_init(meta, initialBox, strat, verb);
    /* initialize power sums */
    if (metadatas_usePowerSums(meta))
        metadatas_set_pwSuDatas( meta, NULL, cacheApp_getDegree(cache), 2, 1, 1, verb );
    
    connCmp_list_init(qRes);
    
    ccluster_algo( qRes, initialBox, eps, cache, meta);
    metadatas_count(meta);
    metadatas_fprint(stdout, meta, eps);
    
    if (verb>=3) {
        connCmp_list_print_for_results(stdout, qRes, meta);
//         connCmp_list_print_for_results(stdout, qRes, 500, 40, meta);
    }
    
    /* feed the results */
//     int nbClus = connCmp_list_get_size(qRes);
    int nbClus = 0;
    compBox_t containingBox;
    compBox_init(containingBox);
    connCmp_list_iterator it = connCmp_list_begin(qRes);
    while (it!=connCmp_list_end() ) {
        
        connCmp_componentBox( containingBox, connCmp_list_elmt(it), metadatas_initBref(meta));
        realRat_set( centerRe[nbClus], compRat_realref(compBox_centerref(containingBox)) );
        realRat_set( centerIm[nbClus], compRat_imagref(compBox_centerref(containingBox)) );
        mults[nbClus] = connCmp_nSols(connCmp_list_elmt(it));
        
        it = connCmp_list_next(it);
        nbClus++;
    }
    compBox_clear(containingBox);
       
    cacheApp_clear(cache);
    strategies_clear(strat);
    metadatas_clear(meta);
    connCmp_list_clear(qRes);
    
    return nbClus;
}

int ccluster_interface_poly_real( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                  const realRat_poly_t poly, 
                                  const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                  const realRat_t eps, 
                                  int st, 
                                  int verb){
    
    /* initial Box */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_3realRat(initialBox, initialBox_cr, initialBox_ci, initialBox_wi);
    /* polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set_realRat_poly(p,poly);
    
    /* call */
    int nbClus = ccluster_interface_poly( centerRe, centerIm, mults, p, initialBox, eps, st, verb);
    
    /* clear */
    compBox_clear(initialBox);
    compRat_poly_clear(p);
    
    return nbClus;
    
}

int ccluster_interface_poly_real_imag( realRat_t * centerRe, realRat_t * centerIm, int * mults,
                                       const realRat_poly_t poly_real, const realRat_poly_t poly_imag, 
                                       const realRat_t initialBox_cr, const realRat_t initialBox_ci, const realRat_t initialBox_wi,
                                       const realRat_t eps, 
                                       int st, 
                                       int verb){
    
    /* initial Box */
    compBox_t initialBox;
    compBox_init(initialBox);
    compBox_set_3realRat(initialBox, initialBox_cr, initialBox_ci, initialBox_wi);
    /* polynomial */
    compRat_poly_t p;
    compRat_poly_init(p);
    compRat_poly_set2_realRat_poly(p, poly_real, poly_imag);
    
    /* call */
    int nbClus = ccluster_interface_poly( centerRe, centerIm, mults, p, initialBox, eps, st, verb);
    
    /* clear */
    compBox_clear(initialBox);
    compRat_poly_clear(p);
    
    return nbClus;
    
}
