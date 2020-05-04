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

#include "ccluster/ccluster.h"

void risolate_compBox_gnuplot(FILE * f, 
                     const compBox_t b){
    
    int nbdigits = 12;
    int prec = 53;
    
    realRat_t xinf, xsup, yinf, ysup;
    realRat_t factor;
    realApp_t xinfa, xsupa, yinfa, ysupa;
    
    realRat_init(factor);
    realRat_init(xinf);
    realRat_init(xsup);
    realRat_init(yinf);
    realRat_init(ysup);
    realApp_init(xinfa);
    realApp_init(xsupa);
    realApp_init(yinfa);
    realApp_init(ysupa);
    
    realRat_set_si(factor, 1, 2);
    realRat_mul(factor, factor, compBox_bwidthref(b));
    realRat_sub(xinf, compRat_realref(compBox_centerref(b)), factor);
    realRat_add(xsup, compRat_realref(compBox_centerref(b)), factor);
    realRat_sub(yinf, compRat_imagref(compBox_centerref(b)), factor);
    realRat_add(ysup, compRat_imagref(compBox_centerref(b)), factor);
    realApp_set_realRat(xinfa, xinf, prec);
    realApp_set_realRat(xsupa, xsup, prec);
    realApp_set_realRat(yinfa, yinf, prec);
    realApp_set_realRat(ysupa, ysup, prec);
    
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
                     
void risolate_connCmp_gnuplot(FILE * f, 
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
    
    connCmp_risolate_componentBox( containingBox, c, metadatas_initBref(meta));
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

void risolate_connCmp_list_gnuplot(FILE * f, 
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
        risolate_connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
    fprintf(f, "e\n");
//     fprintf(f, "%s", Bcommand3);
    
    it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        risolate_connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
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

void risolate_connCmp_list_gnuplot_drawSubdiv(FILE * f, 
                          const connCmp_list_t l, 
                          const compBox_list_t lb,
                          metadatas_t meta){
    
    char preamble[100] = "# Ccluster output for GNUPLOT\n#Pipe it to gnuplot!\n";
//     char Bcommand[100] = "set pointsize 0.3\nplot '-' title 'Computed clusters' with xyerrorbars\n";
    char Bcommand1[100] = "set pointsize 1\n";
    char Bcommand2[1000] = "plot '-' title 'Computed clusters' with circles lc rgb \"#008080\" fs transparent solid 0.15 noborder,\\\n";
    char Bcommand3[1000] = "     '-' u 1:2 title 'centers of clusters' pt 2 lc rgb \"#008080\"";
    char Ecommand[100] = "\npause mouse close\n";
    
    connCmp_list_iterator it;
    
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
    fprintf(f, "set xrange["); realApp_fprintn(f, xinfa, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, ":");realApp_fprintn(f, xsupa, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, "]\n");
    fprintf(f, "set yrange["); realApp_fprintn(f, yinfa, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, ":");realApp_fprintn(f, ysupa, nbdigits, ARB_STR_NO_RADIUS);
    fprintf(f, "]\n");
    
//     fprintf(f, "%s", Bcommand);
    fprintf(f, "%s", Bcommand1);
    fprintf(f, "%s", Bcommand2);
    fprintf(f, "%s", Bcommand3);
    fprintf(f, ",\\\n     '-' title 'initial box' with lines lw 2 lc rgb \"black\"");
    
    /* iterate on cc in qres */
    it = connCmp_list_begin(l);
    while (it!=connCmp_list_end() ) {
        compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( connCmp_list_elmt(it)) );
        while (itb!= compBox_list_end() ){
            fprintf(f, ",\\\n     '-' title '' with lines lw 2 lc rgb \"#008080\"");
            itb = compBox_list_next( itb );
        }
        it = connCmp_list_next(it);
    }
    
    for (int s=0;s<compBox_list_get_size(lb); s++){
        fprintf(f, ",\\\n     '-' title '' with lines lw 1 lc rgb \"black\"");
    }
    
    fprintf(f, "\n");
    
    /* disks */
    it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        risolate_connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
    fprintf(f, "e\n");
//     fprintf(f, "%s", Bcommand3);
    
    /* centers */
    it = connCmp_list_begin(l);
    
    while (it!=connCmp_list_end() ) {
        risolate_connCmp_gnuplot(f, connCmp_list_elmt(it), meta);
        it = connCmp_list_next(it);
        fprintf(f, "\n");
    }
    fprintf(f, "e\n");
    
    /* initial box */
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
    
    fprintf(f, "e\n");
    
    /* boxes in res */
    it = connCmp_list_begin(l);
    while (it!=connCmp_list_end() ) {
        compBox_list_iterator itb = compBox_list_begin( connCmp_boxesref( connCmp_list_elmt(it)) );
        while (itb!= compBox_list_end() ){
            risolate_compBox_gnuplot(f, compBox_list_elmt(itb));
            itb = compBox_list_next( itb );
            fprintf(f, "e\n");
        }
        it = connCmp_list_next(it);
    }
    
    /* excluded boxes */
    connCmp_list_iterator itb = compBox_list_begin(lb);
    
    while (itb!=compBox_list_end() ) {
        risolate_compBox_gnuplot(f, compBox_list_elmt(itb));
        itb = compBox_list_next(itb);
        fprintf(f, "e\n");
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
