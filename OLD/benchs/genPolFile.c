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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "polynomials/compRat_poly.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/app_rat_poly.h"
#include "ccluster/ccluster.h"

void genBernPolFile( FILE * file, slong degree);
void genMignPolFile( FILE * file, slong degree, slong bitsize);
void genMignMulPolFile( FILE * file, slong degree, slong bitsize, slong power);
void genWilkPolFile( FILE * file, slong degree);
void genWilkRatPolFile( FILE * file, slong degree);
void genWilkMulPolFile( FILE * file, slong degree);
void genSpiralPolFile( FILE * file, slong degree, slong prec);
void clustersIterate( compApp_poly_ptr tabres, compApp_poly_ptr tabprec, int i, slong prec);
void genClusterPolFile( FILE * file, int iterations, slong prec);






int main(int argc, char **argv){
    
    if (argc<3){
        printf("usage: %s Bernoulli   degree\n", argv[0]);
        printf("usage: %s Mignotte    degree bitsize\n", argv[0]);
        printf("usage: %s MignotteMul degree bitsize power\n", argv[0]);
        printf("usage: %s Wilkinson   degree\n", argv[0]);
        printf("usage: %s WilkRat     degree\n", argv[0]);
        printf("usage: %s WilkMul     degree\n", argv[0]);
        printf("usage: %s Spiral      degree prec\n", argv[0]);
        printf("usage: %s Cluster     degree prec\n", argv[0]);
        return -1;
    }
    
    char poly[100];
//     char filename[100];
    char bernoulli[] = "Bernoulli\0";
    char mignotte[] = "Mignotte\0";
    char mignotteMul[] = "MignotteMul\0";
    char wilkinson[] = "Wilkinson\0";
    char wilkRat[] = "WilkRat\0";
    char wilkMul[] = "WilkMul\0";
    char spiral[] = "Spiral\0";
    char cluster[] = "Cluster\0";
    int degree = 0;
    int thirdArg = 0;
    int fourthArg = 0;
//     FILE * pFile;
    
    sscanf(argv[1], "%s", poly);
    sscanf(argv[2], "%d", &degree);
    if (argc>=4)
        sscanf(argv[3], "%d", &thirdArg);
    
    if (argc>=5)
        sscanf(argv[4], "%d", &fourthArg);
    
    if (strcmp(poly, bernoulli)==0) {
        genBernPolFile(stdout, degree);
    }
    
    if (strcmp(poly, wilkinson)==0) {
        genWilkPolFile(stdout, degree);
    }
    
    if (strcmp(poly, wilkRat)==0) {
        genWilkRatPolFile(stdout, degree);
    }
    
    if (strcmp(poly, wilkMul)==0) {
        genWilkMulPolFile(stdout, degree);
    }
    
    if (strcmp(poly, mignotte)==0) {
        if (argc<4)
            printf("usage: %s Mignotte degree bitsize\n", argv[0]);
        else{
            genMignPolFile(stdout, degree, thirdArg);
        }
    }
    
    if (strcmp(poly, mignotteMul)==0) {
        if (argc<5)
            printf("usage: %s Mignotte degree bitsize power\n", argv[0]);
        else{
            genMignMulPolFile(stdout, degree, thirdArg, fourthArg);
        }
    }
    
    if (strcmp(poly, spiral)==0) {
        if (argc<4)
            printf("usage: %s Spiral degree prec\n", argv[0]);
        else{
            genSpiralPolFile(stdout, degree, thirdArg);
        }
    }
    
    if (strcmp(poly, cluster)==0) {
        if (argc<4)
            printf("usage: %s Cluster degree prec\n", argv[0]);
        else{
            genClusterPolFile(stdout, degree, thirdArg);
        }
    }
    
}

void genBernPolFile( FILE * file, slong degree){
    
    realRat_poly_t pbern;
    realRat_poly_init(pbern);
    realRat_t coeff;
    realRat_init(coeff);
    
    bernoulli_polynomial( pbern, degree);
    
    fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Rational;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "\n");
    
    for(int i = degree; i>=0; i--){
        realRat_poly_get_coeff_realRat(coeff, pbern, i);
        if (!realRat_is_zero(coeff)) {
            fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
        }
    }
    realRat_poly_clear(pbern);
    realRat_clear(coeff);
}

void genMignPolFile( FILE * file, slong degree, slong bitsize){
    
    realRat_poly_t pmign;
    realRat_poly_init(pmign);
    realRat_t coeff;
    realRat_init(coeff);
    
    mignotte_polynomial(pmign, degree, bitsize);
    
    fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Rational;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "\n");
    
    for(int i = degree; i>=0; i--){
        realRat_poly_get_coeff_realRat(coeff, pmign, i);
        if (!realRat_is_zero(coeff)) {
            fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
        }
    }
    realRat_poly_clear(pmign);
    realRat_clear(coeff);
}

void genMignMulPolFile( FILE * file, slong degree, slong bitsize, slong power){
    realRat_poly_t pmign;
    realRat_poly_init(pmign);
    realRat_t coeff;
    realRat_init(coeff);
    
    mignotte_polynomial(pmign, degree, bitsize);
    realRat_poly_pow(pmign, pmign, (ulong) power);
    
    fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Rational;\n");
    fprintf(file, "Degree = %d;\n", (int) (degree*power));
    fprintf(file, "\n");
    
    for(int i = degree*power; i>=0; i--){
        realRat_poly_get_coeff_realRat(coeff, pmign, i);
        if (!realRat_is_zero(coeff)) {
            fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
        }
    }
    realRat_poly_clear(pmign);
    realRat_clear(coeff);
}

void genWilkPolFile( FILE * file, slong degree){
    
    realRat_poly_t pwilk, ptemp;
    realRat_poly_init(pwilk);
    realRat_poly_init2(ptemp,2);
    realRat_poly_one(pwilk);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(pwilk, pwilk, ptemp);
    }
    
    realRat_t coeff;
    realRat_init(coeff);
    
    fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Integer;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "\n");
    
    for(int i = degree; i>=0; i--){
        realRat_poly_get_coeff_realRat(coeff, pwilk, i);
        if (!realRat_is_zero(coeff)) {
            fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
        }
    }
    realRat_poly_clear(pwilk);
    realRat_poly_clear(ptemp);
    realRat_clear(coeff);
}

void genWilkRatPolFile( FILE * file, slong degree){
    
    realRat_poly_t pwilk, ptemp;
    realRat_poly_init(pwilk);
    realRat_poly_init2(ptemp,2);
    realRat_poly_one(pwilk);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, ((ulong) degree)+1);
        realRat_poly_mul(pwilk, pwilk, ptemp);
    }
    
    realRat_t coeff;
    realRat_init(coeff);
    
    fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Rational;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "\n");
    
    for(int i = degree; i>=0; i--){
        realRat_poly_get_coeff_realRat(coeff, pwilk, i);
        if (!realRat_is_zero(coeff)) {
            fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
        }
    }
    realRat_poly_clear(pwilk);
    realRat_poly_clear(ptemp);
    realRat_clear(coeff);
}

void genWilkMulPolFile( FILE * file, slong degree){
    
    realRat_poly_t pwilk, ptemp;
    realRat_poly_init(pwilk);
    realRat_poly_init2(ptemp,degree+1);
    realRat_poly_one(pwilk);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_zero(ptemp);
        realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_pow(ptemp, ptemp, (ulong) i);
//         realRat_poly_pow(ptemp, ptemp, (ulong) (degree-i)+1);
        realRat_poly_mul(pwilk, pwilk, ptemp);
    }
    
    realRat_t coeff;
    realRat_init(coeff);
    
    int truedeg = (int) ((degree*degree + degree)/2);
    
//     fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Real;\n");
    fprintf(file, "Integer;\n");
    fprintf(file, "Degree = %d;\n", (int) truedeg);
    fprintf(file, "\n");
    
    for(int i = 0; i<=truedeg; i++){
        realRat_poly_get_coeff_realRat(coeff, pwilk, i);
//         if (!realRat_is_zero(coeff)) {
//             fprintf(file, "%d ", i);
            realRat_fprint(file, coeff);
            fprintf(file, "\n");
//         }
    }
    realRat_poly_clear(pwilk);
    realRat_poly_clear(ptemp);
    realRat_clear(coeff);
}

void genSpiralPolFile( FILE * file, slong degree, slong prec){
    
    compApp_poly_t dest;
    compApp_poly_init(dest);
    
    realRat_t modu;
    realRat_t argu;
    compApp_t a_modu;
    compApp_t a_argu;
    compApp_t coeff;
    
    realRat_init(modu);
    realRat_init(argu);
    compApp_init(a_modu);
    compApp_init(a_argu);
    compApp_init(coeff);
    
    compApp_poly_t temp;
    compApp_poly_init2(temp,2);
    compApp_poly_set_coeff_si(temp, 1, 1);
    
    compApp_poly_one(dest);
    slong prectemp = 10000;
    
    for(int i=1; i<=degree; i++){
        realRat_set_si(modu, -i, (ulong) degree);
        realRat_set_si(argu, 4*i, (ulong) degree);
        compApp_set_realRat( a_modu, modu, prectemp);
        compApp_set_realRat( a_argu, argu, prectemp);
        compApp_exp_pi_i( coeff, a_argu, prectemp);
        compApp_mul( coeff, coeff, a_modu, prectemp);
        compApp_poly_set_coeff_compApp(temp, 0, coeff);
        compApp_poly_mul(dest, dest, temp, prectemp);
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
    }
    
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_poly_clear(temp);
    
    
//     fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Complex;\n");
    fprintf(file, "FloatingPoint;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "Precision = %d;\n", (int) prec);
    fprintf(file, "\n");
    
    char tempstr[100*prec];
    char * temp2;
    
    for(int i = 0; i<=degree; i++){
        compApp_set(coeff, compApp_poly_getCoeff(dest, i));
        
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
        
        temp2 = arb_get_str(compApp_realref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s ", tempstr);
        temp2 = arb_get_str(compApp_imagref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s\n", tempstr);
        
    }

    compApp_clear(coeff);    
    compApp_poly_clear(dest);
}

void clustersIterate( compApp_poly_ptr tabres, compApp_poly_ptr tabprec, int i, slong prec){
    // tabres is a table of 3^i compApp_poly
    // tabres is a table of 3^(i-1) compApp_poly
    realRat_t modu;
    realRat_t argu;
    compApp_t a_modu;
    compApp_t a_argu;
    compApp_t coeff;
    
    realRat_init(modu);
    realRat_init(argu);
    compApp_init(a_modu);
    compApp_init(a_argu);
    compApp_init(coeff);
    int indexInTabRes = 0;
//     printf("pow(3,i-1): %d\n", (int) pow(3,i-1));
    for (int j = 0; j<((int) pow(3,i-1)); j++){
//         printf("%d\n", j);
//         realRat_set_si(modu, -1, (ulong) 0x1<<(4*(i-1)));
        realRat_set_si(modu, -1, (ulong) pow(4, 2*(i-1)));
        compApp_set_realRat( a_modu, modu, prec);
        
        realRat_set_si(argu, 2, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 4, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 6, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
    }
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_clear(coeff);    
}

void genClusterPolFile( FILE * file, int iterations, slong prec){
    
    compApp_poly_ptr tabprec;
    tabprec = (compApp_poly_ptr) malloc (1*sizeof(compApp_poly));
    compApp_poly_init2( tabprec, 2);
    compApp_poly_zero(tabprec);
    compApp_poly_set_coeff_si(tabprec, 1, 1);
    int degree = 3;
    
    slong prectemp = 10000;
    
    for (int i=1; i<=iterations; i++) {
        compApp_poly_ptr tabres;
        tabres = (compApp_poly_ptr) malloc (degree*sizeof(compApp_poly));
        for(int j = 0; j<degree; j++) compApp_poly_init2( tabres + j, 2);
        
        clustersIterate( tabres, tabprec, i, prectemp);
        for(int j = 0; j<((int) (degree/3)); j++) compApp_poly_clear( tabprec + j);
        free(tabprec);
        tabprec = tabres;
        degree = degree*3;
    }
    
    degree = (int) degree/3;
    compApp_poly_t dest;
    compApp_poly_init2(dest, degree+1);
    compApp_poly_one(dest);
    for(int j = 0; j<degree; j++) compApp_poly_mul(dest, dest, tabprec +j, prectemp);
    for(int j = 0; j<degree; j++) compApp_poly_clear( tabprec + j);
    free(tabprec);
    
    compApp_t coeff;
    compApp_init(coeff);
    
//     fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Complex;\n");
    fprintf(file, "FloatingPoint;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "Precision = %d;\n", (int) prec);
    fprintf(file, "\n");
    
    char tempstr[100*prec];
    char * temp2;
    
    for(int i = 0; i<=degree; i++){
        compApp_set(coeff, compApp_poly_getCoeff(dest, i));
        
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
        temp2 = arb_get_str(compApp_realref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s ", tempstr);
        temp2 = arb_get_str(compApp_imagref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s\n", tempstr);
        
    }

    compApp_clear(coeff);    
    compApp_poly_clear(dest);
}