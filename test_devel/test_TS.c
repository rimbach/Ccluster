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

#include <stdio.h>
#include "flint/flint.h"
#include "polynomials/compApp_poly.h"
#include "polynomials/compRat_poly.h"
#include "polynomials/app_rat_poly.h"
#include <time.h>

void realRat_poly_wilkinson( realRat_poly_t res, slong degree );

void realRat_poly_wilkRat( realRat_poly_t res, slong degree );

void realRat_poly_wilkFac( realRat_poly_t res, slong degree );

void realRat_poly_wilkFacRat( realRat_poly_t res, slong degree );

void regularGrid_polynomial(realRat_poly_t preg, int res);

void realRat_poly_ones( realRat_poly_t res, slong degree );

void compApp_poly_div_si ( compApp_poly_t res, const compApp_poly_t p, slong y, slong prec);

void precompute( compApp_poly_ptr tab, const compApp_poly_t p, slong prec );

void compApp_poly_taylorShift_home( compApp_poly_t res, compApp_poly_ptr tab, const compApp_t center, slong prec);

void clear_table( compApp_poly_ptr tab );



fmpz ** precompute_binomial_table( fmpz ** binCoeffsTable, slong degree );

void clear_binomial_table( fmpz ** binCoeffsTable, slong degree );

/* powers already initialized */
void compute_powers( compApp_ptr powers, const compRat_t center, slong len, slong prec){
    
//     printf("center: \n"); compRat_print(center); printf("\n\n");
    
    compApp_one( powers+0 );
    compApp_set_compRat( powers+1, center, prec );
//     printf("powers+1: \n"); compApp_printd(powers+1,10); printf("\n\n");
    for (slong i = 2; i<len; i++) {
        compApp_mul_compRat( powers+i, powers + (i-1), center, prec );
//         printf("powers+i: \n"); compApp_printd(powers+i,10); printf("\n\n");
//         printf("center: \n"); compRat_print(center); printf("\n\n");
    }
    
//     for (slong i = 0; i<len; i++){
//         compApp_printd(powers+i, 10);
//         printf(" || ");
//     }
//     printf("\n");
}

// void compute_powersRat( compRat_ptr powers, const compRat_t center, slong len){
//     
// //     printf("center: \n"); compRat_print(center); printf("\n\n");
//     
//     compRat_one( powers+0 );
//     compRat_set( powers+1, center );
// //     printf("powers+1: \n"); compApp_printd(powers+1,10); printf("\n\n");
//     for (slong i = 2; i<len; i++) {
//         compRat_mul( powers+i, powers + (i-1), center);
// //         printf("powers+i: \n"); compApp_printd(powers+i,10); printf("\n\n");
// //         printf("center: \n"); compRat_print(center); printf("\n\n");
//     }
//     
//     for (slong i = 0; i<len; i++){
//         compRat_print(powers+i);
//         printf(" || ");
//     }
//     printf("\n");
// }

void compute_column( compApp_ptr column, 
                     compApp_ptr powers, 
                     fmpz ** binCoeffsTable,
                     const compApp_poly_t p,
                     slong ind_in_block,
                     slong ind_of_block,
                     slong len_of_block,
                     slong degree,
                     slong prec ) {
    
    slong index_of_column = ind_of_block*len_of_block + ind_in_block;
//     printf("index_of_column: %ld\n", index_of_column);
    
    if ( compApp_is_zero( p->coeffs + index_of_column ) ) {
//         printf("index_of_column: %ld; coeff is zero\n", index_of_column);
        return;
    }
    
    compApp_t temp, temp2;
    compApp_init(temp);
    compApp_init(temp2);
    
    slong i = 0;
    slong sym = index_of_column - i;
    
    compApp_mul( temp, powers + ind_in_block, p->coeffs + index_of_column, prec);
    /* first element in column: (p->coeffs + index_of_column)*x^(ind_in_block) */
    compApp_add( column+i, column+i, temp, prec);
    /* last element in column: (p->coeffs + index_of_column)*x^0 */
    if (index_of_column>0) /*not when i==sym==0*/
        compApp_add( column+sym, column+sym, p->coeffs + index_of_column, prec);
    
    slong ind_in_block_save = ind_in_block;
    slong ind_of_block_save = ind_of_block;
    ind_in_block=len_of_block-2;
    i = 1;
    sym = index_of_column - i;
    
    slong stop = (slong) floor(index_of_column/2);
    
    while (ind_of_block>0 && i <= stop) {    
        
        compApp_mul_fmpz( temp2, temp, binCoeffsTable[degree - index_of_column] + i-1, prec);
        compApp_add( column+i, column+i, temp2, prec);
        
        slong sym = index_of_column - i;
        if ( (sym <= (ind_of_block_save*len_of_block)) && (i!=sym) )
            compApp_add( column+sym, column+sym, temp2, prec);
        
//         printf("ind_of_block: %ld, ind_in_block: %ld  ", ind_of_block, ind_in_block);
// //         compApp_printd(column+i, 10);
//         printf("binomial coeff: "); fmpz_print(binCoeffsTable[i-1] + index_of_column-i);
//         if (i<= stop) {
//             printf(" index with same coeff: %ld", sym);
//             printf(" binomial coeff at that place: "); fmpz_print(binCoeffsTable[sym-1] + index_of_column-sym);
// //             if (i==sym) printf(" *");
//             if ( (sym <= (ind_of_block_save*len_of_block)) && (i!=sym) )
//                 printf(" *");
//         }
//         printf("\n");
        
        if (ind_in_block==0){
            ind_in_block = len_of_block - 1;
            ind_of_block--;
        }
        else
            ind_in_block--;
        i++;
    }

    i=ind_of_block_save*len_of_block+1;
    
    ind_in_block = ind_in_block_save-1;
    while (ind_in_block>=1){
//         printf("ind_in_block: %ld \n", ind_in_block);
        compApp_mul( temp, powers + ind_in_block, p->coeffs + index_of_column, prec);
        compApp_mul_fmpz( temp2, temp, binCoeffsTable[degree - index_of_column] + i-1, prec);
        compApp_add( column+i, column+i, temp2, prec);
        
//         printf("ind_of_block: %ld, ind_in_block: %ld  ", ind_of_block, ind_in_block);
//         compApp_printd(column+i, 10);
//         printf("\n");
        
        ind_in_block--;
        i++;
    }
    
    compApp_clear(temp);
    compApp_clear(temp2);
    
//     printf( "column %ld: ", index_of_column );
//     for ( slong i = 0; i < index_of_column+1; i++){
//         compApp_printd(column+i, 10);
//         printf(" || ");
//     }
//     printf("\n");
}

void add_block( compApp_poly_t pres,
                compApp_ptr column, 
                compApp_ptr powers, 
                fmpz ** binCoeffsTable,
                slong ind_in_block,
                slong ind_of_block,
                slong len_of_block,
                slong degree,
                slong prec ) {
    
    /* one block is computed; ind_in_block is 0 */
    slong index_of_column = ind_of_block*len_of_block + ind_in_block;
//     printf("index_of_column: %ld\n", index_of_column);
    for ( slong i = 0; i < index_of_column; i++){
       
        if (ind_of_block>0){
            /* multiply by power ind_of_block*len_of_block-i*/
//             printf("power of x: %ld\n", ind_of_block*len_of_block-i);
            compApp_mul( column + i, column + i, powers + (ind_of_block*len_of_block-i), prec);
        }
        /* index of coeff in pres: i */
//         printf("index of coeff in pres: %ld\n", i);
        compApp_add(pres->coeffs + i, pres->coeffs + i, column + i, prec);
        /* set column +i to zero */
        compApp_zero( column+i );
    }
    for ( slong i = 0; (i < len_of_block) && (ind_of_block*len_of_block + i < degree +1); i++){
        
        compApp_add(pres->coeffs + i + index_of_column, pres->coeffs + i + index_of_column, column + i + index_of_column, prec);
        compApp_zero( column+i + index_of_column);
        
    }
    
//     pres->length=degree+1;
//     printf("pres: \n"); compApp_poly_printd(pres, prec); printf("\n\n");
    
}

void my_Taylor_Shift( compApp_poly_t pres, 
                      const compApp_poly_t p, 
                      const compRat_t center,
                      fmpz ** binCoeffsTable,
                      slong degree, 
                      slong prec){
    
    compApp_poly_fit_length(pres, degree+1);
    
    slong len_of_block = (slong) floor(sqrt( (double) degree+1 ));
    slong num_of_block = (slong) (degree+1)/len_of_block;
    slong ind_in_block = (degree+1) - num_of_block*len_of_block;
    slong ind_of_block = num_of_block;
    slong r = degree+1-len_of_block +1;
    if (ind_in_block==0){
        ind_in_block = len_of_block-1;
        ind_of_block = ind_of_block-1;
    }
    else {
        r = num_of_block*len_of_block +1;
        ind_in_block = ind_in_block-1;
    }

//     printf("len_of_block: %ld, num_of_block: %ld, ind_in_block: %ld, r: %ld\n", len_of_block, num_of_block, ind_in_block, r);
    
    /*compute powers of the point*/
    compApp_ptr powers = (compApp_ptr) ccluster_malloc  (r*sizeof(compApp));
    for (slong i = 0; i<r; i++){
        compApp_init( powers+i );
    }
    compute_powers( powers, center, r, prec);
    
//     compRat_ptr powersRat = (compRat_ptr) ccluster_malloc  (r*sizeof(compRat));
//     for (slong i = 0; i<r; i++){
//         compRat_init( powersRat+i );
//     }
//     compute_powersRat( powersRat, center, r);
   
    /*column temp*/
    compApp_ptr column_temp = (compApp_ptr) ccluster_malloc  ((degree+1)*sizeof(compApp));
    for (slong i = 0; i<degree+1; i++){
        compApp_init( column_temp+i );
        compApp_zero( column_temp+i );
    }
    
    /* begin with highest degree coefficient */
    while ( (ind_of_block>=0)&&(ind_in_block>=0) ) {
//         printf("index_in_block: %ld, index_of_block: %ld\n", ind_in_block,ind_of_block);
        compute_column( column_temp, powers, binCoeffsTable, p, ind_in_block, ind_of_block, len_of_block, degree, prec );
        
        if (ind_in_block==0){
            add_block( pres, column_temp, powers, binCoeffsTable, ind_in_block, ind_of_block, len_of_block, degree, prec );
            ind_in_block = len_of_block -1;
            ind_of_block--;
        }
        else
            ind_in_block--;
    }
    
    /* clear powers */
    for (slong i = 0; i<r; i++)
        compApp_clear(powers+i);
    ccluster_free(powers);
    
//     /* clear powersRat */
//     for (slong i = 0; i<r; i++)
//         compRat_clear(powersRat+i);
//     ccluster_free(powersRat);
    
    /* clear column_temp */
    for (slong i = 0; i<degree+1; i++)
        compApp_clear(column_temp+i);
    ccluster_free(column_temp);
    
    pres->length=degree+1;
//     printf("pres: \n"); compApp_poly_printd(pres, 10); printf("\n\n");
    
}

compApp ** precompute_coeffsTable( compApp ** coeffsTable, fmpz ** binCoeffsTable, const compApp_poly_t p, slong degree, slong prec );

void clear_coeffsTable( compApp ** coeffsTable, slong degree );

void compute_column_index0( compApp_ptr column, 
                     compApp_ptr powers, 
                     compApp ** coeffsTable,
                     const compApp_poly_t p,
                     slong ind_of_block,
                     slong len_of_block,
                     slong degree,
                     slong prec );

void compute_column2( compApp_ptr column, 
                     compApp_ptr powers, 
                     compApp ** coeffsTable,
                     const compApp_poly_t p,
                     slong ind_in_block,
                     slong ind_of_block,
                     slong len_of_block,
                     slong degree,
                     slong prec ) {
    
    slong index_of_column = ind_of_block*len_of_block + ind_in_block;
//     printf("index_of_column: %ld\n", index_of_column);
    
    if ( compApp_is_zero( p->coeffs + index_of_column ) ) {
//         printf("index_of_column: %ld; coeff is zero\n", index_of_column);
        return;
    }
    
    if (ind_in_block==0) {
        compute_column_index0( column, powers, coeffsTable, p, ind_of_block, len_of_block, degree, prec );
        return;
    }
    
    compApp_t temp;
    compApp_init(temp);
    
    slong i = 0;
    slong sym = index_of_column - i;
    
    compApp_mul( temp, powers + ind_in_block, p->coeffs + index_of_column, prec);
    /* first element in column: (p->coeffs + index_of_column)*x^(ind_in_block) */
    compApp_add( column+i, column+i, temp, prec);
    /* last element in column: (p->coeffs + index_of_column)*x^0 */
    if (index_of_column>0) /*not when i==sym==0*/
        compApp_add( column+sym, column+sym, p->coeffs + index_of_column, prec);
    
    i = 1;
    sym = index_of_column - i;
    
    slong stop = (slong) floor(index_of_column/2);
    
    if (ind_of_block>0)
        while (i <= stop) {    
        
            compApp_mul( temp, powers + ind_in_block, coeffsTable[degree - index_of_column] + i-1, prec);
            compApp_add( column+i, column+i, temp, prec);
        
            sym = index_of_column - i;
            if ( (sym <= (ind_of_block*len_of_block)) && (i!=sym) )
                compApp_add( column+sym, column+sym, temp, prec);
        
            i++;
        }

    i=ind_of_block*len_of_block+1;
    ind_in_block--;
    
    while (ind_in_block>=1){
        
        compApp_mul( temp, powers + ind_in_block, coeffsTable[degree - index_of_column] + i-1, prec);
        compApp_add( column+i, column+i, temp, prec);
        
        ind_in_block--;
        i++;
    }
    
    compApp_clear(temp);
    
//     printf( "column %ld: ", index_of_column );
//     for ( slong i = 0; i < index_of_column+1; i++){
//         compApp_printd(column+i, 10);
//         printf(" || ");
//     }
//     printf("\n");
}

void compute_column_index0( compApp_ptr column, 
                     compApp_ptr powers, 
                     compApp ** coeffsTable,
                     const compApp_poly_t p,
                     slong ind_of_block,
                     slong len_of_block,
                     slong degree,
                     slong prec ) {
    
    slong index_of_column = ind_of_block*len_of_block;
//     printf("index_of_column: %ld\n", index_of_column);
    
    slong i = 0;
    slong sym = index_of_column - i;
    
    /* first element in column: (p->coeffs + index_of_column) */
    compApp_add( column+i, column+i, p->coeffs + index_of_column, prec);
    /* last element in column: (p->coeffs + index_of_column) */
    if (index_of_column>0) /*not when i==sym==0*/
        compApp_add( column+sym, column+sym, p->coeffs + index_of_column, prec);
    
    i = 1;
    sym = index_of_column - i;
    
    slong stop = (slong) floor(index_of_column/2);
    
    while (i <= stop) {    
        
        compApp_add( column+i, column+i, coeffsTable[degree - index_of_column] + i-1, prec);
        
        sym = index_of_column - i;
        if ( (sym <= (index_of_column)) && (i!=sym) )
            compApp_add( column+sym, column+sym, coeffsTable[degree - index_of_column] + i-1, prec);
        i++;
    }
    
//     printf( "column %ld: ", index_of_column );
//     for ( slong i = 0; i < index_of_column+1; i++){
//         compApp_printd(column+i, 10);
//         printf(" || ");
//     }
//     printf("\n");
}

void add_block2( compApp_poly_t pres,
                compApp_ptr column, 
                compApp_ptr powers, 
                compApp ** coeffsTable,
                slong ind_of_block,
                slong len_of_block,
                slong degree,
                slong prec ) {
    
    /* one block is computed; ind_in_block is 0 */
    
    slong index_of_column = ind_of_block*len_of_block;
//     printf("index_of_column: %ld\n", index_of_column);
    for ( slong i = 0; i < index_of_column; i++){
       
        if (ind_of_block>0){
            /* multiply by power index_of_column-i*/
            compApp_mul( column + i, column + i, powers + (index_of_column-i), prec);
        }
        /* index of coeff in pres: i */
        compApp_add(pres->coeffs + i, pres->coeffs + i, column + i, prec);
        /* set column +i to zero */
        compApp_zero( column+i );
    }
    for ( slong i = 0; (i < len_of_block) && (index_of_column + i < degree +1); i++){
        
        compApp_add(pres->coeffs + i + index_of_column, pres->coeffs + i + index_of_column, column + i + index_of_column, prec);
        compApp_zero( column+i + index_of_column);
        
    }
    
}

void my_Taylor_Shift2( compApp_poly_t pres, 
                      const compApp_poly_t p, 
                      const compRat_t center,
                      compApp ** coeffsTable,
                      slong degree, 
                      slong prec){
    
    compApp_poly_fit_length(pres, degree+1);
    
    slong len_of_block = (slong) floor(sqrt( (double) degree+1 ));
//     slong len_of_block = 128;
    slong num_of_block = (slong) (degree+1)/len_of_block;
    slong ind_in_block = (degree+1) - num_of_block*len_of_block;
    slong ind_of_block = num_of_block;
    slong r = degree+1-len_of_block +1;
    if (ind_in_block==0){
        ind_in_block = len_of_block-1;
        ind_of_block = ind_of_block-1;
    }
    else {
        r = num_of_block*len_of_block +1;
        ind_in_block = ind_in_block-1;
    }

//     printf("len_of_block: %ld, num_of_block: %ld, ind_in_block: %ld, r: %ld\n", len_of_block, num_of_block, ind_in_block, r);
    
    /*compute powers of the point*/
    compApp_ptr powers = (compApp_ptr) ccluster_malloc  (r*sizeof(compApp));
    for (slong i = 0; i<r; i++){
        compApp_init( powers+i );
    }
    compute_powers( powers, center, r, prec);
   
    /*column temp*/
    compApp_ptr column_temp = (compApp_ptr) ccluster_malloc  ((degree+1)*sizeof(compApp));
    for (slong i = 0; i<degree+1; i++){
        compApp_init( column_temp+i );
        compApp_zero( column_temp+i );
    }
    
    /* begin with highest degree coefficient */
    while ( (ind_of_block>=0)&&(ind_in_block>=0) ) {
//         printf("index_in_block: %ld, index_of_block: %ld\n", ind_in_block,ind_of_block);
        compute_column2( column_temp, powers, coeffsTable, p, ind_in_block, ind_of_block, len_of_block, degree, prec );
        
        if (ind_in_block==0){
            add_block2( pres, column_temp, powers, coeffsTable, ind_of_block, len_of_block, degree, prec );
            ind_in_block = len_of_block -1;
            ind_of_block--;
        }
        else
            ind_in_block--;
    }
    
    /* clear powers */
    for (slong i = 0; i<r; i++)
        compApp_clear(powers+i);
    ccluster_free(powers);
    
    /* clear column_temp */
    for (slong i = 0; i<degree+1; i++)
        compApp_clear(column_temp+i);
    ccluster_free(column_temp);
    
    pres->length=degree+1;
//     printf("pres: \n"); compApp_poly_printd(pres, 10); printf("\n\n");
    
}

int main() {
    
//     slong degree = 10;
    slong degree = 9;
    slong prec = 53;
    int nbtests = 1000;
    
    compApp_poly_t p, pshift1, pshift2, pdiff;
    realRat_poly_t pbern;
    compRat_t center;
    compApp_t c;
    compRat_init(center);
    compApp_init(c);
    compApp_poly_init(p);
    compApp_poly_init(pshift1);
    compApp_poly_init(pshift2);
    compApp_poly_init(pdiff);
    realRat_poly_init(pbern);
    
    compRat_set_sisi(center, 1,1,1,1);
//     compRat_set_sisi(center, 2,1,0,1);
    compApp_set_compRat(c, center, prec);
//     bernoulli_polynomial(pbern, degree);
    realRat_poly_ones(pbern, degree);
    
//     mignotte_polynomial(pbern, degree, 200);
//     realRat_poly_wilkinson(pbern, degree);
//     realRat_poly_wilkRat(pbern, degree);
//     realRat_poly_wilkFac(pbern, degree);
//     realRat_poly_wilkRat(pbern, degree);
//     regularGrid_polynomial(pbern, 7);
    degree = (pbern->length) -1;
    
    compApp_poly_set_realRat_poly( p, pbern, prec);
//     printf("p: \n"); compApp_poly_printd(p, 10); printf("\n\n");
    
    compApp_poly_ptr ptab = (compApp_poly_ptr) malloc ((degree+1)*sizeof(compApp_poly));
    
    
    fmpz ** binCoeffsTable = NULL;
    binCoeffsTable = precompute_binomial_table( binCoeffsTable, degree );
    
    compApp ** coeffsTable = NULL;
    coeffsTable = precompute_coeffsTable( coeffsTable, binCoeffsTable, p, degree, prec );
    
//     my_Taylor_Shift( pshift2, p, center, binCoeffsTable, degree, prec);
    
    my_Taylor_Shift2( pshift2, p, center, coeffsTable, degree, prec);
    
//     clear_binomial_table( binCoeffsTable, degree );
    
    compApp_poly_set(pshift1, p);
    _acb_poly_taylor_shift_convolution(pshift1->coeffs, c, pshift1->length, prec);
//     printf("pshift1: \n"); compApp_poly_printd(pshift1, 10); printf("\n\n");
    
    compApp_poly_sub(pdiff, pshift1, pshift2, prec);
//     printf("pdiff: \n"); compApp_poly_printd(pdiff, 10); printf("\n\n");
    
    /* compute the norm 2 of the difference */
    compApp_t norm2, temp;
    compApp_init(norm2);
    compApp_init(temp);
    compApp_zero(norm2);
    for (slong i=0; i<=degree; i++){
        compApp_mul(temp, (pdiff->coeffs)+i,  (pdiff->coeffs)+i, prec);
        compApp_add(norm2, norm2,  temp, prec);
    }
    compApp_sqrt(norm2, norm2, prec);
        
    printf("norm 2 of the dirrerence: \n"); compApp_printd(norm2, 10); printf("\n\n");
    compApp_clear(norm2);
    compApp_clear(temp);
    
    clock_t ti;
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_set(pshift1, p);
//         _acb_poly_taylor_shift_horner(pshift1->coeffs, c, pshift1->length, prec);
        _acb_poly_taylor_shift_convolution(pshift1->coeffs, c, pshift1->length, prec);
//         _acb_poly_taylor_shift_divconquer(pshift1->coeffs, c, pshift1->length, prec);
//         printf("pshift1: \n"); compApp_poly_printd(pshift1, prec); printf("\n\n");
    }
    ti = clock() - ti;
    printf ("time for %d convo taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++) {
        compApp_poly_set(pshift1, p);
//         _acb_poly_taylor_shift_horner(pshift1->coeffs, c, pshift1->length, prec);
//         _acb_poly_taylor_shift_convolution(pshift1->coeffs, c, pshift1->length, prec);
        _acb_poly_taylor_shift_divconquer(pshift1->coeffs, c, pshift1->length, prec);
//         printf("pshift1: \n"); compApp_poly_printd(pshift1, prec); printf("\n\n");
    }
    ti = clock() - ti;
    printf ("time for %d div conq taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    precompute( ptab, p, prec );
    ti = clock() - ti;
    printf ("time for precomputing derivatives: %f seconds.\n", ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++)
        compApp_poly_taylorShift_home( pshift2, ptab, c, prec);
    ti = clock() - ti;
    printf ("time for %d home taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
//     printf("pshift2: \n"); compApp_poly_printd(pshift2, prec); printf("\n\n");
    
    ti = clock();
    binCoeffsTable = precompute_binomial_table( binCoeffsTable, degree );
    ti = clock() - ti;
    printf ("time for precomputing binomial coeffs: %f seconds.\n", ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++)
        my_Taylor_Shift( pshift2, p, center, binCoeffsTable, degree, prec);
    ti = clock() - ti;
    printf ("time for %d my taylor shifts, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    coeffsTable = precompute_coeffsTable( coeffsTable, binCoeffsTable, p, degree, prec );
    ti = clock() - ti;
    printf ("time for precomputing coeffs: %f seconds.\n", ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    slong len_of_block = (slong) floor(sqrt( (double) degree+1 ));
    slong r = degree+1-len_of_block +1;
    for (int i = 0; i<nbtests; i++){
        /*compute powers of the point*/
        compApp_ptr powers = (compApp_ptr) ccluster_malloc  (r*sizeof(compApp));
        for (slong i = 0; i<r; i++){
            compApp_init( powers+i );
        }
        compute_powers( powers, center, r, prec);
    
        /* clear powers */
        for (slong i = 0; i<r; i++)
            compApp_clear(powers+i);
        ccluster_free(powers);
    }
    ti = clock() - ti;
    printf ("time for %d precomputing %ld powers: %f seconds.\n", nbtests, r, ((float)ti)/CLOCKS_PER_SEC);
    
    ti = clock();
    for (int i = 0; i<nbtests; i++)
        my_Taylor_Shift2( pshift2, p, center, coeffsTable, degree, prec);
    ti = clock() - ti;
    printf ("time for %d my taylor shift 2, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
//     ti = clock();
//     for (int i = 0; i<nbtests; i++) {
//         compApp_poly_evaluate( (pshift2->coeffs), ptab, c, prec);
//         compApp_poly_evaluate( (pshift2->coeffs)+1, ptab+1, c, prec);
//     }
//     ti = clock() - ti;
//     printf ("time for %d rectangular evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
//     
//     ti = clock();
//     for (int i = 0; i<nbtests; i++) {
//         compApp_poly_evaluate( (pshift2->coeffs), ptab+degree/2, c, prec);
//         compApp_poly_evaluate( (pshift2->coeffs)+1, ptab+degree/2 +1, c, prec);
//     }
//     ti = clock() - ti;
//     printf ("time for %d rectangular evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
//     
//     ti = clock();
//     for (int i = 0; i<nbtests; i++) {
//         compApp_poly_evaluate_horner( (pshift2->coeffs), ptab, c, prec);
//         compApp_poly_evaluate_horner( (pshift2->coeffs)+1, ptab+1, c, prec);
//     }
//     ti = clock() - ti;
//     printf ("time for %d horner evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
//     
//     ti = clock();
//     for (int i = 0; i<nbtests; i++) {
//         compApp_poly_evaluate_horner( (pshift2->coeffs), ptab+degree/2, c, prec);
//         compApp_poly_evaluate_horner( (pshift2->coeffs)+1, ptab+degree/2 +1, c, prec);
//     }
//     ti = clock() - ti;
//     printf ("time for %d horner evaluations, degree %d, prec %d: %f seconds.\n", nbtests, (int) degree,(int) prec, ((float)ti)/CLOCKS_PER_SEC);
    
    clear_binomial_table( binCoeffsTable, degree );
    clear_coeffsTable( coeffsTable, degree );
    
    clear_table( ptab );
    
    free(ptab);
    compRat_clear(center);
    compApp_clear(c);
    compApp_poly_clear(p);
    compApp_poly_clear(pshift1);
    compApp_poly_clear(pshift2);
    realRat_poly_clear(pbern);

    
    return 0;
}

/* new home */
fmpz ** precompute_binomial_table( fmpz ** binCoeffsTable, slong degree ){
    
    binCoeffsTable = (fmpz **) ccluster_malloc ( (degree)*sizeof(fmpz *) );
    for (slong i = 0; i < degree; i++)
        binCoeffsTable[i] = _fmpz_vec_init(degree-i);
    /*fill column by column*/
    for (slong j = 0; j < degree; j++)
        for (slong i = degree - j -1; i >=0; i--){
            if (j==0)
                fmpz_set_si(binCoeffsTable[i]+j, degree -i);
            else {
                if (i == (degree-j-1) )
                    fmpz_set_si(binCoeffsTable[i]+j, 1);
                else
                    fmpz_add( binCoeffsTable[i]+j, binCoeffsTable[i+1]+(j), binCoeffsTable[i+1]+(j-1) );
            }
                
        }
    
//     for (slong i = 0; i < degree; i++){
//         for (slong j=0; j < degree-i; j++){
//             printf("%4ld ", fmpz_get_si(binCoeffsTable[i]+j));
//         }
//         printf("\n");
//     }
    
    return binCoeffsTable;
}

compApp ** precompute_coeffsTable( compApp ** coeffsTable, fmpz ** binCoeffsTable, const compApp_poly_t p, slong degree, slong prec ){
    
    coeffsTable = (compApp **) ccluster_malloc ( (degree)*sizeof(compApp *) );
    for (slong i = 0; i < degree; i++){
        coeffsTable[i] = (compApp *) ccluster_malloc ( (degree-i)*sizeof(compApp) );
        for (slong j=0; j < degree-i; j++)
            compApp_init(coeffsTable[i] +j);
    }
    /*fill line by line*/
    for (slong i = 0; i < degree; i++)
        for (slong j = 0; j <degree -i; j++)
            compApp_mul_fmpz( coeffsTable[i] +j, (p->coeffs) + degree -i, binCoeffsTable[i]+j, prec ); 
    
//     for (slong i = 0; i < degree; i++){
//         for (slong j=0; j < degree-i; j++){
//             compApp_printd( coeffsTable[i] +j, 10 ); printf("  "); 
//         }
//         printf("\n");
//     }
    
    return coeffsTable;
}
void clear_coeffsTable( compApp ** coeffsTable, slong degree ){
    for (slong i = 0; i < degree; i++){
        for (slong j = 0; j < degree -i; j++)
            compApp_clear(coeffsTable[i] + j);
        ccluster_free(coeffsTable[i]);
    }
    ccluster_free(coeffsTable);
}

void clear_binomial_table( fmpz ** binCoeffsTable, slong degree ){
    for (slong i = 0; i < degree; i++)
        _fmpz_vec_clear(binCoeffsTable[i], degree-i);
    ccluster_free(binCoeffsTable);
}


/* old home */
void compApp_poly_div_si ( compApp_poly_t res, const compApp_poly_t p, slong y, slong prec) {
    slong len = p->length;
    while (len>0) {
        compApp_div_si( (res->coeffs) + (len-1), (p->coeffs) + (len-1), y, prec);
        len--;
    }
}

void precompute( compApp_poly_ptr tab, const compApp_poly_t p, slong prec ){
    
    slong d = p->length -1;
    slong i;
    compApp_poly_init2( tab, p->length);
    compApp_poly_set( tab, p);
//     printf("0-th derivative of p: \n"); compApp_poly_printd(tab, prec); printf("\n\n");
    for (i=1; i<=d; i++) {
        compApp_poly_init2( tab+i, p->length -i);
        compApp_poly_derivative( tab + i, tab +(i-1), prec);
        compApp_poly_div_si ( tab + i, tab + i, i, prec);
//         printf("%i-th derivative of p: \n", (int) i); compApp_poly_printd(tab + i, prec); printf("\n\n");
    }
    
}

void compApp_poly_taylorShift_home( compApp_poly_t res, compApp_poly_ptr tab, const compApp_t center, slong prec){
    slong lentable = tab->length;
//     printf("len of table: %d\n", (int) lentable);
    compApp_poly_fit_length(res, lentable);
    slong i;
    for (i=0;i<lentable;i++){
//         printf("i: %d\n", (int) i);
//         printf("i-th coeff: "); compApp_printd( (res->coeffs) + i, prec); printf("\n");
        compApp_poly_evaluate( (res->coeffs) + i, tab+i, center, prec);
//         if (i+1<lentable) {
//             compApp_poly_evaluate2( (res->coeffs) + i,(res->coeffs) + i+1, tab+i, center, prec);
//             compApp_div_si( (res->coeffs) + i+1, (res->coeffs) + i+1, i+1, prec);
//             i++;
//         }
//         else {
//             compApp_poly_evaluate( (res->coeffs) + i, tab+i, center, prec);
//         }
    }
    compApp_poly_set_length(res, lentable);
}

void clear_table( compApp_poly_ptr tab ) {
    slong lentable = tab->length;
    slong i;
    for (i=0; i<lentable; i++) compApp_poly_clear(tab + i);
}
/* pols */

void regularGrid_polynomial(realRat_poly_t preg, int res){
    
    realRat_poly_t ptemp, ptemp2;
    realRat_poly_init2(ptemp,2);
    realRat_poly_init2(ptemp2,2);
    realRat_poly_one(preg);
    realRat_poly_zero(ptemp);
    realRat_poly_zero(ptemp2);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_poly_set_coeff_si_ui(ptemp2, 2, 1, 1);
    
    for (int i=0; i<=res; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(preg, preg, ptemp);
        
        for (int j=1; j<=res; j++){
            realRat_poly_set_coeff_si_ui(ptemp2, 1, 2*i, 1);
            realRat_poly_set_coeff_si_ui(ptemp2, 0, i*i+j*j, 1);
            realRat_poly_mul(preg, preg, ptemp2);
        }
        
        if (i>0) {
            realRat_poly_set_coeff_si_ui(ptemp, 0, i, 1);
            realRat_poly_mul(preg, preg, ptemp);
            
            for (int j=1; j<=res; j++){
                realRat_poly_set_coeff_si_ui(ptemp2, 1, -2*i, 1);
                realRat_poly_set_coeff_si_ui(ptemp2, 0, i*i+j*j, 1);
                realRat_poly_mul(preg, preg, ptemp2);
            }
        }
            
    }
    
    realRat_poly_clear(ptemp);
    realRat_poly_clear(ptemp2);
}

void realRat_poly_ones( realRat_poly_t res, slong degree ){
    for (int i=0; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(res, i, 1, 1);
    }
}

void realRat_poly_wilkinson( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(res, res, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkRat( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, ((ulong) degree)+1);
        realRat_poly_mul(res, res, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkFac( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_t faci, mfaci;
    realRat_init(faci);
    realRat_init(mfaci);
    realRat_set_si(faci, 1,1);
    
    for (int i=1; i<=degree; i++){
        realRat_mul_si(mfaci, faci, -1);
        realRat_poly_set_coeff_realRat(ptemp, 0, mfaci);
        realRat_poly_mul(res, res, ptemp);
        realRat_mul_si(faci, faci, (slong) i+1 );
    }
    
    realRat_clear(faci);
    realRat_clear(mfaci);
    realRat_poly_clear(ptemp);
}

void realRat_poly_wilkFacRat( realRat_poly_t res, slong degree ){
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(res);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_t faci, mfaci, facd;
    realRat_init(faci);
    realRat_init(mfaci);
    realRat_init(facd);
    realRat_set_si(faci, 1,1);
    realRat_set_si(facd, 1,1);
    
    for (int i=1; i<degree; i++)
        realRat_mul_si(facd, facd, (slong) i+1 );
    
    for (int i=1; i<=degree; i++){
        realRat_mul_si(mfaci, faci, -1);
        realRat_div(mfaci, mfaci, facd);
        realRat_poly_set_coeff_realRat(ptemp, 0, mfaci);
        realRat_poly_mul(res, res, ptemp);
        realRat_mul_si(faci, faci, (slong) i+1 );
    }
    
    realRat_clear(faci);
    realRat_clear(mfaci);
    realRat_clear(facd);
    realRat_poly_clear(ptemp);
}
