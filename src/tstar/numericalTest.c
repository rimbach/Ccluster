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

#include "tstar/tstar.h"

#ifdef CCLUSTER_EXPERIMENTAL

// int tstar_D0_test(cacheApp_t cache, const compApp_t coeff0, 
//                                     const compApp_t coeff1, 
//                                     const compApp_t coeffn, 
//                                     const compDsk_t d, int depth, slong prec, metadatas_t meta ){
//     slong degree = cacheApp_getDegree(cache);
//     slong i = 0;
//     int contains_zero = 1;
//     compApp_t c, ball, pball, pballsave, pballApp, pballApp1, pballCor, ballmc, ballmctoi;
//     realApp_t sumofabscoeff, abscoeff0, abscoeffi, abscoeffn;
//     realRat_t nrad, fac, scaleFactor;
//     compApp_poly_t pp;
//     
//     /*initialize variables*/
//     compApp_init(c);
//     compApp_init(ball);
//     compApp_init(pball);
//     compApp_init(pballsave);
//     compApp_init(pballApp);
//     compApp_init(pballApp1);
//     compApp_init(pballCor);
//     compApp_init(ballmc);
//     compApp_init(ballmctoi);
//     realRat_init(nrad);
//     realRat_init(fac);
//     compApp_poly_init2(pp, degree+1);
//     
//     realApp_init(abscoeff0);
//     realApp_init(abscoeffi);
//     realApp_init(abscoeffn);
//     realApp_init(sumofabscoeff);
//     realRat_init(scaleFactor);
//     
//     tstar_getApproximation( pp, cache, prec, meta);
//     compApp_zero(pballsave);
// //     compApp_poly_set(pp, pApprox);
//     realRat_set_si(fac,3,4);
//     realRat_mul(nrad, compDsk_radiusref(d), fac);
//     compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
//     compApp_set_compRat( c, compDsk_centerref(d), prec);
//     compApp_sub(ballmc, ball, c, prec);
//     compApp_one(ballmctoi);
//     realRat_set_si(fac,1,1);
//     
//     compApp_set(pballApp, coeff0);
//     compApp_set(pballApp1, coeff1);
//     compApp_abs(abscoeff0, pballApp, prec);
//     compApp_abs(abscoeffi, pballApp1, prec);
//     compApp_abs(abscoeffn, coeffn, prec);
//     realRat_pow_si(scaleFactor, compDsk_radiusref(d), degree);
//     realApp_mul_realRat( abscoeffn, abscoeffn, scaleFactor, prec );
//     realApp_mul_realRat( abscoeffi, abscoeffi, compDsk_radiusref(d), prec );
//     realApp_add(sumofabscoeff, abscoeffi, abscoeffn, prec);
//     realRat_set_si(scaleFactor, 1,1);
//     int dominates = 1;
//     
//     /* order 0 evaluation */
//     tstar_evaluate(pballCor, pp, ball, prec, meta, depth);
//     contains_zero = compApp_contains_zero(pballCor);
//     
// //     if (contains_zero) {
// //     realRat_set_si(scaleFactor, 1,20000);
// //     realRat_mul(scaleFactor, scaleFactor, nrad);
// //     compApp_set_compDsk( ball, compDsk_centerref(d), scaleFactor, prec);
// //     tstar_evaluate(pballsave, pp, ball, prec, meta, depth);
// //     int contains_zero2 = compApp_contains_zero(pballsave);
// //     if (contains_zero2==0) dominates = 0;
// //     compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
// //     compApp_zero(pballsave);
// //     realRat_set_si(scaleFactor, 1,1);
// //     }
//     
//     while (dominates&&(contains_zero!=0)&&( i < degree )) {
// //     while (dominates&&(contains_zero!=0)&&( i < 3 )) {
//         
//         /* actualize pballsave */
//         if (i==0) {}                                            /* p(c) is already in pballApp */
//         else if (i==1) compApp_set(pballApp, pballApp1);        /* pballApp <- pballApp1=p'(c)*/
//         else if (i < degree-1 ) { 
//             tstar_evaluate(pballApp, pp, c, prec, meta, depth); /* pballApp <- (p^(i))(c) */
//             compApp_mul_realRat_in_place(pballApp, fac, prec);  /* pballApp <- pballApp / i! */
//             compApp_abs(abscoeffi, pballApp, prec);
//             realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec ); /* abscoeffi <- |pballApp| *( r^i / i!) */
//             realApp_add(sumofabscoeff, sumofabscoeff, abscoeffi, prec);
//         }
//         else compApp_set(pballApp, coeffn);                     /* pballApp <- coeffn=(p^(degree))(c)*/
//         
//         compApp_mul(pballApp, pballApp, ballmctoi, prec); /* pballApp <- pballApp * (B-c)^i */
//         compApp_add( pballsave, pballsave, pballApp, prec); /* pballsave = sum of (p^(i))(c) * ( (B-c)^i / i! ) */
//         
//         i++; /* the actual order of the evaluation */
//         compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
//         fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
//         realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
//         tstar_getDerivative (pp, cache, prec, i , meta);
//         
//         tstar_evaluate(pballCor, pp, ball, prec, meta, depth); /* pballCor <- (p^(i))(B)*/
//         compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
//         compApp_mul(pballCor, pballCor, ballmctoi, prec);      /* pballCor <- pballCor*(B-c)^i */
//         compApp_add( pball, pballsave, pballCor, prec);
//         
//         contains_zero = compApp_contains_zero(pball);
//         
//         if (i>1) {
//             dominates = realApp_soft_compare( abscoeff0, sumofabscoeff, prec );
//             if (dominates!=1) dominates=0;
//         }
//         
//     }
//     
//     /* clear variables */
//     compApp_clear(c);
//     compApp_clear(ball);
//     compApp_clear(pball);
//     compApp_clear(pballsave);
//     compApp_clear(pballApp);
//     compApp_clear(pballApp1);
//     compApp_clear(pballCor);
//     compApp_clear(ballmc);
//     compApp_clear(ballmctoi);
//     realRat_clear(nrad);
//     realRat_clear(fac);
//     compApp_poly_clear(pp);
//     realApp_clear(abscoeff0);
//     realApp_clear(abscoeffi);
//     realApp_clear(abscoeffn);
//     realApp_clear(sumofabscoeff);
//     realRat_clear(scaleFactor);
//     
//     if ((contains_zero==0)||(dominates)){
// //     if ((contains_zero==0)){    
// //         printf("&&depth: %d, success of D0 discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
//         return 0;
//     }
//     else {
// //         printf("&&depth: %d, fail    of D0 discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
//         return -1;
//     }
// }

void completeTaylorShift(cacheApp_t cache, compApp_poly_t shiftedPol, int nbCoeffComputed, const compDsk_t d, slong prec, metadatas_t meta ){
    
    slong deg = cacheApp_getDegree(cache);
    slong deg2 = deg - (nbCoeffComputed-1);
    
    compApp_poly_set_length(shiftedPol, deg+1);
    
    if (deg2>1) {
        compApp_poly_t pp;
        realRat_t fac, scaleFactor;
        realRat_init(fac);
        realRat_init(scaleFactor);
//         printf("deg: %d, nbCoeffComputed: %d, deg2: %d\n", (int) deg, (int) nbCoeffComputed, (int) deg2);
        compApp_poly_init2(pp, deg2+1);
        
        tstar_getDerivative (pp, cache, prec, nbCoeffComputed , meta);
        tstar_taylor_shift_inplace( pp, d, prec, meta);
        
        realRat_pow_si(scaleFactor, compDsk_radiusref(d), (slong) nbCoeffComputed);
        realRat_set_si(fac, 1,1);
        for(int i=2; i<nbCoeffComputed; i++) 
            fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
        
        for(int i=nbCoeffComputed; i<=deg; i++) {
            
            fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
            if ((i-nbCoeffComputed)>=1)
                fmpz_mul_si( realRat_numref(fac), realRat_numref(fac), ((slong) (i-nbCoeffComputed)));
            
            compApp_mul_realRat( (shiftedPol->coeffs) + i,
                                (pp->coeffs) + (i-nbCoeffComputed),
                                scaleFactor, prec);
            compApp_mul_realRat_in_place( (shiftedPol->coeffs) + i,
                                        fac, prec);
        }
        
        realRat_clear(fac);
        realRat_clear(scaleFactor);
        compApp_poly_clear(pp);
    }
}
    
int tstar_D0_test(cacheApp_t cache, 
//                   const compApp_t coeff0, 
//                   const compApp_t coeff1,
                  compApp_poly_t shiftedPol,
                  int * nbCoeffComputed,
                  const compDsk_t d, int depth, slong prec, metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    realRat_t fac, nrad, scaleFactor;
    realApp_t abscoeff0, abscoeff1, abscoeffi, approx, sumofabscoeff;
    compApp_t c, ball, pballCor, ballmc, ballmctoi;
    compApp_poly_t pp;
    
    realRat_init(fac);
    realRat_init(nrad);
    realRat_init(scaleFactor);
    realApp_init(abscoeff0);
    realApp_init(abscoeff1);
    realApp_init(abscoeffi);
    realApp_init(approx);
    realApp_init(sumofabscoeff);
    compApp_init(c);
    compApp_init(ball);
    compApp_init(pballCor);
    compApp_init(ballmc);
    compApp_init(ballmctoi);
    compApp_poly_init2(pp, degree+1);
    
//     realRat_set_si(fac,3,4);
    realRat_set_si(fac,1,1);
    realRat_mul(nrad, compDsk_radiusref(d), fac);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
    compApp_set_compRat( c, compDsk_centerref(d), prec);
    compApp_sub(ballmc, ball, c, prec);
    compApp_set(ballmctoi, ballmc);
    realRat_set_si(fac,1,2);
    realRat_set(scaleFactor, compDsk_radiusref(d));
    
//     compApp_abs(abscoeff0, coeff0, prec);
//     compApp_abs(abscoeff1, coeff1, prec);
    compApp_abs(abscoeff0, (shiftedPol->coeffs)+0, prec);
    compApp_abs(abscoeff1, (shiftedPol->coeffs)+1, prec);
//     realApp_mul_realRat( abscoeff1, abscoeff1, scaleFactor, prec );
    realApp_set(sumofabscoeff, abscoeff1);
    
    realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
    compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
    
    slong i = 2;
    int dominates = 1;
    int res = 0;
    while ( (res==0)&&dominates && (i<=degree)) {   
        tstar_getDerivative (pp, cache, prec, i , meta);
        tstar_evaluate(pballCor, pp, ball, prec, meta, depth); /* pballCor <- (p^(i))(B)*/
        compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
        compApp_mul(pballCor, pballCor, ballmctoi, prec);      /* pballCor <- pballCor*(B-c)^i */
        compApp_abs(abscoeffi, pballCor, prec);
        realApp_add(approx, sumofabscoeff, abscoeffi, prec);
        
        res = realApp_soft_compare( abscoeff0, approx, prec );
        
        if (!(res==1)){
            res = 0;
//             tstar_evaluate(pballCor, pp, c, prec, meta, depth);    /* pballCor <- (p^(i))(c)*/
//             compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
//             compApp_abs(abscoeffi, pballCor, prec);
//             realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec );
            
            tstar_evaluate((shiftedPol->coeffs)+((int) i), pp, c, prec, meta, depth);
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), fac, prec);
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), scaleFactor, prec);
            compApp_abs(abscoeffi, (shiftedPol->coeffs)+((int) i), prec);
            
            realApp_add(sumofabscoeff, sumofabscoeff, abscoeffi, prec);
            dominates = realApp_soft_compare( abscoeff0, sumofabscoeff, prec );
            if (!(dominates==1)) dominates = 0;
            i++;
            compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
            realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
            fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
            *nbCoeffComputed = (int) i;
        }
        
    }
    
    realRat_clear(fac);
    realRat_clear(nrad);
    realRat_clear(scaleFactor);
    realApp_clear(abscoeff0);
    realApp_clear(abscoeff1);
    realApp_clear(abscoeffi);
    realApp_clear(approx);
    realApp_clear(sumofabscoeff);
    compApp_clear(c);
    compApp_clear(ball);
    compApp_clear(pballCor);
    compApp_clear(ballmc);
    compApp_clear(ballmctoi);
    compApp_poly_clear(pp);
    
    if (res==1) {
//         printf("&&depth: %d, success of D0 discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return 0;
    }
    else {
//         printf("&&depth: %d, fail    of D0 discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return -1;
    }
}
    
int tstar_N1_test(cacheApp_t cache, const compApp_t coeff0, 
                                    const compApp_t coeff1, 
                                    const compDsk_t d, int depth, slong prec, metadatas_t meta ){
    
    slong degree = cacheApp_getDegree(cache);
    realRat_t fac, nrad, scaleFactor;
    realApp_t abscoeff0, abscoeff1, abscoeffi, approx, sumofabscoeff;
    compApp_t c, ball, pballCor, ballmc, ballmctoi;
    compApp_poly_t pp;
    
    realRat_init(fac);
    realRat_init(nrad);
    realRat_init(scaleFactor);
    realApp_init(abscoeff0);
    realApp_init(abscoeff1);
    realApp_init(abscoeffi);
    realApp_init(approx);
    realApp_init(sumofabscoeff);
    compApp_init(c);
    compApp_init(ball);
    compApp_init(pballCor);
    compApp_init(ballmc);
    compApp_init(ballmctoi);
    compApp_poly_init2(pp, degree+1);
    
//     realRat_set_si(fac,3,4);
    realRat_set_si(fac,1,1);
    realRat_mul(nrad, compDsk_radiusref(d), fac);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
    compApp_set_compRat( c, compDsk_centerref(d), prec);
    compApp_sub(ballmc, ball, c, prec);
    compApp_set(ballmctoi, ballmc);
    realRat_set_si(fac,1,2);
    realRat_set(scaleFactor, compDsk_radiusref(d));
    
    compApp_abs(abscoeff0, coeff0, prec);
    compApp_abs(abscoeff1, coeff1, prec);
    realApp_mul_realRat( abscoeff1, abscoeff1, scaleFactor, prec );
    realApp_set(sumofabscoeff, abscoeff0);
    
    realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
    compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
    
    slong i = 2;
    int dominates = 1;
    int res = 0;
    while ( (res==0)&&dominates && (i<=degree)) {   
        tstar_getDerivative (pp, cache, prec, i , meta);
        tstar_evaluate(pballCor, pp, ball, prec, meta, depth); /* pballCor <- (p^(i))(B)*/
        compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
        compApp_mul(pballCor, pballCor, ballmctoi, prec);      /* pballCor <- pballCor*(B-c)^i */
        compApp_abs(abscoeffi, pballCor, prec);
        realApp_add(approx, sumofabscoeff, abscoeffi, prec);
        
        res = realApp_soft_compare( abscoeff1, approx, prec );
        
        if (!(res==1)){
            res = 0;
            tstar_evaluate(pballCor, pp, c, prec, meta, depth);    /* pballCor <- (p^(i))(c)*/
            compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
            compApp_abs(abscoeffi, pballCor, prec);
            realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec );
            realApp_add(sumofabscoeff, sumofabscoeff, abscoeffi, prec);
            dominates = realApp_soft_compare( abscoeff1, sumofabscoeff, prec );
            if (!(dominates==1)) dominates = 0;
            i++;
            compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
            realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
            fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
        }
        
    }
    
    realRat_clear(fac);
    realRat_clear(nrad);
    realRat_clear(scaleFactor);
    realApp_clear(abscoeff0);
    realApp_clear(abscoeff1);
    realApp_clear(abscoeffi);
    realApp_clear(approx);
    realApp_clear(sumofabscoeff);
    compApp_clear(c);
    compApp_clear(ball);
    compApp_clear(pballCor);
    compApp_clear(ballmc);
    compApp_clear(ballmctoi);
    compApp_poly_clear(pp);
    
    if (res==1) {
//         printf("&&depth: %d, success of N1 non-discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return 1;
    }
    else return -1;
}

int tstar_C0_test( cacheApp_t cache, const compDsk_t d, int depth, slong prec, metadatas_t meta ) {

//     slong order = pApprox->length;
    slong i = 0;
    int contains_zero = 1;
    compApp_t c, ball, pball, pballsave, pballApp, pballApp1, pballCor, ballmc, ballmctoi;
    realApp_t abscoeff0, abscoeffi, sumofabscoeff;
    realRat_t nrad, fac, scaleFactor;
    compApp_poly_t pApprox, pp;
    
    /*initialize variables*/
    compApp_init(c);
    compApp_init(ball);
    compApp_init(pball);
    compApp_init(pballsave);
    compApp_init(pballApp);
    compApp_init(pballApp1);
    compApp_init(pballCor);
    compApp_init(ballmc);
    compApp_init(ballmctoi);
    realRat_init(nrad);
    realRat_init(fac);
    compApp_poly_init2(pApprox, cacheApp_getDegree(cache)+1);
    compApp_poly_init2(pp, cacheApp_getDegree(cache)+1);
    
    realApp_init(abscoeff0);
    realApp_init(abscoeffi);
    realApp_init(sumofabscoeff);
    realRat_init(scaleFactor);
    
    tstar_getApproximation( pApprox, cache, prec, meta);
    compApp_zero(pballsave);
    compApp_poly_set(pp, pApprox);
    realRat_set_si(fac,3,4);
    realRat_mul(nrad, compDsk_radiusref(d), fac);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
    compApp_set_compRat( c, compDsk_centerref(d), prec);
    compApp_sub(ballmc, ball, c, prec);
    compApp_one(ballmctoi);
    realRat_set_si(fac,1,1);
    
    realRat_set_si(scaleFactor, 1,1);
    
    /* initial test */
    tstar_evaluate(pballApp, pp, c, prec, meta, depth);
    tstar_getDerivative(pp, cache, prec, 1 , meta);
    tstar_evaluate(pballApp1, pp, c, prec, meta, depth);
    compApp_poly_set(pp, pApprox);
    
    compApp_abs(abscoeff0, pballApp, prec);
    compApp_abs(abscoeffi, pballApp1, prec);
    realApp_mul_realRat( abscoeffi, abscoeffi, compDsk_radiusref(d), prec );
    int dominates = realApp_soft_compare( abscoeff0, abscoeffi, prec );
    if (dominates!=1) dominates=0;
    
    if (dominates) {
        /* order 0 evaluation */
        tstar_evaluate(pballCor, pp, ball, prec, meta, depth);
//         compApp_add( pball, pballsave, pballCor, prec);
        contains_zero = compApp_contains_zero(pballCor);
        
    //     compApp_abs(abscoeff0, pballCor, prec);
        realApp_zero(sumofabscoeff);
        compApp_zero(pballsave);
        
        while (dominates&&(contains_zero!=0)&&( i < cacheApp_getDegree(cache) )) {
//         while (dominates&&(contains_zero!=0)&&( i < 3 )) {
            
            /* actualize pballsave */
            if (i==0) {
            }
            else if (i==1) {
                compApp_set(pballApp, pballApp1);
            }
            else
                tstar_evaluate(pballApp, pp, c, prec, meta, depth); /* already done for i==0 */
            compApp_mul_realRat_in_place(pballApp, fac, prec);
//             if (i==0){
//                 compApp_abs(abscoeff0, pballApp, prec);
//                 realApp_mul_realRat( abscoeff0, abscoeff0, scaleFactor, prec );
//     //             printf("--i: %d, abscoeff: ", (int) i); realApp_printd( abscoeff0, 20); printf("\n");
//             }
//             else {
//                 compApp_abs(abscoeffi, pballApp, prec);
//                 realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec );
//     //             printf("--i: %d, abscoeff: ", (int) i); realApp_printd( abscoeffi, 20); printf("\n");
//             }
            if (i>1) {
                compApp_abs(abscoeffi, pballApp, prec);
                realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec );
            }
            compApp_mul(pballApp, pballApp, ballmctoi, prec);
            compApp_add( pballsave, pballsave, pballApp, prec);
            
            i++; /* the actual order of the evaluation */
            compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
            fmpz_mul_si( realRat_denref(fac), realRat_denref(fac), (slong) i);
            realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
    //         compApp_poly_derivative(pp, pp, prec);
    //         printf("i: %d, pp : ", i); compApp_poly_printd(pp,20); printf("\n");
    //         printf("Ici: \n");
    //         printf("i: %d pp2: ",(int) i); compApp_poly_printd(cacheApp_getDerivative ( cache, prec, i ),20); printf("\n");
    //         printf("La: \n");
//             compApp_poly_set(pp,  );
            tstar_getDerivative (pp, cache, prec, i , meta);
            
    //         printf("Et La: \n");
            tstar_evaluate(pballCor, pp, ball, prec, meta, depth);
            compApp_mul_realRat_in_place(pballCor, fac, prec);
            compApp_mul(pballCor, pballCor, ballmctoi, prec);
            compApp_add( pball, pballsave, pballCor, prec);
            
            contains_zero = compApp_contains_zero(pball);
            
            if (i>1) {
                dominates = realApp_soft_compare( abscoeff0, abscoeffi, prec );
                if (dominates!=1) dominates=0;
                if (dominates==1) {
                    realApp_add(sumofabscoeff, sumofabscoeff, abscoeffi, prec);
                    dominates = realApp_soft_compare( abscoeff0, sumofabscoeff, prec );
                    if (dominates!=1) dominates=0;
                }
            }
            
        }
    }
    
    /* clear variables */
    compApp_clear(c);
    compApp_clear(ball);
    compApp_clear(pball);
    compApp_clear(pballsave);
    compApp_clear(pballApp);
    compApp_clear(pballApp1);
    compApp_clear(pballCor);
    compApp_clear(ballmc);
    compApp_clear(ballmctoi);
    realRat_clear(nrad);
    realRat_clear(fac);
    compApp_poly_clear(pp);
    compApp_poly_clear(pApprox);
    
    realApp_clear(abscoeff0);
    realApp_clear(abscoeffi);
    realApp_clear(sumofabscoeff);
    realRat_clear(scaleFactor);
    
    
    if ((contains_zero==0)||(dominates)){
//     if ((contains_zero==0)){    
//         printf("&&depth: %d, success of C0 predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return 0;
    }
    else {
//         printf("&&depth: %d, fail    of C0 predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return -1;
    }
}

int tstar_numerical_test( compApp_poly_t pApprox, const compDsk_t d, slong prec, metadatas_t meta){
    
    int nbsols = -1;
    
    compApp_t ball, pball, dmc, ppball, pppball;
    compApp_t center, pcenter, ppcenter;
    compApp_poly_t pp;
    realRat_t half, nrad;
    compApp_init(ball);
    compApp_init(pball);
    compApp_init(dmc);
    compApp_init(ppball);
    compApp_init(pppball);
    compApp_init(center);
    compApp_init(pcenter);
    compApp_init(ppcenter);
    compApp_poly_init2(pp, pApprox->length);
    realRat_init(half);
    realRat_init(nrad);
    
    /* order 0 evaluation */
    realRat_set_si(half,3,4);
    realRat_mul(nrad, compDsk_radiusref(d), half);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
//     compApp_set_compDsk( ball, compDsk_centerref(d), compDsk_radiusref(d), prec);
    compApp_poly_evaluate(pball, pApprox, ball, prec);
    if (compApp_contains_zero(pball)==0) {
        nbsols = 0;
    }
//     printf("evaluation at order 0: "); compApp_printd(pball,10); printf("\n");
    
    /* order 2 evaluation */
    if (nbsols==-1) {
        compApp_set_compRat(center, compDsk_centerref(d), prec);
        compApp_sub(dmc, ball, center, prec);
        compApp_poly_evaluate2(pcenter,ppcenter, pApprox, center, prec);
        compApp_poly_derivative(pp,pApprox,prec);
        compApp_poly_derivative(pp,pp,prec);
        compApp_poly_evaluate(pppball, pp, ball, prec);
        realRat_set_si(half,1,2);
        compApp_mul(pball, dmc, dmc, prec);        /*val =                                    (ball-c)^2*/
        compApp_mul_realRat_in_place(pball, half, prec); /*val =                         (1/2)*(ball-c)^2*/
        compApp_mul(pball, pball, pppball, prec); /*val =                         (1/2)*(ball-c)^2*p''(ball)*/
        compApp_mul(ppcenter, ppcenter, dmc, prec); /*valt=        p'(c)*(ball-c) */
        compApp_add(pball, pball,ppcenter, prec); /*val =        p'(c)*(ball-c) + (1/2)*(ball-c)^2*p''(ball)*/
        compApp_add(pball, pball, pcenter, prec); /*val = p(c) + p'(c)*(ball-c) + (1/2)*(ball-c)^2*p''(ball)*/
        if (compApp_contains_zero(pball)==0) {
            nbsols = 0;
        }
    }
    
    /* order 3 evaluation */
    if (nbsols==-1) {
        compApp_t pppcenter;
        compApp_init(pppcenter);
        
        compApp_set_compRat(center, compDsk_centerref(d), prec);
        compApp_sub(dmc, ball, center, prec);
        
        compApp_poly_evaluate2(pcenter,ppcenter, pApprox, center, prec);
        compApp_poly_derivative(pp,pApprox,prec);
        compApp_poly_derivative(pp,pp,prec);
        compApp_poly_evaluate(pppcenter, pp, center, prec);
        compApp_poly_derivative(pp,pp,prec);
        compApp_poly_evaluate(pppball, pp, ball, prec);
        
        compApp_mul(ppcenter, ppcenter, dmc, prec); /* ppcenter = (ball-c)*p'(c) */
        
        compApp_mul(pball, dmc, dmc, prec);         /* pball    = (ball-c)^2*/
        realRat_set_si(half,1,2);
        compApp_mul_realRat_in_place(pppcenter, half, prec); /* pppcenter = (1/2)*p''(c) */
        compApp_mul(pppcenter, pppcenter, pball, prec);      /* pppcenter = (1/2)*(ball-c)^2*p''(c) */
        
        compApp_mul(pball, pball, dmc, prec);       /* pball = (ball-c)^3 */
        realRat_set_si(half,1,6);
        compApp_mul_realRat_in_place(pppball, half, prec);   /* pppball = (1/6)*p'''(ball) */
        compApp_mul(pppball, pppball, pball, prec);          /* pppball = (1/6)*(ball-c)^3*p'''(ball) */
        
        compApp_add(pball, pcenter, ppcenter, prec);         /* pball = p(c) + (ball-c)*p'(c)*/
        compApp_add(pball, pball,   pppcenter, prec);        /* pball = p(c) + (ball-c)*p'(c) + (1/2)*(ball-c)^2*p''(c)*/
        compApp_add(pball, pball,   pppball,  prec); /*pball = p(c) + (ball-c)*p'(c) + (1/2)*(ball-c)^2*p''(c) + (1/6)*(ball-c)^3*p'''(ball)*/
        if (compApp_contains_zero(pball)==0) {
            nbsols = 0;
        }
        compApp_clear(pppcenter);
    }
    
    /* try to show that there is 1 root in the box */
    if ((nbsols==-1)&&CCLUSTER_EXP_NUM_T1(meta)) {
        realApp_t coeff0, coeff1, othercoeffs;
        realApp_init(coeff0);
        realApp_init(coeff1);
        realApp_init(othercoeffs);
        
        compApp_abs( coeff0, pcenter, prec);
        compApp_abs( coeff1, ppcenter, prec);
        compApp_mul(pball, dmc, dmc, prec);        /*val =                                    (ball-c)^2*/
        compApp_mul_realRat_in_place(pball, half, prec); /*val =                         (1/2)*(ball-c)^2*/
        compApp_mul(pball, pball, pppball, prec); /*val =                         (1/2)*(ball-c)^2*p''(ball)*/
        compApp_abs( othercoeffs, pball, prec);
        realApp_add( othercoeffs, coeff0, othercoeffs, prec);
        nbsols = realApp_soft_compare( coeff1, othercoeffs, prec );
        if (nbsols==-2) {
//             printf("NOT ENOUGH PRECISION IN NUMERICAL T1TEST\n");
        }
        if (nbsols==1) nbsols = 1;
        else nbsols = -1;
        
        realApp_clear(coeff0);
        realApp_clear(coeff1);
        realApp_clear(othercoeffs);
    }
    
    compApp_clear(ball);
    compApp_clear(pball);
    compApp_clear(dmc);
    compApp_clear(ppball);
    compApp_clear(pppball);
    compApp_clear(center);
    compApp_clear(pcenter);
    compApp_clear(ppcenter);
    compApp_poly_clear(pp);
    realRat_clear(half);
    realRat_clear(nrad);
    
    return nbsols;
}
#endif