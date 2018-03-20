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

int D0N1_test ( cacheApp_t cache, compApp_poly_t pApprox, const compDsk_t d, int depth, slong prec, metadatas_t meta ){
    
    realRat_t radtothen;
    realApp_t sum, abscoeff0, abscoeff1, abscoeffn;
    compApp_t c;
    compApp_poly_t pshifted;
    slong deg = cacheApp_getDegree(cache);
    int nbCoeffcomputed=2;
    
    realRat_init(radtothen);
    realApp_init(sum);
    compApp_init(c);
    realApp_init(abscoeff0);
    realApp_init(abscoeff1);
    realApp_init(abscoeffn);
    compApp_poly_init2(pshifted, deg+1);
    
    realRat_pow_si(radtothen, compDsk_radiusref(d), deg );
    compApp_set_compRat( c, compDsk_centerref(d), prec);
    tstar_getApproximation( pApprox, cache, prec, meta);
    compApp_mul_realRat((pshifted->coeffs)+((int) deg), (pApprox->coeffs)+((int) deg), radtothen, prec);
    tstar_evaluate((pshifted->coeffs)+0, pApprox, c, prec, meta, depth);
    tstar_getDerivative(pApprox, cache, prec, 1 , meta);
    tstar_evaluate((pshifted->coeffs)+1, pApprox, c, prec, meta, depth);
    compApp_mul_realRat_in_place((pshifted->coeffs)+1, compDsk_radiusref(d), prec);
        
    compApp_abs( abscoeffn, (pshifted->coeffs)+((int) deg), prec);
    compApp_abs( abscoeff0, (pshifted->coeffs)+0, prec);
    compApp_abs( abscoeff1, (pshifted->coeffs)+1, prec);
    realApp_add(sum, abscoeff1, abscoeffn, prec);
    int test = realApp_soft_compare( abscoeff0, sum, prec );
    int restemptemp = -1;

    if ( (test==1) && CCLUSTER_EXP_NUM_T0(meta) )
        restemptemp = tstar_D0_test(cache, pshifted, &nbCoeffcomputed, d, depth, prec, meta );
    
    if ( (test==0) && CCLUSTER_EXP_NUM_T1(meta) ) {
        realApp_add(sum, abscoeff0, abscoeffn, prec);
        test = realApp_soft_compare( abscoeff1, sum, prec );
        if (test==1)
            restemptemp = tstar_N1_test(cache, pshifted, &nbCoeffcomputed, d, depth, prec, meta );
    }
    
    if (restemptemp==-1) {
        completeTaylorShift(cache, pshifted, nbCoeffcomputed, d, prec, meta );
        compApp_poly_set(pApprox, pshifted);
    }
    
    realRat_clear(radtothen);
    realApp_clear(sum);
    compApp_clear(c);
    realApp_clear(abscoeff0);
    realApp_clear(abscoeff1);
    realApp_clear(abscoeffn);
    compApp_poly_clear(pshifted);
    
    return restemptemp;
    
}
        
        
void completeTaylorShift(cacheApp_t cache, compApp_poly_t shiftedPol, int nbCoeffComputed, const compDsk_t d, slong prec, metadatas_t meta ){
    
    slong deg = cacheApp_getDegree(cache);
    slong deg2 = deg - (nbCoeffComputed-1);
    
    compApp_poly_set_length(shiftedPol, deg+1);
    
    if (deg2>1) {
        compApp_poly_t pp;
        realRat_t scaleFactor;
        realRat_init(scaleFactor);
        compApp_poly_init2(pp, deg2+1);
        
        tstar_getDerivative (pp, cache, prec, nbCoeffComputed , meta);
        tstar_taylor_shift_inplace( pp, d, prec, meta);
        
        realRat_pow_si(scaleFactor, compDsk_radiusref(d), (slong) nbCoeffComputed);
        for(int i=2; i<nbCoeffComputed; i++) 
            fmpz_mul_si( realRat_denref(scaleFactor), realRat_denref(scaleFactor), (slong) i);
        
        for(int i=nbCoeffComputed; i<=deg; i++) {
            
            fmpz_mul_si( realRat_denref(scaleFactor), realRat_denref(scaleFactor), (slong) i);
            if ((i-nbCoeffComputed)>=1)
                fmpz_mul_si( realRat_numref(scaleFactor), realRat_numref(scaleFactor), ((slong) (i-nbCoeffComputed)));
            
            compApp_mul_realRat( (shiftedPol->coeffs) + i,
                                (pp->coeffs) + (i-nbCoeffComputed),
                                scaleFactor, prec);
        }
        
        realRat_clear(scaleFactor);
        compApp_poly_clear(pp);
    }
}
    
int tstar_D0_test(cacheApp_t cache, 
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
    
    realRat_set_si(fac,1,1);
    realRat_mul(nrad, compDsk_radiusref(d), fac);
    compApp_set_compDsk( ball, compDsk_centerref(d), nrad, prec);
    compApp_set_compRat( c, compDsk_centerref(d), prec);
    compApp_sub(ballmc, ball, c, prec);
    compApp_set(ballmctoi, ballmc);
    realRat_set_si(fac,1,2);
    realRat_set(scaleFactor, compDsk_radiusref(d));
    
    compApp_abs(abscoeff0, (shiftedPol->coeffs)+0, prec);
    compApp_abs(abscoeff1, (shiftedPol->coeffs)+1, prec);
    realApp_set(sumofabscoeff, abscoeff1);
    
    realRat_mul( scaleFactor, scaleFactor, compDsk_radiusref(d));
    compApp_mul(ballmctoi, ballmctoi, ballmc, prec);
    
    slong i = 2;
    int dominates = 1;
    int res = 0;
    while ( (res==0)&&dominates && (i<=degree)) { 
//     while ( (res==0)&&dominates && (i<=2)) { 
        tstar_getDerivative (pp, cache, prec, i , meta);
        tstar_evaluate(pballCor, pp, ball, prec, meta, depth); /* pballCor <- (p^(i))(B)*/
        compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
        compApp_mul(pballCor, pballCor, ballmctoi, prec);      /* pballCor <- pballCor*(B-c)^i */
        compApp_abs(abscoeffi, pballCor, prec);
        realApp_add(approx, sumofabscoeff, abscoeffi, prec);
        
        res = realApp_soft_compare( abscoeff0, approx, prec );
        
        if (!(res==1)){
            res = 0;
            tstar_evaluate((shiftedPol->coeffs)+((int) i), pp, c, prec, meta, depth); /* i-thcoeff of shifted <- (p^(i))(c)*/
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), fac, prec);  /* i-thcoeff <- i-thcoeff / i!*/
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), scaleFactor, prec); /* i-thcoeff <- i-thcoeff *r^i */
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
    
int tstar_N1_test(cacheApp_t cache, compApp_poly_t shiftedPol,
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
//             tstar_evaluate(pballCor, pp, c, prec, meta, depth);    /* pballCor <- (p^(i))(c)*/
//             compApp_mul_realRat_in_place(pballCor, fac, prec);     /* pballCor <- pballCor / i!*/
//             compApp_abs(abscoeffi, pballCor, prec);
//             realApp_mul_realRat( abscoeffi, abscoeffi, scaleFactor, prec );
            
            tstar_evaluate((shiftedPol->coeffs)+((int) i), pp, c, prec, meta, depth); /* i-thcoeff of shifted <- (p^(i))(c)*/
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), fac, prec);  /* i-thcoeff <- i-thcoeff / i!*/
            compApp_mul_realRat_in_place((shiftedPol->coeffs)+((int) i), scaleFactor, prec); /* i-thcoeff <- i-thcoeff *r^i */
            compApp_abs(abscoeffi, (shiftedPol->coeffs)+((int) i), prec);
            
            realApp_add(sumofabscoeff, sumofabscoeff, abscoeffi, prec);
            dominates = realApp_soft_compare( abscoeff1, sumofabscoeff, prec );
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
//         printf("&&depth: %d, success of N1 non-discarding predicate, order of evaluation: %d, dominates: %d, prec: %d\n", depth, (int) i, dominates, (int) prec);
        return 1;
    }
    else return -1;
}
#endif