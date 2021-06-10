/* ************************************************************************** */
/*  Copyright (C) 2021 Remi Imbach                                            */
/*                                                                            */
/*  This file is part of Ccluster.                                            */
/*                                                                            */
/*  Ccluster is free software: you can redistribute it and/or modify it under */
/*  the terms of the GNU Lesser General Public License (LGPL) as published    */
/*  by the Free Software Foundation; either version 2.1 of the License, or    */
/*  (at your option) any later version.  See <http://www.gnu.org/licenses/>.  */
/* ************************************************************************** */

#include "cauchy_rootRadii/cauchy_rootRadii.h"

// void cauchyRootRadii_root_radius( const compRat_t center,
//                                   realRat_t radInf, /* radInf < r_{d+1-m}(center, p) */
//                                   realRat_t radSup, /* radSup > r_{d+1-m}(center, p) */
//                                   realRat_t relativeError, /* want relativeError*radInf >= radSup */ 
//                                   slong nbOfRoots,
//                                   cacheCauchy_t cacheCau,
//                                   cacheApp_t cache,
//                                   metadatas_t meta ){
//     
//     /* 2^logRelError <= relativeError */
//     ulong logRelError = fmpz_flog_ui( realRat_numref( relativeError ), 2 ) - fmpz_clog_ui( realRat_denref( relativeError ), 2 );
//     slong precForMiddle = 2*logRelError;
//     
//     cauchyTest_res cres;
//     cres.appPrec = CCLUSTER_DEFAULT_PREC;
//     
//     int stop = 0;
//     realRat_t ratio, middle;
//     realRat_init(ratio);
//     realRat_init(middle);
//     
//     realApp_t mid;
//     realApp_init(mid);
//     
//     fmpz_t a,b,exp;
//     fmpz_init(a);
//     fmpz_init(b);
//     fmpz_init(exp);
//     
//     compDsk_t Delta;
//     compDsk_init(Delta);
//     compRat_set( compDsk_centerref(Delta), center );
//     
//     realRat_div(ratio, radSup, radInf);
//     stop = ( realRat_cmp( ratio, relativeError ) <= 0);
//     
//     
//     while (stop==0) {
//         
//         realRat_mul( middle, radInf, radSup );
//         realApp_set_realRat( mid, middle, precForMiddle );
//         realApp_sqrt(mid, mid, precForMiddle);
//         arb_get_interval_fmpz_2exp(a, b, exp, mid);
//         realRat_set_fmpz(middle, a);
//         fmpz_set_si(b, 2);
//         fmpz_pow_fmpz(b, exp);
//         fmpq_mul_fmpz(middle, middle, b);
//         
//         /* call probabilistic root counter */
//         realRat_set( compDsk_radiusref(Delta), middle );
//         cres = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//         if (cres.nbOfSol == -1) { /* there are roots in A(c, middle/(theta), middle*theta ) with theta = 4/3 */
//                                   /* r_{d+1-m}(center, p) > middle/(theta) */
//             realRat_mul( radInf, middle, cacheCauchy_isoRatioref( cacheCau ) );
//         } else if (cres.nbOfSol != nbOfRoots) { /* either the result is correct and r_{d+1-m}(center, p) > middle */
//                                                /* or     the result is not correct and there are roots in A(c, middle/(theta), middle*theta ) */
//                                                /* r_{d+1-m}(center, p) > middle/(theta) */
//             realRat_mul( radInf, middle, cacheCauchy_isoRatioref( cacheCau ) );
//         } else { /* need to check the result with deterministic counting */
//             cres = cauchyTest_deterministic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//             if (cres.nbOfSol == -1) { /* there are roots in A(c, 13/24*middle, 53/24*middle); */
//                                       /* r_{d+1-m}(center, p) > 13/24*middle */
//                 realRat_mul_si(radInf, middle, 13);
//                 realRat_div_ui(radInf, radInf, 24);
//         } else if (cres.nbOfSol != nbOfRoots) { /* necessarily cres.nbOfSol < nbOfRoots because the result is correct */
//                                                 /* r_{d+1-m}(center, p) > middle */
//         
//                 realRat_set( radInf, middle );
//         } else { /* cres.nbOfSol == nbOfRoots */
//                  /* r_{d+1-m}(center, p) < middle */
//                 realRat_set( radSup, middle ); 
//         }
//         /* check stop criterion */
//         realRat_div(ratio, radSup, radInf);
//         stop = ( realRat_cmp( ratio, relativeError ) <= 0);
//     }
//     
//     realRat_clear(ratio);
//     realRat_clear(middle);
//     realApp_clear(mid);
//     fmpz_clear(a);
//     fmpz_clear(b);
//     fmpz_clear(exp);
// }
