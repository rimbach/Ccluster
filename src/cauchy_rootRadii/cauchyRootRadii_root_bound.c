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

void cauchyRootRadii_root_bound( realRat_t upperBound,
                                 cacheCauchy_t cacheCau,
                                 cacheApp_t cache,
                                 metadatas_t meta ) {
    
    compDsk_t Delta;
    compDsk_init(Delta);
    cauchyTest_res cres;
    cres.appPrec = CCLUSTER_DEFAULT_PREC;
    
    compRat_zero(compDsk_centerref(Delta));
    /* initialize double exponential sieve */
    realRat_set_si(compDsk_radiusref(Delta), 2, 1);
    
    cres = cauchyTest_deterministic_verification( Delta, cacheCauchy_degreeref(cacheCau), cacheCauchy_isoRatioref(cacheCau), 
                                                  cache, cacheCau, cres.appPrec, meta, 0);
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("#---cauchy root bound algorithm: \n");
        printf("#-------cauchy test res %d\n", cres.nbOfSol);
        printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
        printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
    }
    
    while (cres.nbOfSol!=cacheCauchy_degreeref(cacheCau)) {
        realRat_mul(compDsk_radiusref(Delta), compDsk_radiusref(Delta), compDsk_radiusref(Delta));
        cres = cauchyTest_deterministic_verification( Delta, cacheCauchy_degreeref(cacheCau), cacheCauchy_isoRatioref(cacheCau), 
                                                      cache, cacheCau, cres.appPrec, meta, 0);
        if (metadatas_getVerbo(meta)>=3) {
            printf("#-------cauchy test res %d\n", cres.nbOfSol);
            printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
            printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
        }
        
    }
    
    
    realRat_set(upperBound, compDsk_radiusref(Delta));
    compDsk_clear(Delta);
    
}

/* OLD Version */
// void cauchyRootRadii_root_bound( realRat_t upperBound,
//                                  cacheCauchy_t cacheCau,
//                                  cacheApp_t cache,
//                                  metadatas_t meta ) {
//     
//     compDsk_t Delta;
//     compDsk_init(Delta);
//     cauchyTest_res cres;
//     cres.appPrec = CCLUSTER_DEFAULT_PREC;
//     
//     compRat_zero(compDsk_centerref(Delta));
//     /* initialize double exponential sieve */
//     realRat_set_si(compDsk_radiusref(Delta), 2, 1);
//     
//     cres = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//     
//     if (metadatas_getVerbo(meta)>=3) {
//         printf("#---cauchy root bound algorithm: \n");
//         printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
//         printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
//     }
//     /* confirm result with deterministic counting test */
//     if (cres.nbOfSol==cacheCauchy_degreeref(cacheCau)) {
//         cres = cauchyTest_deterministic_counting( Delta, cacheCauchy_isoRatioref(cacheCau), cache, cacheCau, cres.appPrec, meta, 0);
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#-------cauchy test res %d\n", cres.nbOfSol);
//             printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
//             printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
//         }
//     }
//     
//     
//     while (cres.nbOfSol!=cacheCauchy_degreeref(cacheCau)) {
//         realRat_mul(compDsk_radiusref(Delta), compDsk_radiusref(Delta), compDsk_radiusref(Delta));
//         cres = cauchyTest_probabilistic_counting( Delta, cache, cacheCau, cres.appPrec, meta, 0);
//         if (metadatas_getVerbo(meta)>=3) {
//             printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
//             printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
//         }
//         /* confirm result with deterministic counting test */
//         if (cres.nbOfSol==cacheCauchy_degreeref(cacheCau)) {
//             cres = cauchyTest_deterministic_counting( Delta, cacheCauchy_isoRatioref(cacheCau), cache, cacheCau, cres.appPrec, meta, 0);
//             if (metadatas_getVerbo(meta)>=3) {
//                 printf("#-------cauchy test res %d\n", cres.nbOfSol);
//                 printf("#------ upperBound: "); realRat_print(compDsk_radiusref(Delta)); printf("\n");
//                 printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
//             }
//         
//         }
//         
//     }
//     
//     
//     realRat_set(upperBound, compDsk_radiusref(Delta));
//     compDsk_clear(Delta);
//     
// }
