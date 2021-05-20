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
    
    compRat_t center;
    realRat_t radInf, radSup;
    compRat_init(center);
    realRat_init(radInf);
    realRat_init(radSup);
    cauchyTest_res cres;
    cres.appPrec = CCLUSTER_DEFAULT_PREC;
    
    compRat_zero(center);
    /* initialize double exponential sieve */
    realRat_set_si(upperBound, 2, 1);
    cres = cauchyTest_probabilistic_counting_test( center, upperBound, cache, cacheCau, cres.appPrec, meta, 0);
    
    if (metadatas_getVerbo(meta)>=3) {
        printf("#---cauchy root bound algorithm: \n");
        printf("#------ upperBound: "); realRat_print(upperBound); printf("\n");
        printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
    }
    /* confirm result with deterministic counting test */
    if (cres.nbOfSol==cacheCauchy_degreeref(cacheCau)) {
        realRat_set(radInf, upperBound);
        realRat_mul_si(radInf, radInf, 4);
        realRat_div_ui(radInf, radInf, 3);
        realRat_mul_si(radSup, upperBound, 6);
        cres = cauchyTest_deterministic_counting_test_for_newton( center,
                                                                  radInf, radSup, cacheCauchy_degreeref(cacheCau),
                                                                  cache, cacheCau, cres.appPrec, meta, 0);
        if (cres.nbOfSol)
                realRat_set(upperBound, radInf);
        if (metadatas_getVerbo(meta)>=3) {
            printf("#------ upperBound: "); realRat_print(upperBound); printf("\n");
            printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
        }
    }
    
    
    while (cres.nbOfSol!=cacheCauchy_degreeref(cacheCau)) {
        realRat_mul(upperBound, upperBound, upperBound);
        cres = cauchyTest_probabilistic_counting_test( center, upperBound, cache, cacheCau, cres.appPrec, meta, 0);
        if (metadatas_getVerbo(meta)>=3) {
            printf("#------ upperBound: "); realRat_print(upperBound); printf("\n");
            printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
        }
        /* confirm result with deterministic counting test */
        if (cres.nbOfSol==cacheCauchy_degreeref(cacheCau)) {
            realRat_set(radInf, upperBound);
            realRat_mul_si(radInf, radInf, 4);
            realRat_div_ui(radInf, radInf, 3);
            realRat_mul_si(radSup, upperBound, 6);
            cres = cauchyTest_deterministic_counting_test_for_newton( center,
                                                                      radInf, radSup, cacheCauchy_degreeref(cacheCau),
                                                                      cache, cacheCau, cres.appPrec, meta, 0);
            if (cres.nbOfSol)
                realRat_set(upperBound, radInf);
            if (metadatas_getVerbo(meta)>=3) {
                printf("#------ upperBound: "); realRat_print(upperBound); printf("\n");
                printf("#------ cres:       %d, %ld\n", cres.nbOfSol, cres.appPrec);
            }
        }
    }
    
    
    
    compRat_clear(center);
    realRat_clear(radInf);
    realRat_clear(radSup);
    
}
