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

#include <string.h>
#include <stdio.h>
#include "geometry/compBox.h"
#include "numbers/realRat.h"
#define GLOBAL_STR_NAME "global"
#define INFINITY_STR_NAME "infinity"

int scan_degree( char * argv, int * target) {
    return sscanf(argv, "%d", target);
}

int scan_bitsize( char * argv, int * target) {
    return sscanf(argv, "%d", target);
}

int scan_power( char * argv, int * target) {
    return sscanf(argv, "%d", target);
}

int scan_initialBox( char * argv, compBox_t target ){
    
    if (strcmp( argv, GLOBAL_STR_NAME ) == 0)
        return 2;
    
    char * tokRe=NULL;
    char * tokIm=NULL;
    char * tokWi=NULL;
    
    tokRe = strtok (argv,",");
    tokIm = strtok (NULL,",");
    tokWi = strtok (NULL,",");
    if ( (tokRe==NULL)||(tokIm==NULL)||(tokWi==NULL) ) {
        printf("error in parsing initial box!\n");
        return 0;
    }
    
    if ( compBox_set_str_pretty(target, tokRe, tokIm, tokWi) == -1 ){
        printf("error in parsing initial box! %s %s %s\n", tokRe, tokIm, tokWi);
        return 0;
    }
    
    if ( realRat_sgn( compBox_bwidthref(target))<=0 ) {
        printf("error: width of initial box should be positive! %s\n", tokWi);
        return 0;
    }
        
    return 1;
//     char * tok=NULL;
//     char temp[1000];
//     char cReN[100]= "0", cImN[100]= "0", widN[100]= "100";
//     char cReD[100]= "1", cImD[100]= "1", widD[100]="1";
//     
//     sscanf(argv, "%s", temp);
//     tok = strtok (temp,",");
//     sprintf(cReN, "%s", tok);
//     tok = strtok (NULL,",");
//     sprintf(cReD, "%s", tok);
//     tok = strtok (NULL,",");
//     sprintf(cImN, "%s", tok);
//     tok = strtok (NULL,",");
//     sprintf(cImD, "%s", tok);
//     tok = strtok (NULL,",");
//     sprintf(widN, "%s", tok);
//     tok = strtok (NULL,",");
//     sprintf(widD, "%s", tok);
//     tok = strtok (NULL,",");
//     
//     if (compBox_set_str(target, cReN, cReD, cImN, cImD, widN, widD, 10)==-1){
//         printf("error in parsing initial box! %s %s %s %s %s %s\n", cReN, cReD, cImN, cImD, widN, widD);
//         return 0;
//     }
//     
//     return 1;
}

int scan_initialInterval( char * argv, compBox_t target ){
    
    if (strcmp( argv, GLOBAL_STR_NAME ) == 0)
        return 2;
    
    char * tokRe=NULL;
//     char * tokIm=NULL;
    char * tokWi=NULL;
    
    tokRe = strtok (argv,",");
//     tokIm = strtok (NULL,",");
    tokWi = strtok (NULL,",");
    if ( (tokRe==NULL)||(tokWi==NULL) ) {
        printf("error in parsing initial interval!\n");
        return 0;
    }
    
    if ( compBox_set_str_pretty(target, tokRe, "0", tokWi) == -1 ){
        printf("error in parsing initial interval! %s %s\n", tokRe, tokWi);
        return 0;
    }
    
    if ( realRat_sgn( compBox_bwidthref(target))<=0 ) {
        printf("error: width of initial interval should be positive! %s\n", tokWi);
        return 0;
    }
        
    return 1;
}

int scan_epsilon( char * argv, realRat_t target ){
    
    if ( (strcmp( argv, INFINITY_STR_NAME ) == 0)
         ||(strcmp( argv, "inf" ) == 0) 
         ||(strcmp( argv, "+inf" ) == 0) ) {
        realRat_set_si(target, 1,0);
        return 2;
    }
    
    if ( realRat_set_str_pretty(target, argv) == 0 ) {
        if (realRat_is_zero(target)) {
            printf("error epsilon should not be 0 ! %s\n", argv);
            return 0;
        }
        if (realRat_sgn( target ) < 0) {
            if ( !(realRat_is_den_one(target)) ) {
            printf("error epsilon should not be either positive, or a negative integer ! %s\n", argv);
            return 0;
            }
            slong p;
            realRat_set_si(target, 2,1);
            sscanf(argv, "%ld", &p);
            realRat_pow_si(target, target, p);
        }
        return 1;
    }
    
    printf("error in parsing epsilon ! %s\n", argv);
    
    return 0;
    
//     char * tok=NULL;
//     char temp[1000];
//     char epsN[100]="1", epsD[100]="100";
//     
//     sscanf(argv, "%s", temp);
//     tok = strtok (temp,",");
//     sprintf(epsN, "%s", tok);
//     tok = strtok (NULL,",");
//     
//     
//     if (tok == NULL){
//         int p;
//         sscanf(epsN, "%d", &p);
//         realRat_set_si(target, 2,1);
//         realRat_pow_si(target, target, p);
//     }
//     else {
//         sprintf(epsD, "%s", tok);
//         if (realRat_set_str(target, epsN, epsD, 10)==-1){
//             printf("error in parsing epsilon ! %s %s\n", epsN, epsD);
//             return 0;
//         }
//     }
//     
//     return 1;
        
}

int scan_output( char * argv, int * target ) {
    
    if ( (strcmp( argv, "g" ) == 0)
         ||(strcmp( argv, "G" ) == 0) 
         ||(strcmp( argv, "Gnuplot" ) == 0)
         ||(strcmp( argv, "gnuplot" ) == 0) ) {
        *target = -2;
        return 1;
    } else if ( (strcmp( argv, "gs" ) == 0)
         ||(strcmp( argv, "GS" ) == 0) 
         ||(strcmp( argv, "GnuplotS" ) == 0)
         ||(strcmp( argv, "gnuplotS" ) == 0) ) {
        *target = -3;
        return 1;
    } else if ( (strcmp( argv, "r" ) == 0)
         ||(strcmp( argv, "R" ) == 0) 
         ||(strcmp( argv, "rational" ) == 0)
         ||(strcmp( argv, "Rational" ) == 0) ) {
        *target = -1;
        return 1;
    }
    else {
        sscanf(argv, "%d", target);
        return 1;
    }
}

int scan_strategy(char * argv, char * target ) {
    return sscanf(argv, "%s", target);
}

int scan_nbthreads(char * argv, int * target ) {
    return sscanf(argv, "%d", target);
}

int scan_verbosity(char * argv, int * target ) {
    return sscanf(argv, "%d", target);
}
