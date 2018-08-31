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

#include "doubRealApp.h"
#include <fenv.h> 

void doubRealApp_swap(doubRealApp_t x, doubRealApp_t y) {
    double t1,t2;
    t1=x->low;
    t2=x->upp;
    x->low=y->low;
    x->upp=y->upp;
    y->low=t1;
    y->upp=t2;
}

void doubRealApp_set_realApp   (doubRealApp_t y, const realApp_t x ) { 
    double rad = mag_get_d( arb_radref(x) );
    y->upp= arf_get_d( arb_midref(x), ARF_RND_CEIL) + rad; 
    y->low= -arf_get_d( arb_midref(x), ARF_RND_FLOOR) + rad;
}

void doubRealApp_get_realApp   (realApp_t y, const doubRealApp_t x ){
    double rad = (x->upp + x->low)/2;
    double midupp = -(x->low) + rad;
    double midlow = x->low - rad;
    rad = rad + (midupp + midlow);
    
    arf_set_d( arb_midref(y), midupp);
    mag_set_d( arb_radref(y), rad);
}

// void doubRealApp_get_mid_doubRealApp(doubRealApp_t y, const doubRealApp_t x){
//     y->mid = x->mid;
//     y->rad = x->rad;
//     if (x->rad == 0.) 
//         y->rad = 0.;
//     else
//         y->rad = y->mid - nextafter(y->mid,-INFINITY);
// }
// 
void doubRealApp_add   (doubRealApp_t res, const doubRealApp_t x,  const doubRealApp_t y) { 
    res->upp = x->upp + y->upp;
    res->low = x->low + y->low;
}

void doubRealApp_sub   (doubRealApp_t res, const doubRealApp_t x,  const doubRealApp_t y) {
    double low = x->low + y->upp;
    res->upp = x->upp + y->low;
    res->low = low;
}
void doubRealApp_neg(doubRealApp_t res, const doubRealApp_t x){
    double low = x->upp;
    res->upp = x->low;
    res->low = low;
}

void doubRealApp_sqr (doubRealApp_t res, const doubRealApp_t x) {
    double upp, low;
    if (x->low <= 0) { /* case low of x >= 0 */
        upp = (x->upp)*(x->upp);
        low = (x->low)*(-x->low);
    }
    else if (x->upp <= 0) { /* case upp of x <= 0 */
        upp = (-x->low)*(-x->low);
        low = (-x->upp)*(x->upp);
    }
    else { /* case low of x<0 and upp of x > 0 */
        low=0.0;
        double tmp = (x->upp)*(x->upp);
        upp = (x->low)*(x->low);
        upp = CCLUSTER_MAX( upp, tmp );
    }
    res->upp = upp;
    res->low = low;
}

void _doubRealApp_mul   (doubRealApp_t res, const doubRealApp_t x,  const doubRealApp_t y) { 
    
    double upp, low;
    if (x->low <= 0) { /* case low of x >= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(y->upp);
            low = (x->low)*(-y->low);
//             fesetround(FE_DOWNWARD);
//             double low2 = (-(x->low))*(-(y->low));
//             low = -low2;
//             fesetround(FE_UPWARD);
//             printf("test: low = %0.23lf\n", low2);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(y->upp);
            low = (x->upp)*(y->low);
//             fesetround(FE_DOWNWARD);
//             double low2 = (x->upp)*(-(y->low));
//             printf("test: low = %0.23lf\n", low2);
//             fesetround(FE_UPWARD);
        }
        else { /* case low of y<0 and upp of y > 0 */
            upp = (x->upp)*(y->upp);
            low = (x->upp)*(y->low);
//             fesetround(FE_DOWNWARD);
//             double low2 = (x->upp)*(-(y->low));
//             printf("test: low = %0.23lf\n", low2);
//             fesetround(FE_UPWARD);
        }
    }
    else if (x->upp <= 0) { /* case upp of x <= 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(-y->low);
            low = (x->low)*(y->upp);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(-y->low);
            low = (-x->upp)*(y->upp);
        }
        else { /* case low of y<0 and upp of y > 0 */
            upp = (-x->low)*(-y->low);
            low = (x->low)*(y->upp);
        }
    }
    else { /* case low of x<0 and upp of x > 0 */
        if (y->low <= 0) { /* case low of y >= 0 */
            upp = (x->upp)*(y->upp);
            low = (x->low)*(y->upp);
        }
        else if (y->upp <= 0){ /* case upp of y <= 0 */
            upp = (-x->low)*(-y->low);
            low = (x->upp)*(y->low);
        }
        else { /* case low of y<0 and upp of y > 0 */
            double tmp = (x->upp)*(y->low);
            low = (x->low)*(y->upp);
            low = CCLUSTER_MAX( low, tmp );
            tmp = (x->upp)*(y->upp);
            upp = (x->low)*(y->low);
            upp = CCLUSTER_MAX( upp, tmp );
        }
    }
    res->upp = upp;
    res->low = low;
}
// 
// void doubRealApp_addmul   (doubRealApp_t res, const doubRealApp_t x,  const doubRealApp_t y) { 
//     double temp = x->mid*y->mid; 
//     double diff = 0.;
//     res->rad= res->rad + CCLUSTER_ABS(x->mid)*(y->rad) + (x->rad)*CCLUSTER_ABS(y->mid) + (x->rad)*(y->rad);
//     if ((!(y->mid==0))&&(!(temp / y->mid == x->mid))){
//         diff = temp - nextafter(temp,-INFINITY);
//         res->rad= res->rad + diff;
//     }
//     temp = res->mid + temp;
//     if (!(temp - x->mid*y->mid == res->mid)){
//         diff = temp - nextafter(temp,-INFINITY);
// //         printf("doubRealApp_add: is exact: %d, value of diff: %.23lf\n", diff==0, diff);
//         res->rad= res->rad + diff;
//     }
//     res->mid = temp;
// }
// 
void doubRealApp_mul_si   (doubRealApp_t res, const doubRealApp_t x,  slong y) { 
    res->upp = y*x->upp;
    res->low = y*x->low;
}

void doubRealApp_mul_ui   (doubRealApp_t res, const doubRealApp_t x,  ulong y) {
    res->upp = y*x->upp;
    res->low = y*x->low;
}

void doubRealApp_inv(doubRealApp_t z, const doubRealApp_t x){
    if (doubRealApp_contains_zero(x)) {
        z->upp = +INFINITY;
        z->low = -INFINITY;
    }
    double temp = 1./x->upp;
    z->upp = 1./(-x->low);
    z->low = -temp;
}

// void doubRealApp_sqrt   (doubRealApp_t res, const doubRealApp_t x) { 
//     res->mid= sqrt(x->mid); 
//     res->rad= CCLUSTER_ABS(x->mid)*(y->rad) + (x->rad)*CCLUSTER_ABS(y->mid) + (x->rad)*(y->rad);
//     if ((!(y->mid==0))&&(!(res->mid / y->mid == x->mid))){
//         double diff = res->mid - nextafter(res->mid,-INFINITY);
// //         printf("doubRealApp_mul: is exact: %d, value of diff: %.23lf\n", diff==0, diff);
//         res->rad= res->rad + diff;
//     }
// }

void doubRealApp_fprint (FILE * file, const doubRealApp_t x){
//     realApp_t conv;
//     realApp_init(conv);
//     doubRealApp_get_realApp(conv, x);
    fprintf(file, "[ %.23lf, %.23lf ]", -(*doubRealApp_lowref(x)), *doubRealApp_uppref(x)); 
//     realApp_fprint(file, conv);
//     fprintf(file, "\n");
//     realApp_clear(conv);
}
