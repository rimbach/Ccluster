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

#include "doubCompApp.h"

doubCompApp_ptr _doubCompApp_vec_init(slong n){
//     slong i;
    doubCompApp_ptr v = (doubCompApp_ptr) flint_malloc(sizeof(doubCompApp) * n);

//     for (i = 0; i < n; i++)
//         doubCompApp_init(v + i);

    return v;
}

void _doubCompApp_vec_clear(doubCompApp_ptr v, slong n){
//     slong i;
//     for (i = 0; i < n; i++)
//         doubCompApp_clear(v + i);
    flint_free(v);
}

void doubCompApp_mul   ( doubCompApp_t z, const doubCompApp_t x, const doubCompApp_t y){ 
#define a doubCompApp_realref(x)
#define b doubCompApp_imagref(x)
#define c doubCompApp_realref(y)
#define d doubCompApp_imagref(y)
#define e doubCompApp_realref(z)
#define f doubCompApp_imagref(z)
    
    if ( doubRealApp_is_zero(b) ) {
        doubRealApp_mul(f, d, a);
        doubRealApp_mul(e, c, a);
    }
    else if ( doubRealApp_is_zero(d) ) {
        doubRealApp_mul(f, b, c);
        doubRealApp_mul(e, a, c);
    }
    else if ( doubRealApp_is_zero(a) ) {
        doubRealApp_mul(e, c, b);
        doubRealApp_mul(f, d, b);
        doubCompApp_mul_onei(z,z);
    }
    else if ( doubRealApp_is_zero(c) ) {
        doubRealApp_mul(e, a, d);
        doubRealApp_mul(f, b, d);
        doubCompApp_mul_onei(z,z);
    }
    else {
        if (x == y) { 
            /*squaring*/
            doubCompApp_t temp;
            doubRealApp_sqr(doubCompApp_realref(temp), a);
            doubRealApp_sqr(doubCompApp_imagref(temp), b);
            doubRealApp_sub(doubCompApp_realref(temp), doubCompApp_realref(temp), doubCompApp_imagref(temp));
            
            doubRealApp_mul(doubCompApp_imagref(temp), a, b);
            doubRealApp_mul_si(doubCompApp_imagref(temp), doubCompApp_imagref(temp), 2);
            
            doubCompApp_set(z, temp);
        }
        else {
            /* Gauss multiplication: e = ac-bd
                                     f = (a+b)(c+d) - ac - bd
                                     classical: f = ad + bc
            */
            doubRealApp_t t,u,v,w;
            
//             doubRealApp_mul(t, a, c);
//             doubRealApp_mul(u, b, d);
//             
//             doubRealApp_add(v, a, b);
//             doubRealApp_add(w, c, d);
//             doubRealApp_mul(f, v, w);
//             doubRealApp_sub(f, f, t);
//             doubRealApp_sub(f, f, u);
//             
//             doubRealApp_sub(e, t, u);
            
            /* classical multiplication: e = ac-bd
                                         f = ad + bc
                                         (a+ib)(c+id)=ac + iad + ibc +iibd
            */
            doubRealApp_mul(t, a, c);
            doubRealApp_mul(u, b, d);
            doubRealApp_mul(v, a, d);
            doubRealApp_mul(w, b, c);
            doubRealApp_sub(e, t, u);
            doubRealApp_add(f, v, w);
        }
    }
#undef a
#undef b
#undef c
#undef d
#undef e
#undef f    
}

void doubCompApp_fprint (FILE * file, const doubCompApp_t x){
    fprintf(file, "( ");
    doubRealApp_fprint (file, doubCompApp_realref(x));
    fprintf(file, " +i ");
    doubRealApp_fprint (file, doubCompApp_imagref(x));
    fprintf(file, " )");
}

/* DEPRECATED */
// #define a doubCompApp_realref(x)
// #define b doubCompApp_imagref(x)
// #define c doubCompApp_realref(y)
// #define d doubCompApp_imagref(y)
// #define e doubCompApp_realref(z)
// #define f doubCompApp_imagref(z)
// #define ar a->rad
// #define br b->rad
// #define cr c->rad
// #define dr d->rad
// void doubCompApp_sqr   ( doubCompApp_t z, const doubCompApp_t x ) {
//     double am, bm, er, fr, em;
//     am = a->mid;
//     bm = b->mid;
//     
//     er = 2*(am*ar + bm*br) + (ar*ar) + (br*br);
//     fr = 2*(bm*ar + am*br + br*br);
//     
// //     e->mid = am*am - bm*bm;
//     em = (bm)*(bm);
//     if ((!(bm==0))&&(!(em / bm == bm))){
//         double diff = em - nextafter(em,-INFINITY);
//         er= er + diff;
//     }
//     e->mid = (am)*(am);
//     if ((!(am==0))&&(!(e->mid / am == am))){
//         double diff = e->mid - nextafter(e->mid,-INFINITY);
//         er= er + diff;
//     }
//     e->mid = e->mid - em;
//     if (!(e->mid + em == e->mid)){
//         double diff = e->mid - nextafter(e->mid,-INFINITY);
//         er= er + diff;
//     }
//     
// //     f->mid = 2*(am)*(bm); 
//     f->mid = 2*(am)*(bm);
//     if ((!(bm==0))&&(!(f->mid / (2*bm) == bm))){
//         double diff = f->mid - nextafter(f->mid,-INFINITY);
//         fr= fr + diff;
//     }
//     
//     e->rad = er;
//     f->rad = fr;
// }
// #undef a
// #undef b
// #undef c
// #undef d
// #undef e
// #undef f
// #undef ar
// #undef br
// #undef cr
// #undef dr
