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

#include "fpci_poly.h"

/* Memory management */
void fpci_poly_init(fpci_poly_t poly){
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void fpci_poly_init2(fpci_poly_t poly, slong len){
    fpci_poly_init(poly);
    fpci_poly_fit_length(poly, len);
}

void fpci_poly_clear(fpci_poly_t poly){
    fpri_free(poly->coeffs);
}

void fpci_poly_fit_length(fpci_poly_t poly, slong len){
    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = fpri_realloc(poly->coeffs, len * sizeof(fpci_struct));
        poly->alloc = len;
    }
}

void _fpci_poly_set_length(fpci_poly_t poly, slong len){
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            fpci_zero(poly->coeffs + i);
    }

    poly->length = len;
}

void _fpci_poly_normalise(fpci_poly_t poly){
    slong i;

    for (i = poly->length - 1;
        (i >= 0) && fpci_is_zero(poly->coeffs + i); i--);

    poly->length = i + 1;
}

void fpci_poly_set (fpci_poly_t dest, const fpci_poly_t src){
    slong i;
    slong len = fpci_poly_length(src);

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_acb_poly (fpci_poly_t dest, const acb_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_acb(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_fpri (fpci_poly_t dest, const fpri_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_fpri(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

void fpci_poly_set_arb_poly (fpci_poly_t dest, const arb_poly_t src){
    slong i;
    slong len = src->length;

    fpci_poly_fit_length(dest, len);
    
    for (i = 0; i < len; i++)
        fpci_set_arb(dest->coeffs + i, src->coeffs + i);
    
    _fpci_poly_set_length(dest, len);
}

/* evaluation */
void _fpci_poly_evaluate2_horner(fpci_t y, fpci_t z, fpci_srcptr poly, slong len, const fpci_t x) {
    if (len == 0) {
        fpci_zero(y);
        fpci_zero(z);
    }
    else if (len == 1) {
        fpci_set(y, poly + 0);
        fpci_zero(z);
    }
    else if (fpci_is_zero(x)) {
        fpci_set(y, poly + 0);
        fpci_set(z, poly + 1);
    }
    else if (len == 2) {
        fpci_mul(y, x, poly + 1);
        fpci_add(y, y, poly + 0);
        fpci_set(z, poly + 1);
    }
    else {
        fpci_t t, u, v;
        slong i;

        fpci_set(u, poly + len - 1);
        fpci_zero(v);

        for (i = len - 2; i >= 0; i--) {
            fpci_mul(t, v, x);
            fpci_add(v, u, t);
            fpci_mul(t, u, x);
            fpci_add(u, t, poly + i);
        }

        fpci_swap(y, u);
        fpci_swap(z, v);
    }
}
void fpci_poly_evaluate2_horner(fpci_t y, fpci_t z, const fpci_poly_t p, const fpci_t x) {
    _fpci_poly_evaluate2_horner(y, z, p->coeffs, p->length, x);
}

void _fpci_poly_evaluate2_rectangulart(fpci_t y, fpci_t z, fpci_srcptr poly, slong len, const fpci_t x){
    
    slong i,j, m, r;
    fpci_ptr xs, xss;
    fpci_t s1, s2, temp1, temp2;
    
    fpci_zero(y);
    fpci_zero(z);
    
    m = n_sqrt(len) + 1;
    r = (len / m) + 1;
    xs  = _fpci_vec_init(m);
    xss = _fpci_vec_init(r);
    
    /* compute the m+1 first powers x^0, x, x^2, ..., x^(m-1) of x */
    fpci_one(xs + 0);
    fpci_set(xs + 1, x);
    for(i=2; i<m; i++)
        fpci_mul(xs+i, xs+(i-1), x);
    /* compute the first powers x^0, x^m, x^(2m), x^(3m), ..., x^((r-1)m) of x */
    fpci_one(xss + 0);
    fpci_mul(xss+1, xs+(m-1), x);
    for(i=2; i < r; i++)
       fpci_mul(xss+i, xss+(i-1), xss+1);
    
    /* loop on r-2 first blocks of m coeffs */
    for (j=0; j<r-2; j++) {
        fpci_set(s1, poly +  j*m    );
        fpci_set(s2, poly + (j*m +1));
        fpci_mul_si(s2, s2, j*m + 1 );
        for (i=1; i<m; i++) {
//             fpci_mul(temp1, poly + (j*m + i), xs + i);
//             fpci_mul(temp2, poly + (j*m + i + 1), xs + i);
//             _fpri_mul(fpci_realref(temp1), fpci_realref(poly + (j*m + i)), fpci_realref(xs + i));
//             _fpri_mul(fpci_imagref(temp1), fpci_realref(poly + (j*m + i)), fpci_imagref(xs + i));    
//             _fpri_mul(fpci_realref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_realref(xs + i));
//             _fpri_mul(fpci_imagref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_imagref(xs + i));
            
//             _fpri_mul(fpci_realref(temp1), fpci_realref(xs + i), fpci_realref(poly + (j*m + i))     );
//             _fpri_mul(fpci_imagref(temp1), fpci_imagref(xs + i), fpci_realref(poly + (j*m + i))     );    
//             _fpri_mul(fpci_realref(temp2), fpci_realref(xs + i), fpci_realref(poly + (j*m + i + 1)) );
//             _fpri_mul(fpci_imagref(temp2), fpci_imagref(xs + i), fpci_realref(poly + (j*m + i + 1)) );

             fpri_srcptr a=fpci_realref(xs + i)               ;
             fpri_srcptr b=fpci_imagref(xs + i)               ;
             fpri_srcptr c=fpci_realref(poly + (j*m + i))     ;
             fpri_srcptr d=fpci_realref(poly + (j*m + i + 1)) ;
             fpri_ptr e=fpci_realref(temp1)                ;
             fpri_ptr f=fpci_imagref(temp1)                ;
             fpri_ptr g=fpci_realref(temp2)                ;
             fpri_ptr h=fpci_imagref(temp2)                ;
            if ( ( (a)->low <= 0 ) && ( (b)->low <= 0 ) ) {
                if ( (c)->low <=0 ){
                    e->low = ( a->low)*(-c->low);
                    e->upp = ( a->upp)*( c->upp);
                    f->low = ( b->low)*(-c->low);
                    f->upp = ( b->upp)*( c->upp);
                } else if (c->upp <= 0){
                    e->low = ( a->upp)*( c->low);
                    e->upp = (-a->low)*( c->upp);
                    f->low = ( b->upp)*( c->low);
                    f->upp = (-b->low)*( c->upp);
                } else {
                    e->low = ( a->upp)*( c->low);
                    e->upp = ( a->upp)*( c->upp);
                    f->low = ( b->upp)*( c->low);
                    f->upp = ( b->upp)*( c->upp);
                } 
                if ( (d)->low <=0 ){
                    g->low = ( a->low)*(-d->low);
                    g->upp = ( a->upp)*( d->upp);
                    h->low = ( b->low)*(-d->low);
                    h->upp = ( b->upp)*( d->upp);
                } else if (d->upp <= 0){
                    g->low = ( a->upp)*( d->low);
                    g->upp = (-a->low)*( d->upp);
                    h->low = ( b->upp)*( d->low);
                    h->upp = (-b->low)*( d->upp);
                } else {
                    g->low = ( a->upp)*( d->low);
                    g->upp = ( a->upp)*( d->upp);
                    h->low = ( b->upp)*( d->low);
                    h->upp = ( b->upp)*( d->upp);
                }
            } else if ( ( (a)->upp <= 0 ) && ( (b)->upp <= 0 ) ) {
                if ( (c)->low <=0 ){
                    e->low = ( a->low)*( c->upp);
                    e->upp = ( a->upp)*(-c->low);
                    f->low = ( b->low)*( c->upp);
                    f->upp = ( b->upp)*(-c->low);
                } else if (c->upp <= 0){
                    e->low = (-a->upp)*( c->upp);
                    e->upp = (-a->low)*(-c->low);
                    f->low = (-b->upp)*( c->upp);
                    f->upp = (-b->low)*(-c->low);
                } else {
                    e->low = ( a->low)*( c->upp);
                    e->upp = (-a->low)*(-c->low);
                    f->low = ( b->low)*( c->upp);
                    f->upp = (-b->low)*(-c->low);
                } 
                if ( (d)->low <=0 ){
                    g->low = ( a->low)*( d->upp);
                    g->upp = ( a->upp)*(-d->low);
                    h->low = ( b->low)*( d->upp);
                    h->upp = ( b->upp)*(-d->low);
                } else if (d->upp <= 0){
                    g->low = (-a->upp)*( d->upp);
                    g->upp = (-a->low)*(-d->low);
                    h->low = (-b->upp)*( d->upp);
                    h->upp = (-b->low)*(-d->low);
                } else {
                    g->low = ( a->low)*( d->upp);
                    g->upp = (-a->low)*(-d->low);
                    h->low = ( b->low)*( d->upp);
                    h->upp = (-b->low)*(-d->low);
                }
            } else {
                _fpri_mul(e, a, c );
                _fpri_mul(f, b, c );    
                _fpri_mul(g, a, d );
                _fpri_mul(h, b, d );
            }
            
            fpci_mul_si( temp2, temp2, j*m + i + 1 );
            fpci_add(s1, s1, temp1);
            fpci_add(s2, s2, temp2);
        }
        fpci_mul(temp1, s1, xss + j);
        fpci_mul(temp2, s2, xss + j);
        fpci_add(y, y, temp1);
        fpci_add(z, z, temp2);
    }
    
    /* (r-1)-th block */
    j=r-2;
    fpci_set(s1, poly +  j*m    );
    fpci_set(s2, poly + (j*m +1));
    fpci_mul_si(s2, s2, j*m + 1 );
    for (i=1; i< m-1; i++) {
//          fpci_mul(temp1, poly + (j*m + i), xs + i);
         _fpri_mul(fpci_realref(temp1), fpci_realref(poly + (j*m + i)), fpci_realref(xs + i));
         _fpri_mul(fpci_imagref(temp1), fpci_realref(poly + (j*m + i)), fpci_imagref(xs + i));
//          fpci_mul(temp2, poly + (j*m + i + 1), xs + i);
         _fpri_mul(fpci_realref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_realref(xs + i));
         _fpri_mul(fpci_imagref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_imagref(xs + i));
         fpci_mul_si( temp2, temp2, j*m + i + 1);
         fpci_add(s1, s1, temp1);
         fpci_add(s2, s2, temp2);
    }
    fpci_mul(temp1, poly + (j*m + m-1), xs + m-1);
    fpci_add(s1, s1, temp1);
//     fpci_addmul(s1, poly + (j*m + m-1), xs + m-1);
    if ((j*m + m)<len) {
        fpci_mul(temp2, poly + (j*m + m), xs + m-1);
        fpci_mul_si( temp2, temp2, j*m + m );
        fpci_add(s2, s2, temp2);
    }
    fpci_mul(temp1, s1, xss + j);
    fpci_mul(temp2, s2, xss + j);
    fpci_add(y, y, temp1);
    fpci_add(z, z, temp2);
    
    /* (r)-th block */
    j=r-1;
    fpci_set(s1, poly +  j*m    );
    if ((j*m + 1)<len) {
        fpci_set(s2, poly + (j*m +1));
        fpci_mul_si(s2, s2, j*m + 1 );
    }
    else 
        fpci_zero(s2);
    for (i=1; (j*m + i) < (len-1); i++) {
//          fpci_mul(temp1, poly + (j*m + i), xs + i);
         _fpri_mul(fpci_realref(temp1), fpci_realref(poly + (j*m + i)), fpci_realref(xs + i));
         _fpri_mul(fpci_imagref(temp1), fpci_realref(poly + (j*m + i)), fpci_imagref(xs + i));
//          fpci_mul(temp2, poly + (j*m + i + 1), xs + i);
         _fpri_mul(fpci_realref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_realref(xs + i));
         _fpri_mul(fpci_imagref(temp2), fpci_realref(poly + (j*m + i + 1)), fpci_imagref(xs + i));
         fpci_mul_si( temp2, temp2, j*m + i + 1 );
         fpci_add(s1, s1, temp1);
         fpci_add(s2, s2, temp2);
    }
    fpci_mul(temp1, poly + (len-1), xs + m-1);
    fpci_add(s1, s1, temp1);
    fpci_mul(temp1, s1, xss + j);
    fpci_mul(temp2, s2, xss + j);
    fpci_add(y, y, temp1);
    fpci_add(z, z, temp2);
    
    
    _fpci_vec_clear( xs,  m);
    _fpci_vec_clear( xss, r);
}

void fpci_poly_evaluate2_rectangulart(fpci_t y, fpci_t z, const fpci_poly_t p, const fpci_t x) {
    _fpci_poly_evaluate2_rectangulart(y, z, p->coeffs, p->length, x);
}

void _fpci_poly_evaluate2_rectangular(fpci_t y, fpci_t z, fpci_srcptr poly, slong len, const fpci_t x){
    
    slong i,j, m, r;
    fpci_ptr xs, xss;
    fpci_t s1, s2, temp;
    
    fpci_zero(y);
    fpci_zero(z);
    
    m = n_sqrt(len) + 1;
    r = (len / m) + 1;
    xs  = _fpci_vec_init(m);
    xss = _fpci_vec_init(r);
    
    /* compute the m+1 first powers x^0, x, x^2, ..., x^(m-1) of x */
    fpci_one(xs + 0);
    fpci_set(xs + 1, x);
    for(i=2; i<m; i++)
        fpci_mul(xs+i, xs+(i-1), x);
    /* compute the first powers x^0, x^m, x^(2m), x^(3m), ..., x^((r-1)m) of x */
    fpci_one(xss + 0);
    fpci_mul(xss+1, xs+(m-1), x);
    for(i=2; i < r; i++)
       fpci_mul(xss+i, xss+(i-1), xss+1);
    
    /* loop on r-2 first blocks of m coeffs */
    for (j=0; j<r-2; j++) {
        fpci_set(s1, poly +  j*m    );
        fpci_set(s2, poly + (j*m +1));
        fpci_mul_si(s2, s2, j*m + 1 );
        for (i=1; i<m; i++) {
            fpci_mul(temp, poly + (j*m + i), xs + i);
            fpci_add(s1, s1, temp);
            fpci_mul(temp, poly + (j*m + i + 1), xs + i);
            fpci_mul_si( temp, temp, j*m + i + 1 );
            fpci_add(s2, s2, temp);
        }
        fpci_mul(temp, s1, xss + j);
        fpci_add(y, y, temp);
        fpci_mul(temp, s2, xss + j);
        fpci_add(z, z, temp);
    }
    
    /* (r-1)-th block */
    j=r-2;
    fpci_set(s1, poly +  j*m    );
    fpci_set(s2, poly + (j*m +1));
    fpci_mul_si(s2, s2, j*m + 1 );
    for (i=1; i< m-1; i++) {
         fpci_mul(temp, poly + (j*m + i), xs + i);
         fpci_add(s1, s1, temp);
         fpci_mul(temp, poly + (j*m + i + 1), xs + i);
         fpci_mul_si( temp, temp, j*m + i + 1);
         fpci_add(s2, s2, temp);
    }
    fpci_mul(temp, poly + (j*m + m-1), xs + m-1);
    fpci_add(s1, s1, temp);
    if ((j*m + m)<len) {
        fpci_mul(temp, poly + (j*m + m), xs + m-1);
        fpci_mul_si( temp, temp, j*m + m );
        fpci_add(s2, s2, temp);
    }
    fpci_mul(temp, s1, xss + j);
    fpci_add(y, y, temp);
    fpci_mul(temp, s2, xss + j);
    fpci_add(z, z, temp);
    
    /* (r)-th block */
    j=r-1;
    fpci_set(s1, poly +  j*m    );
    if ((j*m + 1)<len) {
        fpci_set(s2, poly + (j*m +1));
        fpci_mul_si(s2, s2, j*m + 1 );
    }
    else 
        fpci_zero(s2);
    for (i=1; (j*m + i) < (len-1); i++) {
         fpci_mul(temp, poly + (j*m + i), xs + i);
         fpci_add(s1, s1, temp);
         fpci_mul(temp, poly + (j*m + i + 1), xs + i);
         fpci_mul_si( temp, temp, j*m + i + 1 );
         fpci_add(s2, s2, temp);
    }
    fpci_mul(temp, poly + (len-1), xs + m-1);
    fpci_add(s1, s1, temp);
    fpci_mul(temp, s1, xss + j);
    fpci_add(y, y, temp);
    fpci_mul(temp, s2, xss + j);
    fpci_add(z, z, temp);
    
    
    _fpci_vec_clear( xs,  m);
    _fpci_vec_clear( xss, r);
}

void fpci_poly_evaluate2_rectangular(fpci_t y, fpci_t z, const fpci_poly_t p, const fpci_t x) {
    _fpci_poly_evaluate2_rectangular(y, z, p->coeffs, p->length, x);
}

/* sparse evaluation */
slong fpci_init_sparse_eval( slong ** inNZC, const fpci_poly_t p ){
    slong ind;
    slong nbNZC=0;
    for (ind=0; ind < p->length; ind++)
        if ( ! fpci_is_zero( (p->coeffs) + ind ) )
            nbNZC++;
    *inNZC = (slong *) fpri_malloc (nbNZC*(sizeof(slong)));
    nbNZC=0;
    for (ind=0; ind < p->length; ind++)
        if ( ! fpci_is_zero( (p->coeffs) + ind ) ) {
            (*inNZC)[nbNZC] = ind;
            nbNZC++;
        }
    return nbNZC;
}

void fpci_clear_sparse_eval( slong ** inNZC, slong  nbNZC){
    fpri_free(*inNZC);
    *inNZC = NULL;
}

void fpci_poly_sparse_eval(fpci_t y, fpci_t z, const fpci_poly_t p, slong inNZC[], slong nbNZC, const fpci_t point){
    fpci_set( y, (p->coeffs) + 0);
    fpci_zero( z );
    
    fpci_t x, xp, mon;
//     fpci_init(x);
//     fpci_init(xp);
//     fpci_init(mon);
    
    slong ind=0;
    if (inNZC[ind] == 0)
        ind++;
    fpci_pow_ui( xp, point, inNZC[ind]-1 );
    fpci_mul(x, xp, point);
    
    while ( ind < nbNZC ){
        fpci_mul(mon, (p->coeffs) + inNZC[ind], x);
        fpci_add(y, y, mon);
//         fpci_addmul(y, (p->coeffs) + inNZC[ind], x);
        
        fpci_mul(mon, (p->coeffs) + inNZC[ind], xp);
        fpci_mul_si( mon, mon, inNZC[ind]);
        fpci_add(z, z, mon);
        
        ind++;
        
        if (ind < nbNZC) {
            fpci_pow_ui(mon, point, inNZC[ind] - inNZC[ind-1]);
            fpci_mul( x, x, mon);
            fpci_mul( xp, xp, mon);
        }
    }
    
    fpci_clear(x);
    fpci_clear(xp);
    fpci_clear(mon);
}
 
void fpci_poly_sparse_eval2(fpci_t fval, fpci_t fderval, const fpci_poly_t p, slong inNZC[], slong nbNZC, const fpci_t x){
        
    int log2deg = (int) ceil(log2( (inNZC[nbNZC-1]) - 1 ));
    fpci_ptr pows = (fpci_ptr) fpri_malloc (log2deg*sizeof(fpci_struct));
    fpci_t powp, powd;
//     fpci_init(powp);
//     fpci_init(powd);
    
    fpci_zero( fval );
    fpci_zero( fderval );
    
    fpci_zero(pows+0);
    fpci_set(pows+0, x); 
    
    for (int i=1; i<log2deg; i++){
//         fpci_zero(pows+i);
        fpci_sqr(pows+i, pows+(i-1)); 
    }
    
    for (slong i=0; i<nbNZC; i++){
//         printf("i: %ld, inNZC[i]: %ld\n", i, inNZC[i]);
        if ( inNZC[i]==0 ){
            fpci_zero(powd);
            fpci_one (powp);
        } else if ( inNZC[i]==1 ){ 
            fpci_one(powd);
            fpci_set (powp, x);
        } else if ( inNZC[i]==2 ){
            fpci_set (powd, pows + 0);
            fpci_set (powp, pows + 1);
        } else {
            fpci_one (powd);
            slong pow = inNZC[i] - 1;
            int indmax = (int) ceil(log2( pow ));
//             printf("pow: %ld, ind: %d\n", pow, indmax);
            for(int ind = 0; ind<indmax; ind++){
                if (pow%2)
                    fpci_mul( powd, powd, pows + ind);
                pow=pow>>1;
            }
            fpci_mul(powp, powd, x);
        }
        fpci_mul_si(powd, powd, inNZC[i]);
//         fpci_addmul(fval,    powp, (p->coeffs) + inNZC[i]);
//         fpci_addmul(fderval, powd, (p->coeffs) + inNZC[i]);
        fpci_t temp;
        fpci_mul(temp, powp, (p->coeffs) + inNZC[i]);
        fpci_add(fval, fval, temp);
        fpci_mul(temp, powd, (p->coeffs) + inNZC[i]);
        fpci_add(fderval, fderval, temp);
    }
    
//     for (int i=0; i<log2deg; i++)
//         compApp_clear(pows+i);
    
    fpri_free(pows); 
}
            
/* printing */
void fpci_poly_fprint (FILE * file, const fpci_poly_t x){
    slong i;
    slong len = x->length;
    if (len==0)
        fprintf(file, "degree: 0, coeffs: 0");
    else {
        fprintf(file, "degree: %ld, coeffs: ", len-1);
        for (i = 0; i <len; i++){ 
            fpci_fprint(file, x->coeffs+i);
            if (i<len-1)
                fprintf(file, "\n");
        }
    }
}
