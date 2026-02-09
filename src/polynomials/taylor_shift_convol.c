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

#include "compApp_poly.h"

void _compApp_poly_taylor_shift_conv_pre(compApp_ptr p, realApp_t f, compApp_ptr t, slong len, slong prec){
    
    slong i, n = len - 1;
    
    arb_one(f);
    for (i = 2; i <= n; i++)
    {
        arb_mul_ui(f, f, i, prec);
        acb_mul_arb(p + i, p + i, f, prec);
    }

    _acb_poly_reverse(p, p, len, len);
    
    acb_one(t + n);
    for (i = n; i > 0; i--)
        acb_mul_ui(t + i - 1, t + i, i, prec);
    
}

void compApp_poly_taylor_shift_conv_pre(compApp_poly_t dest, const compApp_poly_t p, realApp_t f, compApp_ptr t, slong prec){
    compApp_poly_set(dest, p);
    _compApp_poly_taylor_shift_conv_pre(dest->coeffs, f, t, dest->length, prec);
}

void compApp_poly_taylor_shift_convolution(compApp_poly_t dest, const compApp_poly_t p, const compApp_t c, slong prec){
    
    realApp_t f;
    acb_ptr t;
    realApp_init(f);
    
    t = _acb_vec_init(p->length);
    
    compApp_poly_taylor_shift_conv_pre(dest, p, f, t, prec);
    _compApp_poly_taylor_shift_convolution(dest->coeffs, f, t, c, dest->length, prec);
    
    realApp_clear(f);
    _acb_vec_clear(t, p->length);
}

void _compApp_poly_taylor_shift_convolution(compApp_ptr p, realApp_t f, compApp_ptr t, const compApp_t c, slong len, slong prec){
    slong i, n = len - 1;
    acb_t d;
    arb_t fp;
    acb_ptr tp, u;
//     acb_ptr u;

    if (acb_is_zero(c) || len <= 1)
        return;

    tp = _acb_vec_init(len);
    u = _acb_vec_init(len);

    arb_init(fp);
    acb_init(d);

//     _compApp_poly_taylor_shift_conv_pre(p, f, t, len, prec);

    if (acb_equal_si(c, -1))
    {
        for (i = 1; i <= n; i += 2)
            acb_neg(tp + i, t + i);
    }
    else if (!acb_is_one(c))
    {
        acb_set(d, c);

        for (i = 1; i <= n; i++)
        {
            acb_mul(tp + i, t + i, d, prec);
            acb_mul(d, d, c, prec);
        }
    }
    else {
        for (i = 1; i <= n; i++)
            acb_set(tp + i, t + i);
    }
        

    _acb_poly_mullow(u, p, len, tp, len, len, prec);

    arb_mul(fp, f, f, prec);

    if (arb_bits(fp) > 0.25 * prec)
    {
        arb_inv(fp, fp, prec);
    }
    else
    {
        for (i = 0; i <= n; i++)
            acb_div_arb(u + i, u + i, fp, prec);

        arb_one(fp);
    }

    for (i = n; i >= 0; i--)
    {
        acb_mul_arb(p + i, u + n - i, fp, prec);
        arb_mul_ui(fp, fp, (i == 0) ? 1 : i, prec);
    }

    _acb_vec_clear(tp, len);
    _acb_vec_clear(u, len);

    arb_clear(fp);
    acb_clear(d);
}

#ifdef ARB_IN_FLINT

void
_arb_poly_taylor_shift_convolution(arb_ptr p, const arb_t c, slong len, slong prec)
{
    slong i, n = len - 1;
    arb_t f, d;
    arb_ptr t, u;

    if (arb_is_zero(c) || len <= 1)
        return;

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);

    arb_init(f);
    arb_init(d);

    arb_one(f);
    for (i = 2; i <= n; i++)
    {
        arb_mul_ui(f, f, i, prec);
        arb_mul(p + i, p + i, f, prec);
    }

    _arb_poly_reverse(p, p, len, len);

    arb_one(t + n);
    for (i = n; i > 0; i--)
        arb_mul_ui(t + i - 1, t + i, i, prec);

    if (arb_equal_si(c, -1))
    {
        for (i = 1; i <= n; i += 2)
            arb_neg(t + i, t + i);
    }
    else if (!arb_is_one(c))
    {
        arb_set(d, c);

        for (i = 1; i <= n; i++)
        {
            arb_mul(t + i, t + i, d, prec);
            arb_mul(d, d, c, prec);
        }
    }

    _arb_poly_mullow(u, p, len, t, len, len, prec);

    arb_mul(f, f, f, prec);

    if (arb_bits(f) > 0.25 * prec)
    {
        arb_inv(f, f, prec);
    }
    else
    {
        for (i = 0; i <= n; i++)
            arb_div(u + i, u + i, f, prec);

        arb_one(f);
    }

    for (i = n; i >= 0; i--)
    {
        arb_mul(p + i, u + n - i, f, prec);
        arb_mul_ui(f, f, (i == 0) ? 1 : i, prec);
    }

    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);

    arb_clear(f);
    arb_clear(d);
}

void
arb_poly_taylor_shift_convolution(arb_poly_t g, const arb_poly_t f,
        const arb_t c, slong prec)
{
    if (f != g)
        arb_poly_set_round(g, f, prec);

    _arb_poly_taylor_shift_convolution(g->coeffs, c, g->length, prec);
}

void
_acb_poly_taylor_shift_convolution(acb_ptr p, const acb_t c, slong len, slong prec)
{
    slong i, n = len - 1;
    acb_t d;
    arb_t f;
    acb_ptr t, u;

    if (acb_is_zero(c) || len <= 1)
        return;

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    arb_init(f);
    acb_init(d);

    arb_one(f);
    for (i = 2; i <= n; i++)
    {
        arb_mul_ui(f, f, i, prec);
        acb_mul_arb(p + i, p + i, f, prec);
    }

    _acb_poly_reverse(p, p, len, len);

    acb_one(t + n);
    for (i = n; i > 0; i--)
        acb_mul_ui(t + i - 1, t + i, i, prec);

    if (acb_equal_si(c, -1))
    {
        for (i = 1; i <= n; i += 2)
            acb_neg(t + i, t + i);
    }
    else if (!acb_is_one(c))
    {
        acb_set(d, c);

        for (i = 1; i <= n; i++)
        {
            acb_mul(t + i, t + i, d, prec);
            acb_mul(d, d, c, prec);
        }
    }

    _acb_poly_mullow(u, p, len, t, len, len, prec);

    arb_mul(f, f, f, prec);

    if (arb_bits(f) > 0.25 * prec)
    {
        arb_inv(f, f, prec);
    }
    else
    {
        for (i = 0; i <= n; i++)
            acb_div_arb(u + i, u + i, f, prec);

        arb_one(f);
    }

    for (i = n; i >= 0; i--)
    {
        acb_mul_arb(p + i, u + n - i, f, prec);
        arb_mul_ui(f, f, (i == 0) ? 1 : i, prec);
    }

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);

    arb_clear(f);
    acb_clear(d);
}

void
acb_poly_taylor_shift_convolution(acb_poly_t g, const acb_poly_t f,
        const acb_t c, slong prec)
{
    if (f != g)
        acb_poly_set_round(g, f, prec);

    _acb_poly_taylor_shift_convolution(g->coeffs, c, g->length, prec);
}
#endif

