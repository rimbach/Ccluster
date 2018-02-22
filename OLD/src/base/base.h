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

#ifndef BASE_H
#define BASE_H

#define CCLUSTER_DEFAULT_PREC 53

#define CCLUSTER_MIN(A,B) (A<=B? A : B)
#define CCLUSTER_MAX(A,B) (A>=B? A : B)

#define CCLUSTER_EXP_NUM_T0(A) ( ((&(A)->strat)->_additionalFlags)&(0x1<<0) )
#define CCLUSTER_EXP_NUM_T1(A) ( ((&(A)->strat)->_additionalFlags)&(0x1<<1) )
#define CCLUSTER_INC_TEST(A) ( ((&(A)->strat)->_additionalFlags)&(0x1<<2) )

#endif