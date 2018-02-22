#
#  Copyright (C) 2018 Remi Imbach
#
#  This file is part of Ccluster.
#
#  Ccluster is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License (LGPL) as published
#  by the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
#

using Nemo

R, x = PolynomialRing(Nemo.QQ, "x")

n = 8 #degree
P = zero(R)
Nemo.bernoulli_cache(n)
for k = 0:n
    coefficient = (Nemo.binom(n,k))*(Nemo.bernoulli(n-k))
    P = P + coefficient*x^k
end
# N = round(Int,5+ceil(log2(1+log2(degree(P)))))
# print("N: $N, 2^$N: $(2^N), N-log2(degree(P)): $(N-log2(degree(P)))\n")


function getAppBern( dest::Ptr{acb_poly}, prec::Int )
    ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
                (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &P, prec)
end

include("../julia/ccluster.jl")

bInit = [fmpq(0,1),fmpq(0,1),fmpq(150,1)]
# bInit = [fmpq(5,4),fmpq(0,1),fmpq(1,100)]
# bInit = [fmpq(0,1),fmpq(0,1),fmpq(5,1)]
# bInit = [0,0,1]
eps = fmpq(1,100)
#eps = fmpq(1,2^(53))
# eps = fmpq(2^(23),1)
#case with adventitious components: n=30, bInit = box(0,0,5), eps = Nemo.fmpq(1,10)
#case with adventitious components: n=70, bInit = box(0,0,10), eps = Nemo.fmpq(1,10)
    
Res = Ccluster(getAppBern, bInit, eps, 15, 0);
plotCcluster(Res, bInit, false)

