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

const pkgdir = realpath(joinpath(dirname(@__FILE__), "."))
const libcclusterdir = joinpath(pkgdir,"..")
const libccluster = joinpath(libcclusterdir, "libccluster")

push!(Libdl.DL_LOAD_PATH, libccluster)
# push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
if is_linux() 
	Libdl.dlopen(libccluster)
end

using Nemo
using PyPlot
using PyCall

include("box.jl")
include("connComp.jl")
include("listConnComp.jl")
include("plotCcluster.jl")

   
function Ccluster( getApprox::Function, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int = 0 )
    
    initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    lccRes = listConnComp()
    
    ccall( (:ccluster_interface_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &lccRes,           getApp_c,   &initBox, &eps,      strat, verbose )
     
    queueResults = []
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),[getCenterRe(tempBO),getCenterIm(tempBO),fmpq(3,4)*getWidth(tempBO)]])
    end
    
    return queueResults
    
end

function Ccluster( getApprox::Function, initBox::box, eps::fmpq, strat::Int, verbose::Int = 0 )
    
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
    const getApp_c = cfunction(getApprox, Void, (Ptr{acb_poly}, Int))
    
    lccRes = listConnComp()
    ccall( (:ccluster_interface_forJulia, :libccluster), 
             Void, (Ptr{listConnComp}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
                    &lccRes,           getApp_c,   &initBox, &eps,      strat, verbose )
     
    queueResults = []
    while !isEmpty(lccRes)
        tempCC = pop(lccRes)
        tempBO = getComponentBox(tempCC,initBox)
        push!(queueResults, [getNbSols(tempCC),tempBO])
    end
    
    return queueResults
    
end

# POLY_GLOBAL = fmpq_poly()
# 
# function GETAPP_GLOBAL( dest::Ptr{acb_poly}, prec::Int )
#     ccall((:acb_poly_set_fmpq_poly, :libarb), Void,
#                 (Ptr{acb_poly}, Ptr{fmpq_poly}, Int), dest, &POLY_GLOBAL, prec)
# end
# 
# function Ccluster( poly::fmpq_poly, initialBox::Array{fmpq,1}, eps::fmpq, strat::Int, verbose::Int = 0 )
#     
#     POLY_GLOBAL = poly
#     
#     initBox::box = box(initialBox[1],initialBox[2],initialBox[3])
#     const getApp_c = cfunction(GETAPP_GLOBAL, Void, (Ptr{acb_poly}, Int))
#     
#     lccRes = listConnComp()
#     
#     ccall( (:ccluster_interface_forJulia, :libccluster), 
#              Void, (Ptr{listConnComp}, Ptr{Void},    Ptr{box}, Ptr{fmpq}, Int,   Int), 
#                     &lccRes,           getApp_c,   &initBox, &eps,      strat, verbose )
#      
#      
#     queueResults::Array{Array{fmpq,1},1} = []
#     while !isEmpty(lccRes)
#         temp = getComponentBox(pop(lccRes),initBox)
#         push!(queueResults, [getCenterRe(temp),getCenterIm(temp),fmpq(3,4)*getWidth(temp)])
#     end
#     
#     return queueResults
#     
# end 
