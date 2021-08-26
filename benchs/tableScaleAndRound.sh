#!/bin/bash

ECHO_ME="/bin/echo -e"
BGREEN="\e[1;32m"
BRED="\e[1;31m"
NORMAL="\e[0m"

# #tables with Bernoulli pol
# DEGREES="128 191"
# DEGREESRISO="383 512"
# $ECHO_ME $BGREEN "Bernoulli Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableBernoulli.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableBernoulli.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableBernoulli.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableBernoulli.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"
# # 
# # #tables with Laguerre pol
# DEGREES="128 191"
# DEGREESRISO="383 512"
# $ECHO_ME $BGREEN "Laguerre Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableLaguerre.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableLaguerre.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableLaguerre.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableLaguerre.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"
# 
# #tables with Mandelbrot pol
# DEGREES="6 7 8"
# DEGREESRISO="8 9 10"
# $ECHO_ME $BGREEN "Mandelbrot Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableMandelbrot.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableMandelbrot.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableMandelbrot.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableMandelbrot.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"
# 
# #table with Mignotte pol
# DEGREES="128"
# BITSIZES="32 256"
# DEGREESRISO="128"
# BITSIZESRISO="32 256"
# EPSCCL="-53"
# EPSMPS="16"
# $ECHO_ME $BGREEN "Mignotte Polynomials" $NORMAL
# $ECHO_ME $BGREEN "    No root radii, No scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO"  --mflag="default" --rep="tableScaleAndRound"
# 
# EPSCCL="-2120"
# EPSMPS="640"
# $ECHO_ME $BGREEN "Mignotte Polynomials" $NORMAL
# $ECHO_ME $BGREEN "    No root radii, No scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="V7" --rep="tableScaleAndRound" --epsilonMPS="$EPSMPS" --epsilonCCL="$EPSCCL"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="onlySubd" --rep="tableScaleAndRound" --epsilonMPS="$EPSMPS" --epsilonCCL="$EPSCCL"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO" --mflag="V8" --rep="tableScaleAndRound" --epsilonMPS="$EPSMPS" --epsilonCCL="$EPSCCL"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableMignotte.sh --purge --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --bitsizes="$BITSIZES" --bitsizesRiso="$BITSIZESRISO"  --mflag="default" --rep="tableScaleAndRound" --epsilonMPS="$EPSMPS" --epsilonCCL="$EPSCCL"

# tables with Nested Clusters pol
# ITERATIONS="3 4"
# RELWIDTHS="16 256"
# NBSOLS="3"
# ITERATIONSRISO="3 4"
# RELWIDTHSRISO="16 256"
# NBSOLSRISO="3"
# $ECHO_ME $BGREEN "Nested Clusters Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableNestedClusters.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --relWidth="$RELWIDTHS" --nbSols="$NBSOLS" --relWidthRiso="$RELWIDTHSRISO" --nbSolsRiso="$NBSOLSRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableNestedClusters.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --relWidth="$RELWIDTHS" --nbSols="$NBSOLS" --relWidthRiso="$RELWIDTHSRISO" --nbSolsRiso="$NBSOLSRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableNestedClusters.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --relWidth="$RELWIDTHS" --nbSols="$NBSOLS" --relWidthRiso="$RELWIDTHSRISO" --nbSolsRiso="$NBSOLSRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableNestedClusters.sh --iterations="$DEGREES" --iterationsRiso="$DEGREESRISO" --relWidth="$RELWIDTHS" --nbSols="$NBSOLS" --relWidthRiso="$RELWIDTHSRISO" --nbSolsRiso="$NBSOLSRISO"  --mflag="default" --rep="tableScaleAndRound"

# #tables with Wilkinson pol
# DEGREES="128 191"
# DEGREESRISO="383 512"
# $ECHO_ME $BGREEN "Wilkinson Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableWilkinson.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableWilkinson.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableWilkinson.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableWilkinson.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"

#tables with Wilkinson rat
# DEGREES="128 191"
# DEGREESRISO="383 512"
# $ECHO_ME $BGREEN "WilkRat Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableWilkRat.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableWilkRat.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableWilkRat.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableWilkRat.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"

# #tables with Wilkinson Mul
# DEGREES="10 15 19"
# DEGREESRISO="5"
# $ECHO_ME $BGREEN "WilkMul Polynomials" $NORMAL
# $ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
# ./tableWilkMul.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
# ./tableWilkMul.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
# ./tableWilkMul.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
# $ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
# ./tableWilkMul.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"

# #tables with Wilkinson Clus
DEGREES="12 16 20 24 28 32"
DEGREESRISO="12 16 20 24 28 32"
$ECHO_ME $BGREEN "WilkClus Polynomials" $NORMAL
$ECHO_ME $BGREEN "                     No root radii, No scale and round" $NORMAL
./tableWilkClus.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V7" --rep="tableScaleAndRound"
$ECHO_ME $BGREEN "    No root radii, scale and round" $NORMAL
./tableWilkClus.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="onlySubd" --rep="tableScaleAndRound"
$ECHO_ME $BGREEN "    Root radii,    No scale and round" $NORMAL
./tableWilkClus.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO" --mflag="V8" --rep="tableScaleAndRound"
$ECHO_ME $BGREEN "    Root radii,    scale and round" $NORMAL
./tableWilkClus.sh --degrees="$DEGREES" --degreesRiso="$DEGREESRISO"  --mflag="default" --rep="tableScaleAndRound"
