#!/bin/bash

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

TRUE=1
FALSE=0
##########################getting arguments
DEGREES=$1
BITSIZES=$2
EPSILONCCL=$3
EPSILONMPS=$4
CCLUSTER=$5
MPSOLVE=$6
MPSOLVE_S=$7

##########################solving parameters
#BINIT="0,1,0,1,100,1"
BINIT="0,1,0,1,1000000000000000000000,1"
# BINIT="0,1,0,1,2,1"
# EPSIL="1,100"
EPSIL=$EPSILONCCL
# STRAT="7" # newton tstarOpt predPrec
# STRAT_DIS="newton tstarOpt predPrec"
#STRAT="15" # newton tstarOpt predPrec stopWhenCompact
#STRAT_DIS="newton tstarOpt predPrec stopWhenCompact"
#STRAT="23" # newton tstarOpt predPrec anticip
#STRAT_DIS="newton tstarOpt predPrec anticip"
STRAT="31" # newton tstarOpt predPrec stopWhenCompact anticip
STRAT_DIS="newton tstarOpt predPrec stopWhenCompact anticip"
##########################solvers
SOLVER_PATH="./"
MPSOLVE_CALL="mpsolve -au -Gi -o"$EPSILONMPS" -j1"
MPSOLVE_CALL_S="mpsolve -as -Gi -o"$EPSILONMPS" -j1"

##########################naming
POLY_NAME="Mignotte"
CCLUSTER_CALL="./bench"$POLY_NAME
GENPOLFILE_CALL="./genPolFile "$POLY_NAME
FILE_PATH="./"$POLY_NAME
FILENAME_TABLE=$FILE_PATH"/"$POLY_NAME"_table"
echo -e "\$d\$ & nb sols & nbclus & tree depth & tree size & time & time \\\\\\\\" > $FILENAME_TABLE

# ############################################ Variables for colours
neutre='\e[0;m'
bleufonce='\e[0;34m'
vertfonce='\e[0;32m'
rouge='\e[0;31m'

############################################## Loop

echo -e $bleufonce --------$POLY_NAME polynomial-------- $neutre
echo  initial box: $BINIT
echo  epsilon:     $EPSIL
echo  strategy:    $STRAT_DIS 
echo -e $bleufonce ----------------------------------- $neutre

for DEGREE in $DEGREES; do
    for BITSIZE in $BITSIZES; do
        if [ $CCLUSTER -eq $TRUE ]; then
            echo  "Clustering roots in bInit with Ccluster..."
            FILENAME_CCLUSTER_OUT=$FILE_PATH"/"$POLY_NAME"_"$DEGREE"_"$BITSIZE"_ccluster.out"
    #         echo -e $FILENAME_CCLUSTER_OUT
            $CCLUSTER_CALL $DEGREE $BITSIZE $BINIT $EPSIL $STRAT "2" > $FILENAME_CCLUSTER_OUT
            TCCLUSTER=$(grep "time:" $FILENAME_CCLUSTER_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
            TDCCLUSTER=$(grep "tree depth:" $FILENAME_CCLUSTER_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
            TSCCLUSTER=$(grep "tree size:" $FILENAME_CCLUSTER_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
            NBCLUSTERS=$(grep "number of clusters:" $FILENAME_CCLUSTER_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
            NBSOLUTION=$(grep "number of solutions:" $FILENAME_CCLUSTER_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
            echo -e time: $vertfonce $TCCLUSTER s $neutre 
            echo -e $vertfonce $NBSOLUTION $neutre found in $vertfonce $NBCLUSTERS $neutre clusters. "\c"
            echo -e tree depth: $vertfonce $TDCCLUSTER $neutre, tree size: $vertfonce $TSCCLUSTER $neutre
            echo log can be found in file $FILENAME_CCLUSTER_OUT
            echo -e $bleufonce ----------------------------------- $neutre
        fi
        if [ $MPSOLVE -eq $TRUE ]; then
            echo "Isolating roots in C with MPSOLVE........."
            FILENAME_MPSOLVE_IN=$FILE_PATH"/"$POLY_NAME"_"$DEGREE"_"$BITSIZE".pol"
            FILENAME_MPSOLVE_OUT=$FILE_PATH"/"$POLY_NAME"_"$DEGREE"_"$BITSIZE"_mpsolve.out"
            $GENPOLFILE_CALL $DEGREE $BITSIZE > $FILENAME_MPSOLVE_IN
            (time $MPSOLVE_CALL $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_OUT) &>> $FILENAME_MPSOLVE_OUT
            TMPSOLVE=$(grep "user" $FILENAME_MPSOLVE_OUT| cut -f2 -d'r' | tr -d ' ')
            echo -e time: $vertfonce $TMPSOLVE $neutre
            echo log can be found in file $FILENAME_MPSOLVE_OUT
            echo -e $bleufonce ----------------------------------- $neutre
        fi
        if [ $MPSOLVE_S -eq $TRUE ]; then    
            echo "Isolating roots in C with MPSOLVE SECULAR........."
            FILENAME_MPSOLVE_IN=$FILE_PATH"/"$POLY_NAME"_"$DEGREE"_"$BITSIZE".pol"
            FILENAME_MPSOLVE_S_OUT=$FILE_PATH"/"$POLY_NAME"_"$DEGREE"_"$BITSIZE"_mpsolve_s.out"
            $GENPOLFILE_CALL $DEGREE $BITSIZE > $FILENAME_MPSOLVE_IN
            (time $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
#             $MPSOLVE_CALL_S -Ogf $FILENAME_MPSOLVE_IN | gnuplot
            TMPSOLVE_S=$(grep "user" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'r' | tr -d ' ')
            echo -e time: $vertfonce $TMPSOLVE_S $neutre
            echo log can be found in file $FILENAME_MPSOLVE_S_OUT
            echo -e $bleufonce ----------------------------------- $neutre
        fi
        echo -e $DEGREE  "&" $NBSOLUTION "&" $NBCLUSTERS "&" $TDCCLUSTER "&" $TSCCLUSTER "&" $TCCLUSTER "&" $TMPSOLVE "&" $TMPSOLVE_S "\\\\\\\\" >> $FILENAME_TABLE
    done
done

cat $FILENAME_TABLE
