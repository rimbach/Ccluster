#!/bin/bash

usage()
{
   echo "Usage: ./table_dense.sh <options> <args>"
   echo "where options are:"
   echo "--purge: erase all files in rep"
#    echo "--purgeCCL: erase ccluster result files"
   echo "--purgeCAU: erase cauchy   result files"
}

##########################getting arguments
while [ "$1" != "" ]; do
   PARAM=`echo $1 | sed 's/=.*//'`
   VALUE=`echo $1 | sed 's/[^=]*//; s/=//'`
   case "$PARAM" in
      -h|--help)
         usage
         exit 0
         ;;
#       --degrees)
#         DEGREES=$VALUE
#         ;;
#       --bitsize)
#         BITSIZE=$VALUE
#         ;;
#       --nbpols)
#         NBPOLS=$VALUE
#         ;;
#       --nbitts)
#         NBITT=$VALUE
#         ;;
#       --sizegrid)
#         SIZEGRID=$VALUE
#         ;;
#       --nbterms)
#         NBTERMS=$VALUE
#         ;;
#       --epsilonCCL)
#         EPSILONCCL=$VALUE
#         ;;
      --epsilonCAU)
        EPSILONCAU=$VALUE
        ;;
      --epsilonMPL)
        EPSILONMPL=$VALUE
        ;;
      --purge)
        PURGE=1
        ;;
#       --purgeCCL)
#         PURGECCL=1
#         ;;
      --purgeCAU)
        PURGECAU=1
        ;;
      --purgeMPL)
        PURGEMPL=1
        ;;
#       --generate)
#         GENERATE=1
#         ;;
      *)
        usage
        exit 1
        ;;
    esac
    shift
done

#default values
if [ -z "$DEGREES" ]; then
   DEGREES="64 128 191"
#    DEGREES="128 191 256 391 512"
   LENP=5
fi

# if [ -z "$BITSIZE" ]; then
#    BITSIZE="8"
# fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$EPSILONCAU" ]; then
   EPSILONCAU="-53"
fi

if [ -z "$EPSILONMPL" ]; then
   EPSILONMPL="16"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

# if [ -z "$PURGECCL" ]; then
#    PURGECCL=0
# fi

if [ -z "$PURGECAU" ]; then
   PURGECAU=0
fi

if [ -z "$PURGEMPL" ]; then
   PURGEMPL=0
fi
# if [ -z "$GENERATE" ]; then
#    GENERATE=0
# fi
# 
# if [ -z "$NBPOLS" ]; then
#    NBPOLS=10
# fi
# 
# if [ -z "$NBITT" ]; then
#    NBITT="8 9 10 11"
# fi
# 
# if [ -z "$SIZEGRID" ]; then
#    SIZEGRID="6 8 10 12 14"
# #    SIZEGRID="5 6 7 8 10 11 12 13 14"
# #    SIZEGRID="10 11 12 13 14"
# fi
# 
# if [ -z "$NBTERMS" ]; then
#    NBTERMS=10
# fi

CCLUSTER_PATH="../../"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
CAUCHY_CALL=$CCLUSTER_PATH"/bin/cauchy"
CAUCHY_MANDELBROT_CALL=$CCLUSTER_PATH"/bin/cauchy_mandelbrot"
CAUCHY_RUNNELS_CALL=$CCLUSTER_PATH"/bin/cauchy_runnels"
CAUCHY_MIGNOTTE_CALL=$CCLUSTER_PATH"/bin/cauchy_mignotte"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
GENRANDPOLFI_CALL=$CCLUSTER_PATH"/bin/genRandPolFile"
CCLUSTER_OPTS="-v 2 -m onlySubd"
CAUCHY_OPTS="-v 2"

MPSOLVE_OPTS="-as -Ga -j1"
MPSOLVE_CALL_S="../../../softs/MPSolve/src/mpsolve/mpsolve"

TEMPTABFILE="temptab_cauchy.txt"
touch $TEMPTABFILE

NBCOLS=10

source ./functions_cauchy.sh

REP="tab_ISSAC2022"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

LIMDEGCAU_WO=3050

DEGREES="64 128 191 256 391"
# DEGREES="128 191"
# POLNAMES="Bernoulli Wilkinson"
POLNAMES="Bernoulli"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg      $NAME $POLNAME $DEG
    run_cauchy_comp   $NAME $POLNAME $DEG $EPSILONCAU
    run_mpsolve       $NAME $POLNAME $DEG $EPSILONMPL
    
    stats_pol_procedural $NAME $DEG
done 
done
echo "\\hline" >> $TEMPTABFILE

cat $TEMPTABFILE
rm -f $TEMPTABFILE
