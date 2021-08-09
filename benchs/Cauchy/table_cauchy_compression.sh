#!/bin/bash

usage()
{
   echo "Usage: ./table_cauchy_compression.sh <options> <args>"
   echo "where options are:"
   echo "--purge: erase all files in rep"
   echo "--purgeCAUWO: erase cauchy without compression result files"
   echo "--purgeCAUWI: erase cauchy with compression    result files"
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
      --epsilonCAU)
        EPSILONCAU=$VALUE
        ;;
      --purge)
        PURGE=1
        ;;
      --purgeCAUWO)
        PURGECAUWO=1
        ;;
      --purgeCAUWI)
        PURGECAUWI=1
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

if [ -z "$EPSILONCAU" ]; then
   EPSILONCAU="-53"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

if [ -z "$PURGECAUWO" ]; then
   PURGECAUWO=0
fi

if [ -z "$PURGECAUWI" ]; then
   PURGECAUWI=0
fi

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

TEMPTABFILE="temptab_cauchy.txt"
touch $TEMPTABFILE

NBCOLS=13

source ./functions_cauchy.sh

REP="tab_cauchy"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

KS="3 4"
POLNAME="nestedClusters"
AS="16 256"
# AS="16 256 4096"
CS="3"
for C in $CS; do
for A in $AS; do
    echo "\\multicolumn{$NBCOLS}{c}{Polynomial with nested clusters, \$a=$A\$, \$c=$C\$} \\\\\\hline" >> $TEMPTABFILE
for K in $KS; do
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$C"_"$A"_"$K
    
    gen_with_c_a_k   $NAME $POLNAME $C $A $K
    run_cauchy_comp  $NAME $POLNAME $EPSILONCAU
    
    stats_pol_comp $NAME $K
done
echo "\\hline" >> $TEMPTABFILE
done
done

DEGREES="64 128 191 256"
POLNAMES="Bernoulli"

for POLNAME in $POLNAMES; do
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_cauchy_comp $NAME $POLNAME $EPSILONCAU
    
    stats_pol_comp $NAME $DEG
done 
echo "\\hline" >> $TEMPTABFILE
done

DEGREES="64 128 191"
BITS="14"
POLNAMES="Mignotte"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials, \$a=$BITS\$} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS
    
    gen_with_deg_bs          $NAME $POLNAME $DEG $BITS
    run_cauchy_mignotte_comp $NAME $POLNAME $DEG $BITS $EPSILONCAU
    
    stats_pol_comp $NAME $DEG
done 
done
echo "\\hline" >> $TEMPTABFILE

cat $TEMPTABFILE
rm -f $TEMPTABFILE
