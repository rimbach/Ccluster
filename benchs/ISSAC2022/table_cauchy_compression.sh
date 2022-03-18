#!/bin/bash

usage()
{
   echo "Usage: ./table_cauchy_compression.sh <options> <args>"
   echo "where options are:"
   echo "--purge: erase all files in rep"
   echo "--purgeCAUWO: redo cauchy without compression "
   echo "--purgeCAUWI: redo cauchy with compression    "
   echo "--purgeMPS:   redo cauchy with compression    "
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
      --nbCorrectDigits)
        NBCORDIG=$VALUE
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
      --purgeMPS)
        PURGEMPS=1
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
   DEGREES="1024 2048"
#    DEGREES="128 191 256 391 512"
   LENP=5
fi

# if [ -z "$BITSIZE" ]; then
#    BITSIZE="8"
# fi

if [ -z "$NBCORDIG" ]; then
   NBCORDIG="5 10 50"
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

if [ -z "$PURGEMPS" ]; then
   PURGEMPS=0
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
CAUCHY_OPTS="-v 2 -f 1"

MPSOLVE_OPTS="-as -Ga -j1"
MPSOLVE_CALL_S="../../../softs/MPSolve/src/mpsolve/mpsolve"

TEMPTABFILE="temptab_cauchy.txt"
touch $TEMPTABFILE

NBCOLS=9

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

DEGREES="1024 1535 2048"
BITS="16"
POLNAMES="Mignotte"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials, \$a=$BITS\$} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    for NBC in $NBCORDIG; do
#         echo $NBC
        NBC2=`convert_nbdigits $NBC`
#         echo $NBC2
        REPNAME=$REP
        NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS
        
        gen_with_deg_bs          $NAME $POLNAME $DEG $BITS
        run_cauchy_comp          $NAME $POLNAME $DEG $NBC2
        run_mpsolve              $NAME $POLNAME $DEG $NBC
        stats_pol_comp           $NAME $DEG          $NBC $NBC2
        
    done
done 
echo "\\hline" >> $TEMPTABFILE
done
echo "\\hline" >> $TEMPTABFILE

NBITT="10 11"
POLNAMES="Mandelbrot"
MPSOLVE_MAX_IT=11

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials} \\\\\\hline" >> $TEMPTABFILE
for DEG in $NBITT; do
#     echo $DEG
    for NBC in $NBCORDIG; do
#         echo $DEG
        NBC2=`convert_nbdigits $NBC`
#         echo $DEG
        REPNAME=$REP
        NAME=$REPNAME"/"$POLNAME"_"$DEG
    
        gen_with_deg                 $NAME $POLNAME $DEG
        run_cauchy_mandelbrot_comp   $NAME $POLNAME $DEG $NBC2
#         if [ $DEG -le $MPSOLVE_MAX_IT ]; then
            run_mpsolve              $NAME $POLNAME $DEG $NBC
#         fi
#         run_mpsolve_man              $NAME $POLNAME $DEG $EPSILONMPL
    
        stats_pol_comp           $NAME $DEG          $NBC $NBC2
    done
done 
echo "\\hline" >> $TEMPTABFILE
done
echo "\\hline" >> $TEMPTABFILE

cat $TEMPTABFILE
rm -f $TEMPTABFILE
