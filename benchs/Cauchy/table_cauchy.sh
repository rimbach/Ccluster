#!/bin/bash

usage()
{
   echo "Usage: ./table_cauchy.sh <options> <args>"
   echo "where options are:"
   echo "--purge: erase all files in rep"
   echo "--purgeCCL: erase ccluster result files"
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
      --epsilonCCL)
        EPSILONCCL=$VALUE
        ;;
      --epsilonCAU)
        EPSILONCAU=$VALUE
        ;;
      --epsilonMPL)
        EPSILONMPL=$VALUE
        ;;
      --purge)
        PURGE=1
        ;;
      --purgeCCL)
        PURGECCL=1
        ;;
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

if [ -z "$PURGECCL" ]; then
   PURGECCL=0
fi

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
#    NBITT="7 8 9"
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

NBCOLS=11

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

LIMDEGCCL=512
LIMINDCCL=9

POLNAME="randomSparse"
# DEGREES="64 128 191 256 391 512 717 1024 1523 2048"
DEGREES="512 717 1024 1523 2048 3051 4096 6099"
NBTERMS="5 10 20"
# BITSIZES="16"
BITSIZES="256"
NBPOLS="2"
for NBT in $NBTERMS; do
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials with $NBT terms} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
for BIT in $BITSIZES; do
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT"_"$NBT
    
    genRand_with_deg_bs_nbterms $NAME $POLNAME $DEG $BIT $NBT $NBPOLS $REPNAME
    for CURIND in `seq 1 $NBPOLS`; do
        NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT"_"$NBT"_"$CURIND
        run_ccluster $NAME $POLNAME $DEG $EPSILONCCL
        run_cauchy   $NAME $POLNAME $DEG $EPSILONCAU
        run_mpsolve  $NAME $POLNAME $DEG $EPSILONMPL
    done
    
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT"_"$NBT
    stats_pol_rand $NAME $DEG $NBPOLS
    
done
done
done

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
    
    gen_with_deg $NAME $POLNAME $DEG
    run_ccluster $NAME $POLNAME $DEG $EPSILONCCL
    run_cauchy   $NAME $POLNAME $DEG $EPSILONCAU
    run_mpsolve  $NAME $POLNAME $DEG $EPSILONMPL
    
    stats_pol $NAME $DEG
done 
done
echo "\\hline" >> $TEMPTABFILE

# DEGREES="64 128 191 256 391 512 717 1024 1523 2048 3051 4096"
DEGREES="512 717 1024 1523 2048 3051 4096 6099 8192"
BITS="14"
POLNAMES="Mignotte"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials, \$a=$BITS\$} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS
    
    gen_with_deg_bs         $NAME $POLNAME $DEG $BITS
    run_ccluster            $NAME $POLNAME $DEG       $EPSILONCCL
#     run_cauchy_mignotte     $NAME $POLNAME $DEG $BITS $EPSILONCAU
    run_cauchy              $NAME $POLNAME $DEG $EPSILONCAU
    run_mpsolve             $NAME $POLNAME $DEG       $EPSILONMPL
    
    stats_pol $NAME $DEG
done 
done
echo "\\hline" >> $TEMPTABFILE

DEGREES="5 6 7 8 9 10 11"
POLNAMES="Mandelbrot"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg            $NAME $POLNAME $DEG
    run_ccluster_with_ind   $NAME $POLNAME $DEG $EPSILONCCL
    run_cauchy_mandelbrot   $NAME $POLNAME $DEG $EPSILONCAU
    run_mpsolve             $NAME $POLNAME $DEG $EPSILONMPL
    
    stats_pol $NAME $DEG
done 
done
echo "\\hline" >> $TEMPTABFILE

DEGREES="5 6 7 8 9 10 11 12"
POLNAMES="Runnels"

for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
    echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials} \\\\\\hline" >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg            $NAME $POLNAME $DEG
    run_ccluster_with_ind   $NAME $POLNAME $DEG $EPSILONCCL
    run_cauchy_runnels      $NAME $POLNAME $DEG $EPSILONCAU
    run_mpsolve             $NAME $POLNAME $DEG $EPSILONMPL
    
    stats_pol $NAME $DEG
done 
done

cat $TEMPTABFILE
rm -f $TEMPTABFILE

# DEGREES="256"
# BITS="14"
# EPSILONS="-1000 -10000"
# POLNAMES="Mignotte"
#     POLNAME="Mignotte"
# echo "\\multicolumn{$NBCOLS}{c}{$POLNAME polynomials, \$a=$BITS\$} \\\\\\hline" >> $TEMPTABFILE
# # for EPS in $EPSILONS; do
#     DEG="256"
#     EPS="-2048"
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS"_"$EPS
#     
#     gen_with_deg_bs         $NAME $POLNAME $DEG $BITS
#     run_ccluster            $NAME $POLNAME $DEG       $EPS
#     run_cauchy_mignotte     $NAME $POLNAME $DEG $BITS $EPS
#     
#     stats_pol_mig $NAME $DEG $EPS
#     
#     DEG="391"
#     EPS="-4096"
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS"_"$EPS
#     
#     gen_with_deg_bs         $NAME $POLNAME $DEG $BITS
#     run_ccluster            $NAME $POLNAME $DEG       $EPS
#     run_cauchy_mignotte     $NAME $POLNAME $DEG $BITS $EPS
#     
#     stats_pol_mig $NAME $DEG $EPS
#     
#     DEG="512"
#     EPS="-4096"
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITS"_"$EPS
#     
#     gen_with_deg_bs         $NAME $POLNAME $DEG $BITS
#     run_ccluster            $NAME $POLNAME $DEG       $EPS
#     run_cauchy_mignotte     $NAME $POLNAME $DEG $BITS $EPS
#     
#     stats_pol_mig $NAME $DEG $EPS
# # done
# echo "\\hline" >> $TEMPTABFILE
# 
# cat $TEMPTABFILE
# rm -f $TEMPTABFILE
