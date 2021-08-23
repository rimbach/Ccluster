#!/bin/bash

POLNAME="Mignotte"

usage()
{
   echo "Usage: ./table$POLNAME <options> <args>"
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
      --verbose)
        VERBOSE=1
        ;;
      --degrees)
        DEGREES=$VALUE
        ;;
      --bitsizes)
        BITSIZES=$VALUE
        ;;
      --degreesRiso)
        DEGREESRISO=$VALUE
        ;;
      --bitsizesRiso)
        BITSIZESRISO=$VALUE
        ;;
      --epsilonCCL)
        EPSILONCCL=$VALUE
        ;;
      --epsilonMPS)
        EPSILONMPS=$VALUE
        ;;
      --blocal)
        BLOCAL="$VALUE"
        ;;
      --blocalRiso)
        BLOCALRISO="$VALUE"
        ;;
      --purge)
        PURGE=1
        ;;
      --purgeMPSOLVE)
        PURGEMPSOLVE=1
        ;;
      --purgeDSC)
        PURGEDSC=1
        ;;
      --purgeCCLLOCAL)
        PURGECCLLOCAL=1
        ;;
      --purgeCCLGLOBAL)
        PURGECCLGLOBAL=1
        ;;
      --purgeRISOLOCAL)
        PURGECCLLOCAL=1
        ;;
      --purgeRISOGLOBAL)
        PURGERISOGLOBAL=1
        ;;
      --mflag)
        MFLAG="$VALUE"
        ;;
      --rep)
        REP="$VALUE"
        ;;
      *)
        usage
        exit 1
        ;;
    esac
    shift
done

#default values
if [ -z "$VERBOSE" ]; then
   VERBOSE=0
fi
ECHO=""
if [ $VERBOSE -eq 0 ]; then
     ECHO="echo -e \c "
else
     ECHO="echo "
fi

if [ -z "$DEGREES" ]; then
   DEGREES="64 128"
fi

if [ -z "$BITSIZES" ]; then
   BITSIZES="8"
fi

if [ -z "$DEGREESRISO" ]; then
   DEGREESRISO="64 128"
fi

if [ -z "$BITSIZESRISO" ]; then
   BITSIZESRISO="8"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$EPSILONMPS" ]; then
   EPSILONMPS="16"
fi

if [ -z "$BLOCAL" ]; then
   BLOCAL="0/1,0/1,2/1"
fi

if [ -z "$BLOCALRISO" ]; then
   BLOCALRISO="0/1,2/1"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

if [ -z "$PURGEMPSOLVE" ]; then
   PURGEMPSOLVE=0
fi

if [ -z "$PURGEDSC" ]; then
   PURGEDSC=0
fi

if [ -z "$PURGECCLLOCAL" ]; then
   PURGECCLLOCAL=0
fi

if [ -z "$PURGECCLGLOBAL" ]; then
   PURGECCLGLOBAL=0
fi

if [ -z "$PURGERISOLOCAL" ]; then
   PURGERISOLOCAL=0
fi

if [ -z "$PURGERISOGLOBAL" ]; then
   PURGERISOGLOBAL=0
fi

if [ -z "$MFLAG" ]; then
   MFLAG="default"
fi

if [ -z "$REP" ]; then
   REP="table$POLNAME"
fi

##########################solvers
CCLUSTER_PATH="../"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
RISOLATE_CALL=$CCLUSTER_PATH"/bin/risolate"
GENPOLFILE_CALL=$CCLUSTER_PATH"/bin/genPolFile"

MPSOLVE_CALL_S="mpsolve -as -Ga -o"$EPSILONMPS" -j1"

ANEWDSC_PATH="/work/softs"
ANEWDSC_CALL=$ANEWDSC_PATH"/test_descartes_linux64"

source functions.sh

init_rep $REP

TEMPTABFILE="temptabfileMign.txt"
touch $TEMPTABFILE

for DEG in $DEGREES; do
    for BIT in $BITSIZES; do
        FILENAME=$REP"/"$POLNAME"_"$DEG"_"$BIT
        gen_with_deg_bs $FILENAME $POLNAME $DEG $BIT
        run_ccluster_local_global $FILENAME $POLNAME $DEG
        run_mpsolve $FILENAME $POLNAME $DEG 
        stats_pol_ccluster_l_g_mpsolve $FILENAME $DEG
    done
done

echo $HEAD_TABLE
echo $FIRST_LINE
echo $SECOND_LINE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE

TEMPTABFILE="temptabfileMign.txt"
touch $TEMPTABFILE

for DEG in $DEGREESRISO; do
    for BIT in $BITSIZESRISO; do
        FILENAME=$REP"/"$POLNAME"_"$DEG"_"$BIT
        gen_with_deg_bs $FILENAME $POLNAME $DEG $BIT
        run_risolate_local_global $FILENAME $POLNAME $DEG
        run_aNewDsc $FILENAME $POLNAME $DEG "0"
        stats_pol_risolate_l_g_anewdsc $FILENAME $DEG
    done
done

echo $HEAD_TABLE_RISO
echo $FIRST_LINE_RISO
echo $SECOND_LINE_RISO
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE
