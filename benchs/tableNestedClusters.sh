#!/bin/bash

POLNAME="nestedClusters"

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
      --iterations)
        DEGREES=$VALUE
        ;;
      --relWidth)
        RELWIDTH=$VALUE
        ;;
      --nbSols)
        NBSOLS=$VALUE
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
      --purge)
        PURGE=1
        ;;
      --purgeMPSOLVE)
        PURGEMPSOLVE=1
        ;;
      --purgeCCLLOCAL)
        PURGECCLLOCAL=1
        ;;
      --purgeCCLGLOBAL)
        PURGECCLGLOBAL=1
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
   DEGREES="2 3"
fi

if [ -z "$RELWIDTH" ]; then
   RELWIDTH="16"
fi

if [ -z "$NBSOLS" ]; then
   NBSOLS="3"
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

if [ -z "$PURGE" ]; then
   PURGE=0
fi

if [ -z "$PURGEMPSOLVE" ]; then
   PURGEMPSOLVE=0
fi

if [ -z "$PURGECCLLOCAL" ]; then
   PURGECCLLOCAL=0
fi

if [ -z "$PURGECCLGLOBAL" ]; then
   PURGECCLGLOBAL=0
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
GENPOLFILE_CALL=$CCLUSTER_PATH"/bin/genPolFile"

MPSOLVE_CALL_S="mpsolve -as -Ga -o"$EPSILONMPS" -j1"

source functions.sh

init_rep $REP
TEMPTABFILE="temptabfileNes.txt"

touch $TEMPTABFILE

for NS in $NBSOLS; do
    for RW in $RELWIDTH; do
        for DEG in $DEGREES; do
            FILENAME=$REP"/"$POLNAME"_"$NS"_"$RW"_"$DEG
            gen_with_deg_rw_ns $FILENAME $POLNAME $DEG $RW $NS
            run_ccluster_local_global $FILENAME $POLNAME $DEG
            run_mpsolve $FILENAME $POLNAME $DEG 
            stats_pol_ccluster_l_g_mpsolve $FILENAME $DEG
        done
    done
done

echo $HEAD_TABLE
echo $FIRST_LINE
echo $SECOND_LINE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE
