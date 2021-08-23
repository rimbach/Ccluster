#!/bin/bash

usage()
{
   echo "Usage: ./tableParallel <options> <args>"
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
      --epsilonCCL)
        EPSILONCCL=$VALUE
        ;;
      --nbthreads)
        NBTHREADS="$VALUE"
        ;;
      --purge)
        PURGE=1
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
   DEGREES="128 256"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$NBTHREADS" ]; then
   NBTHREADS="2 4"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

if [ -z "$MFLAG" ]; then
   MFLAG="default"
fi

if [ -z "$REP" ]; then
   REP="tableParallel"
fi

##########################solvers
CCLUSTER_PATH="../"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
GENPOLFILE_CALL=$CCLUSTER_PATH"/bin/genPolFile"

POLNAME="Bernoulli"

source functions.sh

init_rep $REP

TEMPTABFILE="temptabfileParall.txt"
touch $TEMPTABFILE

for DEG in $DEGREES; do
    LINE_TAB=$DEG

    FILENAME=$REP"/"$POLNAME"_"$DEG
    gen_with_deg $FILENAME $POLNAME $DEG
    
    $ECHO  "Clustering roots for $POLNAME, degree $DEG, global, 1 threads"
    (/usr/bin/time -p $CCLUSTER_CALL $FILENAME".ccl" -e $EPSILONCCL -m "onlySubd" -v 2 > tempparall.txt ) &>> tempparall.txt
#     cat tempparall.txt
    TIMEREF=$(grep "real" tempparall.txt| cut -f2 -d' '| tr -d ' ')
    LINE_TAB=$LINE_TAB" & "`format_time $TIMEREF`
        
    for NBT in $NBTHREADS; do
    
        $ECHO  "Clustering roots for $POLNAME, degree $DEG, global, $NBT threads"
        (/usr/bin/time -p $CCLUSTER_CALL $FILENAME".ccl" -e $EPSILONCCL -m "onlySubd" -v 2 -j $NBT > tempparall.txt ) &>> tempparall.txt
        TIME=$(grep "real" tempparall.txt| cut -f2 -d' '| tr -d ' ')
        LINE_TAB=$LINE_TAB" & "`ratio_time $TIMEREF $TIME`
#         echo `ratio_time $TIMEREF $TIME`
    done
    LINE_TAB=$LINE_TAB" \\\\"
    
    echo $LINE_TAB >> $TEMPTABFILE
done

# echo $HEAD_TABLE
# echo $FIRST_LINE
# echo $SECOND_LINE
cat $TEMPTABFILE
# echo $TAIL_TAB

rm -rf $TEMPTABFILE
rm -rf tempparall.txt
