#!/bin/bash

usage()
{
   echo "Usage: ./tableRealCoeffs <options> <args>"
}

format_time()
{
    TIME1=$1
    TIME2=$1
    TIME1=`echo $TIME1 | cut -f1 -d'.'`
    TIME2=`echo $TIME2 | cut -f2 -d'.'`
    STIME1=${#TIME1}
    STIME1=$(( 3 - $STIME1 ))
    STIME2=$(( 3 < $STIME1 ? 3 : $STIME1 ))
    STIME2=$(( 0 > $STIME2 ? 0 : $STIME2 ))
    if [ $STIME2 -eq 0 ]; then
        echo $TIME1
    else
        TIME2=`echo $TIME2 | cut -c-$( echo $STIME2)`
        echo $TIME1"."$TIME2
    fi
} 

ratio_time()
{
    NUM=$1
    DEN=$2
    RATIO=`echo $NUM/$DEN|bc -l`
    echo `format_time $RATIO`
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
      --degrees)
        DEGREES=$VALUE
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
      --bglobal)
        BGLOBAL="$VALUE"
        ;;
      --mpsolve)
        MPSOLVE="$VALUE"
        ;;
      --stop-when-compact)
        STOPWHENCOMPACT=1
        ;;
      --anticipate)
        ANTICIPATE=1
        ;;
      --purge)
        PURGE=1
        ;;
      *)
        usage
        exit 1
        ;;
    esac
    shift
done

#default values
if [ -z "$DEGREES" ]; then
   DEGREES="64 128 256"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$EPSILONMPS" ]; then
   EPSILONMPS="16"
fi

if [ -z "$BLOCAL" ]; then
   BLOCAL="0,1,0,1,2,1"
fi

if [ -z "$BGLOBAL" ]; then
   BGLOBAL="0,1,0,1,300,1"
fi

if [ -z "$MPSOLVE" ]; then
   MPSOLVE=1
fi

if [ -z "$STOPWHENCOMPACT" ]; then
   STOPWHENCOMPACT=0
fi

if [ -z "$ANTICIPATE" ]; then
   ANTICIPATE=0
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

VNFLAG=23
VRFLAG=55

REP="tableRealCoeffs"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

make_line()
{
    LINE=""
    FILE1=$1
    FILE2=$2
#     N3V3=$(grep -m 1 "total number GR:" $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     N3V4=$(grep -m 1 "total number GR:" $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TVN=$(grep "total time:"                  $FILE1 | cut -f2 -d':' | cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TVR=$(grep "total time:"                  $FILE2 | cut -f2 -d':' | cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTERN=$(grep "tree depth:"          $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTERN=$(grep "tree size:"           $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERSN=$(grep "number of clusters:"  $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTIONN=$(grep "number of solutions:" $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTERR=$(grep "tree depth:"          $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTERR=$(grep "tree size:"           $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERSR=$(grep "number of clusters:"  $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTIONR=$(grep "number of solutions:" $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    LINE=$LINE"&( $TDCCLUSTERN, $TSCCLUSTERN )"
    LINE=$LINE"&( $NBCLUSTERSN, $NBSOLUTIONN )"
    LINE=$LINE"& `format_time $TVN`"
    LINE=$LINE"&( $TDCCLUSTERR, $TSCCLUSTERR )"
    LINE=$LINE"&( $NBCLUSTERSR, $NBSOLUTIONR )"
    LINE=$LINE"& `ratio_time $TVN $TVR`"
}

TEMPTABFILE="temptabfile.txt"
touch $TEMPTABFILE

POL_NAME="Bernoulli"
CCL_CALL="../test/bernoulli"

for DEG in $DEGREES; do
    LINE_TAB="Bernoulli, \$d=$DEG\$"  
    FILE1=$REP"/"$POL_NAME"_"$DEG".out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_real.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG output in "$FILE1 > /dev/stderr
        $CCL_CALL $DEG $BGLOBAL $EPSILONCCL $VNFLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG with real coeffs output in "$FILE2 > /dev/stderr
        $CCL_CALL $DEG $BGLOBAL $EPSILONCCL $VRFLAG "3" > $FILE2
    fi
      make_line $FILE1 $FILE2
      echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

POL_NAME="Mignotte"
CCL_CALL="../test/mignotte"

for DEG in $DEGREES; do
    BS="8"
    LINE_TAB="Mignotte, bs=$BS, \$d=$DEG\$"  
    FILE1=$REP"/"$POL_NAME"_"$DEG".out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_real.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG output in "$FILE1 > /dev/stderr
        $CCL_CALL $DEG $BS $BGLOBAL $EPSILONCCL $VNFLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG with real coeffs output in "$FILE2 > /dev/stderr
        $CCL_CALL $DEG $BS $BGLOBAL $EPSILONCCL $VRFLAG "3" > $FILE2
    fi
      make_line $FILE1 $FILE2
      echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

HEAD_TABLE="\begin{tabular}{l||c|c|c||c|c|c||}"
FIRST_LINE="           & \multicolumn{3}{|c||}{\texttt{Ccluster}}"
FIRST_LINE=$FIRST_LINE"& \multicolumn{3}{|c||}{\texttt{Ccluster} real coeffs}\\\\"
SECOND_LINE="          & (depth,size) & (\#Clus,\#Sols) & t (s)"
SECOND_LINE=$SECOND_LINE"& (depth,size) & (\#Clus,\#Sols) & ratio\\\\\\hline"
TAIL_TAB="\end{tabular}"

echo $HEAD_TABLE
echo $FIRST_LINE
echo $SECOND_LINE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE
