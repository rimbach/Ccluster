#!/bin/bash

usage()
{
   echo "Usage: ./tableParallel <options> <args>"
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
      --bglobal)
        BGLOBAL="$VALUE"
        ;;
      --nbthreads)
        NBTHREADS="$VALUE"
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
   DEGREES="128 256"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$BGLOBAL" ]; then
   BGLOBAL="global"
fi

if [ -z "$NBTHREADS" ]; then
   NBTHREADS="2 4"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

STRATFLAG="V5"

VERBOFLAG=2

REP="tableParallel"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

##########################solvers
CCLUSTER_PATH="../"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
GENPOLFILE_CALL=$CCLUSTER_PATH"/bin/genPolFile"

# MPSOLVE_CALL="mpsolve -au -Gi -o"$EPSILONMPS" -j1"
# MPSOLVE_CALL_S="mpsolve -as -Gi -o"$EPSILONMPS" -j1"

##########################naming
POL_NAME="Bernoulli"

# 
# HEAD_TABLE="\begin{tabular}{l||c|c|c||c|c|c||c|c|}"
# FIRST_LINE="     & \multicolumn{3}{|c||}{\texttt{Ccluster} local (\$[-1,1]^2\$)}"
# FIRST_LINE=$FIRST_LINE"& \multicolumn{3}{|c||}{\texttt{Ccluster} global (\$[-75,75]^2\$)}"
# FIRST_LINE=$FIRST_LINE"& \texttt{unisolve} & \texttt{secsolve} \\\\\\hline"
# SECOND_LINE="d& (\#Sols:\#Clus) & (depth:size) & \$\tau_\ell\$ (s)"
# SECOND_LINE=$SECOND_LINE"& (\#Sols:\#Clus) & (depth:size) & \$\tau_g\$ (s)"
# SECOND_LINE=$SECOND_LINE"& \$\tau_u\$ (s) & \$\tau_s\$ (s)   \\\\\\hline\hline"
# TAIL_TAB="\end{tabular}"

TEMPTABFILE="temptabfileParall.txt"
touch $TEMPTABFILE

for DEG in $DEGREES; do
    LINE_TAB=$DEG
    
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    $GENPOLFILE_CALL $POL_NAME $DEG $FILEIN -f 1
    
    echo  "Clustering roots for $POL_NAME, degree $DEG, global, 1 threads"
    (/usr/bin/time -p $CCLUSTER_CALL $FILEIN -e $EPSILONCCL -m $STRATFLAG -v $VERBOFLAG > tempparall.txt ) &>> tempparall.txt
#     cat tempparall.txt
    TIMEREF=$(grep "real" tempparall.txt| cut -f2 -d' '| tr -d ' ')
    LINE_TAB=$LINE_TAB" & "`format_time $TIMEREF`
        
    for NBT in $NBTHREADS; do
    
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, $NBT threads"
        (/usr/bin/time -p $CCLUSTER_CALL $FILEIN -e $EPSILONCCL -m $STRATFLAG -v $VERBOFLAG -j $NBT > tempparall.txt ) &>> tempparall.txt
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
