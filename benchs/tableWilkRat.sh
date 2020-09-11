#!/bin/bash

usage()
{
   echo "Usage: ./tableWilkRat <options> <args>"
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
   DEGREES="64 128"
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
   BGLOBAL="global"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

VFLAG="V5"

REP="tableWilkRat"

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

MPSOLVE_CALL="mpsolve -au -Gi -o"$EPSILONMPS" -j1"
MPSOLVE_CALL_S="mpsolve -as -Gi -o"$EPSILONMPS" -j1"

##########################naming
POL_NAME="WilkRat"

HEAD_TABLE="\begin{tabular}{l||c|c|c||c|c|c||c|c|}"
FIRST_LINE="     & \multicolumn{3}{|c||}{\texttt{Ccluster} local (\$[-1,1]^2\$)}"
FIRST_LINE=$FIRST_LINE"& \multicolumn{3}{|c||}{\texttt{Ccluster} global (\$[-75,75]^2\$)}"
FIRST_LINE=$FIRST_LINE"& \texttt{unisolve} & \texttt{secsolve} \\\\\\hline"
SECOND_LINE="d& (\#Sols:\#Clus) & (depth:size) & \$\tau_\ell\$ (s)"
SECOND_LINE=$SECOND_LINE"& (\#Sols:\#Clus) & (depth:size) & \$\tau_g\$ (s)"
SECOND_LINE=$SECOND_LINE"& \$\tau_u\$ (s) & \$\tau_s\$ (s)   \\\\\\hline\hline"
TAIL_TAB="\end{tabular}"

TEMPTABFILE="temptabfileWilk.txt"
touch $TEMPTABFILE

for DEG in $DEGREES; do
    LINE_TAB=$DEG
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    FILENAME_CCLUSTER_L_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local.out"
    FILENAME_CCLUSTER_G_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_U_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_u.out"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    $GENPOLFILE_CALL $POL_NAME $DEG $FILEIN -f 1
    $GENPOLFILE_CALL $POL_NAME $DEG $FILENAME_MPSOLVE_IN -f 2
    
    if [ ! -e $FILENAME_CCLUSTER_L_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_OUT > /dev/stderr
        $CCLUSTER_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAG -v 3 > $FILENAME_CCLUSTER_L_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_OUT > /dev/stderr
        $CCLUSTER_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAG -v 3 > $FILENAME_CCLUSTER_G_OUT
    fi
     
    if [ ! -e $FILENAME_MPSOLVE_U_OUT ]; then
        if [ $(( $DEG < 300 ? 0 : 1 )) -eq 0 ]; then
            echo "Isolating roots in C with MPSOLVE UNISOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_U_OUT) &>> $FILENAME_MPSOLVE_U_OUT
        else
            echo "real    \$>1000\$" > $FILENAME_MPSOLVE_U_OUT
        fi
    fi
    if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
        echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
        (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
    fi
    TCCLUSTER_L=$(grep "time:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G=$(grep "time:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G=$(grep "tree depth:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G=$(grep "tree size:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G=$(grep "number of clusters:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_U_OUT| cut -f2 -d'l' | tr -d ' ')
    TMPSOLVE_S=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&("$NBSOLUTION_L":"$NBCLUSTERS_L")&("$TDCCLUSTER_L":"$TSCCLUSTER_L")&`format_time $TCCLUSTER_L`"
    LINE_TAB=$LINE_TAB"&("$NBSOLUTION_G":"$NBCLUSTERS_G")&("$TDCCLUSTER_G":"$TSCCLUSTER_G")&`format_time $TCCLUSTER_G`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE"&"$TMPSOLVE_S
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
done

echo $HEAD_TABLE
echo $FIRST_LINE
echo $SECOND_LINE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE
