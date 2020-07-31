#!/bin/bash

usage()
{
   echo "Usage: ./tableMACIS <options> <args>"
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
   DEGREES="64"
   DEGREES="64 128"
#    DEGREES="64 128 191"
#    DEGREES="64 128 191 256"
#    DEGREES="64 128 191 256 383"
fi

if [ -z "$DEGMAND" ]; then
   DEGMAND="7 8 9"
fi

if [ -z "$DEGRUNN" ]; then
   DEGRUNN="8 9 10"
fi

if [ -z "$BITSIZE" ]; then
   BITSIZE="14"
#    DEGREES="64 128 191 256 383"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$EPSILONMPS" ]; then
   EPSILONMPS="16"
fi

if [ -z "$BLOCAL" ]; then
   BLOCAL="0,1,0,1,1,1"
fi

if [ -z "$BLOCALM" ]; then
   BLOCALM="-1,1,0,1,1,2"
fi

if [ -z "$BGLOBAL" ]; then
   BGLOBAL="global"
fi

if [ -z "$MPSOLVE" ]; then
   MPSOLVE=1
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

VFLAGV4="V4"
VFLAGV5="default"
VFLAGV6="V6"

REP="tableMACIS"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

##########################solvers
SOLVER_PATH="./"
MPSOLVE_CALL_S="mpsolve -as -Gi -o"$EPSILONMPS" -j1"

HEAD_TABLE="\begin{tabular}{l||c|c|c||c|c|c||c|c|}"
FIRST_LINE="     & \multicolumn{3}{|c||}{\texttt{Ccluster} local (\$[-1,1]^2\$)}"
FIRST_LINE=$FIRST_LINE"& \multicolumn{3}{|c||}{\texttt{Ccluster} global (\$[-75,75]^2\$)}"
FIRST_LINE=$FIRST_LINE"& \texttt{unisolve} & \texttt{secsolve} \\\\\\hline"
SECOND_LINE="d& (\#Sols:\#Clus) & (depth:size) & \$\tau_\ell\$ (s)"
SECOND_LINE=$SECOND_LINE"& (\#Sols:\#Clus) & (depth:size) & \$\tau_g\$ (s)"
SECOND_LINE=$SECOND_LINE"& \$\tau_u\$ (s) & \$\tau_s\$ (s)   \\\\\\hline\hline"
TAIL_TAB="\end{tabular}"

TEMPTABFILE="temptabMACIS.txt"
touch $TEMPTABFILE

TEMPTABFILE2="temptabMACIS2.txt"
touch $TEMPTABFILE2

CCL_CALL="../../bin/ccluster"
GPL_CALL="../../bin/genPolFile "

##########################naming
POL_NAME="Bernoulli"

for DEG in $DEGREES; do
#     LINE_TAB=$POL_NAME", \$d="$DEG"\$"
#     LINE_TAB2=$POL_NAME", \$d="$DEG"\$"
    LINE_TAB="\$\Ber{"$DEG"}\$"
    LINE_TAB2="\$\Ber{"$DEG"}\$"
    
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V5.out"
    FILENAME_CCLUSTER_G_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V5.out"
    FILENAME_CCLUSTER_L_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    $GPL_CALL $POL_NAME $DEG $FILEIN
    $GPL_CALL $POL_NAME $DEG $FILENAME_MPSOLVE_IN -f 2
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_L_V5_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_G_V5_OUT
    fi
    
    if [ ! -e $FILENAME_CCLUSTER_L_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V6_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_L_V6_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V6_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_G_V6_OUT
    fi

    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:"                 $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:"                 $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V5=$(grep "time:"                 $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V5=$(grep "time:"                 $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V6=$(grep "time:"                 $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V6=$(grep "time:"                 $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_V6"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_L_V4 $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_V6"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    LINE_TAB2=$LINE_TAB2"&("$NBCLUSTERS_G_V4", "$NBSOLUTION_G_V4")"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V4", "$TSCCLUSTER_G_V4")&`format_time $TCCLUSTER_G_V4`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V5", "$TSCCLUSTER_G_V5")&`format_time $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V6", "$TSCCLUSTER_G_V6")&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V5 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    echo $LINE_TAB2 >> $TEMPTABFILE2
    
done

##########################naming
POL_NAME="Mignotte"

for DEG in $DEGREES; do
#     LINE_TAB=$POL_NAME", \$a="$BITSIZE"\$, \$d="$DEG"\$"
#     LINE_TAB2=$POL_NAME", \$a="$BITSIZE"\$, \$d="$DEG"\$"
    LINE_TAB="\$\Mig{"$BITSIZE"}{"$DEG"}\$"
    LINE_TAB2="\$\Mig{"$BITSIZE"}{"$DEG"}\$"
    
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V5.out"
    FILENAME_CCLUSTER_G_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V5.out"
    FILENAME_CCLUSTER_L_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    $GPL_CALL $POL_NAME $DEG $FILEIN -b 14
    $GPL_CALL $POL_NAME $DEG $FILENAME_MPSOLVE_IN -b 14 -f 2
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_L_V5_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_G_V5_OUT
    fi
    
    if [ ! -e $FILENAME_CCLUSTER_L_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V6_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_L_V6_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V6_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_G_V6_OUT
    fi
    
    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:"                 $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:"                 $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V5=$(grep "time:"                 $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V5=$(grep "time:"                 $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V6=$(grep "time:"                 $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V6=$(grep "time:"                 $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_V6"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_L_V4 $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_V6"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    LINE_TAB2=$LINE_TAB2"&("$NBCLUSTERS_G_V4", "$NBSOLUTION_G_V4")"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V4", "$TSCCLUSTER_G_V4")&`format_time $TCCLUSTER_G_V4`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V5", "$TSCCLUSTER_G_V5")&`format_time $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V6", "$TSCCLUSTER_G_V6")&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V5 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    echo $LINE_TAB2 >> $TEMPTABFILE2
    
done

##########################naming
POL_NAME="Mandelbrot"
CCL_CALLV6="../../bin/MACIS19/ccluster_mandelbrot"

for DEG in $DEGMAND; do
#     LINE_TAB=$POL_NAME", \$d="$DEG"\$"
#     LINE_TAB2=$POL_NAME", \$d="$DEG"\$"
    LINE_TAB="\$\Man{"$DEG"}\$"
    LINE_TAB2="\$\Man{"$DEG"}\$"
    
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V5.out"
    FILENAME_CCLUSTER_G_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V5.out"
    FILENAME_CCLUSTER_L_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    $GPL_CALL $POL_NAME $DEG $FILEIN -b 14
    $GPL_CALL $POL_NAME $DEG $FILENAME_MPSOLVE_IN -b 14 -f 2
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_L_V5_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_G_V5_OUT
    fi
    
    if [ ! -e $FILENAME_CCLUSTER_L_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V6_OUT > /dev/stderr
        $CCL_CALLV6 $DEG -d $BLOCAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_L_V6_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V6_OUT > /dev/stderr
        $CCL_CALLV6 $DEG -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_G_V6_OUT
    fi
    
    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:"                 $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:"                 $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V5=$(grep "time:"                 $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V5=$(grep "time:"                 $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V6=$(grep "time:"                 $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V6=$(grep "time:"                 $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_V6"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_L_V4 $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_V6"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    LINE_TAB2=$LINE_TAB2"&("$NBCLUSTERS_G_V4", "$NBSOLUTION_G_V4")"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V4", "$TSCCLUSTER_G_V4")&`format_time $TCCLUSTER_G_V4`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V5", "$TSCCLUSTER_G_V5")&`format_time $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V6", "$TSCCLUSTER_G_V6")&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V5 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    echo $LINE_TAB2 >> $TEMPTABFILE2
    
done

##########################naming
POL_NAME="Runnels"
CCL_CALLV6="../../bin/MACIS19/ccluster_runnels"

for DEG in $DEGRUNN; do
#     LINE_TAB=$POL_NAME", \$d="$DEG"\$"
#     LINE_TAB2=$POL_NAME", \$d="$DEG"\$"
    LINE_TAB="\$\Run{"$DEG"}\$"
    LINE_TAB2="\$\Run{"$DEG"}\$"
    
    FILEIN=$REP"/"$POL_NAME"_"$DEG".ccl"
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V5.out"
    FILENAME_CCLUSTER_G_V5_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V5.out"
    FILENAME_CCLUSTER_L_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_V6_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    $GPL_CALL $POL_NAME $DEG $FILEIN -b 14
    $GPL_CALL $POL_NAME $DEG $FILENAME_MPSOLVE_IN -b 14 -f 2
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV4 -v "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BLOCAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_L_V5_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V5_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V5_OUT > /dev/stderr
        $CCL_CALL $FILEIN -d $BGLOBAL -e $EPSILONCCL -m $VFLAGV5 -v "3" > $FILENAME_CCLUSTER_G_V5_OUT
    fi
    
    if [ ! -e $FILENAME_CCLUSTER_L_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V6_OUT > /dev/stderr
        $CCL_CALLV6 $DEG -d $BLOCAL -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_L_V6_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V6_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V6_OUT > /dev/stderr
        $CCL_CALLV6 $DEG -e $EPSILONCCL -m $VFLAGV6 -v "3" > $FILENAME_CCLUSTER_G_V6_OUT
    fi
    
    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:"                 $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:"                 $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V4=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V4=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V4=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V5=$(grep "time:"                 $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V5=$(grep "time:"                 $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V5=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V5=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V5=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V5=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V5_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_V6=$(grep "time:"                 $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_L_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_L_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_L_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_V6=$(grep "time:"                 $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_V6=$(grep "tree depth:"          $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_V6=$(grep "tree size:"           $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V6=$(grep "number of clusters:"  $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_V6=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_V6_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_V6"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_L_V4 $TCCLUSTER_L_V6`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_V6"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    LINE_TAB2=$LINE_TAB2"&("$NBCLUSTERS_G_V4", "$NBSOLUTION_G_V4")"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V4", "$TSCCLUSTER_G_V4")&`format_time $TCCLUSTER_G_V4`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V5", "$TSCCLUSTER_G_V5")&`format_time $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V5`"
    LINE_TAB2=$LINE_TAB2"&("$TDCCLUSTER_G_V6", "$TSCCLUSTER_G_V6")&`format_time $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V5 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"&`ratio_time $TCCLUSTER_G_V4 $TCCLUSTER_G_V6`"
    LINE_TAB2=$LINE_TAB2"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    echo $LINE_TAB2 >> $TEMPTABFILE2
    
done

# echo $HEAD_TABLE
# echo $FIRST_LINE
# echo $SECOND_LINE
cat $TEMPTABFILE
# echo $TAIL_TAB

cat $TEMPTABFILE2

rm -rf $TEMPTABFILE
rm -rf $TEMPTABFILE2
