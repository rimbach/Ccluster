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
   DEGREES="64 128 191"
   DEGREES="64 128 191 256"
   DEGREES="64 128 191 256 383"
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

if [ -z "$BGLOBAL" ]; then
   BGLOBAL="0,1,0,1,1000,1"
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

VFLAG=7
if [ $STOPWHENCOMPACT -eq 1 ]; then
    VFLAG=$(( $VFLAG + 8 ))
fi
if [ $ANTICIPATE -eq 1 ]; then
    VFLAG=$(( $VFLAG + 16 ))
fi

VFLAGV4="V4"
VFLAGT="test"

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

##########################naming
POL_NAME="Bernoulli"
CCLUSTER_CALL="../test/bernoulli"
GENPOLFILE_CALL="./genPolFile "$POL_NAME

for DEG in $DEGREES; do
    LINE_TAB=$DEG
    
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_T_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_T_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BLOCAL $EPSILONCCL $VFLAGV4 "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BGLOBAL $EPSILONCCL $VFLAGV4 "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_T_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_T_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BLOCAL $EPSILONCCL $VFLAGT "3" > $FILENAME_CCLUSTER_L_T_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_T_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_T_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BGLOBAL $EPSILONCCL $VFLAGT "3" > $FILENAME_CCLUSTER_G_T_OUT
    fi
    
    $GENPOLFILE_CALL $DEG > $FILENAME_MPSOLVE_IN
    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_G=$(grep "tree depth:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_G=$(grep "tree size:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_G=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_T=$(grep "time:" $FILENAME_CCLUSTER_L_T_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_T=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_T_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_T=$(grep "time:" $FILENAME_CCLUSTER_G_T_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_G=$(grep "tree depth:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_G=$(grep "tree size:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_T=$(grep "number of clusters:" $FILENAME_CCLUSTER_G_T_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_G=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_T"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_T`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_T"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_T`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    
done

##########################naming
POL_NAME="Mignotte"
CCLUSTER_CALL="../test/mignottePS"
GENPOLFILE_CALL="./genPolFile "$POL_NAME

for DEG in $DEGREES; do
    LINE_TAB=$DEG
    
    FILENAME_CCLUSTER_L_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_V4.out"
    FILENAME_CCLUSTER_G_V4_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_V4.out"
    FILENAME_CCLUSTER_L_T_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_local_test.out"
    FILENAME_CCLUSTER_G_T_OUT=$REP"/"$POL_NAME"_"$DEG"_ccluster_global_test.out"
    FILENAME_MPSOLVE_IN=$REP"/"$POL_NAME"_"$DEG".pol"
    FILENAME_MPSOLVE_S_OUT=$REP"/"$POL_NAME"_"$DEG"_mpsolve_s.out"
    
    if [ ! -e $FILENAME_CCLUSTER_L_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_V4_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BITSIZE $BLOCAL $EPSILONCCL $VFLAGV4 "3" > $FILENAME_CCLUSTER_L_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_V4_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_V4_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BITSIZE $BGLOBAL $EPSILONCCL $VFLAGV4 "3" > $FILENAME_CCLUSTER_G_V4_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_L_T_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, local, output in "$FILENAME_CCLUSTER_L_T_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BITSIZE $BLOCAL $EPSILONCCL $VFLAGT "3" > $FILENAME_CCLUSTER_L_T_OUT
    fi
    if [ ! -e $FILENAME_CCLUSTER_G_T_OUT ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG, global, output in "$FILENAME_CCLUSTER_G_T_OUT > /dev/stderr
        $CCLUSTER_CALL $DEG $BITSIZE $BGLOBAL $EPSILONCCL $VFLAGT "3" > $FILENAME_CCLUSTER_G_T_OUT
    fi
    
    $GENPOLFILE_CALL $DEG $BITSIZE > $FILENAME_MPSOLVE_IN
    if [ $MPSOLVE -eq 1 ]; then
        if [ ! -e $FILENAME_MPSOLVE_S_OUT ]; then
            echo "Isolating roots in C with MPSOLVE SECSOLVE........." > /dev/stderr
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $FILENAME_MPSOLVE_IN > $FILENAME_MPSOLVE_S_OUT) &>> $FILENAME_MPSOLVE_S_OUT
        fi
    fi
    
    TCCLUSTER_L_V4=$(grep "time:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_V4=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G_V4=$(grep "time:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_G=$(grep "tree depth:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_G=$(grep "tree size:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_V4=$(grep "number of clusters:" $FILENAME_CCLUSTER_G_V4_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_G=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_L_T=$(grep "time:" $FILENAME_CCLUSTER_L_T_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_L_T=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_T_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TCCLUSTER_G_T=$(grep "time:" $FILENAME_CCLUSTER_G_T_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_G=$(grep "tree depth:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_G=$(grep "tree size:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_T=$(grep "number of clusters:" $FILENAME_CCLUSTER_G_T_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_G=$(grep "number of solutions:" $FILENAME_CCLUSTER_G_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE=$(grep "real" $FILENAME_MPSOLVE_S_OUT| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_L_T"&`format_time $TCCLUSTER_L_V4`"
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_L_T`"
    LINE_TAB=$LINE_TAB"&"$NBCLUSTERS_G_T"&`format_time $TCCLUSTER_G_V4`" 
    LINE_TAB=$LINE_TAB"&`format_time $TCCLUSTER_G_T`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
    echo $LINE_TAB >> $TEMPTABFILE
    
done

# echo $HEAD_TABLE
# echo $FIRST_LINE
# echo $SECOND_LINE
cat $TEMPTABFILE
# echo $TAIL_TAB

rm -rf $TEMPTABFILE
