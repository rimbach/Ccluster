#!/bin/bash

echo "EPSILON AS AN ESCAPE BOUND IS A DEPRECATED OPTION"
exit -1;

# bernoulli: deg 64, eps = -53 -530 -5300
# mignotte : deg 64, bitsize = 14, eps = -53, -530, -5300
# #            deg 128, bitsize = 14, eps = -53, -530, -5300
# cluster  : deg 81, eps = -2, -4, -16

# wilkinson: deg 64, eps = -53, -530, -5300
# wilkMul  : deg 66, eps = -53, -530, -5300
# spiral   : deg 64, eps = -53, -530, -5300
#mignotteGen: deg64, eps = -53, -530, -5300
usage()
{
   echo "Usage: ./tableV4escape <options> <args>"
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
    GT1=`echo $NUM'>'$DEN|bc -l`
    RATIO=""
    if [ $GT1 -eq 0 ]; then
        RATIO=0`echo $NUM/$DEN|bc -l`
    else
        RATIO=`echo $NUM/$DEN|bc -l`
    fi
    echo `format_time $RATIO`
}
   
##########################getting arguments
# DEGREE=$1
# EPSILONCCL=$2
# BITSIZE=$3
# STOPWHENCOMPACT=$4
# BOX=$5

while [ "$1" != "" ]; do
   PARAM=`echo $1 | sed 's/=.*//'`
   VALUE=`echo $1 | sed 's/[^=]*//; s/=//'`
   case "$PARAM" in
      -h|--help)
         usage
         exit 0
         ;;
      --degree)
        DEGREE=$VALUE
        ;;
      --epsilon)
        EPSILONCCL=$VALUE
        ;;
      --eps1)
        EPS1=$VALUE
        ;;
      --eps2)
        EPS2=$VALUE
        ;;
      --eps3)
        EPS3=$VALUE
        ;;
      --bitsize)
        BITSIZE=$VALUE
        ;;
      --nbSols)
        NBSOLS=$VALUE
        ;;
      --pow)
        POW=$VALUE
        ;;
      --itts)
        ITTS=$VALUE
        ;;
#       --stop-when-compact)
#         STOPWHENCOMPACT=1
#         ;;
      --box)
        BOX="$VALUE"
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
if [ -z "$DEGREE" ]; then
   DEGREE="64"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$EPS1" ]; then
   EPS1="-53"
fi

if [ -z "$EPS2" ]; then
   EPS2="-530"
fi

if [ -z "$EPS3" ]; then
   EPS3="-5300"
fi

if [ -z "$BITSIZE" ]; then
   BITSIZE="14"
fi

if [ -z "$NBSOLS" ]; then
   NBSOLS="11"
fi

if [ -z "$POW" ]; then
   POW="3"
fi

if [ -z "$ITTS" ]; then
   ITTS="4"
fi

# if [ -z "$STOPWHENCOMPACT" ]; then
#    STOPWHENCOMPACT=1
# fi

if [ -z "$BOX" ]; then
   BOX="0,1,0,1,100,1"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

# echo $DEGREE
# echo $EPSILONCCL
# echo $BITSIZE
# echo $STOPWHENCOMPACT
# echo $BOX
##########################constants
TRUE=1
FALSE=0
REP="tableV4esc_epsilon"
V4FLAG=23
V4PFLAG=31
# if [ $STOPWHENCOMPACT -eq 1 ]; then
#     V3FLAG=$(( $V3FLAG + 8 ))
#     V4FLAG=$(( $V4FLAG + 8 ))
# fi

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

COL_V41_TABLE=" V4', \$\epsilon = 2^{$EPS1}\$ "
COL_V42_TABLE=" V4', \$\epsilon = 2^{$EPS2}\$ "
COL_V43_TABLE=" V4', \$\epsilon = 2^{$EPS3}\$ "

HEAD_TABLE="\begin{tabular}{l||c|c||c|c||c|c||}"
FIRST_LINE_TABLE="&\multicolumn{2}{c||}{$COL_V41_TABLE}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{2}{c||}{$COL_V42_TABLE}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{2}{c||}{$COL_V43_TABLE}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"\\\\\\hline"
SECOND_LINE_TABLE="& (n1, TD, n4) & t53 "
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& (n1, TD, n4) & t530/t53 "
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& (n1, TD, n4) & t5300/t53 "
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"\\\\\\hline\\hline"
TAIL_TAB="\end{tabular}"

TEMPTABFILE="temptabfile.txt"

touch $TEMPTABFILE

#-----------------------------------------------bernoulli-----------------------------------------
POL_NAME="Bernoulli"
CCL_CALL="../test/bernoulli"

LINE_TAB="$POL_NAME, \$d=$DEGREE\$"  
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------mignotte-----------------------------------------
POL_NAME="Mignotte"
CCL_CALL="../test/mignotte"

LINE_TAB="$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$"

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"
CCL_CALL="../test/wilkinson"

LINE_TAB="$POL_NAME, \$d=$DEGREE\$"

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"
CCL_CALL="../test/spiral"

LINE_TAB="$POL_NAME, \$d=$DEGREE\$"

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $DEGREE $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------WilkMul-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
DEGREE=$(( $DEGREE / 2 ))

POL_NAME="WilkMul"
CCL_CALL="../test/wilkMul"

LINE_TAB="$POL_NAME, \$nbSols=$NBSOLS\$, \$d=$DEGREE\$"

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $NBSOLS $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $NBSOLS $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $NBSOLS $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE
DEGREE=$DEGSAVE

#-----------------------------------------------MignotteGen-----------------------------------------
POL_NAME="MignotteGen"
CCL_CALL="../test/mignotteGen"

LINE_TAB="$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$, \$k=$POW\$"

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $POW $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $POW $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $DEGREE $BITSIZE $POW $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Cluster-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( 3 ** $ITTS ))

POL_NAME="NestedClusters"
CCL_CALL="../test/cluster"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    $CCL_CALL $ITTS $BOX $EPS1 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    $CCL_CALL $ITTS $BOX $EPS2 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    $CCL_CALL $ITTS $BOX $EPS3 $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

LINE_TAB=$LINE_TAB"\\\\\\hline"
echo $LINE_TAB >> $TEMPTABFILE
DEGREE=$DEGSAVE
# 
echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE
