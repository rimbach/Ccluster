#!/bin/bash
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
#     echo $TIME1
#     echo $TIME2
    STIME1=${#TIME1}
#     echo $STIME1
    STIME1=$(( 3 - $STIME1 ))
#     echo $STIME1
    STIME2=$(( 3 < $STIME1 ? 3 : $STIME1 ))
#     echo $STIME1
    STIME2=$(( 0 > $STIME2 ? 0 : $STIME2 ))
#     echo $STIME2
    if [ $STIME2 -eq 0 ]; then
        echo $TIME1
    else
#         echo $TIME2 | cut -c-$( echo $STIME2 )
        TIME2=`echo $TIME2 | cut -c-$( echo $STIME2)`
        echo $TIME1"."$TIME2
    fi
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
#    STOPWHENCOMPACT=0
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
REP="tableV4escape"
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

HEAD_TABLE="\begin{tabular}{c||"
FIRST_LINE_TABLE="    "
SECOND_LINE_TABLE="   "
TAIL_TAB="\end{tabular}"
# LINE_V3_TABLE=" (V3) "
LINE_V41_TABLE=" (V4), \$\epsilon = 2^{$EPS1}\$ "
LINE_V42_TABLE=" (V4), \$\epsilon = 2^{$EPS2}\$ "
LINE_V43_TABLE=" (V4), \$\epsilon = 2^{$EPS3}\$ "
LINE_V4P_TABLE=" (V4'), \$\epsilon = 2^{$EPSILONCCL}\$ "

#-----------------------------------------------bernoulli-----------------------------------------
POL_NAME="Bernoulli"

HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

#-----------------------------------------------mignotte-----------------------------------------
POL_NAME="Mignotte"

HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

#-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

HEAD_TABLE=$HEAD_TABLE"}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"\\\\\hline"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"\\\\\hline\hline"
LINE_V4P_TABLE=$LINE_V4P_TABLE"\\\\\hline"
LINE_V41_TABLE=$LINE_V41_TABLE"\\\\\hline"
LINE_V42_TABLE=$LINE_V42_TABLE"\\\\\hline"
LINE_V43_TABLE=$LINE_V43_TABLE"\\\\\hline"


echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
echo $LINE_V4P_TABLE
echo $LINE_V41_TABLE
echo $LINE_V42_TABLE
echo $LINE_V43_TABLE
echo $TAIL_TAB

HEAD_TABLE="\begin{tabular}{c||"
FIRST_LINE_TABLE="    "
SECOND_LINE_TABLE="   "
TAIL_TAB="\end{tabular}"
# LINE_V3_TABLE=" (V3) "
LINE_V41_TABLE=" (V4), \$\epsilon = 2^{$EPS1}\$ "
LINE_V42_TABLE=" (V4), \$\epsilon = 2^{$EPS2}\$ "
LINE_V43_TABLE=" (V4), \$\epsilon = 2^{$EPS3}\$ "
LINE_V4P_TABLE=" (V4'), \$\epsilon = 2^{$EPSILONCCL}\$ "

#-----------------------------------------------WilkMul-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
DEGREE=$(( $DEGREE / 2 ))

POL_NAME="WilkMul"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$nbSols=$NBSOLS\$, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

DEGREE=$DEGSAVE

#-----------------------------------------------MignotteGen-----------------------------------------
POL_NAME="MignotteGen"

HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$, \$k=$POW\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

#-----------------------------------------------Cluster-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( 3 ** $ITTS ))

POL_NAME="Cluster"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& TSize & Tdepth & n4 & time"  

if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4p eps=2^$EPSILONCCL output in "$REP"/"$POL_NAME"_"$DEGREE"_v4p.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPSILONCCL $V4PFLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
fi
if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
fi

TSV4P=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV4P=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V4P=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV4P=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v4p.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V4P_TABLE=$LINE_V4P_TABLE" & $TSV4P & $TDV4P & $N4V4P & `format_time $TV4P` s "
LINE_V41_TABLE=$LINE_V41_TABLE" & $TSV41 & $TDV41 & $N4V41 & `format_time $TV41` s "
LINE_V42_TABLE=$LINE_V42_TABLE" & $TSV42 & $TDV42 & $N4V42 & `format_time $TV42` s "
LINE_V43_TABLE=$LINE_V43_TABLE" & $TSV43 & $TDV43 & $N4V43 & `format_time $TV43` s "

DEGREE=$DEGSAVE

HEAD_TABLE=$HEAD_TABLE"}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"\\\\\hline"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"\\\\\hline\hline"
LINE_V4P_TABLE=$LINE_V4P_TABLE"\\\\\hline"
LINE_V41_TABLE=$LINE_V41_TABLE"\\\\\hline"
LINE_V42_TABLE=$LINE_V42_TABLE"\\\\\hline"
LINE_V43_TABLE=$LINE_V43_TABLE"\\\\\hline"


echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
echo $LINE_V4P_TABLE
echo $LINE_V41_TABLE
echo $LINE_V42_TABLE
echo $LINE_V43_TABLE
echo $TAIL_TAB

