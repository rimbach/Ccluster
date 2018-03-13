#!/bin/bash

usage()
{
   echo "Usage: ./tableV1V2V3 <options> <args>"
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
      --stop-when-compact)
        STOPWHENCOMPACT=1
        ;;
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

if [ -z "$STOPWHENCOMPACT" ]; then
   STOPWHENCOMPACT=0
fi

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
REP="tableV1V2V3"
V1FLAG=5
V2FLAG=519
V3FLAG=7
if [ $STOPWHENCOMPACT -eq 1 ]; then
    V1FLAG=$(( $V1FLAG + 8 ))
    V2FLAG=$(( $V2FLAG + 8 ))
    V3FLAG=$(( $V3FLAG + 8 ))
fi

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

# echo $V1FLAG
# echo $V2FLAG
# echo $V3FLAG

HEAD_TABLE="\begin{tabular}{c||"
FIRST_LINE_TABLE="    "
SECOND_LINE_TABLE="   "
TAIL_TAB="\end{tabular}"

LINE_V1_TABLE=" (V1) "
LINE_V2_TABLE=" (V2) "
LINE_V3_TABLE=" (V3) "

#-----------------------------------------------bernoulli-----------------------------------------
POL_NAME="Bernoulli"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  


if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchBernoulli $DEGREE $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

#-----------------------------------------------Mignotte-----------------------------------------
POL_NAME="Mignotte"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  

if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchMignotte $DEGREE $BITSIZE $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  

if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchWilkinson $DEGREE $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

#-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGREE\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  

if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchSpiral $DEGREE $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "


HEAD_TABLE=$HEAD_TABLE"}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"\\\\\hline"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"\\\\\hline\hline"
LINE_V1_TABLE=$LINE_V1_TABLE"\\\\\hline"
LINE_V2_TABLE=$LINE_V2_TABLE"\\\\\hline"
LINE_V3_TABLE=$LINE_V3_TABLE"\\\\\hline"

echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
echo $LINE_V1_TABLE
echo $LINE_V2_TABLE
echo $LINE_V3_TABLE
echo $TAIL_TAB


HEAD_TABLE="\begin{tabular}{c||"
FIRST_LINE_TABLE="    "
SECOND_LINE_TABLE="   "
TAIL_TAB="\end{tabular}"

LINE_V1_TABLE=" (V1) "
LINE_V2_TABLE=" (V2) "
LINE_V3_TABLE=" (V3) "

#-----------------------------------------------WilkMul-----------------------------------------
DEGTEMP=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
DEGTEMP=$(( $DEGTEMP / 2 ))

POL_NAME="WilkMul"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$nbSols=$NBSOLS\$, \$d=$DEGTEMP\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  


if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchWilkMul $NBSOLS $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

#-----------------------------------------------MignotteGen-----------------------------------------

POL_NAME="MignotteGen"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$, \$k=$POW\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  

if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

#-----------------------------------------------Cluster-----------------------------------------
DEGTEMP=$(( 3 ** $ITTS ))

POL_NAME="Cluster"
HEAD_TABLE=$HEAD_TABLE"c|c|c|c||"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"&\multicolumn{4}{|c||}{$POL_NAME, \$d=$DEGTEMP\$}"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"& n1 & n2 & n3 & time"  

if [ ! -e $REP"/"$POL_NAME"_v1.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$REP"/"$POL_NAME"_v1.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPSILONCCL $V1FLAG "3" > $REP"/"$POL_NAME"_v1.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v2.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$REP"/"$POL_NAME"_v2.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPSILONCCL $V2FLAG "3" > $REP"/"$POL_NAME"_v2.out"
fi
if [ ! -e $REP"/"$POL_NAME"_v3.out" ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$REP"/"$POL_NAME"_v3.out" > /dev/stderr
    ./benchCluster $ITTS $BOX $EPSILONCCL $V3FLAG "3" > $REP"/"$POL_NAME"_v3.out"
fi
N1V1=$(grep "tree size:" $REP"/"$POL_NAME"_v1.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V2=$(grep "tree size:" $REP"/"$POL_NAME"_v2.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N1V3=$(grep "tree size:" $REP"/"$POL_NAME"_v3.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N2V1=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V2=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N2V3=$(grep -m 1 "conclusion" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N3V1=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V2=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N3V3=$(grep -m 1 "total number GR:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV1=$(grep "total time:" $REP"/"$POL_NAME"_v1.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV2=$(grep "total time:" $REP"/"$POL_NAME"_v2.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV3=$(grep "total time:" $REP"/"$POL_NAME"_v3.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_V1_TABLE=$LINE_V1_TABLE" & $N1V1 & $N2V1 & $N3V1 & `format_time $TV1` s "
LINE_V2_TABLE=$LINE_V2_TABLE" & $N1V2 & $N2V2 & $N3V2 & `format_time $TV2` s "
LINE_V3_TABLE=$LINE_V3_TABLE" & $N1V3 & $N2V3 & $N3V3 & `format_time $TV3` s "

HEAD_TABLE=$HEAD_TABLE"}"
FIRST_LINE_TABLE=$FIRST_LINE_TABLE"\\\\\hline"
SECOND_LINE_TABLE=$SECOND_LINE_TABLE"\\\\\hline\hline"
LINE_V1_TABLE=$LINE_V1_TABLE"\\\\\hline"
LINE_V2_TABLE=$LINE_V2_TABLE"\\\\\hline"
LINE_V3_TABLE=$LINE_V3_TABLE"\\\\\hline"

echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
echo $LINE_V1_TABLE
echo $LINE_V2_TABLE
echo $LINE_V3_TABLE
echo $TAIL_TAB


