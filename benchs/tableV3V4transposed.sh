#!/bin/bash

usage()
{
   echo "Usage: ./tableV3V4 <options> <args>"
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

ratio_time()
{
    NUM=$1
    DEN=$2
#     RATIO=0`echo $NUM/$DEN|bc -l`
    RATIO=`echo $NUM/$DEN|bc -l`
#     echo $RATIO
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
   DEGREE="64 128 256"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$BITSIZE" ]; then
   BITSIZE="14"
fi

if [ -z "$NBSOLS" ]; then
   NBSOLS="11 12 13"
fi

if [ -z "$POW" ]; then
   POW="3"
fi

if [ -z "$ITTS" ]; then
   ITTS="3 4 5"
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
REP="tableV3V4"
V3FLAG=7
V4FLAG=23
if [ $STOPWHENCOMPACT -eq 1 ]; then
    V3FLAG=$(( $V3FLAG + 8 ))
    V4FLAG=$(( $V4FLAG + 8 ))
fi

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

HEAD_TABLE="\begin{tabular}{l||c|c||c|c||}"
FIRST_LINE_TABLE="     &\multicolumn{2}{c||}{ V3 }&\multicolumn{2}{c||}{ V4 }\\\\\\hline"
# SECOND_LINE_TABLE="     & (n1, n2, n3) & tV3 & (n1, n2, n3) & tV4/tV3\\\\\\hline\\hline "
SECOND_LINE_TABLE="     & n3 & tV3 & n3 & tV3/tV4\\\\\\hline\\hline "
TAIL_TAB="\end{tabular}"

LINES_TAB=""
TEMPTABFILE="temptabfile.txt"
rm -rf $TEMPTABFILE
touch $TEMPTABFILE

make_line()
{
    LINE=""
    FILE1=$1
    FILE2=$2
#     N1V3=$(grep "tree size:" $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     N1V4=$(grep "tree size:" $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     N2V3=$(grep -m 1 "conclusion" $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     N2V4=$(grep -m 1 "conclusion" $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    N3V3=$(grep -m 1 "total number GR:" $FILE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    N3V4=$(grep -m 1 "total number GR:" $FILE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TV3=$(grep "total time:" $FILE1 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TV4=$(grep "total time:" $FILE2 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     LINE=$LINE"& ($N1V3,$N2V3,$N3V3) & `format_time $TV3`"
#     LINE=$LINE"& ($N1V4,$N2V4,$N3V4) & `ratio_time $TV4 $TV3`"
    LINE=$LINE"& $N3V3 & `format_time $TV3`"
    LINE=$LINE"& $N3V4 & `ratio_time $TV3 $TV4`"
}

#-----------------------------------------------bernoulli-----------------------------------------

POL_NAME="Bernoulli"

for DEG in $DEGREE; do
    LINE_TAB="Bernoulli, \$d=$DEG\$"  
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchBernoulli $DEG $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchBernoulli $DEG $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
      make_line $FILE1 $FILE2
      echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

#-----------------------------------------------Mignotte-----------------------------------------
POL_NAME="Mignotte"

for DEG in $DEGREE; do
    LINE_TAB="Mignotte, \$d=$DEG\$, \$\sigma=$BITSIZE\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchMignotte $DEG $BITSIZE $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchMignotte $DEG $BITSIZE $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done
# 
#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"

for DEG in $DEGREE; do
    LINE_TAB="Wilkinson, \$d=$DEG\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchWilkinson $DEG $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchWilkinson $DEG $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

# #-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"

for DEG in $DEGREE; do
    LINE_TAB="Spiral, \$d=$DEG\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchSpiral $DEG $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchSpiral $DEG $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done


# #-----------------------------------------------WilkMul-----------------------------------------
POL_NAME="WilkMul"

for NBSOL in $NBSOLS; do
    
    DEG=$(( $NBSOL * $(( $NBSOL + 1 )) ))
    DEG=$(( $DEG / 2 ))

    LINE_TAB="Wilkinson Mult., \$nbSols=$NBSOL\$, \$d=$DEG\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, $NBSOL sols version V3 output in "$FILE1 > /dev/stderr
        ./benchWilkMul $NBSOL $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, $NBSOL sols version V4 output in "$FILE2 > /dev/stderr
        ./benchWilkMul $NBSOL $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

#-----------------------------------------------MignotteGen-----------------------------------------
POL_NAME="MignotteGen"

for DEG in $DEGREE; do
    LINE_TAB="Mignotte Mult., \$d=$DEG\$, \$\sigma=$BITSIZE\$, \$k=$POW\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchMignotteGen $DEG $BITSIZE $POW $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchMignotteGen $DEG $BITSIZE $POW $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done
# 
# #-----------------------------------------------Cluster-----------------------------------------
POL_NAME="NestedClustrs"

for ITT in $ITTS; do
    
    DEG=$(( 3 ** $ITT ))

    LINE_TAB="Nested Clusters, \$d=$DEG\$"
    FILE1=$REP"/"$POL_NAME"_"$DEG"_v3.out"
    FILE2=$REP"/"$POL_NAME"_"$DEG"_v4.out"
    if [ ! -e $FILE1 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V3 output in "$FILE1 > /dev/stderr
        ./benchCluster $ITT $BOX $EPSILONCCL $V3FLAG "3" > $FILE1
    fi
    if [ ! -e $FILE2 ]; then
        echo  "Clustering roots for $POL_NAME, degree $DEG version V4 output in "$FILE2 > /dev/stderr
        ./benchCluster $ITT $BOX $EPSILONCCL $V4FLAG "3" > $FILE2
    fi
    make_line $FILE1 $FILE2
    echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
done

echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
cat $TEMPTABFILE
echo $TAIL_TAB