#!/bin/bash

usage()
{
   echo "Usage: ./tableV1V3 <options> <args>"
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

if [ -z "$BOX" ]; then
   BOX="0,1,0,1,100,1"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

##########################constants
TRUE=1
FALSE=0
REP="tableV1V3"
V1FLAG="V1"
V3FLAG="V3"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

HEAD_TABLE="\begin{tabular}{l||c|c||c|c||}"
FIRST_LINE_TABLE="     &\multicolumn{2}{c||}{ V1 }&\multicolumn{2}{c||}{ V3 }\\\\\\hline"
SECOND_LINE_TABLE="     & (n1, n2, n3) & tV1 & (n1, n2, n3) & tV3/tV1\\\\\\hline\\hline "
TAIL_TAB="\end{tabular}"

LINES_TAB=""
TEMPTABFILE="temptabfile.txt"
touch $TEMPTABFILE

make_line()
{
    LINE=""
    FILEOUTV1=$1
    FILEOUTV2=$2
    N1V1=$(grep "tree size:" $FILEOUTV1| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    N1V3=$(grep "tree size:" $FILEOUTV3| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    N2V1=$(grep -m 1 "conclusion" $FILEOUTV1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    N2V3=$(grep -m 1 "conclusion" $FILEOUTV3 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    N3V1=$(grep -m 1 "total number GR:" $FILEOUTV1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    N3V3=$(grep -m 1 "total number GR:" $FILEOUTV3 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

    TV1=$(grep "total time:" $FILEOUTV1 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TV3=$(grep "total time:" $FILEOUTV3 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

    LINE_TAB=$LINE_TAB"& ($N1V1,$N2V1,$N3V1) & `format_time $TV1`"
    LINE_TAB=$LINE_TAB"& ($N1V3,$N2V3,$N3V3) & `ratio_time $TV1 $TV3`"
}

CCL_CALL="../../bin/ccluster"
GPL_CALL="../../bin/genPolFile"

#-----------------------------------------------bernoulli-----------------------------------------
POL_NAME="Bernoulli"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"
CCL_CALLV1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 "  
CCL_CALLV3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi

make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------Mignotte-----------------------------------------
POL_NAME="Mignotte"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"
CCL_CALLV1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 "
CCL_CALLV3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN -b $BITSIZE

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------MignotteGen-----------------------------------------
POL_NAME="MignotteGen"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"
CCL_CALLV1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 "  
CCL_CALLV3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN -b $BITSIZE -p $POW

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"
CCL_CALLV1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 " 
CCL_CALLV3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------WilkMul-----------------------------------------
DEGTEMP=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
DEGTEMP=$(( $DEGTEMP / 2 ))

POL_NAME="WilkMul"
LINE_TAB="$POL_NAME, \$d=$DEGTEMP\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"
CCL_CALLV1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 " 
CCL_CALLV3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $NBSOLS $FILEIN

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGTEMP version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGTEMP version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"

CCL_CALLV1=$CCL_CALL"_spiral "$DEGREE" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 "   
CCL_CALLV3=$CCL_CALL"_spiral "$DEGREE" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

#-----------------------------------------------Cluster-----------------------------------------
DEGTEMP=$(( 3 ** $ITTS ))

POL_NAME="Cluster"
LINE_TAB="$POL_NAME, \$d=$DEGTEMP\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTV1=$REP"/"$POL_NAME"_v1.out"
FILEOUTV3=$REP"/"$POL_NAME"_v3.out"

CCL_CALLV1=$CCL_CALL"_nested "$ITTS" -d "$BOX" -e "$EPSILONCCL" -m "$V1FLAG" -v 3 "      
CCL_CALLV3=$CCL_CALL"_nested "$ITTS" -d "$BOX" -e "$EPSILONCCL" -m "$V3FLAG" -v 3 "

if [ ! -e $FILEOUTV1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGTEMP version V1 output in "$FILEOUTV1 > /dev/stderr
    $CCL_CALLV1 > $FILEOUTV1
fi
if [ ! -e $FILEOUTV3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGTEMP version V3 output in "$FILEOUTV3 > /dev/stderr
    $CCL_CALLV3 > $FILEOUTV3
fi
make_line $FILEOUTV1 $FILEOUTV3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE


echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
cat $TEMPTABFILE
echo $TAIL_TAB

rm -rf $TEMPTABFILE

