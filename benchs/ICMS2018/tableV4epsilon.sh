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
REP="tableV4epsilon"
V4FLAG="V4"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

COL_V41_TABLE=" V4, \$\epsilon = 2^{$EPS1}\$ "
COL_V42_TABLE=" V4, \$\epsilon = 2^{$EPS2}\$ "
COL_V43_TABLE=" V4, \$\epsilon = 2^{$EPS3}\$ "

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
rm -rf $TEMPTABFILE
touch $TEMPTABFILE

make_line()
{
    LINE=""
    FILEOUTE1=$1
    FILEOUTE2=$2
    FILEOUTE3=$3
    
TSV41=$(grep "tree size:" $FILEOUTE1| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV42=$(grep "tree size:" $FILEOUTE2| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TSV43=$(grep "tree size:" $FILEOUTE3| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TDV41=$(grep -m 1 "tree depth:" $FILEOUTE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV42=$(grep -m 1 "tree depth:" $FILEOUTE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
TDV43=$(grep -m 1 "tree depth:" $FILEOUTE3 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

N4V41=$(grep -m 1 "total number NE:" $FILEOUTE1 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V42=$(grep -m 1 "total number NE:" $FILEOUTE2 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
N4V43=$(grep -m 1 "total number NE:" $FILEOUTE3 | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

TV41=$(grep "total time:" $FILEOUTE1 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV42=$(grep "total time:" $FILEOUTE2 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
TV43=$(grep "total time:" $FILEOUTE3 | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"

}

CCL_CALL="../../bin/ccluster"
GPL_CALL="../../bin/genPolFile"

#-----------------------------------------------bernoulli-----------------------------------------
POL_NAME="Bernoulli"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

# #-----------------------------------------------bernoulli-----------------------------------------
# POL_NAME="Bernoulli"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$"  
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchBernoulli $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchBernoulli $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchBernoulli $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Mignotte-----------------------------------------
POL_NAME="Mignotte"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN -b $BITSIZE

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

# #-----------------------------------------------mignotte-----------------------------------------
# POL_NAME="Mignotte"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$"
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchMignotte $DEGREE $BITSIZE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchMignotte $DEGREE $BITSIZE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchMignotte $DEGREE $BITSIZE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------MignotteGen-----------------------------------------
POL_NAME="MignotteGen"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN -b $BITSIZE -p $POW

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

# #-----------------------------------------------MignotteGen-----------------------------------------
# POL_NAME="MignotteGen"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$, $\sigma=$BITSIZE\$, \$k=$POW\$"
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE, bitsize $BITSIZE, pow $POW version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchMignotteGen $DEGREE $BITSIZE $POW $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Wilkinson-----------------------------------------
POL_NAME="Wilkinson"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $DEGREE $FILEIN

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

# #-----------------------------------------------Wilkinson-----------------------------------------
# POL_NAME="Wilkinson"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$"
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchWilkinson $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchWilkinson $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchWilkinson $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------WilkMul-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
DEGREE=$(( $DEGREE / 2 ))

POL_NAME="WilkMul"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL" "$FILEIN" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

#generate input file
$GPL_CALL $POL_NAME $NBSOLS $FILEIN

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
DEGREE=$DEGSAVE

#-----------------------------------------------WilkMul-----------------------------------------
# DEGSAVE=$DEGREE
# DEGREE=$(( $NBSOLS * $(( $NBSOLS + 1 )) ))
# DEGREE=$(( $DEGREE / 2 ))
# 
# POL_NAME="WilkMul"
# LINE_TAB="$POL_NAME, \$nbSols=$NBSOLS\$, \$d=$DEGREE\$"
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchWilkMul $NBSOLS $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchWilkMul $NBSOLS $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchWilkMul $NBSOLS $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE
# DEGREE=$DEGSAVE

#-----------------------------------------------Spiral-----------------------------------------
POL_NAME="Spiral"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL"_spiral "$DEGREE" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL"_spiral "$DEGREE" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL"_spiral "$DEGREE" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE

# #-----------------------------------------------Spiral-----------------------------------------
# POL_NAME="Spiral"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$"
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchSpiral $DEGREE $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchSpiral $DEGREE $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, degree $DEGREE version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchSpiral $DEGREE $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE

#-----------------------------------------------Spiral-----------------------------------------
DEGSAVE=$DEGREE
DEGREE=$(( 3 ** $ITTS ))

POL_NAME="Spiral"
LINE_TAB="$POL_NAME, \$d=$DEGREE\$" 
FILEIN=$REP"/"$POL_NAME"_"$DEGREE".ccl"
FILEOUTE1=$REP"/"$POL_NAME"_e1.out"
FILEOUTE2=$REP"/"$POL_NAME"_e2.out"
FILEOUTE3=$REP"/"$POL_NAME"_e3.out"
CCL_CALLE1=$CCL_CALL"_nested "$ITTS" -d "$BOX" -e "$EPS1" -m "$V4FLAG" -v 3 " 
CCL_CALLE2=$CCL_CALL"_nested "$ITTS" -d "$BOX" -e "$EPS2" -m "$V4FLAG" -v 3 "
CCL_CALLE3=$CCL_CALL"_nested "$ITTS" -d "$BOX" -e "$EPS3" -m "$V4FLAG" -v 3 "

if [ ! -e $FILEOUTE1 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V1 output in "$FILEOUTE1 > /dev/stderr
    $CCL_CALLE1 > $FILEOUTE1
fi
if [ ! -e $FILEOUTE2 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V2 output in "$FILEOUTE2 > /dev/stderr
    $CCL_CALLE2 > $FILEOUTE2
fi
if [ ! -e $FILEOUTE3 ]; then
    echo  "Clustering roots for $POL_NAME, degree $DEGREE version V3 output in "$FILEOUTE3 > /dev/stderr
    $CCL_CALLE3 > $FILEOUTE3
fi

make_line $FILEOUTE1 $FILEOUTE2 $FILEOUTE3
echo $LINE_TAB $LINE"\\\\\\hline">> $TEMPTABFILE
DEGREE=$DEGSAVE

# #-----------------------------------------------Cluster-----------------------------------------
# DEGSAVE=$DEGREE
# DEGREE=$(( 3 ** $ITTS ))
# 
# POL_NAME="Cluster"
# LINE_TAB="$POL_NAME, \$d=$DEGREE\$"  
# 
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v41.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS1 output in "$REP"/"$POL_NAME"_"$DEGREE"_v41.out" > /dev/stderr
#     ./benchCluster $ITTS $BOX $EPS1 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v41.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v42.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS2 output in "$REP"/"$POL_NAME"_"$DEGREE"_v42.out" > /dev/stderr
#     ./benchCluster $ITTS $BOX $EPS2 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v42.out"
# fi
# if [ ! -e $REP"/"$POL_NAME"_"$DEGREE"_v43.out" ]; then
#     echo  "Clustering roots for $POL_NAME, $NBSOLS sols version V4 eps=2^$EPS3 output in "$REP"/"$POL_NAME"_"$DEGREE"_v43.out" > /dev/stderr
#     ./benchCluster $ITTS $BOX $EPS3 $V4FLAG "3" > $REP"/"$POL_NAME"_"$DEGREE"_v43.out"
# fi
# 
# TSV41=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV42=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TSV43=$(grep "tree size:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out"| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TDV41=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV42=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# TDV43=$(grep -m 1 "tree depth:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# N4V41=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V42=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# N4V43=$(grep -m 1 "total number NE:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
# 
# TV41=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v41.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV42=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v42.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# TV43=$(grep "total time:" $REP"/"$POL_NAME"_"$DEGREE"_v43.out" | cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
# LINE_TAB=$LINE_TAB"& ($TSV41,$TDV41,$N4V41) & `format_time $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV42,$TDV42,$N4V42) & `ratio_time $TV42 $TV41`"
# LINE_TAB=$LINE_TAB"& ($TSV43,$TDV43,$N4V43) & `ratio_time $TV43 $TV41`"
# 
# LINE_TAB=$LINE_TAB"\\\\\\hline"
# echo $LINE_TAB >> $TEMPTABFILE
# DEGREE=$DEGSAVE
# 
echo $HEAD_TABLE
echo $FIRST_LINE_TABLE
echo $SECOND_LINE_TABLE
cat $TEMPTABFILE
echo $TAIL_TAB

