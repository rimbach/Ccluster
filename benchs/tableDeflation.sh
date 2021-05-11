#!/bin/bash

usage()
{
   echo "Usage: ./tableDeflation.sh <options> <args>"
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
      --bitsize)
        BITSIZE=$VALUE
        ;;
      --nbpols)
        NBPOLS=$VALUE
        ;;
      --nbitts)
        NBITT=$VALUE
        ;;
      --sizegrid)
        SIZEGRID=$VALUE
        ;;
      --nbterms)
        NBTERMS=$VALUE
        ;;
      --epsilonCCL)
        EPSILONCCL=$VALUE
        ;;
      --purge)
        PURGE=1
        ;;
      --purgeRiso)
        PURGERISO=1
        ;;
      --purgeRisoDef)
        PURGERISODE=1
        ;;
      --purgeCclu)
        PURGECCLU=1
        ;;
      --purgeCcluDef)
        PURGECCLUDE=1
        ;;
      --purgeDSC)
        PURGEDSC=1
        ;;
        --purgeMPS)
        PURGEMPS=1
        ;;
      --generate)
        GENERATE=1
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
   DEGREES="128 191 256 391 512"
   LENP=5
fi

if [ -z "$BITSIZE" ]; then
   BITSIZE="8"
fi

if [ -z "$EPSILONCCL" ]; then
   EPSILONCCL="-53"
fi

if [ -z "$PURGE" ]; then
   PURGE=0
fi

if [ -z "$PURGERISO" ]; then
   PURGERISO=0
fi

if [ -z "$PURGERISODE" ]; then
   PURGERISODE=0
fi

if [ -z "$PURGECCLU" ]; then
   PURGECCLU=0
fi

if [ -z "$PURGECCLUDE" ]; then
   PURGECCLUDE=0
fi

if [ -z "$PURGEDSC" ]; then
   PURGEDSC=0
fi

if [ -z "$PURGEMPS" ]; then
   PURGEMPS=0
fi

if [ -z "$GENERATE" ]; then
   GENERATE=0
fi

if [ -z "$NBPOLS" ]; then
#    NBPOLS=20
   NBPOLS=10
#    NBPOLS=1
fi

# if [ -z "$NBITT" ]; then
#    NBITT="7 8 9"
# fi

# if [ -z "$SIZEGRID" ]; then
# #    SIZEGRID="5 6"
# #    SIZEGRID="5 6 7 8 10 11 12 13 14"
# #    SIZEGRID="6 7 8 9 10 11 12 13 14"
#    SIZEGRID="6 8 10 12 14 16"
# #    SIZEGRID="6 8"
# fi

# if [ -z "$NBTERMS" ]; then
#    NBTERMS=10
# fi

CCLUSTER_PATH="../"
RISOLATE_CALL=$CCLUSTER_PATH"/bin/risolate"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
GENRANDPOLFI_CALL=$CCLUSTER_PATH"/bin/genRandPolFile"
RISOLATE_OPTS="-v 2"
CCLUSTER_OPTS="-v 2"

MPSOLVE_CALL="mpsolve -as -Ga -j1"

ANEWDSC_PATH="../../softs"
# ANEWDSC_CALL=$ANEWDSC_PATH"/test_descartes_linux64 --sqrfree 1"
ANEWDSC_CALL=$ANEWDSC_PATH"/test_descartes_linux64"

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
        echo $TIME1"."
    else
        TIME2=`echo $TIME2 | cut -c-$( echo $STIME2)`
        echo $TIME1"."$TIME2
    fi
}

format_numb()
{
#     printf "%$2d" $1
    NUMB1=$1
    SNUMB1=${#NUMB1}
#     echo $SNUMB1
    SDIFF1=$(( $2 - $SNUMB1 ))
#     echo $SDIFF1
    if [ $SDIFF1 -le 0 ]; then
        echo $NUMB1
    else 
        RES=""
        while [ $SDIFF1 -gt 0 ]
        do
            RES=$RES" "
            SDIFF1=$(( $SDIFF1-1 ))
        done
        RES="$RES$1"
        echo "$RES"
    fi
}

ratio_time()
{
    NUM=$1
    DEN=$2
    RATIO=`echo $NUM/$DEN|bc -l`
    echo `format_time $RATIO`
}

percent_time()
{
    NUM=$1
    DEN=$2
    RATIO=`echo 100*$NUM/$DEN|bc -l`
    echo `format_time $RATIO`
}

gen_with_deg(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".ccl"
    NAME_IN2=$NAME".mpl"
    NAME_IN3=$NAME".dsc"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1
    fi
    
    if [ ! -e $NAME_IN3 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN3
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN3 -f 3
    fi
    
}

gen_with_deg_bs(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    BS=$4
    NAME_IN=$NAME".ccl"
    NAME_IN2=$NAME".mpl"
    NAME_IN3=$NAME".dsc"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS
    fi
    
    if [ ! -e $NAME_IN2 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN2
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2 -b $BS
    fi
    
    if [ ! -e $NAME_IN3 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN3
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN3 -f 3 -b $BS
    fi
    
}

run_risolate()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_riso"
    NAME_OUTDE=$NAME".out_riso_de"
    
    if [ $PURGERISO -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, bitsize $BIT, global, without deflation, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m V7"
            $CALL > $NAME_OUT
    fi
    
    if [ $PURGERISODE -eq 1 ]; then
        rm -f $NAME_OUTDE
    fi
  
    if [ ! -e $NAME_OUTDE ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, bitsize $BIT, global, with deflation, output in " $NAME_OUTDE
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m onlySubd"
            $CALL > $NAME_OUTDE
    fi
    
}

run_ccluster()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    EPS=$5
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_cclu"
    NAME_OUTDE=$NAME".out_cclu_de"
    
    if [ $PURGECCLU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, bitsize $BIT, epsilon $EPS, global, without deflation, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CCLUSTER_CALL $NAME_IN -e $EPS -v 2 -m V7"
            $CALL > $NAME_OUT
    fi
    
    if [ $PURGECCLUDE -eq 1 ]; then
        rm -f $NAME_OUTDE
    fi
  
    if [ ! -e $NAME_OUTDE ]; then
            echo  "Clustering roots for $POLNAME degree $DEG,bitsize $BIT, epsilon $EPS, global, with deflation, output in " $NAME_OUTDE
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
            CALL="$CCLUSTER_CALL $NAME_IN -e $EPS -v 2 -m onlySubd"
            $CALL > $NAME_OUTDE
    fi
    
}

run_aNewDsc()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    SQUAREFREE=$5
    NAME_IN=$NAME".dsc"
    NAME_OUT=$NAME".out_dsc"
    
    if [ $PURGEDSC -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    CALL=$ANEWDSC_CALL
    if [ $SQUAREFREE -eq 1 ]; then
        CALL=$CALL" --sqrfree 1"
    fi
    
#     echo $CALL
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, bitsize $BIT, with aNewDsc, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $CALL $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

run_MPsolve()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    EPS=$5
    NAME_IN=$NAME".mpl"
    NAME_OUT=$NAME".out_mpl"
    
    if [ $PURGEMPS -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Approximating roots for $POLNAME degree $DEG, bitsize $BIT, with MPSolve, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL "-o"$EPS $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

TEMPTABFILE1="temptab_deflation1.txt"
touch $TEMPTABFILE1
TEMPTABFILE2="temptab_deflation2.txt"
touch $TEMPTABFILE2

stats_pol()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTDE=$NAME".out_riso_de"
    NAME_OUTDSC=$NAME".out_dsc"
    DEG=$2
    
    BITSI=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NSOLS=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBEXT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBCOT=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_NSOLS=$(grep "number of real roots:"       $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_TTIME=$(grep "total time:"                   $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TSIZE=$(grep "tree size:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TDEPT=$(grep "tree depth:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBEXT=$(grep "total number DT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBCOT=$(grep "total number VT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

    DSC_NSOLS=$(grep "Number of roots:" $NAME_OUTDSC| cut -f2 -d':'| tr -d ' ')
#     DSC_TSIZE=$(grep "TREESIZE=" $NAME_OUTDSC| cut -f2 -d'='| tr -d ' ')
    DSC_TTIME=$(grep "real" $NAME_OUTDSC| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $NSOLS $LENP` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DE_NSOLS $LENP` & `format_time $DE_TTIME` & `percent_time $DE_TTIME $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_NSOLS $LENP` & `format_time $DSC_TTIME`\\\\"
    
    echo $LINE_TAB1 >> $TEMPTABFILE1
}

stats_pol_only_risolate()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTDE=$NAME".out_riso_de"
    DEG=$2
    
    BITSI=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NSOLS=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBEXT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBCOT=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_NSOLS=$(grep "number of real roots:"       $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_TTIME=$(grep "total time:"                   $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TSIZE=$(grep "tree size:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TDEPT=$(grep "tree depth:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBEXT=$(grep "total number DT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBCOT=$(grep "total number VT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $NSOLS $LENP` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DE_NSOLS $LENP` & `format_time $DE_TTIME` & `percent_time $DE_TTIME $TTIME`\\\\"
    
    echo $LINE_TAB1 >> $TEMPTABFILE1
}

stats_pol_ccluster()
{
    NAME=$1
    NAME_OUT=$NAME".out_cclu"
    NAME_OUTDE=$NAME".out_cclu_de"
    NAME_OUTMPL=$NAME".out_mpl"
    DEG=$2
    
    BITSI=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NCLUS=$(grep "number of clusters:"           $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NSOLS=$(grep "number of solutions:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBEXT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBCOT=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_CLUS=$(grep "number of clusters:"        $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_NSOLS=$(grep "number of solutions:"       $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    DE_TTIME=$(grep "total time:"                   $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TSIZE=$(grep "tree size:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_TDEPT=$(grep "tree depth:"                    $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBEXT=$(grep "total number DT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     DE_NBCOT=$(grep "total number VT:"              $NAME_OUTDE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

#     MPL_NSOLS=$(grep "Number of roots:" $NAME_OUTMPL| cut -f2 -d':'| tr -d ' ')
    MPL_TTIME=$(grep "real" $NAME_OUTMPL| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $NCLUS $LENP`  & `format_numb $NSOLS $LENP` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DE_NCLUS $LENP` & `format_numb $DE_NSOLS $LENP` & `format_time $DE_TTIME` & `percent_time $DE_TTIME $TTIME`"
    LINE_TAB1=$LINE_TAB1" &`format_time $MPL_TTIME`\\\\"
    
    echo $LINE_TAB1 >> $TEMPTABFILE2
}

REP="tableDeflation"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

DEGREES="256 391"
BITSIZES="256 512 1024 2048"
POLNAME="Mignotte"

echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2

for DEG in $DEGREES; do
    for BIT in $BITSIZES; do
    
        REPNAME=$REP
        NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT
    
        gen_with_deg_bs $NAME $POLNAME $DEG $BIT
        run_risolate $NAME $POLNAME $DEG $BIT
        run_aNewDsc  $NAME $POLNAME $DEG $BIT 1
    
        stats_pol $NAME $DEG
        echo "$POLNAME & $DEG & $BIT & `format_numb $NSOLS` & `format_numb $DE_NSOLS` & `format_numb $DSC_NSOLS` "
    done 
done

NBSOLS="10 15 17"
POLNAME="WilkMul"

echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2

for NBSOL in $NBSOLS; do
    
        REPNAME=$REP
        NAME=$REPNAME"/"$POLNAME"_"$NBSOL
    
        gen_with_deg $NAME $POLNAME $NBSOL
        run_risolate $NAME $POLNAME $NBSOL "?"
#         run_aNewDsc  $NAME $POLNAME $DEG 1
        
        run_ccluster $NAME $POLNAME $NBSOL "?" "-100000"
        run_MPsolve  $NAME $POLNAME $NBSOL "?" "30190"
    
        stats_pol_only_risolate $NAME $NBSOL
        stats_pol_ccluster $NAME $NBSOL
        echo "$POLNAME & $NBSOL & `format_numb $NSOLS` & `format_numb $DE_NSOLS` "
done

cat $TEMPTABFILE1
rm -f $TEMPTABFILE1
cat $TEMPTABFILE2
rm -f $TEMPTABFILE2
