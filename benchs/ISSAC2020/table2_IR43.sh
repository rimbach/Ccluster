#!/bin/bash

usage()
{
   echo "Usage: ./table2_IR43.sh <options> <args>"
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
   DEGREES="64 128 191"
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

if [ -z "$GENERATE" ]; then
   GENERATE=0
fi

if [ -z "$NBPOLS" ]; then
   NBPOLS=100
fi

if [ -z "$NBITT" ]; then
   NBITT="7 8 9"
fi

if [ -z "$SIZEGRID" ]; then
   SIZEGRID="5 6"
fi

if [ -z "$NBTERMS" ]; then
   NBTERMS=10
fi

CCLUSTER_PATH="../../"
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
CCLUSTER_EXPE=$CCLUSTER_PATH"/bin/ISSAC20/ccluster_issac20"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
CCLUSTER_OPTS="-v 2"

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
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2
    fi
    
}

gen_with_deg_bs(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    BS=$4
    NAME_IN=$NAME".ccl"
    NAME_IN2=$NAME".mpl"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2 -b $BS
    fi
    
}

gen_with_deg_bs_nbterms(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    BS=$4
    NBT=$5
    NAME_IN=$NAME".ccl"
    NAME_IN2=$NAME".mpl"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS -n $NBT
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2 -b $BS -n $NBT
    fi
    
}

run_ccluster()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out"
    NAME_OUTE=$NAME".out_expe"
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CCLUSTER_CALL $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUTE ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, experimental, output in " $NAME_OUTE
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
            CALL="$CCLUSTER_EXPE $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUTE
    fi
    
}

run_ccluster_mandelbrot()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out"
    NAME_OUTE=$NAME".out_expe"
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CCLUSTER_CALL $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUTE ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, experimental, output in " $NAME_OUTE
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
#             CALL="$CCLUSTER_EXPE $NAME_IN -v 2 -m default"
#             $CALL > $NAME_OUTE
              CALL="$CCLUSTER_PATH/bin/ISSAC20/ccluster_issac20_mandelbrot $DEG -v 2"
              $CALL > $NAME_OUTE
    fi
    
}

stats_pol()
{
    NAME=$1
    NAME_OUT=$NAME".out"
    NAME_OUTE=$NAME".out_expe"
    
    TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_NBFAILS=$(grep "failure:"                 $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TTIME=$(grep "total time:"                   $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TSIZE=$(grep "tree size:"                    $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TDEPT=$(grep "tree depth:"                    $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     EXPE_TINEVAL=$(grep "time in Evaluation:"                $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TINPOST=$(grep "time in Ps counting tests:"                $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    LINE_TAB=$TSIZE" & "$TDEPT" & `format_time $TTIME`"
    LINE_TAB=$LINE_TAB" & "$EXPE_NBFAILS" & "$EXPE_TSIZE" & "$EXPE_TDEPT" & `format_time $EXPE_TTIME` & `percent_time $EXPE_TTIME $TTIME`"
#     LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINEVAL $EXPE_TTIME` & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    echo $LINE_TAB
}

stats_pol_rand()
{
    NAME=$1
    NAME_OUT=$NAME".out"
    NAME_OUTE=$NAME".out_expe"
    
    TSIZE_T=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_T=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_T=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_NBFAILS_T=$(grep "failure:"                 $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TTIME_T=$(grep "total time:"                   $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TSIZE_T=$(grep "tree size:"                    $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TDEPT_T=$(grep "tree depth:"                    $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     EXPE_TINEVAL_T=$(grep "time in Evaluation:"                $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    EXPE_TINPOST_T=$(grep "time in Ps counting tests:"                $NAME_OUTE| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TSIZE=`echo $TSIZE+$TSIZE_T|bc -l`
    TDEPT=`echo $TDEPT+$TDEPT_T|bc -l`
    TTIME=`echo $TTIME+$TTIME_T|bc -l`
    EXPE_NBFAILS=`echo $EXPE_NBFAILS+$EXPE_NBFAILS_T|bc -l`
    EXPE_TTIME=`echo $EXPE_TTIME+$EXPE_TTIME_T|bc -l`
    EXPE_TSIZE=`echo $EXPE_TSIZE+$EXPE_TSIZE_T|bc -l`
    EXPE_TDEPT=`echo $EXPE_TDEPT+$EXPE_TDEPT_T|bc -l`
    EXPE_TINPOST=`echo $EXPE_TINPOST+$EXPE_TINPOST_T|bc -l`
}

REP="table2_IR43"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

TEMPTABFILE="temptab2_IR43.txt"
touch $TEMPTABFILE

POLNAME="randomDense"

#solve random polynomials with ccluster
echo $POLNAME >> $TEMPTABFILE
for DEG in $DEGREES; do

    REPNAME=$REP"/"$POLNAME"_"$DEG
    FILENAME=$POLNAME"_"$DEG
    
    if [ ! -d $REPNAME ]; then
#         echo $REPNAME
         mkdir -p $REPNAME
    fi
    
    TSIZE=0
    TDEPT=0
    TTIME=0
    EXPE_NBFAILS=0
    EXPE_TTIME=0
    EXPE_TSIZE=0
    EXPE_TDEPT=0
    EXPE_TINPOST=0
    
    for CURIND in `seq 1 $NBPOLS`; do
        NAME=$REPNAME"/"$FILENAME"_"$CURIND
        gen_with_deg_bs $NAME $POLNAME $DEG $BITSIZE
        run_ccluster $NAME $POLNAME $DEG
#         cp $REP"/rand_"$DEG"/rand_"$DEG"_"$CURIND".ccl" $NAME."ccl"
#         cp $REP"/rand_"$DEG"/rand_"$DEG"_"$CURIND".out" $NAME."out"
#         cp $REP"/rand_"$DEG"/rand_"$DEG"_"$CURIND".out0" $NAME."out0"
#         cp $REP"/rand_"$DEG"/rand_"$DEG"_"$CURIND".out1" $NAME."out1"
#         cp $REP"/rand_"$DEG"/rand_"$DEG"_"$CURIND".out2" $NAME."out2"
        stats_pol_rand $NAME
    done
    
    LINE_TAB=" "$DEG" & "$TSIZE" & "$TDEPT" & `format_time $TTIME`"
    LINE_TAB=$LINE_TAB" & "$EXPE_NBFAILS" & "$EXPE_TSIZE" & "$EXPE_TDEPT" & `format_time $EXPE_TTIME` & `percent_time $EXPE_TTIME $TTIME`"
#     LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINEVAL $EXPE_TTIME` & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    LINE_TAB=$LINE_TAB"\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
    
done
# # 
#Other polynomials
# DEGREES="64 128 191 256 383 512"
# DEGREES="64"
POLNAMES="Bernoulli Chebyshev1 Legendre"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_ccluster $NAME $POLNAME $DEG
#     gen_and_run_ccluster $NAME
    
    LINE_TAB=" "$DEG" & "
    LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
    
#     stats_pol $NAME
done 
done
# 


POLNAME="RegularGrid"

echo $POLNAME >> $TEMPTABFILE
for SIZ in $SIZEGRID; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$SIZ
    
    gen_with_deg $NAME $POLNAME $SIZ
    run_ccluster $NAME $POLNAME $SIZ
#     gen_and_run_ccluster $NAME
    
    LINE_TAB=" "$SIZ" & "
    LINE_TAB=$LINE_TAB"`stats_pol $REPNAME"/"$POLNAME"_"$SIZ`\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
    
#     stats_pol $REPNAME"/"$POLNAME"_"$DEG
done

POLNAME="randomSparse"

#solve random sparse polynomials with ccluster
echo $POLNAME >> $TEMPTABFILE
for DEG in $DEGREES; do

    REPNAME=$REP"/"$POLNAME"_"$DEG
    FILENAME=$POLNAME"_"$DEG
    
    if [ ! -d $REPNAME ]; then
#         echo $REPNAME
         mkdir -p $REPNAME
    fi
    
    TSIZE=0
    TDEPT=0
    TTIME=0
    EXPE_NBFAILS=0
    EXPE_TTIME=0
    EXPE_TSIZE=0
    EXPE_TDEPT=0
    EXPE_TINEVAL=0
    EXPE_TINPOST=0
    
#     REP2="tableISSAC1sparse_IR43"
    
    for CURIND in `seq 1 $NBPOLS`; do
        NAME=$REPNAME"/"$FILENAME"_"$CURIND
        gen_with_deg_bs_nbterms $NAME $POLNAME $DEG $BITSIZE $NBTERMS
        run_ccluster $NAME $POLNAME $DEG
#         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".ccl" $NAME."ccl"
#         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out" $NAME."out"
#         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out0" $NAME."out0"
#         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out1" $NAME."out1"
#         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out2" $NAME."out2"

        stats_pol_rand $NAME
    done
    
    LINE_TAB=" "$DEG" & "$TSIZE" & "$TDEPT" & `format_time $TTIME`"
    LINE_TAB=$LINE_TAB" & "$EXPE_NBFAILS" & "$EXPE_TSIZE" & "$EXPE_TDEPT" & `format_time $EXPE_TTIME` & `percent_time $EXPE_TTIME $TTIME`"
    LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINEVAL $EXPE_TTIME` & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
    LINE_TAB=$LINE_TAB"\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
    
#     LINE_TAB=" "$DEG" & "$NUMDISTESTS" & `percent_time $TTIMEINDIST $TTIMEINCCLU`"
#     LINE_TAB=$LINE_TAB" & "$NUMTN0" & "$NUMFP0" & `ratio_time $TTIME0 $TTIMEINDIST`"
#     LINE_TAB=$LINE_TAB" & "$NUMTN1" & "$NUMFP1" & `ratio_time $TTIME1 $TTIMEINDIST`"
#     LINE_TAB=$LINE_TAB" & "$NUMTN2" & "$NUMFP2" & `ratio_time $TTIME2 $TTIMEINDIST`"
#     LINE_TAB=$LINE_TAB"\\\\"
# #     echo $LINE_TAB
#     echo $LINE_TAB >> $TEMPTABFILE
done

#Other polynomials
# DEGREES="64 128 191 256 383 512"
# DEGREES="64"
POLNAMES="Mignotte"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg_bs $NAME $POLNAME $DEG $BITSIZE
    run_ccluster $NAME $POLNAME $DEG
#     gen_and_run_ccluster $NAME
    
    LINE_TAB=" "$DEG" & "
    LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
#     stats_pol $NAME
done 
done

#Procedural polynomials

POLNAMES="Mandelbrot"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE
for DEG in $NBITT; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
#     run_ccluster $NAME $POLNAME $DEG
    run_ccluster_mandelbrot $NAME $POLNAME $DEG
#     gen_and_run_ccluster $NAME
    
    LINE_TAB=" "$DEG" & "
    LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
    echo $LINE_TAB >> $TEMPTABFILE
#     stats_pol $NAME
done 
done


cat $TEMPTABFILE
rm -f $TEMPTABFILE
