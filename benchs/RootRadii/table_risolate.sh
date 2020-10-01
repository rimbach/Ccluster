#!/bin/bash

usage()
{
   echo "Usage: ./table_risolate.sh <options> <args>"
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
#    DEGREES="64 128 191"
#    DEGREES="128 191 256 391 512 791 1024"
   DEGREES="128 191 256 391 512 791"
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

if [ -z "$GENERATE" ]; then
   GENERATE=0
fi

if [ -z "$NBPOLS" ]; then
   NBPOLS=10
fi

if [ -z "$NBITT" ]; then
   NBITT="7 8 9"
fi

if [ -z "$SIZEGRID" ]; then
#    SIZEGRID="5 6"
   SIZEGRID="5 6 7 8 10 11 12 13 14"
fi

if [ -z "$NBTERMS" ]; then
   NBTERMS=10
fi

# CCLUSTER_PATH="../../"
# CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
# CCLUSTER_EXPE=$CCLUSTER_PATH"/bin/ISSAC20/ccluster_issac20"
# GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
# CCLUSTER_OPTS="-v 2"

CCLUSTER_PATH="../../"
RISOLATE_CALL=$CCLUSTER_PATH"/bin/risolate"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
GENRANDPOLFI_CALL=$CCLUSTER_PATH"/bin/genRandPolFile"
RISOLATE_OPTS="-v 2"

ANEWDSC_PATH="../../../softs"
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
    
    if [ ! -e $NAME_IN3 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN3
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN3 -f 3 -b $BS
    fi
    
}

genRand_with_deg_bs(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    BS=$4
    NBPOLS=$5
    LOC=$6
    NAME_IN=$NAME"_*.ccl"
    NAME_IN2=$NAME"_*.mpl"
    NAME_IN3=$NAME"_*.dsc"
    
#     if [ ! -e $NAME_IN ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 1 -b $BS -p $NBPOLS -l $LOC
#     fi
    
#     if [ ! -e $NAME_IN3 ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN3
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 3 -b $BS -p $NBPOLS -l $LOC
#     fi
    
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

run_risolate()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUTRR ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, root radii, output in " $NAME_OUTRR
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m V7"
            $CALL > $NAME_OUTRR
    fi
    
}

run_aNewDsc()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".dsc"
    NAME_OUT=$NAME".out_dsc"
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, with aNewDsc, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $ANEWDSC_CALL $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

TEMPTABFILE1="temptab_risolate1.txt"
touch $TEMPTABFILE1
TEMPTABFILE2="temptab_risolate2.txt"
touch $TEMPTABFILE2



# run_ccluster_mandelbrot()
# {
#     NAME=$1
#     POLNAME=$2
#     DEG=$3
#     NAME_IN=$NAME".ccl"
#     NAME_OUT=$NAME".out"
#     NAME_OUTE=$NAME".out_expe"
#     
#     if [ ! -e $NAME_OUT ]; then
#             echo  "Clustering roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
# #             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
#             CALL="$CCLUSTER_CALL $NAME_IN -v 2 -m default"
#             $CALL > $NAME_OUT
#     fi
#     
#     if [ ! -e $NAME_OUTE ]; then
#             echo  "Clustering roots for $POLNAME degree $DEG, global, experimental, output in " $NAME_OUTE
# #             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
# #             CALL="$CCLUSTER_EXPE $NAME_IN -v 2 -m default"
# #             $CALL > $NAME_OUTE
#               CALL="$CCLUSTER_PATH/bin/ISSAC20/ccluster_issac20_mandelbrot $DEG -v 2"
#               $CALL > $NAME_OUTE
#     fi
#     
# }

stats_pol()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    NAME_OUTDSC=$NAME".out_dsc"
    DEG=$2
    
    BITSI=$(grep "bitsize of input polynomial:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NSOLS=$(grep "number of solutions:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate
#     NBTZT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBTST=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate RR
    RR_NSOLS=$(grep "number of solutions:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_NBVTS=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TIVTS=$(grep "total time spent in tests VT:" $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PRECN=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PPREC=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_NBGRA=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBGRR=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_TINGR=$(grep "time in Graeffe iterations:"    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TINRR=$(grep "time in computing root radii:"   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

    DSC_NSOLS=$(grep "Number of roots:" $NAME_OUTDSC| cut -f2 -d':'| tr -d ' ')
    DSC_TSIZE=$(grep "TREESIZE=" $NAME_OUTDSC| cut -f2 -d'='| tr -d ' ')
    DSC_TTIME=$(grep "real" $NAME_OUTDSC| cut -f2 -d'l' | tr -d ' ')
    
#     LINE_TAB="`format_numb $NSOLS $LENP` & `format_numb $TSIZE $LENP` & `format_numb $TDEPT 2` & `format_time $TTIME`"
# #     LINE_TAB=$LINE_TAB" & `format_numb $NBTZT $LENP` & `format_numb $NBTST $LENP`"
#     LINE_TAB=$LINE_TAB" & `format_numb $RR_NSOLS $LENP` & `format_numb $RR_TSIZE $LENP` & `format_numb $RR_TDEPT 2`"
# #     LINE_TAB=$LINE_TAB" & `format_numb $RR_NBVTS $LENP` &  `format_time $RR_TIVTS`"
#     LINE_TAB=$LINE_TAB" & `format_numb $RR_PRECN $LENP` "
#     LINE_TAB=$LINE_TAB" & `format_time $RR_TINRR` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
#     
#     LINE_TAB=$LINE_TAB" & `format_numb $DSC_NSOLS $LENP` & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`"
    
    LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_numb $NSOLS $LENP` & `format_numb $TSIZE $LENP` & `format_numb $TDEPT 2` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $RR_TSIZE $LENP` & `format_numb $RR_TDEPT 2` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`\\\\"
    
    LINE_TAB2="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_time $RR_TTIME` & `percent_time $RR_TINRR $RR_TTIME`"
    LINE_TAB2=$LINE_TAB2" & `format_numb $RR_PPREC $LENP` & `format_numb $RR_PRECN $LENP`"
    LINE_TAB2=$LINE_TAB2" & `format_numb $RR_NBGRA $LENP` & `format_numb $RR_NBGRR $LENP`"
    LINE_TAB2=$LINE_TAB2" & `percent_time $RR_TINGR $RR_TTIME`\\\\"  
    
    echo $LINE_TAB1 >> $TEMPTABFILE1
    echo $LINE_TAB2 >> $TEMPTABFILE2
}

stats_pol_rand()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    NAME_OUTDSC=$NAME".out_dsc"
    
    BITSI_T=$(grep "bitsize of input polynomial:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NSOLS_T=$(grep "number of solutions:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE_T=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_T=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_T=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
        # details on risolate #
#     NBTZT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBTST_T=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate RR #
    RR_NSOLS_T=$(grep "number of solutions:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME_T=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE_T=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT_T=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_NBVTS=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TIVTS=$(grep "total time spent in tests VT:" $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PRECN_T=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PPREC_T=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_NBGRA_T=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBGRR_T=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_TINGR_T=$(grep "time in Graeffe iterations:"    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TINRR_T=$(grep "time in computing root radii:"   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

    DSC_NSOLS_T=$(grep "Number of roots:" $NAME_OUTDSC| cut -f2 -d':'| tr -d ' ')
    DSC_TSIZE_T=$(grep "TREESIZE=" $NAME_OUTDSC| cut -f2 -d'='| tr -d ' ')
    DSC_TTIME_T=$(grep "real" $NAME_OUTDSC| cut -f2 -d'l' | tr -d ' ')
    
    BITSI=`echo $BITSI+$BITSI_T|bc -l`
    NSOLS=`echo $NSOLS+$NSOLS_T|bc -l`
    TSIZE=`echo $TSIZE+$TSIZE_T|bc -l`
    TDEPT=`echo $TDEPT+$TDEPT_T|bc -l`
    TTIME=`echo $TTIME+$TTIME_T|bc -l`
#     NBTZT=`echo $NBTZT+$NBTZT_T|bc -l`
#     NBTST=`echo $NBTST+$NBTST_T|bc -l`
    RR_NSOLS=`echo $RR_NSOLS +$RR_NSOLS_T|bc -l`
    RR_TTIME=`echo $RR_TTIME +$RR_TTIME_T|bc -l`
    RR_TSIZE=`echo $RR_TSIZE +$RR_TSIZE_T|bc -l`
    RR_TDEPT=`echo $RR_TDEPT +$RR_TDEPT_T|bc -l`
    RR_PRECN=`echo $RR_PRECN +$RR_PRECN_T|bc -l`
    RR_PPREC=`echo $RR_PPREC +$RR_PPREC_T|bc -l`
    RR_NBGRA=`echo $RR_NBGRA +$RR_NBGRA_T|bc -l`
    RR_NBGRR=`echo $RR_NBGRR +$RR_NBGRR_T|bc -l`
    RR_TINGR=`echo $RR_TINGR +$RR_TINGR_T|bc -l`
    RR_TINRR=`echo $RR_TINRR +$RR_TINRR_T|bc -l`
    DSC_NSOLS=`echo $DSC_NSOLS+$DSC_NSOLS_T|bc -l`
    DSC_TSIZE=`echo $DSC_TSIZE+$DSC_TSIZE_T|bc -l`
    DSC_TTIME=`echo $DSC_TTIME+$DSC_TTIME_T|bc -l`
    
}

REP="tab_risolate"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

POLNAME="randomDense"

#solve random polynomials with ccluster
echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2
NBPOLS=1
for DEG in $DEGREES; do

#     REPNAME=$REP"/"$POLNAME"_"$DEG
#     FILENAME=$POLNAME"_"$DEG
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITSIZE
    
#     if [ ! -d $REPNAME ]; then
# #         echo $REPNAME
#          mkdir -p $REPNAME
#     fi
    
    BITSI=0
    NSOLS=0
    TSIZE=0
    TDEPT=0
    TTIME=0
#     NBTZT=0
#     NBTST=0
    RR_NSOLS=0
    RR_TTIME=0
    RR_TSIZE=0
    RR_TDEPT=0
    RR_PRECN=0
    RR_PPREC=0
    RR_NBGRA=0
    RR_NBGRR=0
    RR_TINGR=0
    RR_TINRR=0
    DSC_NSOLS=0
    DSC_TSIZE=0
    DSC_TTIME=0
    
    genRand_with_deg_bs $NAME $POLNAME $DEG $BITSIZE $NBPOLS $REPNAME
    for CURIND in `seq 1 $NBPOLS`; do
        NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITSIZE"_"$CURIND
        run_risolate $NAME $POLNAME $DEG
        run_aNewDsc  $NAME $POLNAME $DEG
        stats_pol_rand $NAME
    done
    
    BITSI=`echo     $BITSI    /$NBPOLS     |bc -l`
    NSOLS=`echo     $NSOLS    /$NBPOLS     |bc -l`
    TSIZE=`echo     $TSIZE    /$NBPOLS     |bc -l`
    TDEPT=`echo     $TDEPT    /$NBPOLS     |bc -l`
    TTIME=`echo     $TTIME    /$NBPOLS     |bc -l`
#     NBTZT=`echo     $NBTZT    /$NBPOLS     |bc -l`
#     NBTST=`echo     $NBTST    /$NBPOLS     |bc -l`
    RR_NSOLS=`echo  $RR_NSOLS /$NBPOLS     |bc -l`
    RR_TTIME=`echo  $RR_TTIME /$NBPOLS     |bc -l`
    RR_TSIZE=`echo  $RR_TSIZE /$NBPOLS     |bc -l`
    RR_TDEPT=`echo  $RR_TDEPT /$NBPOLS     |bc -l`
    RR_PRECN=`echo  $RR_PRECN /$NBPOLS     |bc -l`
    RR_PPREC=`echo  $RR_PPREC /$NBPOLS     |bc -l`
    RR_NBGRA=`echo  $RR_NBGRA /$NBPOLS     |bc -l`
    RR_NBGRR=`echo  $RR_NBGRR /$NBPOLS     |bc -l`
    RR_TINGR=`echo  $RR_TINGR /$NBPOLS     |bc -l`
    RR_TINRR=`echo  $RR_TINRR /$NBPOLS     |bc -l`
    DSC_NSOLS=`echo $DSC_NSOLS/$NBPOLS     |bc -l`
    DSC_TSIZE=`echo $DSC_TSIZE/$NBPOLS     |bc -l`
    DSC_TTIME=`echo $DSC_TTIME/$NBPOLS     |bc -l`
    
#     LINE_TAB=" `format_numb $DEG $LENP` & `format_time $NSOLS `"
#     LINE_TAB=$LINE_TAB" &`format_time $NSOLS ` & `format_time $TSIZE` & `format_time $TDEPT` & `format_time $TTIME`"
#     
#     LINE_TAB=$LINE_TAB" & `format_time $RR_NSOLS` & `format_time $RR_TSIZE` & `format_time $RR_TDEPT`"
# #     LINE_TAB=$LINE_TAB" & `format_numb $RR_NBVTS` &  `format_time $RR_TIVTS`"
#     LINE_TAB=$LINE_TAB" & `format_time $RR_PRECN` "
#     LINE_TAB=$LINE_TAB" & `format_time $RR_TINRR` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
#     
#     LINE_TAB=$LINE_TAB" & `format_time $DSC_NSOLS` & `format_time $DSC_TSIZE` &`format_time $DSC_TTIME`"
#     echo $LINE_TAB >> $TEMPTABFILE
    
    LINE_TAB1="$DEG & `format_time $NSOLS` & `format_time $BITSI` & `format_time $TSIZE` & `format_time $TDEPT` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_time $RR_TSIZE` & `format_time $RR_TDEPT` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`\\\\"
    
#     LINE_TAB2="$DEG & `format_time $RR_TTIME` & `format_time $RR_PRECN` & `format_time $RR_TINRR` & `percent_time $RR_TINRR $RR_TTIME`\\\\" 
    LINE_TAB2="`format_numb $DEG $LENP` & `format_time $BITSI` & `format_time $RR_TTIME` & `percent_time $RR_TINRR $RR_TTIME`"
    LINE_TAB2=$LINE_TAB2" & `format_time $RR_PPREC` & `format_time $RR_PRECN`"
    LINE_TAB2=$LINE_TAB2" & `format_time $RR_NBGRA` & `format_time $RR_NBGRR`"
    LINE_TAB2=$LINE_TAB2" & `percent_time $RR_TINGR $RR_TTIME`\\\\"   
    
    echo $LINE_TAB1 >> $TEMPTABFILE1
    echo $LINE_TAB2 >> $TEMPTABFILE2
    
done

# POLNAME="randomDense"
# 
# #solve random polynomials with ccluster
# echo $POLNAME >> $TEMPTABFILE1
# echo $POLNAME >> $TEMPTABFILE2
# DEGF=256
# BITSIZES="10 50 100 500 1000 5000"
# for BIT in $BITSIZES; do
# 
#     REPNAME=$REP"/"$POLNAME"_"$DEGF"_"$BIT
#     FILENAME=$POLNAME"_"$DEG"_"$BIT
#     
#     if [ ! -d $REPNAME ]; then
# #         echo $REPNAME
#          mkdir -p $REPNAME
#     fi
#     
#     BITSI=0
#     NSOLS=0
#     TSIZE=0
#     TDEPT=0
#     TTIME=0
# #     NBTZT=0
# #     NBTST=0
#     RR_NSOLS=0
#     RR_TTIME=0
#     RR_TSIZE=0
#     RR_TDEPT=0
#     RR_PRECN=0
#     RR_PPREC=0
#     RR_NBGRA=0
#     RR_NBGRR=0
#     RR_TINGR=0
#     RR_TINRR=0
#     DSC_NSOLS=0
#     DSC_TSIZE=0
#     DSC_TTIME=0
#     
#     for CURIND in `seq 1 $NBPOLS`; do
#         NAME=$REPNAME"/"$FILENAME"_"$CURIND
#         gen_with_deg_bs $NAME $POLNAME $DEGF $BIT
#         run_risolate $NAME $POLNAME $DEGF
#         run_aNewDsc  $NAME $POLNAME $DEGF
#         stats_pol_rand $NAME
#     done
#     
#     BITSI=`echo     $BITSI    /$NBPOLS     |bc -l`
#     NSOLS=`echo     $NSOLS    /$NBPOLS     |bc -l`
#     TSIZE=`echo     $TSIZE    /$NBPOLS     |bc -l`
#     TDEPT=`echo     $TDEPT    /$NBPOLS     |bc -l`
#     TTIME=`echo     $TTIME    /$NBPOLS     |bc -l`
# #     NBTZT=`echo     $NBTZT    /$NBPOLS     |bc -l`
# #     NBTST=`echo     $NBTST    /$NBPOLS     |bc -l`
#     RR_NSOLS=`echo  $RR_NSOLS /$NBPOLS     |bc -l`
#     RR_TTIME=`echo  $RR_TTIME /$NBPOLS     |bc -l`
#     RR_TSIZE=`echo  $RR_TSIZE /$NBPOLS     |bc -l`
#     RR_TDEPT=`echo  $RR_TDEPT /$NBPOLS     |bc -l`
#     RR_PRECN=`echo  $RR_PRECN /$NBPOLS     |bc -l`
#     RR_PPREC=`echo  $RR_PPREC /$NBPOLS     |bc -l`
#     RR_NBGRA=`echo  $RR_NBGRA /$NBPOLS     |bc -l`
#     RR_NBGRR=`echo  $RR_NBGRR /$NBPOLS     |bc -l`
#     RR_TINGR=`echo  $RR_TINGR /$NBPOLS     |bc -l`
#     RR_TINRR=`echo  $RR_TINRR /$NBPOLS     |bc -l`
#     DSC_NSOLS=`echo $DSC_NSOLS/$NBPOLS     |bc -l`
#     DSC_TSIZE=`echo $DSC_TSIZE/$NBPOLS     |bc -l`
#     DSC_TTIME=`echo $DSC_TTIME/$NBPOLS     |bc -l`
#     
# #     LINE_TAB=" `format_numb $DEG $LENP` & `format_time $NSOLS `"
# #     LINE_TAB=$LINE_TAB" &`format_time $NSOLS ` & `format_time $TSIZE` & `format_time $TDEPT` & `format_time $TTIME`"
# #     
# #     LINE_TAB=$LINE_TAB" & `format_time $RR_NSOLS` & `format_time $RR_TSIZE` & `format_time $RR_TDEPT`"
# # #     LINE_TAB=$LINE_TAB" & `format_numb $RR_NBVTS` &  `format_time $RR_TIVTS`"
# #     LINE_TAB=$LINE_TAB" & `format_time $RR_PRECN` "
# #     LINE_TAB=$LINE_TAB" & `format_time $RR_TINRR` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
# #     
# #     LINE_TAB=$LINE_TAB" & `format_time $DSC_NSOLS` & `format_time $DSC_TSIZE` &`format_time $DSC_TTIME`"
# #     echo $LINE_TAB >> $TEMPTABFILE
#     
#     LINE_TAB1="$DEG & `format_time $NSOLS` & `format_time $BITSI` & `format_time $TSIZE` & `format_time $TDEPT` & `format_time $TTIME`"
#     LINE_TAB1=$LINE_TAB1" & `format_time $RR_TSIZE` & `format_time $RR_TDEPT` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
#     LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`\\\\"
#     
# #     LINE_TAB2="$DEG & `format_time $RR_TTIME` & `format_time $RR_PRECN` & `format_time $RR_TINRR` & `percent_time $RR_TINRR $RR_TTIME`\\\\" 
#     LINE_TAB2="`format_numb $DEG $LENP` & `format_time $BITSI` & `format_time $RR_TTIME` & `percent_time $RR_TINRR $RR_TTIME`"
#     LINE_TAB2=$LINE_TAB2" & `format_time $RR_PPREC` & `format_time $RR_PRECN`"
#     LINE_TAB2=$LINE_TAB2" & `format_time $RR_NBGRA` & `format_time $RR_NBGRR`"
#     LINE_TAB2=$LINE_TAB2" & `percent_time $RR_TINGR $RR_TTIME`\\\\"   
#     
#     echo $LINE_TAB1 >> $TEMPTABFILE1
#     echo $LINE_TAB2 >> $TEMPTABFILE2
#     
# done


# # 
#Other polynomials
# DEGREES="64 128 191 256 383 512"
# DEGREES="64"
POLNAMES="Bernoulli Chebyshev1 Legendre Wilkinson"
# POLNAMES="Bernoulli"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE1
    echo $POLNAME >> $TEMPTABFILE2
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_risolate $NAME $POLNAME $DEG
    run_aNewDsc  $NAME $POLNAME $DEG
#     gen_and_run_ccluster $NAME
    
    stats_pol $NAME $DEG
#     LINE_TAB=" `format_numb $DEG $LENP` & "
#     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
    
#     stats_pol $NAME
done 
done
# 


POLNAME="RegularGrid"

echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2
for SIZ in $SIZEGRID; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$SIZ
    
    gen_with_deg $NAME $POLNAME $SIZ
    run_risolate $NAME $POLNAME $SIZ
    run_aNewDsc  $NAME $POLNAME $SIZ
#     gen_and_run_ccluster $NAME
    
    DEG=$(( 2*$SIZ+1 ))
    DEG=$(( $DEG*$DEG ))
    stats_pol $NAME $DEG
    
#     LINE_TAB=" `format_numb $SIZ $LENP` & "
#     LINE_TAB=$LINE_TAB"`stats_pol $REPNAME"/"$POLNAME"_"$SIZ`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
    
    
#     stats_pol $REPNAME"/"$POLNAME"_"$DEG
done
# 
# POLNAME="randomSparse"
# 
# #solve random sparse polynomials with ccluster
# echo $POLNAME >> $TEMPTABFILE
# for DEG in $DEGREES; do
# 
#     REPNAME=$REP"/"$POLNAME"_"$DEG
#     FILENAME=$POLNAME"_"$DEG
#     
#     if [ ! -d $REPNAME ]; then
# #         echo $REPNAME
#          mkdir -p $REPNAME
#     fi
#     
#     TSIZE=0
#     TDEPT=0
#     TTIME=0
#     EXPE_NBFAILS=0
#     EXPE_TTIME=0
#     EXPE_TSIZE=0
#     EXPE_TDEPT=0
#     EXPE_TINEVAL=0
#     EXPE_TINPOST=0
#     
# #     REP2="tableISSAC1sparse_IR43"
#     
#     for CURIND in `seq 1 $NBPOLS`; do
#         NAME=$REPNAME"/"$FILENAME"_"$CURIND
#         gen_with_deg_bs_nbterms $NAME $POLNAME $DEG $BITSIZE $NBTERMS
#         run_ccluster $NAME $POLNAME $DEG
# #         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".ccl" $NAME."ccl"
# #         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out" $NAME."out"
# #         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out0" $NAME."out0"
# #         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out1" $NAME."out1"
# #         cp $REP2"/randsparse_"$DEG"/randsparse_"$DEG"_"$CURIND".out2" $NAME."out2"
# 
#         stats_pol_rand $NAME
#     done
#     
#     LINE_TAB=" "$DEG" & "$TSIZE" & "$TDEPT" & `format_time $TTIME`"
#     LINE_TAB=$LINE_TAB" & "$EXPE_NBFAILS" & "$EXPE_TSIZE" & "$EXPE_TDEPT" & `format_time $EXPE_TTIME` & `percent_time $EXPE_TTIME $TTIME`"
#     LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINEVAL $EXPE_TTIME` & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
#     LINE_TAB=$LINE_TAB" & `percent_time $EXPE_TINPOST $EXPE_TTIME`"
#     LINE_TAB=$LINE_TAB"\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
#     
# #     LINE_TAB=" "$DEG" & "$NUMDISTESTS" & `percent_time $TTIMEINDIST $TTIMEINCCLU`"
# #     LINE_TAB=$LINE_TAB" & "$NUMTN0" & "$NUMFP0" & `ratio_time $TTIME0 $TTIMEINDIST`"
# #     LINE_TAB=$LINE_TAB" & "$NUMTN1" & "$NUMFP1" & `ratio_time $TTIME1 $TTIMEINDIST`"
# #     LINE_TAB=$LINE_TAB" & "$NUMTN2" & "$NUMFP2" & `ratio_time $TTIME2 $TTIMEINDIST`"
# #     LINE_TAB=$LINE_TAB"\\\\"
# # #     echo $LINE_TAB
# #     echo $LINE_TAB >> $TEMPTABFILE
# done
# 
# #Other polynomials
# # DEGREES="64 128 191 256 383 512"
# # DEGREES="64"
POLNAMES="Mignotte"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE1
    echo $POLNAME >> $TEMPTABFILE2
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg_bs $NAME $POLNAME $DEG $BITSIZE
    run_risolate $NAME $POLNAME $DEG
    run_aNewDsc  $NAME $POLNAME $DEG
#     gen_and_run_ccluster $NAME
    
    stats_pol $NAME $DEG
    
#     LINE_TAB=" `format_numb $DEG $LENP` & "
#     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
#     stats_pol $NAME
done 
done
# 
# #Procedural polynomials
# 
# POLNAMES="Mandelbrot"
# 
# for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE
# for DEG in $NBITT; do
#     
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG
#     
#     gen_with_deg $NAME $POLNAME $DEG
# #     run_ccluster $NAME $POLNAME $DEG
#     run_ccluster_mandelbrot $NAME $POLNAME $DEG
# #     gen_and_run_ccluster $NAME
#     
#     LINE_TAB=" "$DEG" & "
#     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
# #     stats_pol $NAME
# done 
# done


cat $TEMPTABFILE1
rm -f $TEMPTABFILE1
cat $TEMPTABFILE2
rm -f $TEMPTABFILE2
