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
      --purgeRiso)
        PURGERISO=1
        ;;
      --purgeRR)
        PURGERR=1
        ;;
      --purgeDSC)
        PURGEDSC=1
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
#    DEGREES="64 128 191 256 391 512 791"
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

if [ -z "$PURGERR" ]; then
   PURGERR=0
fi

if [ -z "$PURGEDSC" ]; then
   PURGEDSC=0
fi

if [ -z "$GENERATE" ]; then
   GENERATE=0
fi

if [ -z "$NBPOLS" ]; then
#    NBPOLS=20
   NBPOLS=10
#    NBPOLS=1
fi

if [ -z "$NBITT" ]; then
   NBITT="7 8 9"
fi

if [ -z "$SIZEGRID" ]; then
#    SIZEGRID="5 6"
#    SIZEGRID="5 6 7 8 10 11 12 13 14"
#    SIZEGRID="6 7 8 9 10 11 12 13 14"
   SIZEGRID="6 8 10 12 14 16"
#    SIZEGRID="6 8"
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
    NAME_IN=$NAME"_nbp.ccl"
    NAME_IN2=$NAME"_nbp.mpl"
    NAME_IN3=$NAME"_nbp.dsc"
    NAME_IN_MAX=$NAME"_"$NBPOLS".ccl"
    NAME_IN3_MAX=$NAME"_"$NBPOLS".dsc"
    
    if [ ! -e $NAME_IN_MAX ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 1 -b $BS -p $NBPOLS -l $LOC
    fi
    
    if [ ! -e $NAME_IN3_MAX ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN3
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 3 -b $BS -p $NBPOLS -l $LOC
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

run_risolate()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    
    if [ $PURGERISO -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m onlySubd"
            $CALL > $NAME_OUT
    fi
    
    if [ $PURGERR -eq 1 ]; then
        rm -f $NAME_OUTRR
    fi
  
    if [ ! -e $NAME_OUTRR ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, global, root radii, output in " $NAME_OUTRR
#             ./ccluster $NAME_IN "global" $EPSILONCCL "test" 2 > $NAME_OUT0
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUTRR
    fi
    
}

run_aNewDsc()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    SQUAREFREE=$4
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
            echo  "Isolating real roots for $POLNAME degree $DEG, with aNewDsc, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $CALL $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

TEMPTABFILE1="temptab_risolate1.txt"
touch $TEMPTABFILE1
TEMPTABFILE2="temptab_risolate2.txt"
touch $TEMPTABFILE2
TEMPTABFILE3="temptab_risolate3.txt"
touch $TEMPTABFILE3


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
    
#     BITSI=$(grep "bitsize of input polynomial:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    BITSI=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NSOLS=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCOT=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate
#     NBTZT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBTST=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate RR
    RR_NSOLS=$(grep "number of real roots:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBEXT=$(grep "total number DT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBCOT=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
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
    
    LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_numb $NSOLS $LENP`  & `format_numb $RR_NSOLS $LENP` & `format_numb $TSIZE $LENP` & `format_numb $TDEPT 2` & `format_time $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $RR_TSIZE $LENP` & `format_numb $RR_TDEPT 2` & `format_time $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
    LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`\\\\"
    
    LINE_TAB2="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_time $RR_TTIME` & `percent_time $RR_TINRR $RR_TTIME`"
    LINE_TAB2=$LINE_TAB2" & `format_numb $RR_PPREC $LENP` & `format_numb $RR_PRECN $LENP`"
    LINE_TAB2=$LINE_TAB2" & `format_numb $RR_NBGRA $LENP` & `format_numb $RR_NBGRR $LENP`"
    LINE_TAB2=$LINE_TAB2" & `percent_time $RR_TINGR $RR_TTIME`\\\\"  
    
    LINE_TAB3="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_numb $NSOLS $LENP`"
    LINE_TAB3=$LINE_TAB3" & `format_time $TTIME`     & `format_numb $NBEXT $LENP`,`format_numb $NBCOT $LENP`"
    LINE_TAB3=$LINE_TAB3" & `format_time $RR_TTIME`  & `format_numb $RR_NBEXT $LENP`,`format_numb $RR_NBCOT $LENP`"
#     LINE_TAB3=$LINE_TAB3" & `format_numb $RR_PRECN $LENP` "
    LINE_TAB3=$LINE_TAB3" & `percent_time $RR_TINRR $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
    LINE_TAB3=$LINE_TAB3" & `format_time $DSC_TTIME`\\\\"
    
#     echo $LINE_TAB1
    echo $LINE_TAB1 >> $TEMPTABFILE1
    echo $LINE_TAB2 >> $TEMPTABFILE2
    echo $LINE_TAB3 >> $TEMPTABFILE3
}

stats_pol_rand()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    NAME_OUTDSC=$NAME".out_dsc"
#     DEG=$2
    
    BITSI_T=$(grep "bitsize of input polynomial:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NSOLS_T=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TSIZE_T=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_T=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCOT_T=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_T=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
        # details on risolate #
#     NBTZT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBTST_T=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    # details on risolate RR #
    RR_NSOLS_T=$(grep "number of real roots:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME_T=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE_T=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT_T=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBEXT_T=$(grep "total number DT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBCOT_T=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
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
    
#     BITSI=`echo $BITSI+$BITSI_T|bc -l`
    NSOLS=`echo $NSOLS+$NSOLS_T|bc -l`
    TSIZE=`echo $TSIZE+$TSIZE_T|bc -l`
    TDEPT=`echo $TDEPT+$TDEPT_T|bc -l`
    NBEXT=`echo $NBEXT+$NBEXT_T|bc -l`
    NBCOT=`echo $NBCOT+$NBCOT_T|bc -l`
    TTIME=`echo $TTIME+$TTIME_T|bc -l`
    TTIME_SQ=`echo $TTIME_SQ+$TTIME_T^2|bc -l`
#     NBTZT=`echo $NBTZT+$NBTZT_T|bc -l`
#     NBTST=`echo $NBTST+$NBTST_T|bc -l`
    RR_NSOLS=`echo $RR_NSOLS +$RR_NSOLS_T|bc -l`
    RR_TTIME=`echo $RR_TTIME +$RR_TTIME_T|bc -l`
    RR_TTIME_SQ=`echo $RR_TTIME_SQ+$RR_TTIME_T^2|bc -l`
    RR_TSIZE=`echo $RR_TSIZE +$RR_TSIZE_T|bc -l`
    RR_TDEPT=`echo $RR_TDEPT +$RR_TDEPT_T|bc -l`
    RR_NBCOT=`echo $RR_NBCOT +$RR_NBCOT_T|bc -l`
    RR_NBEXT=`echo $RR_NBEXT +$RR_NBEXT_T|bc -l`
    RR_PRECN=`echo $RR_PRECN +$RR_PRECN_T|bc -l`
    RR_PPREC=`echo $RR_PPREC +$RR_PPREC_T|bc -l`
    RR_NBGRA=`echo $RR_NBGRA +$RR_NBGRA_T|bc -l`
    RR_NBGRR=`echo $RR_NBGRR +$RR_NBGRR_T|bc -l`
    RR_TINGR=`echo $RR_TINGR +$RR_TINGR_T|bc -l`
    RR_TINRR=`echo $RR_TINRR +$RR_TINRR_T|bc -l`
    DSC_NSOLS=`echo $DSC_NSOLS+$DSC_NSOLS_T|bc -l`
    DSC_TSIZE=`echo $DSC_TSIZE+$DSC_TSIZE_T|bc -l`
    DSC_TTIME=`echo $DSC_TTIME+$DSC_TTIME_T|bc -l`
    DSC_TTIME_SQ=`echo $DSC_TTIME_SQ+$DSC_TTIME_T^2|bc -l`
#     echo $DSC_TTIME_SQ
    
#     LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $BITSI_T $LENP` & `format_numb $NSOLS_T $LENP` & `format_numb $RR_NSOLS_T $LENP` & `format_numb $TSIZE_T $LENP` & `format_numb $TDEPT_T 2` & `format_time $TTIME_T`"
#     LINE_TAB1=$LINE_TAB1" & `format_numb $RR_TSIZE_T $LENP` & `format_numb $RR_TDEPT_T 2` & `format_time $RR_TTIME_T` & `percent_time $RR_TTIME_T $TTIME_T`"
#     LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE_T $LENP` &`format_time $DSC_TTIME_T`\\\\"
#     
#     echo $LINE_TAB1
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
# 
# #solve random polynomials with risolate
echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2
echo $POLNAME >> $TEMPTABFILE3
# 
# DEGREES="256 391 512"
DEGREES="16 32 64 256 391 512"
# BITSIZES="16 8192 16384 32768"
BITSIZES="16 8192 16384 32768 65536"
for DEG in $DEGREES; do
    for BIT in $BITSIZES; do
            REPNAME=$REP
            NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT
            
            NSOLS=0
            TSIZE=0
            TDEPT=0
            NBEXT=0
            NBCOT=0
            TTIME=0
            TTIME_SQ=0
            RR_NSOLS=0
            RR_TTIME=0
            RR_TTIME_SQ=0
            RR_TSIZE=0
            RR_TDEPT=0
            RR_NBEXT=0
            RR_NBCOT=0
            RR_PRECN=0
            RR_PPREC=0
            RR_NBGRA=0
            RR_NBGRR=0
            RR_TINGR=0
            RR_TINRR=0
            DSC_NSOLS=0
            DSC_TSIZE=0
            DSC_TTIME=0
            DSC_TTIME_SQ=0
            
            genRand_with_deg_bs $NAME $POLNAME $DEG $BIT $NBPOLS $REPNAME
            for CURIND in `seq 1 $NBPOLS`; do
                NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT"_"$CURIND
                run_risolate $NAME $POLNAME $DEG
        #         run_mpsolve  $NAME $POLNAME $DEG
                run_aNewDsc  $NAME $POLNAME $DEG 0
                stats_pol_rand $NAME
            done
            
            NSOLS=`echo     $NSOLS    /$NBPOLS     |bc -l`
            TSIZE=`echo     $TSIZE    /$NBPOLS     |bc -l`
            TDEPT=`echo     $TDEPT    /$NBPOLS     |bc -l`
            NBEXT=`echo     $NBEXT    /$NBPOLS     |bc -l`
            NBCOT=`echo     $NBCOT    /$NBPOLS     |bc -l`
            TTIME=`echo     $TTIME    /$NBPOLS     |bc -l`
            TTIME_SQ=`echo "sqrt("$TTIME_SQ"/"$NBPOLS "-" $TTIME"^2)"   |bc -l`
            RR_NSOLS=`echo  $RR_NSOLS /$NBPOLS     |bc -l`
            RR_TTIME=`echo  $RR_TTIME /$NBPOLS     |bc -l`
            RR_TTIME_SQ=`echo "sqrt("$RR_TTIME_SQ"/"$NBPOLS "-" $RR_TTIME"^2)"   |bc -l`
            RR_TSIZE=`echo  $RR_TSIZE /$NBPOLS     |bc -l`
            RR_TDEPT=`echo  $RR_TDEPT /$NBPOLS     |bc -l`
            RR_NBEXT=`echo  $RR_NBEXT /$NBPOLS     |bc -l`
            RR_NBCOT=`echo  $RR_NBCOT /$NBPOLS     |bc -l`
            RR_PRECN=`echo  $RR_PRECN /$NBPOLS     |bc -l`
            RR_PPREC=`echo  $RR_PPREC /$NBPOLS     |bc -l`
            RR_NBGRA=`echo  $RR_NBGRA /$NBPOLS     |bc -l`
            RR_NBGRR=`echo  $RR_NBGRR /$NBPOLS     |bc -l`
            RR_TINGR=`echo  $RR_TINGR /$NBPOLS     |bc -l`
            RR_TINRR=`echo  $RR_TINRR /$NBPOLS     |bc -l`
            DSC_NSOLS=`echo $DSC_NSOLS/$NBPOLS     |bc -l`
            DSC_TSIZE=`echo $DSC_TSIZE/$NBPOLS     |bc -l`
            DSC_TTIME=`echo $DSC_TTIME/$NBPOLS     |bc -l`
            DSC_TTIME_SQ=`echo "sqrt("$DSC_TTIME_SQ"/"$NBPOLS "-" $DSC_TTIME"^2)"   |bc -l`
            
            echo "$POLNAME & $DEG & $BIT & `format_time $NSOLS` & `format_time $RR_NSOLS` & `format_time $DSC_NSOLS` "
            
            LINE_TAB3="`format_numb $DEG $LENP` & `format_numb $BIT $LENP` & `format_time $NSOLS`"
            LINE_TAB3=$LINE_TAB3" & `format_time $TTIME` (`format_time $TTIME_SQ`) & `format_time $NBEXT`,`format_time $NBCOT`"
            LINE_TAB3=$LINE_TAB3" & `format_time $RR_TTIME` (`format_time $RR_TTIME_SQ`)  & `format_time $RR_NBEXT`,`format_time $RR_NBCOT`"
        #     LINE_TAB3=$LINE_TAB3" & `format_time $RR_PRECN` "
            LINE_TAB3=$LINE_TAB3" & `percent_time $RR_TINRR $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
            LINE_TAB3=$LINE_TAB3" & `format_time $DSC_TTIME` (`format_time $DSC_TTIME_SQ`)\\\\"
            
#             echo $LINE_TAB1 >> $TEMPTABFILE1
#             echo $LINE_TAB2 >> $TEMPTABFILE2
            echo $LINE_TAB3 >> $TEMPTABFILE3
    done
done

#Other polynomials
# DEGREES="64 128 191 256 383 512"
# DEGREES="64"
# DEGREES="128 191 256 391 512 791 1024"
# # POLNAMES="Bernoulli Chebyshev1 Legendre Wilkinson"
# POLNAMES="Bernoulli Wilkinson"
# # POLNAMES="Bernoulli Chebyshev1 Legendre"

# DEGREES="128 196 256"
DEGREES="256 391 512 791 1024"
# DEGREES="256 391"
# POLNAMES="Bernoulli Chebyshev1 Legendre Wilkinson"
POLNAMES="Bernoulli Wilkinson"
# POLNAMES="Bernoulli Chebyshev1 Legendre"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE1
    echo $POLNAME >> $TEMPTABFILE2
    echo $POLNAME >> $TEMPTABFILE3
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_risolate $NAME $POLNAME $DEG
    run_aNewDsc  $NAME $POLNAME $DEG 1
#     gen_and_run_ccluster $NAME
    
    stats_pol $NAME $DEG
    echo "$POLNAME & $DEG & $BIT & `format_numb $NSOLS` & `format_numb $RR_NSOLS` & `format_numb $DSC_NSOLS` "
#     LINE_TAB=" `format_numb $DEG $LENP` & "
#     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
    
#     stats_pol $NAME
done 
done
# # 
# 
# 
POLNAME="RegularGrid"
SIZEGRID="8 10 12 14 16"
# SIZEGRID="6 8"

echo $POLNAME >> $TEMPTABFILE1
echo $POLNAME >> $TEMPTABFILE2
echo $POLNAME >> $TEMPTABFILE3
for SIZ in $SIZEGRID; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$SIZ
    
    gen_with_deg $NAME $POLNAME $SIZ
    run_risolate $NAME $POLNAME $SIZ
    run_aNewDsc  $NAME $POLNAME $SIZ 1
#     gen_and_run_ccluster $NAME
    
    DEG=$(( 2*$SIZ+1 ))
    DEG=$(( $DEG*$DEG ))
    stats_pol $NAME $DEG
    echo "$POLNAME & $DEG & $BIT & `format_numb $NSOLS` & `format_numb $RR_NSOLS` & `format_numb $DSC_NSOLS` "
#     LINE_TAB=" `format_numb $SIZ $LENP` & "
#     LINE_TAB=$LINE_TAB"`stats_pol $REPNAME"/"$POLNAME"_"$SIZ`\\\\"
#     echo $LINE_TAB >> $TEMPTABFILE
    
    
#     stats_pol $REPNAME"/"$POLNAME"_"$DEG
done

# DEGREES="128 191 256 391 512 791 1024"
# POLNAMES="Mignotte"
# 
# for POLNAME in $POLNAMES; do
#     echo $POLNAME >> $TEMPTABFILE1
#     echo $POLNAME >> $TEMPTABFILE2
#     echo $POLNAME >> $TEMPTABFILE3
# for DEG in $DEGREES; do
#     
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG
#     
#     gen_with_deg_bs $NAME $POLNAME $DEG $BITSIZE
#     run_risolate $NAME $POLNAME $DEG
#     run_aNewDsc  $NAME $POLNAME $DEG 1
# #     gen_and_run_ccluster $NAME
#     
#     stats_pol $NAME $DEG
#     echo "$POLNAME & $DEG & $BIT & `format_numb $NSOLS` & `format_numb $RR_NSOLS` & `format_numb $DSC_NSOLS` "
# #     LINE_TAB=" `format_numb $DEG $LENP` & "
# #     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
# #     echo $LINE_TAB >> $TEMPTABFILE
# #     stats_pol $NAME
# done 
# done
# 
# POLNAMES="Mignotte"
# # with increasing bit-size
# 
# # DEGF=512
# DEGF=256
# BITSIZES="128 256 512 1024 2048"
# # BITSIZES="128 256 512 1024 2048 4096 8192 16384 32768 65536"
# for BIT in $BITSIZES; do
#     
#     DEG=$DEGF
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT
#     
#     gen_with_deg_bs $NAME $POLNAME $DEG $BIT
#     run_risolate $NAME $POLNAME $DEG
#     run_aNewDsc  $NAME $POLNAME $DEG 1
# #     gen_and_run_ccluster $NAME
#     
#     stats_pol $NAME $DEG
#     echo "$POLNAME & $DEG & $BIT & `format_numb $NSOLS` & `format_numb $RR_NSOLS` & `format_numb $DSC_NSOLS` "
# #     LINE_TAB=" `format_numb $DEG $LENP` & "
# #     LINE_TAB=$LINE_TAB"`stats_pol $NAME`\\\\"
# #     echo $LINE_TAB >> $TEMPTABFILE
# #     stats_pol $NAME
# done

# 
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


# cat $TEMPTABFILE1
rm -f $TEMPTABFILE1
# cat $TEMPTABFILE2
rm -f $TEMPTABFILE2
cat $TEMPTABFILE3
rm -f $TEMPTABFILE3
