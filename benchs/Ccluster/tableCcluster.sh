#!/bin/bash

usage()
{
   echo "Usage: ./tableCcluster.sh <options> <args>"
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
      --purgeCclu)
        PURGECCLU=1
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
   DEGREES="64 128 191 256"
#    DEGREES="128 191 256 391 512 791 1024"
#    DEGREES="128 191 256 391 512 791"
#    DEGREES="128 191 256 391 512"
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

if [ -z "$PURGECCLU" ]; then
   PURGECCLU=0
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
   SIZEGRID="5 6"
#    SIZEGRID="5 6 7 8 10 11 12 13 14"
#    SIZEGRID="10 11 12 13 14"
#    SIZEGRID="12 13 14 15"
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
CCLUSTER_CALL=$CCLUSTER_PATH"/bin/ccluster"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
GENRANDPOLFI_CALL=$CCLUSTER_PATH"/bin/genRandPolFile"
CCLUSTER_OPTS="-v 2"

MPSOLVE_CALL="mpsolve -au -Gi -o"$EPSILONMPS" -j1"
MPSOLVE_CALL_S="mpsolve -as -Gi -o"$EPSILONMPS" -j1"

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
            RES=$RES"_"
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
    
    if [ ! -e $NAME_IN2 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN2
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
    NAME_IN3=$NAME".dsc"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS
    fi
    
    if [ ! -e $NAME_IN2 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN2
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2 -b $BS
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
    NAME_IN2_MAX=$NAME"_"$NBPOLS".mpl"
    NAME_IN3_MAX=$NAME"_"$NBPOLS".dsc"
    
    if [ ! -e $NAME_IN_MAX ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 1 -b $BS -p $NBPOLS -l $LOC
    fi
    
    if [ ! -e $NAME_IN2_MAX ]; then
            echo  "Generating $NBPOLS files for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN2
            $GENRANDPOLFI_CALL $POLNAME $DEG -f 2 -b $BS -p $NBPOLS -l $LOC
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
    NAME_OUT=$NAME".out_cclu"
    
    if [ $PURGECCLU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CCLUSTER_CALL $NAME_IN -v 2 -m default"
            $CALL > $NAME_OUT
    fi
    
}

run_mpsolve()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    NAME_IN=$NAME".mpl"
    NAME_OUT=$NAME".out_mpl"
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Approximating complex roots for $POLNAME degree $DEG, with mpsolve, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $MPSOLVE_CALL_S $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

TEMPTABFILE1="temptab_ccluster1.txt"
touch $TEMPTABFILE1



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
    NAME_OUT=$NAME".out_cclu"
    NAME_OUTMPL=$NAME".out_mpl"
    DEG=$2
    
#     TCCLUSTER_L=$(grep "time:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDCCLUSTER_L=$(grep "tree depth:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     TSCCLUSTER_L=$(grep "tree size:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBCLUSTERS_L=$(grep "number of clusters:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
#     NBSOLUTION_L=$(grep "number of solutions:" $FILENAME_CCLUSTER_L_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TCCLUSTER_G=$( grep "time:"                $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G=$(grep "tree depth:"          $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G=$(grep "tree size:"           $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G=$(grep "number of clusters:"  $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G=$(grep "number of solutions:" $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE_S=$(grep "real"                   $NAME_OUTMPL| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB=$DEG" & ("$NBSOLUTION_G":"$NBCLUSTERS_G")&("$TDCCLUSTER_G":"$TSCCLUSTER_G")&`format_time $TCCLUSTER_G`"
    LINE_TAB=$LINE_TAB"&"$TMPSOLVE_S
    LINE_TAB=$LINE_TAB"\\\\\\hline"
    
#     LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $NSOLS $LENP`  & `format_numb $TSIZE $LENP` & `format_numb $TDEPT 2` & `format_time $TTIME`"
#     LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE $LENP` &`format_time $DSC_TTIME`\\\\"
    
    echo $LINE_TAB >> $TEMPTABFILE1
}

stats_pol_rand()
{
    NAME=$1
    NAME_OUT=$NAME".out_cclu"
    NAME_OUTMPL=$NAME".out_mpl"
    DEG=$2
    
    TCCLUSTER_G_T=$( grep "time:"                $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDCCLUSTER_G_T=$(grep "tree depth:"          $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    TSCCLUSTER_G_T=$(grep "tree size:"           $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBCLUSTERS_G_T=$(grep "number of clusters:"  $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    NBSOLUTION_G_T=$(grep "number of solutions:" $NAME_OUT| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
    
    TMPSOLVE_S_T=$(grep "real"                   $NAME_OUTMPL| cut -f2 -d'l' | tr -d ' ')
    
     TCCLUSTER_G=`echo  $TCCLUSTER_G+$TCCLUSTER_G_T |bc -l`
    TDCCLUSTER_G=`echo $TDCCLUSTER_G+$TDCCLUSTER_G_T|bc -l`
    TSCCLUSTER_G=`echo $TSCCLUSTER_G+$TSCCLUSTER_G_T|bc -l`
    NBCLUSTERS_G=`echo $NBCLUSTERS_G+$NBCLUSTERS_G_T|bc -l`
    NBSOLUTION_G=`echo $NBSOLUTION_G+$NBSOLUTION_G_T|bc -l`
      TMPSOLVE_S=`echo   $TMPSOLVE_S+$TMPSOLVE_S_T  |bc -l`
    
    
#     LINE_TAB1="`format_numb $DEG $LENP` & `format_numb $NSOLS_T $LENP` & `format_numb $TSIZE_T $LENP` & `format_numb $TDEPT_T 2` & `format_time $TTIME_T`"
#     LINE_TAB1=$LINE_TAB1" & `format_numb $DSC_TSIZE_T $LENP` &`format_time $DSC_TTIME_T`\\\\"
    
#     echo $LINE_TAB1
}

REP="tabCcluster"

if [ -d "$REP" ]; then
  if [ $PURGE -eq 1 ]; then
    rm -rf $REP
    mkdir $REP
  fi
else
  mkdir $REP
fi

# POLNAME="randomDense"
# 
# #solve random polynomials with ccluster
# echo $POLNAME >> $TEMPTABFILE1
# 
# for DEG in $DEGREES; do
#                                 
# #     REPNAME=$REP"/"$POLNAME"_"$DEG
# #     FILENAME=$POLNAME"_"$DEG
#     
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITSIZE
#     
# #     if [ ! -d $REPNAME ]; then
# # #         echo $REPNAME
# #          mkdir -p $REPNAME
# #     fi
#     
#      TCCLUSTER_G=0
#     TDCCLUSTER_G=0
#     TSCCLUSTER_G=0
#     NBCLUSTERS_G=0
#     NBSOLUTION_G=0
#       TMPSOLVE_S=0
#     
#     genRand_with_deg_bs $NAME $POLNAME $DEG $BITSIZE $NBPOLS $REPNAME
#     for CURIND in `seq 1 $NBPOLS`; do
#         NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BITSIZE"_"$CURIND
#         run_ccluster $NAME $POLNAME $DEG
#         run_mpsolve  $NAME $POLNAME $DEG
#         stats_pol_rand $NAME $DEG
#     done
#     
#      TCCLUSTER_G=`echo     $TCCLUSTER_G     /$NBPOLS     |bc -l`
#     TDCCLUSTER_G=`echo     $TDCCLUSTER_G    /$NBPOLS     |bc -l`
#     TSCCLUSTER_G=`echo     $TSCCLUSTER_G    /$NBPOLS     |bc -l`
#     NBCLUSTERS_G=`echo     $NBCLUSTERS_G    /$NBPOLS     |bc -l`
#     NBSOLUTION_G=`echo     $NBSOLUTION_G    /$NBPOLS     |bc -l`
#       TMPSOLVE_S=`echo     $TMPSOLVE_S      /$NBPOLS     |bc -l`
#     
#     
#     LINE_TAB=$DEG" & ("`format_time $NBSOLUTION_G`":"`format_time $NBCLUSTERS_G`")&("`format_time $TDCCLUSTER_G`":"`format_time $TSCCLUSTER_G`")&`format_time $TCCLUSTER_G`"
#     LINE_TAB=$LINE_TAB"&"`format_time $TMPSOLVE_S`
#     LINE_TAB=$LINE_TAB"\\\\\\hline"
#     
#     echo $LINE_TAB >> $TEMPTABFILE1
#     
# done
# # 
# POLNAME="randomDense"
# 
# #solve random polynomials with ccluster
# echo $POLNAME >> $TEMPTABFILE1
# DEGF=256
# # BITSIZES="50 100 500 1000 5000"
# BITSIZES="50 100 500"
# for BIT in $BITSIZES; do
# 
#     REPNAME=$REP
#     NAME=$REPNAME"/"$POLNAME"_"$DEGF"_"$BIT
#     
# #     if [ ! -d $REPNAME ]; then
# # #         echo $REPNAME
# #          mkdir -p $REPNAME
# #     fi
#     
#      TCCLUSTER_G=0
#     TDCCLUSTER_G=0
#     TSCCLUSTER_G=0
#     NBCLUSTERS_G=0
#     NBSOLUTION_G=0
#       TMPSOLVE_S=0
#     
#     genRand_with_deg_bs $NAME $POLNAME $DEGF $BIT $NBPOLS $REPNAME
#     for CURIND in `seq 1 $NBPOLS`; do
#         NAME=$REPNAME"/"$POLNAME"_"$DEGF"_"$BIT"_"$CURIND
#         run_ccluster $NAME $POLNAME $DEGF
#         run_mpsolve  $NAME $POLNAME $DEGF
#         stats_pol_rand $NAME $DEGF
#     done
#     
#      TCCLUSTER_G=`echo     $TCCLUSTER_G     /$NBPOLS     |bc -l`
#     TDCCLUSTER_G=`echo     $TDCCLUSTER_G    /$NBPOLS     |bc -l`
#     TSCCLUSTER_G=`echo     $TSCCLUSTER_G    /$NBPOLS     |bc -l`
#     NBCLUSTERS_G=`echo     $NBCLUSTERS_G    /$NBPOLS     |bc -l`
#     NBSOLUTION_G=`echo     $NBSOLUTION_G    /$NBPOLS     |bc -l`
#       TMPSOLVE_S=`echo     $TMPSOLVE_S      /$NBPOLS     |bc -l`
#     
#     
#     LINE_TAB=$DEG" & ("`format_time $NBSOLUTION_G`":"`format_time $NBCLUSTERS_G`")&("`format_time $TDCCLUSTER_G`":"`format_time $TSCCLUSTER_G`")&`format_time $TCCLUSTER_G`"
#     LINE_TAB=$LINE_TAB"&"`format_time $TMPSOLVE_S`
#     LINE_TAB=$LINE_TAB"\\\\\\hline"
#     
#     echo $LINE_TAB >> $TEMPTABFILE1
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
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_ccluster $NAME $POLNAME $DEG
    run_mpsolve  $NAME $POLNAME $DEG
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
for SIZ in $SIZEGRID; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$SIZ
    
    gen_with_deg $NAME $POLNAME $SIZ
    run_ccluster $NAME $POLNAME $SIZ
    run_mpsolve  $NAME $POLNAME $SIZ
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
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg_bs $NAME $POLNAME $DEG $BITSIZE
    run_ccluster $NAME $POLNAME $DEG
    run_mpsolve  $NAME $POLNAME $DEG
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
