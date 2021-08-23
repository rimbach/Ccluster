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
#             RES=$RES"_"
            RES=$RES""
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
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1
    fi
    
}

gen_with_deg_bs(){

    NAME=$1
    POLNAME=$2
    DEG=$3
    BS=$4
    NAME_IN=$NAME".ccl"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS
    fi
    
}

gen_with_c_a_k(){

    NAME=$1
    POLNAME=$2
    C=$3
    A=$4
    K=$5
    NAME_IN=$NAME".ccl"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME c=$C, a=$A, k=$K, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $K $NAME_IN -f 1 -c $C -a $A
    fi
    
}

run_ccluster()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    EPS=$4
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_ccl"
    
    if [ $PURGECCL -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering complex roots for $POLNAME degree $DEG, global, default, output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CCLUSTER_CALL $NAME_IN -e $EPS $CCLUSTER_OPTS"
#             echo $CALL
            $CALL > $NAME_OUT
    fi
    
}

run_cauchy()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    EPS=$4
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_cau"
    
    if [ $PURGECAU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering complex roots for $POLNAME degree $DEG, with cauchy output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CAUCHY_CALL $NAME_IN -e $EPS $CAUCHY_OPTS"
#             echo $CALL
            $CALL > $NAME_OUT
    fi
    
}

run_cauchy_comp()
{
    NAME=$1
    POLNAME=$2
    EPS=$3
    NAME_IN=$NAME".ccl"
    NAME_OUT_WO=$NAME".out_cau_wo"
    NAME_OUT_WI=$NAME".out_cau"
    
    if [ $PURGECAUWO -eq 1 ]; then
        rm -f $NAME_OUT_WO
    fi
    
    if [ ! -e $NAME_OUT_WO ]; then
            echo  "Clustering complex roots for $NAME with cauchy  without compression output in " $NAME_OUT_WO
            CALL="$CAUCHY_CALL $NAME_IN -e $EPS -m \"C1\" $CAUCHY_OPTS"
            $CALL > $NAME_OUT_WO
    fi
    
    if [ $PURGECAUWI -eq 1 ]; then
        rm -f $NAME_OUT_WI
    fi
    
    if [ ! -e $NAME_OUT_WI ]; then
            echo  "Clustering complex roots for $NAME with cauchy  with    compression output in " $NAME_OUT_WI
#             CALL="$CAUCHY_CALL $NAME_IN -e $EPS -m \"C2\" $CAUCHY_OPTS"
            CALL="$CAUCHY_CALL $NAME_IN -e $EPS $CAUCHY_OPTS"
#             echo $CALL
            $CALL > $NAME_OUT_WI
    fi
}

run_cauchy_mignotte_comp()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    EPS=$5
    NAME_IN=$NAME".ccl"
    NAME_OUT_WO=$NAME".out_cau_wo"
    NAME_OUT_WI=$NAME".out_cau"
    
    if [ $PURGECAUWO -eq 1 ]; then
        rm -f $NAME_OUT_WO
    fi
    
    if [ ! -e $NAME_OUT_WO ]; then
            echo  "Clustering complex roots for $NAME with cauchy  without compression output in " $NAME_OUT_WO
            CALL="$CAUCHY_MIGNOTTE_CALL $DEG $BIT -e $EPS -m \"C1\" $CAUCHY_OPTS"
            $CALL > $NAME_OUT_WO
    fi
    
    if [ $PURGECAUWI -eq 1 ]; then
        rm -f $NAME_OUT_WI
    fi
    
    if [ ! -e $NAME_OUT_WI ]; then
            echo  "Clustering complex roots for $NAME with cauchy  with    compression output in " $NAME_OUT_WI
#             CALL="$CAUCHY_CALL $NAME_IN -e $EPS -m \"C2\" $CAUCHY_OPTS"
            CALL="$CAUCHY_MIGNOTTE_CALL $DEG $BIT -e $EPS $CAUCHY_OPTS"
#             echo $CALL
            $CALL > $NAME_OUT_WI
    fi
}

run_cauchy_mandelbrot()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    EPS=$4
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_cau"
    
    if [ $PURGECAU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering complex roots for $POLNAME degree $DEG, with cauchy output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CAUCHY_MANDELBROT_CALL $DEG -e $EPS $CAUCHY_OPTS"
            $CALL > $NAME_OUT
    fi
    
}

run_cauchy_runnels()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    EPS=$4
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_cau"
    
    if [ $PURGECAU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering complex roots for $POLNAME degree $DEG, with cauchy output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CAUCHY_RUNNELS_CALL $DEG -e $EPS $CAUCHY_OPTS"
            $CALL > $NAME_OUT
    fi
    
}

run_cauchy_mignotte()
{
    NAME=$1
    POLNAME=$2
    DEG=$3
    BIT=$4
    EPS=$5
    NAME_IN=$NAME".ccl"
    NAME_OUT=$NAME".out_cau"
    
    if [ $PURGECAU -eq 1 ]; then
        rm -f $NAME_OUT
    fi
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Clustering complex roots for $POLNAME degree $DEG, with cauchy output in " $NAME_OUT
#             ./ccluster $NAME_IN "global" $EPSILONCCL "default" 2 > $NAME_OUT
            CALL="$CAUCHY_MIGNOTTE_CALL $DEG $BIT -e $EPS $CAUCHY_OPTS"
            $CALL > $NAME_OUT
    fi
    
}

stats_pol()
{
    NAME=$1
    NAME_OUTCCL=$NAME".out_ccl"
    NAME_OUTCAU=$NAME".out_cau"
    DEG=$2
    
    NBSOLS=$(grep "number of solutions"              $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE_CCL=$(grep "tree size:"                    $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_CCL=$(grep "tree depth:"                   $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_CCL=$(grep "total number DT:"              $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_CCL=$(grep "total time:"                   $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TSIZE_CAU=$(grep "tree size:"                    $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_CAU=$(grep "tree depth:"                   $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_CAU=$(grep "total number ET:"              $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_CAU=$(grep "total time:"                   $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TTIME_CCL=`format_time $TTIME_CCL`
    TTIME_CAU=`format_time $TTIME_CAU`
#     echo $TTIME_CCL
    COLORCCL="\\coblue{"
    COLORCAU="\\cored{"
    TEST=`echo "$TTIME_CCL > $TTIME_CAU" | bc -l`
    if [ $TEST -eq 1 ]; then
        COLORCCL="\\cored{"
        COLORCAU="\\coblue{"
    fi
    
    K=$DEG
    D=$NBSOLS
    TEST=`echo "$K == $D" | bc -l`
    if [ $TEST -eq 1 ]; then
        K=" "
    fi
    
    LINE_TAB="$K & $D"
    LINE_TAB=$LINE_TAB" & `format_numb $TSIZE_CCL $LENP` & `format_numb $TDEPT_CCL 2`"
    LINE_TAB=$LINE_TAB" & $COLORCCL$TTIME_CCL}           & `format_numb $NBEXT_CCL $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $TSIZE_CAU $LENP` & `format_numb $TDEPT_CAU 2`"
    LINE_TAB=$LINE_TAB" & $COLORCAU$TTIME_CAU}       & `format_numb $NBEXT_CAU $LENP`"
    LINE_TAB=$LINE_TAB"\\\\"
    
#     echo $LINE_TAB
    echo $LINE_TAB >> $TEMPTABFILE
}

stats_pol_comp()
{
    NAME=$1
    NAME_OUT_WO=$NAME".out_cau_wo"
    NAME_OUT_WI=$NAME".out_cau"
    K=$2
    
    NBSOLS=$(grep "number of solutions"              $NAME_OUT_WO| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_WO=$(grep "total number ET:"              $NAME_OUT_WO| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_WO=$(grep "total time:"                   $NAME_OUT_WO| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBNIT_WO=$(grep "total number NE:"              $NAME_OUT_WO| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBNFA_WO=$(grep "number of fails:"              $NAME_OUT_WO| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    NBEXT_WI=$(grep "total number ET:"              $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_WI=$(grep "total time:"                   $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBNIT_WI=$(grep "total number NE:"              $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBNFA_WI=$(grep "number of fails:"              $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    COMPR_WI=$(grep "time spent in certi. RR algo:"  $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    COMPT_WI=$(grep "total time spent in compression:" $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCOM_WI=$(grep "total number for clus of >1 root:" $NAME_OUT_WI| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    D=$NBSOLS
    TEST=`echo "$K == $D" | bc -l`
    if [ $TEST -eq 1 ]; then
        K=" "
    fi
    
    TTIME_WO=`format_time $TTIME_WO`
    TTIME_WI=`format_time $TTIME_WI`
#     echo $TTIME_CCL
    COLORWO="\\coblue{"
    COLORWI="\\cored{"
    TEST=`echo "$TTIME_WO > $TTIME_WI" | bc -l`
    if [ $TEST -eq 1 ]; then
        COLORWO="\\cored{"
        COLORWI="\\coblue{"
    fi
    
    LINE_TAB="$K & $D"
#     LINE_TAB=$LINE_TAB" & `format_numb $TSIZE_CCL $LENP` & `format_numb $TDEPT_CCL 2`"
    LINE_TAB=$LINE_TAB" & $COLORWO$TTIME_WO}           & `format_numb $NBEXT_WO $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $NBNIT_WO $LENP` & `format_numb $NBNFA_WO $LENP`"
    LINE_TAB=$LINE_TAB" & $COLORWI$TTIME_WI}       & `format_numb $NBEXT_WI $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $NBNIT_WI $LENP` & `format_numb $NBNFA_WI $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $NBCOM_WI $LENP`"
    LINE_TAB=$LINE_TAB" & `format_time $COMPR_WI`  &  `format_time $COMPT_WI`"
    LINE_TAB=$LINE_TAB"\\\\"
    
#     echo $LINE_TAB
    echo $LINE_TAB >> $TEMPTABFILE
}

stats_pol_mig()
{
    NAME=$1
    NAME_OUTCCL=$NAME".out_ccl"
    NAME_OUTCAU=$NAME".out_cau"
    DEG=$2
    EPS=$3
    
    NBSOLS=$(grep "number of solutions"              $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE_CCL=$(grep "tree size:"                    $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_CCL=$(grep "tree depth:"                   $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_CCL=$(grep "total number DT:"              $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_CCL=$(grep "total time:"                   $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCLU_CCL=$(grep "number of clusters"            $NAME_OUTCCL| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TSIZE_CAU=$(grep "tree size:"                    $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_CAU=$(grep "tree depth:"                   $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_CAU=$(grep "total number ET:"              $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_CAU=$(grep "total time:"                   $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCLU_CAU=$(grep "number of clusters"            $NAME_OUTCAU| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    LINE_TAB="`format_numb $DEG $LENP` & `format_numb $NBSOLS $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $EPS $LENP` & `format_numb $NBCLU_CAU $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $TSIZE_CCL $LENP` & `format_numb $TDEPT_CCL 2`"
    LINE_TAB=$LINE_TAB" & `format_time $TTIME_CCL`       & `format_numb $NBEXT_CCL $LENP`"
    LINE_TAB=$LINE_TAB" & `format_numb $TSIZE_CAU $LENP` & `format_numb $TDEPT_CAU 2`"
    LINE_TAB=$LINE_TAB" & `format_time $TTIME_CAU`       & `format_numb $NBEXT_CAU $LENP`"
    LINE_TAB=$LINE_TAB"\\\\"
    
#     echo $LINE_TAB
    echo $LINE_TAB >> $TEMPTABFILE
}

# stats_pol_rand()
# {
#     NAME=$1
#     NAME_OUT=$NAME".out_ccl"
#     NAME_OUTRR=$NAME".out_ccl_rr"
#     NAME_OUTMPSOLVE=$NAME".out_mpl"
# #     NAME_OUTDSC=$NAME".out_dsc"
# #     DEG=$2
#     
# #     BITSI_T=$(grep "bitsize of input polynomial:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NSOLS_T=$(grep "number of solutions:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TSIZE_T=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TDEPT_T=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     NBEXT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     TTIME_T=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#         # details on risolate #
# #     NBTZT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# #     NBTST_T=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     # details on risolate RR #
#     RR_NSOLS_T=$(grep "number of solutions:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TTIME_T=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TSIZE_T=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TDEPT_T=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_NBEXT_T=$(grep "total number DT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# #     RR_NBVTS=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# #     RR_TIVTS=$(grep "total time spent in tests VT:" $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_PRECN_T=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_PPREC_T=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
#     RR_NBGRA_T=$(grep "number of Graeffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_NBGRR_T=$(grep "number of Graeffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
#     RR_TINGR_T=$(grep "time in Graeffe iterations:"    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
#     RR_TINRR_T=$(grep "time in computing root radii:"   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
# 
#     TMPSOLVE_S_T=$(grep "real" $NAME_OUTMPSOLVE| cut -f2 -d'l' | tr -d ' ')
# #     DSC_NSOLS_T=$(grep "Number of roots:" $NAME_OUTDSC| cut -f2 -d':'| tr -d ' ')
# #     DSC_TSIZE_T=$(grep "TREESIZE=" $NAME_OUTDSC| cut -f2 -d'='| tr -d ' ')
# #     DSC_TTIME_T=$(grep "real" $NAME_OUTDSC| cut -f2 -d'l' | tr -d ' ')
# 
# #     BITSI=`echo $BITSI+$BITSI_T|bc -l`
#     NSOLS=`echo $NSOLS+$NSOLS_T|bc -l`
#     TSIZE=`echo $TSIZE+$TSIZE_T|bc -l`
#     TDEPT=`echo $TDEPT+$TDEPT_T|bc -l`
#     TTIME=`echo $TTIME+$TTIME_T|bc -l`
#     TTIME_SQ=`echo $TTIME_SQ+$TTIME_T^2|bc -l`
#     NBEXT=`echo $NBEXT+$NBEXT_T|bc -l`
#     RR_NSOLS=`echo $RR_NSOLS +$RR_NSOLS_T|bc -l`
#     RR_TTIME=`echo $RR_TTIME +$RR_TTIME_T|bc -l`
#     RR_TTIME_SQ=`echo $RR_TTIME_SQ+$RR_TTIME_T^2|bc -l`
#     RR_TSIZE=`echo $RR_TSIZE +$RR_TSIZE_T|bc -l`
#     RR_TDEPT=`echo $RR_TDEPT +$RR_TDEPT_T|bc -l`
#     RR_NBEXT=`echo $RR_NBEXT +$RR_NBEXT_T|bc -l`
#     RR_PRECN=`echo $RR_PRECN +$RR_PRECN_T|bc -l`
#     RR_PPREC=`echo $RR_PPREC +$RR_PPREC_T|bc -l`
#     RR_NBGRA=`echo $RR_NBGRA +$RR_NBGRA_T|bc -l`
#     RR_NBGRR=`echo $RR_NBGRR +$RR_NBGRR_T|bc -l`
#     RR_TINGR=`echo $RR_TINGR +$RR_TINGR_T|bc -l`
#     RR_TINRR=`echo $RR_TINRR +$RR_TINRR_T|bc -l`
#     TMPSOLVE_S=`echo $TMPSOLVE_S +$TMPSOLVE_S_T|bc -l`
#     TMPSOLVE_S_SQ=`echo $TMPSOLVE_S_SQ+$TMPSOLVE_S_T^2|bc -l`
# }
