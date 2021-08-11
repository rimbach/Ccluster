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

if [ -z "$PURGERR" ]; then
   PURGERR=0
fi

if [ -z "$PURGEDSC" ]; then
   PURGEDSC=0
fi

if [ -z "$NBPOLS" ]; then
   NBPOLS=10
fi

CCLUSTER_PATH="../../"
RISOLATE_CALL=$CCLUSTER_PATH"/bin/risolate"
GENPOLFI_CALL=$CCLUSTER_PATH"/bin/genPolFile"
GENRANDPOLFI_CALL=$CCLUSTER_PATH"/bin/genRandPolFile"
RISOLATE_OPTS="-v 2"

ANEWDSC_PATH="../../../softs"
ANEWDSC_CALL=$ANEWDSC_PATH"/test_descartes_linux64"

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
            CALL="$RISOLATE_CALL $NAME_IN -v 2 -m onlySubd"
            $CALL > $NAME_OUT
    fi
    
    if [ $PURGERR -eq 1 ]; then
        rm -f $NAME_OUTRR
    fi
  
    if [ ! -e $NAME_OUTRR ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, global, root radii, output in " $NAME_OUTRR
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
    
    if [ ! -e $NAME_OUT ]; then
            echo  "Isolating real roots for $POLNAME degree $DEG, with aNewDsc, output in " $NAME_OUT
            (/usr/bin/time -f "real\t%e" $CALL $NAME_IN > $NAME_OUT) &>> $NAME_OUT
    fi
    
}

TEMPTABFILE="temptab_risolate.txt"
touch $TEMPTABFILE

stats_pol()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    NAME_OUTDSC=$NAME".out_dsc"
    DEG=$2
    
    BITSI=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NSOLS=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TSIZE=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT=$(grep "tree depth:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCOT=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NSOLS=$(grep "number of real roots:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBEXT=$(grep "total number DT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBCOT=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PRECN=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_PPREC=$(grep "precision required/predicted:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_NBGRA=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBGRR=$(grep "number of Greaffe Iterations:"  $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f2 -d'|' | tr -d ' ')
    RR_TINGR=$(grep "time in Graeffe iterations:"    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TINRR=$(grep "time in computing root radii:"   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')

    DSC_NSOLS=$(grep "Number of roots:" $NAME_OUTDSC| cut -f2 -d':'| tr -d ' ')
    DSC_TSIZE=$(grep "TREESIZE=" $NAME_OUTDSC| cut -f2 -d'='| tr -d ' ')
    DSC_TTIME=$(grep "real" $NAME_OUTDSC| cut -f2 -d'l' | tr -d ' ')
    
    LINE_TAB3="`format_numb $DEG $LENP` & `format_numb $BITSI $LENP` & `format_numb $NSOLS $LENP`"
    LINE_TAB3=$LINE_TAB3" & `format_time $TTIME`     & `format_numb $NBEXT $LENP`,`format_numb $NBCOT $LENP`"
    LINE_TAB3=$LINE_TAB3" & `format_time $RR_TTIME`  & `format_numb $RR_NBEXT $LENP`,`format_numb $RR_NBCOT $LENP`"
    LINE_TAB3=$LINE_TAB3" & `percent_time $RR_TINRR $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
    LINE_TAB3=$LINE_TAB3" & `format_time $DSC_TTIME`\\\\"
    
    echo $LINE_TAB3 >> $TEMPTABFILE
}

stats_pol_rand()
{
    NAME=$1
    NAME_OUT=$NAME".out_riso"
    NAME_OUTRR=$NAME".out_riso_rr"
    NAME_OUTDSC=$NAME".out_dsc"
    BITSI_T=$(grep "bitsize: "                     $NAME_OUT| tr -d ' ' | cut -f3 -d':' | cut -f1 -d'|')
    NSOLS_T=$(grep "number of real roots:"          $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    
    TSIZE_T=$(grep "tree size:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TDEPT_T=$(grep "tree depth:"                    $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBEXT_T=$(grep "total number DT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    NBCOT_T=$(grep "total number VT:"              $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    TTIME_T=$(grep "total time:"                   $NAME_OUT| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NSOLS_T=$(grep "number of real roots:"       $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TTIME_T=$(grep "total time:"                   $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TSIZE_T=$(grep "tree size:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_TDEPT_T=$(grep "tree depth:"                    $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBEXT_T=$(grep "total number DT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
    RR_NBCOT_T=$(grep "total number VT:"              $NAME_OUTRR| cut -f2 -d':'| cut -f1 -d's' | cut -f1 -d'|' | tr -d ' ')
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
    
}

source ./functions.sh
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
echo $POLNAME >> $TEMPTABFILE
DEGREES="256 391 512"
BITSIZES="8192 16384 32768 65536"
# DEGREES="256"
# BITSIZES="8192"
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
            
            LINE_TAB3="`format_numb $DEG $LENP` & `format_numb $BIT $LENP` & `format_time $NSOLS`"
            LINE_TAB3=$LINE_TAB3" & `format_time $TTIME` (`format_time $TTIME_SQ`) & `format_time $NBEXT`,`format_time $NBCOT`"
            LINE_TAB3=$LINE_TAB3" & `format_time $RR_TTIME` (`format_time $RR_TTIME_SQ`)  & `format_time $RR_NBEXT`,`format_time $RR_NBCOT`"
            LINE_TAB3=$LINE_TAB3" & `percent_time $RR_TINRR $RR_TTIME` & `percent_time $RR_TTIME $TTIME`"
            LINE_TAB3=$LINE_TAB3" & `format_time $DSC_TTIME` (`format_time $DSC_TTIME_SQ`)\\\\"
            
            echo $LINE_TAB3 >> $TEMPTABFILE
    done
done

DEGREES="256 391 512 791 1024"
# DEGREES="256"
POLNAMES="Bernoulli Wilkinson"

for POLNAME in $POLNAMES; do
    echo $POLNAME >> $TEMPTABFILE
for DEG in $DEGREES; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG
    
    gen_with_deg $NAME $POLNAME $DEG
    run_risolate $NAME $POLNAME $DEG
    run_aNewDsc  $NAME $POLNAME $DEG 1
    
    stats_pol $NAME $DEG
done 
done
# # 
# 
# 
POLNAME="RegularGrid"
SIZEGRID="8 10 12 14 16"
# SIZEGRID="8"

echo $POLNAME >> $TEMPTABFILE
for SIZ in $SIZEGRID; do
    
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$SIZ
    
    gen_with_deg $NAME $POLNAME $SIZ
    run_risolate $NAME $POLNAME $SIZ
    run_aNewDsc  $NAME $POLNAME $SIZ 1
    
    DEG=$(( 2*$SIZ+1 ))
    DEG=$(( $DEG*$DEG ))
    stats_pol $NAME $DEG
done

# 
POLNAME="Mignotte"
# with increasing bit-size

echo $POLNAME >> $TEMPTABFILE

DEGF=512
BITSIZES="127 255 511 1023 2047"
# BITSIZES="127"
for BIT in $BITSIZES; do
    
    DEG=$DEGF
    REPNAME=$REP
    NAME=$REPNAME"/"$POLNAME"_"$DEG"_"$BIT
    
    gen_with_deg_bs $NAME $POLNAME $DEG $BIT
    run_risolate $NAME $POLNAME $DEG
    run_aNewDsc  $NAME $POLNAME $DEG 1
    
    stats_pol $NAME $DEG
done

cat $TEMPTABFILE
rm -f $TEMPTABFILE
