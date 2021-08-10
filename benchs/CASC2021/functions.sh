#!/bin/bash

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
    
    if [ ! -e $NAME_IN2 ]; then
            echo  "Generating file for $POLNAME degree $DEG, pol in " $NAME_IN2
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2
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
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 3 -b $BS
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
    NAME_IN3=$NAME".dsc"
    
    if [ ! -e $NAME_IN ]; then
            echo  "Generating file for $POLNAME degree $DEG bitsize $BS, pol in " $NAME_IN
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN -f 1 -b $BS -n $NBT
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN2 -f 2 -b $BS -n $NBT
            $GENPOLFI_CALL $POLNAME $DEG $NAME_IN3 -f 3 -b $BS -n $NBT
    fi
    
}
