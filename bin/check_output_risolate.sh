#!/bin/bash

BGREEN="\e[1;32m"
BRED="\e[1;31m"
NORMAL="\e[0m"

absolute_path(){
   dirlist="$1"
   retval=""
   for dir in $dirlist; do
      case $dir in
        /*) dir=$dir;;
        *) dir=$PWD/$dir;;
      esac
      retval=$retval" "$dir
   done
   echo $retval
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | sed 's/=.*//'`
    VALUE=`echo $1 | sed 's/[^=]*//; s/=//'`
    case "$PARAM" in
        --nbSols)
          NBSOLS=$VALUE
          ;;
        --nbClus)
          NBCLUS=$VALUE
          ;;
        --file)
          FILEIN=$(absolute_path "$VALUE")
    esac
    shift
done

NBCLUSF=$(grep "number of distinct real roots:"             $FILEIN| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')
NBSOLSF=$(grep "number of real roots:"            $FILEIN| cut -f2 -d':' | cut -f1 -d'|' | tr -d ' ')

# echo "number of clusters : " $NBCLUS " " $NBCLUSF
# echo "number of solutions: " $NBSOLS " " $NBSOLSF

if [ "$NBCLUS" = "$NBCLUSF" ] && [ "$NBSOLS" = "$NBSOLSF" ]; then
    echo -e $BGREEN "test PASSED" $NORMAL
else  
    echo -e $BRED "test FAILED: should have found " $NBSOLS " roots in " $NBCLUS "clusters"  $NORMAL
fi
