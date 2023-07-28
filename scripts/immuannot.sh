#!/bin/bash
# Set some default values:
CONTIG=unset
REFDIR=unset

# optional parameters
OUTPREF=immuannot-out
OVERLAP=0.9
DIFF=0.03
THREAD=3

SCRIPTPATH=$(dirname $0)

usage()
{
  echo "
  Usage: bash ${SCRIPTPATH}/immuannot.sh  [OPTIONS] value
                           [ -c | --contig  target assembly (.fa, .fa.gz)       ] 
                           [ -r | --refdir  references                          ] 
                           [ -o | --outpref output prefix (optional)            ] 
                           [ -t | --thread  num of thread (optional, default 3) ] 
                           [ --overlaprate  OVERLAP (optional, default 0.9)     ] 
                           [ --diff         DIFF (optional, default 0.03)       ] 
                           "
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n immuannot \
    -o c:r:o:t: \
    --long contig:,refdir:,outpref:,thread:,overlaprate:,diff:\
    -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
  exit 1
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -c | --contig)  CONTIG="$2"    ; shift 2 ;;
    -r | --refdir)  REFDIR="$2"    ; shift 2 ;;
    -o | --outpref) OUTPREF="$2"   ; shift 2 ;;
    -t | --thread)  THREAD="$2"    ; shift 2 ;;
    --overlap)  OVERLAP="$2"    ; shift 2 ;;
    --diff)  DIFF="$2"    ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1."
       usage
       ;;
  esac
done

if [ $CONTIG == unset ]; then
  echo "Error: target contig seq is required."
  usage
fi

if [ $REFDIR == unset ]; then
  echo "Error: reference data set is required."
  usage
fi

if [ $OUTPREF == unset ]; then
  echo "Error: output prefix is required."
  usage
else
  mkdir -p $OUTPREF
  if [ $? != "0" ]; then
    echo "${OUTPREF} (-o) is not valid."
    usage
  fi
fi

start=`date +%s.%N`

echo "########################################"
echo "##Welcome###############################"
echo "########################################"
echo ""
echo "Starting time: ${start}"
echo "#####parameters:########################"
echo "CONTIG(-c)           : $CONTIG"
echo "REFDIR(-r)           : $REFDIR "
echo "OUTPREF(-o)          : $OUTPREF"
echo "THREAD(-t)           : $THREAD"
echo "OVERLAP(--overlap)   : $OVERLAP"
echo "DIFF(--diff)         : $DIFF"
#echo "Parameters remaining are: $@"




if true; then
  echo ""
  echo "#### processing IPD HLA/Kir genes#############"
  bash ${SCRIPTPATH}/annot.IPD.sh ${CONTIG} ${REFDIR} ${OUTPREF} \
    ${THREAD} ${OVERLAP} ${DIFF}
  if [ $? != "0" ]; then
    echo ERROR: Failed to annotate IPD HLA/Kir genes
    exit 1
  fi
fi

if true; then
  echo ""
  echo "#### processing C4 genes###############"
  bash ${SCRIPTPATH}/annot.C4.sh ${CONTIG} ${REFDIR} ${OUTPREF} ${THREAD}
  if [ $? != "0" ]; then
    echo ERROR: Failed to annotate C4 genes
    exit 1
  fi
fi

##combine the results
if true; then
  echo ""
  echo "#### combine IPD and C4 ###############"
  bash ${SCRIPTPATH}/annot.combine.sh ${REFDIR} ${OUTPREF}
  if [ $? != "0" ]; then
    echo ERROR: Failed to annotate combine IPD and C4
    exit 1
  fi
fi


end=`date +%s.%N`
echo "Ending time: ${end}"
echo "####$0 done"

