ctgfa=$1
refdir=$2
outpref=$3
mm2thread=$4

scriptpath=$(dirname $0)

exonfa=${refdir}/C4.exon.fa.gz
geneinfo=${refdir}/gene.names.txt.gz

# check input files:
if [ ! -f ${exonfa} ]; then
  echo "${exonfa} does not exist"
  exit 1
fi

if [ ! -f ${geneinfo} ]; then
  echo "${geneinfo} does not exist"
  exit 1
fi


exonpaf=${outpref}/mm2.c4.exon.paf.gz
mm2options="-t ${mm2thread} -C5 --cs -cx splice:hq"
minimap2 ${mm2options} ${ctgfa} ${exonfa} | gzip -c > ${exonpaf}


outgtf=${outpref}/tmp.c4.gtf.gz
python3 ${scriptpath}/callC4Allele.py \
  ${exonpaf} ${exonfa} ${ctgfa} ${geneinfo} | gzip -c > ${outgtf} 
