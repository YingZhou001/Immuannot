refdir=$1
outpref=$2

cleantmp=true

ipdgtf=${outpref}/tmp.ipd.gtf.gz
c4gtf=${outpref}/tmp.c4.gtf.gz
geneinfo=${refdir}/gene.names.txt.gz
blacklist=${refdir}/gene.black.list

scriptpath=$(dirname $0)

python3 ${scriptpath}/combine.py ${ipdgtf} ${c4gtf} ${geneinfo} \
  ${blacklist} | gzip -c > ${outpref}.gtf.gz

if ${cleantmp} ; then
  rm ${outpref}/tmp.ipd.gtf.gz 
  rm ${outpref}/tmp.c4.gtf.gz
  rm ${outpref}/tmp.newgene.txt
  rm ${outpref}/tmp.allele.csv
  rm ${outpref}/tmp.gene.csv
  rm ${outpref}/tmp.ipd.cds.fa.gz
fi
