#! /bin/bash

## genotyping HLA alleles

tctg=$1
refdir=$2
outpref=$3
mm2thread=$4

overlap=$5
diff=$6

scriptpath=$(dirname $0)

# refs
eledat=${refdir}/alleles.csv.gz
genfa=${refdir}/gen.fa.gz
geneinfo=${refdir}/gene.names.txt.gz
allelerecord=${refdir}/alleles.csv.gz

#check input
for file in ${eledat} ${genfa} ${geneinfo} ${allelerecord}
do
  if [ ! -f ${file} ] ; then
    echo "${file} does not exist"
    echo "Check the input directory: ${refdir}"
    exit 1
  fi
done

mm2paf=${outpref}/mm2.ipd.gen.paf.gz

if true ; then
  # Run minimap2 for the gene seq data
  options="-t ${mm2thread} -cx asm5 --cs --end-bonus=10"
  minimap2 ${options} ${tctg} ${genfa} | gzip -c > ${mm2paf}
fi

## will generate *csv, *cds.fa, *cds.prot, and *.newgene.txt
searchTemplate="python3 ${scriptpath}/searchTemplate.py"
${searchTemplate} ${mm2paf} ${eledat} ${genfa} \
  ${outpref} ${overlap} ${diff}


tcdsctg=${outpref}/cds.fa.gz
genelist=${outpref}/tmp.newgene.txt
calltmpfile=${outpref}/tmp.gene.csv

cdspaf=${outpref}/mm2.ipd.cds.paf.gz
qcdsseq=${outpref}/tmp.ipd.cds.fa.gz

if true ; then
  read -a genes <<< `cat ${genelist} | tr '\n' ' '`
  echo -n '' | gzip -c > ${cdspaf}
  echo -n '' | gzip -c > ${qcdsseq}
  for gene in ${genes[@]}
  do
    qctg=${refdir}/CDSseq/${gene}.fa.gz
    if [ ! -f ${qctg} ] ; then
      echo "${qctg} does not exist"; exit 1
    fi
    zcat ${qctg} | gzip >> ${qcdsseq}
  done
  option="-t ${mm2thread} -c --cs --end-bonus=10"
  minimap2 ${option} ${tcdsctg} ${qcdsseq} | gzip >> ${cdspaf}
fi

calltmpfile2=${outpref}/tmp.allele.csv
callnewallele="python3 ${scriptpath}/callIPDallele.py"
${callnewallele} ${cdspaf} ${qcdsseq} ${calltmpfile} > ${calltmpfile2}

##annotation to gtf file
annot="python3 ${scriptpath}/annotIPDallele.py"
gtfout=${outpref}/tmp.ipd.gtf.gz
${annot} ${calltmpfile2} ${allelerecord} ${geneinfo}| gzip -c > ${gtfout}
