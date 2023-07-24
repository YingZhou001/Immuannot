#! /bin/bash

# Extract gene seq and CDS seq from IDP-HLA data base
# from the downloaded hla.dat

hladat=$1
kirdat=$2
c4dat=$3
refdir=$4

scriptpath=$(dirname $0)

easyhla="python3 ${scriptpath}/easyipd/easyIPD.py"
extractC4="python3 ${scriptpath}/extractC4Seq.py"

## switches
formatC4=false
formatIPDgeneStruct=false
generateGeneInformation=false
extractIPDgeneSeq=false
extractIPDCDSseq=true
simpleStatistics=false


if ! test -d ${refdir}; then
  mkdir -p ${refdir}
fi

## prepare C4 input data
if ${formatC4} ; then
  ${extractC4} ${c4dat} ${refdir}/C4
fi

## extract gene elements for selected alleles
gcsv=${refdir}/alleles.csv.gz
if ${formatIPDgeneStruct} ; then
  cat ${hladat} | ${easyhla} --elements | gzip -c > ${gcsv}
  cat ${kirdat} | ${easyhla} --elements | gzip  >> ${gcsv}
fi

## extract all gene names and add index
if ${generateGeneInformation} ; then
  datapref="IA"
  read -a genes <<< `zcat  ${gcsv} | cut -f1,2 | sort | uniq \
    | sed 's/gene=//g' | sed 's/\tversion=/,/g' | tr '\n' ' '`
  echo -n '' | gzip -c > ${refdir}/gene.names.txt.gz
  i=0
  for g in ${genes[@]}
  do
    i=$((i + 1))
    if [ ${i} -lt 10 ] ; then
      s="00000"${i}
    elif [ ${i} -lt 100 ] ; then
      s="0000"${i}
    else
      s=${i}
    fi
      echo ${g},${datapref}G${s},${datapref}T${s} \
        | gzip >> ${refdir}/gene.names.txt.gz
  done
  i=$((i + 1))
  s="0000"${i}
  echo "C4A,RefSeq,${datapref}G${s},${datapref}T${s}" \
    | gzip >> ${refdir}/gene.names.txt.gz
  i=$((i + 1))
  s="0000"${i}
  echo "C4B,RefSeq,${datapref}G${s},${datapref}T${s}" \
    | gzip >> ${refdir}/gene.names.txt.gz
  i=$((i + 1))
  s="0000"${i}
  echo "C4X,RefSeq,${datapref}G${s},${datapref}T${s}" \
    | gzip >> ${refdir}/gene.names.txt.gz
fi

##extract complete allele seqs with introns
gfa=${refdir}/gen.fa.gz
if ${extractIPDgeneSeq} ; then
  zgrep -v "partial" ${gcsv} | grep -v "intron=NA" \
    | cut -f 4 | sed "s/allele=//g" > ${refdir}/alleles.gene.txt
  zgrep -v "partial" ${gcsv} | grep "intron=NA" | grep "nexon=1" \
    | cut -f 4 | sed "s/allele=//g" >> ${refdir}/alleles.gene.txt
  cat ${hladat} ${kirdat} \
    | ${easyhla} --geneSeq --alleleF ${refdir}/alleles.gene.txt \
    | gzip -c > ${gfa}
fi

## extract cds seqs
cdsdir=${refdir}/CDSseq
if ${extractIPDCDSseq} ; then
  if ! test -d ${cdsdir}; then
    mkdir -p ${cdsdir}
  fi
  read -d '\r' -a genes <<< `zcat ${refdir}/gene.names.txt.gz | grep -v "^C4" \
    | cut -f1 -d','`
  for g in ${genes[@]}
  do
    echo ${g}
    zcat  ${gcsv} | grep -w "gene=${g}" | cut -f3 | sed "s/ID=//g" > tmp.cds.allele
    cfa=${cdsdir}/${g}.fa.gz
    cat ${hladat} ${kirdat} \
      | ${easyhla} --cdsSeq --alleleF tmp.cds.allele | gzip -c > ${cfa}
    rm tmp.cds.allele
  done
fi

## statistics
if ${simpleStatistics} ; then

  zcat ${refdir}/gene.names.txt.gz | cut -f1 -d',' > tmp.g
  zgrep 'goodTemplate=True' ${gcsv} | \
    cut -f1 | sed "s/gene=//g" > tmp.tlist #good template
  cut ${refdir}/alleles.gene.txt -f1 -d '*' | sed "s/>//g" > tmp.glist
  zcat ${gcsv} | cut -f1 | sed "s/gene=//g" > tmp.clist

  read -a genes <<< `zcat ${refdir}/gene.names.txt.gz | grep -v "^C4" \
    | cut -f1 -d','`

  echo "# complete: not partial, include introns" > ${refdir}/genSeqcounts.txt
  echo "# UTR: from UTR to UTR, good template" >> ${refdir}/genSeqcounts.txt
  echo "# updated date: $(date +'%F')" >> ${refdir}/genSeqcounts.txt
  echo "gene,complete,CDS,UTR" >> ${refdir}/genSeqcounts.txt

  echo -n '' > ${refdir}/gene.black.list
  for g in ${genes[@]}
  do
    count1=$(grep -w ${g} tmp.glist | wc -l)
    count2=$(grep -w ${g} tmp.clist | wc -l)
    count3=$(grep -w ${g} tmp.tlist | wc -l)
    if [ ${count1} == 0 ]; then
      echo ${g} >> ${refdir}/gene.black.list
    fi
    echo ${g},${count1},${count2},${count3} >> ${refdir}/genSeqcounts.txt
  done

  rm tmp.glist tmp.clist tmp.tlist tmp.g

fi
