
## dowbload the latest kir.dat from https://github.com/ANHIG/IPDKIR
wget https://raw.githubusercontent.com/ANHIG/IPDKIR/refs/heads/Latest/kir.dat
## download the latest hla.dat from https://github.com/ANHIG/IMGTHLA
wget https://github.com/ANHIG/IMGTHLA/raw/refs/heads/Latest/hla.dat.zip

mv kir.dat kir.2.13.0.dat
unzip hla.dat.zip
mv hla.dat hla.3.58.0.dat

# output
refdir=Data-$(date +%Y%b%d)
hla=hla.3.58.0.dat
kir=kir.2.13.0.dat
c4=c4a.gb

mkdir -p ${refdir}

tooldir=../../tools # you need to change this to your own path to the script
bash ${tooldir}/predata.sh ${hla} ${kir} ${c4} ${refdir}

cat curated/DRB3.ele.csv | gzip >> ${refdir}/alleles.csv.gz
cat curated/DRB3.new.fa | gzip >> ${refdir}/gen.fa.gz
cat curated/DRB3.new.cds.fa | gzip >> ${refdir}/CDSseq/HLA-DRB3.fa.gz
