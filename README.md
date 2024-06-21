# Immuannot pipeline introduction

Immuannot is a pipeline built on [minimap2](https://github.com/lh3/minimap2) 
to detect and annotate immunological genes for human genome assembly. 
By taking advantage of gene sequences from 
[IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/),
[IPD-KIR](https://www.ebi.ac.uk/ipd/kir/), and
[RefSeq](accession number NG_011638.1)
as references,
it is able to annotate HLA and KIR allele at full precision (if exists in the
reference data set) and to report novel allele by locating new mutations that do not
exist in the reference allele set.

Last update date : 04/09/2024

Content
--------

- [Detection Strategy](#detection-strategy)
- [Gene Coverage](#gene-coverage)
- [Installation](#installation)
    - [Requirement](#requirement)
    - [Download](#download)
- [Inputs and Outputs](#inputs-and-outputs)
- [A Running Example](#a-running-example)
- [Limitations](#limitations)


# Detection Strategy


<img src=figs/pipeline-simple.png width="550" height="700" />

An __HLA__ gene is detected if any reference allele sequences of that gene are well mapped 
to the target assembly (overlapping rate > 90%). 
In a region with multiple mapping allele sequences, the one with the smallest 
edit-distance and longest matching length is chosen as the template allele, 
which is futher used for gene structure annotation.

If the edit-distance of the chosen template allele is zero, which means a perfect
matching, the contig then is reproted to carry that allele. 
Otherwise, CDS is extracted for calling the allele type.

__C4__ gene is detected through split alignment. Taking the exons from
one gene sequence as the query, the boundaries of each exon, 
the length of the 9th intron, and key pipetides in the 26th exon can 
be estimated through aligning the query to the target contig. 
Those information are used for gene structure annotation and gene typing. 

Because each reference gene sequence could be mapped to different regions
of the target genome, copy number of a particular gene 
is reported naturally with the number of mapping clusters.

See our manuscript for more details:

**Full resolution HLA and KIR genes annotation for human genome assemblies** (doi: 10.1101/gr. 278985.124)

# Gene coverage

## __HLA__ genes: 

HLA-A,HLA-B,HLA-C,HLA-DMA,HLA-DMB,HLA-DOA,HLA-DOB,HLA-DPA1,HLA-DPA2,HLA-DPB1,HLA-DPB2,HLA-DQA1,HLA-DQA2,HLA-DQB1,HLA-DQB2,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,HLA-DRB5,HLA-E,HLA-F,HLA-G,HLA-HFE,HLA-H,HLA-J,HLA-K,HLA-L,HLA-N,HLA-P,HLA-S,HLA-T,HLA-U,HLA-V,HLA-W,HLA-Y,MICA,MICB,TAP1,TAP2,C4A,C4B

## __KIR__ genes: 

KIR2DL1,KIR2DL2,KIR2DL3,KIR2DL4,KIR2DL5A,KIR2DL5B,KIR2DP1,KIR2DS1,KIR2DS2,KIR2DS3,KIR2DS4,KIR2DS5,KIR3DL1,KIR3DL2,KIR3DL3,KIR3DP1,KIR3DS1

[\[top\]](#content)

# Installation

## Requirement

This pipeline is developped and desgined under linux environment, following
programs are required to be pre-intalled and added in the [system searching PATH](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/):

* minimap2 (test version 2.17-r941)
* python3 (test version 3.9.13)
* bash (test version 4.4.20)

## Download

Download the reference data set 'Data-2024Feb02.tar.gz' from zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10948964.svg)](https://doi.org/10.5281/zenodo.10948964)





```bash
git clone https://github.com/YingZhou001/Immuannot.git
cd Immuannot
# put Data-2024Feb02.tar.gz here and 
tar xvf Data-2024Feb02.tar.gz
```

Testing :

```bash
>> bash scripts/immuannot.sh
Error: target contig seq is required.

  Usage: bash scripts/immuannot.sh  [OPTION] value
                           [ -c | --contig  target assembly (.fa, .fa.gz)       ]
                           [ -r | --refdir  references                          ]
                           [ -o | --outpref output prefix (optional)            ]
                           [ -t | --thread  num of thread (optional, default 3) ]
                           [ --overlaprate  OVERLAP (optional, default 0.9)     ]
                           [ --diff         DIFF (optional, default 0.03)       ]
```


[\[top\]](#content)

# Inputs and Outputs

Immuannot takes genome assembly (-c/--contig, fasta/fastq format) and gene reference
sequences data directory (-r/--refdir along this software) as inputs. 
It is better to set your own output prefix, otherwise all the results will be
prefixed with "immuannot-out".

Optically, the number of threads ( -t/--thread), 
the proportion of the query sequence mapped to the target (--overlaprate),
and the differential ratio cutoff of the matched region (--diff) 
can all be customized as need.

The output annotation is in gzip-compressed [gtf](https://useast.ensembl.org/info/website/upload/gff.html) format.

In the attribute column, new allele would have a 'new' tag at the rightest field 
of the 'consensus' call, it is usually accompanied with one or more most similar
'alleles' calls. 
CDS sequence is also examined during the template search, typically, 
'template\_warning' would give warnings as "no-start\_codon", "no-stop\_codon", 
"incomplete\_CDS", and "inframe\_stop" for gene structure annotation.

Intermediate results are also saved in the output folder.

[\[top\]](#content)

# A Running example

A running example is included along this pipeline. 
The target file 'test.fa.gz' include two contigs, one for HLA region 
and the other for KIR region.

User can test the pipeline by running the script bellow (a few mins):

```bash
## check the file path before running, take a few minutes
ctg=example/test.fa.gz
script=scripts/immuannot.sh
refdir=Data-2024Feb02
outpref=test-run
bash ${script} -c ${ctg} -r ${refdir} -o ${outpref}
```

Immuannot would output file "test-run.gtf.gz" for annotation and a folder named
"test-run" for intermediate results.

[\[top\]](#content)

# Limitations

Because Immuannot is mainly based on gene sequence alignment, except for C4 gene, a novel gene may not be reported if it is largely different from the reference data sequences.

Immuannot takes one set of genome as the input for each run. Annotating multiple-sample sequences in one run could lead to lots of missing calls.

[\[top\]](#content)


# Todo list

 * Combine split-alignment and gene sequence aligment to get a better annotation for large novel indels
 * ~~use zenodo link to replace data larger than 1Mb~~
