

# Pipeline introduction

Immuannot is a pipeline building on [minimap2](https://github.com/lh3/minimap2) 
to annotate immunological genes for human genome assembly. 
By taking adavntage of gene sequences from 
[IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/),
[IPD-KIR](https://www.ebi.ac.uk/ipd/kir/), and
[National library of medicine](https://www.ncbi.nlm.nih.gov/gene/720),
it is able to detect HLA and KIR alleles at full precision (if exist in the
reference data set) or report novel alleles by locating the postions of
mutations.

Last update date : \Today

# Detection Strategy

# Gene coverage

* __HLA__ genes: HLA-A,HLA-B,HLA-C,HLA-DMA,HLA-DMB,HLA-DOA,HLA-DOB,HLA-DPA1,HLA-DPA2,HLA-DPB1,HLA-DPB2,HLA-DQA1,HLA-DQA2,HLA-DQB1,HLA-DQB2,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,HLA-DRB5,HLA-E,HLA-F,HLA-G,HLA-HFE,HLA-H,HLA-J,HLA-K,HLA-L,HLA-N,HLA-P,HLA-S,HLA-T,HLA-U,HLA-V,HLA-W,HLA-Y, MICA,MICB,TAP1,TAP2,C4A,C4B
* __KIR__ genes: KIR2DL1,KIR2DL2,KIR2DL3,KIR2DL4,KIR2DL5A,KIR2DL5B,KIR2DP1,KIR2DS1,KIR2DS2,KIR2DS3,KIR2DS4,KIR2DS5,KIR3DL1,KIR3DL2,KIR3DL3,KIR3DP1,KIR3DS1

# Installation

## Requirement

This pipeline is developped and desgined under linux environment, following
programs are required to be pre-intalled:

* minimap2 (test version 2.17-r941)
* python3 (test version 3.9.13)
* bash (test version 4.4.20)

## Download and install

(TBD)

# Inputs and Outputs

# A Running example

# Limitations

# Lisence
