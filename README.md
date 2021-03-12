Table of Contents
=================

   * [Introduction](#introduction)
   * [Notes for use](#notes-for-use)
   * [Change Log](#change-log)


# Introduction

## Family History informed Genetic Risk score (FHiGR score)
This project began as a way to create a unified score that combines genetic risk and family history—-a family history informed genetic risk score (FHiGRS). 
Step 1. Using weights from Genome Wide Association Study (GWAS) summary statistics you can estimate a GRS in your study.  
Step 2. Using self-reported family history you can better stratfiy samples with FHiGR score.  

In evaluating this score, we determined models with family history and GRS separately would be optimal. The FHiGRS moniker has persisted as the project evolved. These scripts allow for generating GRS from weights, creating data visualizations, and evaluating statistical models. They are generalized where noted. Deprecated scripts are noted below for reference.

## Acknowledgements
Code base authored by Dr. Sarah Graham, Kimmo Pääkkönen, Dr. Ida Surakka, and Brooke Wolford

Genome-wide polygenic score (GPS) weights from [Khera, A. V. et al. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. Nature Genetics 50, 1219 (2018).](https://doi.org/10.1038/s41588-018-0183-z)   
Weights can be downloaded [here](http://www.broadcvdi.org/informational/data). Weights can also be downloaded at the [PGS Catalog](pgscatalog.org).

# Notes for use

### :warning: This branch is under development

### Dependencies 
bcftools, qctool v2, split, python2, python3
python packages: argparse, glob, math, multiprocessing, numpy, os, pandas, random, scipy, string, subprocess, sys, time, tempfile, memory-profiler
R packages: optparse, data.table, and more

### File format descriptions 

* Weight file
SNPID EFFECT_ALLELE WEIGHT CHR POS A1 A2
header noted with #

* Variant Call Format (VCF)
CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT {GT:DS per sample}
header noted with #

* GEN (from BGEN v1.2)
CHR SNPID RSID POS A1 A2 {Genotype trios per sample}
with trios homozygous A1, heterozygous, homozygous A2
chromosome number has leading 0 for 1-9
associated .sample file should be the order of the genotype trios

* Phenotype file
Covariates (e.g. sex, birthyear, genetic principal components)
Stratification column (e.g. family history) should be coded as 1/0, NA allowed
Phenotype column for binary trait should be coded as 1/0, NA allowed

# Scripts

### processWeights.py
Takes a weights file and writes chunked weight and region files for parallelizing GRS estimation. 

`python processWeights.py -w weights.txt -c 3 -o output_prefix -cc 3 -pc 4 -rc 5 -ac 6 -k --num_chunk 10000`

### processRegionWeights.py
Takes a weights file and refers to a file with chromosomal start and end coordinates which were used for creating chunked genetic data files. Writes region and weights files for specified marker chunks and also a config file matching the weight and region files to corresponding genetic information (e.g. 20 weight files match with one VCF which is one of 20 VCFs for a given chromosome). Not generalized.

`python processRegionWeights.py -w weights.txt -cc 3 -pc 4 -rc 5 -ac 6 -o output_prefix`

### processDosages.py
Estimates GRS using weights file and regions file from above and a VCF.

`python processDosages.py  -r regions.txt -w weights.txt -v <vcf> -cc 3 -pc 4 -rc 5 -ac 6 -o output_prefix -i samples.txt`

### sumGRS.py 
Sum across the chromosomal chunks for one GRS per sample. Option to return inverse normalized score as well. Also checks files have the same number of samples. Written for Python 2.7.14.  

`python sumGRS.py -s CHR22.AtrialFibrillation_PRS_LDpred_rho0.003_v3.sample -c AtrialFibrillation_PRS_LDpred_rho0.003_v3.config -o HUNT.AtrialFibrillation_PRS_LDpred_rho0.003_v3.GRS -i`

### sumGRS.R
Sum across chunks of samples and scores to create one GRS per sample. Inverse normalizes score. Takes a file list with each file name on a new line and an output prefix.

 `Rscript ~/FHiGR_score/sumGRS.R file_list.txt output_prefix`
 
### checkDim.R
Checks number of rows and columns for each file in a provided file. Useful when running massively parallelized/chunked job arrays.

### FHiGR_raincloud.R
Creates visualizations of GRS distirbutions. Requires R_rainclouds.R

### FHiGR_prevalence.R
Generalized to create visualations of disease prevalence binned by GRS quantiles and stratified by a provided variable. Requires helperFunctions.R.

### FHiGR_ROC.R
Create ROC curves for FHiGRS and GRS.

### FHiGR_model.R
Generalized to check various association models using family history and genetic risk as predictors.

# :no_entry_sign: Deprecated

### estimateGRS.py
Revised from Dr. Sarah Graham. Uses subprocess module to parallelize tabix calls. High memory and time requirements make this inadvisable for large datasets.

### VCFtoDOSEforGRS.py 
Revised from Kimmo Pääkkönen & Dr. Ida Surakka. Despite the name, this file converts .vcf to .dose AND .bgen to .gen. This is the beginning, preparatory step for calculating GRS.
Written for Python 2.7.14. This script has low memory requirements, e.g. <1 G per chromosome with ~25K markers from the weights file. Parallization into k chunks is limited by number of cpus for Python's subprocess module to use.

* .bgen to .gen  
`python2 ~/scripts/VCFtoDOSEforGRS.py -f AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt -c 3 -p 4 -o UKBB.CHR3.AtrialFibrillation_PRS_LDpred_rho0.003 -b ukb_imp_chr3_v3.bgen -s ukb24460_imp_chr3_v3_s487395.sample -cn 3 -k 133`

* .vcf to .dose  
`python VCFtoDOSEforGRS.py -f AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt -c 3 -p 4 -o CHR14.AtrialFibrillation_PRS_LDpred_rho0.003_v3 -v CHR14.HRC_WGS.vcf.gz -cn 14 -k 25`

### DOSEtoGRS.py 
Revised from Kimmo Pääkkönen & Dr. Ida Surakka. Despite the name, this file converts .dose AND .gen to GRS. Written for Python 3.6.4. This script has high memory requirements.

* gen to GRS  

* dose to GRS  
`python3 DOSEtoGRS.py --input_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.22.dose --sample_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.sample --snp_weights_fn AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt --chromosome_no 17 --output_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.22.gen --verbose`

### FHiGR_score.R
Calculate odds ratios and visualize differences in disease prevalence across genetic risk score quantiles and stratum.

``
mkdir figures  
cut=`seq 0.8 0.01 0.99 | tr "\n" "," | sed 's/,$//g'` ; Rscript FHiGR_score.R -f MI_pheno_GRS.txt -s 16 -p 15 -g 3 --maintitle 'Family History of MI informed CAD GRS' --xlabel 'CAD GRS' -ylabel 'MI Prevalence' --legend 'MI Family History'  --cut_points $cut -o figures/MI
``

# Change Log

* Updates March 2021 for R version 4.0.4 and package updates 
