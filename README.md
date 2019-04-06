Table of Contents
=================

   * [Introduction](#introduction)
   * [Notes for use](#notes-for-use)
   * [Change Log](#change-log)


# Introduction

## Family History informed Genetic Risk score (FHiGR score)


Step 1. Using weights from Genome Wide Association Study (GWAS) summary statistics you can calclate a GRS in your study.  
Step 2. Using self-reported family history you can better stratfiy samples with FHiGR score.  

## Acknowledgements
Code base authored by Kimmo Pääkkönen, Dr. Ida Surakka, and Brooke Wolford

Genome-wide polygenic score (GPS) weights from [Khera, A. V. et al. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. Nature Genetics 50, 1219 (2018).](10.1038/s41588-018-0183-z)  
Weights can be downloaded [here](http://www.broadcvdi.org/informational/data).

# Notes for use

## :warning: This branch is under development

## Dependencies 
bcftools, qctool v2, split, python2, and python3  
python packages: argparse, glob, math, multiprocessing, numpy, os, pandas, random, scipy, string, subprocess, sys, time, tempfile  
R packages:  

## VCFtoDOSEforGRS.py 
Despite the name, this file converts .vcf to .dose or .bgen to .gen. This is the beginning, preparatory step for calculating GRS.
Written for Python 2.7.14. This script has low memory requirements, e.g. <1 G per chromosome with ~25K markers from the weights file. Parallization into k chunks is limited by number of cpus for Python's subprocess module to use.

* bgen to gen  
`python2 ~/scripts/VCFtoDOSEforGRS.py -f AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt -c 3 -p 4 -o UKBB.CHR3.AtrialFibrillation_PRS_LDpred_rho0.003 -b ukb_imp_chr3_v3.bgen -s ukb24460_imp_chr3_v3_s487395.sample -cn 3 -k 133`

* vcf to dose  
`python VCFtoDOSEforGRS.py -f AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt -c 3 -p 4 -o CHR14.AtrialFibrillation_PRS_LDpred_rho0.003_v3 -v CHR14.HRC_WGS.vcf.gz -cn 14 -k 25`

## DOSEtoGRS.py 
Despite the name, this file converts .dose or .gen to GRS. Written for Python 3.6.4. This script has high memory requirements.

* gen to GRS  

* dose to GRS  
`python3 DOSEtoGRS.py --input_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.22.dose --sample_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.sample --snp_weights_fn AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt --chromosome_no 17 --output_fn CHR17.AtrialFibrillation_PRS_LDpred_rho0.003_v3.22.gen --verbose`

## sumGRS.py 
Sum across the chromosomal chunks for one GRS per sample. Option to return inverse normalized score as well. Written for Python 2.7.14.  

`python sumGRS.py -s CHR22.AtrialFibrillation_PRS_LDpred_rho0.003_v3.sample -c AtrialFibrillation_PRS_LDpred_rho0.003_v3.config -o HUNT.AtrialFibrillation_PRS_LDpred_rho0.003_v3.GRS -i`

# Change Log
