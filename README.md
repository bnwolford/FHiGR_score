Table of Contents
=================

   * [Introduction](#introduction)
   * [Notes for use](#notes-for-use)
   * [Change Log](#change-log)


# Introduction

## FHiGR score
Family history informed genetic risk score

## Acknowledgements
Code base authored by Kimmo Pääkkönen, Dr. Ida Surakka, and Brooke Wolford

Genome-wide polygenic score (GPS) weights from [Khera, A. V. et al. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. Nature Genetics 50, 1219 (2018).](10.1038/s41588-018-0183-z)
Weights can be downloaded [here](http://www.broadcvdi.org/informational/data)

# Notes for use

## :warning: This branch is under development

## Dependencies 
bcftools, qctool v2, split, python2, and python3 

## VCFtoDOSEforGRS.py 
Despite the name, this file converts .vcf to .dose or .bgen to .gen. This is the beginning, preparatory step for calculating GRS.
Written for Python 2.7.14.

`python2 ~/scripts/VCFtoDOSEforGRS.py -f AtrialFibrillation_PRS_LDpred_rho0.003_v3.txt -c 3 -p 4 -o UKBB.CHR3.AtrialFibrillation_PRS_LDpred_rho0.003 -b ukb_imp_chr3_v3.bgen -s ukb24460_imp_chr3_v3_s487395.sample -cn 3 -k 133`

## DOSEtoGRS.py 
Despite the name, this file converts .dose or .gen to GRS. 

## sumGRS.py 

# Change Log
