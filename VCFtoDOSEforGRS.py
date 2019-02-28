#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2019 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# Revised from Dr. Ida Surakka pipeline
# Formatting script to prepare to calculate genome wide polygenic risk scores using weights form LDpred
#==============================================================================

############################
##### IMPORT MODULES #######
###########################
import os
import sys
import subprocess
import argparse
import math
import string
import random
import multiprocessing as mp
from tempfile import NamedTemporaryFile

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",help="File with variants and weights",type=str,required=True)
    parser.add_argument("-c","--chr",help="0-based column with chromosome",type=int,default=1)
    parser.add_argument("-p","--pos",help="0-based column with end position of vairant",type=int,default=3)
    parser.add_argument("-v","--vcf",help="VCF with genetic data",type=str,required=True)
    parser.add_argument("-k","--chunk",help="Chunk each .dose file into X chunks for parallelization",default=0)
    parser.add_argument("-o","--output",help="output prefix",type=str,required=True)
    args=parser.parse_args()
    return args

###############################
######## SUB ROUTINES  #########
###############################

def open_zip(f):
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderrr, "Opening gzipped file %s\n" % f
    elif f == "-":
        sys.exit("Cannot read file from stdin\n")
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
    return command

def multiallelic(dict):
    keys=dict.keys()
    marker=[]
    if len(keys) == 2: #handle tri allelic SNP
        m1=dict[keys[0]]
        m2=dict[keys[1]]
        maf1=m1[7].split(";")[1]
        maf2=m2[7].split(";")[1]
        #take the marker with the larger maf
        if float(maf1.split("=")[1]) < float(maf2.split("=")[1]): 
            marker=dict[keys[1]]
        elif float(maf1.split("=")[1]) > float(maf2.split("=")[1]):
            marker=dict[keys[0]]
        else:
            sys.stderr.write("Failure to handle multiallelic SNP\n")
    elif len(keys) == 1: 
        marker=dict[keys[0]]
    else:
        sys.stderr.write("Failure to handle potential multiallelic SNP\n")
    return marker
                    

#return dosage in terms of effect allele from GWAS
def checkAllele(e,r,a,dos):
    #sys.stderr.write("Checking effect allele\n")
    #print(e+"\t"+a+"\t"+r+"\t"+dos)
    if r == flipStrand(a) and a == flipStrand(r): #ambiguous genotype
        return 0 
    d=float(dos)
    if a==e:
        return d #dosage from vcf in terms of alternate allele 
    elif r==e:
        return 2-d #make dosage from vcf in terms of reference allele
    elif r!=e and a!=e:
        new_r=flipStrand(r)
        new_a=flipStrand(a)
        if new_r==e:
            return 2-d 
        elif new_a==e:
            return d
        else:
            sys.stderr.write("Situation not accounted for\n")
            return None
    else:
        sys.stderr.write("Situation not accounted for\n")
        return None
    

def readSamples(vcf):
    sys.stderr.write("Reading in sample names from vcf\n")
    vcf=vcf.replace("$","01") #assumes sample list is the same in all chromosome VCFs

    #intialize containers
    ddict={}
    cdict={}
    sampleList=[]

    #read samples from vcf
    proc=subprocess.Popen(["bcftools","query","-l",vcf],stdout=subprocess.PIPE)
    output=proc.stdout.read()
    for sampleID in output.split("\n"):
        if len(sampleID) > 0: #make sure sampleID is real, not new line
            ddict[sampleID]=0 #make key in dictionary for data
            cdict[sampleID]=0 #make key in dictionary for snp counts
            sampleList.append(sampleID) #make ordered list of samples
    return ddict,cdict,sampleList

def flipStrand(allele):
    """
    Flip to complementary allele
    """
    if allele=="A":
        return "T"
    elif allele=="T":
        return "A"
    elif allele=="C":
        return "G"
    elif allele=="G":
        return "C"
    else:
        sys.stderr.write("Allele is %s\n" % allele)
        return None

def readWeights(f,c,p):
    """
    Read file with weights from LDpred. Write into temporary bed file. 
    """
    #initialize temporary file
    sys.stderr.write("Writing temporary file for marker names in bed format\n")
    marker_bed = NamedTemporaryFile(delete=True)
    with open(marker_bed.name, 'w') as tmp:
        
        #open file with weights
        command=open_zip(f)
        counter=0
        with command as f:
            for line in f:
                ls=line.rstrip()
                if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit 
                    lineList=ls.split() #assumes whitespace delimiter, space or tab
                    counter+=1
                    tmp.write("\t".join([lineList[c],str(int(lineList[p])-1),lineList[p]]))
                    tmp.write("\n")
        f.close()
    print >> sys.stderr, "Number of markers to pull from VCF is %d\n" % counter
    return(marker_bed,counter) #return tmp file object
            
def callQuery(vcf,tmp,out,chunk,counter):
    """
    Turn VCF into a .dose format for specific list of markers. Assumes path to bcftools.
    """
    #chunking option, make temporary files for the chunuks and call bcftools on each and then close them
    if chunk > 0:
        print >> sys.stderr, "Performing bcftools query to pull markers frorm VCF %s and write to %d .dose files\n" % (vcf,chunk)
        splitPrefix=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        chunkString="/".join(["l",str(chunk)]) #l/N
        subprocess.call(["split",tmp.name,splitPrefix,"-n",chunkString,"-d","--additional-suffix=.bed"])
        tmpFileList=[]
        outFileList=[]
        for i in range(chunk):
            tmpFileList.append(''.join([splitPrefix,"%02d" % i,".bed"]))
            outFileList.append(".".join([out,"%02d" % i,"dose"]))
        #TO DO: multiprocessing or threading to speed this up
        for j in range(len(tmpFileList)):
            subprocess.call(["bcftools","query",vcf,"-R",tmpFileList[j],"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outFileList[j]])
            subprocess.call(["rm",tmpFileList[j]])

    #no chunking
    elif chunk==0:
        outName=".".join([out,"dose"])
        print >> sys.stderr, "Performing bcftools query to pull markers frorm VCF %s and write to %s\n" % (vcf,outName)
        subprocess.call(["bcftools","query",vcf,"-R",tmp.name,"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outName])
    return 0


    
#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args = get_settings()

    #makes bed file of markrers frorm weight file
    tmp_obj,counter=readWeights(args.file,args.chr, args.pos)

    #writes .dose file from VCF, extracts markers from weight file
    #may chunk each .dose file into X chunks
    callQuery(args.vcf,tmp_obj,args.output,int(args.chunk),counter)

    

    
main()
