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
    parser.add_argument("-v","--vcf",help="VCF with genetic data",type=str)
    parser.add_argument("-b","--bgen",help="BGEN with  genetic data",type=str)
    parser.add_argument("-vc","--vcf_chrom",help="Chromosome number of the VCF provided",type=int,default=0)
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
    

def readSamples(vcf,out):
    """
    Make sample file that corresponds to the vcf. To be used later with .dose files. Has header FID and IID.
    """
    outName=".".join([out,"sample"])
    print >> sys.stderr,"Reading in sample names from vcf to %s\n" % outName
    #open sample file for writing
    sampleFile=open(outName,"w")
    sampleFile.write("\t".join(["FID","IID"]))
    sampleFile.write("\n")
    
    #read samples from vcf
    proc=subprocess.Popen(["bcftools","query","-l",vcf],stdout=subprocess.PIPE)
    output=proc.stdout.read()
    
    for sampleID in output.split("\n"):
        if len(sampleID) > 0: #make sure sampleID is real, not new line
            sampleFile.write("\t".join([sampleID,sampleID]))
            sampleFile.write("\n")

    sampleFile.close()
            
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

def readWeights(f,c,p,v):
    """
    Read file with weights from LDpred. Write into temporary bed file. 
    """
    #initialize temporary file
    sys.stderr.write("Writing temporary file for marker names in bed format\n")
    marker_bed = NamedTemporaryFile(delete=True,suffix=".bed")
    with open(marker_bed.name, 'w') as tmp:
        
        #open file with weights
        command=open_zip(f)
        counter=0
        with command as f:
            for line in f:
                ls=line.rstrip()
                if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit 
                    lineList=ls.split() #assumes whitespace delimiter, space or tab
                    if v > 0: #vcf is just for one chromosome so we can ignore weights from other chrom
                        if v==int(lineList[c]): 
                            counter+=1
                            tmp.write("\t".join([lineList[c],str(int(lineList[p])-1),lineList[p]]))
                            tmp.write("\n")
                    else: #vcf is for all chromosomes 
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
        chunkString="/".join(["l",str(chunk)]) #l/N, chunking parameter
        subprocess.call(["split",tmp.name,splitPrefix,"-n",chunkString,"-d","--additional-suffix=.bed"])
        tmpFileList=[]
        outFileList=[]
        for i in range(chunk):
            tmpFileList.append(''.join([splitPrefix,"%02d" % i,".bed"])) #what temporary bed files have the chunked markers 
            outFileList.append(".".join([out,"%02d" % i,"dose"]))  #what will the final .dose files be called
        cmds_list=[]
        procs_list=[]
        procs2_list=[]
        rm_list=[]
        for j in range(len(tmpFileList)):
            cmds_list.append(["bcftools","query",vcf,"-R",tmpFileList[j],"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outFileList[j]])
            rm_list.append(["rm",tmpFileList[j]])
        for k in range(len(cmds_list)):
            procs_list.append(subprocess.Popen(cmds_list[k],stdout=subprocess.PIPE,stderr=subprocess.PIPE))
        for proc in procs_list: #run bcftools parallel
            proc.wait()
        #clean up
        print >> sys.stderr, "Removing temporary .bed files for chunking\n"
        for r in range(len(rm_list)):
            procs2_list.append(subprocess.Popen(rm_list[r],stdout=subprocess.PIPE,stderr=subprocess.PIPE))
        for proc2 in procs2_list: #run rm of tmp bed files parallel
            proc2.wait()

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

    #VCF data
    if args.bgen == None:
    
        #makes bed file of markrers from weight file
        tmp_obj,counter=readWeights(args.file,args.chr, args.pos,args.vcf_chrom)

        #make sample file from VCF
        readSamples(args.vcf,args.output)
    
        #writes .dose file from VCF, extracts markers from weight file
        #may chunk each .dose file into X chunks
        callQuery(args.vcf,tmp_obj,args.output,int(args.chunk),counter)

    else if  args.vcf == None:

        x=5
        #BGEN TO DOSE  IN CHUNKS

    else: 
        print>> sys.stderr,"Must provide --vcf or --bgen\nb"

    
main()
