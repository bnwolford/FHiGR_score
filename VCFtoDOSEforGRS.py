#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2019 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# Revised from Dr. Ida Surakka pipeline
# Step 1: Convert bgen to gen or vcf to dose using a weights file from LDpred with the variants of interest for the genome wide polygenic risk score
#Requires bcftools or qctool_v2 as well as unix's split function
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
    parser.add_argument("-c","--chr",help="0-based column with chromosome in file",type=int,default=1)
    parser.add_argument("-p","--pos",help="0-based column with end position of variant in file",type=int,default=3)
    parser.add_argument("-v","--vcf",help="VCF with genetic data. Will convert to .dose.",type=str)
    parser.add_argument("-b","--bgen",help="BGEN with  genetic data. Will convert to .gen",type=str)
    parser.add_argument("-s","--sample",help="BGEN sample file.",type=str)
    parser.add_argument("-cn","--chr_num",help="Chromosome number of the BGEN or VCF provided if for a single chromosome (not recommended to provide entire genome genetic file, can have leading 0)",type=int,default=0)
    parser.add_argument("-k","--chunk",help="Chunk each .dose file into X chunks for parallelization",default=0)
    parser.add_argument("-o","--output",help="output prefix",type=str,required=True)
    parser.add_argument("--qctool",help="qctool path",type=str,default="/net/snowwhite/home/bwolford/qctool/build/release/qctool_v2.0.1")
    parser.add_argument("--bcftools",help="bcftools path",type=str,default="/usr/local/bin/bcftools")
    parser.add_argument("--split",help="split path",type=str,default="/usr/bin/split")
    parser.add_argument("-u","--cpu",help="Number of CPU cores to utilize for multiprocessing",default=8)
    args=parser.parse_args()
    print >> sys.stderr, "%s\n" % args
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

def readWeights(f,c,p,n):
    """
    Read file with weights from LDpred. Write into temporary bed file. Used to subset VCF provided.
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
                    if n > 0: #vcf is just for one chromosome so we can ignore variants from other chromosome
                        if n==int(lineList[c]): 
                            counter+=1
                            tmp.write("\t".join([lineList[c],str(int(lineList[p])-1),lineList[p]]))
                            tmp.write("\n")
                    else: #vcf is for all chromosomes 
                        counter+=1
                        tmp.write("\t".join([lineList[c],str(int(lineList[p])-1),lineList[p]]))
                        tmp.write("\n")
        f.close()
    tmp.close()
    print >> sys.stderr, "Number of markers to pull from VCF is %d\n" % counter
    return(marker_bed,counter) #return tmp file object
            
def callQuery(vcf,tmp,out,chunk,counter,bcftools,split,cpu):
    """
    Turn VCF into a .dose format for specific list of markers. Assumes path to bcftools.
    """
    #chunking option: make temporary files for the chunuks and call bcftools on each and then close them
    if chunk > 0:
        splitPrefix=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        chunkString="/".join(["l",str(chunk)]) #l/N, chunking parameter
        print >> sys.stderr, "Performing bcftools query to pull markers frorm VCF %s and write to %d .dose files. Temporary .bed files have prefix %s\n" % (vcf,chunk,splitPrefix)
        subprocess.call([split,tmp.name,splitPrefix,"-a","3","-n",chunkString,"-d","--additional-suffix=.bed"])

        tmpFileList=[]
        outFileList=[]
        for i in range(chunk):
            tmpFileList.append(''.join([splitPrefix,"%03d" % i,".bed"])) #what temporary bed files have the chunked markers 
            outFileList.append(".".join([out,"%03d" % i,"dose"]))  #what will the final .dose files be called
        cmds_list=[]
        procs_list=[]
        procs2_list=[]
        rm_list=[]
        for j in range(len(tmpFileList)):
            cmds_list.append(["bcftools","query",vcf,"-R",tmpFileList[j],"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outFileList[j]])
            rm_list.append(["rm",tmpFileList[j]])

        print >> sys.stderr, "Using multiprocessing pool functionality with %d cpus and %d chunks" % (cpu,chunk)
        pool = mp.Pool(cpu) #use user defined number of cpus, user should also specify for job scheduler
        results_list=pool.map(cmd_executor,cmds_list)
        if (sum(results_list)) != 0:
            print >> sys.stderr, "At least one parallelized bcftools query failed\n"
        print >> sys.stderr, "Removing temporary .bed files for chunking\n"
        results_list=pool.map(cmd_executor,rm_list)
            
    #no chunking
    elif chunk==0:
        outName=".".join([out,"dose"])
        print >> sys.stderr, "Performing bcftools query to pull markers frorm VCF %s and write to %s\n" % (vcf,outName)
        subprocess.call([bcftools,"query",vcf,"-R",tmp.name,"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outName])
    return 0

#called with multiprocessing pool when chunk > 5
def cmd_executor(cmd):
    pid=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    pid.wait()
    return 0


def readWeightsForBgen(f,c,p,n):
    """
    Read file with weights from LDpred. Write into temporary text file. Used to subset BGEN provided.
    """
    #initialize temporary file
    print >> sys.stderr, "Writing temporary file for marker names in chr:pos space delimited file\n"
    marker_file = NamedTemporaryFile(delete=True,suffix=".txt")
    with open(marker_file.name, 'w') as tmp:
        #open file with weights
        command=open_zip(f)
        counter=0
        with command as f:
            for line in f:
                ls=line.rstrip()
                if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit
                    lineList=ls.split() #assumes whitespace delimiter, space or tab
                    if n > 0: #bgen is just for one chromosome so we can ignore variants from other chrom
                        if n==int(lineList[c]):
                            counter+=1
                            chrom="0"+lineList[c] if len(lineList[c]) < 2 else lineList[c] #bgen -incl-positions requires CC:pos
                            tmp.write(":".join([chrom,lineList[p]]))
                            tmp.write(" ")
                    else: #bgen is for all chromosomes
                        counter+=1
                        chrom="0"+lineList[c] if len(lineList[c]) < 2 else lineList[c] #bgen -incl-positions requires CC:pos
                        tmp.write("\t".join([chrom,lineList[p]]))
                        tmp.write(" ")
        f.close()
    tmp.close()
    print >> sys.stderr, "Number of markers to pull from BGEN is %d\n" % counter
    return(marker_file,counter) #return tmp file object
                        

def bgenToGen(bgen,tmp,out,chunk,counter,qctool,split,cpu):
    """
    Turn .bgen into a .gen format for specific list of markers. Chunks into given number of chunks and runs qctool command parallely.
    """
    #chunking option, make temporary files for the chunks and call qctool on each and then close them
    if chunk > 0:
        splitPrefix=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        chunkString="/".join(["l",str(chunk)]) #l/N, chunking parameter
        print >> sys.stderr, "Performing qctool query to pull markers from BGEN %s and write to %d .gen files. Temporary .txt files have prefix %s\n" % (bgen,chunk,splitPrefix)
        subprocess.call([split,tmp.name,splitPrefix,"-a","3","-t"," ","-n",chunkString,"-d","--additional-suffix=.txt"])
        tmpFileList=[]
        outFileList=[]
        for i in range(chunk):
            tmpFileList.append(''.join([splitPrefix,"%03d" % i,".txt"])) #what temporary text files have the chunked markers
            outFileList.append(".".join([out,"%03d" % i,"gen"]))  #what will the final .gen files be called
        cmds_list=[]
        procs_list=[]
        procs2_list=[]
        rm_list=[]
        for j in range(len(tmpFileList)):
            cmds_list.append([qctool,"-g",bgen,"-incl-positions",tmpFileList[j],"-og",outFileList[j]])
            rm_list.append(["rm",tmpFileList[j]])

        print >> sys.stderr, "Using multiprocessing pool functionality with %d cpus and %d chunks\n" % (cpu,chunk)
        pool = mp.Pool(cpu) #use user defined cpus, should also submit this to job scheduler
        results_list=pool.map(cmd_executor,cmds_list)
        if (sum(results_list)) != 0:
            print >> sys.stderr, "At least one parallelized bcftools query failed\n"
        print >> sys.stderr, "Removing temporary .txt files for chunking\n"
        results_list=pool.map(cmd_executor,rm_list)
            
    #no chunking
    elif chunk==0:
        outName=".".join([out,"gen"])
        print >> sys.stderr, "Performing qctool query to pull markers from BGEN %s and write to %s\n" % (bgen,outName) 
        subprocess.call([qctool,"-g",bgen,"-incl-positions",tmp.name,"-og",outName])
    return 0

def bgen_sample(sf,out):
    outName=".".join([out,"sample"])
    print >> sys.stderr,"Reading in sample names from .sample file belonging to bgen and writing to %s\n" % outName
    #open new sample file for writing
    sampleFile=open(outName,"w")
    
    #open sample file from bgen
    command=open_zip(sf)
    with command as f:
        for line in f:
            ls=line.rstrip()
            lineList=ls.split()
            sampleFile.write("\t".join([lineList[0],lineList[1]])) #UKBB .sample file is ID_1, ID_2
            sampleFile.write("\n")
    f.close()
    sampleFile.close()
    return(0)
    
#########################
########## MAIN #########
#########################

def main():
    
    #get arguments
    args = get_settings()

    if int(args.chunk) > 1000:
        print >> sys.stderr, "UNIX's split command will exhaust suffixes. Recode this line or choose a smaller chunk size\n"
    
    #VCF data
    if args.bgen == None:
        
        #makes bed file of markers from weight file
        tmp_obj,counter=readWeights(args.file,args.chr, args.pos,args.chr_num)

        #make sample file from VCF
        readSamples(args.vcf,args.output)
    
        #writes .dose file from VCF, extracts markers from weight file
        #may chunk each .dose file into X chunks
        callQuery(args.vcf,tmp_obj,args.output,int(args.chunk),counter,args.bcftools,args.split,args.cpu)

    #BGEN data 
    elif  args.vcf == None:

        #check for sample file
        if args.sample == None:
            print >> sys.stderr, "Must provide sample file for .bgen\n"
        else:
            bgen_sample(args.sample,args.output)
        #To do:write out sample file with just FID and IID 

        #makes text file of markers from weight file
        tmp_obj,counter=readWeightsForBgen(args.file,args.chr,args.pos,args.chr_num)

        #To do: code split by hand, suffixes exhausted, make warning about split/cbuunk
        #writes .gen from .bgen, extracts only markers from weight file
        #may chunk each .bgen file into X chunks
        bgenToGen(args.bgen,tmp_obj,args.output,int(args.chunk),counter,args.qctool,args.split,args.cpu)

    else: 
        print>> sys.stderr,"Must provide --vcf or --bgen but not both\n"

    
main()
