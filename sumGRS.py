#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2019 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# Adapted from Dr. Ida Surakka's Chunk_compiler.R
# Merge weighted sums across chromosomes and chunks to calculate one genome wide polygenic risk score per individual
#=============================================================================

############################
##### IMPORT MODULES #######
###########################
import os
import sys
import argparse
import glob

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--sample_file",help="File with FID and IID of the samples you expect to calcualte GRS for",required=True,type=str)
    parser.add_argument("-c","--config",help="Config file with new line delimited list of all the file names to merge",required=True,type=str)
    parser.add_argument("-o","--output",help="Output prefix. Will have .txt appended",required=True,type=str)
    parser.set_defaults(chrom=True)
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

def check(config):
    file_list=[]
    command=open_zip(config)
    with command as f:
        for line in f:
            ls=line.rstrip()
            file_list.append(ls)
        print >> sys.stderr, "Config lists %s files\n" % len(file_list)
        for fn in file_list:
            ok=os.path.isfile(fn)
        if not ok:
            print >> sys.stderr, "Expected %s to exist but it does not\n" % fn
    return(file_list)

def read_samples(sample):
    """
    Reads a sample file with FID and IID. Assumes header. This is the samples that you expect to have GRS values for in each chunk.
    """
    sample_list=[]
    command=open_zip(sample)
    count=0
    with command as f:
        for line in f:
            if count==0: #skip 1 line header, most likely FID IID 
                count+=1
                next
            else:
                ls=line.rstrip()
                sample_list.append(".".join(ls)) #put FID.IID values in a list
    return sample_list
    
def merge(fl,sl):
    """
    Looks at all .gen files matching the string given. Returns dictionary with information and potentially flags any missing chromosomes or chunks.
    """
    ddict={}
    for gen in fl:
        command=open_zip(gen)
        with command as f:
            for line in f:
                ls=line.rstrip() #expects FID, IID, weighted sum from DOSEtoGRS.py
                ids=".".join(ls[2:])
                ddict[ids][gen]=ls[3]
        print ddict
    #check that IDs match those in samples and that each ID has X values summed to create the total
    
    return ddict

def output(o,d):

    #sum across the nested dictionary to get 1 value per ID
    for ids in d.keys():
        for GRS in d[ids].keys():
            print GRS
            
    outname=".".join([o,"txt"])
    out_file=open(outname,"w")
    out_file.write("test")
    return 0

    
#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args = get_settings()

    #check that all the chunk/chr files exist
    file_list=check(args.config)

    #read in samples that are expected
    sample_list=read_samples(args.sample_file)

    #merge data
    data=merge(file_list,sample_list)

    #write output
    output(args.output,data)

    

    
main()
