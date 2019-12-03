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
import numpy as np
import pandas as pd
import scipy.stats as ss

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--plink_file",help="File with header and FID and IID of the samples you expect to calcualte GRS for",type=str)
    parser.add_argument("-s","--sample_file",help="File with no header and one column of samples that you expect to calcualte GRS for",type=str)
    parser.add_argument("-c","--config",help="Config file with new line delimited list of all the file names to merge",required=True,type=str)
    parser.add_argument("-o","--output",help="Output prefix. Will have .txt appended",required=True,type=str)
    parser.add_argument("-i","--invNorm",help="Will print additional column with inverse normalized score",action='store_true')
    parser.add_argument("--id_column",help="0-based column with ID in the score file. IDs should match those from plink_file or sample_file. If not provided, assumes FID, IID, score",type=int)
    parser.add_argument("--score_column",help="0-based column with score in the score file. If not provided, assumes FID, IID, score",type=int)
    parser.add_argument("--header",help="Flag if score file has header.",action='store_true',default=False)
    parser.set_defaults(chrom=True,invNorm=False)
    args=parser.parse_args()
    return args

###############################
######## SUB ROUTINES  #########
###############################

def open_zip(f):
    exists=os.path.isfile(f)
    if exists:
        if ".gz" in f:
            command=gzip.open(f,"rt")
            print >> sys.stderrr, "Opening gzipped file %s\n" % f
        elif f == "-":
            sys.exit("Cannot read file from stdin\n")
        else:
            command=open(f,"rt")
            print >> sys.stderr, "Opening file %s\n" % f
        return command
    else:
        print >> sys.stderr, "%s does not exist\n" %f
        return 0

def check(config):
    """
    Check that all the .gen files in the config file exist 
    """
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
            elif ok:
                if os.path.getsize(fn) == 0:
                    print >> sys.stderr, "%s is an empty file\n" % fn
    return(file_list)


def read_plink(sample):
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
                lineList=ls.split("\t") #assume tab seperator in .sample file from VCFtoDOSEforGRS.py
                sample_list.append(".".join(lineList[0:2])) #put FID.IID values in a list
    return sample_list


def read_sample(sample):
    """
    Reads a sample file without header and just one column of sample IDs that you expect to have GRS values for in each chunk.
    """
    sample_list=[]
    command=open_zip(sample)
    with command as f:
        for line in f:
            ls=line.rstrip()
            sample_list.append(ls)
    return sample_list

def merge(fl,sl):
    """
    Looks at all .gen files matching the string given. Returns dictionary with information and potentially flags any missing chromosomes or chunks.
    """
    ddict={} #initialize dictionary
    for sample in sl:
        ddict[sample]=[] #initialize list

    for gen in fl:
        command=open_zip(gen)
        count=0
        if command!=0: #if file exists
            with command as f:
                for line in f:
                    if count==0: #skip 1 line header
                        count+=1
                        next
                    else:
                        ls=line.rstrip() #expects FID, IID, weighted sum from DOSEtoGRS.py
                        lineList=ls.split(" ")
                        ids=".".join(lineList[0:2])
                        ddict[ids].append(lineList[2])
                        
    return(ddict)


def merge_custom(fl,sl,id_col,score_col,header):
    """
    Looks at all files matching the string given. Returns dictionary with the id and score information.
    """
    ddict={} #initialize dictionary
    for sample in sl:
        ddict[sample]=[] #initialize list

    for score_file in fl:
        command=open_zip(score_file)
        count=0
        if command!=0:
            with command as f:
                for line in f:
                    if header is True and count==0: #skip 1 line header
                        count+=1
                        next
                    else:
                        ls=line.rstrip()
                        lineList=ls.split()
                        sample_id=lineList[id_col]
                        ddict[sample_id].append(lineList[score_col])
    return(ddict)
                        



def output(o,d,list_size,inorm):

    outname=".".join([o,"txt"])
    out_file=open(outname,"w")
    if inorm:
        #initialize lists for inorm
        GRS_list=[]
        id_list=[]
        out_file.write("\t".join(["IID","FID","GPS","invNorm_GPS\n"]))
    else:
        out_file.write("\t".join(["IID","FID","GPS\n"]))

    #sum across the nested dictionary to get 1 value per ID
    for ids in d.keys():
        if len(d[ids]) != list_size: #check that all sub sums are represented 
            print >> sys.stderr, "%s does not have expected number (%d) of sub-chunks to sum\n"  % (ids,list_size)
        GRS=sum(float(sub_sum) for sub_sum in d[ids])
        #if we want to inverse normalize just save data to ordered lists
        if inorm:
            GRS_list.append(GRS)
            id_list.append(ids)
        #otherwise write to file
        else:
            if "." in ids[0]: #if this script created FID.IID
                id_1,id_2=ids.split(".")
                out_file.write("\t".join([id_1,id_2,str(GRS)])+"\n")
            else:
                 out_file.write("\t".join([ids,ids,str(GRS)])+"\n")
    #now we inverse normalize and write out data using lists 
    if inorm:
        trans_GRS_list=rank_INT(GRS_list)
        for i in range(len(trans_GRS_list)):
            if "." in ids[0]: #if this script created FID.IID
                id_1,id_2=id_list[i].split(".")
                out_file.write("\t".join([id_1,id_2,str(GRS_list[i]),str(trans_GRS_list[i])])+"\n")
            else:
                out_file.write("\t".join([id_list[i],id_list[i],str(GRS_list[i]),str(trans_GRS_list[i])])+"\n")
    return 0


#Code adapted from https://github.com/edm1/rank-based-INT/blob/master/rank_based_inverse_normal_transformation.py
def rank_INT(GRS_list,c=3.0/8,stochastic=True):
    """
    Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.
        Args:
            param1 (pandas.Series):   Series of values to transform
            param2 (Optional[float]): Constand parameter (Bloms constant)
            param3 (Optional[bool]):  Whether to randomise rank of ties
        
        Returns:
            pandas.Series
    """
    #set seed
    np.random.seed(1234)
    
    GRS_series=pd.Series(GRS_list)
    orig_index=GRS_series.index
    GRS_series=GRS_series.loc[~pd.isnull(GRS_series)] #drop NAN

    #Get ranks
    #Ties are determined by their position in the series 
    if stochastic==True:
        GRS_series=GRS_series.loc[np.random.permutation(GRS_series.index)]
        rank=ss.rankdata(GRS_series,method="ordinal")
    #Ties are averaged
    else:
        rank=ss.rankdata(GRS_series,method="average")

    #Convert numpy array back to series
    rank=pd.Series(rank,index=GRS_series.index)

    #Convert rank to normal distribution
    transformed=rank.apply(rank_to_normal,c=c,n=len(rank))

    return transformed[orig_index]


def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2*c + 1)
    return ss.norm.ppf(x)
    
    

#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args = get_settings()

    #check that all the chunk/chr files exist
    file_list=check(args.config)

    #read in samples that are expected from either a PLINK type file or sample list like from bcftools query -l
    if (args.sample_file is None and args.plink_file is None) or (args.sample_file is not None and args.plink_file is not None):
        sys.exit("--sample_file or --plink_file is required")
    elif args.plink_file is not None:
        sample_list=read_plink(args.plink_file)
    elif args.sample_file is not None:
        sample_list=read_sample(args.sample_file)

    #merge data assuming FID, IID, score or a custom score results file 
    if (args.score_column is not None and args.id_column is not None):
        data=merge_custom(file_list,sample_list,args.id_column,args.score_column,args.header) #provide custom columns
    else:
        data=merge(file_list,sample_list) #assumes score results file with FID, IID, score
    

    #write output
    output(args.output,data,len(file_list),args.invNorm)

    

    
main()
