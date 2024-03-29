#!/usr/bin/env python

#===============================================================================
# Copyright (c) 2019 Brooke Wolford
# Revised from calculate_risk_scores.py from Dr. Sarah Graham
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#==============================================================================

############################
##### IMPORT MODULES #######
###########################

from __future__ import division
import argparse
import subprocess
import gzip
import io
import numpy as np
from glob import glob
import sys
from tempfile import NamedTemporaryFile
import math
from collections import OrderedDict
from collections import Counter 
import os
import multiprocessing as mp
#import pysam
from functools import partial
import time
import signal
import datetime
now = datetime.datetime.now()
print ("Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))
print(sys.version)
for module in sys.modules:
    try:
        print(module,sys.modules[module].__version__)
    except:
        pass
    
#This python script can be used to process weights files into region .bed files for tabix 
# Takes a fairly generalized weights file
# - Format your weight file as chr:pos, effect allele, weight OR provide 0-based column numbers for chr, pos, chr:pos, effect allele, weight
# - Provide number of lines that have ## or text header to ignore when reading in file
###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser(description='Enter files to use for PRS calculation')
    parser.add_argument('-w', '--weight_file',help="Must be sorted by position. Columns and headers are customizable with arguments. Genom wide")
    parser.add_argument("-f",'--chunk_file',default="/net/dumbo/home/larsf/UKB500/VCF/Chunks_250k.txt",help="File with relevant information about chunked VCFs so we can make weight files and region files corresponding to the chunked VCFs. First 3 columns must be chr, start, end.")
    parser.add_argument("-xc","--prefix_col",help="0-based column with prefix for the chunk from the --chunk_file",default=7)
    parser.add_argument("-cc","--chrom_col",help="0-based column with chromosome in file",type=int)
    parser.add_argument("-pc","--pos_col",help="0-based column with end position of variant in weight file",type=int)
    parser.add_argument("-dc","--coord_col",help="0-based column with chromosome:position:ref:alt of variant in weight file", type=int)
    parser.add_argument("-ec","--ea_col",help="0-based column with effect allele in weight file",type=int,default=1)
    parser.add_argument("-rc","--ref_col",help="0-based column with reference allele in weight file",type=int)
    parser.add_argument("-ac","--alt_col",help="0-based column with alternate allele in weight file",type=int)
    parser.add_argument("-wc","--weight_col",help="0-based column with weight",type=int,default=2)
    parser.add_argument("-l","--header_lines",help="Number of header lines in weight file to skip",type=int)
    parser.add_argument("-n","--num_chunk",help="Number of markers from weight file to run at a time",default=1000,type=int)
    parser.add_argument('-o', '--output_prefix',type=str,default="results")
    parser.add_argument("--split",help="split path",type=str,default="/usr/bin/split")
    args=parser.parse_args()

    ## check if marker information is adequately provided
    if args.coord_col is None: #no coordinate offered
        check_list=[]
        coordinate_columns=[args.chrom_col,args.pos_col,args.ref_col,args.alt_col]
        for f in coordinate_columns:
            check_list.append(f!=None) #how many are empty
        if sum(check_list)==0:
            sys.exit("Need coordinate column for chr:pos:ref:alt or each of these pieces of information individually\n")
        elif sum(check_list)!=4: 
            sys.exit("Need ALL four columns for chromosome, position, reference, alternate\n")

    print >> sys.stderr, "%s\n"  % args

    return args

###############################
######## FUNCTIONS  #########
###############################

#open files
def open_zip(f):
    if f is None:
                 print >> sys.stderr, "File was not given as parameter\n"
    elif ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderrr, "Opening gzipped file %s\n" % f
    elif f == "-":
        command=sys.stdin()
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
    return command
                            
def read_weights(weight_file,chrom,pos,ref,alt,coord,ea,weight,header_lines):
    """
    Read file with weights into dictionary.
    """
    weight_dict=OrderedDict()
    command=open_zip(weight_file)
    counter=0
    with command as f:
        for line in f:
            counter+=1
            if (header_lines is not None and counter > header_lines) or (header_lines is None and line[0]!="#"): #handle header
                ls=line.rstrip()
                lineList=ls.split() #assumes whitespace delimiter, space or tab
                if chrom is not None: #because of argument check function we can trust this means we are making our own coordinate with chrom, pos, ref, al1t
                    coordinate=":".join([str(lineList[chrom]),str(lineList[pos]),str(lineList[ref]),str(lineList[alt])])
                    weight_dict[coordinate]=lineList
                elif coord is not None:
                    if len(lineList[coord].split(":"))!=4:
                        sys.exit("Coordinate must have 4 components chr:pos:ref:alt\n")
                    weight_dict[lineList[coord]]=(lineList[ea],float(lineList[weight]))
                else: #earlier checks mean we should never hit this
                    sys.exit("Error1\n")

    return weight_dict

def read_chunks(chunk_file,prefix_col):
    """Read in chunk file"""
    chunk_list=[]
    command=open_zip(chunk_file)
    with command as f:
        for line in f:
            ls=line.rstrip()
            if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit, also ignores sex chrom
                lineList=ls.split() #assumes whitespace delimiter, space or tab
                chunk_list.append(":".join([lineList[0],lineList[1],lineList[2],lineList[prefix_col]])) #chr:start:end:prefix

    return(chunk_list)

def  make_regions(weight_dict,chunk_list,output_prefix,num_chunk):
    """ Make a region files using weights and known chunks """
    total=len(weight_dict.keys())
    print(weight_dict.keys())
    if total==0:
        sys.stderr.write("Weights file is empty\n")
        sys.exit()
            #open config file
    config="_".join([output_prefix,"config.txt"])
    cfile=open(config,'w+')
    for j in range(len(chunk_list)):
        mega_chunk=chunk_list[j]
        chunk_chrom,chunk_start,chunk_end,prefix=mega_chunk.split(":")
        vcf="".join(["/net/dumbo/home/larsf/UKB500/VCF/ukb24460_imp_",prefix,"_v3_s486743.vcf.gz"]) #need to generalize
        mini_chunk_counter=0
        region="_".join([output_prefix,prefix,".".join([str(mini_chunk_counter),"regions.txt"])]) #new region
        weight="_".join([output_prefix,prefix,".".join([str(mini_chunk_counter),"weights.txt"])]) #new weight
        cfile.write("\t".join([region,weight,vcf,chunk_chrom])+"\n") #write in config
        rfile=open(region,'w')
        wfile=open(weight,'w')
        marker_counter=0

        for entry in weight_dict.keys(): #this is an ordered dictionary
            print(marker_counter)
            print(mini_chunk_counter)
            print(entry,chunk_chrom,chunk_start,chunk_end)
            chrom,pos,a1,a2=entry.split(":")
            if int(chrom)==int(chunk_chrom) and int(pos)>=int(chunk_start) and int(pos)<=int(chunk_end): #still in mega chunk
                rfile.write("\t".join([chrom,pos]) + "\n")
                wfile.write("\t".join(weight_dict[entry]) + "\n")
                print("Write %s" % entry)
                del weight_dict[entry] #delete from dictionary
                marker_counter+=1
                if marker_counter==num_chunk+1: #end of a mini chunk, start new region and write files within the same mega_chunk
                    rfile.close()
                    wfile.close()
                    marker_counter=0
                    mini_chunk_counter+=1
                    break
            elif (int(chrom)!=int(chunk_chrom) or int(pos) > int(chunk_end)) or marker_counter==total: #end of mega chunk
                print("%s not in %s" % (entry,mega_chunk))
                j+=1
                break

    cfile.close()
    clean_up(config)
    return 0


def clean_up(config):
    print("Clean up directory\n")
    keep=[]
    command=open_zip(config)
    with command as f:
        for line in f:
            ls=line.rstrip()
            lineList=ls.split()
            if os.stat(lineList[0]).st_size==0 | os.stat(lineList[1]).st_size==0:
                os.remove(lineList[0])
                os.remove(lineList[1])
            else:
                keep.append(line)
    f.close()
    os.remove(config)
    cfile=open(config,'w+')
    print("Writing new config\n")
    for entry in keep:
        cfile.write(entry)
    cfile.close()
                
                
                          
            
#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args=get_settings()

    print >> sys.stderr, "Not generalized yet. The VCF path printing to config file is hardcoded\n"
    
    #create dictionary of weights per variant
    weight_dict=read_weights(args.weight_file,args.chrom_col,args.pos_col,args.ref_col,args.alt_col,args.coord_col,args.ea_col,args.weight_col,args.header_lines)

    #chr     start   end     N       chunk   bgen    sample  prefi
    chunk_list=read_chunks(args.chunk_file,args.prefix_col)

    #write regions and config file
    make_regions(weight_dict,chunk_list,args.output_prefix,args.num_chunk)
    
##### Call main 
if __name__ == "__main__":
        main()
