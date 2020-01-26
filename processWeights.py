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
    parser.add_argument('-w', '--weight_file',help="Must be sorted by position. Columns and headers are customizable with arugments.")
    parser.add_argument("-cc","--chrom_col",help="0-based column with chromosome in file",type=int)
    parser.add_argument("-pc","--pos_col",help="0-based column with end position of variant in weight file",type=int)
    parser.add_argument("-dc","--coord_col",help="0-based column with chromosome:position:ref:alt of variant in weight file", type=int)
    parser.add_argument("-ec","--ea_col",help="0-based column with effect allele in weight file",type=int,default=1)
    parser.add_argument("-rc","--ref_col",help="0-based column with reference allele in weight file",type=int)
    parser.add_argument("-ac","--alt_col",help="0-based column with alternate allele in weight file",type=int)
    parser.add_argument("-wc","--weight_col",help="0-based column with weight",type=int,default=2)
    parser.add_argument("-l","--header_lines",help="Number of header lines in weight file to skip",default=16)
    parser.add_argument("-k","--chunk",help="Split weights file into -n number of markers. If this is not used and many markers are in -w, the process is quite memory intensive",action="store_true")
    parser.add_argument("-n","--num_chunk",help="Number of markers from weight file to run at a time",default=1000,type=int)
    parser.add_argument("-c","--chrom",help="Provide a chromosome number", type=int)
    parser.add_argument('-o', '--output_prefix',type=str,default="results")
    parser.add_argument("--split",help="split path",type=str,default="/usr/bin/split")
    args=parser.parse_args()

    ## catches people using X chromosome VCF but only if doing on a per chromosome basis
    if str(args.chrom)=="X":
        sys.exit("This script currently only handles autosomes\n")
    
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
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderrr, "Opening gzipped file %s\n" % f
    elif f == "-":
        command=sys.stdin()
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
    return command
                            
def read_weights(weight_file,chrom,pos,ref,alt,coord,ea,weight,chrom_num):
    """
    Read file with weights into dictionary and regions files for tabix.
    """
    weight_dict=OrderedDict()
    command=open_zip(weight_file)
    counter=0
    with command as f:
        for line in f:
            ls=line.rstrip()
            if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit
                lineList=ls.split() #assumes whitespace delimiter, space or tab
                if chrom is not None: #because of argument check function we can trust this means we are making our own coordinate with chrom, pos, ref, al1t
                    if chrom_num is not None: #only save info for chromosome of interes
                        if int(lineList[chrom])==chrom_num:
                            coordinate=":".join([str(lineList[chrom]),str(lineList[pos]),str(lineList[ref]),str(lineList[alt])])
                            weight_dict[coordinate]=(lineList[ea],float(lineList[weight]))
                    else:
                        coordinate=":".join([str(lineList[chrom]),str(lineList[pos]),str(lineList[ref]),str(lineList[alt])])
                        weight_dict[coordinate]=(lineList[ea],float(lineList[weight]))
                elif coord is not None:
                    if len(lineList[coord].split(":"))!=4:
                        sys.exit("Coordinate must have 4 components chr:pos:ref:alt\n")
                    if chrom_num is not None:
                        if lineList[coord].split(":")[0]==chrom_num: #only save info for chromosome of interest
                           weight_dict[lineList[coord]]=(lineList[ea],float(lineList[weight]))
                else:  #dont need an else condition because of argument check
                    continue
    return weight_dict


#write out regions file to use with tabix, gets coordinates from weights file
def make_regions_file(weight_dict, output_name,chunk,chunkTF,chrom_num):
    number_markers=len(weight_dict.keys())
    if chrom_num is None:
        chrom_num="all"
    if chunkTF is False:
        num_files=1
        chunk=number_markers
        sys.stderr.write("--chunk flag should be provided if weight markers are to be chunked into multiple tabix processes. --num_chunk is being ignored if provided. \n")
    elif int(math.ceil(number_markers / chunk)) == 1: #chunk value >= than markers
        num_files=1
        chunk=number_markers
        sys.stderr.write("--num_chunk parameter is greater than or equal to the number of markers. May want to consider a more appropriate parameter.\n")
    else:
        num_files=int(math.ceil(number_markers / chunk))
    sys.stderr.write("Writing %d files for marker regions\n" % num_files)
    file_list=[]
    for i in range(num_files):
        output="_".join(["".join([output_name,"chrom",str(chrom_num)]),"".join(["chunk",str(i)])])
        file_list.append(output)
        with open(output, 'w') as ofile:
            for coord in weight_dict.keys()[i*chunk:(i+1)*chunk]:
                chrom,pos,ref,alt=coord.split(":")
                ofile.write(chrom+"\t"+pos + "\n")
    return file_list,number_markers



#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args=get_settings()


    if args.chrom is not None:
         print >> sys.stderr, "Region file(s) only contain chromosome %s\n"  %  args.chrom
        
    #create dictionary of weights per variant
    weight_dict=read_weights(args.weight_file,args.chrom_col,args.pos_col,args.ref_col,args.alt_col,args.coord_col,args.ea_col,args.weight_col,args.chrom)
    
    #Write out regions file for tabix using dictionary and chunk parameter 
    #regions_output_name = "Regions_" + str(args.vcf_chrom) + "_" +  args.weight_file.split("/")[-1]
    file_names,num_markers=make_regions_file(weight_dict,args.output_prefix,args.num_chunk,args.chunk,args.chrom)
    print(file_names)
    print(num_markers)
##### Call main 
if __name__ == "__main__":
        main()
