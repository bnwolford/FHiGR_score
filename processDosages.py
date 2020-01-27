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
    
#This python script can be used to calculate genetic risk scores from dosages in VCF file and a score file with weights
#Takes a fairly generalized weights file
# - Format your weight file as chr:pos, effect allele, weight OR provide 0-based column numbers for chr, pos, chr:pos, effect allele, weight
# - Provide number of lines that have ## or text header to ignore when reading in file
###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser(description='Enter files to use for PRS calculation')
    parser.add_argument('-r','--region_file',help="Region file already made for tabix from processWeights.py",type=str)
    parser.add_argument('-w', '--weight_file',help="Must be sorted by position. Columns and headers are customizable with arugments.")
    parser.add_argument("-cc","--chrom_col",help="0-based column with chromosome in file",type=int)
    parser.add_argument("-pc","--pos_col",help="0-based column with end position of variant in weight file",type=int)
    parser.add_argument("-dc","--coord_col",help="0-based column with chromosome:position:ref:alt of variant in weight file", type=int)
    parser.add_argument("-ec","--ea_col",help="0-based column with effect allele in weight file",type=int,default=1)
    parser.add_argument("-rc","--ref_col",help="0-based column with reference allele in weight file",type=int)
    parser.add_argument("-ac","--alt_col",help="0-based column with alternate allele in weight file",type=int)
    parser.add_argument("-wc","--weight_col",help="0-based column with weight",type=int,default=2)
    parser.add_argument("-l","--header_lines",help="Number of header lines in weight file to skip",default=16)
    parser.add_argument('-v', '--single_vcf',help="Path to VCF. Expects GT:DS")
    parser.add_argument("-c","--vcf_chrom",help="Provide a chromosome number  of VCF of multi VCFs for efficiency",type=int)
    parser.add_argument('-i', '--id_file', help="File with sampleIDs in the same order as the VCF",type=str)
    parser.add_argument('-o', '--output_prefix',type=str,default="results")
    parser.add_argument("--split",help="split path",type=str,default="/usr/bin/split")
    parser.add_argument("--tabix",help="bcftools path",type=str,default="/usr/local/bin/tabix")
    parser.add_argument("-u","--cpu",help="Number of CPU cores to utilize for multiprocessing",default=8)
    args=parser.parse_args()

    ## catches people using X chromosome VCF but only if doing on a per chromosome basis
    if str(args.vcf_chrom)=="X":
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
    if f is None:
         print >> sys.stderr, "File was not given as parameter\n"
    elif ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderr, "Opening gzipped file %s\n" % f
    elif f == "-":
        command=sys.stdin()
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
    return command
                            

#flatten nested lists
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)
            
#create dictionary of weight per variant
def read_weights(weight_file,chrom,pos,ref,alt,coord,ea,weight,vcf_chrom):
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
                if chrom is not None: #because of argument check function we can trust this means we are making our own coordinate with chrom, pos, ref, alt
                    if vcf_chrom is not None: #only save info for chromosome of interes
                        if int(lineList[chrom])==vcf_chrom:
                            coordinate=":".join([str(lineList[chrom]),str(lineList[pos]),str(lineList[ref]),str(lineList[alt])])
                            weight_dict[coordinate]=(lineList[ea],float(lineList[weight]))
                    else:
                        coordinate=":".join([str(lineList[chrom]),str(lineList[pos]),str(lineList[ref]),str(lineList[alt])])
                        weight_dict[coordinate]=(lineList[ea],float(lineList[weight]))
                elif coord is not None:
                    if len(lineList[coord].split(":"))!=4:
                        sys.exit("Coordinate must have 4 components chr:pos:ref:alt\n")
                    if vcf_chrom is not None:
                        if lineList[coord].split(":")[0]==vcf_chrom: #only save info for chromosome of interest
                           weight_dict[lineList[coord]]=(lineList[ea],float(lineList[weight]))
                else:  #dont need an else condition because of argument check
                    continue
    return weight_dict

def getDosage(region_file,tabix_path,vcf,cpu,weight_dict,sample_id,output):
    print >> sys.stderr, "Calling tabix on %s to subset markers from  %s\n" % (vcf,region_file)
    #cmd=[tabix_path, '-R',region_file , vcf]
    cmd=["/usr/local/bin/bcftools","query","-R",region_file,vcf,"-f","%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%DS\t]\n"]
    marker_count=0
    test_results=[]
    try:
        f = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1)
        with f.stdout:
            for line in iter(f.stdout.readline,b''):
                ls=line.split()
                ls[-1]=ls[-1].rstrip()
                #print(ls)
                if ls[0] != "#":
                    coord=ls[0] + ":" + ls[1] + ":" + ls[3] + ":" + ls[4]
                    alt_coord=ls[0] + ":" + ls[1] + ":" + ls[4] + ":" + ls[3] #flip ref and alt in case weight file is in that order
                
                    if coord in weight_dict:
                        #test_results.append(flatten((weight_dict[coord][0], weight_dict[coord][1], ls[0:5], [value.split(":")[1] for value in ls[9:]])))
                        test_results.append(flatten((weight_dict[coord][0],weight_dict[coord][1], ls))) 
                    elif alt_coord in weight_dict:
                        #test_results.append(flatten((weight_dict[alt_coord][0], weight_dict[alt_coord][1], ls[0:5],[value.split(":")[1] for value in ls[9:]])))
                        test_results.append(flatten((weight_dict[alt_coord][0], weight_dict[alt_coord][1],ls)))
        test_results=np.asarray(test_results)
        #print(test_results)
        #print(len(test_results))
        #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
        #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....
        marker_count+=len(test_results)
        if (np.shape(test_results)[0] == 1):  #if just 1 row (i.e. 1 marker)
            #print(np.shape(test_results))
            test_results = np.vstack((test_results, np.zeros(np.shape(test_results))))
        elif (np.shape(test_results)[0]==0): #if no rows (i.e. 0 markers)
            #assumes length of sample ID is samples we are pulling from VCF when test_results is successful
            test_results = np.vstack((np.zeros(len(sample_id)),np.zeros(len(sample_id))))
            #Assumes DS in VCF is in terms of the alternate allele
            #Where effect allele from risk score formula matches alternative allele (column 6), multiply directly
        matching_effect_allele = test_results[np.where(test_results[:,0] == test_results[:,6])]
        #[:, np.newaxis] this is needed to do the multipication element wise (first column * all dosages in row)
        matching_effect_allele = matching_effect_allele[:,1].astype(float)[:, np.newaxis] * matching_effect_allele[:,7:].astype(float)

        #Where effect allele matches reference allele, take 2-dosage, then multiply (so flip dosage to be for alternative allele)
        # To Do: check before subtracting from 2 to flip it because currently assumes autosome VCF only
        matching_reference_allele = test_results[np.where(test_results[:,0] == test_results[:,5])]
        matching_reference_allele = matching_reference_allele[:,1].astype(float)[:, np.newaxis] * (2 - matching_reference_allele[:,7:].astype(float))

        #Sum down columns
        dosage_scores_sum = np.sum(matching_reference_allele, axis=0) + np.sum(matching_effect_allele, axis=0)
        sample_score_dict = {sample_id[x]: score for x, score in enumerate(dosage_scores_sum)}
        sample_score_dict["count"]=marker_count #record number of markers from weights that are present in the VCF

        f.wait()
        
    except KeyboardInterrupt:
        print >> sys.stderr, "Caught KeyboardInterrupt, terminating multiprocesses\n"
        sys.exit("Exiting program\n")

    #write output file
    outputname=output + "_" + "scores.txt"
    print >> sys.stderr, "Writing output file %s \n" % (outputname)
    with open(outputname, 'w') as out:
        out.write("%s\t%s\n" % ("individual", "score"))
        for x in range(len(sample_id)):
            out.write("%s\t%.8f\n" % (sample_id[x], sample_score_dict[sample_id[x]]))

##handle situation if there are empty weights to avoid weird key error problem
def empty_weights(sample_id,output):
    outputname=output + "_" + "scores.txt"
    print >> sys.stderr, "Writing output file %s \n" % (outputname)
    with open(outputname, 'w') as out:
        out.write("%s\t%s\n" % ("individual", "score"))
        for x in range(len(sample_id)):
            out.write("%s\t%.8f\n" % (sample_id[x], 0))
    out.close()
    sys.exit()
    
#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args=get_settings()

    if args.vcf_chrom is not None:
         print >> sys.stderr, "Only considering chromosome %s\n"  %  args.vcf_chrom
        
    #open ID file
    #Assumes order in VCF is the same across everything provided
    with open(args.id_file) as f:
        sample_id = [line.rstrip() for line in f]

    #create dictionary of weights per variant
    weight_dict=read_weights(args.weight_file,args.chrom_col,args.pos_col,args.ref_col,args.alt_col,args.coord_col,args.ea_col,args.weight_col,args.vcf_chrom)
    if len(weight_dict.keys())==0:
        print >> sys.stderr, "No markers in this weight file %s\n" % args.weight_file
        empty_weights(sample_id,args.output_prefix)
    
    #Calculate weighted dosages per individual, Make sure to check allele, print output 
    getDosage(args.region_file,args.tabix,args.single_vcf,args.cpu,weight_dict,sample_id,args.output_prefix)
        
        
##### Call main 
if __name__ == "__main__":
        main()
