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

print(sys.version)

#This python script can be used to calculate genetic risk scores from dosages in VCF file and a score file with weights
#Takes a fairly generalized weights file
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
    parser.add_argument("-k","--chunk",help="Number of markers from weight file to run at a time",default=1000,type=int)
    parser.add_argument('-m', '--multi_vcf', nargs='*')
    parser.add_argument('-v', '--single_vcf')
    parser.add_argument("-c","--vcf_chrom",help="Provide a chromosome number  of VCF of multi VCFs for efficiency",type=int)
    parser.add_argument('-i', '--id_file', default="/net/fantasia/home/sarahgra/Collaborator_projects/PRS_prediction_Cristen/MGI_sample_IDs")
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

    ## check VCFs
    if args.single_vcf is None and args.multi_vcf is None:
         sys.exit("Need a path to VCF with * for multiple VCF if needed\n")
            
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


#write out regions file to use with tabix, gets coordinates from weights file
def make_regions_file(weight_dict, output_name,chunk):
    number_markers=len(weight_dict.keys())
    if int(math.ceil(number_markers / chunk)) == 1: #chunk value >= than markers
        num_files=1
        chunk=number_markers
        sys.stderr.write("--chunk parameter is greater than or equal to the number of markers. May want to consider a more appropriate chunk parameter. Writing one marker per region file\n")
    else:
        num_files=int(math.ceil(number_markers / chunk))
    tmpFileList=[]
    sys.stderr.write("Writing %d temporary files for marker regions\n" % num_files)

    for i in range(num_files):
        regions = NamedTemporaryFile(delete=False)
        tmpFileList.append(regions.name)
        with open(regions.name, 'w') as tmp:
            for coord in weight_dict.keys()[i*chunk:(i+1)*chunk]:
                chrom,pos,ref,alt=coord.split(":")
                tmp.write(chrom+"\t"+pos + "\n")
    print(tmpFileList)
    return tmpFileList,number_markers

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def getDosage(tmpFileNames,tabix_path,vcf_list,cpu,weight_dict,sample_id,output):
    cmd_list=[] #initialize command list
    for vcf in vcf_list: #loop over vcf(s)
        for region in tmpFileNames: #loop over chunked region files
            cmd_list.append([tabix_path, '-R',region , vcf]) #make list of lists for commands
    pool = mp.Pool(cpu,init_worker) #use user defined number of cpus, user should also specify this value for job scheduler
    pfunc=partial(process_function,weight_dict=weight_dict,sample_id=sample_id) #set weight_dict as a standing variable for the process_function
    results_list=pool.map_async(pfunc,cmd_list)
    try:
        print >> sys.stderr, "Using multiprocessing pool functionality with %d cpus and %d processes to perform tabix on %d VCF(s)" % (cpu,len(cmd_list),len(vcf_list))
        time.sleep(10)
    except KeyboardInterrupt:
        print >> sys.stderr, "Caught KeyboardInterrupt, terminating multiprocesses\n"
        pool.terminate()
        pool.join()
        sys.exit("Exiting program\n")
    else:
        pool.close()
        pool.join()
        print >> sys.stderr, "Normal termination of multiprocesses upon completion\n"
    #merge all the dosages 
    sys.stderr.write("Merging per sample scores across chunked regions and VCF(s)\n")
    c=Counter() #initialize counter
    for dictionaries in results_list.get():
        c.update(dictionaries) #sum across all the dictionaries 
    print >> sys.stderr, "%d variants were in the region file(s) and %d were ultimately found in the VCF(s)"  % (len(weight_dict),c["count"])
    #write output file
    outputname=output + "_" + "scores.txt"
    print >> sys.stderr, "Writing output file %s \n" % (outputname)
    with open(outputname, 'w') as out:
        out.write("%s\t%s\n" % ("individual", "score"))
        for x in range(len(sample_id)):
            out.write("%s\t%.8f\n" % (sample_id[x], c[sample_id[x]]))

        
#function to multiprocess
def process_function(cmd,weight_dict,sample_id):
    #Create dictionary to keep track of total scores per person, set initial value to zero
    sample_score_dict = {x:0 for x in sample_id}
    marker_count=0
    f = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    f = f.communicate()[0]
    test_results=[]
    if f!="": #if tabix query finds at  least one of  the  risk variants in that chunk
        for line in f.rstrip().split("\n"):
            if line.split()[0][0] != "#":
                coord=line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4]
                alt_coord=line.split()[0] + ":" + line.split()[1] + ":" + line.split()[4] + ":" + line.split()[3] #flip ref and alt in case weight file is in that order 
                if coord in weight_dict:
                    test_results.append(flatten((weight_dict[coord][0], weight_dict[coord][1], line.rstrip().split()[0:5], [value.split(":")[1] for value in line.rstrip().split()[9:]])))
                elif alt_coord in weight_dict:
                     test_results.append(flatten((weight_dict[alt_coord][0], weight_dict[alt_coord][1], line.rstrip().split()[0:5],[value.split(":")[1] for value in line.rstrip().split()[9:]])))

       # test_results = [flatten((weight_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][0], weight_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][1], line.rstrip().split()[0:5], [value.split(":")[1] for value in line.rstrip().split()[9:]])) for line in f.rstrip().split("\n") if line.split()[0][0] != "#" if (line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4]) in weight_dict]
    test_results=np.asarray(test_results)

    #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
    #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....
    marker_count+=len(test_results)
    if (np.shape(test_results)[0] == 1):  #if just 1 row (i.e. 1 marker)
        #print(np.shape(test_results))
        test_results = np.vstack((test_results, np.zeros(np.shape(test_results))))
        #print(np.shape(test_results))
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
    for x in range(len(dosage_scores_sum)):
        sample_score_dict[sample_id[x]] = sample_score_dict[sample_id[x]] + dosage_scores_sum[x]
        
    sample_score_dict["count"]=marker_count #record number of markers from weights that are present in the VCF 

    return(sample_score_dict)


#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args=get_settings()

    ## Handle multi or single VCF
    if args.multi_vcf is not None:
        file_list=args.multi_vcf
    elif args.single_vcf is not None:
        file_list = [args.single_vcf]

    if args.vcf_chrom is not None:
         print >> sys.stderr, "Region file(s) only contain chromosome %s\n"  %  args.vcf_chrom
        
    #open ID file
    #Assumes order in VCF is the same across everything provided
    with open(args.id_file) as f:
        sample_id = [line.rstrip() for line in f]

    #create dictionary of weights per variant
    weight_dict=read_weights(args.weight_file,args.chrom_col,args.pos_col,args.ref_col,args.alt_col,args.coord_col,args.ea_col,args.weight_col,args.vcf_chrom)
    
    #Write out regions file for tabix using dictionary and chunk parameter 
    #regions_output_name = "Regions_" + str(args.vcf_chrom) + "_" +  args.weight_file.split("/")[-1]
    tmpFileNames,num_markers=make_regions_file(weight_dict,args.output_prefix,args.chunk)

    ## Handle multi or single VCF
    if args.multi_vcf is not None:
        vcf_list=args.multi_vcf
    elif args.single_vcf is not None:
        vcf_list = [args.single_vcf]

    #Calculate weighted dosages per individual, Make sure to check allele, print output 
    getDosage(tmpFileNames,args.tabix,vcf_list,args.cpu,weight_dict,sample_id,args.output_prefix)
        
        
    #record the number of markers actually found and included because risk marker list may differ from variants in data of interest
 #   variant_count=0
    #loop over VCFs
  #  for file_x in file_list:
   #     print >> sys.stderr, "Now reading in: %s" % file_x
    #    f = Popen(['tabix', '-R', regions_output_name, file_x], stdout=PIPE)
     #   f = f.communicate()[0]

      #  print(sample_score_dict)
       # print(weight_file_dict)
       # print(f)
        #If the tabix command finds one of the risk variants in that chunk:
        #To do: record the number of markers actually found and included because risk marker list may differ from variants in data of interest
        #if f != '':
         #   variant_count+=1 #add to marker count
            #This will return the effect allele and effect for each line if variant is in the risk score and not a comment line, and then the dosage line
          #  test_results = np.asarray([flatten((weight_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][0], weight_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][1], line.rstrip().split()[0:5], [value.split(":")[1] for value in line.rstrip().split()[9:]])) for line in f.rstrip().split("\n") if line.split()[0][0] != "#" if (line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4]) in weight_file_dict])
            #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
            #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....
           # print(test_results)
            #if (np.shape(test_results)[0] == 1):
            #    test_results = np.vstack(test_results, np.zeros(np.shape(test_results))
            #Assumes DS in VCF is in terms of the alternate allele 
            #Where effect allele from risk score formula matches alternative allele, multiply directly
            #matching_effect_allele = test_results[np.where(test_results[:,0] == test_results[:,6])]
            #[:, np.newaxis] this is needed to do the multipication element wise (first column * all dosages in row)
            #matching_effect_allele = matching_effect_allele[:,1].astype(float)[:, np.newaxis] * matching_effect_allele[:,7:].astype(float)
            
            #Where effect allele matches reference allele, take 2-dosage, then multiply (so flip dosage to be for alternative allele)
            # To Do: check before subtracting from 2 to flip it because currently assumes autosome VCF only
            #matching_reference_allele = test_results[np.where(test_results[:,0] == test_results[:,5])]
            #matching_reference_allele = matching_reference_allele[:,1].astype(float)[:, np.newaxis] * (2 - matching_reference_allele[:,7:].astype(float))
            
            #Sum down columns
            #dosage_scores_sum = np.sum(matching_reference_allele, axis=0) + np.sum(matching_effect_allele, axis=0)
            #for x in range(len(dosage_scores_sum)):
            #    sample_score_dict[sample_id[x]] = sample_score_dict[sample_id[x]] + dosage_scores_sum[x]

#    print >> sys.stderr, "%d variants were in the region file and %d were ultimately found in the VCF(s)"  % (region_count,variant_count)
    
 #   with open(args.output_file, 'w') as out:
  #      out.write("%s\t%s\n" % ("individual", "score"))
   #     for x in range(len(sample_id)):
    #        out.write("%s\t%.8f\n" % (sample_id[x], sample_score_dict[sample_id[x]]))

    #delete temporary region files
    for filename in tmpFileNames:
        os.remove(filename)
    
##### Call main 
if __name__ == "__main__":
        main()
