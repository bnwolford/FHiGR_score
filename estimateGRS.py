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
from subprocess import Popen, PIPE
import gzip
import io
import numpy as np
from glob import glob
import sys

#This python script can be used to calculate genetic risk scores from dosages in VCF file and a score file with weights
#Takes a fairly generalized weights file
# - Format your weight file as chr:pos, effect allele, weight OR provide 0-based column numbers for chr, pos, chr:pos, effect allele, weight
# - Provide number of lines that have ## or text header to ignore when reading in file
###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser(description='Enter files to use for PRS calculation')
    parser.add_argument('-w', '--weight_file')
    parser.add_argument("-r","--chrom_col",help="0-based column with chromosome in file",type=int)
    parser.add_argument("-p","--pos_col",help="0-based column with end position of variant in file",type=int)
    parser.add_argument("-d","--coord_col",help="0-based column with chromosome:position of variant in weight file", type=int,default=0)
    parser.add_argument("-a","--ea_col",help="0-based column with effect allele",type=int,default=1)
    parser.add_argument("-t","--weight_col",help="0-based column with weight",type=int,default=2)
    parser.add_argument("-l","--header_lines",help="Number of header lines in weight file to skip",default=16)
    parser.add_argument('-m', '--multi_vcf', nargs='*')
    parser.add_argument('-v', '--single_vcf')
    parser.add_argument("-c","--chrom",help="Provide a chromosome number  of VCF of multi VCFs for efficiency",type=int)
    parser.add_argument('-i', '--id_file', default="/net/fantasia/home/sarahgra/Collaborator_projects/PRS_prediction_Cristen/MGI_sample_IDs")
    parser.add_argument('-o', '--output_file', default="polygenic_risk_score_results.txt")
    args=parser.parse_args()
    print >> sys.stderr, "%s\n"  % args
    return args



###############################
######## FUNCTIONS  #########
###############################

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
def make_weights_dict(weight_file,chrom,pos,coord,ea,weight,vcf_chrom):
    #dictionary format variant:(effectallele, effect)
    weight_dict={}
    with open(weight_file) as f:
        for line in f:
            ls=line.rstrip()
            if ls[0].isdigit(): #assumes we ignore header lines not starting with a digit
                lineList=ls.split() #assumes whitespace delimiter, space or tab
                if chrom is not None and pos is not None: #chr and pos are separate  in weight file
                    if vcf_chrom is not None: #only save info for chromosome of interest
                        if int(lineList[chrom])==vcf_chrom:
                            coordinate=":".join([str(lineList[chrom]),str(lineList[pos])])
                            weight_dict[coordinate]=(lineList[ea],float(lineList[weight]))
                else: #chr and pos are in coordinate form in weight file
                    if vcf_chrom is not None:  #only save info for chromosome of interest
                        if lineList[coord].split(":")[0]==vcf_chrom:
                            weight_dict[lineList[coord]]=(lineList[ea],float(lineList[weight]))
    return weight_dict

#write out regions file to use with tabix, gets coordinates from weights file
def make_regions_file(weight_file,chrom,pos,coord,ea,weight,vcf_chrom,output_name):
    with open(weight_file) as f:
        if chrom is not None and pos is not None: #chr and pos are separate  in weight file
            regions = np.genfromtxt(f, usecols=(chrom,pos), names=("Chr", "Pos"), dtype=None, skip_header=16)
        else:   #coordinate in weight file
             regions = np.genfromtxt(f, delimiter=":", usecols=(coord,coord+1), names=("Chr", "Pos"), dtype=None, skip_header=16) 
    if vcf_chrom is not None:
        regions=regions[regions["Chr"] == int(vcf_chrom)] #filter score file to just the chromosome the VCF corresponds to
    regions = np.sort(regions, order=["Chr", "Pos"])
    region_count=regions.size
    print >> sys.stderr, "Now writing regions file: %s" % output_name
    with open(output_name, 'w') as out:
        np.savetxt(output_name, regions, delimiter="\t", fmt='%d')
    return region_count

#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args=get_settings()

    #create dictionary of weights per variant and write out regions file for tabix
    weight_file_dict=make_weights_dict(args.weight_file,args.chrom_col,args.pos_col,args.coord_col,args.ea_col,args.weight_col,args.chrom)

    #open ID file
    with open(args.id_file) as f:
        sample_id = [line.rstrip() for line in f]

    #Create dictionary to keep track of total scores per person, set initial value to zero
    sample_score_dict = {x:0 for x in sample_id}    

    #Write out regions file for tabix, getting regions from score file
    regions_output_name = "Regions_" + str(args.chrom) + args.weight_file.split("/")[-1]
    region_count=make_regions_file(args.weight_file,args.chrom_col,args.pos_col,args.coord_col,args.ea_col,args.weight_col,args.chrom,regions_output_name)
        
    #Calculate weighted dosages per individual, Make sure to check allele

    ## Handle multi or single VCF
    if args.multi_vcf is not None:
        file_list = glob(args.multi_vcf)
    elif args.single_vcf is not None:
        file_list = [args.single_vcf]

    #record the number of markers actually found and included because risk marker list may differ from variants in data of interest
    variant_count=0
    #loop over VCFs
    for file_x in file_list:
        print >> sys.stderr, "Now reading in: %s" % file_x
        f = Popen(['tabix', '-R', regions_output_name, file_x], stdout=PIPE)
        f = f.communicate()[0]
        
        #If the tabix command finds one of the risk variants in that chunk:
        #To do: record the number of markers actually found and included because risk marker list may differ from variants in data of interest
        if f != '':
            variant_count+=1 #add to marker count
            #This will return the effect allele and effect for each line if variant is in the risk score and not a comment line, and then the dosage line
            test_results = np.asarray([flatten((score_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][0], weight_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][1], line.rstrip().split()[0:5], [value.split(":")[1] for value in line.rstrip().split()[9:]])) for line in f.rstrip().split("\n") if line.split()[0][0] != "#" if (line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4]) in weight_file_dict])
            #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
            #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....

            #if (np.shape(test_results)[0] == 1):
            #    test_results = np.vstack(test_results, np.zeros(np.shape(test_results))
            #Assumes DS in VCF is in terms of the alternate allele 
            #Where effect allele from risk score formula matches alternative allele, multiply directly
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

    print >> sys.stderr, "%d variants were in the region file and %d were ultimately found in the VCF(s)"  % (region_count,variant_count)
    
    with open(args.output_file, 'w') as out:
        out.write("%s\t%s\n" % ("individual", "score"))
        for x in range(len(sample_id)):
            out.write("%s\t%.8f\n" % (sample_id[x], sample_score_dict[sample_id[x]]))

##### Call main 
if __name__ == "__main__":
        main()
