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

from __future__ import division
import argparse
from subprocess import Popen, PIPE
import gzip
import io
import numpy as np
from glob import glob

#This python script can be used to calculate genetic risk scores from MGI dosages in VCF file
#"EUR_p_lt_5e-8_r2_0.2_plink.clumped_scores_per_variant.txt"
parser = argparse.ArgumentParser(description='Enter files to use for PRS calculation')
parser.add_argument('-s', '--score_file', default="scores_per_variant.txt")
parser.add_argument('-t', '--test_dosages', nargs='*', default="/net/hunt/disk2/MGI_data_2/MGI_all_data_HPI-3498/Willer_HPI3498_HUM00094409_mgi_HRC.chr*.vcf.gz")
parser.add_argument('-i', '--id_file', default="/net/fantasia/home/sarahgra/Collaborator_projects/PRS_prediction_Cristen/MGI_sample_IDs")
parser.add_argument('-o', '--output_file', default="polygenic_risk_score_results.txt")
args=parser.parse_args()

#Create dictionary of scores per variant
#dictionary format variant:(effectallele, effect)
with open(args.score_file) as f:
    score_file_dict = {line.split()[0]:(line.split()[1],float(line.split()[2])) for line in f}

with open(args.id_file) as f:
    sample_id = [line.rstrip() for line in f]

#Create dictionary to keep track of total scores per person, set initial value to zero
sample_score_dict = {x:0 for x in sample_id}    

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

#Write out regions file for tabix 
with open(args.score_file) as f:
    regions = np.genfromtxt(f, delimiter=":", usecols=(0,1), names=("Chr", "Pos"), dtype=None)

regions = np.sort(regions, order=["Chr", "Pos"])

regions_output_name = "Regions_" + args.score_file
with open(regions_output_name, 'w') as out:
    np.savetxt(regions_output_name, regions, delimiter="\t", fmt='%d')
        


#Calculate dosages per individual
#Make sure to check allele
file_list = glob(args.test_dosages)

for file_x in file_list:
    print("Now reading in: %s" % file_x)
    f = Popen(['tabix', '-R', regions_output_name, file_x], stdout=PIPE)
    f = f.communicate()[0]

    #If the tabix command finds one of the risk variants in that chunk:
    if f != '':
        #This will return the effect allele and effect for each line if variant is in the risk score and not a comment line, and then the dosage line
        test_results = np.asarray([flatten((score_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][0], score_file_dict[(line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4])][1], line.rstrip().split()[0:5], [value.split(":")[1] for value in line.rstrip().split()[9:]])) for line in f.rstrip().split("\n") if line.split()[0][0] != "#" if (line.split()[0] + ":" + line.split()[1] + ":" + line.split()[3] + ":" + line.split()[4]) in score_file_dict])
        #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
        #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....

        #if (np.shape(test_results)[0] == 1):
        #    test_results = np.vstack(test_results, np.zeros(np.shape(test_results))
        #Where effect allele from risk score formula matches alternative allele, multiply directly
        matching_effect_allele = test_results[np.where(test_results[:,0] == test_results[:,6])]
        #[:, np.newaxis] this is needed to do the multipication element wise (first column * all dosages in row)
        matching_effect_allele = matching_effect_allele[:,1].astype(float)[:, np.newaxis] * matching_effect_allele[:,7:].astype(float)
            
        #Where effect allele matches reference allele, take 2-dosage, then multiply (so flip dosage to be for alternative allele)
        matching_reference_allele = test_results[np.where(test_results[:,0] == test_results[:,5])]
        matching_reference_allele = matching_reference_allele[:,1].astype(float)[:, np.newaxis] * (2 - matching_reference_allele[:,7:].astype(float))
            
        #Sum down columns
        dosage_scores_sum = np.sum(matching_reference_allele, axis=0) + np.sum(matching_effect_allele, axis=0)
        for x in range(len(dosage_scores_sum)):
            sample_score_dict[sample_id[x]] = sample_score_dict[sample_id[x]] + dosage_scores_sum[x]

with open(args.output_file, 'w') as out:
    out.write("%s\t%s\n" % ("individual", "score"))
    for x in range(len(sample_id)):
        out.write("%s\t%.8f\n" % (sample_id[x], sample_score_dict[sample_id[x]]))
