#!/usr/bin/python3

### Written by Dr. Ida Surakka in 2018
### Revised by Brooke Wolford in March 2019


'''
Created on 24.1.2018
Class to read in .dose file and write out dosage file
@author: isurakka
'''
import argparse
from collections import OrderedDict
import os
import sys

__version__= "1.0"

class SNP(object):
    """
    Object to contain data from single .dose line
    Will also keep track of dosages and weights if needed
    """
    def __init__(self, line, allow_indels):
        self.rsid = None
        self.chrom = None
        self.pos = None
        self.a1 = None
        self.a2 = None
        self.dosages = None
        self.weight = None
        self.reversed = False
        self.allow_indels = False
        ok = self.readGenFileLine(line)
        if not ok:
            raise ValueError(line)
        
    def readGenFileLine(self, line):
        """
        Reads single line of genotype
        """
        words = line.split()
        
        # first set the position that holds the first allele, element no 4
        first_allele = 3
            
        if self.allow_indels == False:
            # ignore line if the allele is insertion or deletion
            if len(words[first_allele]) > 1 or len(words[first_allele+1]) > 1:
                return False
        
        self.chrom = words[1]
        # remove leading zero from chrom if necessary
        if self.chrom[0] == "0":
            self.chrom = self.chrom[1:]
        self.rsid = words[0]
        self.pos = words[2]
        self.a1 = words[first_allele]
        self.a2 = words[first_allele + 1]
        # start column of the genotype data
        geno_start = first_allele + 2
        # each snp line contains n genotypes
        # genotypes consists of a single number; dosage for the A2
        # count number of genotypes in the snp line
        no_genotypes = int(len(words) - geno_start)
        self.dosages = [0] * no_genotypes
        # calculates and stores all the dosages for the snp
        for geno_no in range(no_genotypes):
            # calculate dosage for single genotype
            geno_index = geno_no + geno_start
            self.dosages[geno_no] = float(words[geno_index])
        return True
    
    def getRSID(self):
        return self.rsid
    
    def getChrom(self):
        return self.chrom
    
    def getUniquePosID(self):
        """
        Returns position ID which is unique for all positions -
        chrom:pos_a1_a2
        """
        return self.chrom + ":" + self.pos + "_" + self.a1 + "_" + self.a2
    
    def getPosition(self):
        return self.pos
    
    def getAlleles(self):
        return self.a1, self.a2
    
    def getDosages(self):
        return self.dosages
    
    def getDosage(self, i):
        return self.dosages[i]
    
    def setWeight(self, weight, reversed):
        if self.weight != None:
#            print ("Snp %s %s:%s" % (self.rsid, self.chrom, self.pos))
#            print ("Setting weight second time; original weight %s" % self.weight)
#            print ("New weight: %s" % weight)
#            print ("Exiting program")
            sys.exit()
        self.weight = weight
        self.reversed = reversed
    
    def getWeightedDosage(self, i):
        if self.weight == None:
            return None
        if not self.reversed:
#            print ("calculating weight for snp %s" % self.rsid)
#            print ("dosage is %s" % self.dosages[i])
#            print ("weights is %s" % self.weight)
            return self.dosages[i] * self.weight
        return (2.0 - self.dosages[i]) * self.weight
    
    def getUsedAllele(self):
        if self.reversed:
            return self.a1
        return self.a2
    
    def getDosageLine(self):
        line = "%s %s %s %s %s" % (self.chrom, self.rsid, self.pos, self.a1, self.a2)
        for dosage in self.dosages:
            line += " %f" % dosage
        return line

    def __repr__(self):
        s = "%s %s %s %s %s" % (self.rsid, self.chrom, self.pos, self.a1, self.a2)
        return s

class Impute2Dosage(object):
    '''
    classdocs
    '''
    def __init__(self, args = sys.argv[1:]):
        '''
        Constructor gets the command line arguments and feeds them to argparse
        self.snps is an ordered dictionary of SNP objects and will contain all snps from the .dose file
        '''
        self.args = parseArguments(args)
        self.snps = OrderedDict()
        
    def checkInputFiles(self):
        files_to_check = (self.args.input_fn, self.args.sample_fn, self.args.snp_weights_fn)
        for fn in files_to_check:
            if fn != None:
                if not os.path.isfile(fn):
                    print ("%s is not a file!\n" % fn)
                    return False
        return True 
                    
    def readGenFile(self):
        """
        Reads the .dose file. If there are insertions or deletions, those lines 
        are skipped from analysis but they are written into file 
        inserts_and_deletions.txt
        """
        if self.args.sample_fn == None:
            if self.args.verbose:
                print ("Opening file %s for output\n" % self.args.output_fn)
            out_f = open(self.args.output_fn, "w")
        error_f = None
        with open(self.args.input_fn) as gen_f:
            i = 0
            for line in gen_f:
                ignore_line = False
                try:
                    snp = SNP(line, self.args.allow_indels)
                except ValueError:
                    # the snp was insertion or deletion - write the line to 
                    # file inserts_and_deletions.txt
                    ignore_line = True
 #                   if error_f == None:
 #                       error_f = open("inserts_and_deletions.txt", "w")
 #                   error_f.write("%s" % line)
                if not ignore_line:
                    if self.args.sample_fn == None:
                        # no sample file, just print out the dosage file
                        out_f.write("%s\n" % snp.getDosageLine())
                    else:
                        # there is a sample file, so the SNP information has to be stored
                        snp_id = snp.getUniquePosID()
                        self.snps[snp_id] = snp
                if i % 1000 == 0:
                    print ("\rReading row number %d in file %s\n" % (i, self.args.input_fn), end="")
                    sys.stdout.flush()
                i += 1
        if error_f != None:
            error_f.close()
        if self.args.verbose:
            print ("File %s read\n" % self.args.input_fn)
        if self.args.sample_fn == None:
            out_f.close()
            print ("File %s written\n" % self.args.output_fn)
        return
            
    def readSampleFile(self):
        """
        Reads the sample file. Skips first header row, change that if they are used.
        """
        samples = []
        i = 0
        with open(self.args.sample_fn) as sample_file:
            for line in sample_file:
                # skip first header row
                if i not in [0]:
                    # data row, check the number of elements
                    words = line.split()
                    if len(words) < 2:
                        # too short data row, give error and stop
                        print("Sample file %s has row with too few elements\n" % self.args.sample_fn)
                        print("Every row has to have FID & IID\n")
                        print("%s\n" % line)
                        sys.exit()
                    else:
                        # ok data row - store only 2 first elements          
                        samples.append(words[:2])
                i += 1
        return samples

    def readSNPWeights(self):
        """
        Reads the weights file, finds matching SNP and sets the weight for it. Assumes header, either text or #.
        """

        with open(self.args.snp_weights_fn) as weights_f:
            for line in weights_f:
                if line[0].isdigit(): #expects first value will be a number, not text or #
                    found_match = False
                    words = line.split()
                    chrom = words[3]
                    # remove leading zero from chrom if necessary
                    if chrom[0] == "0":
                        chrom = chrom[1:]
                    if chrom == self.args.chromosome_no or self.args.chromosome_no == None:
                        rsid = words[0]
                        allele = words[1]
                        weight = float(words[2])
                        pos = words[4]
                        #print ("Reading rsid %s: %s" % (rsid, allele))
                        for snp_id in self.snps.keys():
                            # the snp is identified by chromosome, position and allele
                            if (chrom == self.snps[snp_id].getChrom() and
                                pos == self.snps[snp_id].getPosition()):
                                # chromosome and position match, check the allele
                                found_match = self.checkAlleleMatch(allele.upper(), weight, snp_id)
                                if found_match:
                                    break
                        if not found_match and self.args.verbose:
                            print ("No match found for entry: %s\n" % line)
     
    def checkAlleleMatch(self, allele, weight, snp_id):
        """
        Check that the alleles match, or in case of alleles being I and D,
        check that entry is indel
        """
        found_match = False
        a1,a2 = self.snps[snp_id].getAlleles()
        reversed = False
        if allele.upper() not in ["I", "D"]:
            # just allele strings, check for match
            if allele == a1 or allele == a2:
                # got match
                if allele == a1:
                    reversed = True                                      
                found_match = True
        else:
            # check for indel
            if len(a1) > 1 or len(a2) > 1:
                # match found, determine reversed
                # indel, set longer to I (insert) and shorter to D (deletion)
                if allele.upper() == "I":
                    # insertion
                    if len(a1) > len(a2):
                        reversed = True
                else:
                    # deletion
                    if len(a1) < len(a2):
                        reversed = True
                found_match = True
        if found_match:
            self.snps[snp_id].setWeight(weight, reversed)
        return found_match
                            
    def printSampleFile(self, samples):
        """
        Prints sample file with dosages as entries
        """
        out_f = open(self.args.output_fn, "w")
        header = "FID IID"
        for snp_id in self.snps.keys():
            used_allele = self.snps[snp_id].getUsedAllele()
            header += " %s_%s" % (snp_id, used_allele)
        out_f.write("%s\n" % header)
        i = 0
        for sample in samples:
            row = "%s %s" % (sample[0], sample[1])
            for snp_id in self.snps.keys():
                # print the asked dosage from the array
                row += " %s" % self.snps[snp_id].getDosage(i)
            out_f.write("%s\n" % row)
            if i % 100 == 0:
                print ("Writing row number %d to file %s\n" % (i, self.args.output_fn))
                sys.stdout.flush()
            i += 1
        print ("Printed sample file %s\n" % self.args.output_fn)
        
    def printWeightedSampleFile(self, samples):
        """
        Prints sample file with weighted dosages as entries
        If weighted dosage cannot be calculated for some snps, their id's will be printed to a file
        failed_snp_ids.txt
        """
        out_f = open(self.args.output_fn, "w")
        header = "FID IID"
        for snp_id in self.snps.keys():
            used_allele = self.snps[snp_id].getUsedAllele()
#            header += " %s_%s" % (snp_id, used_allele)
        out_f.write("%s Summed_weight\n" % header)
        i = 0
        failed_snp_ids = []
        for sample in samples:
            row = "%s %s" % (sample[0], sample[1])
            weight_sum = 0
            for snp_id in self.snps.keys():
                # print the asked dosage from the array
                weighted_dosage = self.snps[snp_id].getWeightedDosage(i)
                if weighted_dosage == None:
                    if snp_id not in failed_snp_ids:
                        failed_snp_ids.append(snp_id)
#                    row += " NA"
                else:
#                    row += " %s" % weighted_dosage
                    weight_sum += weighted_dosage
            out_f.write("%s %f\n" % (row, weight_sum))
            if i % 10000 == 0:
                print ("Writing row number %d to file %s\n" % (i, self.args.output_fn))
                sys.stdout.flush()
            i += 1
#        if len(failed_snp_ids) > 0:
#            # if some snip id:s failed, print them into separate file
#            out_f = open("failed_snp_ids.txt", "w")
#            out_f.write("Failed snp ids:\n")
#            for snp_id in failed_snp_ids:
#                out_f.write("%s\n" % snp_id)
#            out_f.close()
        print ("Printed weighted sample file %s\n" % self.args.output_fn)
        
    def main(self):
#         if self.args.snp_weights_fn != None and self.args.chromosome_no == None:
#             print ("Chromosome number (--chromosome_no) is obligatory when weighted dosages are calculated")
#             return
        if self.args.snp_weights_fn != None and self.args.sample_fn == None:
            print("Cannot calculate SNP weights without sample file (--sample_fn)\n")
            return
        ok = self.checkInputFiles()
        if not ok:
            return
        self.readGenFile()
        if self.args.sample_fn != None:
            samples = self.readSampleFile()
            if self.args.snp_weights_fn == None:
                self.printSampleFile(samples)
            else:
                self.readSNPWeights()
                self.printWeightedSampleFile(samples)
                
        
def parseArguments(args):
    """
    Handles the argument parsing
    """
    parser = argparse.ArgumentParser(prog = "impute2dosage",
                                     description = "Converts .dose file into a weighted dosages file if both sample file and SNP weights (and chromosome number) are provided.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_fn", help = "The .dose file name")
    parser.add_argument("--chromosome_no", help = "The number of the chromosome being analysed", default = None)
    parser.add_argument("--output_fn", help = "Output file name, default is dosages.gen", default = "dosages.gen")
    parser.add_argument("--sample_fn", help = "Sample file name", default = None)
    parser.add_argument("--snp_weights_fn", help = "File with SNP weights. Also used for checking the snps by chromosome, position and allele. If this file is given, output file contains weighted dosages",
                        default = None)
    #parser.add_argument("--use_allele1", help = "If set, the first allele is used in calculations instead of the second",
    #                    action = "store_true", default = False)
    parser.add_argument("--allow_indels", help ="Allow matching indels in the gen file",
                        default = False, action = "store_true")
    parser.add_argument("--verbose", help="If true, prints possible error messages during calculations. Default is False.",
                        action = "store_true", default = False)
    parser.add_argument("--version", help = "Shows the version number of the program",
                        action = 'version', version = "%(prog)s {version}".format(version=__version__))
    return parser.parse_args(args)
        
if __name__ == '__main__':
    app = Impute2Dosage()
    app.main()
    
