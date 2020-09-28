#!/usr/bin/env python3

#===============================================================================
# Copyright (c) 2020 Brooke Wolford
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
#============================================================================

############################
##### IMPORT MODULES #######
###########################
import os
import sys
import sqlite3
import rsidx
import argparse
import gzip

###########################
##### PARSE ARGUMENTS ####
###########################
def get_settings():
    parser = argparse.ArgumentParser(description="Convert weights files with just rsIDs to ones with chromosomes and coordinates too\n")
    parser.add_argument("-v","--vcf",help="bgzipped and tabixed vcf file",type=str,required=True)
    parser.add_argument("-r","--rsidx",help="rsidx file",type=str)
    parser.add_argument("-w","--weights",help="weights file. assumes one header. will skip lines with #",type=str,required=True)
    parser.add_argument("-c","--col",help="0-based column in weights file that rsID is in",type=int,default=0)
    parser.add_argument("-p","--prefix",help="Prefix for output file, including path",type=str,required=True)
    args=parser.parse_args()
    return args

###############################
######## SUB ROUTINES  #########
###############################

def open_zip(f):
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print("Opening gzipped file %s\n" % f,file=sys.stderr)
    elif f == "-":
        command=sys.stdin()
    else:
        command=open(f,"rt")
        print("Opening file %s\n" % f,file=sys.stderr)
    return command

def index(vcf):
    with sqlite3.connect('myidx.db') as dbconn, open(vcf, 'r') as vcffh:
        rsidx.index.index(dbconn, vcffh)

#            rsidx index 00-All.vcf.gz 00-All.vcf.rsidx
    
def search(rsidlist,vcf,index):
    print("Querying markers from weights file in VCF\n",file=sys.stderr)
    in_len=len(rsidlist)
    rsid_dict={}
    with sqlite3.connect(index) as dbconn:
        for line in rsidx.search.search(rsidlist, dbconn, vcf):
            ls=line.rstrip()
            lineList=ls.split("\t")
            rsid_dict[lineList[2]]=lineList[:5]  #assumes VCF is chr, pos, rsID, REF, ALT
    out_len=len(rsid_dict.keys())
    if in_len!=out_len:
        diff=int(in_len)-int(out_len)
        print("Not all rsIDs from weights file could be found in the VCF. Missing %d of %d\n" % (diff,in_len),file=sys.stderr)
    else:
        print("All %d rsIDs from weights file could be found in the VCF.\n" % in_len,file=sys.stderr)
    return rsid_dict
               

def rsid_from_weights(weights,col):
    print("Getting rsIDs from weights file %s\n" %weights, file=sys.stderr)
    command=open_zip(weights)
    rsid_list=[]
    header_count=0
    with command as f:
        for line in f:
            ls=line.rstrip()
            if ls[0] != "#":
                if header_count==0: #skip first header after any lines with #
                    header_count+=1
                    next
                else:
                    lineList=ls.split()
                    rsid_list.append(lineList[col])
    return rsid_list


def merge(weights,col,rsid_dict,prefix):
    command=open_zip(weights)
    header_count=0
    output=prefix + "_reformat.txt"
    print("Writing new weights file %s\n" %output, file=sys.stderr)
    with open(output,"w") as o:
        with command as f:
            for line in f:
                ls=line.rstrip()
                if ls[0] == "#":
                    o.write(ls+"\n")
                elif header_count==0:
                    lineList=ls.split()
                    o.write("\t".join(lineList+["CHR","POS","REF","ALT"])+"\n") #write header
                    header_count+=1
                else:
                    lineList=ls.split()
                    try:
                        from_vcf=rsid_dict[lineList[col]] #look up from dictionary from vcf
                        ## handle occurence of multiple alt alleles by printing each potential entry as a newline
                        for alt_allele in from_vcf[4].split(","):
                            o.write("\t".join(lineList+from_vcf[0:2]+[alt_allele])+"\n")
                    except KeyError:
                        o.write("\t".join(lineList+["NA"]*4)+"\n") #rsid not in VCF
        f.close()
    o.close()
    os.system("gzip "+output) #ystem call to gzip
    return 0

#########################
########## MAIN #########
#########################

def main():

    #get arguments
    args = get_settings()
    print(args)
    
    #github package https://github.com/bioforensics/rsidx

    if not os.path.exists(args.vcf):
        sys.exit("VCF does not exist\n")
    if not os.path.exists(args.rsidx):
        sys.exit("RSIDX does not exist\n")
    #index(args.vcf)

    #get rsids from weights file
    rsid_list=rsid_from_weights(args.weights,args.col)
    
    #search vcf
    rsid_dict=search(rsid_list,args.vcf,args.rsidx)

    #merge new info with weights file
    merge(args.weights,args.col,rsid_dict,args.prefix)

#call main
if __name__ == "__main__":
      main()
111
