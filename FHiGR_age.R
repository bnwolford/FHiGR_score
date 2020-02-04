#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################
library(dplyr)
library(ggplot2)
library(data.table)
#library(ggridges)
library(gridExtra)
library(RNOmni)
library(optparse)

print(Sys.time())
print(sessionInfo())

###########################################################
################### Read Command Line Parameters ##########
###########################################################

optionList <- list(
  make_option(c("-f","--file"),type="character",help="File with sample IDs, phenotypes, self reported family history, and GRS. Expeects header. White space delimited."),
  make_option(c("-s","--stratum_col"),type="numeric",help="1-based column with stratum (e.g. family history). Must be on binary scale 1/0 with 1 being affirmative. NAs ok."),
  make_option(c("-p","--pheno_col"),type="numeric",help="1-based column with phenotype (e.g. disease status). Must be on binary scale 1/0 with 1 being case. NAs ok."),
  make_option(c("-g","--grs_col"),type="numeric",help="1-based column with GRS, not inverse normalized."),
  make_option(c("-a","--age_col"),type="numeric",help="1-based column with participation or enrollment age e.g. age of self reported family history"),
  make_option(c("-b","--birthYear_col"),type="numeric",help="1-based column with birthyear"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--codeDir",type="character",default="/FHiGR_score/",help="Directory for repository for sourcing other code in code base [default=/FHiGR_score/]")
)

parser <- OptionParser(
  usage="%prog --file --stratum_col --pheno_col --grs_col --age_col --output --digits --header --maintitle --codeDir",
  option_list=optionList
)
arguments <- parse_args(parser, positional_arguments=TRUE)
print(arguments$options)

##print warnings when they happen
options(warn=1)

#check arguments without defaults
file<-arguments$options$file
if (length(file)==0){
  warning("File is required -f or --file")
}
strat_col<-arguments$options$stratum_col
if (length(strat_col)==0){
  warning("Stratum col is required -s or --stratum_col")
}
pheno_col<-arguments$options$pheno_col
if (length(pheno_col)==0){
  warning("Pheno col is required -p or --pheno_col")
}
grs_col<-arguments$options$grs_col
if (length(grs_col)==0){
  warning("GRS col is required -g or --grs_col")
}
age_col<-arguments$options$age_col
if (length(age_col)==0){
    warning("Age col is required -a or --age_col")
}
birthYear_col<-arguments$options$birthYear_col
if (length(birthYear_col)==0){
    warning("Birthyear col is required -b or --birthYear_col")
}
out<-arguments$options$output
main<-arguments$options$maintitle
dig<-arguments$option$digits
header<-arguments$options$header


###########################################################
#################### FUNCTIONS #############################
###########################################################




###########################################################
#################### MAIN #################################
###########################################################

##read data 
dat<-fread(file,header=header)
print(paste("Data dimensions are:",dim(dat)[1],dim(dat)[2]))

print(table(dat[[strat_col]],useNA="always"))

## age at time of self report for everyone                                      
pdf_fn<-paste(sep=".",out,"age_at_SR.pdf")
pdf(file=pdf_fn,height=6,width=6,useDingbats=FALSE)
ggplot(dat,aes(x=get(names(dat)[age_col]))) + geom_density(fill="black",alpha=0.75) + theme_bw() +
    labs(x="Enrollment age") +
    theme(title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()


##age at time of self report and family history for 1/0
pdf_fn<-paste(sep=".",out,"age_at_SR_stratify.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[!is.na(dat[[strat_col]])],aes(x=get(names(dat)[age_col]),fill=factor(get(names(dat)[strat_col])))) + geom_density(alpha=0.7) + theme_bw() +
    scale_fill_manual(values=c("goldenrod","dark blue"),name="Family History",labels=c("Negative","Positive")) +
        labs(x="Enrollment Age") +  theme(legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

## birth year distribution stratified
pdf_fn<-paste(sep=".",out,"birthYear_stratify.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[!is.na(dat[[strat_col]])],aes(x=get(names(dat)[birthYear_col]),fill=factor(get(names(dat)[strat_col])))) + geom_density(alpha=0.7) + theme_bw() + scale_fill_manual(values=c("goldenrod","dark blue"),name="Family History",labels=c("Negative","Positive")) +  labs(x="BirthYear") +  theme(legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

##age at time of self report and family history for NA (age may also be missing if FH is missing)
pdf_fn<-paste(sep=".",out,"age_at_SR_stratifyNA.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[is.na(dat[[strat_col]])],aes(x=get(names(dat)[age_col]))) + geom_density(alpha=0.7,fill="grey") + theme_bw() +
       labs(x="Enrollment Age") +  theme(legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

## birth yaer distirbution for samples with NA for stratum
pdf_fn<-paste(sep=".",out,"birthYear_stratifyNA.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[is.na(dat[[strat_col]])],aes(x=get(names(dat)[birthYear_col]))) + geom_density(alpha=0.7,fill="grey") + theme_bw() +
    labs(x="Birth Year") +  theme(legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()


### to do histogram?



