

#https://github.com/RainCloudPlots/RainCloudPlots/tree/master/tutorial_R

#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################

print(Sys.time())
print(sessionInfo())
library(optparse)
library(data.table)
library(ROCR)
library(dplyr)

###########################################################
################### Read Command Line Parameters ##########
###########################################################

optionList <- list(
  make_option(c("-f","--file"),type="character",help="File with sample IDs, phenotypes, self reported family history, and GRS. Expeects header. White space delimited."),
  make_option(c("-s","--stratum_col"),type="numeric",help="1-based column with stratum (e.g. family history). Must be on binary scale 1/0 with 1 being affirmative. NAs ok."),
  make_option(c("-p","--pheno_col"),type="numeric",help="1-based column with phenotype (e.g. disease status). Must be on binary scale 1/0 with 1 being case. NAs ok."),
  make_option(c("-g","--grs_col"),type="numeric",help="1-based column with GRS, not inverse normalized."),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-i","--invNorm"),type="logical",default=FALSE,help="Inverse normalize GRS for entire population [default=FALSE]"),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--xlabel",type="character",default="GRS",help="X-axis label [default='']"),
  make_option("--ylabel",type="character",default="Prevalence",help="Y-axis label [default='']")
)

parser <- OptionParser(
  usage="%prog --file --stratum_col --pheno_col --grs_col",
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
  warning("Stratum col is required -sc or --stratum_col")
}
pheno_col<-arguments$options$pheno_col
if (length(pheno_col)==0){
  warning("Pheno col is required -pc or --pheno_col")
}
grs_col<-arguments$options$grs_col
if (length(grs_col)==0){
  warning("GRS col is required -gc or --grs_col")
}
out<-arguments$options$output
main<-arguments$options$maintitle
xlab<-arguments$options$xlabel
ylab<-arguments$options$ylabel
dig<-arguments$option$digits
invNorm<-arguments$options$invNorm
header<-arguments$options$header

###########################################################
#################### FUNTIONS #############################
###########################################################

## all data
make_pred_obj<-function(df,value){
    q<-quantile(df[[grs_col]],value)
    labels<-df[[pheno_col]] #phenotype labels
    pred<-ifelse(df[[grs_col]]<=q,0,1)
    return(list(labels,pred))
}

#return ROC plot xy pair which are false positive rate and true positive rate
ROC_pair<-function(obj){
    m<-as.matrix(table(obj))
    neg_predictive<-m[1,1]/sum(m[1,])
    pos_predictive<-m[2,2]/sum(m[2,])
    specificity<-m[1,1]/sum(m[,1])
    sensitivity<-m[2,2]/sum(m[,2]) #true positive rate
    tpr<-sensitivity
    fpr<-m[2,1]/sum(m[2,])
    accuracy<-(m[1,1] + m[2,2])/sum(m)

    ##data.frame(false__pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy)
    return(data.frame(fpr,tpr))
}

###########################################################
########################## MAIN ###########################
###########################################################

##read data 
df<-fread(file,header=header)
#print(dim(df))

subset<-df[!is.na(df[[strat_col]])] #remove if NA for stratum
print(dim(subset))
#ggplot(subset, aes_string(x=names(subset)[grs_col])) + geom_density()

if (invNorm==TRUE){
  ## inverse rank normalize GRS in thepopulation
  subset$invNormGRS<-rankNorm(subset[[grs_col]])
  grs_col<-which(names(subset)=="invNormGRS") #new GRS col
}

cutpts<-(1:99)/100
pred_obj<-lapply(cutpts,make_pred_obj,df=subset)
pairs<-lapply(pred_obj,ROC_pair)
roc_df<-bind_rows(pairs)
roc_df$cutpts<-cutpts

pdf_fn<-paste(sep=".",out,"ROC.pdf")
pdf(file=pdf_fn,height=6,width=6)
ggplot(roc_df,aes(x=fpr,y=tpr)) + theme_bw() + geom_point()
dev.off()
