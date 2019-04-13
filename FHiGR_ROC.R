

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
library(ggplot2)

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
  make_option("--maintitle", type="character", default="ROC",help="Plot title [default='']"),
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
xlabel<-arguments$options$xlabel
ylabel<-arguments$options$ylabel
dig<-arguments$option$digits
invNorm<-arguments$options$invNorm
header<-arguments$options$header

###########################################################
#################### FUNTIONS #############################
###########################################################

## identify top and bottom of distribution (e.g. prediction) for ROC from all data
make_pred_obj<-function(df,value,grs_col,pheno_col){
    q<-quantile(df[[grs_col]],value)
    labels<-df[[pheno_col]] #phenotype labels
    pred<-ifelse(df[[grs_col]]<=q,0,1) #top of distribution or bottom, considered as prediction
    return(list(pred,labels))
}

## identify top and bottom of distribution (e.g. prediction) for ROC using top and FH=1 versus all other distribution
make_pred_obj_strat<-function(df,value,grs_col,pheno_col,strat_col,qfirst=FALSE){
  labels<-df[[pheno_col]] #phenotype labels
  if (qfirst==TRUE){
    q<-quantile(df[[grs_col]],value) #quantile in all data
    #top of bottom of distribution as binary variable, bottom as reference group
    df$dist<-2
    df$dist<-ifelse(df[[grs_col]]<=q,0,1) #1 if in top, 0 in bottom
    df[df$dist==2]$dist<-NA
  } else {
    df$dist<-2
    for (s in c(1,2)){
      q<-quantile(df[df[[strat_col]]==(s-1)][[grs_col]],value,na.rm=TRUE)
      df[df[[strat_col]]==(s-1),'dist']<-ifelse(df[df[[strat_col]]==(s-1)][[grs_col]]<=q,0,1) #1 if in top, 0 in bottom
    } 
    df[df$dist==2,'dist']<-NA #if stratum was missing 1/0
  }
  #collapse FH=0 and bottom distirbution FH=1 into 1 
  df$FHIGR<-2
  df[df[[strat_col]]==1 & df$dist==1,'FHIGR']<-1 #positive fam hx, top dist
  df[df[[strat_col]]==1 & df$dist==0,'FHIGR']<-0 #positive fam hx, bottom dist
  df[df[[strat_col]]==0,'FHIGR']<-0 #negative fam hx
  df[df$FHIGR==2,'FHIGR']<-NA
  df$dist<-df$FHIGR #replace I(dist) with I(dist & stratum=1) 
  dist_col<-which(names(df)=="dist")
  return(list(df[[dist_col]],labels))
}
  

#return ROC plot xy pair which are false positive rate and true positive rate
ROC_pair<-function(obj){
    m<-as.matrix(table(obj))
    sensitivity<-m[2,2]/sum(m[,2]) #true positive rate
    specificity<-m[1,1]/sum(m[,1]) 
    tpr<-sensitivity
    fpr<-1-specificity
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

## make cut in entire distribution
cutpts<-(1:99)/100
pred_obj<-lapply(cutpts,make_pred_obj,df=subset,pheno_col=pheno_col,grs_col=grs_col)
pairs<-lapply(pred_obj,ROC_pair)
roc_df<-bind_rows(pairs)
roc_df$cutpts<-cutpts
roc_df$method<-"standard"

all<-roc_df

## make cut in F=1 
logical_list<-c("TRUE","FALSE")
label_list<-c("quantileFirst","stratifyFirst")
for (l in c(1,2)){ #do for each division logic
  roc_df<-NULL
  pred_obj<-lapply(cutpts, make_pred_obj_strat,df=subset,pheno_col,grs_col=grs_col,strat_col=strat_col,qfirst=logical_list[l])
  pairs<-lapply(pred_obj,ROC_pair)
  roc_df<-rbind(roc_df,cbind(bind_rows(pairs),cutpts,method="FHiGRS"))
  roc_df<-rbind(all,roc_df)
  
  pdf_fn<-paste(sep=".",out,label_list[l],"ROC.pdf")
  pdf(file=pdf_fn,height=4,width=6)
  print(ggplot(roc_df,aes(x=fpr,y=tpr,color=method)) + theme_bw() + geom_point() +
    coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + scale_color_manual(values=c("grey","darkblue")) + 
    labs(title=main,xlab="False Positive Rate",ylab="True Positive Rate"))
  dev.off()
  
}
                  




