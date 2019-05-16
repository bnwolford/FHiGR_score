

#https://github.com/RainCloudPlots/RainCloudPlots/tree/master/tutorial_R

#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################

print(Sys.time())
print(sessionInfo())
library(optparse)
library(data.table)
library(ROCR) #ROC package
library(dplyr)
library(ggplot2)
library(RNOmni)

###########################################################
################### Read Command Line Parameters ##########
###########################################################

optionList <- list(
  make_option(c("-f","--file"),type="character",help="File with sample IDs, phenotypes, self reported family history, and GRS. Expeects header. White space delimited."),
  make_option(c("-s","--stratum_col"),type="numeric",help="1-based column with stratum (e.g. family history). Must be on binary scale 1/0 with 1 being affirmative. NAs ok."),
  make_option(c("-p","--pheno_col"),type="numeric",help="1-based column with phenotype (e.g. disease status). Must be on binary scale 1/0 with 1 being case. NAs ok."),
  make_option(c("-q","--quantile"),type="numeric",help="Number of quantiles for calculating FHiGRS [default=20]",default=20),
  make_option(c("-g","--grs_col"),type="numeric",help="1-based column with GRS, not inverse normalized."),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-v","--covariates"),type="character",default="",help="Comma separated list of 1-based column incides for model covariates"),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="ROC",help="Plot title [default='']"),
  make_option("--codeDir",type="character",default="/FHiGRS_score/",help="Directory for repository for sourcing other code in code base [default=/FHiGRS_score/]")
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
dig<-arguments$option$digits
invNorm<-arguments$options$invNorm
covar<-as.numeric(strsplit(arguments$options$covariates,",")[[1]])
header<-arguments$options$header
quantile<-arguments$options$quantile
##source relevant code from code base
source(paste0(arguments$options$codeDir,"helperFunctions.R")) ##will be used to calculate FHiGRS


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
    ##data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy)
    return(data.frame(fpr,tpr))
}

## make ROC curve and calculate AUC given a dataframe with predictions and labels
ROCR_package<-function(d){
    p<-prediction(d$pred,d$label) #ROCR prediction object
    roc.perf=performance(p,measure="tpr",x.measure="fpr")
    acc.perf=performance(p,measure="acc")
    #plot(acc.perf)
    auc.perf=performance(p,measure="auc")
    data<-data.frame(x=unlist(roc.perf@x.values),y=unlist(roc.perf@y.values),auc=unlist(auc.perf@y.values))
    return(data)
}

###########################################################
########################## MAIN ###########################
###########################################################

##read data 
df<-fread(file,header=header)
#print(dim(df))

subset<-df[!is.na(df[[strat_col]])] #remove if NA for stratum
print(dim(subset))

## make cut in entire distribution
cutpts<-(1:99)/100
pred_obj<-lapply(cutpts,make_pred_obj,df=subset,pheno_col=pheno_col,grs_col=grs_col)
pairs<-lapply(pred_obj,ROC_pair)
roc_df<-bind_rows(pairs)
roc_df$cutpts<-cutpts
roc_df$method<-"GRS"
all<-roc_df

##estimate FHiGRS (functions from helperFunctions.R)
sobj<-prev_per_quantile_stratum(qtile=quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=FALSE)
for (j in c(1,2)){ #across stratum
    list_length<-length(sobj$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-sobj$prev[j,][-list_length]
    se<-sobj$se[j,][-list_length]
    n<-sobj$n[j,][-list_length]
    tiles<-sobj$tiles[-1]
    lower_tile<-sobj$tiles[1:list_length-1]
    upper_tile<-sobj$tiles[2:list_length]
    percents<-names(tiles)
    bins<-rep(quantile,list_length-1)
    strat<-rep(j-1,list_length-1)
    if (j==1) {
        sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
    } else {
        sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
    }
}
qsub<-estimate_FHiGRS(sdf,subset,strat_col,grs_col)
fhigrs_col<-which(names(qsub)=="FHIGRS")

## inverse normalize GRS
qsub$invNormGRS<-rankNorm(qsub[[grs_col]])
grs_col<-which(names(qsub)=="invNormGRS")

## make cut in FHiGRS distribution
fpred_obj<-lapply(cutpts,make_pred_obj,df=qsub,pheno_col=pheno_col,grs_col=fhigrs_col)
pairs<-lapply(fpred_obj,ROC_pair)
roc_df<-bind_rows(pairs)
roc_df$cutpts<-cutpts
roc_df$method<-"FHiGRS"
fhigrs<-roc_df

## make cut in F=1 stratum
logical_list<-c("TRUE","FALSE")
label_list<-c("quantileFirst","stratifyFirst")
for (l in c(1,2)){ #do for each division logic
  roc_df<-NULL
  pred_obj<-lapply(cutpts, make_pred_obj_strat,df=qsub,pheno_col,grs_col=grs_col,strat_col=strat_col,qfirst=logical_list[l])
  pairs<-lapply(pred_obj,ROC_pair)
  roc_df<-rbind(roc_df,cbind(bind_rows(pairs),cutpts,method="GRS|FH=1"))
  roc_df<-rbind(all,roc_df) #ROC for all GRS
  roc_df<-rbind(fhigrs,roc_df) #ROC for FHiGRS
  #plot
  pdf_fn<-paste(sep=".",out,label_list[l],"ROC.pdf")
  pdf(file=pdf_fn,height=4,width=5,useDingbats=FALSE)
  print(ggplot(roc_df,aes(x=fpr,y=tpr,color=method)) + theme_bw() + geom_point(alpha=0.8) +
        coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
        scale_color_manual(values=c("seagreen4","darkblue","grey")) +
        labs(title=main,x="False Positive Rate",y="True Positive Rate") +
        geom_abline(slope=1,intercept=0,linetype="dashed",color="black"))
  dev.off()
  
}
                  

##### ROCR curve and AUC for GRS, FHIGRS, GRS+FH without covariates

##additive model predicted values 
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(strat_col,grs_col)],collapse="*"),
                          sep = ""))
glm.obj<-glm(formula=formula,data=qsub,family="binomial")
add_fitted<-fitted(glm.obj)

##use functions from ROCR package using the user defined function ROCR_package
fhigrs_df<-data.frame(pred=qsub[[fhigrs_col]],label=qsub[[pheno_col]])
grs_df<-data.frame(pred=qsub[[grs_col]],label=qsub[[pheno_col]])
add_df<-data.frame(pred=add_fitted,label=qsub[[pheno_col]])

grs_roc<-ROCR_package(grs_df)
fhigrs_roc<-ROCR_package(fhigrs_df)
add_roc<-ROCR_package(add_df)

##label models
grs_roc$method<-"GRS"
fhigrs_roc$method<-"FHiGRS"
add_roc$method<-"GRS+FH"
roc_df<-rbind(grs_roc,fhigrs_roc,add_roc)

##pull out AUC
grs_auc<-unique(roc_df[roc_df$method=="GRS",]$auc)
fhigrs_auc<-unique(roc_df[roc_df$method=="FHiGRS",]$auc)
add_auc<-unique(roc_df[roc_df$method=="GRS+FH",]$auc)

#plot
pdf_fn<-paste(sep=".",out,"ROC.pdf")
pdf(file=pdf_fn,height=3,width=4,useDingbats=FALSE)
print(ggplot(roc_df,aes(x=x,y=y,color=method)) + theme_bw() +geom_line() +
      coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
      scale_color_manual(values=c("darkblue","grey","seagreen4"),name="Score") +
      labs(title=main,x="False Positive Rate",y="True Positive Rate") +
      annotate("text",x=0.8,y=0, label=paste0("GRS AUC  ",format(grs_auc,digits=dig,format="f")),color="darkgrey",size=2) + 
      annotate("text",x=0.8,y=0.05,label=paste0("FHiGRS AUC  ",format(fhigrs_auc,digits=dig,format="f")),color="darkblue",size=2) +
      annotate("text",x=0.8,y=0.1,label=paste0("GRS + FH AUC ",format(add_auc,digits=dig,format="f")),color="seagreen4",size=2) +
      geom_abline(slope=1,intercept=0,linetype="dashed",color="black"))
dev.off()


##### ROCR curve and AUC for GRS, FHIGRS, GRS+FH with covariates

## model for GRS
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(covar,grs_col)], collapse = "+"),
                          sep = ""))
grs<-glm(formula=formula,data=qsub,family="binomial")
grs_fitted<-fitted(grs)

## model for FHIGRS
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(covar,fhigrs_col)], collapse = "+"),
                          sep = ""))
fhigrs<-glm(formula=formula,data=qsub,family="binomial")
fhigrs_fitted<-fitted(fhigrs)

## model for family history + GRS with interaction term
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(covar)], collapse = "+"), "+" ,
                          paste(colnames(qsub)[c(strat_col,grs_col)],collapse="*"),
                          sep = ""))
int<-glm(formula=formula,data=qsub,family="binomial")
int_fitted<-fitted(int)

## model for additive GRS + family history without interaction term
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(covar,strat_col,grs_col)], collapse = "+"),
                          sep = ""))
add<-glm(formula=formula,data=qsub,family="binomial")
add_fitted<-fitted(add)

## model for family history
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[c(covar,strat_col)], collapse = "+"),
                          sep = ""))
fh<-glm(formula=formula,data=qsub,family="binomial")
fh_fitted<-fitted(fh)

## intercept only model
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~","1"))
null<-glm(formula=formula,data=qsub,family="binomial")
null_fitted<-fitted(null)

## reduced model
formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                          paste(colnames(qsub)[covar],collapse="+"),
                          sep=""))
red<-glm(formula=formula,data=qsub,family="binomial")
red_fitted<-fitted(red)


##use functions from ROCR package using the user defined function ROCR_package
fhigrs_df<-data.frame(pred=fhigrs_fitted,label=qsub[[pheno_col]])
grs_df<-data.frame(pred=grs_fitted,label=qsub[[pheno_col]])
add_df<-data.frame(pred=add_fitted,label=qsub[[pheno_col]])
int_df<-data.frame(pred=int_fitted,label=qsub[[pheno_col]])
fh_df<-data.frame(pred=fh_fitted,label=qsub[[pheno_col]])
null_df<-data.frame(pred=null_fitted,label=qsub[[pheno_col]])
red_df<-data.frame(pred=red_fitted,label=qsub[[pheno_col]])
                       
grs_roc<-ROCR_package(grs_df)
fhigrs_roc<-ROCR_package(fhigrs_df)
add_roc<-ROCR_package(add_df)
int_roc<-ROCR_package(int_df)
fh_roc<-ROCR_package(fh_df)
null_roc<-ROCR_package(null_df)
red_roc<-ROCR_package(red_df)

##label models
grs_roc$method<-"GRS"
fhigrs_roc$method<-"FHiGRS"
add_roc$method<-"GRS+FH"
int_roc$method<-"GRS*FH"
fh_roc$method<-"FH"
null_roc$method<-"null"
red_roc$method<-"reduced"
roc_df<-rbind(grs_roc,fhigrs_roc,add_roc,int_roc,fh_roc,null_roc,red_roc)

##pull out AUC
grs_auc<-unique(roc_df[roc_df$method=="GRS",]$auc)
fhigrs_auc<-unique(roc_df[roc_df$method=="FHiGRS",]$auc)
add_auc<-unique(roc_df[roc_df$method=="GRS+FH",]$auc)
int_auc<-unique(roc_df[roc_df$method=="GRS*FH",]$auc)
fh_auc<-unique(roc_df[roc_df$method=="FH",]$auc)
null_auc<-unique(roc_df[roc_df$method="null",]$auc)
red_auc<-unique(roc_df[roc_df$method=="reduced",]$auc)


##plot
pdf_fn<-paste(sep=".",out,"withcovar_ROC.pdf")
pdf(file=pdf_fn,height=3,width=4,useDingbats=FALSE)
print(ggplot(roc_df,aes(x=x,y=y,color=method)) + theme_bw() +geom_line() +
      coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
      labs(title=main,x="False Positive Rate",y="True Positive Rate") +
      annotate("text",x=0.8,y=0, label=paste0("GRS AUC  ",format(grs_auc,digits=dig,format="f")),color="darkgrey",size=2) +
      annotate("text",x=0.8,y=0.05,label=paste0("FHiGRS AUC  ",format(fhigrs_auc,digits=dig,format="f")),color="darkblue",size=2) +
      annotate("text",x=0.8,y=0.1,label=paste0("GRS + FH AUC ",format(add_auc,digits=dig,format="f")),color="seagreen4",size=2) +
      annotate("text",x=0.8,y=0.15,label=paste0("GRS*FH AUC",format(int_auc,digits=dig,format="f")),color="purple",size=2) +
      annotate("text",x=0.8,y=0.2,label=paste0("Covariates AUC",format(red_auc,digits=dig,format="f")),color="red",size=2) +
      annotate("text",x=0.8,y=0.25,label=paste0("FH AUC",format(fh_auc,digits=dig,format="f")),color=
      geom_abline(slope=1,intercept=0,linetype="dashed",color="black"))
dev.off()
