#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
#library(ggridges)
library(gridExtra)
library(RNOmni)
library(optparse)
library(ResourceSelection)


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
  make_option(c("-c","--cut_points"),type="character",help="Comma separated list of percentile values in decimal format with which to compare top and bottom of distribution [default=0.8,0.9,0.95,0,0.99,0.995]",default="0.8,0.9,0.95,0.99,.995"),
  make_option(c("-q","--quantiles"),type="character",help="Comma separated list of q-quantiles for binning the GRS distribution (e.g. 4,10,20 gives quartile, decile, ventile) [default=4,5,10,20,100]",default="4,5,10,20,100"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-v","--covariates"),type="character",default="",help="Comma separated list of 1-based column incides for model covariates"),
  make_option(c("-i","--invNorm"),type="logical",default=FALSE,help="Inverse normalize GRS for entire population [default=FALSE]"),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--xlabel",type="character",default="GRS",help="X-axis label [default='']"),
  make_option("--ylabel",type="character",default="Prevalence",help="Y-axis label [default='']"),
  make_option("--legend",type="character",default="Binary stratum",help="Legend title which is stratum [default='Binary stratum']"),
  make_option("--codeDir",type="character",default="/FHiGRS_score/",help="Directory for repository for sourcing other code in code base [default=/FHiGRS_score/]")
)

parser <- OptionParser(
  usage="%prog --file --stratum_col --pheno_col --grs_col --cut_points --quantile",
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
quantiles<-as.numeric(strsplit(arguments$options$quantiles,",")[[1]])
cutpts<-as.numeric(strsplit(arguments$options$cut_points,",")[[1]])
out<-arguments$options$output
main<-arguments$options$maintitle
xlabel<-arguments$options$xlabel
ylabel<-arguments$options$ylabel
legend<-arguments$options$legend
dig<-arguments$option$digits
covar<-as.numeric(strsplit(arguments$options$covariates,",")[[1]])
invNorm<-arguments$options$invNorm
header<-arguments$options$header
##source relevant code from code base
#source(paste0(arguments$options$codeDir,"helperFunctions.R")) ##will be used to calculate FHiGRS
#TO do: get helper functions working and pare this script down

###########################################################
#################### FUNTIONS #############################
###########################################################

##quantile is how many divisions we should make of dataset (e.g. 20-quantile is ventiles for 20 groups)
##returns object with prevalences and standard error for distribution 
##split into provided number of quantiles and then by provided stratum if qfirst==TRUE
##TO DO: checking for binary variable isn't working as expected
prev_per_quantile_stratum<-function(df,GRS_col,prev_col,strat_col,qtile,qfirst=FALSE){
  if (!sum(unique(df[[strat_col]])==c(0,1))==2) {
    print("Stratum column must be a binary variable. Expects 0 and 1.")
  }
  if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
    print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
  }
  if (sum(qtile)<2*length(qtile)){ #check qtile
    print("q-quantiles should be number of divisions for data set and must be greater than 1")
  } 
  #initialize data structures
  p<-(100/qtile)/100
  index<-c(seq(from=0,to=1,by=p)*100)
  prevalences<-matrix(NA,2,qtile+1) #initialize prevalence matrix
  ns<-matrix(NA,2,qtile+1) #initialize count matrix
  ses<-matrix(NA,2,qtile+1)#initialize se matrix
  
  if (qfirst==TRUE){ #calculate quantiles before stratification
    tiles<-quantile(df[[GRS_col]],seq(from=0,to=1,by=p)) #quantile values
    for (i in 1:length(index)-1) {
      for (r in c(1,2)){ #iterate over stratum
        prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] <= tiles[i+1] & df[[strat_col]]==(r-1)][[prev_col]]
        prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ns[r,i]<-length(prev_list)
        ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
      }
    }
  } else { #stratify before calculating quantiles
    for (r in c(1,2)){ #iterate over stratum
      sdf<-df[df[[strat_col]]==(r-1),] #make data set for one stratum
      tiles<-quantile(sdf[[GRS_col]],seq(from=0,to=1,by=p)) #quantile values
      for (i in 1:length(index)-1) {
        prev_list<-sdf[sdf[[GRS_col]] > tiles[i] & sdf[[GRS_col]] <= tiles[i+1]][[prev_col]]
        prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ns[r,i]<-length(prev_list)
        ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
      } 
    }
  }
  #create object 
  pqs<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
  class(pqs)<-"prev_quantile_stratum_obj"
  return(pqs)
}


###### model function
## use logistic regression to model disease ~ GRS + covariates
## write out table of gooodness of fit
model<-function(df,grs_col,fhigrs_col,strat_col,pheno_col,covar,out){
  mobj<-list() #initialize object
  ## inverse rank normalize GRS in the population because FHIGRS is inv normal scale
  df$invNormGRS<-rankNorm(df[[grs_col]])
  grs_col<-which(names(df)=="invNormGRS") #new GRS col
  
  #model for GRS
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  gdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"GRS")) #standard GRS only
  mobj[["GRS"]]<-gdf
  grs<-glm.obj
  
  ## model for FHIGRS
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                              paste(colnames(df)[c(covar,fhigrs_col)], collapse = "+"),
                              sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FHIGRS"))
  mobj[["FHiGRS"]]<-fdf
  fhigrs<-glm.obj
  
  ## model for family history + GRS with interaction term
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar)], collapse = "+"), "+" ,
                            paste(colnames(df)[c(strat_col,grs_col)],collapse="*"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  idf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH","GRS","FH*GRS"))
  mobj[["interaction"]]<-idf
  int<-glm.obj
  
  ## model for additive GRS + family history + GRS without interaction term
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  adf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH","GRS"))
  mobj[["additive"]]<-adf
  add<-glm.obj
  
  ## model for family history 
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH"))
  mobj[["family"]]<-fdf
  fh<-glm.obj
  
  ## intercept only model 
  formula<-as.formula(paste(colnames(df)[pheno_col], "~","1"))
  null.obj<-glm(formula=formula,data=df,family="binomial")
  
  ## saturated model
  #a <- factor(1:nrow(df))
  #fit <- glm(df[[pheno_col]]~a,family=binomial)
  
  ## compare goodness of model fit 
  model_list<-c("grs","fhigrs","int","add","fh")
  fit<- matrix(ncol=6, nrow=length(model_list))
  for (i in 1:length(model_list)){
    fit[i,1]<-get(model_list[i])$deviance
    fit[i,2]<-get(model_list[i])$null.deviance
    fit[i,3]<-get(model_list[i])$aic
    fit[i,4]<-hoslem.test(df[[pheno_col]], fitted(get(model_list[i])))$p.value
    
    #Cox-Snell pseudo R2 and max adjusted R2
    l0=as.numeric(logLik(null.obj))
    l1=as.numeric(logLik(get(model_list[i])))
    R2_cox<-1-(exp(l0 - l1))^{2/nrow(df)}
    R2_max<-1-exp(l0)^{2/nrow(df)}
    max_adj_R2<-R2_cox/R2_max
    fit[i,5]<-R2_cox
    fit[i,6]<-max_adj_R2
  }
  fit_df<-data.frame(fit)
  names(fit_df)<-c("deviance","null_deviance","aic","hoslem","cox_snell_r2","max_adj_r2")
  fit_df$model<-c("GRS","FHiGRS","GRS*FH","GRS+FH","FH")
  file_n<-paste(sep=".",out,"modelGOF.txt")
  write.table(format(fit_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
  
  ## anova comparing to additive model
  anova_df<-rbind(
    data.frame(dev=anova(add,grs,test="LRT")$Deviance[2],pval=anova(add,grs,test="LRT")$`Pr(>Chi)`[2],method="GRS"),
    data.frame(dev=anova(add,fh,test="LRT")$Deviance[2],pval=anova(add,fh,test="LRT")$`Pr(>Chi)`[2],method="FH"),
    data.frame(dev=anova(add,fhigrs,test="LRT")$Deviance[2],pval=anova(add,fhigrs,test="LRT")$`Pr(>Chi)`[2],method="FHiGRS"))
  file_n<-paste(sep=".",out,"compareToAdditive.txt")
  write.table(format(anova_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
  
  
  class(mobj)<-"model_obj"
  return(mobj)

}


## use logistic regression to model disease ~ I(top of GRS distribution) + covariates
model_indicator<-function(df,value,fhigrs_col,grs_col,strat_col,pheno_col,covar,qfirst=FALSE){
  if (value < 0.5 | value >= 1){
    print("Percentile for dividing GRS/FHiGRS distribution must be >= 0.5 < 1 ")
  }
  mobj<-list() #initialize object
  
  ##model I(GRS distribution) for all data
  #to do: inv normalize GRS?
  q<-quantile(df[[grs_col]],value) #quantile in all data
  #top of bottom of distribution as binary variable, bottom as reference group
  df$dist<-2
  df$dist<-ifelse(df[[grs_col]]>q,1,0) #1 if in top, 0 in bottom
  df[df$dist==2]$dist<-NA
  dist_col<-which(names(df)=="dist")
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,dist_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  gdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")) #all
  m<-matrix(table(df$dist,df[[pheno_col]]),byrow=FALSE,nrow=2)
  gdf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  gdf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["GRS"]]<-gdf
  
  ##if we want to stratify before we decide on quantiles, remake I(GRS distribution) column, otherwise use previous
  if (qfirst!=TRUE){
    df$dist<-2
    for (s in c(1,2)){
      q<-quantile(df[df[[strat_col]]==(s-1)][[grs_col]],value,na.rm=TRUE)
      df[df[[strat_col]]==(s-1)]$dist<-ifelse(df[df[[strat_col]]==(s-1)][[grs_col]]>q,1,0) #1 if in top, 0 in bottom
    } 
    df[df$dist==2]$dist<-NA
  }
  dist_col<-which(names(df)=="dist")
  
  ## model I(dist) per stratum
  for (s in c(1,2)){ #iterate over stratum
    formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                              paste(colnames(df)[c(covar,dist_col)], collapse = "+"),
                              sep = ""))
    glm.obj<-glm(formula=formula,data=df[df[[strat_col]]==(s-1)],family="binomial") #subsets data frame to stratum
    df_name<-paste0("stratum",s-1)
    m<-matrix(table(df[df[[strat_col]]==(s-1)]$dist,df[df[[strat_col]]==(s-1)][[pheno_col]])/nrow(df[df[[strat_col]]==(s-1)]),byrow=FALSE,nrow=2)
    top<-m[2,2]/sum(m[2,])  #prevalence in top of distribution
    bottom<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
    assign(df_name,cbind(data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")),"top_prev"=top,"bottom_prev"=bottom))
    mobj[[df_name]]<-get(df_name)
  }
  
  ## model using I(dist & stratum=1)
  df$conditional<-2
  df[df[[strat_col]]==1 & df$dist==1]$conditional<-1
  df[df[[strat_col]]==1 & df$dist==0]$conditional<-0
  df[df[[strat_col]]==0]$conditional<-0
  df[df$conditional==2]$conditional<-NA
  df$dist<-df$conditional #replace I(GRS dist) with I(GRS dist & stratum=1)
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,dist_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  cdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist"))
  m<-matrix(table(df$dist,df[[pheno_col]]),byrow=FALSE,nrow=2)
  cdf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  cdf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["conditional"]]<-cdf
  
  
  ##model using I(FHiGRS dist)
  q<-quantile(df[[fhigrs_col]],value) #quantile in all data
  #top of bottom of distribution as binary variable, bottom as reference group
  df$fhigrs<-2
  df$fhigrs<-ifelse(df[[fhigrs_col]]>q,1,0) #1 if in top, 0 in bottom
  df[df$fhigrs==2]$dist<-NA
  df$dist<-df$fhigrs #replace I(GRS dist) with I(FHiGRS dist)
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,dist_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")) #all
  m<-matrix(table(df$dist,df[[pheno_col]]),byrow=FALSE,nrow=2)
  fdf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  fdf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["FHIGRS"]]<-fdf
  
  ##model using I(FH) NOT REALLY A CUT POINT JUST POSITIVE OR NEGATIVE
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  ydf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")) 
  m<-matrix(table(df[[strat_col]],df[[pheno_col]]),byrow=FALSE,nrow=2)
  ydf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  ydf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["family"]]<-ydf  
  
  class(mobj)<-"model_obj"
  return(mobj)
}


## screen N people how many do we catch and how many do we miss given screening strategy (prioritize by family history or no)
#OR and 95% confidence intervals for 2 by 2 contingency tables of cases/control in top/bottom distribution for various scenarios 
clinical_impact<-function(df,value,grs_col,fhigrs_col,pheno_col,strat_col,N=10000,qfirst=FALSE){

  ## inverse rank normalize standard GRS in the population because FHIGRS is inv normal scale
  df$invNormGRS<-rankNorm(df[[grs_col]])
  grs_col<-which(names(df)=="invNormGRS") #new GRS col
  counts<-list(rep(NA,2))
  score_list<-c(fhigrs_col,grs_col)
  for (k in c(1,2)){
    q<-quantile(df[[score_list[k]]],value) #take quantile in score of interest
    df$dist<-2
    df$dist<-ifelse(df[[score_list[k]]]>q,1,0) #1 if in top, 0 in bottom
    df[df$dist==2]$dist<-NA
    dist_col<-which(names(df)=="dist")
    counts[[k]]<-as.matrix(table(df[[pheno_col]],df[[dist_col]])) #2 by 2 table of binary phenotype and binary top/bottom of distribution
  }
  
  #first row of matrix is counts from bottom of distribution at given cut point
  #second row of matrix is counts from top of distirbution at given cut point
  #first column of matrix is counts of controls
  #second column of matrix is counts of cases 
  
  #FHIGRS
  screen<-c(counts[[1]][2,1],counts[[1]][2,2]) #control,case of top distribution
  no_screen<-c(sum(counts[[1]][1,1],counts[[1]][,1]),sum(counts[[1]][1,2],counts[[1]][,2])) #control,case of bottom distribution
  m<-matrix(c(no_screen,screen),nrow=2,ncol=2,byrow=TRUE)
  m_frac<-m/sum(m)
  scenario1<-m_frac*N
  
  ##GRS
  screen<-c(counts[[2]][2,1],counts[[2]][2,2]) #control,case of top distribution
  no_screen<-c(sum(counts[[2]][1,1],counts[[2]][,1]),sum(counts[[2]][1,2],counts[[2]][,2])) #control,case of bottom distribution
  grsm<-matrix(c(no_screen,screen),nrow=2,ncol=2,byrow=TRUE)
  grsm_frac<-grsm/sum(grsm)
  scenario2<-grsm_frac*N #no screen is top row, screen is bottom row, control is first column, case is second column
  
  #FHIGRS relative to GRS
  false_pos<-(scenario1-scenario2)[2,1]
  false_neg<-(scenario1-scenario2)[1,2]
  
  ## scenario 1 FHiGRS
  neg_predictive<-m[1,1]/sum(m[1,])
  pos_predictive<-m[2,2]/sum(m[2,])
  specificity<-m[1,1]/sum(m[,1])
  sensitivity<-m[2,2]/sum(m[,2])
  accuracy<-(m[1,1] + m[2,2])/sum(m) 
  top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  OR<-(m[2,2]/m[2,1])/(m[1,2]/m[1,1])  #case/control top distribution over case/control bottom distribution
  SE<-sqrt(sum(1/m)) #log odds scale
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  scenario1_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy, OR=OR,SE=SE, UB=UB, LB=LB,top_prev,bottom_prev,scenario="FHiGRS",cutpt=value)
  
  ## scenario 2 GRS
  m<-grsm
  neg_predictive<-m[1,1]/sum(m[1,])
  pos_predictive<-m[2,2]/sum(m[2,])
  specificity<-m[1,1]/sum(m[,1])
  sensitivity<-m[2,2]/sum(m[,2])
  accuracy<-(m[1,1] + m[2,2])/sum(m)
  top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  OR<-(m[2,2]/m[2,1])/(m[1,2]/m[1,1])  #case/control top distribution over case/control bottom distribution
  SE<-sqrt(sum(1/m)) #log odds scale
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  scenario2_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy,  OR=OR,SE=SE, UB=UB, LB=LB, top_prev, bottom_prev, scenario="GRS",cutpt=value)
  
  #scenario 3 GRS stratified by FH
  #top dist is > quantile and FH=1, reference group is < quantile and FH=1 +  all of FH=0
  index<-c(value*100,100)
  prevalences<-matrix(NA,2,length(index)) #initialize prevalence matrix
  counts<-list(matrix(NA,length(index),length(index)),matrix(NA,length(index),length(index))) #list of matrices, one matrix per stratum
  ses<-matrix(NA,2,length(index))#initialize se matrix
  if (qfirst==TRUE) { #make quantiles before stratifying data
    tiles<-quantile(df[[grs_col]],c(0,value,1))
    for (i in 1:length(index)) {
      for (r in c(1,2)){ #iterate over stratum
        prev_list<-df[df[[grs_col]] > tiles[i] & df[[grs_col]] <= tiles[i+1] & df[[strat_col]]==(r-1)][[pheno_col]]
        prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ses[r,i]<-sqrt((prevalences[i]*(1-prevalences[i]))/length(prev_list)) #what is SE for this prevalence
        counts[[r]][i,]<-as.vector(table(prev_list)) #counts for OR
    }}
  } else {
      for (r in c(1,2)){ #iterate over stratum
        sdf<-df[df[[strat_col]]==(r-1),] #make data set for one stratum
        tiles<-quantile(sdf[[grs_col]],c(0,value,1)) #quantile values
        for (i in 1:length(index)) {
          prev_list<-sdf[sdf[[grs_col]] > tiles[i] & sdf[[grs_col]] <= tiles[i+1]][[pheno_col]]
          prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
          ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
          counts[[r]][i,]<-as.vector(table(prev_list))
        }}}
  num<-counts[[2]][2,2]/counts[[2]][2,1] #case/control of top distribution with stratum=1
  denom<-sum(counts[[2]][1,2],counts[[1]][,2])/sum(counts[[2]][1,1],counts[[1]][,1]) #case/control of top distribution with stratum=0 + bottom distribution stratum=0|1
  OR<-num/denom
  SE<-sqrt(sum(1/counts[[2]][2,2],1/counts[[2]][2,1],1/sum(counts[[2]][1,2],counts[[1]][,2]),1/sum(counts[[2]][1,1],counts[[1]][,1])))
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  top_prev<-counts[[2]][2,2]/sum(counts[[2]][2,]) #cases in top dist over all in top
  bottom_prev<-sum(counts[[1]][,2],counts[[2]][1,2])/(sum(counts[[2]][1,])+sum(counts[[1]])) #cases in bottom dist over all in bottom
  
  #FH+GRS relative to GRS
  strat<-matrix(c(sum(counts[[2]][1,2],counts[[1]][,2]),sum(counts[[2]][1,1],counts[[1]][,1]),counts[[2]][2,2],counts[[2]][2,1]),nrow=2,ncol=2,byrow=TRUE)
  stratfrac<-strat/sum(strat)
  scenario3<-stratfrac*N #no screen is top row, screen is bottom row, control is first column, case is second column
  false_pos<-(scenario3-scenario2)[2,1]
  false_neg<-(scenario3-scenario2)[1,2]
  
  scenario3_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy, OR=OR,SE=SE, UB=UB, LB=LB, top_prev, bottom_prev, scenario="GRS|FH",cutpt=value)
  
  #scenario 4 just family history
    m<-table(df[[strat_col]],df[[pheno_col]])
  neg_predictive<-m[1,1]/sum(m[1,])
  pos_predictive<-m[2,2]/sum(m[2,])
  specificity<-m[1,1]/sum(m[,1])
  sensitivity<-m[2,2]/sum(m[,2])
  accuracy<-(m[1,1] + m[2,2])/sum(m)
  top_prev<-m[2,2]/sum(m[2,]) #prevalence in family history positive
  bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in family history negative
  OR<-(m[2,2]/m[2,1])/(m[1,2]/m[1,1])  #case/control top distribution over case/control bottom distribution
  SE<-sqrt(sum(1/m)) #log odds scale
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  
  #FH relative to GRS
  fhfrac<-m/sum(m)
  scenario4<-fhfrac*N #no screen is top row, screen is bottom row, control is first column, case is second column
  false_pos<-(scenario4-scenario2)[2,1]
  false_neg<-(scenario4-scenario2)[1,2]
  
  scenario4_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy, OR=OR,SE=SE, UB=UB, LB=LB, top_prev, bottom_prev, scenario="FH",cutpt=0)

  return(rbind(scenario1_df,scenario2_df,scenario3_df,scenario4_df))
}

##function to estimate FHiGRS
#takes data frame created from prev_per_quantile_stratum object
estimate_FHiGRS<-function(prev_df,main_df,strat_col,grs_col){
  ##put q-quantiles in order of prevalence
  ##to do: use prevalence from reference population
  prev_df_order<-prev_df[order(prev_df$prev),]
  prev_df_order$rank<-seq(100,nrow(prev_df)*100,100) #make rankings in scale of 100
  
  ##assign samples to groups
  main_df$rank<-as.numeric(1) #initialize column for FHiGRS
  for (s in c(1,2)){
    prev_df_order_strat<-prev_df_order[prev_df_order$stratum==(s-1),]
    prev_df_order_strat[1,'lower_tile']<- -100 #condition to put samples equiv to minimum in bottom bin
    for (j in 1:nrow(prev_df_order_strat)){
      u<-prev_df_order_strat[j,'upper_tile']
      l<-prev_df_order_strat[j,'lower_tile']
      rank<-prev_df_order_strat[j,'rank']
      main_df[main_df[[strat_col]]==(s-1) & main_df[[grs_col]]>l & main_df[[grs_col]]<=u,'rank']<-rank
    }
  }
  main_df$grs_rank<-main_df[[grs_col]] + main_df$rank #add rank to GRS
  main_df$FHIGRS<-qnorm((rank(main_df$grs_rank,na.last="keep")-0.5)/sum(!is.na(main_df$grs_rank))) #new FHIGRS
  return(main_df)
}

  
###########################################################
#################### MAIN #################################
###########################################################

##read data 
dat<-fread(file,header=header)
print(paste("Data dimensions are:",dim(dat)[1],dim(dat)[2]))
## To DO: check column assumptions 

##subset to data with stratum available
subset<-dat[!is.na(dat[[strat_col]])] #remove if NA for stratum
print(paste("Data dimensions after removing samples with NA stratum:",dim(subset)[1],dim(subset)[2]))

##stratify then calculate quantiles
sobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=FALSE)
#print(sobj)

size<-length(sobj) #number of quantiles being tested, size of obj
for (i in 1:size){ #across q-quantiles
  for (j in c(1,2)){ #across stratum
    list_length<-length(sobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-sobj[[i]]$prev[j,][-list_length]
    se<-sobj[[i]]$se[j,][-list_length]
    n<-sobj[[i]]$n[j,][-list_length]
    tiles<-sobj[[i]]$tiles[-1]
    lower_tile<-sobj[[i]]$tiles[1:list_length-1]
    upper_tile<-sobj[[i]]$tiles[2:list_length]
    percents<-names(tiles)
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(j-1,list_length-1)
    if (j==1) {
      sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
    } else {
      sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
    }
  }
  
  #estimate FHiGRS
  qsub<-estimate_FHiGRS(sdf,subset,strat_col,grs_col)
  fhigrs_col<-which(names(qsub)=="FHIGRS")
  
  ## compare GRS and FHIGRS with logistic regression and covariates
  mobj<-model(df=qsub,grs_col=grs_col,fhigrs_col=fhigrs_col,pheno_col=pheno_col,strat_col=strat_col,covar=covar,out=paste(sep=".",out,quantiles[i]))
  model_list<-c("GRS","FHiGRS","interaction","additive","family") #model being tested 
  for (l in 1:length(mobj)){ #loop over list from model function
    mobj[[names(mobj[l])]]$model<-names(mobj)[l]
    mobj[[l]]$OR<-exp(mobj[[l]][,'Estimate'])
    mobj[[l]]$LB<-exp(mobj[[l]][,'Estimate']-mobj[[l]][,'Std..Error'])
    mobj[[l]]$UB<-exp(mobj[[l]][,'Estimate']+mobj[[l]][,'Std..Error'])
    mobj[[l]]$model<-model_list[l]
    mobj[[l]]$pred<-row.names(mobj[[l]])
  }
  d<-do.call("rbind",mobj)
  d$qtile=quantiles[i]
  
  #pull out important covariate/model combinations
  dsub<-d[d$pred %in% c("FHIGRS","GRS","FH","FH*GRS"),]
  dsub$model<-as.factor(dsub$model)
  dsub$model<-factor(dsub$model,levels=c("GRS", "family", "additive","FHiGRS","interaction"))
  dsub$model<-revalue(dsub$model,c("family"="FamilyHistory"))

  if (i==1) {
      cmodel<-dsub    
  } else {
      cmodel<-rbind(cmodel,dsub)
  }
  
  ## compare GRS, FHIGRS, GRS+FH, and FH with 2 by 2 contingency tables, also clinical impact and false negatives/positives/accuracy 
  logical_list<-c(TRUE,FALSE)
  label_list<-c("quantileFirst","stratifyFirst")
  for (qlogic in c(1,2)){ #over 2 logic options for creating quantiles 
    #2 by 2 tables
    obj<-lapply(cutpts,clinical_impact,df=qsub,pheno_col=pheno_col,grs_col=grs_col,fhigrs_col=fhigrs_col,strat_col=strat_col,qfirst=logical_list[qlogic])
    if (qlogic==1 & i==1) {
      clin_df<-data.frame(bind_rows(obj),qfirst=label_list[qlogic],qtile=quantiles[i])
    } else {
      clin_df<-rbind(clin_df,data.frame(bind_rows(obj),qfirst=label_list[qlogic],qtile=quantiles[i]))
    }
  }
  
  ## compare I(distribution of GRS and FHIGRS) with logistic regression and covariates
  for (qlogic in c(1,2)){ #over 2 logic options for creating quantiles 
    mobj<-lapply(cutpts,model_indicator,df=qsub,covar=covar,pheno_col=pheno_col,grs_col=grs_col,fhigrs_col=fhigrs_col,strat_col=strat_col,qfirst=logical_list[qlogic])
    for (l in 1:length(mobj)){ #loop over list from model function, one per cutpt value
      name_list<-names(mobj[[l]])
      for (n in 1:length(name_list)){ #across the models tested 
        colnames<-lapply(mobj[[l]][name_list[n]],names)[[1]] #pull column names from data frame
        z<-data.frame(mobj[[l]][name_list[n]])['dist',] #pull out summ stats for variable of interest
        names(z)<-colnames
        OR=exp(z[['Estimate']])
        LB=exp(z[['Estimate']]-1.96*z[['Std..Error']])
        UB=exp(z[['Estimate']]+1.96*z[['Std..Error']])
        pval=z[['Pr...z..']]
        if (qlogic==1 & n==1 & i==1 & l==1) { #inititalize df for first q-quantile, first cutpoint, first quantile logic
         model_df<-cbind(z,qtile=quantiles[i],cutpt=cutpts[l],qfirst=label_list[qlogic],name=name_list[n],OR=OR,LB=LB,UB=UB,pval=pval)
        } else {
          model_df<-rbind(model_df,cbind(z,qtile=quantiles[i],cutpt=cutpts[l],qfirst=label_list[qlogic],name=name_list[n],OR=OR,LB=LB,UB=UB,pval=pval))
        }
      }
    }
  }

  #plot of false positives and false negatives in stratified versus regular screening schemes, convert decimal cutpoints to whole numbers
  pdf_fn<-paste(sep=".",out,"clinicalImpact.pdf")
  nudge_factor<- diff(range(clin_df$false_pos))/10 #if the x axis scale is kind of small we don't need to nudge labels too far from points
  n<-length(quantiles)
  if (n<=5){
    scale<-100*as.vector(cutpts) #use all cutpts as breaks in legend
  } else {
    scale<-100*(c(cutpts[1],cutpts[floor(n/2)],cutpts[n])) #use max, middle, and min as cutpts in legend 
  }
  pdf(file=pdf_fn,height=4,width=6,useDingbats=FALSE)
  print(ggplot(clin_df,aes(x=false_pos,y=false_neg,label=cutpt*100,size=100*cutpt)) + geom_point(alpha=0.7,color="orchid4") +   
          theme_bw() + geom_text(nudge_x=nudge_factor,color="black",size=5) + 
          labs(title="FHiGR score versus standard GRS",x="False Positive",y="False Negative") + guides(color=FALSE) +
          geom_vline(xintercept=0,linetype="dashed",color="black",alpha=0.5) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.5) +
          scale_size_continuous(name="Cut Point",breaks=scale))
  dev.off()
 
}

####### Write tables and make plots (plots have fewer things visualized than exist in the tables)


model_df$name<-revalue(model_df$name, c("stratum0"="FamilyHistory-", "stratum1"="FamilyHistory+","conditional"="ConditionalFamilyHistory+","family"="FamilyHistory","FHIGRS"="FHiGRS"))
model_df_sub<-model_df[model_df$name=="GRS" | model_df$name=="FHiGRS" | model_df$name=="FamilyHistory",] # remove conditional family history since its a smaller group of people and not equal comparison
#write table comparing logistic regression with indicator variable for top/bottom distribution across # of bins to divide data for FHIGRS
file_n<-paste(sep=".",out," modelIndicator.compareScores.txt")
write.table(format(model_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
##plot comparison, make line for family history OR since it is not a threshold
by(model_df_sub,model_df_sub$qtile,
   function(x){
     name=unique(x$qtile)
     pdf_fn<-paste(sep=".",out,name,"modelIndicator.compareScores.pdf")
     pdf(file=pdf_fn,height=4,width=6,useDingbats=FALSE)
     print(ggplot(x[x$name!="FamilyHistory",],aes(x=OR,y=as.factor(cutpt),color=name))  + geom_point() + theme_bw() + facet_wrap(~qfirst) + 
             geom_errorbarh(aes(xmin=x[x$name!="FamilyHistory",]$LB, xmax=x[x$name!="FamilyHistory",]$UB)) +
             labs(title=paste(main,"Model with covariates and high risk threshold"),y="Threshold for High Risk Group",x="Odds Ratio between High Risk & Reference Group") + scale_color_manual(values=c("grey","orchid4"),name="") +
             geom_vline(linetype="dashed",color="black",xintercept=1,alpha=0.7) + geom_vline(color="red",xintercept=x[x$name=="FamilyHistory",]$OR,alpha=0.2))
     dev.off()
   })

#write table comparing scores from logistic regression across # of bins to divide data for FHIGRS
file_n<-paste(sep=".",out," modelContinuous.compareScores.txt")
write.table(format(cmodel,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
print(head(cmodel))
##plot comparison
by(cmodel, cmodel$qtile,
   function(x){
     name=unique(x$qtile)
     pdf_fn<-paste(sep=".",out,name,"modelContinuous.compareScores.pdf")
     pdf(file=pdf_fn,height=6,width=6,useDingbats=FALSE)
     print(ggplot(x,aes(x=OR,y=as.factor(pred),color=model)) + facet_wrap(~model,ncol=1,scales="free_y") + geom_point() + theme_bw() + geom_errorbarh(aes(xmin=x$LB,xmax=x$UB)) + 
           labs(title=paste(main,"Model with covariates and continuous variables"),y="Predictor from model",x="Odds Ratio") +
           scale_color_manual(values=c("grey","darkblue","goldenrod3","orchid4","seagreen4"),name="") +
           geom_vline(linetype="dashed",color="black",xintercept=1,alpha=0.7))
     dev.off()
   })


#write table comparing scores from 2 by 2 contingency tables for top/bottom distribution across # of bins to divide data for FHIGRS
file_n<-paste(sep=".",out,"table.compareScores.txt")
write.table(format(clin_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
clin_df$scenario<-factor(clin_df$scenario,levels(clin_df$scenario)[c(2,1,3,4)])
print(head(clin_df))
#plot comparison 
by(clin_df, clin_df$qtile,
   function(x){
       name=unique(x$qtile)
       pdf_fn<-paste(sep=".",out,name,"table.compareScores.pdf")
       pdf(file=pdf_fn,height=4,width=6,useDingbats=FALSE)
       print(ggplot(x[x$scenario!="FH" & x$scenario!="GRS|FH",],aes(y=as.factor(cutpt),x=OR,color=scenario)) + facet_wrap(~qfirst) + geom_point() + theme_bw() +
             geom_errorbarh(aes(xmin=x[x$scenario!="FH" & x$scenario!="GRS|FH",]$LB,xmax=x[x$scenario!="FH" & x$scenario!="GRS|FH",]$UB)) +
             labs(title=paste(main, "Reduced model with high risk threshold"),y="Threshold for High Risk Group",x="Odds Ratio between High Risk & Reference Group") + scale_color_manual(values=c("grey","orchid4"),name="") +
             geom_vline(linetype="dashed",color="black",xintercept=1,alpha=0.7) +
             geom_vline(color="red",xintercept=x[x$scenario=="FH",]$OR,alpha=0.2)
             )
    dev.off()
})




