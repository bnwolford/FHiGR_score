#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RNOmni)
library(optparse)
library(ResourceSelection)
library(Hmisc)
library(corrplot)
library(GenomicRanges)
library(ltm)
library(pROC)
library(RColorBrewer)
library(DescTools) 
library(rcompanion)
library(LaCroixColoR)
library(cvAUC)
library(grid)

print(Sys.time())
print(sessionInfo())
###########################################################
################### Read Command Line Parameters ##########
###########################################################

optionList <- list(
  make_option(c("-f","--file"),type="character",help="File with sample IDs, phenotypes, self reported family history, and GRS. Expeects header. White space delimited.",default=NULL),
  make_option(c("-s","--stratum_col"),type="numeric",help="1-based column with stratum (e.g. family history). Must be on binary scale 1/0 with 1 being affirmative. NAs ok.",default=NULL),
  make_option(c("-p","--pheno_col"),type="numeric",help="1-based column with phenotype (e.g. disease status). Must be on binary scale 1/0 with 1 being case. NAs ok.",default=NULL),
  make_option(c("-g","--grs_col"),type="numeric",help="1-based column with GRS, not inverse normalized.",default=NULL),
  make_option(c("-c","--cut_points"),type="character",help="Comma separated list of percentile values in decimal format with which to compare top and bottom of distribution [default=0.8,0.9,0.95,0,0.99,0.995]",default="0.8,0.9,0.95,0.99,.995"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files",default=NULL),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-i","--invNorm"),type="logical",default=FALSE,help="Inverse normalize GRS for entire population [default=FALSE]"),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--xlabel",type="character",default="GRS",help="X-axis label [default='']"),
  make_option("--ylabel",type="character",default="Prevalence",help="Y-axis label [default='']"),
  make_option("--legend",type="character",default="Binary stratum",help="Legend title which is stratum [default='Binary stratum']"),
  make_option("--codeDir",type="character",default="/FHiGRS_score/",help="Directory for repository for sourcing other code in code base [default=/FHiGRS_score/]"),
  make_option("--sex_col",type="numeric",help="1-based colummn with sex",default=NULL),
  make_option(c("-b","--birthYear_col"),type="numeric",help="1-based column with birthyear",default=NULL),
  make_option(c("-a","--age_col"),type="numeric",help="1-based column with participation or enrollment age e.g. age of self reported family history",default=NULL),
  make_option(c("--binary_covar_cols"),type="character",help="Comma separated string of 1-based column numbers with binary covariates for models",default=NULL),
  make_option(c("--binary_covar_labels"),type="character",help="Comma separated string of of 1-based column labels with binary covariates for models",default=NULL),      
  make_option(c("--quant_covar_cols"),type="character",help="Comma separated string of of 1-based column numbers with quantitative covariates for models",default=NULL),
  make_option(c("--quant_covar_labels"),type="character",help="Comma separated string of of 1-based column labels with quantitative covariates for models",default=NULL)    
  
)

parser <- OptionParser(
  usage="%prog",
  option_list=optionList
)
arguments <- parse_args(parser, positional_arguments=TRUE)
print(arguments$options)

##print warnings when they happen
options(warn=1)

#check arguments without defaults
check<-c("file","stratum_col","sex_col","birthYear_col","age_col","output","pheno_col","grs_col","stratum_col")
for (c in 1:length(check)){
  if (is.null(arguments$options[[check[c]]])){
    warning(paste(sep=" ","Missing argument",check[c]))
  }
}

cutpt<-arguments$options$cut_points
file<-arguments$options$file
strat_col<-arguments$options$stratum_col
pheno_col<-arguments$options$pheno_col
sex<-arguments$options$sex_col
grs<-arguments$options$grs_col
birthYear<-arguments$options$birthYear_col
part_age<-arguments$options$age_col
out<-arguments$options$output

#covariates
if (!is.null(arguments$options$binary_covar_cols)){
  bcovar<-as.numeric(strsplit(arguments$options$binary_covar_cols,",")[[1]])
  print(bcovar)
}
if (!is.null(arguments$options$binary_covar_labels)){
  bcovar_lab<-strsplit(arguments$options$binary_covar_labels,",")[[1]]
  print(bcovar_lab)
}
if (!is.null(arguments$options$quant_covar_cols)){
  qcovar<-as.numeric(strsplit(arguments$options$quant_covar_cols,",")[[1]])
  print(qcovar)
}
if (!is.null(arguments$options$quant_covar_cols)){
  qcovar_lab<-strsplit(arguments$options$quant_covar_labels,",")[[1]]
  print(qcovar_lab)
}

#to do: censor, young

###########################################################
###################### Functions ##########################
###########################################################



#calculate odds ratio for top of dist and positive family history
odds_ratio<-function(df,qtile=0.95,prev_col,grs_col,strat_col){
  q<-quantile(df[[grs_col]],qtile) #take quantile in score of interest
  df$dist<-2
  df$dist<-ifelse(df[[grs_col]]>q & df[[strat_col]]==1,1,0) #1 if in top, 0 in bottom
  df[df$dist==2]$dist<-NA
  dist_col<-which(names(df)=="dist")
  m<-as.matrix(table(df[[prev_col]],df[[dist_col]])) #2 by 2 table of binary phenotype and binary top/bottom of distribution
  print(m)
  OR<-(m[2,2]/m[1,2])/(m[2,1]/m[1,1])  #case/control top distribution over case/control bottom distribution
  SE<-sqrt(sum(1/m)) #log odds scale
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  return(c(OR,LB,UB))
}

##### test models with outcome not indicator variable 
model<-function(df,outcome,pred,pred_labels,out,name){
  pred<-unlist(pred)
  pred_labels<-unlist(pred_labels)
  mobj<-list() #initialize object
  
  #model for covar only
  formula<-as.formula(paste(colnames(df)[outcome],"~",
                            paste(colnames(df)[c(pred)],collapse="+"),
                            sep=""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  glm_df<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",pred_labels)) #this may fail if some of the covariates aren't in the data frame, To do: make a check
  
  
  r<-ROC(df,formula,outcome,out,paste(sep="_",name)) #plot ROC
  auroc<-data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2])
  
  fit<- matrix(ncol=8, nrow=1)
  fit[1,1]<-glm.obj$deviance
  fit[1,2]<-glm.obj$null.deviance
  fit[1,3]<-glm.obj$aic
  fit[1,4]<-hoslem.test(df[[outcome]],fitted(glm.obj))$p.value
  
  #Cox-Snell pseudo R2 and max adjusted R2
  ## intercept only model 
  formula<-as.formula(paste(colnames(df)[pheno_col], "~","1"))
  null.obj<-glm(formula=formula,data=df,family="binomial")
  l0=as.numeric(logLik(null.obj))
  l1=as.numeric(logLik(glm.obj))
  R2_cox<-1-(exp(l0 - l1))^{2/nrow(df)}
  R2_max<-1-exp(l0)^{2/nrow(df)}
  max_adj_R2<-R2_cox/R2_max
  fit[1,5]<-R2_cox
  fit[1,6]<-max_adj_R2
  fit[1,7]<-PseudoR2(glm.obj,which="all")[["Nagelkerke"]]
  fit[1,8]<-BrierScore(glm.obj,scaled=FALSE) #what is the range when unscaled?!
  fit_df<-data.frame(fit)
  names(fit_df)<-c("deviance","null_deviance","aic","hoslem","cox_snell_r2","max_adj_r2","Nagelkerke","Brier")
  fit_df$cvAUC<-auroc$cvAUC
  fit_df$AUC_LB<-auroc$LB
  fit_df$AUC_UB<-auroc$UB
  fit_df$dev_test<-1-pchisq(fit_df$deviance,glm.obj$df.residual)
  
  #mobj$glm<-glm.obj
  mobj$glm<-glm_df
  mobj$fit<-fit_df
  class(mobj)<-"model_obj"
  return(mobj)
  
}


##### ROC plot function that calculates 95% CI
ROC<-function(qsub,formula,pheno_col,out,name,V=5){
  set.seed(1234)
  .cvFolds <- function(Y, V){  #Create CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  .doFit <- function(v, folds, data,formula){  #Train/test glm for each fold
    fit <- glm(formula, data=data[-folds[[v]],], family=binomial)
    pred <- predict(fit, newdata=data[folds[[v]],], type="response")
    return(pred)
  }
  folds <- .cvFolds(Y=qsub[[pheno_col]], V=V)  #Create folds
  predictions_list <- sapply(seq(V), .doFit, folds=folds, data=qsub,formula=formula)  #CV train/predict
  predictions<-unlist(predictions_list)
  predictions[unlist(folds)] <- predictions  #Re-order pred values
  # Get CV AUC and confidence interval
  ci_cv <- ci.cvAUC(predictions=predictions, labels=qsub[[pheno_col]], folds=folds, confidence=0.95)
  labels_list<-list()
  for (v in seq(V)) {   labels_list[v]<-list(qsub[unlist(folds[v])][[pheno_col]]) }
  cv<-cvAUC(predictions=predictions_list,labels=labels_list)
  
  pdf(file=paste(sep="_",out,name,"ROC_all_folds.pdf"),height=5,width=6,useDingbats = )
  plot(cv$perf,col=1:V,main=paste0(V,"-fold Cross Validated ROC"))
  abline(a=0,b=1,col="red")
  text(0.75,0.5-0.05,paste0(V,"-fold AUROC"))
  for (v in seq(V)) { text(0.75,0.5-((v+1)*0.05),format(cv$fold.AUC[v],digits=3))}
  text(0.75,0.5-(V+2)*0.05,paste0("cvAUROC=",format(ci_cv$cvAUC,digits=3)," [", format(ci_cv$ci[1],digits=3),",", format(ci_cv$ci[2],digits=3),"]"))
  dev.off()
  
  return(ci_cv)
}


#Rscript FHiGR_model.R --file ~/Documents/MyStuff/UM/Research/PRS_famHX/pheno_files/UKBB_CoronaryArteryDisease_PRS_LDpred_rho0.001_allchr.scores_pheno.txt 
#--stratum_col 34 --pheno_col 21 --sex_col 7 --birthYear_col 8 --age_col 9 --grs_col 3 --output test 
#--quant_covar_cols 9,10,11,12 --quant_covar_labels PC1,PC2,PC3,PC4


###########################################################
###################### Main ##########################
###########################################################

dat<-fread(file)

#to do: censor, young

# stratum is 1 for yes, 0 for no, NA for unknown/missing
# to do: check binary variables
subset<-dat[!is.na(dat[[strat_col]])] #remove if NA for stratum
print(paste("Data dimensions after removing samples with NA stratum:",dim(subset)[1],dim(subset)[2]))

subset<-subset[!is.na(subset[[pheno_col]])] #remove if NA for pheno
print(paste("Data dimensions after removing samples with NA phenotype:", dim(subset)[1],dim(subset)[2]))

qsub<-subset #rename for ease of already written code
  ##censor
##  if (censor==TRUE){
  #  if (young==TRUE){
   #   qsub<-qsub[qsub[[part_age]]<60,]
    #  out<-paste(sep="_",out,"young","censor")
    #} else if (young==FALSE) {
    #  qsub<-qsub[qsub[[part_age]]>=60,]
    #  out<-paste(sep="_",out,"old","censor")
  #  }}
  
  
qsub$old<-qsub[[part_age]]>50 #change age if needed
qsub$age<-2019-qsub[[birthYear]]
qsub$age_std<-scale(qsub$age)
qsub$part_age_sq<-qsub[[part_age]]^2
part_age_sq<-which(names(qsub)=="part_age_sq")
qsub$part_age_std<-scale(qsub[[part_age]]) #standardize so comparable to binary and inverse normalized GRS
qsub$part_age_sq_std<-scale(qsub[[part_age_sq]])
qsub$birthYear_std<-scale(qsub[[birthYear]])
qsub$invNormGRS<-rankNorm(qsub[[grs]]) #inverse normalize GRS
part_age_sq_std<-which(names(qsub)=="part_age_sq_std")
part_age_std<-which(names(qsub)=="part_age_std")
birthYear_std<-which(names(qsub)=="birthYear_std")
age_std<-which(names(qsub)=="age_std")

  
#plot distribution of GRS with sex and FH and age
pdf_fn<-paste(sep=".",out,"GRS_distribution.pdf")
pdf(file=pdf_fn,height=6,width=6,useDingbats=FALSE)
x<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=factor(get(names(qsub)[pheno_col]))))  +  geom_density(alpha=0.5) + theme_bw() + 
  scale_fill_manual(values=c("purple","green3"),name="Trait Status",labels=c("Control","Case")) + labs(x="GRS")
y<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=factor(get(names(qsub)[strat_col])))) + geom_density(alpha=0.5) + theme_bw() + 
  scale_fill_manual(values=c("goldenrod3","darkblue"),name="Family History",labels=c("Negative","Positive")) + labs(x="GRS")
z<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=old)) + geom_density(alpha=0.5) + theme_bw() + 
    scale_fill_manual(values=c("light blue","grey"),name="Age",labels=c("<=50",">50")) + labs(x="GRS")
grid.arrange(x,y,z)
dev.off()
  
#establish color palette early 
pamp<-lacroix_palette("Pamplemousse",type = "discrete", n=6) #6 elements
pamp2<-c(pamp[1:6],"#8B4789") #7 elements 
  
### TO do correlations

###### Test logistic regression
predictor_list<-list(strat_col,
                  grs,
                  c(strat_col,bcovar,qcovar,sex),
                  c(grs,bcovar,qcovar,sex),
                  c(bcovar,qcovar,sex),
                  c(strat_col,grs),
                  c(strat_col,grs,bcovar,qcovar,sex))
print(str(predictor_list))
predictor_list_labels<-list("Family History",
                         "GRS",
                         c(bcovar_lab,qcovar_lab,"sex"),
                         c("FH",bcovar_lab,qcovar_lab,"sex"),
                         c("GRS",bcovar_lab,qcovar_lab,"sex"),
                         c("Family History","GRS"),
                         c("Family History","GRS",bcovar_lab,qcovar_lab,"sex"))
                         
print(str(predictor_list_labels))
model_name<-list("Family History Only",
              "GRS Only",
              "Covariates Only",
              "Family History + Covariates",
              "GRS + Covariates",
              "Additive",
              "Additive + Covariates")

model_objects<-list()
for (i in length(predictor_list)){
  obj<-model(qsub,pheno_col,predictor_list[i],predictor_list_labels[i],out,model_name[i])
  print(str(obj))
}
#setNames(model_objects,model_name)
#print(str(model_objects))

  
 