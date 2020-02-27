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

#https://medium.com/@outside2SDs/an-overview-of-correlation-measures-between-categorical-and-continuous-variables-4c7f85610365

###########################################################
###################### Functions ##########################
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
  ts<-matrix(NA,2,qtile+1) #initialize tiles matrix
  if (qfirst==TRUE){ #calculate quantiles before stratification
    tiles<-quantile(df[[GRS_col]],seq(from=0,to=1,by=p)) #quantile values
    for (i in 1:length(index)-1) {
      for (r in c(1,2)){ #iterate over stratum
        percents<-names(tiles)
        ts[r,]<-tiles
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
      ts[r,]<-tiles
      percents<-names(tiles)
      for (i in 1:length(index)-1) {
        prev_list<-sdf[sdf[[GRS_col]] > tiles[i] & sdf[[GRS_col]] <= tiles[i+1]][[prev_col]]
        prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ns[r,i]<-length(prev_list)
        ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
      } 
    }
  }
  #create object 
  pqs<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=ts,percents=percents)
  class(pqs)<-"prev_quantile_stratum_obj"
  return(pqs)
}

#calcualte odds ratio for top of dist and positive family history
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
model<-function(df,grs_col,strat_col,pheno_col,covar,out,name,FH_int){
  mobj<-list() #initialize object
  auroc<-data.frame() #initialize AUROC list 
  
  ## inverse rank normalize GRS in the population of interest because FHIGRS is inv normal scale
  df$invNormGRS<-rankNorm(df[[grs_col]])
  grs_col<-which(names(df)=="invNormGRS") #new GRS col
  
  #model for covar only
  formula<-as.formula(paste(colnames(df)[pheno_col],"~",
                            paste(colnames(df)[c(covar)],collapse="+"),
                            sep=""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  print(summary(glm.obj))
  cdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar))
  mobj[["covar_only"]]<-cdf
  conly<-glm.obj
  #auroc<-append(auroc,my_ROC_plot(df,formula,out,paste(sep="_",name,"covar_only")))
  r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"covar_only"))
  auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  
  #model for GRS
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  gdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"GRS")) #standard GRS only
  mobj[["GRS"]]<-gdf
  grs<-glm.obj
  #auroc<-append(auroc,my_ROC_plot(df,formula,out,paste(sep="_",name,"GRS")))
  r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"covar_only"))
  auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  
  
  ## model for family history + GRS with interaction term
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar)], collapse = "+"), "+" ,
                            paste(colnames(df)[c(strat_col,grs_col)],collapse="*"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  idf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH","GRS","FH*GRS"))
  mobj[["interaction"]]<-idf
  int<-glm.obj
  #auroc<-append(auroc,my_ROC_plot(df,formula,out,paste(sep="_",name,"interaction")))
  r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"covar_only"))
  auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  
  ## model for additive GRS + family history
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  adf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH","GRS"))
  mobj[["additive"]]<-adf
  add<-glm.obj
  #auroc<-append(auroc,my_ROC_plot(df,formula,out,paste(sep="_",name,"additive")))
  r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"covar_only"))
  auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  
  ## model for family history 
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH"))
  mobj[["family"]]<-fdf
  fh<-glm.obj
  #auroc<-append(auroc,my_ROC_plot(df,formula,out,paste(sep="_",name,"FH")))
  r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"covar_only"))
  auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  
  ##model  for interaction with  family history
  #if the FH_int variable is in the covar list we need to remove it to maintain expected order of the variables in the glm model output 
  if (!is.null(FH_int)){
    covar2<-covar
    if (FH_int %in% covar){
      covar2<-covar2[-(covar2==FH_int)]
    }
    if (strat_col %in% covar){
      covar2<-covar2[-(covar2==FH)]
    }
    formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                              paste(sep="+",paste(colnames(df)[c(covar2,strat_col,FH_int)], collapse = "+"),paste(colnames(df)[c(strat_col,FH_int)],collapse="*")),
                              sep = ""))
    glm.obj<-glm(formula=formula,data=df,family="binomial")
    #colnames(df)[FH_int]<-"FH_interaction_term"
    
    fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar2,"FH","FH_int","FH*FH_int")) #generalized so FH_int can be anything
    mobj[["FH_interaction"]]<-fdf
    fh_int<-glm.obj
    r<-ROC_v2(df,formula,pheno_col,out,paste(sep="_",name,"FH_interaction"))
    auroc<-rbind(auroc,data.frame(cvAUC=r$cvAUC,se=r$se,LB=r$ci[1],UB=r$ci[2]))
  }
  
  ## saturated model
  #a <- factor(1:nrow(df))
  #fit <- glm(df[[pheno_col]]~a,family=binomial)
  ## compare goodness of model fit 
  if (!is.null(FH_int)){
    model_list<-c("conly","grs","int","add","fh","fh_int")
  } else { model_list<-c("conly","grs","int","add","fh")
  }
  fit<- matrix(ncol=8, nrow=length(model_list))
  for (i in 1:length(model_list)){
    fit[i,1]<-get(model_list[i])$deviance
    fit[i,2]<-get(model_list[i])$null.deviance
    fit[i,3]<-get(model_list[i])$aic
    fit[i,4]<-hoslem.test(df[[pheno_col]], fitted(get(model_list[i])))$p.value
    
    #Cox-Snell pseudo R2 and max adjusted R2
    ## intercept only model 
    formula<-as.formula(paste(colnames(df)[pheno_col], "~","1"))
    null.obj<-glm(formula=formula,data=df,family="binomial")
    l0=as.numeric(logLik(null.obj))
    l1=as.numeric(logLik(get(model_list[i])))
    R2_cox<-1-(exp(l0 - l1))^{2/nrow(df)}
    R2_max<-1-exp(l0)^{2/nrow(df)}
    max_adj_R2<-R2_cox/R2_max
    fit[i,5]<-R2_cox
    fit[i,6]<-max_adj_R2
    fit[i,7]<-PseudoR2(get(model_list[i]),which="all")[["Nagelkerke"]]
    fit[i,8]<-BrierScore(get(model_list[i]),scaled=FALSE) #what is the range when unscaled?!
  }
  fit_df<-data.frame(fit)
  names(fit_df)<-c("deviance","null_deviance","aic","hoslem","cox_snell_r2","max_adj_r2","Nagelkerke","Brier")
  fit_df$cvAUC<-auroc$cvAUC
  fit_df$AUC_LB<-auroc$LB
  fit_df$AUC_UB<-auroc$UB
  fit_df$dev_test<-1-pchisq(fit_df$deviance,get(model_list[i])$df.residual)
  if (!is.null(FH_int)){
    fit_df$model<-c("covar_only","GRS","GRS*FH","GRS+FH","FH","FH_int")
  } else {
    fit_df$model<-c("covar_only","GRS","GRS*FH","GRS+FH","FH")
  }
  
  file_n<-paste(sep=".",out,name,"modelGOF.txt")
  write.table(format(fit_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
  
  ## anova comparing to additive model
  anova_df<-rbind(
    data.frame(dev=anova(add,grs,test="LRT")$Deviance[2],pval=anova(add,grs,test="LRT")$`Pr(>Chi)`[2],method="GRS"),
    data.frame(dev=anova(add,fh,test="LRT")$Deviance[2],pval=anova(add,fh,test="LRT")$`Pr(>Chi)`[2],method="FH"))
  file_n<-paste(sep=".",out,name,"compareToAdditive.txt")
  write.table(format(anova_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
  
  class(mobj)<-"model_obj"
  return(mobj)
  
}

#indicator variable used for outcome Y/N in top X-tile
model_indicator<-function(df,value,grs_col,strat_col,pheno_col,covar,qfirst=FALSE){
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
  #gdf$dev<-summary(glm.obj)$null.deviance-summary(glm.obj)$deviance
  gdf$nagelkerke<-PseudoR2(glm.obj,which="all")[["Nagelkerke"]]
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
    glm.obj<-glm(formula=formula,data=df[df[[strat_col]]==(s-1),],family="binomial") #subsets data frame to stratum
    df_name<-paste0("stratum",s-1)
    m<-matrix(table(df[df[[strat_col]]==(s-1)]$dist,df[df[[strat_col]]==(s-1)][[pheno_col]])/nrow(df[df[[strat_col]]==(s-1)]),byrow=FALSE,nrow=2)
    top<-m[2,2]/sum(m[2,])  #prevalence in top of distribution
    bottom<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
    assign(df_name,cbind(data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")), 
                         "top_prev"=top,"bottom_prev"=bottom,nagelkerke=PseudoR2(glm.obj,which="all")[["Nagelkerke"]]))
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
  cdf$nagelkerke<-PseudoR2(glm.obj,which="all")[["Nagelkerke"]]
  mobj[["conditional"]]<-cdf
  
  
  ##model using I(FH) NOT REALLY A CUT POINT JUST ZERO or ONE
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  ydf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")) 
  m<-matrix(table(df[[strat_col]],df[[pheno_col]]),byrow=FALSE,nrow=2)
  ydf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  ydf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  ydf$nagelkerke<-PseudoR2(glm.obj,which="all")[["Nagelkerke"]]
  mobj[["family"]]<-ydf  
  
  class(mobj)<-"model_obj"
  return(mobj)
}


######### ROC 
my_ROC_plot<-function(qsub,formula,out,name){
  
  k<-5 #number of folds
  
  # set seed and shuffle again so we get the same order as before, this time we are  going to plot all together
  #Randomly shuffle the data after setting seed
  set.seed(1234)
  shufData<-qsub[sample(nrow(qsub)),]
  folds <- cut(seq(1,nrow(qsub)),breaks=k,labels=FALSE)
  
  #pick k colors for plots
  #you can pick different palette here: https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
  colors<-brewer.pal(k, "Dark2")	
  
  #initialize list for saving data
  roc_list<-list()				 
  
  #Perform k fold cross validation and plot all together
  pdf(file=paste(sep="_",out,name,"ROC_all_folds.pdf"),height=5,width=6,useDingbats = )
  for(i in 1:k){
    #Segment your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- shufData[testIndexes, ]
    trainData <- shufData[-testIndexes, ]
    #fit model on training data
    glm.obj<-glm(formula=formula,data=trainData,family="binomial")
    
    #use on testing data
    testPreds<-predict(glm.obj,testData)
    if (i==1){				  
      roc.obj<-roc(testData[[pheno]], testPreds, plot=TRUE,ci=TRUE,col=colors[i],print.auc=TRUE,print.auc.col=colors[i],print.auc.y=0.5-(i*0.05))
    } else {
      roc.obj<-roc(testData[[pheno]], testPreds, plot=TRUE,ci=TRUE,add=TRUE,col=colors[i],print.auc=TRUE,print.auc.col=colors[i],print.auc.y=0.5-(i*0.05))
    }
    roc_list[[i]]<-roc.obj	#save ROC  object			
  }
  avgAUC<-mean(unlist(lapply(roc_list,auc)))
  text(0.9,1,paste(sep=" ","mean AUC:",format(avgAUC,digits=3)))
  dev.off()
  #parse the object
  avgAUC<-mean(unlist(lapply(roc_list,auc)))
  return(avgAUC)
}

##### ROC plot function that calculates 95% CI
ROC_v2<-function(qsub,formula,pheno_col,out,name,V=5){
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




###########################################################
###################### Main ##########################
###########################################################



############# T2D in UKBB
t2d<-fread("pheno_files/UKBB_Type2Diabetes_PRS_LDpred_rho0.01_allchr.scores_pheno.txt")
batch<-6 #use genotyping array instead of actual batch, too many batches
part_age<-91
fh<-69
birthYear<-8
pcs<-c(9,10,11,12)
sex<-7
pheno<-26
grs<-3
out<-"UKBB.T2D"

############# CAD in UKBB
cad<-fread("pheno_files/UKBB_CoronaryArteryDisease_PRS_LDpred_rho0.001_allchr.scores_pheno.txt")
batch<-6 #use genotyping array instead of actual batch, too many batches
part_age<-91
fh<-34
birthYear<-8
pcs<-c(9,10,11,12)
sex<-7
pheno<-21
grs<-3
out<-"UKBB.CAD"

##variables
cutpt<-0.9
quantiles<-c(20) ##  To do: test with more 
censor<-FALSE
young<-FALSE
#dat<-t2d
dat<-cad
out<-out
strat_col<-fh
pheno_col<-pheno
grs_col<-grs
dig<-3

#fh is 1 for yes, 0 for no, NA for missing/unknown
# fh is 1 for yes, 0 for no, NA for unknown/missing

#########################################################


subset<-dat[!is.na(dat[[strat_col]])] #remove if NA for stratum
print(paste("Data dimensions after removing samples with NA stratum:",dim(subset)[1],dim(subset)[2]))

subset<-subset[!is.na(subset[[pheno_col]])] #remove if NA for pheno
print(paste("Data dimensions after removing samples with NA phenotype:", dim(subset)[1],dim(subset)[2]))

##stratify then calculate quantiles
sobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs,prev_col=pheno,strat_col=fh,qfirst=FALSE)
##calculate quantiles then stratify
qobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=TRUE)

size<-length(quantiles) #number of quantiles being tested, size of obj
for (i in 1:size){ #across q-quantiles
  for (j in c(1,2)){ #across stratum
    list_length<-length(sobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-sobj[[i]]$prev[j,][-list_length]
    se<-sobj[[i]]$se[j,][-list_length]
    n<-sobj[[i]]$n[j,][-list_length]
    tiles<-sobj[[i]]$tiles[j,][-1]
    lower_tile<-sobj[[i]]$tiles[j,][1:list_length-1]
    upper_tile<-sobj[[i]]$tiles[j,][2:list_length]
    percents<-sobj[[i]]$percents[-1]
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(j-1,list_length-1)
    if (j==1) {
      sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
    } else {
      sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
    }
  }
  
  for (j in c(1,2)){ #across stratum
    list_length<-length(qobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-qobj[[i]]$prev[j,][-list_length]
    se<-qobj[[i]]$se[j,][-list_length]
    n<-qobj[[i]]$n[j,][-list_length]
    tiles<-qobj[[i]]$tiles[j,][-1]
    lower_tile<-qobj[[i]]$tiles[j,][1:list_length-1]
    upper_tile<-qobj[[i]]$tiles[j,][2:list_length]
    percents<-qobj[[i]]$percents[-1]
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(j-1,list_length-1)
    if (i==1 & j==1) {
      qdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
    } else {
      qdf<-rbind(qdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
    }
  }
  
  qsub<-subset
  
  ##censor
  if (censor==TRUE){
    if (young==TRUE){
      qsub<-qsub[qsub[[part_age]]<60,]
      out<-paste(sep="_",out,"young","censor")
    } else if (young==FALSE) {
      qsub<-qsub[qsub[[part_age]]>=60,]
      out<-paste(sep="_",out,"old","censor")
    }}
  
  
  qsub$old<-qsub[[part_age]]>50 #change age if needed
  qsub$age<-2019-qsub[[birthYear]]
  qsub$age_std<-scale(qsub$age)
  qsub$part_age_sq<-qsub[[part_age]]^2
  part_age_sq<-which(names(qsub)=="part_age_sq")
  qsub$part_age_std<-scale(qsub[[part_age]]) #standardize so comparable to binary and inverse normalized GRS
  qsub$part_age_sq_std<-scale(qsub[[part_age_sq]])
  qsub$birthYear_std<-scale(qsub[[birthYear]])
  part_age_sq_std<-which(names(qsub)=="part_age_sq_std")
  part_age_std<-which(names(qsub)=="part_age_std")
  birthYear_std<-which(names(qsub)=="birthYear_std")
  age_std<-which(names(qsub)=="age_std")
  
  #plot distribution of GRS with sex and FH and age
  pdf_fn<-paste(sep=".",out,quantiles[i],"GRS_distribution.pdf")
  pdf(file=pdf_fn,height=6,width=6,useDingbats=FALSE)
  x<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=factor(get(names(qsub)[pheno_col]))))  +  geom_density(alpha=0.5) + theme_bw() + 
    scale_fill_manual(values=c("purple","green3"),name="Trait Status",labels=c("Control","Case")) + labs(x="GRS")
  y<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=factor(get(names(qsub)[fh])))) + geom_density(alpha=0.5) + theme_bw() + 
    scale_fill_manual(values=c("goldenrod3","darkblue"),name="Family History",labels=c("Negative","Positive")) + labs(x="GRS")
  z<-ggplot(qsub,aes(x=get(names(qsub)[grs]),fill=old)) + geom_density(alpha=0.5) + theme_bw() + 
    scale_fill_manual(values=c("light blue","grey"),name="Age",labels=c("<=50",">50")) + labs(x="GRS")
  grid.arrange(x,y,z)
  dev.off()
  
  #establish color palette early 
  pamp<-lacroix_palette("Pamplemousse",type = "discrete", n=6) #6 elements
  pamp2<-c(pamp[1:6],"#8B4789") #7 elements 
  
  #correlation matrix between covariates, need to handle the covarites given for model and birthyear/sex
  var_list<-c(pheno,batch,birthYear,sex,pcs,strat_col,grs,part_age)
  var_names<-c("Pheno","Batch","birthYear","Sex","PC1","PC2","PC3","PC4","FamHx","GRS","EnrollmentAge")
  mydata.rcorr<-subset(qsub,select=var_list) %>% as.matrix() %>% rcorr(type="pearson") ### Pearson ok?
  row.names(mydata.rcorr$r)<-var_names
  colnames(mydata.rcorr$r)<-var_names
  pdf_fn<-paste(sep=".",out,quantiles[i],"correlation.pdf")
  col1<-colorRampPalette(c(pamp[[1]],"#FFFFFF",pamp[[6]]))
  pdf(file=pdf_fn,height=8,width=8,useDingbats=FALSE)
  par(mar=c(1.25,1.25,1.25,1.25))
  print(corrplot.mixed(mydata.rcorr$r,upper="ellipse",lower.col="black",upper.col=col1(100),tl.pos="lt",tl.col="black"))
  dev.off()
  
  #### correlation between age of participation and FH
  #visual
  ggplot(qsub,aes(y=get(names(qsub)[part_age]))) + geom_boxplot() + theme_bw() + labs(y="Age at time of self report Family History")
  ggplot(qsub,aes(x=get(names(qsub)[part_age]),y=get(names(qsub)[strat_col]))) + geom_point(alpha=0.2) + theme_bw()
  #logreg is enrollment age predictive  of family history
  formula<-as.formula(paste(colnames(qsub)[strat_col], "~",
                            paste(colnames(qsub)[part_age_std], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=qsub,family="binomial")
  print(summary(glm.obj))
  #point biserial correlation
  biserial.cor(qsub[[part_age]],qsub[[fh]]) 
  biserial.cor(qsub[[birthYear]],qsub[[fh]]) 
  biserial.cor(qsub[[grs]],qsub[[fh]]) 
  
  #### leave one variable out at a time and test participation age instead of birthyear, standard model without cutpoint 
  tcdf<-data.frame() #continuous
  idf<-data.frame() #indicator
  edf<-data.frame() #continuous extended
  test_covar<-function(qsub,covar,out,name,int_term){
    #change these if other covars are included
    covar_indices<-as.character(c(birthYear_std, batch, pcs,sex, part_age_std,part_age_sq_std)) #standardized age variables 
    covar_names<-c("birthYear","batch","PC1", "PC2", "PC3", "PC4","sex","participation_age","participation_age_squared")
    
    #contious variable as predictor
    if (length(covar)>0){
      mobj<-model(df=qsub,grs_col=grs,pheno_col=pheno,strat_col=fh,covar=covar,out,name,int_term) #last variable is the column number for the variable we want to test interaction wth family history
      if (!is.null(int_term)){
        model_list<-c("covar_only","GRS","interaction","additive","family","FH_interaction") #model being tested 
        pred_list<-c("GRS","FH","FH*GRS","Int")
        model_list_reorder<-c("GRS", "family", "additive","FH_interaction")
      } else {
        model_list<-c("covar_only","GRS","interaction","additive","family")
        pred_list<-c("GRS","FH","FH*GRS")
        model_list_reorder<-c("GRS", "family", "additive")
      }
      
      for (l in 1:length(mobj)){ #loop over list from model function
        mobj[[names(mobj[l])]]$model<-names(mobj)[l]
        mobj[[l]]$OR<-exp(mobj[[l]][,'Estimate'])
        mobj[[l]]$LB<-exp(mobj[[l]][,'Estimate']-(1.96*mobj[[l]][,'Std..Error']))
        mobj[[l]]$UB<-exp(mobj[[l]][,'Estimate']+(1.96*mobj[[l]][,'Std..Error']))
        mobj[[l]]$model<-model_list[l]
        mobj[[l]]$pred<-row.names(mobj[[l]])
      }
      d<-do.call("rbind",mobj)
      d$pred<-mapvalues(d$pred,from=covar_indices,to=covar_names) #name covariate indices, throws error but should be ok
      d$qtile=quantiles[i]
      dsub<-d[d$pred %in% pred_list,] #we can ignore pcs, batch, sex, etc.
      dsub$model<-as.factor(dsub$model)
      dsub$model<-factor(dsub$model,levels=model_list_reorder) #reorder levels
      dsub$model<-revalue(dsub$model,c("family"="FamilyHistory"))
    } else {
      dsub<-NA
      d<-NA
    }
    #indicator for high or low in distirbution used as predictor
    cutpts<-c(0.8,0.9,0.95,0.99)
    iobj<-lapply(cutpts,model_indicator,df=qsub,covar=covar,pheno_col=pheno,grs_col=grs,strat_col=fh,qfirst=FALSE)
    for (l in 1:length(iobj)){ #loop over list from model function, one per cutpt value
      name_list<-names(iobj[[l]])
      for (n in 1:length(name_list)){ #across the models tested 
        colnames<-lapply(iobj[[l]][name_list[n]],names)[[1]] #pull column names from data frame
        z<-data.frame(iobj[[l]][name_list[n]])['dist',] #pull out summ stats for variable of interest
        names(z)<-colnames
        OR=exp(z[['Estimate']])
        LB=exp(z[['Estimate']]-1.96*z[['Std..Error']])
        UB=exp(z[['Estimate']]+1.96*z[['Std..Error']])
        pval=z[['Pr...z..']]
        if (n==1 & i==1 & l==1) { #inititalize df for first q-quantile, first cutpoint, first quantile logic
          model_df<-cbind(z,qtile=quantiles[i],cutpt=cutpts[l],qfirst="stratifyFirst",name=name_list[n],OR=OR,LB=LB,UB=UB,pval=pval)
        } else {
          model_df<-rbind(model_df,cbind(z,qtile=quantiles[i],cutpt=cutpts[l],qfirst="stratifyFirst",name=name_list[n],OR=OR,LB=LB,UB=UB,pval=pval))
        }
      }}
    obj<-list() #initialize object
    obj[["continuous"]]<-dsub
    obj[["continuous_extended"]]<-d
    obj[["indicator"]]<-model_df
    return(obj)
  }
  
  ##reset birthyear_std to be the column for age_std so that it is on the same direction as participaation (baseline) age
  birthYear_std<-age_std
  
  covar<-c(birthYear_std,batch,pcs,sex) #birthyear + batch + pcs + sex + [FH, GRS, FH+GRS, FH*GRS, FHiGRS, null]
  name<-"birthYear_batch_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std) #testing interaction between participation age and FH
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(part_age_std,batch,pcs,sex) #use participation age instead of birthyear, highly negative correlated 
  name<-"partAge_batch_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(part_age_std,part_age_sq_std,batch,pcs,sex) #use participation age and age squared variable
  name<-"partAge_partAgeSq_batch_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(pcs,sex,batch) #pcs + sex + batch (no age or birthyear avariable)
  name<-"batch_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(birthYear_std,pcs,sex) #birthyear + pcs + sex (no batch)
  name<-"birthYear_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(part_age_std,pcs,sex) #participaton age + pcs + sex (no batch)
  name<-"partAge_PCs_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(birthYear_std,sex,batch) #birthyear + sex + batch (no pcs)
  name<-"batch_birthYear_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(part_age_std,sex,batch) #participation age + sex + batch (no pcs)
  name<-"batch_partAge_sex"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c(part_age_std,birthYear_std,pcs,sex,batch) #participation age and birthyear with other covarites
  name<-"batch_partAge_sex_PCs_birthYear"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  covar<-c()
  name<-"none"
  tmp<-test_covar(qsub,covar,out,name,part_age_std)
  tmp$indicator$covar<-name
  idf<-rbind(idf,tmp$indicator)
  
  pdf(file=paste0(out,"_LOCO_model_continuous.pdf"),height=8,width=12,useDingbats = FALSE)
  print(ggplot(tcdf,aes(x=OR,y=pred,color=pred)) + geom_point() + theme_bw() + facet_wrap(~covar+model) + 
          geom_errorbarh(aes(xmin=tcdf$LB,xmax=tcdf$UB)))
  dev.off()
  
  edf$model<-as.factor(edf$model)
  edf$model<-factor(edf$model,levels(edf$model)[c(2,3,6,5,1,7,4)]) #covar, family, GRS, FHIGRS, additive, interaction, FH_interaction
  names(edf)[4]<-"pvalue"
  edf<-edf[edf$pred!="Int",]
  edf[edf$pvalue==0,]$pvalue<-1e-324
  edf$neglog10<- -log10(edf$pvalue)
  pdf(file=paste0(out,"_all_model_continuous.pdf"),height=20,width=20,useDingbats = FALSE)
  print(ggplot(edf,aes(x=OR,y=pred,color=pred)) + geom_point() + theme_bw() + facet_wrap(~covar+model,ncol=6) + geom_vline(xintercept=1,linetype="dashed",color="red") +
          geom_errorbarh(aes(xmin=edf$LB,xmax=edf$UB)) + coord_cartesian(xlim=c(0,2.5)))
  dev.off()
  
  #pvalues
  pdf(file=paste0(out,"_all_model_continuous_tile.pdf"),height=6,width=12,useDingbats = FALSE)
  print(ggplot(edf,aes(x=model,y=pred,fill=neglog10)) + geom_tile() + theme_bw() + scale_fill_gradient(high="red",low="grey") + facet_wrap(~covar))
  dev.off()
  
  #goodness of fit
  gof_df<-data.frame()
  files<-Sys.glob(paste0(out,"*.modelGOF.txt")) #make sure there aren't extraneous modelGOF.txt files 
  for (fi in 1:length(files)){
    file_df<-read.table(file=files[fi],header=T)
    file_df$covar<-gsub(paste0(out,"."),"",gsub(".modelGOF.txt","",files[fi]))
    print(head(file_df))
    gof_df<-rbind(gof_df,file_df)
    #file.remove(files[fi]) #delete file so it doesn't get pulled back up later
  }
  gof_df$model<-revalue(gof_df$model,c("FH"='family',"GRS+FH"="additive","GRS*FH"="FH interaction with GRS","FH_int"="FH interaction with enrollment age"))
  gof_df$model<-factor(gof_df$model,levels(gof_df$model)[c(1,2,4,7,6,3)])
  write.table(gof_df,file=paste0(out,".all.model.GOF.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  
  pdf(file=paste0(out,"_all_model_continuous_gof.pdf"),height=6,width=8,useDingbats = FALSE)
  print(ggplot(gof_df,aes(x=Brier,y=Nagelkerke,color=model)) + geom_point() + theme_bw() + facet_wrap(~covar))
  dev.off()
  
  pdf(file=paste0(out,"_all_model_continuous_gof_OR.pdf"),height=6,width=20,useDingbats = FALSE)
  x<-ggplot(gof_df,aes(x=model,y=Nagelkerke)) + geom_point() + theme_bw() + facet_wrap(~covar,nrow=1) + theme(axis.text.x=element_text(angle=45))
  y<-ggplot(edf,aes(x=model,y=OR,color=pred)) + geom_point() + theme_bw() + facet_wrap(~covar,nrow=1) + geom_vline(xintercept=1,linetype="dashed",color="red") +
    geom_errorbar(aes(ymin=edf$LB,ymax=edf$UB)) + coord_cartesian(ylim=c(0,3)) +  theme(axis.text.x=element_text(angle=45),legend.position="bottom")
  grid.arrange(x,y,nrow=2)  
  dev.off()
  
  #remove predictors that are just controlling for confounders 
  edf_sub<-edf[grep("PC",edf$pred,invert=TRUE),]
  edf_sub<-edf_sub[edf_sub$pred!="batch" & edf_sub$pred!="sex",]
  #Nagelkerke's R
  pdf(file=paste0(out,"_all_model_continuous_gof_OR_simple.pdf"),height=6,width=20,useDingbats = FALSE)
  x<-ggplot(gof_df,aes(x=model,y=Nagelkerke)) + geom_point() + theme_bw() + facet_wrap(~covar,nrow=1) + theme(axis.text.x=element_text(angle=45))
  y<-ggplot(edf_sub,aes(x=model,y=OR,color=pred)) + geom_point(alpha=0.5) + theme_bw() + facet_wrap(~covar,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="red",alpha=0.25) +
    geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB)) + coord_cartesian(ylim=c(0,3)) +  theme(axis.text.x=element_text(angle=45),legend.position="bottom") + scale_color_brewer(palette = "Set1")
  grid.arrange(x,y,nrow=2)  
  dev.off()
  
  #AUROC 
  pdf(file=paste0(out,"_all_model_continuous_auroc_OR_simple.pdf"),height=6,width=20,useDingbats = FALSE)
  x<-ggplot(gof_df,aes(x=model,y=cvAUC)) + geom_point() + theme_bw() + facet_wrap(~covar,nrow=1) + theme(axis.text.x=element_text(angle=45))
  y<-ggplot(edf_sub,aes(x=model,y=OR,color=pred)) + geom_point(alpha=0.5) + theme_bw() + facet_wrap(~covar,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="red",alpha=0.25) +
    geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB)) + coord_cartesian(ylim=c(0,3)) +  theme(axis.text.x=element_text(angle=45),legend.position="bottom") + scale_color_brewer(palette = "Set1")
  grid.arrange(x,y,nrow=2)  
  dev.off()
  
  #### clean plot for poster
  #look at fewer model
  gof_df_sub<-gof_df[gof_df$covar %in% c("batch_PCs","birthYear_batch_PCs","partAge_batch_PCs", "batch_partAge_sex_PCs_birthYear"),]
  edf_sub2$covar<-factor(edf_sub2$covar,levels(as.factor(edf_sub2$covar))[c(2,3,4,1)])
  gof_df_sub$covar<-factor(gof_df_sub$covar,levels(as.factor(gof_df_sub$covar))[c(2,3,4,1)])
  edf_sub2$covar<-revalue(edf_sub2$covar,c("birthYear_batch_PCs"="age + batch + \nPCs + sex","partAge_batch_PCs"="enrollment age + batch + \nPCS ", "batch_PCs"="batch + PCS", 
                                           "batch_partAge_PCs_birthYear"="enrollment age + age + \nbatch + sex + PCs"))
  gof_df_sub$covar<-revalue(gof_df_sub$covar,c("birthYear_batch_PCs"="age + batch + \nPCs","partAge_batch_PCs"="enrollment age + batch + \nPCS ", "batch_PCs"="batch + PCS ", 
                                               "batch_partAge_PCs_birthYear"="enrollment age + age + \nbatch  + PCs"))
  
  gof_df_sub$model<-revalue(gof_df_sub$model,c("family"="family history","covar_only"="covariates only"))
  edf_sub2$model<-revalue(edf_sub2$model,c("family"="family history","covar_only"="covariates only"))
  edf_sub2$model<-revalue(edf_sub2$model,c("interaction"="FH interaction with GRS","FH_interaction"="FH interaction with enrollment age"))
  
  edf_sub2$pred<-revalue(edf_sub2$pred,c("birthYear"="age","FH"="family history","participation_age"="enrollment age","FH*GRS"="FHxGRS interaction term","FH_int"="enrollment age", "FH*FH_int"="FH x enrollment age interaction term"))
  edf_sub2$pred<-factor(edf_sub2$pred,levels(as.factor(edf_sub2$pred))[c(1,3,2,7,5,6,4)])
  pdf(file=paste0(out,"_all_model_continuous_auroc_OR_simple_poster.pdf"),height=5.5,width=9,useDingbats = FALSE)
  x<-ggplot(gof_df_sub,aes(x=model,y=cvAUC)) + geom_point(size=2,alpha=0.5) + theme_bw() + facet_wrap(~covar,nrow=1) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.background =element_rect(fill="white")) +
    labs(x="",y="5-fold CV AUROC") + scale_y_continuous(labels=function(x) sprintf("%.1f", x)) + geom_errorbar(aes(ymin=gof_df_sub$AUC_LB,ymax=gof_df_sub$AUC_UB),width=0.4) +
    geom_hline(linetype="dashed",color="red",yintercept=0.5)
  y<-ggplot(edf_sub2,aes(x=model,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~covar,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
    geom_errorbar(aes(ymin=edf_sub2$LB,ymax=edf_sub2$UB),width=0.4) +  theme(strip.background =element_rect(fill="white"),axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank()) +  
    scale_color_manual(values=pamp2) + labs(y="Odds Ratio",x="Predictors") + scale_y_continuous(labels=function(x) sprintf("%.1f", x)) +
    guides(col = guide_legend(nrow = 1))
  grid.arrange(x,y,sub=textGrob(paste0("N=",nrow(qsub))),nrow=3,heights=c(0.36,0.59,0.05)) 
  dev.off()
  
  ##Bhramar plot for poster
  gof_df_roc<-unique(gof_df[c("cvAUC","AUC_LB","AUC_UB","model","covar")])
  #compare model with and without baseline age
  roc<-gof_df_roc[gof_df_roc$model=="covar_only"&gof_df_roc$covar=="batch_PCs",]
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="family"&gof_df_roc$covar=="batch_PCs",])
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="family"&gof_df_roc$covar=="partAge_batch_PCs",])
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="GRS"&gof_df_roc$covar=="batch_PCs",])
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="GRS"&gof_df_roc$covar=="partAge_batch_PCs",])
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="additive"&gof_df_roc$covar=="batch_PCs",])
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="additive"&gof_df_roc$covar=="partAge_batch_PCs",])
  roc$labels<-c("covariates only","covariates + family history","covariates + family history\n + enrollment age","covariates + GRS","covariates + GRS\n + enrollment age", "covariates + family history + GRS","covariates + family history + GRS\n + enrollment age")
  roc$labels<-factor(roc$labels,levels(as.factor(roc$labels))[c(7,1,2,5,6,3,4)])
  roc$group<-as.character(c(0,1,1,2,2,3,3))
  col<-c("black",rep(pamp[2],2),rep(pamp[4],2),rep(pamp[5],2))
  names(col)<-as.character(roc$group)
  pdf(file=paste0(out,"_auroc_model_selection.pdf"),height=3.5,width=4.5,useDingbats = FALSE)
  ggplot(roc,aes(x=labels,y=cvAUC,color=group)) + geom_point(size=2,alpha=0.5) + theme_bw() + labs(x="Predictors",y="5-fold CV AUROC")+
    geom_errorbar(aes(ymin=roc$AUC_LB,ymax=roc$AUC_UB),width=0.3) + theme(axis.text.x=element_text(size=10,angle=45,hjust=1),legend.position="none") +
    scale_color_manual(values=col) + geom_hline(linetype="dashed",color="red",yintercept=0.5)
  dev.off()
  #compare interactions with additive
  roc<-gof_df_roc[gof_df_roc$model=="family" & gof_df_roc$covar=="partAge_batch_PCs",] #covar + partAge + FH
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="FH interaction with enrollment age"&gof_df_roc$covar=="batch_PCs",]) #covar + partAge + FH + partAge*FH
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="additive" & gof_df_roc$covar=="batch_PCs",]) #covar + FH + GRS
  roc<-rbind(roc,gof_df_roc[gof_df_roc$model=="FH interaction with GRS" & gof_df_roc$covar=="batch_PCs",]) #covar + FH + GRS + FH*GRS
  roc$labels<-c("covariates + family history\n + enrollmentAge","covariates + family history + enrollmentAge\n + FH*enrollmentAge","covariates + family history\n + GRS","covariates + family history + GRS\n + FH*GRS")
  roc$group<-as.character(c(0,0,1,1))
  col<-c(rep(pamp2[7],2),rep(pamp[5],2))
  names(col)<-as.character(roc$group)
  pdf(file=paste0(out,"_auroc_model_selection_interactions.pdf"),height=4.5,width=5,useDingbats = FALSE)
  ggplot(roc,aes(x=labels,y=cvAUC,color=group)) + geom_point(size=2,alpha=0.5) + theme_bw() + labs(x="Predictors",y="5-fold CV AUROC")+
    geom_errorbar(aes(ymin=roc$AUC_LB,ymax=roc$AUC_UB),width=0.3) + theme(axis.text.x=element_text(size=10,angle=45,hjust=1),legend.position="none") +
    scale_color_manual(values=col) + geom_hline(linetype="dashed",color="red",yintercept=0.5)
  dev.off()
  
  
  ### make smaller plot with sele cted points
  
  #look at HUNT 1,2,3 as a covariate to explain participation age
  #titrate age for old vs young
  #look at other clinical covariates
  #look at T2D
  #calculate in UKBB
  ### To Do: look at this in T2D, birthyeaar and  participation age  relaly different? age instead of birth  year?
  sub<-idf[idf$name=="FHIGRS"|idf$name=="GRS",]
  pdf(file=paste0(out,"_LOCO_model_indicator.pdf"),height=8,width=8,useDingbats = FALSE)
  print(ggplot(sub,aes(x=OR,y=as.factor(cutpt),color=name)) + geom_point() + theme_bw() + facet_wrap(~covar) + 
          geom_errorbarh(aes(xmin=sub$LB,xmax=sub$UB),height=0.2) + ylab("Threshold for High Risk Group") + xlab("Odds Ratio"))
  dev.off()
  
  ###make plots for the talk from idf
  pdf(file=paste0(out,"_model_progression.pdf"),height=8,width=12,useDingbats = FALSE)
  x<-ggplot(sub[sub$covar=="none",],aes(x=OR,y=as.factor(cutpt),color=name)) + geom_point() + theme_bw() + 
    geom_errorbarh(aes(xmin=sub[sub$covar=="none",]$LB,xmax=sub[sub$covar=="none",]$UB)) + ylab("Threshold for High Risk Group") + xlab("Odds Ratio")
  y<-ggplot(sub[sub$covar=="birthYear_batch_PCs",],aes(x=OR,y=as.factor(cutpt),color=name)) + geom_point() + theme_bw() +
    geom_errorbarh(aes(xmin=sub[sub$covar=="birthYear_batch_PCs",]$LB,xmax=sub[sub$covar=="birthYear_batch_PCs",]$UB)) + ylab("Threshold for High Risk Group") + xlab("Odds Ratio")
  z<-ggplot(sub[sub$covar=="batch_PCs",],aes(x=OR,y=as.factor(cutpt),color=name)) + geom_point() + theme_bw() + 
    geom_errorbarh(aes(xmin=sub[sub$covar=="batch_PCs",]$LB,xmax=sub[sub$covar=="batch_PCs",]$UB)) + ylab("Threshold for High Risk Group") + xlab("Odds Ratio") 
  grid.arrange(x,y,z,ncol=3)
  dev.off()
  
  sub2<-sub[sub$covar=="none" | sub$covar=="birthYear_batch_PCs" | sub$covar=="batch_PCs",]
  sub2$covar<-factor(sub2$covar, levels = c("none","birthYear_batch_PCs","batch_PCs"))
  pdf(file=paste0(out,"_model_progression_facet.pdf"),height=4,width=8,useDingbats = FALSE)
  ggplot(sub2,aes(x=OR,y=as.factor(cutpt),color=name)) + geom_point() + theme_bw() + 
    geom_errorbarh(aes(xmin=sub2$LB,xmax=sub2$UB),height=0.2) + ylab("Threshold for High Risk Group") + xlab("Odds Ratio") +
    facet_wrap(~covar) + theme(legend.position="left",axis.text.x=element_text(size=10,angle=30,hjust=0.5),axis.title=element_text(size=12),legend.text=element_text(size=12)) + 
    scale_color_manual(values=lacroix_palette("Pamplemousse",type = "discrete", n=2)) + geom_vline(xintercept=1,linetype="dashed",color="black")
  dev.off()
  
  #get the colors to match the previous plots
  sub2$covar<-revalue(sub2$covar,c("birthYear_batch_PCs"="MI ~ I(top X-tile) + batch + \nPCs  + birth year", 
                                   "batch_PCs"="MI ~ I(top X-tile) + \nbatch + PCS", "none"="MI ~ I(top X-tile)"))
  pdf(file=paste0(out,"_model_progression_facet_v2.pdf"),height=3,width=7,useDingbats = FALSE)
  ggplot(sub2,aes(y=OR,x=as.factor(cutpt),color=name)) + geom_point() + theme_bw() + 
    geom_errorbar(aes(ymin=sub2$LB,ymax=sub2$UB),width=0.2) + xlab("Threshold for High Risk Group") + ylab("Odds Ratio") +
    facet_wrap(~covar) + theme(legend.position="left",axis.text.x=element_text(size=10,angle=30,hjust=0.5),axis.title=element_text(size=12),legend.text=element_text(size=12)) + 
    scale_color_manual(values=c(pamp[4],pamp[5]),name="Score") + geom_hline(yintercept=1,linetype="dashed",color="black") +
    scale_y_continuous(labels=function(x) sprintf("%.1f", x)) +  theme(strip.background =element_rect(fill="white"),legend.position="bottom")
  dev.off()
  
  
  #####Prevalence across age bins
  #using subset and part_age not  standaardized
  sdf<-data.frame()
  #subset$bin<-cut(subset[[part_age]],3)
  #subset$bin<-cut_number(subset[[part_age]],5)
  #subset$bin<-cut(subset[[part_age]],breaks=c(20,40,60,80,120))
  subset$bin<-cut(subset[[part_age]],breaks=c(35,40,45,50,55,60,65,70))
  age_bins<-unique(subset$bin)
  prop<-c()
  count<-c()
  for (b in 1:length(age_bins)) {
    print(age_bins[b])
    bin_sub<-subset[subset$bin==age_bins[b],]
    prop[b]<-sum(bin_sub[[fh]])/nrow(bin_sub) #how many within each age group have a self reported relative with disease
    count[b]<-nrow(bin_sub)
    print(odds_ratio(bin_sub,grs_col=grs,prev_col=pheno,strat_col=fh))
    sobj<-lapply(quantiles,prev_per_quantile_stratum,df=bin_sub,GRS_col=grs,prev_col=pheno,strat_col=fh,qfirst=FALSE)
    i<-1
    print(str(sobj))
    for (j in c(1,2)){ #across stratum
      list_length<-length(sobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
      prev<-sobj[[i]]$prev[j,][-list_length]
      se<-sobj[[i]]$se[j,][-list_length]
      n<-sobj[[i]]$n[j,][-list_length]
      tiles<-sobj[[i]]$tiles[j,][-1]
      lower_tile<-sobj[[i]]$tiles[j,][1:list_length-1]
      upper_tile<-sobj[[i]]$tiles[j,][2:list_length]
      percents<-sobj[[i]]$percents[-1]
      bins<-rep(quantiles[i],list_length-1)
      strat<-rep(j-1,list_length-1)
      if (j==1 & b==1) {
        sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL,bin=age_bins[b],propFH=prop[b],count=count[b])
      } else {
        sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL,bin=age_bins[b],propFH=prop[b],count=count[b]))
      }
    }
  }
  sdf$frac<-as.numeric(sub("%","",sdf$percents)) #convert factor percentages to numeric
  sdf<-sdf[sdf$frac!=1.00,]
  sdf$ub<-sdf$prev+(1.96*sdf$se)
  sdf$lb<-sdf$prev-(1.96*sdf$se)
  ggplot(sdf,aes(x=frac,y=prev)) + facet_wrap(~bin) + geom_point()
  if (unique(sdf$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=sdf$frac}
  ymax=1
  sdf$label<-paste(sep="\n",paste0("Enrollment age ",sdf$bin),paste0(format(sdf$propFH*100,digits=3),sep="% with family history"),paste0("N=",sdf$count))
  pdf(file=paste(sep="_",out,unique(sdf$q),"prev_by_part_age.pdf"),height=5,width=10,useDingbats = FALSE)
  ggplot(sdf,aes(x=frac,y=prev,color=as.factor(stratum))) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=sdf$lb,ymax=sdf$ub)) +
    scale_color_manual(values=c("goldenrod3","darkblue"),name="Family History",labels=c("Negative","Positive")) + xlab("GRS quantile") + ylab("Prevalence")  +
    scale_x_continuous(breaks=breaks) +  facet_wrap(~label,nrow=1) + theme(legend.position="bottom",legend.text=element_text(size=14),title=element_text(size=14),axis.title=element_text(size=14),axis.text.x=element_text(size=6,angle=45)) +
    theme(strip.background =element_rect(fill="white"))
  dev.off()
  
} ###end of qtile loop



#TO DO: int_term isn't working so just made it NULL

####### models across the ages 
age_analysis<-function(df,grs_col,strat_col,pheno_col,out,part_age,part_age_std,cut_age=60){
  #df$age_bins<-cut_number(df[[part_age]],4)
  df$age_bins<-cut(df[[part_age]],breaks=c(35,40,45,50,55,60,65,70)) #decade bins
  #df$age_bins<-cut(df[[part_age]],breaks=c(15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)) #five years 
  #df$age_bins<-cut(df[[part_age]],breaks=c(20,40,60,80,100)) #20 year bins
  #df$age_bins<-cut(df[[part_age]],15) #roughly 5 year, this is failing at model indicator function figure out why
  #df$age_bins<-cut(df[[part_age]],breaks=seq(40,80)) #check cross over age by looking at every age
  df<-df[!is.na(df$age_bins),]
  bins<-unique(df$age_bins)
  tcdf<-data.frame() #continuous
  idf<-data.frame() #indicator
  edf<-data.frame() #continuous extended
  #covar<-c(batch,pcs,sex,part_age_std)
  #name<-"batch_PCs_sex_partAge"
  covar<-c(batch,pcs,sex)
  name<-"batch_PCs_sex"
  #by bins
  for (b in 1:length(bins)){
    age_sub<-df[df$age_bins==bins[b],] #age subset
    print(bins[b])
    print(nrow(age_sub))
    if (nrow(age_sub) > 1){
      tmp<-test_covar(age_sub,covar,out,name,int_term=NULL)
      tmp$continuous$covar<-name
      tmp$indicator$covar<-name
      tmp$continuous_extended$covar<-name
      tmp$continuous$age_bin<-bins[b]
      tmp$indicator$age_bin<-bins[b]
      tmp$continuous_extended$age_bin<-bins[b]
      tmp$continuous_extended$count<-nrow(age_sub)
      tmp$continuous_extended$FH_count<-table(age_sub[[strat_col]])[2]
      tmp$continuous_extended$pheno_count<-table(age_sub[[pheno_col]])[2]
      tcdf<-rbind(tcdf,tmp$continuous)
      idf<-rbind(idf,tmp$indicator)
      edf<-rbind(edf,tmp$continuous_extended)
      print(edf)
    }
  }
  #old
  tmp<-test_covar(df[df[[part_age]]>cut_age,],covar,out,paste(sep="_",name,"old"),int_term=NULL)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tmp$continuous$age_bin<-"old"
  tmp$indicator$age_bin<-"old"
  tmp$continuous_extended$age_bin<-"old"
  tmp$continuous_extended$count<-nrow(qsub[qsub[[part_age]]>cut_age,])
  tmp$continuous_extended$FH_count<-table(qsub[qsub[[part_age]]>cut_age,][[strat_col]])[2]
  tmp$continuous_extended$pheno_count<-table(qsub[qsub[[part_age]]>cut_age,][[pheno_col]])[2]
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  #young
  tmp<-test_covar(df[df[[part_age]]<=cut_age,],covar,out,paste0(name,"young"),int_term=NULL)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tmp$continuous$age_bin<-"young"
  tmp$indicator$age_bin<-"young"
  tmp$continuous_extended$age_bin<-"young"
  tmp$continuous_extended$count<-nrow(qsub[qsub[[part_age]]<=cut_age,])
  tmp$continuous_extended$FH_count<-table(qsub[qsub[[part_age]]<=cut_age,][[strat_col]])[2]
  tmp$continuous_extended$pheno_count<-table(qsub[qsub[[part_age]]<=cut_age,][[pheno_col]])[2]
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  return(list(tcdf=tcdf,idf=idf,edf=edf))
  
}

##### looking at models across age bins using age_analysis function

cut_age<-60
aobj<-age_analysis(qsub,grs_col,strat_col,pheno_col,out,part_age,part_age_std,cut_age)
#if model_df not found we need to reset i to  1, we should put this  in the  big loop over the q quantiles
edf_sub<-aobj$edf[(aobj$edf$pred=="FH"|aobj$edf$pred=="GRS") & aobj$edf$model!="interaction" & aobj$edf$age_bin!="old" & aobj$edf$age_bin!="young",]
edf_sub$model<-as.factor(edf_sub$model)
edf_sub$prop_FH<- edf_sub$FH_count/edf_sub$count #proportion family history
levels(edf_sub$model)<-c("family history + GRS","family history only","GRS only")
edf_sub$xlab<-paste(sep="\n",edf_sub$age_bin,paste0("(",format(edf_sub$prop_FH*100,digits=3),"% of ", edf_sub$count,")")) ##make column with label
#currently not using because order was off
xlabel_list<-c()
for (x in 1:nrow(unique(edf_sub[c("age_bin","prop_FH","count")]))){
  xlabel_list[x]<-paste(sep="\n",unique(edf_sub[c("age_bin","prop_FH","count")])[x,1],paste0("(",format(unique(edf_sub[c("age_bin","prop_FH","count")])[x,2]*100,digits=3),"% of ", unique(edf_sub[c("age_bin","prop_FH","count")])[x,3],")"))
}

pdf(file=paste0(out,"_age.pdf"),height=5,width=10,useDingbats = FALSE)
ggplot(edf_sub,aes(x=xlab,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB)) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank()) +  
  scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="[Enrollment Age Bin]\n(% Family History of n)") + scale_y_continuous(labels=function(x) sprintf("%.1f", x))
dev.off()

pdf(file=paste0(out,"_age_smooth.pdf"),height=6,width=10,useDingbats = FALSE)
edf_sub$age<-sapply(strsplit(as.character(edf_sub$age_bin),",",fixed=TRUE), '[', 1)
edf_sub$age<-as.numeric(sapply(strsplit(as.character(edf_sub$age),"(",fixed=TRUE), '[', 2))
ggplot(edf_sub,aes(x=age,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB)) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank()) +  
  scale_color_manual(values=c(pamp[2],pamp[5],pamp[6])) + labs(y="Odds Ratio",x="[Enrollment Age Bin]\n(% Family History of n)") +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x))  + geom_smooth(alpha=.25,se=FALSE)
dev.off()

##cross over age in additive model
GRS_smooth <- predict(loess(OR~age,edf_sub[edf_sub$model=="family history + GRS" & edf_sub$pred=="GRS",]), seq(20,80,2))
FH_smooth<-predict(loess(OR~age,edf_sub[edf_sub$model=="family history + GRS" & edf_sub$pred=="FH",]), seq(20,80,2))
cross<-data.frame(GRS_smooth,FH_smooth, age=seq(20,80,2))
ggplot(cross,aes(x=GRS_smooth,y=FH_smooth,label=age)) + geom_point() + geom_text() + geom_abline(slope=1,intercept=0,linetype="dashed",color="red") 
#To do: need to optimize  to find where x=y


edf_sub<-edf_sub[edf_sub$model=="GRS only" | edf_sub$model=="family history + GRS" | edf_sub$model=="family history only",]
edf_sub$pred<-revalue(edf_sub$pred,c("FH"="family history"))
pdf(file=paste0(out,"_age_small.pdf"),height=4,width=8,useDingbats = FALSE)
ggplot(edf_sub,aes(x=xlab,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB),width=0.4) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank(),strip.background =element_rect(fill="white")) +  
  scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="[Enrollment Age Bin]\n(% Positive Family History of n)") + scale_y_continuous(labels=function(x) sprintf("%.1f", x)) 
dev.off()

if (length(unique(edf_sub$age_bin)) > 10) {
  edf_sub<-edf_sub[edf_sub$model=="family history + GRS",]
  edf_sub$age<-sapply(strsplit(as.character(edf_sub$age_bin),",",fixed=TRUE), '[', 1)
  edf_sub$age<-as.numeric(sapply(strsplit(as.character(edf_sub$age),"(",fixed=TRUE), '[', 2))
  pdf(file=paste0(out,"_many_ages.pdf"),height=4,width=9,useDingbats = FALSE)
  print(ggplot(edf_sub,aes(x=as.numeric(age),y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
          geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB),width=0.4) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank(),strip.background =element_rect(fill="white")) +  
          scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="[Baseline Age Bin]\n(% Positive Family History)") + scale_y_continuous(labels=function(x) sprintf("%.1f", x)))
  dev.off()
}

##plot sex stratified
sex_edf_sub<-aobj$sex_edf[(aobj$sex_edf$pred=="FH"|aobj$sex_edf$pred=="GRS") & aobj$sex_edf$model!="interaction" & aobj$sex_edf$age_bin!="old" & aobj$sex_edf$age_bin!="young",]
sex_edf_sub$prop_FH<- sex_edf_sub$FH_count/sex_edf_sub$count #proportion family history
sex_edf_sub$xlab<-paste(sep="\n",sex_edf_sub$age_bin,paste0("(",format(sex_edf_sub$prop_FH*100,digits=3),"% of ", sex_edf_sub$count,")")) ##make column with label
sex_edf_sub$sex<-as.factor(as.character(sex_edf_sub$sex))
sex_edf_sub$sex<-revalue(sex_edf_sub$sex,c("2"="Female","1"="Male"))
sex_edf_sub$model<-as.factor(sex_edf_sub$model)
levels(sex_edf_sub$model)<-c("family history + GRS","family history only","FHiGRS","GRS only")
pdf(file=paste0(out,"_age_small_sex.pdf"),height=8,width=8,useDingbats = FALSE)
ggplot(sex_edf_sub,aes(x=xlab,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~sex+model,nrow=2,scales="free_x") + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=sex_edf_sub$LB,ymax=sex_edf_sub$UB),width=0.4) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank(),strip.background =element_rect(fill="white")) +  
  scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="[Enrollment Age Bin]\n(% Positive Family History of n)") + scale_y_continuous(labels=function(x) sprintf("%.1f", x)) 
dev.off()

edf_sub$model<-as.factor(edf_sub$model)
edf_sub$prop_FH<- edf_sub$FH_count/edf_sub$count #proportion family history
levels(edf_sub$model)<-c("family history + GRS","family history only","FHiGRS","GRS only")

edf_sub<-edf_sub[edf_sub$model=="GRS only" | edf_sub$model=="family history + GRS" | edf_sub$model=="family history only",]
edf_sub$pred<-revalue(edf_sub$pred,c("FH"="family history"))

##old vs young
edf_sub<-aobj$edf[(aobj$edf$pred=="FH"|aobj$edf$pred=="GRS") & aobj$edf$model!="interaction" & (aobj$edf$age_bin=="old" | aobj$edf$age_bin=="young"),]
edf_sub$pred<-revalue(edf_sub$pred,c("FH"="family history"))
edf_sub$prop_FH<- edf_sub$FH_count/edf_sub$count #proportion family history
xlabel_list<-c()
for (x in 1:nrow(edf_sub[c("age_bin","prop_FH")])){
  xlabel_list[x]<-ifelse(edf_sub["age_bin"][x,1]==as.factor("young"),paste0("<=",cut_age,"yrs"),paste0(">",cut_age,"yrs"))
}
edf_sub$label<-xlabel_list
edf_sub$model<-revalue(edf_sub$model,c("family"="family history","additive"="family history + GRS"))
edf_sub$model<-factor(edf_sub$model,levels(as.factor(edf_sub$model))[c(2,1,3)])
pdf(file=paste0(out,"_old_vs_young_small.pdf"),height=3,width=6,useDingbats = FALSE)
ggplot(edf_sub,aes(x=label,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB),width=0.3) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank(),strip.background =element_rect(fill="white")) +  
  scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="Enrollment Age\n(% Positive Family History)") + scale_y_continuous(labels=function(x) sprintf("%.1f", x))
dev.off()


#look at cumulative 




###### BMI 
bmi_analysis<-function(df,grs_col,fhigrs_col,strat_col,pheno_col,out,part_age,part_age_std,bmi,cut_bmi=30){
  df$bmi_bins<-cut(df[[bmi]],breaks=c(10,20,25,30,35,40,45,60)) #bmi bins
  df<-df[!is.na(df$bmi_bins),]
  bins<-unique(df$bmi_bins)
  tcdf<-data.frame() #continuous
  idf<-data.frame() #indicator
  edf<-data.frame() #continuous extended
  covar<-c(batch,pcs,sex,part_age_std)
  name<-"batch_sex_PCs_partAge"
  #by bins
  for (b in 1:length(bins)){
    bmi_sub<-df[df$bmi_bins==bins[b],] #age subset
    print(bins[b])
    print(nrow(bmi_sub))
    tmp<-test_covar(bmi_sub,covar,out,name,int_term=NULL)
    tmp$continuous$covar<-name
    tmp$indicator$covar<-name
    tmp$continuous_extended$covar<-name
    tmp$continuous$bmi_bin<-bins[b]
    tmp$indicator$bmi_bin<-bins[b]
    tmp$continuous_extended$age_bin<-bins[b]
    tmp$continuous_extended$count<-nrow(bmi_sub)
    tmp$continuous_extended$FH_count<-table(bmi_sub[[strat_col]])[2]
    tmp$continuous_extended$pheno_count<-table(bmi_sub[[pheno_col]])[2]
    tcdf<-rbind(tcdf,tmp$continuous)
    idf<-rbind(idf,tmp$indicator)
    edf<-rbind(edf,tmp$continuous_extended)
  }
  #obese
  tmp<-test_covar(df[df[[bmi]]>cut_bmi,],covar,out,paste(sep="_",name,"obese"),int_term=NULL)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tmp$continuous$age_bin<-"obese"
  tmp$indicator$age_bin<-"obese"
  tmp$continuous_extended$age_bin<-"obese"
  tmp$continuous_extended$count<-nrow(qsub[qsub[[bmi]]>cut_bmi,])
  tmp$continuous_extended$FH_count<-table(qsub[qsub[[bmi]]>cut_bmi,][[strat_col]])[2]
  tmp$continuous_extended$pheno_count<-table(qsub[qsub[[bmi]]>cut_bmi,][[pheno_col]])[2]
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  #non obese
  tmp<-test_covar(df[df[[bmi]]<=cut_bmi,],covar,out,paste0(name,"non-obese"),int_term=NULL)
  tmp$continuous$covar<-name
  tmp$indicator$covar<-name
  tmp$continuous_extended$covar<-name
  tmp$continuous$age_bin<-"non_obese"
  tmp$indicator$age_bin<-"non_obese"
  tmp$continuous_extended$age_bin<-"non_obese"
  tmp$continuous_extended$count<-nrow(qsub[qsub[[bmi]]<=cut_bmi,])
  tmp$continuous_extended$FH_count<-table(qsub[qsub[[bmi]]<=cut_bmi,][[strat_col]])[2]
  tmp$continuous_extended$pheno_count<-table(qsub[qsub[[bmi]]<=cut_bmi,][[pheno_col]])[2]
  tcdf<-rbind(tcdf,tmp$continuous)
  idf<-rbind(idf,tmp$indicator)
  edf<-rbind(edf,tmp$continuous_extended)
  
  return(list(tcdf=tcdf,idf=idf,edf=edf))
}
bobj<-bmi_analysis(qsub,grs_col,fhigrs_col,strat_col,pheno_col,out,part_age,part_age_std,bmi)
edf_sub<-bobj$edf[(bobj$edf$pred=="FH"|bobj$edf$pred=="GRS") & bobj$edf$model!="interaction",]
pdf(file=paste0(out,"_bmi.pdf"),height=5,width=8,useDingbats = FALSE)
ggplot(edf_sub,aes(x=bmi_bin,y=OR,color=pred)) + geom_point(alpha=0.5,size=2) + theme_bw() + facet_wrap(~model,nrow=1) + geom_hline(yintercept=1,linetype="dashed",color="black",alpha=0.25) +
  geom_errorbar(aes(ymin=edf_sub$LB,ymax=edf_sub$UB),width=0.3) +  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom",legend.title=element_blank(),strip.background =element_rect(fill="white")) +  
  scale_color_manual(values=c(pamp[2],pamp[6])) + labs(y="Odds Ratio",x="BMI") + scale_y_continuous(labels=function(x) sprintf("%.1f", x))
dev.off()