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
    make_option("--legend",type="character",default="Binary stratum",help="Legend title which is stratum [default='Binary stratum']")
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
xlab<-arguments$options$xlabel
ylab<-arguments$options$ylabel
legend<-arguments$options$legend
dig<-arguments$option$digits
covar<-as.numeric(strsplit(arguments$options$covariates,",")[[1]])
invNorm<-arguments$options$invNorm
header<-arguments$options$header

###########################################################
#################### FUNTIONS #############################
###########################################################

tile_info<-function(num_list){
    val<-c(2,3,4,5,6,7,8,9,10,12,16,20,100,1000)
    names<-c("median","tertiles","quartiles","quantiles","sextiles","septiles","octiles","deciles","dodeciles","hexadeciles","ventiles","percentiles")
    name_list<-c()
    for (i in 1:length(num_list)){
        if (length(which(val==num_list[i]))==0){ #use generic quantile label 
            name_list[i]<-"quantile"
        } else {
            name_list[i]<-names[which(val==num_list[i])] #use specific quantile label 
        }
    }
    return(name_list)
}

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


##quantile is how many divisions we should make of dataset (e.g. 20-quantile is ventiles for 20 groups)
##returns object with prevalences and standard error for distribution 
#does not stratify
prev_per_quantile<-function(df,GRS_col,prev_col,qtile){
  if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
    print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
  }
  if (sum(qtile)<2*length(qtile)){ #check qtile
    print("q-quantiles should be number of divisions for data set and must be greater than 1")
  }
  #initialize data structures
  p<-(100/qtile)/100
  index<-c(seq(from=0,to=1,by=p)*100)
  prevalences<-rep(NA,qtile+1) #initialize prevalence vector
  ns<-rep(NA,qtile+1) #initialize count vector
  ses<-rep(NA,qtile+1)#initialize se vector
  tiles<-quantile(df[[GRS_col]],seq(from=0,to=1,by=p)) #quantile values
  for (i in 1:length(index)-1) {
      prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] <= tiles[i+1]][[prev_col]]
      prevalences[i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
      ns[i]<-length(prev_list)
      ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/length(prev_list)) #what is SE for this prevalence
    }
  #create object 
  pq<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
  class(pq)<-"prev_quantile_obj"
  return(pq)
}


## plot prevalence per GRS quantile bin, color by stratum if stratum==TRUE
plotting<-function(dat,out,qtiles,stratum=FALSE,main,xlab,ylab,legend,ymax=1){
  #dat<-dat[complete.cases(dat),] #remove the last bins which are NA
  dat$frac<-as.numeric(sub("%","",dat$percents)) #convert factor percentages to numeric 
  dat<-dat[dat$frac!=1.00,] 
  dat$ub<-dat$prev+(1.96*dat$se)
  dat$lb<-dat$prev-(1.96*dat$se)
  
  if (ymax==1){ # if ymax is not given to function, just plot to scale of data, otherwise, script can be given custom ymax to match scale of another plot
    ymax=max(df$prev)
  }
  
  if (stratum==TRUE){ #stratify
    by(dat, dat$q, #number of q-quantiles (e.g. break data into 4 bins, 10 bins, etc.)
      function (x) {
        name=unique(x$q)
        if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
        pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5)
        print(ggplot(x,aes(x=frac,y=prev,color=as.factor(stratum))) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=x$lb,ymax=x$ub)) + 
          scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +labs(title=main) + xlab(xlab) + ylab(ylab)  + 
          scale_x_continuous(breaks=breaks) + coord_cartesian(ylim=c(0,ymax))) 
        dev.off()
      }
    )
    
  } else { #all data
    by(dat, dat$q, #number of q-quantiles (e.g. break data into 4 bins, 10 bins, etc.)
       function (x) {
         name=unique(x$q)
         if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
         pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5)
         print(ggplot(x,aes(x=frac,y=prev)) + geom_point(color="grey") + geom_errorbar(aes(ymin=x$lb,ymax=x$ub),color="grey")  + 
                 theme_bw() + labs(title=main) + xlab(xlab) + ylab(ylab) + scale_x_continuous(breaks=breaks) + coord_cartesian(ylim=c(0,ymax))) 
         dev.off()
       }
    )
  }
}  


    
#compare prevalence in top/bottom of a distribution give a cut point value 
top_quantile_stratum<-function(df,GRS_col,prev_col,strat_col,qfirst=FALSE,value){
  if (value < 0.5 | value >= 1){
    print("Percentile for dividing GRS distribution must be >= 0.5 < 1 ")
  }
  ##odds ratio for top vs bottom of distribution by stratum
  ##columns are bottom/top of distribution, rows are stratum 
  index<-c(value*100,100)
  prevalences<-matrix(NA,2,length(index)) #initialize prevalence matrix
  counts<-list(matrix(NA,length(index),length(index)),matrix(NA,length(index),length(index))) #list of matrices, one matrix per stratum
  ses<-matrix(NA,2,length(index))#initialize se matrix
  if (qfirst==TRUE) { #make quantiles before stratifying data
    tiles<-quantile(df[[GRS_col]],c(0,value,1))
    for (i in 1:length(index)) {
      for (r in c(1,2)){ #iterate over stratum
        prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] <= tiles[i+1] & df[[strat_col]]==(r-1)][[prev_col]]
        prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ses[r,i]<-sqrt((prevalences[i]*(1-prevalences[i]))/length(prev_list)) #what is SE for this prevalence
        counts[[r]][i,]<-as.vector(table(prev_list)) #counts for OR
      } 
    }  
  } else{ #stratifying before calculating quantiles
      for (r in c(1,2)){ #iterate over stratum
        sdf<-df[df[[strat_col]]==(r-1),] #make data set for one stratum
        tiles<-quantile(sdf[[GRS_col]],c(0,value,1)) #quantile values
        for (i in 1:length(index)) {
          prev_list<-sdf[sdf[[GRS_col]] > tiles[i] & sdf[[GRS_col]] <= tiles[i+1]][[prev_col]]
          prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
          ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
          counts[[r]][i,]<-as.vector(table(prev_list))
        }
      }
  }
  #odds ratio and se
  OR_list<-rep(NA,length(index))
  SE_list<-rep(NA,length(index))
  UB_list<-rep(NA,length(index))
  LB_list<-rep(NA,length(index))
  for (r in c(1,2)) { #iterate over stratum
    m<-counts[[r]]
    OR_list[r]<-(m[2,2]/m[2,1])/(m[1,2]/m[1,1])  #case/control top distribution over case/control bottom distribution
    SE_list[r]<-sqrt(sum(1/m)) #log odds scale
    LB_list[r]<-exp(log(OR_list[r])-1.96*SE_list[r])
    UB_list[r]<-exp(log(OR_list[r])+1.96*SE_list[r])
  }
  ##odds ratio for top of stratum=1 versus everything else 
  num<-counts[[2]][2,2]/counts[[2]][2,1] #case/control of top distribution with stratum=1
  denom<-sum(counts[[2]][1,2],counts[[1]][,2])/sum(counts[[2]][1,1],counts[[1]][,1]) #case/control of top distribution with stratum=0 + bottom distribution stratum=0|1
  stratum_OR<-num/denom
  stratum_SE<-sqrt(sum(1/counts[[2]][2,2],1/counts[[2]][2,1],1/sum(counts[[2]][1,2],counts[[1]][,2]),1/sum(counts[[2]][1,1],counts[[1]][,1])))
  stratum_LB<-exp(log(stratum_OR)-1.96*stratum_SE)
  stratum_UB<-exp(log(stratum_OR)+1.96*stratum_SE)
  tqs<-list(prev=prevalences,se=ses,i=index,c=counts,OR=OR_list,SE=SE_list,LB=LB_list,UB=UB_list,
            stratum_OR=stratum_OR,stratum_SE=stratum_SE,stratum_LB=stratum_LB,stratum_UB=stratum_UB)
  class(tqs)<-"top_quantile_stratum_obj"
  return(tqs)
}

##compare prevalence in top/bottom of a distribution given a cut point value 
top_quantile<-function(df,GRS_col,prev_col,value){
  if (value < 0.5 | value >= 1){
    print("Percentile for dividing GRS distribution must be >= 0.5 < 1 ")
  }
  ##odds ratio for top vs bottom of distribution
  index<-c(value*100,100)
  prevalences<-rep(NA,length(index)) #initialize prevalence vector
  counts<-matrix(NA,length(index),length(index)) #initialize counts matrix
  ses<-rep(NA,length(index)) #initialize se vector
  tiles<-quantile(df[[GRS_col]],c(0,value,1))
  for (i in 1:length(index)) {
      prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] <= tiles[i+1]][[prev_col]]
      prevalences[i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
      ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/length(prev_list)) #what is SE for this prevalence
      counts[i,]<-as.vector(table(prev_list)) #counts for OR
  }
  OR<-(counts[2,2]/counts[2,1])/(counts[1,2]/counts[1,1])  #case/control top distribution over case/control bottom distribution
  SE<-sqrt(sum(1/counts)) #log odds scale
  LB<-exp(log(OR)-1.96*SE)
  UB<-exp(log(OR)+1.96*SE)
  tq<-list(OR=OR,SE=SE,LB=LB,UB=UB,c=counts)
  class(tq)<-"top_quantile_obj" #What is the odds ratio and 95% CI when you compare top and bottom of GRS distribution at a given cut point
  return(tq)
}
  
## use logistic regression to model disease ~ I(top of GRS distirbution) + covariates
model<-function(df,grs_col,strat_col,pheno_col,covar,value,qfirst=FALSE){
  if (value < 0.5 | value >= 1){
    print("Percentile for dividing GRS distribution must be >= 0.5 < 1 ")
  }
  mobj<-list() #initialize object
  
  ##model I(dist) for all data
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
  adf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist")) #all
  m<-matrix(table(df$dist,df[[pheno_col]]),byrow=FALSE,nrow=2)
  adf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  adf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["standard"]]<-adf
  
  ##if we want to stratify before we decide on quantiles, remake I(dist) column, otherwise use previous
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
  df$FHIGR<-2
  df[df[[strat_col]]==1 & df$dist==1]$FHIGR<-1
  df[df[[strat_col]]==1 & df$dist==0]$FHIGR<-0
  df[df[[strat_col]]==0]$FHIGR<-0
  df[df$FHIGR==2]$FHIGR<-NA
  df$dist<-df$FHIGR #replace I(dist) with I(dist & stratum=1)
  
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,dist_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  fhidf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"dist"))
  m<-matrix(table(df$dist,df[[pheno_col]]),byrow=FALSE,nrow=2)
  fhidf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  fhidf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["FHIGR"]]<-fhidf
  class(mobj)<-"model_obj"
  return(mobj)
}

## screen N people how many do we catch and how many do we miss given screening strategy (prioritize by family history or no)
## expects "top_quantile_stratum_obj"
clinical_impact<-function(obj,N=10000){
  if (class(obj)!="top_quantile_stratum_obj"){
    stop("clinical impact function requires top_quantile_stratum_obj from top_quantile_stratum function\n")
  }
  counts<-obj$c # 2 matrices, first is s=0, second is s=1
  #first row of matrix is counts from bottom of distribution at given cut point
  #second row of matrix is counts from top of distirbution at given cut point
  #first column of matrix is counts of controls
  #second column of matrix is counts of cases 
  
  #scenario 1, prioritize screening based on stratum=1
  screen<-c(counts[[2]][2,1],counts[[2]][2,2]) #control,case of top distribution with stratum=1
  no_screen<-c(sum(counts[[2]][1,1],counts[[1]][,1]),sum(counts[[2]][1,2],counts[[1]][,2])) #control,case of top distribution with stratum=0 + bottom distribution stratum=0|1
  m<-matrix(c(no_screen,screen),nrow=2,ncol=2,byrow=TRUE)
  mfrac<-m/sum(m)
  scenario1<-mfrac*N
  
  #scenario 2, screen top percentle of GRS
  all<-counts[[1]]+counts[[2]]
  allfrac<-all/sum(all)
  scenario2<-allfrac*N #no screen is top row, screen is bottom row, control is first column, case is second column
  
  false_pos<-(scenario1-scenario2)[2,1]
  false_neg<-(scenario1-scenario2)[1,2]
  return(c(false_pos,false_neg))
}

###########################################################
########################## MAIN ###########################
###########################################################

##read data 
df<-fread(file,header=header)
#print(dim(df))
## To DO: check column assumptions 

##subset to data with stratum available
subset<-df[!is.na(df[[strat_col]])] #remove if NA for stratum
print(dim(subset))
size<-length(quantiles) #number of quantiles being tested, size of obj
#ggplot(subset, aes_string(x=names(subset)[grs_col])) + geom_density()

if (invNorm==TRUE){
  ## inverse rank normalize GRS in thepopulation
  subset$invNormGRS<-rankNorm(subset[[grs_col]])
  grs_col<-which(names(subset)=="invNormGRS") #new GRS col
}

############# Prevalences versus GRS############# 
## uses prev_per_quantile_stratum and prev_per_quantile functions

##make quantiles then stratify
qobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=TRUE)
#print(obj)
for (i in 1:size){ #across q-quantiles
    for (j in c(1,2)){ #across stratum
        list_length<-length(qobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
        prev<-qobj[[i]]$prev[j,][-list_length]
        se<-qobj[[i]]$se[j,][-list_length]
        n<-qobj[[i]]$n[j,][-list_length]
        tiles<-qobj[[i]]$tiles[-1]
        percents<-names(tiles)
        bins<-rep(quantiles[i],list_length-1)
        strat<-rep(j-1,list_length-1)
        if (i==1 & j==1) {
          qdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
        } else {
          qdf<-rbind(qdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
        }
    }
}

##stratify then calculate quantiles
sobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=FALSE)
#print(obj)
for (i in 1:size){ #across q-quantiles
  for (j in c(1,2)){ #across stratum
    list_length<-length(sobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-sobj[[i]]$prev[j,][-list_length]
    se<-sobj[[i]]$se[j,][-list_length]
    n<-sobj[[i]]$n[j,][-list_length]
    tiles<-sobj[[i]]$tiles[-1]
    percents<-names(tiles)
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(j-1,list_length-1)
    if (i==1 & j==1) {
      sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
    } else {
      sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
    }
  }
}

ymax<-max(max(qdf$prev+(1.96*qdf$se)),max(sdf$prev+(1.96*sdf$se)))
file_name<-paste(sep=".",out,"quantileFirst.txt")
write.table(format(qdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(df,paste(sep="_",out,"quantileFirst"),quantiles,TRUE,main,xlab,ylab,legend,ymax)

file_name<-paste(sep=".",out,"stratifyFirst.txt")
write.table(format(sdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(df,paste(sep="_",out,"stratifyFirst"),quantiles,TRUE,main,xlab,ylab,legend,ymax)

##calculate quantiles for all data, no stratification
all_obj<-lapply(quantiles,prev_per_quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col)
##print(all_obj)
for (i in 1:size){ #across q-quantiles
    list_length<-length(all_obj[[i]]$prev)  #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-all_obj[[i]]$prev[-list_length]
    se<-all_obj[[i]]$se[-list_length]
    n<-all_obj[[i]]$n[-list_length]
    tiles<-all_obj[[i]]$tiles[-1]
    percents<-names(tiles)
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(NA,list_length-1)
    if (i==1){
      all_df<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
    } else {
      all_df<-rbind(all_df,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
    }
}
file_name<-paste(sep=".",out,"noStratify.txt")
write.table(format(all_df,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(all_df,paste(sep="_",out,"all"),quantiles,FALSE,main,xlab,ylab,legend,ymax)  

################## Odds Ratios by model ############# 
## uses model function

if (1 %in% cutpts){ #code doesn't work with 
  cutpts<-cutpts[-which(cutpts==1)]
}
n<-length(cutpts)
logical_list<-c(TRUE,FALSE)
label_list<-c("quantileFirst","stratifyFirst")
obj<-list()
for (l in c(1,2)){ #quantile first or no
  obj[[l]]<-lapply(cutpts,model,df=subset,grs_col=grs_col,pheno_col=pheno_col,strat_col=strat_col,covar=covar,qfirst=logical_list[l])
  for (c in 1:n) { #per cut point
    for (j in 1:length(obj[[l]][[c]])){ #4 data frames come out of model function
      if (j==1 & l==1 & c==1){ #start data frame if beginning of loops
          d<-data.frame(OR=exp(obj[[l]][[c]][[j]]['dist','Estimate']),
             LB=exp(obj[[l]][[c]][[j]]['dist','Estimate']-1.96*obj[[l]][[l]][[j]]['dist','Std..Error']),
             UB=exp(obj[[l]][[c]][[j]]['dist','Estimate']+1.96*obj[[l]][[l]][[j]]['dist','Std..Error']),
             pval=obj[[l]][[c]][[j]]['dist','Pr...z..'],
             top_prev<-unique(obj[[l]][[c]][[j]]['top_prev']),
             bottom_prev<-unique(obj[[l]][[c]][[j]]['bottom_prev']),
             label=label_list[l],
             name=names(obj[[l]][[c]])[j],
             cutpt=cutpts[c])
    } else {
      d<-rbind(d,data.frame(OR=exp(obj[[l]][[c]][[j]]['dist','Estimate']),
                            LB=exp(obj[[l]][[c]][[j]]['dist','Estimate']-1.96*obj[[l]][[l]][[j]]['dist','Std..Error']),  
                            UB=exp(obj[[l]][[c]][[j]]['dist','Estimate']+1.96*obj[[l]][[l]][[j]]['dist','Std..Error']),
                            pval=obj[[l]][[c]][[j]]['dist','Pr...z..'],
                            top_prev<-unique(obj[[l]][[c]][[j]]['top_prev']),
                            bottom_prev<-unique(obj[[l]][[c]][[j]]['bottom_prev']),
                            label=label_list[l],
                            name=names(obj[[l]][[c]])[j],
                            cutpt=cutpts[c]))
    }}}
    sub<-d[(d$name=="standard"|d$name=="FHIGR") & d$label==label_list[l],]
    pdf_fn<-paste(sep=".",out,"model.compareOR",label_list[l],"pdf")
    pdf(file=pdf_fn,height=4,width=6)
    print(ggplot(sub,aes(x=cutpt,y=OR,color=as.factor(name))) + geom_point(alpha=0.7)+   geom_errorbar(aes(ymin=sub$LB,ymax=sub$UB)) +
          theme_bw() + scale_color_manual(values=c("grey","darkblue"),name="") + geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7) +
          labs(title=main,x="Cut Point",y="Odds Ratio"))
    dev.off()
    
    sub<-d[(d$name=="stratum0"| d$name=="stratum1") & d$label==label_list[l],]
    pdf_fn<-paste(sep=".",out,"model.OR",label_list[l],"pdf")
    pdf(file=pdf_fn,height=4,width=6)
    print(ggplot(sub,aes(x=cutpt,y=OR,color=as.factor(name))) + geom_point(alpha=0.7)+   geom_errorbar(aes(ymin=sub$LB,ymax=sub$UB)) +
            theme_bw() + scale_color_manual(values=c("darkblue","goldenrod3"),name=legend) + geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7) +
            labs(title=main,x="Cut Point",y="Odds Ratio"))
    dev.off()
}
fn<-paste(sep=".",out,"modelOR.txt")
write.table(format(d,digits=dig),fn,quote=FALSE,col.names=T,row.names=F,sep="\t")


################## Odds Ratios by contingency tables  #############
## uses top_quantile_stratum, top_quantile, clinical_impact functions
#cutpts<-seq(0.8,1,by=.01)
if (1.0 %in% cutpts){ #code doesn't work with 
  cutpts<-cutpts[-which(cutpts==1)]
}
n<-length(cutpts)
logical_list<-c(TRUE,FALSE)
label_list<-c("quantileFirst","stratifyFirst")
for (l in c(1,2)){ #do for each division logic
  obj<-lapply(cutpts,top_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=logical_list[l])
  allobj<-lapply(cutpts,top_quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col)
  for (i in 1:n){ #across cut points
    for (j in c(1,2)){ #across stratum
      SE<-obj[[i]]$SE[j]
      OR<-obj[[i]]$OR[j]
      LB<-obj[[i]]$LB[j]
      UB<-obj[[i]]$UB[j]
      strat<-j-1
      if (i==1 & j==1) {
        ORdf<-data.frame(OR=OR,SE=SE,LB=LB,UB=UB,stratum=strat,cutpt=cutpts[i],row.names=NULL)
      } else {
        ORdf<-rbind(ORdf,data.frame(OR=OR,SE=SE,LB=LB,UB=UB,stratum=strat,cutpt=cutpts[i],row.names=NULL))
      }
    }
    if (i==1){
      stratumORdf<-data.frame(stratOR=obj[[i]]$stratum_OR,stratSE=obj[[i]]$stratum_SE,stratLB=obj[[i]]$stratum_LB,stratUB=obj[[i]]$stratum_UB,cutpt=cutpts[i],row.names=NULL)
      allORdf<-data.frame(OR=allobj[[i]]$OR,SE=allobj[[i]]$SE,LB=allobj[[i]]$LB,UB=allobj[[i]]$UB,cutpt=cutpts[i],row.names=NULL)
    }else {
      stratumORdf<-rbind(stratumORdf,data.frame(stratOR=obj[[i]]$stratum_OR,stratSE=obj[[i]]$stratum_SE,stratLB=obj[[i]]$stratum_LB,stratUB=obj[[i]]$stratum_UB,cutpt=cutpts[i],row.names=NULL))
      allORdf<-rbind(allORdf,data.frame(OR=allobj[[i]]$OR,SE=allobj[[i]]$SE,LB=allobj[[i]]$LB,UB=allobj[[i]]$UB,cutpt=cutpts[i],row.names=NULL))
    }
  }
  
  #clinical impact of prioritized screening group by stratum compared to just top of GRS distribution
  clin<-lapply(obj,clinical_impact)
  clin_df <- data.frame(matrix(unlist(clin), nrow=length(clin), byrow=T))
  names(clin_df)<-c("falsepos","falseneg")
  clin_df$cutpt<-cutpts
  clin_df$method<-label_list[l]
  file_name<-paste(sep=".",out,"clinical_impact.txt")
  if (l==1){
    write.table(format(clin_df,digits=dig),file=file_name,col.names=TRUE,quote=FALSE,row.names=FALSE,sep="\t",append=FALSE)
  }else {
    write.table(format(clin_df,digits=dig),file=file_name,col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t",append=TRUE)
  }
  #plot of sensitivty and specificity, convert decimal cutpoints to whole numbers
  pdf_fn<-paste(sep=".",out,"clinical_impact",label_list[l],"pdf")
  nudge_factor<- diff(range(clin_df$falsepos))/10 #if the x axis scale is kind of small we don't need to nudge labels too far from points
  if (n<=5){
    scale<-100*as.vector(cutpts) #use all cutpts as breaks in legend
  } else {
    scale<-100*(c(cutpts[1],cutpts[floor(n/2)],cutpts[n])) #use max, middle, and min as cutpts in legend 
  }
  pdf(file=pdf_fn,height=4,width=6)
  print(ggplot(clin_df,aes(x=falsepos,y=falseneg,label=cutpt*100,size=100*cutpt,color=as.factor(method))) + geom_point(alpha=0.7)+   
          theme_bw() + scale_color_manual(values=c("darkorchid4")) + geom_text(nudge_x=nudge_factor,color="black",size=5) + 
          labs(title="Sensitivity and Specificty of FHiGR relative to standard GRS",x="False Positive",y="False Negative") + guides(color=FALSE) +
          geom_vline(xintercept=0,linetype="dashed",color="black",alpha=0.5) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.5) +
          scale_size_continuous(name="Cut Point",breaks=scale))
  dev.off()
  
  ##odds ratio for each stratum
  ##write file
  file_name<-paste(sep=".",out,"OR",label_list[l],"txt")
  write.table(format(ORdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
  ##plot
  pdf_fn<-paste(sep=".",out,"OR",label_list[l],"pdf")
  pdf(file=pdf_fn,height=4,width=6)
  print(ggplot(ORdf,aes(x=cutpt,y=OR,color=as.factor(stratum))) + geom_point(alpha=0.7)+   geom_errorbar(aes(ymin=ORdf$LB,ymax=ORdf$UB)) +
    theme_bw() + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) + geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7) +
    labs(title=main,x="Cut Point",y="Odds Ratio"))
  dev.off()
  
  ##odds ratio for stratum=1 top distribution vs everything else
  ##write file
  file_name<-paste(sep=".",out,"stratumOR",label_list[l],"txt")
  write.table(format(stratumORdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
  ##plot
  pdf_fn<-paste(sep=".",out,"stratumOR",label_list[l],"pdf")
  pdf(file=pdf_fn,height=4,width=6)
  print(ggplot(stratumORdf,aes(x=cutpt,y=stratOR)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=stratumORdf$stratLB,ymax=stratumORdf$stratUB)) +
          labs(title=main,x="Cut Point",y="Odds Ratio for Stratum=1 & Top of Distribution Compared to Remaining") +
          geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7))
  dev.off()
  
  ##combine allOR and stratumOR dataframes for comparison
  allORdf$indicator<-as.factor(0)
  stratumORdf$indicator<-as.factor(1)
  names(stratumORdf)<-names(allORdf)
  comparedf<-rbind(stratumORdf,allORdf)
  pdf_fn<-paste(sep=".",out,"compareOR",label_list[l],"pdf")
  pdf(file=pdf_fn,height=4,width=6)
  print(ggplot(comparedf,aes(x=cutpt,y=OR,color=indicator)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=comparedf$LB,ymax=comparedf$UB)) +
          labs(title=main,x="Cut Point",y="Odds Ratio") + scale_color_manual(values=c("darkblue","grey"),name="Conditional on Stratum")+ 
          geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7))
  dev.off()
}


##odds ratio for top vs bottom of distribution for entire dataset
##write file
file_name<-paste(sep=".",out,"allOR.txt")
write.table(format(allORdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")



