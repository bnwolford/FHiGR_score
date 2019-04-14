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
xlabel<-arguments$options$xlabel
ylabel<-arguments$options$ylabel
legend<-arguments$options$legend
dig<-arguments$option$digits
covar<-as.numeric(strsplit(arguments$options$covariates,",")[[1]])
invNorm<-arguments$options$invNorm
header<-arguments$options$header

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
model<-function(df,grs_col,fhigrs_col,strat_col,pheno_col,covar){
  mobj<-list() #initialize object
  ## inverse rank normalize GRS in the population because FHIGRS is inv normal scale
  df$invNormGRS<-rankNorm(df[[grs_col]])
  grs_col<-which(names(df)=="invNormGRS") #new GRS col
  
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  gdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"GRS")) #standard GRS 
  m<-matrix(table(df[[grs_col]],df[[pheno_col]]),byrow=FALSE,nrow=2)
  gdf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  gdf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["GRS"]]<-gdf
  
  
  ## model for FHIGRS
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                              paste(colnames(df)[c(covar,fhigrs_col)], collapse = "+"),
                              sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  m<-matrix(table(df[[fhigrs_col]],df[[pheno_col]]),byrow=FALSE,nrow=2)
  fdf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FHIGRS"))
  fdf$bottom_prev<-m[1,2]/sum(m[1,]) #prevalence in bottom of distribution
  fdf$top_prev<-m[2,2]/sum(m[2,]) #prevalence in top of distribution
  mobj[["FHIGRS"]]<-fdf
    
    
  ## model for family history + GRS, additive
  formula<-as.formula(paste(colnames(df)[pheno_col], "~",
                            paste(colnames(df)[c(covar,strat_col,grs_col)], collapse = "+"),
                            sep = ""))
  glm.obj<-glm(formula=formula,data=df,family="binomial")
  adf<-data.frame(summary(glm.obj)$coefficients,row.names=c("Int",covar,"FH","GRS"))
  adf$bottom_prev<-NA
  adf$top_prev<-NA
  mobj[["additive"]]<-adf
  class(mobj)<-"model_obj"
  return(mobj)
  #TO DO: compare model fit 
}


## screen N people how many do we catch and how many do we miss given screening strategy (prioritize by family history or no)
## expects "top_quantile_stratum_obj"
clinical_impact<-function(df,value,grs_col,fhigrs_col,N=10000){
  
  ## inverse rank normalize GRS in the population because FHIGRS is inv normal scale
  df$invNormGRS<-rankNorm(df[[grs_col]])
  grs_col<-which(names(df)=="invNormGRS") #new GRS col
  counts<-list(rep(NA,2))
  score_list<-c(fhigrs_col,grs_col)
  for (k in c(1,2)){
    q<-quantile(df[[score_list[k]]],value)
    df$dist<-2
    df$dist<-ifelse(df[[score_list[k]]]>q,1,0) #1 if in top, 0 in bottom
    df[df$dist==2]$dist<-NA
    dist_col<-which(names(df)=="dist")
    counts[[k]]<-as.matrix(table(df[[pheno_col]],df[[dist_col]]))
  }
  
  #first row of matrix is counts from bottom of distribution at given cut point
  #second row of matrix is counts from top of distirbution at given cut point
  #first column of matrix is counts of controls
  #second column of matrix is counts of cases 
  
  #FHIGRS
  screen<-c(counts[[2]][2,1],counts[[2]][2,2]) #control,case of top distribution with stratum=1
  no_screen<-c(sum(counts[[2]][1,1],counts[[1]][,1]),sum(counts[[2]][1,2],counts[[1]][,2])) #control,case of top distribution with stratum=0 + bottom distribution stratum=0|1
  m<-matrix(c(no_screen,screen),nrow=2,ncol=2,byrow=TRUE)
  mfrac<-m/sum(m)
  scenario1<-mfrac*N
  
  ##GRS
  all<-counts[[1]]+counts[[2]]
  allfrac<-all/sum(all)
  scenario2<-allfrac*N #no screen is top row, screen is bottom row, control is first column, case is second column
  
  false_pos<-(scenario1-scenario2)[2,1]
  false_neg<-(scenario1-scenario2)[1,2]
  
  ## scenario 1 FHiGRS
  neg_predictive<-m[1,1]/sum(m[1,])
  pos_predictive<-m[2,2]/sum(m[2,])
  specificity<-m[1,1]/sum(m[,1])
  sensitivity<-m[2,2]/sum(m[,2])
  accuracy<-(m[1,1] + m[2,2])/sum(m) 
  
  scenario1_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy,scenario="FHiGR",cutpt=value)
  
  ## scenario 2 GRS
  m<-all
  neg_predictive<-m[1,1]/sum(m[1,])
  pos_predictive<-m[2,2]/sum(m[2,])
  specificity<-m[1,1]/sum(m[,1])
  sensitivity<-m[2,2]/sum(m[,2])
  accuracy<-(m[1,1] + m[2,2])/sum(m)
  
  scenario2_df <- data.frame(false_pos,false_neg, pos_predictive, neg_predictive, sensitivity, specificity, accuracy,scenario="GRS",cutpt=value)
  
  return(rbind(scenario1_df,scenario2_df))
}


###########################################################
#################### MAIN #################################
###########################################################

##read data 
dat<-fread(file,header=header)
#print(dim(df))
## To DO: check column assumptions 

##subset to data with stratum available
subset<-dat[!is.na(dat[[strat_col]])] #remove if NA for stratum
print(dim(subset))

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

  #put q-quantiles in order of prevalence 
  #To do: prevalences from a reference population
  sdf_order<-sdf[order(sdf$prev),]
  sdf_order$rank<-seq(100,(list_length-1)*200,100)

  #assign samples to groups
  #copy subset dataframe so we don't alter original data during for loop over q-quantiles
  qsub<-subset
  qsub$rank<-as.numeric(1) #initialize 
  for (s in c(1,2)){
    sdf_order_strat<-sdf_order[sdf_order$stratum==(s-1),]
    sdf_order_strat[1,'lower_tile']<- -100 #condition to put samples equiv to minimum in bottom bin
    for (j in 1:nrow(sdf_order_strat)){
      u<-sdf_order_strat[j,'upper_tile']
      l<-sdf_order_strat[j,'lower_tile']
      rank<-sdf_order_strat[j,'rank']
      qsub[qsub[[strat_col]]==(s-1) & qsub[[grs_col]]>l & qsub[[grs_col]]<=u,'rank']<-rank
    }
  }
  qsub$grs_rank<-qsub[[grs_col]]+qsub$rank #add rank to GRS
  qsub$FHIGRS<-qnorm((rank(qsub$grs_rank,na.last="keep")-0.5)/sum(!is.na(qsub$grs_rank))) #new FHIGRS

  ##ggplot(subset,aes(x=FHIGRS)) + geom_density()
  
  ### FHiGR dotplot
  qsub[[strat_col]]<-as.factor(qsub[[strat_col]])
  levels(qsub[[strat_col]])<-c("FH-","FH+")
  stratum<-names(qsub)[[strat_col]]
  pdf_fn<-paste(sep=".",out,quantiles[i],"FHiGRS.pdf")
  png_fn<-paste(sep=".",out,quantiles[i],"FHiGRS.png")
  ##make pdf
  pdf(file=pdf_fn,height=5,width=6,useDingbats=FALSE)
  print(ggplot(qsub,aes(x=FHIGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/35,dotsize=0.5) + 
    scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) + coord_cartesian(ylim=c(0,0.4))+
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + labs(title=main,ylab="Density",xlab=xlabel) + theme_bw() )
  dev.off()
  ##make png
  png(file=png_fn,height=1250,width=1500,res=200)
  print(ggplot(qsub,aes(x=FHIGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/35,dotsize=0.5) + 
          scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) + coord_cartesian(ylim=c(0,0.4))+
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + labs(title=main,ylab="Density",xlab=xlabel) + theme_bw())
  dev.off()
  
  ### GRS dotplot
  pdf_fn<-paste(sep=".",out,quantiles[i],"GRS.pdf")
  png_fn<-paste(sep=".",out,quantiles[i],"GRS.png")
  qsub$invNormGRS<-rankNorm(qsub[[grs_col]])
  pdf(file=pdf_fn,height=5,width=6,useDingbats=FALSE)
  print(ggplot(qsub,aes(x=invNormGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/30,dotsize=0.5) + 
          scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + labs(title=main,ylab="Density",xlab=xlabel) + theme_bw())
  dev.off()
  ##make png
  png(file=png_fn,height=1250,width=1500,res=200)
  print(ggplot(qsub,aes(x=invNormGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/35,dotsize=0.5) + 
          scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + labs(title=main,ylab="Density",xlab=xlabel) + theme_bw())
  dev.off()
  
  ## compare GRS and FHIGRS with logistic regression
  fhigrs_col<-which(names(qsub)=="FHIGRS")
  qsub[[strat_col]]<-as.numeric(qsub[[strat_col]])-1 #turn strat back to value
  mobj<-model(df=qsub,grs_col=grs_col,fhigrs_col=fhigrs_col,pheno_col=pheno_col,strat_col=strat_col,covar=covar)
  print(mobj)
  rowname_list<-c("GRS","FHIGRS","GRS")
  score_list<-c("GRS","FHIGRS","GRS+FH")
  for (l in c(1,2,3)){
    if (i==1 & l==1){
    d<-data.frame(OR=exp(mobj[[l]][rowname_list[l],'Estimate']),
                  LB=exp(mobj[[l]][rowname_list[l],'Estimate']-1.96*mobj[[l]][rowname_list[l],'Std..Error']),
                  UB=exp(mobj[[l]][rowname_list[l],'Estimate']+1.96*mobj[[l]][rowname_list[l],'Std..Error']),
                  pval=mobj[[l]][score_list[l],'Pr...z..'],
                  top_prev<-unique(mobj[[l]]['top_prev']),
                  bottom_prev<-unique(mobj[[l]]['bottom_prev']),
                  score=score_list[l],
                  qtile=quantiles[i])
    } else {
      d<-rbind(d,data.frame(OR=exp(mobj[[l]][rowname_list[l],'Estimate']),
                            LB=exp(mobj[[l]][rowname_list[l],'Estimate']-1.96*mobj[[l]][rowname_list[l],'Std..Error']),
                            UB=exp(mobj[[l]][rowname_list[l],'Estimate']+1.96*mobj[[l]][rowname_list[l],'Std..Error']),
                            pval=mobj[[l]][score_list[l],'Pr...z..'],
                            top_prev<-unique(mobj[[l]]['top_prev']),
                            bottom_prev<-unique(mobj[[l]]['bottom_prev']),
                            score=score_list[l],
                            qtile=quantiles[i]))
      
    }
  }
  
  
  ## quantify clinical impact
  obj<-lapply(cutpts,clinical_impact,df=qsub,grs_col=grs_col,fhigrs_col=fhigrs_col)
  clin_df<-bind_rows(obj)
  file_n<-paste(sep=".",out,"clinicalImpact.txt")
  write.table(format(clin_df,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
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
#write table 
file_n<-paste(sep=".",out,"compareScores.txt")
write.table(format(d,digits=dig),file=file_n,quote=FALSE,row.names=FALSE,sep="\t")
##plot comparison
pdf_fn<-paste(sep=".",out,"compare.pdf")
pdf(file=pdf_fn,height=4,width=6,useDingbats=FALSE)
print(ggplot(d,aes(x=qtile,y=OR,color=score)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=d$LB,ymax=d$UB)) + scale_x_continuous(breaks=quantiles) + 
        labs(title=main,x="Number of Quantiles in which Prevalence is Estimated",y="Odds Ratio for Score") + scale_color_manual(values=c("grey","darkblue","orchid4"),name="") +
        geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7))
dev.off()





