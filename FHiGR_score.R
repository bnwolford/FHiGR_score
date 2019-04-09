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
    make_option(c("-c","--cut_points"),type="character",help="Comma separated list of percentile values in decimal format with which to compare top and bottom of distribution",default="0.8,0.9,0.95,0.975,0.99"),
    make_option(c("-q","--quantiles"),type="character",help="Comma separated list of q-quantiles for binning the GRS distribution (e.g. 4,10,20 gives quartile, decile, ventile).",default="4,5,10,20,100"),
    make_option(c("-o","--output"),type="character",help="Prefix for output files",default="FHiGR"),
    make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables",default=3),
    make_option("--maintitle", type="character", default="",help="Plot title"),
    make_option("--xlabel",type="character",default="",help="X-axis label"),
    make_option("--ylabel",type="character",default="",help="Y-axis label"),
    make_option("--legend",type="character",default="",help="Legend title which is stratum")
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
#split into provided number of quantiles and then by provided stratum if qfirst==TRUE
prev_per_quantile_stratum<-function(df,GRS_col,prev_col,strat_col,qtile,qfirst=FALSE){
    if (!sum(unique(df[[strat_col]])==c(0,1))==2) {
        print("Stratum column must be a binary variable. Expects 0 and 1.")
    }
    if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
        print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
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
plotting<-function(dat,out,qtiles,stratum=FALSE,main,xlab,ylab,legend){
  #dat<-dat[complete.cases(dat),] #remove the last bins which are NA
  dat$frac<-as.numeric(sub("%","",dat$percents)) #convert factor percentages to numeric 
  dat<-dat[dat$frac!=1.00,] 
  dat$ub<-dat$prev+(1.96*dat$se)
  dat$lb<-dat$prev-1.96*dat$se
  if (stratum==TRUE){
    by(dat, dat$q, 
      function (x) {
        name=unique(x$q)
        if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
        pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5)
        print(ggplot(x,aes(x=frac,y=prev,color=as.factor(stratum))) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=x$lb,ymax=x$ub)) + 
          scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +labs(title=main) + xlab(xlab) + ylab(ylab)  + 
          scale_x_continuous(breaks=breaks)) 
        dev.off()
      }
    )
  } else {
    by(dat, dat$q, 
       function (x) {
         name=unique(x$q)
         if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
         pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5)
         print(ggplot(x,aes(x=frac,y=prev)) + geom_point(color="grey") + geom_errorbar(aes(ymin=x$lb,ymax=x$ub),color="grey")  + 
                 theme_bw() + labs(title=main) + xlab(xlab) + ylab(ylab) + scale_x_continuous(breaks=breaks))
         dev.off()
       }
    )
  }
  return(0)
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
  for (r in c(1,2)) {
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
  tq<-list(OR=OR,SE=SE,LB=LB,UB=UB)
  class(tq)<-"top_quantile_obj"
  return(tq)
}
  

###########################################################
########################## MAIN ###########################
###########################################################

##read data 
df<-fread(file,header=T)
#print(dim(df))
## To DO: check column assumptions 

##subset to data with stratum available
subset<-df[!is.na(df[[strat_col]])]
print(dim(subset))
size<-length(quantiles) #number of quantiles being tested, size of obj

############# Prevalences versus GRS

##make quantiles then stratify
obj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=TRUE)
#print(obj)
for (i in 1:size){ #across q-quantiles
    for (j in c(1,2)){ #across stratum
        list_length<-length(obj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
        prev<-obj[[i]]$prev[j,][-list_length]
        se<-obj[[i]]$se[j,][-list_length]
        n<-obj[[i]]$n[j,][-list_length]
        tiles<-obj[[i]]$tiles[-1]
        percents<-names(tiles)
        bins<-rep(quantiles[i],list_length-1)
        strat<-rep(j-1,list_length-1)
        if (i==1 & j==1) {
          df<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
        } else {
          df<-rbind(df,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
        }
    }
}
file_name<-paste(sep=".",out,"quantileFirst.txt")
write.table(format(df,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(df,paste(sep="_",out,"quantileFirst"),quantiles,TRUE,main,xlab,ylab,legend)

##stratify then calculate quantiles
obj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=FALSE)
#print(obj)
for (i in 1:size){ #across q-quantiles
  for (j in c(1,2)){ #across stratum
    list_length<-length(obj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-obj[[i]]$prev[j,][-list_length]
    se<-obj[[i]]$se[j,][-list_length]
    n<-obj[[i]]$n[j,][-list_length]
    tiles<-obj[[i]]$tiles[-1]
    percents<-names(tiles)
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(j-1,list_length-1)
    if (i==1 & j==1) {
      df<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
    } else {
      df<-rbind(df,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
    }
  }
}
file_name<-paste(sep=".",out,"stratifyFirst.txt")
write.table(format(df,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(df,paste(sep="_",out,"stratifyFirst"),quantiles,TRUE,main,xlab,ylab,legend)

##calculate quantiles for all data, no stratification
all_obj<-lapply(quantiles,prev_per_quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col)
print(all_obj)
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
plotting(all_df,paste(sep="_",out,"all"),quantiles,FALSE,main,xlab,ylab,legend)  


################## Odds Ratios 

#cutpts<-seq(0.8,1,by=.01)
if (1 %in% cutpts){ #code doesn't work with 
  cutpts<-cutpts[-which(cutpts==1)]
}
n<-length(cutpts)
logical_list<-c(TRUE,FALSE)
label_list<-c("quantileFirst","stratifyFirst")
for (l in c(1,2)){
  obj<-lapply(cutpts,top_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=logical_list[l])
  allobj<-lapply(cutpts,top_quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col)
  print(allobj)
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
  print(pdf_fn)
  pdf(file=pdf_fn,height=6,width=8)
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
  pdf(file=pdf_fn,height=6,width=8)
  print(ggplot(comparedf,aes(x=cutpt,y=OR,color=indicator)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=comparedf$LB,ymax=comparedf$UB)) +
          labs(title=main,x="Cut Point",y="Odds Ratio") + scale_color_manual(values=c("darkblue","grey"),name="Conditional on Stratum")+ 
          geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7))
  dev.off()
}



##odds ratio for top vs bottom of distribution for entire dataset
##write file
file_name<-paste(sep=".",out,"allOR.txt")
write.table(allORdf,file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
##plot
#pdf_fn<-paste(sep=".","allOR.pdf")
#pdf(file=pdf_fn,height=6,width=8)
#print(ggplot(allORdf,aes(x=cutpt,y=OR)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=allORdf$LB,ymax=allORdf$UB)) +
        #labs(title=main,x="Cut Point",y="Odds Ratio for Stratum=1 & Top of Distribution Compared to Remaining") + 
        #geom_hline(linetype="dashed",color="black",yintercept=1,alpha=0.7))
#dev.off()



