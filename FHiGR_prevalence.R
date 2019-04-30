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
  make_option(c("-q","--quantiles"),type="character",help="Comma separated list of q-quantiles for binning the GRS distribution (e.g. 4, 10,20 gives quartile, decile, ventile) [default=4,5,10,20,100]",default="4,5,10,20,100"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--xlabel",type="character",default="GRS",help="X-axis label [default='GRS']"),
  make_option("--ylabel",type="character",default="Prevalence",help="Y-axis label [default='Prevalence']"),
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
out<-arguments$options$output
main<-arguments$options$maintitle
xlabel<-arguments$options$xlabel
ylabel<-arguments$options$ylabel
legend<-arguments$options$legend
dig<-arguments$option$digits
header<-arguments$options$header
##source relevant code from code base
source(paste0(arguments$options$codeDir,"helperFunctions.R")) ##will be used to calculate FHiGRS



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

  ##fix any NaN which arise due to no cases in a bin
  prevalences[is.nan(prevalences)]<-0
  ses[is.nan(ses)]<-0
  ns[is.nan(ns)]<-0
  ##create object 
  pqs<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
  class(pqs)<-"prev_quantile_stratum_obj"
  return(pqs)
}

##quantile is how many divisions we should make of dataset (e.g. 20-quantile is ventiles for 20 groups)
##returns object with prevalences and standard error for distribution
##does not stratify
prev_per_quantile<-function(df,GRS_col,prev_col,qtile){
    if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
        print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
    }
    if (sum(qtile)<2*length(qtile)){ #check qtile
        print("q-quantiles should be number of divisions for data set and must be greater than 1")
    }
    ##initialize data structures
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
    ##create object
    pq<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
    class(pq)<-"prev_quantile_obj"
    return(pq)
}



## plot prevalence per GRS quantile bin, color by stratum if stratum==TRUE
plotting<-function(dat,out,qtiles,stratum=FALSE,main,xlab,ylab,legend,ymax=1){
    #dat has columns  prev          se    n    tiles  q stratum percents lower_tile upper_tile
    dat$frac<-as.numeric(sub("%","",dat$percents)) #convert factor percentages to numeric
    dat<-dat[dat$frac!=1.00,]
    dat$ub<-dat$prev+(1.96*dat$se)
    dat$lb<-dat$prev-(1.96*dat$se)

    if (ymax==1){ # if ymax is not given to function, just plot to scale of data, otherwise, script can be given custom ymax to match scale of another plot
        ymax=max(dat$prev)
    }

    if (stratum==TRUE){ #stratify
        dat$stratum<-as.factor(dat$stratum)
        levels(dat$stratum)<-c("Negative","Positive") #change labels from 0/1
        by(dat, dat$q, #number of q-quantiles (e.g. break data into 4 bins, 10 bins, etc.)
           function (x) {
               name=unique(x$q)
               if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
               pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5,useDingbats=FALSE)
               print(ggplot(x,aes(x=frac,y=prev,color=as.factor(stratum))) + geom_point() + theme_bw() + geom_errorbar(aes(ymin=x$lb,ymax=x$ub)) +
                     scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +labs(title=main) + xlab(xlab) + ylab(ylab)  +
                     scale_x_continuous(breaks=breaks) + coord_cartesian(ylim=c(0,ymax)))
               dev.off()
           }
        )
    } else { #all data,  plotting hacks to keep plot size the same as stratum
        by(dat, dat$q, #number of q-quantiles (e.g. break data into 4 bins, 10 bins, etc.)
           function (x) {
               name=unique(x$q)
               if (unique(x$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=x$frac}
               pdf(file=paste(sep=".",out,name,"pdf"),height=5,width=5,useDingbats=FALSE)
               print(ggplot(x,aes(x=frac,y=prev,color=as.factor(1))) + geom_point() +
                     scale_color_manual(values=c("grey"),guide=guide_legend(override.aes=list(color="white")),name=legend) +
                     geom_errorbar(aes(ymin=x$lb,ymax=x$ub),color="grey")  +
                     theme_bw() + labs(title=main) + xlab(xlab) + ylab(ylab) + scale_x_continuous(breaks=breaks) +
                     coord_cartesian(ylim=c(0,ymax)) +
                     theme(legend.text=element_text(color = "white"), legend.title = element_text(color = "white"), legend.key = element_rect(fill = "white")))
               dev.off()
           }
       )
    }
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

############# Prevalences versus GRS#############
## uses prev_per_quantile_stratum and prev_per_quantile functions

##make quantiles then stratify
size<-length(quantiles) #number of quantiles being tested, size of obj
qobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=TRUE)
for (i in 1:size){ #across q-quantiles
    for (j in c(1,2)){ #across stratum
        list_length<-length(qobj[[i]]$prev[j,]) #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
        prev<-qobj[[i]]$prev[j,][-list_length]
        se<-qobj[[i]]$se[j,][-list_length]
        n<-qobj[[i]]$n[j,][-list_length]
        tiles<-qobj[[i]]$tiles[-1]
        lower_tile<-qobj[[i]]$tiles[1:list_length-1]
        upper_tile<-qobj[[i]]$tiles[2:list_length]
        percents<-names(tiles)
        bins<-rep(quantiles[i],list_length-1)
        strat<-rep(j-1,list_length-1)
        if (i==1 & j==1) {
            qdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
        } else {
            qdf<-rbind(qdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
        }
    }
}

##stratify then calculate quantiles
sobj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col,qfirst=FALSE)
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
        if (i==1 & j==1) {
            sdf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL)
        } else {
            sdf<-rbind(sdf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,lower_tile=lower_tile,upper_tile=upper_tile,row.names=NULL))
        }
    }
}

ymax<-max(max(qdf$prev+(1.96*qdf$se)),max(sdf$prev+(1.96*sdf$se)))
file_name<-paste(sep=".",out,"quantileFirst.txt")
write.table(format(qdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(qdf,paste(sep="_",out,"quantileFirst"),quantiles,TRUE,main,xlabel,ylabel,legend,ymax)

file_name<-paste(sep=".",out,"stratifyFirst.txt")
write.table(format(sdf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(sdf,paste(sep="_",out,"stratifyFirst"),quantiles,TRUE,main,xlabel,ylabel,legend,ymax)

##calculate quantiles for all data, no stratification
aobj<-lapply(quantiles,prev_per_quantile,df=subset,GRS_col=grs_col,prev_col=pheno_col)
for (i in 1:size){ #across q-quantiles
    list_length<-length(aobj[[i]]$prev)  #need to remove the 0th percentile so the prevalence aligns with correct nth percentile
    prev<-aobj[[i]]$prev[-list_length]
    se<-aobj[[i]]$se[-list_length]
    n<-aobj[[i]]$n[-list_length]
    tiles<-aobj[[i]]$tiles[-1]
    percents<-names(tiles)
    bins<-rep(quantiles[i],list_length-1)
    strat<-rep(NA,list_length-1)
    if (i==1){
        adf<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
    } else {
        adf<-rbind(adf,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
    }
}
file_name<-paste(sep=".",out,"noStratify.txt")
write.table(format(adf,digits=dig),file=file_name,quote=FALSE,row.names=FALSE,sep="\t")
plotting(adf,paste(sep="_",out,"all"),quantiles,FALSE,main,xlabel,ylabel,legend,ymax)

################## Estimate FHiGRS and plot prevalence per quantile  ###################

for (q in 1:length(quantiles)){
    fdf<-estimate_FHiGRS(sdf[sdf$q==quantiles[q],],subset,strat_col,grs_col) #use sdf from prev_quantile_per_stratum object to estimate FHiGRS
    fhigrs_col<-which(names(fdf)=="FHIGRS")
    
    logic_list<-c("TRUE","FALSE")
    label_list<-c("quantileFirst","stratifyFirst")
    #recalculate prevalences per quantile bin this time using FHiGRS instead of GRS
    for (i in c(1,2)){ #across stratify logic
        for (j in c(1,2)) { #across stratum
            fobj<-prev_per_quantile_stratum(qtile=quantiles[q],df=fdf,GRS_col=fhigrs_col,strat_col=strat_col,prev_col=pheno_col,qfirst=logic_list[i])
            list_length<-length(fobj$prev[j,])
            prev<-fobj$prev[j,][-list_length]
            se<-fobj$se[j,][-list_length]
            n<-fobj$n[j,][-list_length]
            tiles<-fobj$tiles[-1]
            percents<-names(tiles)
            bins<-rep(quantiles[q],list_length-1)
            strat<-rep(j-1,list_length-1)
            if (j==1) {
                fobj_df<-data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL)
            } else {
                fobj_df<-rbind(fobj_df,data.frame(prev=prev,se=se,n=n,tiles=tiles,q=bins,stratum=strat,percents=percents,row.names=NULL))
            }
        }

        ##make prevalence plot 
        name<-as.character(unique(fobj_df$q))
        fobj_df$frac<-as.numeric(sub("%","",fobj_df$percents)) #convert factor percentages to numeric
        fobj_df<-fobj_df[fobj_df$frac!=1.00,]
        fobj_df$ub<-fobj_df$prev+(1.96*fobj_df$se)
        fobj_df$lb<-fobj_df$prev-(1.96*fobj_df$se)
        ymax<-max(fobj_df$prev) #this may need to be changed
        fobj_df$stratum<-as.factor(fobj_df$stratum)
        levels(fobj_df$stratum)<-c("Negative","Positive")
        if (unique(fobj_df$q) > 10) {breaks=c(0,10,20,30,40,50,60,70,80,90,100)} else {breaks=fobj_df$frac}
        pdf(file=paste(sep=".",out,"FHiGRS",label_list[i],name,"pdf"),height=5,width=5,useDingbats=FALSE)
        print(ggplot(fobj_df,aes(x=frac,y=prev,color=stratum)) + geom_point() + theme_bw() +
            geom_errorbar(aes(ymin=fobj_df$lb,ymax=fobj_df$ub)) +
            scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
            labs(title=main) + xlab("FHiGR Score") + ylab(ylabel)  +
            scale_x_continuous(breaks=breaks) + coord_cartesian(ylim=c(0,ymax)))
        dev.off()
    }
    ##make dot plot 
    fdf[[strat_col]]<-as.factor(fdf[[strat_col]])
    levels(fdf[[strat_col]])<-c("Negative","Positive") #change labels from 0/1
    stratum<-names(fdf)[[strat_col]]
    pdf_fn<-paste(sep=".",out,quantiles[q],"FHiGRS.dotplot.pdf")
    png_fn<-paste(sep=".",out,quantiles[q],"FHiGRS.dotplot.png")
    ##make pdf
    pdf(file=pdf_fn,height=5,width=6,useDingbats=FALSE)
    print(ggplot(fdf,aes(x=FHIGRS,color=get(stratum),fill=get(stratum)))  +
          geom_dotplot(method="histodot",binwidth=1/38,dotsize=0.5) +
          scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          labs(title=main,ylab="Density",xlab="FHiGR Score") + theme_bw() + scale_y_continuous(NULL, breaks = NULL))
    dev.off()
    ##make png
    png(file=png_fn,height=1000,width=1200,res=200)
    print(ggplot(fdf,aes(x=FHIGRS,color=get(stratum),fill=get(stratum)))  +
          geom_dotplot(method="histodot",binwidth=1/38,dotsize=0.5) +
          scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          labs(title=main,ylab="Density",xlab="FHiGR Score") + theme_bw()  + scale_y_continuous(NULL, breaks = NULL))
    dev.off()
}

### GRS dotplot
pdf_fn<-paste(sep=".",out,"GRS.dotplot.pdf")
png_fn<-paste(sep=".",out,"GRS.dotplot.png")
subset$invNormGRS<-rankNorm(subset[[grs_col]])
subset[[strat_col]]<-as.factor(subset[[strat_col]])
levels(subset[[strat_col]])<-c("Negative","Positive") #change labels from 0/1
stratum<-names(subset)[[strat_col]]
##make pdf
pdf(file=pdf_fn,height=5,width=6,useDingbats=FALSE)
print(ggplot(subset,aes(x=invNormGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/33,dotsize=0.5) +
      scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      labs(title=main,ylab="Density",xlab=xlabel) + theme_bw()  + scale_y_continuous(NULL, breaks = NULL))
dev.off()
##make png
png(file=png_fn,height=1250,width=1500,res=200)
print(ggplot(subset,aes(x=invNormGRS,color=get(stratum),fill=get(stratum)))  +  geom_dotplot(method="histodot",binwidth=1/25,dotsize=0.5) +
      scale_fill_manual(values=c("goldenrod3","darkblue"),name=legend) + scale_color_manual(values=c("goldenrod3","darkblue"),name=legend) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      labs(title=main,ylab="Density",xlab=xlabel) + theme_bw()  + scale_y_continuous(NULL, breaks = NULL))
dev.off()


##To do: fix size of dotplots, weird function of binwidth 
