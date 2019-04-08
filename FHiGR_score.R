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
    make_option(c("-c","--cut_points"),type="character",help="Comma separated list of percentile values with which to compare top and bottom of distribution",default="0.8,0.9,0.95,0.975,0.99"),
    make_option(c("-q","--quantiles"),type="character",help="Comma separated list of q-quantiles for binning the GRS distribution (e.g. 4,10,20 gives quartile, decile, ventile).",default="4,5,10,20")
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


##quantile is how many divisions we should make of dataset
##returns object with prevalences and standard error for distribution split into provided number of quantiles
prev_per_quantile<-function(df,GRS_col,prev_col,qtile){
    if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
        print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
    }
    tiles<-quantile(df[[GRS_col]],seq(from=0,to=1,by=qtile/100)) #quantile values
    index<-seq(qtile,100,by=qtile)
    prevalences<-rep(NA,times=length(tiles)-1) #initialize prevalence list
    ses<-rep(NA,times=length(tiles)-1) #initialize se list
    ns<-rep(NA,times=length(tiles)-1)
    for (i in 1:length(prevalences)) {
        for (r in c(0,1))
            prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] < tiles[i+1]][[prev_col]]
        prevalences[i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
        ns[i]<-length(prev <- list)
        ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/length(prev_list)) #what is SE for this prevalence
    }
    pq<-list(prev=prevalences,se=ses,i=index,tiles=tiles)
    class(pq)<-"prev_quantile_obj"
    return(pq)
}


##quantile is how many divisions we should make of dataset
##returns object with prevalences and standard error for distribution split into provided number of quantiles and by provided stratum
prev_per_quantile_stratum<-function(df,GRS_col,prev_col,strat_col,qtile){
    if (!sum(unique(df[[strat_col]])==c(0,1))==2) {
        print("Stratum column must be a binary variable. Expects 0 and 1.")
    }
    if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
        print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
    }
    tiles<-quantile(df[[GRS_col]],seq(from=0,to=1,by=qtile/100)) #quantile values
    index<-seq(qtile,100,by=qtile)
    prevalences<-matrix(NA,2,length(tiles)-1) #initialize prevalence matrix
    ns<-matrix(NA,2,length(tiles)-1) #initialize count matrix
    ses<-matrix(NA,2,length(tiles)-1)#initialize se matrix

    for (i in 1:length(index)) {
        for (r in c(1,2)){ #iterate over stratum
            prev_list<-df[df[[GRS_col]] > tiles[i] & df[[GRS_col]] < tiles[i+1] & df[[strat_col]]==(r-1)][[prev_col]]
            prevalences[r,i]<-sum(prev_list)/length(prev_list) #how many affected in given quantile
            ns[r,i]<-length(prev_list)
            ses[r,i]<-sqrt((prevalences[r,i]*(1-prevalences[r,i]))/length(prev_list)) #what is SE for this prevalence
        }
    }
    pqs<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
    class(pqs)<-"prev_quantile_stratum_obj"
    return(pqs)
}


###########################################################
########################## MAIN ###########################
###########################################################

##read data 
df<-fread(file,header=T)
print(dim(df))
## To DO: check column assumptions 

##subset to data with stratum available
subset<-df[!is.na(df[[strat_col]])]
print(dim(subset))


#make quantiles then stratify
size<-length(quantiles) #number of quantiles being tested, size of obj
obj<-lapply(quantiles,prev_per_quantile_stratum,df=subset,GRS_col=grs_col,prev_col=pheno_col,strat_col=strat_col)
for (i in 1:size){ #across q-quantiles
    for (j in c(1,2)){ #across stratum
        prev<-obj[[i]]$prev[j,]
        se<-obj[[i]]$se[j,]
        n<-obj[[i]]$n[j,]
        tiles<-obj[[i]]$tiles
        bins<-rep(quantiles[i],length(obj[[i]]$i))

        d<-data.frame(prev=prev,se=se,n=n,tiles=tiles,bins=bins)
        print(d)
             
    }
}

print(obj)

                                        #make plots for each of the bins
                                        #fix language of quantiles bins cut pts etc
#calc ORS at cut poitns 
