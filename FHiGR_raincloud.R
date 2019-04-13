

#https://github.com/RainCloudPlots/RainCloudPlots/tree/master/tutorial_R

#!/usr/bin/Rscript

###########################################################
###################### Libraries ##########################
###########################################################

print(Sys.time())
print(sessionInfo())
library(optparse)
library(data.table)
source("/net/snowwhite/home/bwolford/FHiGR_score/R_rainclouds.R") #loads several other librariesß

###########################################################
################### Read Command Line Parameters ##########
###########################################################

optionList <- list(
  make_option(c("-f","--file"),type="character",help="File with sample IDs, phenotypes, self reported family history, and GRS. Expeects header. White space delimited."),
  make_option(c("-s","--stratum_col"),type="numeric",help="1-based column with stratum (e.g. family history). Must be on binary scale 1/0 with 1 being affirmative. NAs ok."),
  make_option(c("-p","--pheno_col"),type="numeric",help="1-based column with phenotype (e.g. disease status). Must be on binary scale 1/0 with 1 being case. NAs ok."),
  make_option(c("-g","--grs_col"),type="numeric",help="1-based column with GRS, not inverse normalized."),
  make_option(c("-c","--cut_points"),type="character",help="Comma separated list of percentile values in decimal format with which to compare top and bottom of distribution [default=0.8,0.9,0.95,0,0.99,0.995]",default="0.8,0.9,0.95,0.99,.995"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-i","--invNorm"),type="logical",default=FALSE,help="Inverse normalize GRS for entire population [default=FALSE]"),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--xlabel",type="character",default="GRS",help="X-axis label [default='']")
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
cutpts<-as.numeric(strsplit(arguments$options$cut_points,",")[[1]])
out<-arguments$options$output
main<-arguments$options$maintitle
xlab<-arguments$options$xlabel
ylab<-arguments$options$ylabel
dig<-arguments$option$digits
invNorm<-arguments$options$invNorm
header<-arguments$options$header

###########################################################
#################### FUNTIONS #############################
###########################################################

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

## stratum
names(subset)[[strat_col]]<-"group"
names(subset)[[grs_col]]<-"score"
subset$group<-as.factor(subset$group)
xmin<-min(subset$score)
xmax<-max(subset$score)
pdf_fn<-paste(sep=".",out,"stratum_distribution_raincloud.pdf")
pdf(file=pdf_fn,height=6,width=6)
ggplot(subset,aes(x=group,y=score,fill=group,color=group)) + 
  geom_flat_violin(alpha=0.2, position = position_nudge(x = .25, y = 0),adjust =1, trim = FALSE)+
  geom_point(alpha=0.2, position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(group)+0.25, y = score),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  coord_flip(ylim=c(xmin,xmax))+ theme_cowplot() + labs(y=ylab,x=xlab,title=main) + guides(fill = FALSE, colour = FALSE) +
    scale_color_manual(values=c("goldenrod","darkblue")) + scale_fill_manual(values=c("goldenrod","darkblue"))
dev.off()


## stratum + entire population

names(subset)[[strat_col]]<-"group"
names(subset)[[grs_col]]<-"score"
subset$group<-as.factor(subset$group)
all<-subset
all$group<-as.factor(2)
df<-rbind(subset,all)

pdf_fn<-paste(sep=".",out,"stratum_population_distribution_raincloud.pdf")
pdf(file=pdf_fn,height=6,width=6)
ggplot(df,aes(x=group,y=score,fill=group,color=group)) + 
  geom_flat_violin(alpha=0.2, position = position_nudge(x = .20, y = 0),adjust =1, trim = FALSE)+
  geom_point(alpha=0.2, position = position_jitter(width = .15), size = .20)+   
  geom_boxplot(aes(x = as.numeric(group)+0.20, y = score),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  coord_flip()+ theme_cowplot() + labs(y=ylab,x=xlab,title=main) + guides(fill = FALSE, colour = FALSE) +
  scale_color_manual(values=c("goldenrod","darkblue","grey")) + scale_fill_manual(values=c("goldenrod","darkblue","grey")) +
  scale_x_discrete(labels=c("Family History Negative","Family History Positive","Population")) + background_grid(major="xy")
dev.off()

levels(df$group)<-c("Family History Negative","Family History Positive","Population")
pdf_fn<-paste(sep=".",out,"compare_distribution.pdf")
pdf(file=pdf_fn,height=4,width=4)
ggplot(df,aes(x=score,fill=group,color=group)) + geom_density(alpha=0.2) + theme_cowplot() + scale_color_manual(values=c("goldenrod","darkblue","grey")) + 
  scale_fill_manual(values=c("goldenrod","darkblue","grey")) +  labs(y=ylab,x=xlab,title=main) + guides(fill = FALSE, colour = FALSE)
dev.off()

