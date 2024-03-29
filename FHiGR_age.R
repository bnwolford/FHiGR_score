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
library(patchwork)
library(LaCroixColoR)

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
  make_option(c("-a","--age_col"),type="numeric",help="1-based column with participation or enrollment age e.g. age of self reported family history"),
  make_option(c("-b","--birthYear_col"),type="numeric",help="1-based column with birthyear"),
  make_option(c("-o","--output"),type="character",help="Prefix for output files [defualt=FHiGR]",default="FHiGR"),
  make_option(c("-d","--digits"),type="numeric",help="Number of decimal digits to print in tables [default=3]",default=3),
  make_option(c("-r","--header"),type="logical",default=FALSE,help="If phenotype file has a header [default=FALSE]"),
  make_option("--maintitle", type="character", default="",help="Plot title [default='']"),
  make_option("--legend",type="character",default="",help="Legend"),
  make_option("--codeDir",type="character",default="/FHiGR_score/",help="Directory for repository for sourcing other code in code base [default=/FHiGR_score/]")
)

parser <- OptionParser(
  usage="%prog --file --stratum_col --pheno_col --grs_col --age_col --output --digits --header --maintitle --codeDir",
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
  warning("Stratum col is required -s or --stratum_col")
}
pheno_col<-arguments$options$pheno_col
if (length(pheno_col)==0){
  warning("Pheno col is required -p or --pheno_col")
}
grs_col<-arguments$options$grs_col
if (length(grs_col)==0){
  warning("GRS col is required -g or --grs_col")
}
age_col<-arguments$options$age_col
if (length(age_col)==0){
    warning("Age col is required -a or --age_col")
}
birthYear_col<-arguments$options$birthYear_col
if (length(birthYear_col)==0){
    warning("Birthyear col is required -b or --birthYear_col")
}
out<-arguments$options$output
main<-arguments$options$maintitle
legend<-arguments$options$legend
dig<-arguments$option$digits
header<-arguments$options$header

#colors
pamp<-lacroix_palette("Pamplemousse",type = "discrete", n=6) #6 elements
pamp2<-c(pamp[1:6],"#8B4789") #7 elements 
passion<-lacroix_palette("PassionFruit",type="discrete",n=6) #6 elements


###########################################################
#################### FUNCTIONS #############################
###########################################################

#calculates 95% CI from a proportion
CI<-function(prop,n){
    se<-sqrt((prop*(1-prop))/n) #standard error
    ub<-prop+(1.96*se)
    lb<-prop-(1.96*se)
    return(c(n,se,lb,ub))
}


###########################################################
#################### MAIN #################################
###########################################################

##read data 
dat<-fread(file,header=header)
print(paste("Data dimensions are:",dim(dat)[1],dim(dat)[2]))


##make suer strat and pheno columns are integers
dat[[strat_col]]<-as.integer(dat[[strat_col]])
dat[[pheno_col]]<-as.integer(dat[[pheno_col]])


#print(table(dat[[strat_col]],useNA="always"))
n<-nrow(dat)
strat0<-nrow(dat[dat[[strat_col]]==0,])
strat1<-nrow(dat[dat[[strat_col]]==1,])
stratNA<-nrow(dat[is.na(dat[[strat_col]]),])

## age at time of self report for everyone                                      
pdf_fn<-paste(sep=".",out,"age_all.pdf")
pdf(file=pdf_fn,height=6,width=6,useDingbats=FALSE)
ggplot(dat,aes(x=get(names(dat)[age_col]))) + geom_density(fill="black",alpha=0.5) + theme_bw() +
    labs(x="Enrollment age", caption=bquote("N="~.(n)))+
    theme(plot.caption=element_text(hjust=0.5),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()


dat$strat<-as.factor(dat[[strat_col]]) #make factor

##dynamically put the smaller sample size as the top color 
if (strat1>strat0){ #plot bigger sample on bottom (first)
    dat$strat<-relevel(dat$strat,"1")
    colors<-c("dark blue","goldenrod")
    labels<-c("Positive","Negative")
} else {
    dat$strat<-relevel(dat$strat,"0")
    colors<-c("goldenrod","dark blue")
    labels<-c("Negative","Positive")
}
##age at time of self report and family history for 1/0
pdf_fn<-paste(sep=".",out,"age_stratify.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[!is.na(dat$strat)],aes(x=get(names(dat)[age_col]),fill=strat)) + geom_density(alpha=0.7) + theme_bw() +
    scale_fill_manual(values=colors,name=legend,labels=labels) +
    labs(x="Enrollment Age", caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)),title=main)+
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

## birth year distribution stratified
pdf_fn<-paste(sep=".",out,"birthYear_stratify.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[!is.na(dat$strat)],aes(x=get(names(dat)[birthYear_col]),fill=strat)) + geom_density(alpha=0.7) + theme_bw() + scale_fill_manual(values=colors,name=legend,labels=labels) +
    labs(x="BirthYear",caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)),title=main)+
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

##age at time of self report and family history for NA (age may also be missing if FH is missing)
pdf_fn<-paste(sep=".",out,"age_stratifyNA.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[is.na(dat[[strat_col]])],aes(x=get(names(dat)[age_col]))) + geom_density(alpha=0.7,fill="grey") + theme_bw() +
    labs(x="Enrollment Age",caption=bquote(N['NA']~"="~.(stratNA)),tile=main) +
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()


## birth yaer distirbution for samples with NA for stratum
pdf_fn<-paste(sep=".",out,"birthYear_stratifyNA.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat[is.na(dat[[strat_col]])],aes(x=get(names(dat)[birthYear_col]))) + geom_density(alpha=0.7,fill="grey") + theme_bw() +
    labs(x="Birth Year",caption=bquote(N['NA']~"="~.(stratNA)),title=main) + 
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10))
dev.off()

## age at time of self report for 0/1/NA
dat$facet<-ifelse(is.na(dat[[strat_col]]),1,0) #make indicator variable for NA or not
pdf_fn<-paste(sep=".",out,"age_hist.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat,aes(x=get(names(dat)[age_col]),fill=strat)) + geom_histogram(alpha=0.7,binwidth=1,position="identity")  + theme_bw() +
    labs(x="Enrollment Age", caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)~","~N['NA']~"="~.(stratNA)),title=main)+ facet_wrap(~facet) +
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10),strip.text=element_blank()) +
    scale_fill_manual(na.value="grey",values=colors,name=legend,labels=labels)
dev.off()

##proportion of positive strata per age bin 
dat$bin<-cut(dat[[age_col]],breaks=c(15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95)) #age bin breaks
age_bins<-unique(dat$bin)
age_df<-data.frame(matrix(NA,nrow=length(age_bins),ncol=32))
for (i in 1:length(age_bins)){
    sub<-dat[dat$bin==age_bins[i],] #subset to age
    if (nrow(sub) > 30) { #if more than 30 samples in the bin we can analyze it
        prop<-nrow(sub[sub[[strat_col]]==1,])/nrow(sub) #proportion of sample  in bin with strata 1 with 0/1/NA as denominator
        dat_woNA<-dat[!is.na(dat[[strat_col]]),]
        sub_woNA<-sub[!is.na(sub[[strat_col]]),]
        n_NA<-nrow(sub[is.na(sub[[strat_col]]),]) #number of NA samples in bin
        prop_woNA<-nrow(sub_woNA[sub_woNA[[strat_col]]==1,])/nrow(sub_woNA) #proportion of sample in bin with strata 1 using 0/1 as denominator
        prop_strat1<-nrow(sub[sub[[strat_col]]==1,])/nrow(dat[dat[[strat_col]]==1,]) #proportion of total stratum 1 that is in this age bin
        prop_strat1_woNA<-nrow(sub_woNA[sub_woNA[[strat_col]]==1,])/nrow(dat_woNA[dat_woNA[[strat_col]]==1,]) #proportion of total stratum 1 that is in this age bin (excluding NA from denominator) 
        prop_sample<-nrow(sub)/nrow(dat) #proportion of sample in bin  
        prop_sample_woNA<-nrow(sub_woNA)/nrow(dat_woNA) #proportion of sample (without NA) in bin
        ##CI function gives n, se, lb, ub
        age_df[i,]<-c(prop,CI(prop,nrow(sub)),
                      prop_woNA,CI(prop_woNA,nrow(sub_woNA)),
                      prop_strat1,CI(prop_strat1,nrow(dat)), #not sure about N for CI here
                      prop_strat1_woNA,CI(prop_strat1_woNA,nrow(dat_woNA)), #not sure about N for CI here
                      prop_sample,CI(prop_sample,nrow(dat)),
                      prop_sample_woNA,CI(prop_sample_woNA,nrow(dat_woNA)),
                      n_NA,
                      as.character(age_bins[i])) 

    } else { #empty bin
        age_df[i,]<-c(rep(NA,31),as.character(age_bins[i]))
    }
}
names(age_df)<-c("prop","n","se","lb","ub",
                 "prop_woNA","n_woNA","se_woNA","lb_woNA","ub_woNA",
                 "prop_strat1","n_strat1","se_strat1","lb_strat1","ub_strat1",
                 "prop_strat1_woNA","n_strat1_woNA","se_strat1_woNA","lb_strat1_woNA","ub_strat1_woNA",
                 "prop_sample","n_sample","se_sample","lb_sample","ub_sample",
                 "prop_sample_woNA","n_sample_woNA","se_sample_woNA","lb_sample_woNA","ub_sample_woNA",
                 "n_NA","age_bin")
age_df[c(1:31)]<-sapply(age_df[c(1:31)],as.numeric) #make values numeric
print(age_df)

##plot
pdf_fn<-paste(sep=".",out,"prop_strata_by_age.pdf")
pdf(file=pdf_fn,height=5,width=8,useDingbats=FALSE)
p1<-ggplot(age_df,aes(x=age_bin,y=prop)) + geom_point(aes(size=n)) + theme_bw() + geom_errorbar(aes(ymin=age_df$lb,ymax=age_df$ub)) +
    xlab("Enrollment Age Bin") + ylab("Prevalence of Positive Family History\n within bin") +
    theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + scale_size_continuous(name="N") + coord_cartesian(ylim=c(0,0.6)) +
     guides(size= guide_legend(nrow = 2))
p2<-ggplot(age_df,aes(x=age_bin,y=prop_woNA)) + geom_point(aes(size=n_woNA)) + theme_bw() + geom_errorbar(aes(ymin=age_df$lb_woNA,ymax=age_df$ub_woNA)) +
    xlab("Enrollment Age Bin") + ylab("Prevalence of Positive Family History\n within bin without NA samples") + 
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_size_continuous(name="N") + coord_cartesian(ylim=c(0,0.6)) +
    guides(size= guide_legend(nrow = 2))
p1+p2+plot_annotation(title=main)
dev.off()

tmp1<-age_df[,c("prop_woNA","lb_woNA","ub_woNA","age_bin")]
names(tmp1)<-c("proportion","lb","ub","age_bin")
tmp1$label<-"+FH/bin"
tmp2<-age_df[,c("prop_sample_woNA","lb_sample_woNA","ub_sample_woNA","age_bin")]
names(tmp2)<-c("proportion","lb","ub","age_bin")
tmp2$label<-"bin/sample"
prop_df<-rbind(tmp1,tmp2)
prop_df<-prop_df[complete.cases(prop_df),]

pdf_fn<-paste(sep=".",out,"prop_strata_by_age_prop_in_bin.pdf")
pdf(file=pdf_fn,height=5,width=5,useDingbats=FALSE)
ggplot(prop_df,aes(x=age_bin,y=proportion,color=label)) + geom_point() + theme_bw() +
  theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(ymin=prop_df$lb,ymax=prop_df$ub)) + scale_color_manual(name="Proportion",values=passion) +
  coord_cartesian(ylim=c(0,0.75)) + labs(x="Enrollment Age Bin",title=main)
dev.off()

#dual axis plot with bar for count of NA samples
age_df<-age_df[order(age_df$age_bin),]
age_df<-age_df[complete.cases(age_df),]
age_df$group<-1:nrow(age_df)
pdf_fn<-paste(sep=".",out,"prop_strata_by_age_NA_bar.pdf")
pdf(file=pdf_fn,height=5,width=5,useDingbats=FALSE)
par(mar=c(5.1,5.1,2.1,5.1))
bp<-barplot(age_df$n_NA,yaxt="n",ylim=c(0,max(age_df$n_NA+1000)),xlab="")
axis(side=1,at=bp,labels=age_df$age_bin,las=3,xlab="")
mtext("Age Bin",side=1,line=3.5)
axis(side=4)
mtext("Number of NA samples",side=4,line=3)
par(new=TRUE)
plot(bp,age_df$prop_woNA,cex=0.8,pch=19,xaxt='n',ylab="Proportion of Positive Family History within Bin",ylim=c(0,1),xlab="",xlim=c(min(bp)-0.5,max(bp)+0.5))
arrows(bp, age_df$lb_woNA, bp, age_df$ub_woNA, length=0.05, angle=90, code=3)
title(main)
dev.off()


##cumulative proportion 
age_df$cumprop<-cumsum(age_df$prop_strat1)
age_df$cumprop_woNA<-cumsum(age_df$prop_strat1_woNA)
#print(age_df)
pdf_fn<-paste(sep=".",out,"cumulprop_strata_by_age.pdf")
pdf(file=pdf_fn,height=5,width=8,useDingbats=FALSE)
p3<-ggplot(age_df,aes(x=age_bin,y=cumprop))  + geom_point(aes(size=n)) + theme_bw() +
   xlab("Enrollment Age Bin") + ylab("Cumulative Proportion of Positive Family History\n within age bin") +
   theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + scale_size_continuous(name="N") +
   guides(size= guide_legend(nrow = 2)) + geom_line() 
p4<-ggplot(age_df,aes(x=age_bin,y=cumprop_woNA)) + geom_point(aes(size=n_woNA)) +  theme_bw() +
    xlab("Enrollment Age Bin") + ylab("Cumulative Proportion of Positive Family History\n within age bin without NA samples") + 
    theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + scale_size_continuous(name="N") +
    guides(size= guide_legend(nrow = 2)) + geom_line()
p3+p4 +plot_annotation(title=main)
dev.off()

######### CURRENT AGE ###########
## estimated current age for 0/1/NA (NA facet is not on the same scale though)
dat$current_age<-2020-dat[[birthYear_col]]
pdf_fn<-paste(sep=".",out,"currentAge.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat,aes(x=current_age,fill=strat)) + geom_density(alpha=0.7,position="identity")  + theme_bw()  +
    labs(x="Estimated Current Age", caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)~","~N['NA']~"="~.(stratNA)),title=main) + facet_wrap(~facet) +
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10),strip.text=element_blank()) + scale_fill_manual(na.value="grey",values=colors,name=legend,labels=labels) 
dev.off()

pdf_fn<-paste(sep=".",out,"currentAge_hist.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat,aes(x=current_age,fill=strat)) + geom_histogram(alpha=0.7,binwidth=1,position="identity")  + theme_bw() +
    labs(x="Estimated Current Age", caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)~","~N['NA']~"="~.(stratNA)),title=main) + facet_wrap(~facet) +
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10),strip.text=element_blank()) + scale_fill_manual(na.value="grey",values=colors,name=legend,labels=labels) 
dev.off()

pdf_fn<-paste(sep=".",out,"currentAge_area.pdf")
pdf(file=pdf_fn,height=6,width=8,useDingbats=FALSE)
ggplot(dat,aes(x=current_age,fill=strat)) + theme_bw() +
    geom_area(aes(y = ..count.., fill = strat, group = strat), position="identity",stat = "bin",alpha=0.7) +
    labs(x="Estimated Current Age", caption=bquote(N[negative]~"="~.(strat0)~","~N[positive]~"="~.(strat1)~","~N['NA']~"="~.(stratNA)),title=main) + facet_wrap(~facet) +
    theme(plot.caption=element_text(hjust=0.5),legend.text=element_text(size=15),title=element_text(size=15),axis.title=element_text(size=15),axis.text.x=element_text(size=10),strip.text =element_blank()) + scale_fill_manual(na.value="grey",values=colors,name=legend,labels=labels)
dev.off()


    
### glm to predict age using strata
glm.obj<-glm(get(names(dat)[age_col])~get(names(dat)[strat_col]),data=dat)
print(summary(glm.obj))

