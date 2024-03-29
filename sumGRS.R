#!/usr/bin/Rscript

#### Revised from Dr.Sarah Graham

### first argument is config file of file names
### second argument is output file prefix
### writes file with individual, GRS, and invNorm_GRS 

library(data.table)
library(dplyr)
library(RNOmni)

args=commandArgs(trailingOnly=TRUE)
config<-args[1]
out<-args[2]

config_df<-fread(config,header=FALSE)
file_list<-config_df$V1

data_list <- lapply(file_list, function(x)fread(x)) #read in all files

data_list <- lapply(data_list, function(x) unique(x, by="individual")) #unique samples

data_all <- Reduce(function(...) left_join(...,by=c("individual"="individual")), data_list) #join into 1

data_all$GRS <- rowSums(data_all[,2:ncol(data_all)],na.rm=TRUE)  #row sums

data_all <- data_all[,c("individual", "GRS")] #subset to id and sum
data_all$invNorm_GRS <-rankNorm(data_all$GRS) #inverse normalize 

write.table(data_all, paste(sep=".",out,"total.txt"), col.names=TRUE, row.names=FALSE,sep="\t",quote=FALSE)




