#!/usr/bin/Rscript

######### Libraries #########

library(optparse)
library(data.table)

################### Read Command Line Parameters #############################
optionList <- list(
    make_option(c("-f", "--file"), type="character", help="File with file paths per line to check"),
    make_option(c("-r","--row"),type="integer",help="Expected number of rows"),
    make_option(c("-c","--col"),type="integer",help="Expected numbere of columns"),
    make_option(c("-o","--out"),type="character",help="Output file for bad files")
    )

parser <- OptionParser(usage="%prog -f file -c column -r row -o output", option_list=optionList)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
    optionStrings <- character()
    for(item in parserObj@options) {
        optionStrings <- append(optionStrings, c(item@short_flag, item@long_flag))
    }
    optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- parse_args(parser, positional_arguments=TRUE)
print(arguments$options)
############################# MAIN ###########################

#warn if existing file
if (file.exists(arguments$options$out)){
    stop(sprintf("%s already exists. Please delete because this script appends to it\n",arguments$options$out))
}

#read in the new line separatetd list of expected files 
files<-fread(arguments$options$file,header=F)
for (i in 1:nrow(files)){
    fn<-files$V1[i]
    if (file.exists(fn)){
        df<-fread(fn,header=F)
        dimensions<-dim(df)
        ##check dimensions
        if ((arguments$options$row != dimensions[1]) | (arguments$options$col != dimensions[2])){
            write.table(x=fn,file=arguments$options$out,append=TRUE,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
        }
    } else { #file doesn't exist at all
        write.table(x=fn,file=arguments$options$out,append=TRUE,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
    }
}


