#!/usr/bin/Rscript

###########################################################
###################### Helper Functions  ##################
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


##function to estimate FHiGRS          
#takes data frame created from prev_per_quantile_stratum object
estimate_FHiGRS<-function(prev_df,main_df,strat_col,grs_col){
    ##put q-quantiles in order of prevalence
    ##to do: use prevalence from reference population
    prev_df_order<-prev_df[order(prev_df$prev),]
    prev_df_order$rank<-seq(100,nrow(prev_df)*100,100) #make rankings in scale of 100

    ##assign samples to groups
    main_df$rank<-as.numeric(1) #initialize column for FHiGRS
    for (s in c(1,2)){
        prev_df_order_strat<-prev_df_order[prev_df_order$stratum==(s-1),]
        prev_df_order_strat[1,'lower_tile']<- -100 #condition to put samples equiv to minimum in bottom bin
        for (j in 1:nrow(prev_df_order_strat)){
            u<-prev_df_order_strat[j,'upper_tile']
            l<-prev_df_order_strat[j,'lower_tile']
            rank<-prev_df_order_strat[j,'rank']
            main_df[main_df[[strat_col]]==(s-1) & main_df[[grs_col]]>l & main_df[[grs_col]]<=u,'rank']<-rank
        }
    }
    main_df$grs_rank<-main_df[[grs_col]] + main_df$rank #add rank to GRS
    main_df$FHIGRS<-qnorm((rank(main_df$grs_rank,na.last="keep")-0.5)/sum(!is.na(main_df$grs_rank))) #new FHIGRS
    return(main_df)
}


