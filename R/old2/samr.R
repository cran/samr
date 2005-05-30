

samr.const.red.color <- 3
samr.const.green.color <- 10
samr.const.black.color <- 1

samr.const.quantitative.response <- 0
samr.const.twoclass.unpaired.response <- 1
samr.const.survival.response <- 2
samr.const.multiclass.response <- 3
samr.const.oneclass.response <- 4
samr.const.twoclass.paired.response <- 5
samr.const.twoclass.unpaired.timecourse.response <- 6
samr.const.oneclass.timecourse.response <- 7
samr.const.twoclass.paired.timecourse.response <- 8



samr <- function(data, params, resp.type=1, s0=NULL, s0.perc=NULL, nperms=100, center.arrays=FALSE, testStatistic=c("standard","wilcoxon"), time.summary.type=c("slope","signed.area"),  xl.mode=c("regular","firsttime","next20","lasttime"), xl.time=NULL,  xl.prevfit=NULL){


##SAM method. copyright june 2000: Goss, Tibshirani and Chu.
## coded by r tibshirani; 
## y is response measure: 1,2 for two twoclass groups,
## y=1,1,1 for onesample problem,
## -1, 1,2-,2 for paired groups, with -1 paired with 1 etc
## or survival time, or categorical for unordered groups 1,2,3,..
## quantitative for continuous ordered y#
## new resp.type- timearea
##
## params is a list containing user-defined parameters; at this point just has
##   min.foldchange (NULL means no fold change criterion, o/w should be> 1)
##
## s0 is the exchangeability factor; you can specify
## s0 as an actual value
## s0.perc, the percentile of sd values to use for s0
##  or if both s0 and s0.perc are null (the default), then s0 is automatically estimated





## returns
##  evo= expected order statistics (length p=# of genes)
##  tt=numer/(sd+s0)  test statistics on original data (and the ingredients)
##  ttstar0= p by nperms matrix of test statistics on permted data
##  ttstar= ttstar0 with columns sorted (largest values in row 1)
##  also returns permuted values: foldchange.star,  ystar, sdstar, statusstar (for survival data)
##  in xl.mode firsttime or next 20, the function  also  returns  x , which   data$x, except for time-course data,
##   where it is computed from in this function
#  from time summaries of data$x. However, this quantity is deleted the last time the function is called,
#   as it is very large and not needed further
        



  this.call=match.call()  
xl.mode=match.arg(xl.mode)


  
  
if(xl.mode=="regular" | xl.mode=="firsttime"){
  
# initialize some things (harmlessly), just so that xl.mode will work correctly

  
x=NULL
ttstar0=NULL
evo=NULL
ystar=NULL
ystar.keep=NULL
sdstar.keep=NULL
statusstar=NULL
status.keep=NULL
sdstar=NULL
pi0=NULL
stand.contrasts=NULL
stand.contrasts.95=NULL
foldchange=NULL
foldchange.star=NULL
perms=NULL
permsy=NULL
##


testStatistic <- match.arg(testStatistic)
time.summary.type <- match.arg(time.summary.type)
  
  x=data$x

 y=data$y
 argy=y
 resp.type.arg=resp.type
  
  if(resp.type<0 | resp.type>8){
     stop("Error in response type specification")
   }
  
are.blocks.specified=FALSE


# center columns of  array data if requested

if(center.arrays){
 x<-scale(x,center=apply(x,2,median),scale=FALSE)
}

# make sure -1,1,, etc are non-character values  coming from Excel
 if(resp.type==samr.const.twoclass.paired.response){y=as.numeric(y)}

# check if there are blocks for 2 class unpaired case
if(resp.type==samr.const.twoclass.unpaired.response){
    if(substring(y[1],2,6)=="Block"){
       junk=parse.block.labels.for.2classes(y)
       y=junk$y; blocky=junk$blocky
       are.blocks.specified=TRUE
    }}

# parse and summarize, if time course data

if(resp.type==samr.const.twoclass.unpaired.timecourse.response | 
 resp.type==samr.const.twoclass.paired.timecourse.response |
resp.type==samr.const.oneclass.timecourse.response){
    junk=parse.time.labels.and.summarize.data(x,y, resp.type, time.summary.type)
       y=junk$y;
       x=junk$x;
  }

# if the data is timecourse, we have already summarized the time aspect.
# Thus we change the resp.type to the appropriate non-time-course type. Note that the original value
#  of resp.type was saved above in resp.type.arg
  
if(resp.type==samr.const.twoclass.unpaired.timecourse.response){ resp.type=samr.const.twoclass.unpaired.response}
if(resp.type==samr.const.twoclass.paired.timecourse.response){ resp.type=samr.const.twoclass.paired.response}
if(resp.type==samr.const.oneclass.timecourse.response){ resp.type=samr.const.oneclass.response}


  
stand.contrasts=NULL
stand.contrasts.95=NULL

  if(resp.type==samr.const.survival.response){status=data$status}


# do a thorough error  checking of the response data
  
 check.format(y,resp.type=resp.type,status=status)
  
  

  n <- nrow(x)
  ny <- length(y)
  sd <- NULL;numer <- NULL



  ystar.keep <- matrix(0,ncol=nperms,nrow=length(y))
  sdstar.keep <- matrix(0,ncol=nperms,nrow=nrow(x))

  status.keep <- NULL
  if(resp.type==samr.const.survival.response){status.keep <- ystar.keep}



# for wilcoxon, s0 is not needed
if(testStatistic=="wilcoxon"){s0=0;s0.perc=-1}

# estimate s0 if necessary
  if(is.null(s0)){

     if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="standard"){junk <- ttest.func(x,y);numer <- junk$numer;sd <- junk$sd}

 if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="wilcoxon"){junk <- wilcoxon.func(x,y);numer <- junk$numer;sd <- junk$sd}

     if(resp.type==samr.const.oneclass.response){junk <- onesample.ttest.func(x,y);numer <- junk$numer;sd <- junk$sd}
     if(resp.type==samr.const.twoclass.paired.response){junk <- paired.ttest.func(x,y); numer <- junk$numer ; sd <- junk$sd}

     if(resp.type==samr.const.survival.response){junk <- cox.func(x,y,status);numer <- junk$numer;sd <- junk$sd}
     if(resp.type==samr.const.multiclass.response){junk <- multiclass.func(x,y); numer <- junk$numer;sd <- junk$sd}
     if(resp.type==samr.const.quantitative.response){junk <- quantitative.func(x,y);numer <- junk$numer;sd <- junk$sd}
     

     if(!is.null(s0.perc)){
if((s0.perc != -1 & s0.perc < 0) | s0.perc > 100){
       stop("Illegal value for s0.perc: must be between 0 and 100, or equal
to (-1) (meaning that s0 should be set to zero)")
       }
           if(s0.perc== -1){s0=0}
            if(s0.perc>=0) {s0 <- quantile(junk$sd,s0.perc/100)}
          }
     if(is.null(s0.perc)){   
      s0=est.s0(junk$tt,junk$sd)$s0.hat
      s0.perc=100*sum(junk$sd<s0)/length(junk$sd)
  }
}
 

# compute test statistics on original data


  if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="standard"){tt <- ttest.func(x,y,s0=s0)$tt}
  if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="wilcoxon"){tt <- wilcoxon.func(x,y,s0=s0)$tt}
  if(resp.type==samr.const.oneclass.response){tt <- onesample.ttest.func(x,y,s0=s0)$tt}
  if(resp.type==samr.const.twoclass.paired.response){tt <- paired.ttest.func(x,y,s0=s0)$tt}
  if(resp.type==samr.const.survival.response){tt <- cox.func(x,y,status,s0=s0)$tt}
  if(resp.type==samr.const.multiclass.response){
   junk2 <- multiclass.func(x,y,s0=s0)
   tt=junk2$tt
   stand.contrasts=junk2$stand.contrasts
}

  if(resp.type==samr.const.quantitative.response){tt <- quantitative.func(x,y,s0=s0)$tt}

  

# construct matrix of permutations 


if(resp.type==samr.const.quantitative.response |  resp.type==samr.const.multiclass.response |  resp.type==samr.const.survival.response){

   junk<- getperms(y, nperms)
   perms=junk$perms;all.perms.flag=junk$all.perms.flag; nperms.act=junk$ nperms.act
 }

if(resp.type==samr.const.twoclass.unpaired.response){
       if(are.blocks.specified){
       junk=compute.block.perms(y,blocky, nperms)
       permsy=matrix(junk$permsy,ncol=length(y))
   all.perms.flag=junk$all.perms.flag; nperms.act=junk$nperms.act
    }
    else{junk<- getperms(y, nperms)
   permsy= matrix(y[junk$perms],ncol=length(y)) ;all.perms.flag=junk$all.perms.flag; nperms.act=junk$nperms.act
 } }

  
if(resp.type==samr.const.oneclass.response){

    allii= 0:((2^length(y))-1)
    nperms.act=2^length(y)
    all.perms.flag=1
    if((2^length(y))>nperms){
       allii=sample(allii,size=nperms)
       nperms.act=nperms
       all.perms.flag=0
    }
       
     permsy=matrix(NA,nrow=nperms.act,ncol=length(y))
     k=0 
     for(i in  allii ){
       junk=integer.base.b(i,b=2)
       if(length(junk)<length(y)){
          junk=c(rep(0,length(y)-length(junk)), junk)
      }
       k=k+1
       permsy[k,]=y*(2*junk-1)
     }
  }



 if(resp.type==samr.const.twoclass.paired.response){

 junk=compute.block.perms(y,abs(y), nperms)
 permsy=junk$permsy;all.perms.flag=junk$all.perms.flag; nperms.act=junk$nperms.act
}




# compute test statistics on permuted  data


  

  ttstar <- matrix(0,nrow=nrow(x),ncol=nperms)

  foldchange.star=NULL

  if(params$min.foldchange>0){foldchange.star <- matrix(0,nrow=nrow(x),ncol=nperms)}

if(resp.type==samr.const.multiclass.response){
stand.contrasts.star=array(NA,c(nrow(x),length(table(y)),nperms))
}

  # end of if(xltime=="regular" etc
}

  
if(xl.mode=="next20" |  xl.mode=="lasttime"){

  # get stuff from prevfit

 evo= xl.prevfit$evo
 tt= xl.prevfit$tt
 numer=xl.prevfit$numer
 sd=xl.prevfit$sd
 ttstar= xl.prevfit$ttstar
 ttstar0= xl.prevfit$ttstar0
 n= xl.prevfit$n 
pi0= xl.prevfit$pi0 
foldchange= xl.prevfit$foldchange 
 y= xl.prevfit$y
 x=xl.prevfit$x
 argy= xl.prevfit$argy 
 testStatistic= xl.prevfit$testStatistic 
 foldchange.star= xl.prevfit$foldchange.star 
 s0= xl.prevfit$s0 
 s0.perc= xl.prevfit$s0.perc 
 ystar= xl.prevfit$ystar
 ystar.keep=xl.prevfit$ystar.keep
 resp.type= xl.prevfit$resp.type 
 resp.type.arg= xl.prevfit$resp.type.arg 
statusstar= xl.prevfit$statusstar
 status.keep=xl.prevfit$status.keep
 sdstar= xl.prevfit$sdstar
  sdstar.keep= xl.prevfit$sdstar.keep
 resp.type= xl.prevfit$resp.type 
stand.contrasts= xl.prevfit$stand.contrasts 
stand.contrasts.95= xl.prevfit$stand.contrasts.95
 perms=xl.prevfit$perms
  permsy=xl.prevfit$permsy
 nperms= xl.prevfit$nperms 
 nperms.act= xl.prevfit$nperms.act 
 all.perms.flag= xl.prevfit$all.perms.flag 

}



  
if(xl.mode=="regular"){
    first=1;last=nperms.act
  }
  if(xl.mode=="firsttime"){
    first=1;last=1
  }
  if(xl.mode=="next20"){
    first=xl.time; last= min(xl.time+19, nperms.act-1)
  }
  if(xl.mode=="lasttime"){
    first=nperms.act;last=nperms.act
  }



  for(b in first:last){

    cat(c("perm=",b),fill=T)
    xstar <- x


    if(resp.type==samr.const.oneclass.response){

         ystar=permsy[b,]

         ttstar[,b] <- onesample.ttest.func(xstar,ystar,s0=s0)$tt

    }




    if(resp.type==samr.const.twoclass.paired.response){

      ystar=permsy[b,]
      
      ttstar[,b] <- paired.ttest.func(xstar,ystar,s0=s0)$tt

      ystar.keep[,b] <- ystar

      if(params$min.foldchange>0){foldchange.star[,b]=foldchange.paired(xstar,ystar,data$logged2)}



    }



    if(resp.type==samr.const.twoclass.unpaired.response){

       ystar=permsy[b,]
       
if(testStatistic=="standard"){
      junk <- ttest.func(xstar,ystar,s0=s0)
}
if(testStatistic=="wilcoxon"){
      junk <- wilcoxon.func(xstar,ystar,s0=s0)
}



      ttstar[,b] <- junk$tt

      ystar.keep[,b] <- ystar

      sdstar.keep[,b] <- junk$sd

      if(params$min.foldchange>0){foldchange.star[,b]=foldchange.twoclass(xstar,ystar,data$logged2)}



    }



    if(resp.type==samr.const.survival.response){

      o <- perms[b,]

      ttstar[,b] <- cox.func(xstar,y[o],status=status[o],s0=s0)$tt

      ystar.keep[,b] <- y[o]

      status.keep[,b] <- status[o]

    }

    if(resp.type==samr.const.multiclass.response){

      ystar= y[perms[b,]]

      junk <- multiclass.func(xstar,ystar,s0=s0)

      ttstar[,b] <- junk$tt

      ystar.keep[,b] <- ystar

      sdstar.keep[,b] <- junk$sd
     stand.contrasts.star[,,b]=junk$stand.contrasts
    }

    if(resp.type==samr.const.quantitative.response){

        ystar= y[perms[b,]]

      junk <- quantitative.func(xstar,ystar,s0=s0)

      ttstar[,b] <- junk$tt

      ystar.keep[,b] <- ystar

      sdstar.keep[,b] <- junk$sd

    }
# end of xl.mode=="regular" | xl.mode=="one.time"
  }

  

# sort columns of statistics from permuted samples, and compute expected order statistics

  
if(xl.mode=="regular" | xl.mode=="lasttime"){
  ttstar0 <- ttstar

  for(j in 1:ncol(ttstar)) {

    ttstar[, j] <- -1 * sort(-1 * ttstar[, j])

  }

  for(i in 1:nrow(ttstar)) {

    ttstar[i,  ] <- sort(ttstar[i,  ])

  }

  evo <- apply(ttstar, 1, mean)

  evo <- evo[length(evo):1]

  ystar <- ystar.keep;

  statusstar <- status.keep

  sdstar <- sdstar.keep



# estimation of pi0= prop of null genes

  
if(resp.type!=samr.const.multiclass.response){
  qq<-quantile(ttstar,c(.25,.75))
}

if(resp.type==samr.const.multiclass.response){
  qq<-quantile(ttstar,c(0,.50))
}


  pi0 <- sum(tt>qq[1] & tt< qq[2])/(.5*length(tt))



 # compute fold changes, when applicable

  

  foldchange=NULL

  if(params$min.foldchange>0){

    if(resp.type==samr.const.twoclass.unpaired.response){ foldchange=foldchange.twoclass(data$x,data$y,data$logged2)}

    if(resp.type==samr.const.twoclass.paired.response){ foldchange=foldchange.paired(data$x,data$y,data$logged2)}

  }

  
  
stand.contrasts.95=NULL

if(resp.type==samr.const.multiclass.response){
  stand.contrasts.95=quantile(stand.contrasts.star,c(.025,.975))
 } 

  # th last time through, we delete  x, since it is very big and is not needed further
x=NULL
}

  return(list(evo=evo,tt=tt,ttstar=ttstar,ttstar0=ttstar0,n=n,pi0=pi0,foldchange=foldchange, y=y,argy=argy, testStatistic=testStatistic,
              foldchange.star=foldchange.star,  s0=s0, s0.perc=s0.perc,  numer=numer,sd=sd, ystar=ystar, ystar.keep=ystar.keep, sdstar.keep=sdstar.keep, status.keep=status.keep, resp.type=resp.type, resp.type.arg=resp.type.arg, x=x,

              statusstar=statusstar, sdstar=sdstar, resp.type=resp.type,stand.contrasts=stand.contrasts,
stand.contrasts.95=stand.contrasts.95,nperms=nperms, nperms.act=nperms.act, perms=perms, permsy=permsy, all.perms.flag=all.perms.flag, call=this.call))

}
