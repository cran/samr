

## Just for memory sake...

##samr.const.quantitative.response <- 0

##samr.const.twoclass.unpaired.response <- 1

##samr.const.survival.response <- 2

##samr.const.multiclass.response <- 3

##samr.const.oneclass.response <- 4

##samr.const.twoclass.paired.response <- 5


##samr.const.twoclass.unpaired.timecourse.response <- 6
##samr.const.oneclass.timecourse.response <- 7
##samr.const.twoclass.paired.timecourse.response <- 8




##SAM R associated functions



## individual functions for each response type



ttest.func <- function(x,y,s0=0){

  n1 <- sum(y==1)
  n2 <- sum(y==2)
  
  p <- nrow(x)
  m1 <- rowMeans(x[,y==1])
  m2 <- rowMeans(x[,y==2])


  sd <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m1) + (n1-1) * varr(x[, y==1], meanx=m2) )*

             (1/n1+1/n2)/(n1+n2-2) )

  numer <-  m2 - m1

  dif.obs <- (numer)/(sd + s0)
  return(list(tt=dif.obs,numer=numer, sd=sd))
}

  
wilcoxon.func <- function(x,y,s0=0){
  
  n1 <- sum(y==1)
  n2 <- sum(y==2)
p=nrow(x)
  
r2= rowSums(t(apply(x,1,rank))[,y==2])

numer=r2- (n2/2)*(n2+1) -(n1*n2)/2

sd=sqrt(n1*n2*(n1+n2+1)/12)

tt= (numer)/(sd + s0)

  return(list(tt=tt,numer=numer, sd=rep(sd,p)))
}





onesample.ttest.func <- function(x,y,s0=0){

  n <- length(y)

  

  x <- x*matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=T)

  m <- rowMeans(x)

  

  

  

  sd <- sqrt( varr(x, meanx=m)/n )

  dif.obs <- m/(sd + s0)

  

  

  return(list(tt=dif.obs, numer=m,sd=sd))

}



paired.ttest.func <- function(x,y,s0=0,useden=T){

  nc <- ncol(x)/2

  o <- 1:nc

  o1 <- rep(0,ncol(x)/2);o2 <- o1

  for(j in 1:nc){o1[j] <- (1:ncol(x))[y==-o[j]]}

  for(j in 1:nc){o2[j] <- (1:ncol(x))[y==o[j]]}

  

  d <- x[,o2]-x[,o1]

  su <- x[,o2]+x[,o1]

  

  den <- 1

  sd <- NULL

  

  if(is.matrix(d)){ 

    m <-  rowMeans(d)

  }

  

  if(!is.matrix(d)) {m <- mean(d)}

  

  if(useden){

    if(is.matrix(d)){ sd <- sqrt(varr(d, meanx=m)/nc)}

    if(!is.matrix(d)){sd <- sqrt(var(d)/nc)}

    

    

    den <- sd+s0

  }



  

  dif.obs <- m/den

  

  return(list(tt=dif.obs, numer=m, sd=sd))

}



cox.func <- function(x,y,status,s0=0){

  scor <- coxscor(x,y, status)$scor

  sd <- sqrt(coxvar(x,y, status))

  tt <- scor/(sd+s0)

  

  return(list(tt=tt, numer=scor, sd=sd))

}



multiclass.func <- function(x,y,s0=0){

  ##assumes y is coded 1,2...

  nn <- table(y)

  m <- matrix(0,nrow=nrow(x),ncol=length(nn))

  v <- m

  for(j in 1:length(nn)){

    m[,j] <- rowMeans(x[,y==j])

    v[,j] <- (nn[j]-1)*varr(x[,y==j], meanx=m[,j])

  }

  mbar <- rowMeans(x)

  mm <- m-matrix(mbar,nrow=length(mbar),ncol=length(nn))

  fac <- (sum(nn)/prod(nn))

  scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=T)*mm*mm,1,sum)))

  sd <- sqrt(apply(v,1,sum)*(1/sum(nn-1))*sum(1/nn))

  tt <- scor/(sd+s0)

  mm.stand=t(scale(t(mm),center=FALSE,scale=sd))

  return(list(tt=tt, numer=scor, sd=sd,stand.contrasts=mm.stand))

}





quantitative.func <- function(x,y,s0=0){

  yy <- y-mean(y)

  temp <- x%*%yy

  xx <- t(scale(t(x),scale=F))

  sxx <- apply(xx^2,1,sum)

  scor <- temp/sxx

  b0hat <- mean(y)-scor*apply(x,1,mean)

  yhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+x*matrix(scor,nrow=nrow(x),ncol=ncol(x))

  ty <- matrix(y,nrow=nrow(yhat),ncol=ncol(yhat),byrow=T)

  sigma <- sqrt(apply( (ty-yhat)^2,1,sum)/(ncol(yhat)-2))

  sd <- sigma/sqrt(sxx)

  tt <- scor/(sd+s0)

  

  return(list(tt=tt, numer=scor, sd=sd))

}



timearea.func <- function(x, y, s0 = 0) {

  n <- ncol(x)

  xx <- 0.5 * (x[, 2:n] + x[, 1:(n - 1)]) * matrix(diff(y), nrow = nrow(x), ncol = n - 1, byrow = T)

  numer <- rowMeans(xx)

  sd <- sqrt(varr(xx, meanx=numer)/n)

  tt <- numer/sqrt(sd + s0)

  return(list(tt = tt, numer = numer, sd = sd))

}



#########################################



samr.detec.slab <- function(a, del, params) {

  ## find genes above and below the slab of half-width del

  

  n <- length(a$tt)

  tt <- a$tt

  evo <- a$evo

  numer <- a$tt*(a$sd+a$s0)

  tag <- order(tt)

  pup <- NULL

  foldchange.cond.up=rep(T,length(evo))

  foldchange.cond.lo=rep(T,length(evo))

  

  if(!is.null(a$foldchange[1]) & (params$min.foldchange>0)){

    foldchange.cond.up= a$foldchange >= params$min.foldchange

    foldchange.cond.lo= a$foldchange <= 1/params$min.foldchange

  }

  o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0 & foldchange.cond.up  ]

  if(length(o1) > 0) {

    o1 <- o1[1]

    o11 <- o1:n

    o111 <- rep(F, n)

    o111[tag][o11] <- T

    pup <- (1:n)[o111]

  }

  plow <- NULL

  o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0 & foldchange.cond.lo ]

  if(length(o2) > 0) {

    o2 <- o2[length(o2)]

    o22 <- 1:o2

    o222 <- rep(F, n)

    o222[tag][o22] <- T

    plow <- (1:n)[o222]

  }

  

  return(list(plow=plow, pup=pup))

}



f55 <- function(aa) {

  length(aa$pl) + length(aa$pu)

}



samr.compute.delta.table <- function(a, params, dels=NULL, nvals=50) {
                                        # computes delta table, starting with samr object "a", for nvals values of delta

  dels=seq(0,max(abs(sort(a$tt)-a$evo)), length=nvals)

  col=matrix(1,nrow=length(a$evo),ncol=nvals)

  ttstar0 <- a$ttstar0
  tt <- a$tt
  n <- a$n
  evo <- a$evo
  nsim <- ncol(ttstar0)
  res1 <- NULL



    foldchange.cond.up=matrix(T,nrow=nrow(a$ttstar),ncol=ncol(a$ttstar))
    foldchange.cond.lo=matrix(T,nrow=nrow(a$ttstar),ncol=ncol(a$ttstar))

    if(!is.null(a$foldchange[1]) & (params$min.foldchange>0)){

      foldchange.cond.up= a$foldchange.star >= params$min.foldchange

      foldchange.cond.lo= a$foldchange.star <= 1/params$min.foldchange
    }

  
  fr(ii in 1:length(dels)) {
    
     
     ttt <- samr.detec.slab(a,dels[ii], params)
    cutup <- 10e9
    if(length(ttt$pup>0)){ cutup <- min(a$tt[ttt$pup])}
    cutlow <- -10e9
    if(length(ttt$plow)>0){cutlow <- max(a$tt[ttt$plow])}


    resup=a$ttstar0 > cutup & foldchange.cond.up
    reslow=a$ttstar0 < cutlow & foldchange.cond.lo



  
     errup<-colSums(resup) 
    errlow <- colSums(reslow)
    s <- sqrt(var(errup)/nsim + var(errlow)/nsim)
    gmed <- median(errup + errlow)
    g90=quantile(errup + errlow, .90)

    cat(c(dels[ii], errlow + errup), fill = T)

    cat("", fill = T)
    g2 <- f55(ttt)
    res1 <- rbind(res1, c(a$pi0*gmed, a$pi0*g90, g2, a$pi0*gmed/g2, a$pi0*g90/g2, cutlow,cutup))

  }
       
  res1 <- cbind(dels, res1)
  dimnames(res1) <- list(NULL, c("delta", "# med false pos",
                                 "90th perc false pos", "# called", "median FDR", "90th perc FDR","cutlo","cuthi"))

  return(res1)
}



detec.horiz <- function(a,cutlow,cutup, params){

  ## find genes above or below horizontal cutpoints

  dobs <- a$tt

  n <- length(dobs)

  foldchange.cond.up=rep(T,n)

  foldchange.cond.lo=rep(T,n)

  

  if(!is.null(a$foldchange[1]) & (params$min.foldchange>0)){

    foldchange.cond.up= a$foldchange >= params$min.foldchange

    foldchange.cond.lo= a$foldchange <= 1/params$min.foldchange

  }

  

  pup <- (1:n)[dobs> cutup & foldchange.cond.up]

  plow <- (1:n)[dobs< cutlow & foldchange.cond.lo]

  

  return(list(plow=plow,pup=pup))

  

}





samr.plot <- function(a, del, params) {

                                        ## make observed-expected plot

                                        ##

                                        ## takes foldchange into account too



  

  LARGE=10e9

  b <-  samr.detec.slab(a, del, params)

  

  bb <- c(b$pup,b$plow)

  b1= LARGE

  b0=-LARGE

  

  if(!is.null(b$pup)){b1 <- min(a$tt[b$pup])}

  if(!is.null(b$plow)){b0 <- max(a$tt[b$plow])}

  

  c1 <- (1:a$n)[sort(a$tt)>=b1]

  c0 <- (1:a$n)[sort(a$tt)<=b0]

  c2 <- c(c0,c1)

  

  foldchange.cond.up=rep(T,length(a$evo))

  foldchange.cond.lo=rep(T,length(a$evo))

  

  if(!is.null(a$foldchange[1]) & (params$min.foldchange>0)){

    foldchange.cond.up= a$foldchange >= params$min.foldchange

    foldchange.cond.lo= a$foldchange <= 1/params$min.foldchange

  }

  

  col=rep(1,length(a$evo))

  col[b$plow]=3

  col[b$pup]=2

  if(!is.null(a$foldchange[1]) & (params$min.foldchange>0)){

    col[!foldchange.cond.lo & !foldchange.cond.up]=1

  }

  col.ordered=col[order(a$tt)]

  


    ylims <- range(a$tt)

    xlims <- range(a$evo)

    plot(a$evo,sort(a$tt),xlab="expected score", ylab="observed score",ylim=ylims, 

         xlim=xlims, type="n")

    points(a$evo,sort(a$tt),col=col.ordered)

    

    

    abline(0,1)

    abline(del,1,lty=2)

    abline(-del,1,lty=2)


}





stand <- function(x, use.for.stand = 1:ncol(x)) {

  n <- ncol(x)

  xx <- sign(x) * (abs(x))^{

    1/3

  }

  stn <- apply(xx[, use.for.stand], 1, mean)

  xxx <- xx - xx

  for(j in 1:n) {

    a <- lsfit(stn, xx[, j])

    xxx[, j] <- (xx[, j] - a$coef[1])/a$coef[2]

  }

  x3 <- xxx^3

  return(x3)

}





localfdr <- function(a,  params, perc=.01, df=10) {

  

  ## estimates local fdr at score "d", using SAM object "a"

  ## "d" can be a vector of d scores

  ## returns estimate of symmetric fdr

  ## to use: first run SAM in Splus, and then pass the resulting fit object to

  ## localfdr
  ## NOTE: at most 20 of the perms are used to estimate the fdr (for speed sake)

 nperms.to.use=min(20,ncol(a$ttstar))

  d=seq(min(a$tt),max(a$tt),length=100)

  ndscore <- length(d)

  dvector <- rep(NA,ndscore)

  

  ind.foldchange=rep(T,length(a$tt))

  if(!is.null(a$foldchange[1]) & params$min.foldchange>0){

    ind.foldchange= (a$foldchange>=  params$min.foldchange) | (a$foldchange<=  params$min.foldchange)

  }

  

  for(i in 1:ndscore)

    {
      pi0<-a$pi0
      r <- sum(a$tt<d[i])
      r22 <- max(r-length(a$tt)*perc/2, 1)
      dlow.sym <- sort(a$tt)[r22]

      if(d[i]<0)
        {
          r2 <- max(r-length(a$tt)*perc, 1)
          dlow <- sort(a$tt)[r2]
          dup <- d[i]
        }

      r22 <- min(r+length(a$tt)*perc/2, length(a$tt))
      dup.sym <- sort(a$tt)[r22]

      if(d[i]>0)
        {
          r2 <- min(r+length(a$tt)*perc, length(a$tt))
          dup <- sort(a$tt)[r2]
          dlow <- d[i]

        }
      o <- a$tt>dlow & a$tt< dup & ind.foldchange
      oo <- a$tt>dlow.sym & a$tt< dup.sym & ind.foldchange

      nsim <- ncol(a$ttstar)
      fdr <- rep(NA,nsim)
      fdr2 <- fdr


          if(!is.null(a$foldchange[1]) &  params$min.foldchange>0){
             temp=as.vector(a$foldchange.star[,1:nperms.to.use])
            ind.foldchange=(temp >=  params$min.foldchange) | (temp <=  params$min.foldchange)
          }

          

         temp=as.vector(a$ttstar0[,1:nperms.to.use])
          fdr <- pi0*(temp>dlow & temp<dup &ind.foldchange)

          fdr2 <- pi0*(temp > dlow.sym & temp <dup.sym & ind.foldchange)


      fdr <- (sum(fdr)/nperms.to.use)/sum(o)
      fdr.sym <- (sum(fdr2)/nperms.to.use)/sum(oo)
      dlow.sym <- dlow.sym
      dup.sym <-dup.sym
      dlow <- dlow
      dup <- dup

      dvector[i] <- fdr.sym

    }

  

om=!is.na(dvector) & (dvector!=Inf)
aa=smooth.spline(d[om], dvector[om], df=df)

  

  return(list(smooth.object=aa,perc=perc,df=df))

}



predict.localfdr= function(smooth.object,d){
 yhat=predict(smooth.object,d)$y
yhat=pmin(yhat,1)
yhat=pmax(yhat,0)

  return(yhat)

}



samr.compute.siggenes.table=function(a,del, data, delta.table,  params, all.genes=FALSE){



  ## computes significant genes table, starting with samr object "a" and "delta.table"

  ##  for a  **single** value del
# if all.genes is true, all genes are printed (and value of del is ignored)

if(!all.genes){
  sig=samr.detec.slab(a, del, params)
}

if(all.genes){
  p=length(a$tt)
  pup=(1:p)[a$tt>=0]
  plo=(1:p)[a$tt<0]
  sig=list(pup=pup, plo=plo)
}

  aa=localfdr(a, params)

  if(length(sig$pup)>0){fdr.up=predict.localfdr(aa$smooth.object,a$tt[sig$pup])}

 if(length(sig$plo)>0){ fdr.lo=predict.localfdr(aa$smooth.object,a$tt[sig$plo])}

  qvalues=qvalue.func(a,sig, delta.table)

  
done=FALSE
# two class unpaired or paired with fold change specified
  if((a$resp.type==samr.const.twoclass.unpaired.response | a$resp.type==samr.const.twoclass.paired.response) & (params$min.foldchange>0)){


    res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],a$tt[sig$pup],a$numer[sig$pup],a$sd[sig$pup],a$ffoldchange[sig$pup],qvalues$qvalue.up, fdr.up)

    dimnames(res.up)=list(NULL,c("Row","Gene Name","Gene ID", "Score(d)", "Numerator(r)","Denominator(s+s0)", "Fold Change", "q-value","localfdr"))

    res.lo=cbind(sig$plo+1,data$genenames[sig$plo],data$geneid[sig$plo],a$tt[sig$plo],a$numer[sig$plo],a$sd[sig$plo],   a$foldchange[sig$plo],qvalues$qvalue.lo, fdr.lo)

done=TRUE
  }


# multiclass
 if(a$resp.type==samr.const.multiclass.response){


    res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],a$tt[sig$pup],a$numer[sig$pup],a$sd[sig$pup],a$stand.contrasts[sig$pup,],qvalues$qvalue.up, fdr.up)

collabs.contrast=paste("contrast-",as.character(1:ncol(a$stand.contrasts)),sep="")
    dimnames(res.up)=list(NULL,c("Row","Gene Name","Gene ID", "Score(d)", "Numerator(r)","Denominator(s+s0)",collabs.contrast,  "q-value","localfdr"))


res.lo=NULL
done=TRUE
  }





#all other cases


if(!done){
    res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],a$tt[sig$pup],a$numer[sig$pup],a$sd[sig$pup], 
      a$foldchange[sig$pup],qvalues$qvalue.up, fdr.up)

    dimnames(res.up)=list(NULL,c("Row","Gene Name","Gene ID", "Score(d)", "Numerator(r)","Denominator(s+s0)","q-value","localfdr"))

    res.lo=cbind(sig$plo+1,data$genenames[sig$plo],data$geneid[sig$plo],a$tt[sig$plo],a$numer[sig$plo],a$sd[sig$plo], a$foldchange[sig$plo],qvalues$qvalue.lo, fdr.lo)

done=TRUE
  }

  
  

  o1=order(-a$tt[sig$pup])

  res.up=res.up[o1,]

if(!is.null(res.lo)){
  o2=order(a$tt[sig$plo])
  res.lo=res.lo[o2,]
  dimnames(res.lo)=dimnames(res.up)
}


  
color.ind.for.multi=NULL
if(a$resp.type==samr.const.multiclass.response){
   color.ind.for.multi=1*(a$stand.contrasts[sig$pup,]>a$stand.contrasts.95[2]) + (-1)*(a$stand.contrasts[sig$pup,]<a$stand.contrasts.95[1])
}
  return(list(genes.up=res.up,genes.lo=res.lo, color.ind.for.multi=color.ind.for.multi))
}



qvalue.func=function(a,sig, delta.table){

  

  LARGE=10e9

  qvalue.up=rep(NA,length(sig$pup))

  o1=sig$pup

  cutup=delta.table[,8]

  FDR=delta.table[,5]

  

  ii=0

  for(i in o1){

    o= abs(cutup-a$tt[i])

    o[is.na(o)]=LARGE

    oo=(1:length(o))[o==min(o)]

    oo=oo[length(oo)]

    ii=ii+1

    qvalue.up[ii]= FDR[oo]

  }

  

  qvalue.lo=rep(NA,length(sig$plo))

  o2=sig$plo

  cutlo=delta.table[,7]

  

  ii=0

  for(i in o2){

    o= abs(cutlo-a$tt[i])

    o[is.na(o)]=LARGE

    oo=(1:length(o))[o==min(o)]

    oo=oo[length(oo)]

    ii=ii+1

    qvalue.lo[ii]= FDR[oo]

  }

  

  return(list(qvalue.lo=qvalue.lo,qvalue.up=qvalue.up))

}







foldchange.twoclass=function(x,y, logged2){

  

  if(logged2){x=2^x}

  n1 <- sum(y==1)

  n2 <- sum(y==2)

  

  p <- nrow(x)

  m1 <- rowMeans(x[,y==1])

  m2 <- rowMeans(x[,y==2])

  return(m2/m1)

}







foldchange.paired=function(x,y,logged2){

  

  if(logged2){x=2^x}

  

  nc <- ncl(x)/2

  o <- 1:nc

  o1 <- rep(0,ncol(x)/2);o2 <- o1

  for(j in 1:nc){o1[j] <- (1:ncol(x))[y==-o[j]]}

  for(j in 1:nc){o2[j] <- (1:ncol(x))[y==o[j]]}

  

  d <- x[,o2]/x[,o1]

  

  

  m <- rowMeans(d)

  return(m)

}



est.s0<-function(tt,sd,s0.perc=seq(0,.99,len=40)){

  

  ## estimate s0 (exchangeability) factor for denominator.

  ## returns the actual estimate s0 (not a percentile)

  

  a<-cut(sd,breaks=quantile(sd,seq(0,1,len=101)),labels=F)

  a[is.na(a)]<-1

  cv.sd<-rep(0,length(s0.perc))

  

  for(j in 1:length(s0.perc)){

    w<-quantile(sd,s0.perc[j])

    w[j==1]<-0

    tt2<-tt*sd/(sd+w)

    sds<-rep(0,100)

    

    for(i in 1:100){

      sds[i]<-mad(tt2[a==i])

    }

    

    cv.sd[j]<-sqrt(var(sds))/mean(sds)

  }

  o=(1:length(s0.perc))[cv.sd==min(cv.sd)]

  

  s0.hat=min(sd[o])

  return(list(s0.perc=s0.perc,cv.sd=cv.sd, s0.hat= s0.hat))

}





coxscor <- function(x, y, ic, offset = rep(0., length(y))) {

  ## computes cox scor function for rows of nx by n matrix  x

  ## first put everything in time order

  n <- length(y)

  nx <- nrow(x)

  yy <- y + (ic == 0.) * (1e-05)

  otag <- order(yy)

  y <- y[otag]

  ic <- ic[otag]

  x <- x[, otag, drop = F]

  ##compute  unique failure times, d=# of deaths at each failure time, 

  ##dd= expanded version of d to length n, s=sum of covariates at each

  ## failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time

  offset <- offset[otag]

  a <- coxstuff(x, y, ic, offset = offset)

  nf <- a$nf

  fail.times <- a$fail.times

  s <- a$s

  d <- a$d

  dd <- a$dd

  nn <- a$nn

  nno <- a$nno

  w <- rep(0., nx)

  for(i in (1.:nf)) {

    w <- w + s[, i]

    oo<- (1.:n)[y >= fail.times[i]]

    r<-rowSums(x[, oo, drop = F] * exp(offset[oo]))

    w<- w - (d[i]/nno[i])*r 

  }

  return(list(scor = w, coxstuff.obj = a))

}







coxvar <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL){

  ## computes information elements (var) for cox

  ## x is nx by n matrix of expression  values

  nx <- nrow(x)

  n <- length(y)

  yy <- y + (ic == 0.) * (1e-06)

  otag <- order(yy)

  y <- y[otag]

  ic <- ic[otag]

  x <- x[, otag, drop = F]

  offset <- offset[otag]

  if(is.null(coxstuff.obj)) {

    coxstuff.obj <- coxstuff(x, y, ic, offset = offset)

  }

  nf <- coxstuff.obj$nf

  fail.times <- coxstuff.obj$fail.times

  s <- coxstuff.obj$s

  d <- coxstuff.obj$d

  dd <- coxstuff.obj$dd

  nn <- coxstuff.obj$nn

  nno <- coxstuff.obj$nno



  x2<- x^2

  oo <- (1.:n)[y >= fail.times[1] ]

  sx<-(1/nno[1])*rowSums(x[, oo] * exp(offset[oo]))

  s<-(1/nno[1])*rowSums(x2[, oo] * exp(offset[oo]))

  w <-  d[1] * (s - sx * sx)





  for(i in 2.:nf) {

    oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]

    sx<-(1/nno[i])*(nno[i-1]*sx-rowSums(x[, oo,drop=F] * exp(offset[oo])))

    s<-(1/nno[i])*(nno[i-1]*s-rowSums(x2[, oo,drop=F] * exp(offset[oo])))

    w <- w + d[i] * (s - sx * sx)

  }

  return(w)

}









coxstuff<- function(x, y, ic, offset = rep(0., length(y))) {

  fail.times <- unique(y[ic == 1.])

  nf <- length(fail.times)

  n <- length(y)

  nn <- rep(0., nf)

  nno <- rep(0., nf)

  for(i in 1.:nf) {

    nn[i] <- sum(y >= fail.times[i])

    nno[i] <- sum(exp(offset)[y >= fail.times[i]])

  }

  s <- matrix(0., ncol = nf, nrow = nrow(x))

  d <- rep(0., nf)

  ##expand d out to a vector of length n

  for(i in 1.:nf) {

    o <- (1.:n)[(y == fail.times[i]) & (ic == 1.)]

    d[i] <- length(o)

  }

  oo <- match(y, fail.times)

  oo[ic==0]<-NA

  oo[is.na(oo)]<- max(oo[!is.na(oo)])+1

  s<-t(rowsum(t(x),oo))

  if(ncol(s)> nf){s<-s[,-ncol(s)]}

  dd <- rep(0., n)

  for(j in 1.:nf) {

    dd[(y == fail.times[j]) & (ic == 1.)] <- d[j]

  }

  return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))

}







samr.missrate <- function(a, del, delta.table, quant=NULL){

if(is.null(quant)){
if(a$resp.type!=samr.const.multiclass.response){
  quant=c(0, .05,.10, .15, .20, .25, .75, .80, .85,.90,.95, 1.0 )
  }
if(a$resp.type==samr.const.multiclass.response){
  quant=c( .75, .80, .85,.90,.95, 1.0)
  }

}


  ## estimate miss rate from sam object "a"

 o= abs(delta.table[,1]-del)
 oo=(1:nrow(delta.table))[o==min(o)]


cut.lo=delta.table[oo,7]
cut.up=delta.table[oo,8]

ooo=a$tt> cut.lo & a$tt < cut.up
  cuts=quantile(a$tt[ooo],quant)

  ncuts <- length(cuts)



  ngenes <- rep(NA,ncuts)

  ngenes0 <- rep(NA,ncuts)

  ngenes2 <- rep(NA,ncuts)







  missrate <- rep(NA,ncuts)



  nperm=ncol(a$ttstar)





  for(j in 1:(ncuts-1)){

    ngenes2[j] <- sum(a$tt >cuts[j] & a$tt< cuts[j+1])

    ngenes0[j] <- sum(a$ttstar >cuts[j] & a$ttstar< cuts[j+1])/nperm

    missrate[j] <- (ngenes2[j]- a$pi0*ngenes0[j])/ngenes2[j]

    missrate[j] <- max(missrate[j],0)

  }
cuts=round(cuts,3)
res=matrix(NA,ncol=3,nrow=ncuts-1)
missrate=round(missrate,4)
for(i in 1:(ncuts-1)){
res[i,1]=paste(as.character(quant[i]),as.character(quant[i+1]), sep=" -> ")
res[i,2]=paste(as.character(cuts[i]),as.character(cuts[i+1]), sep=" -> ")
res[i,3]=missrate[i]
}

dimnames(res)=list(NULL,c("Quantiles","Cutpoints", "Miss Rate"))

  return(res)

}





varr <- function(x, meanx=NULL){

  n <- ncol(x)

  p <- nrow(x)

  Y <-matrix(1,nrow=n,ncol=1)

  if(is.null(meanx)){   meanx <- rowMeans(x)}

  ans<- rep(1, p)

  xdif <- x - meanx %*% t(Y)

  ans <- (xdif^2) %*% rep(1/(n - 1), n)

  ans <- drop(ans)

  return(ans)

}











samr.options <- list(debug=TRUE, #whether to turn on debugging or not

                    err.file=ifelse(.Platform$OS.type=="windows", "C:/samrtrace.txt", "samrtrace.txt"),

                    image.file=ifelse(.Platform$OS.type=="windows", "C:/samrimage.Rdata", "samrimage.Rdata"))





##

## Our error handler

##



.error.trace <- function() {

 err.message <- geterrmessage()

 if (!is.null(samr.options$image.file)) {

   save.image(samr.options$image.file)

 }

 if (!is.null(samr.options$err.file)) {

   sink(samr.options$err.file)

   print(err.message)

   traceback()

   sink()

 }

 winDialog(type="ok", message=err.message)

}



##

## Upon loading, if we are in a windows environment, we use the windows

## dialog mechanism to display errors. Useful for debugging COM apps

##

.onLoad <- function(lib, pkg) {

 if ( .Platform$OS.type == "windows") {

#    options(error=function() winDialog(type="ok", message=geterrmessage()))

   options(error=samr.xl.error.trace)

 }



}



##

## Upon unload, we set things back the way they were...

##

.onUnload <- function(libpath){

 if ( .Platform$OS.type == "windows") {

   options(error=NULL)

 }

}





samr.xl.build.data <- function(x, y, geneid, genenames, logged2) {

  return(list(x=x,

              y=y,

              geneid=geneid,

              genenames=genenames,

              logged2=logged2))

}






insert.value<-function(vec,newval,pos) {
 if(pos == 1) return(c(newval,vec))
 lvec<-length(vec)
 if(pos > lvec) return(c(vec,newval))
 return(c(vec[1:pos-1],newval,vec[pos:lvec]))
}

permute<-function(elem) {
 if(!missing(elem)) {
  if(length(elem) == 2) return(matrix(c(elem,elem[2],elem[1]),nrow=2))
  last.matrix<-permute(elem[-1])
  dim.last<-dim(last.matrix)
  new.matrix<-matrix(0,nrow=dim.last[1]*(dim.last[2]+1),ncol=dim.last[2]+1)
  for(row in 1:(dim.last[1])) {
   for(col in 1:(dim.last[2]+1))
    new.matrix[row+(col-1)*dim.last[1],]<-insert.value(last.matrix[row,],elem[1]
,col)
  }
  return(new.matrix)
 }
 else cat("Usage: permute(elem)\n\twhere elem is a vector\n")
} 


integer.base.b <-
function(x, b=2){
        xi <- as.integer(x)
        if(xi==0){return(0)}

        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b 
}


  
compute.block.perms=function(y,blocky,nperms){

# y are the data (eg class label 1 vs 2; or -1,1, -2,2 for paired data)
# blocky are the block labels (abs(y) for paired daatr)
    ny=length(y)
     nblocks=length(unique(blocky))
     tab=table(blocky)
     total.nperms=prod(factorial(tab))

# block.perms is a list of all possible permutations
     block.perms=vector("list",nblocks)

# first enumerate all perms, when possible
     if(total.nperms<=nperms){
        all.perms.flag=1
        nperms.act=total.nperms
         for(i in 1:nblocks){
            block.perms[[i]]=permute(y[blocky==i])
          }
        
      kk=0:(factorial(max(tab))^nblocks-1)

#the rows of  the matrix outerm runs through the "outer product"
# first we assume that all blocks have max(tab) members; then we remove rows of outerm that
#  are illegal (ie when a block has fewer members)

      outerm=matrix(0,nrow=length(kk),ncol=nblocks)
     for(i in 1:length(kk)){
	     kkkk=integer.base.b(kk[i],b=factorial(max(tab)))

          if(length(kkkk)>nblocks){kkkk=kkkk[(length(kkkk)-nblocks+1):length(kkkk)]}
       
       outerm[i,(nblocks-length(kkkk)+1):nblocks]=kkkk
       }
       outerm=outerm+1

# now remove rows that are illegal perms
ind=rep(TRUE,nrow(outerm))
   for(j in 1:ncol(outerm)){
      ind=ind&outerm[,j]<=factorial(tab[j])
   }

outerm=outerm[ind,]
# finally, construct permutation matrix from outer product
    permsy=matrix(NA,nrow=total.nperms,ncol=ny)
    for(i in 1:total.nperms){
       junk=NULL
       for(j in 1:nblocks){
         junk=c(junk,block.perms[[j]][outerm[i,j],])
       }
     permsy[i,]=junk
    }}
  

# next handle case when there are too many perms to enumerate
     if(total.nperms>nperms){
       all.perms.flag=0
       nperms.act=nperms
  permsy=NULL

          block.perms=vector("list",nblocks)
         for(j in 1:nblocks){
             block.perms[[j]]=permute(y[blocky==j])
           }
         for(j in 1:nblocks){
             o=sample(1:nrow(block.perms[[j]]),size=nperms,replace=T)
             permsy=cbind(permsy,block.perms[[j]][o,])
          }
   
            
     }




return(list(permsy=permsy,all.perms.flag=all.perms.flag, nperms.act=nperms.act))
  }




getperms=function(y, nperms){
total.perms=factorial(length(y))
   if(total.perms<= nperms){
         perms=permute(1:length(y))
         all.perms.flag=1
         nperms.act=total.perms
}

   if(total.perms> nperms){
         perms=matrix(NA,nrow=nperms,ncol=length(y))
         for(i in 1:nperms){
            perms[i,]=sample(1:length(y), size=length(y))
          }
         all.perms.flag=0
         nperms.act=nperms
       }
 

return(list(perms=perms, all.perms.flag=all.perms.flag, nperms.act=nperms.act))
}

parse.block.labels.for.2classes=function(y){
#this only  works for 2 class case- having form jBlockn, where j=1 or 2
n=length(y)
y.act=rep(NA,n)
 blocky=rep(NA,n)
for(i in 1:n){
   blocky[i]=as.numeric(substring(y[i],7,nchar(y[i])))
   y.act[i]=as.numeric(substring(y[i],1,1))
}
return(list(y.act=y.act,blocky=blocky))
}

parse.time.labels.and.summarize.data=function(x,y,resp.type, time.summary.type){
# parse time labels, and summarize time data for each person, via a slope or area
# does some error checking too
  
  n=length(y)

  

  y.act=rep(NA,n)
  timey=rep(NA,n)
  person.id=rep(NA,n)
  k=1
  end.flag=FALSE
  person.id[1]=1
  if(substring(y[1],nchar(y[1])-4,nchar(y[1]))!="Start"){ 
    stop("Error in format of  time course data")
  }

  for(i in 1:n){
    cat(i)
    j=1
    while(substring(y[i],j,j)!="T"){j=j+1}
    end.of.y= j-1
    y.act[i]=as.numeric(substring(y[i],1,end.of.y))
    timey[i]=substring(y[i],end.of.y +5 ,nchar(y[i]))
    if(nchar(timey[i])>3 & substring(timey[i],nchar(timey[i])-2,nchar(timey[i]))=="End"){
      end.flag=TRUE
      timey[i]=substring(timey[i],1,nchar(timey[i])-3)
    }
    if(nchar(timey[i])>3 & substring(timey[i],nchar(timey[i])-4,nchar(timey[i]))=="Start"){  timey[i]=substring(timey[i],1,nchar(timey[i])-5)}

    if( i<n & !end.flag){ person.id[i+1]=k}
    if( i<n & end.flag){k=k+1; person.id[i+1]=k}

    end.flag=FALSE                
  }
  timey=as.numeric(timey)

                                        

 # do a check that the format was correct
  tt=table(person.id, y.act)
  junk=function(x){sum(x!=0)}
  if(sum(apply(tt,1,junk)!=1)>0){
    stop("Error in format of  time course data")
  }
  
  npeople=length(unique(person.id))

  newx=matrix(NA,nrow=nrow(x),ncol=npeople)

  for(j in 1:npeople){
    jj=person.id==j
    tim=timey[jj]
    xc=t(scale(t(x[,jj]),center=TRUE,scale=FALSE))
    if(time.summary.type=="slope"){
      newx[,j]=quantitative.func(xc,tim-mean(tim))$tt
    }
    if(time.summary.type=="signed.area"){
      newx[,j]=timearea.func(x[,jj],tim)$numer
    }
  }

  y.unique=y.act[!duplicated(person.id)]
  return(list(y=y.unique,  x=newx))
}


check.format=function(y, resp.type, status=NULL){

                                        # here i do some format checks for the input data$y
                                        # note that checks for time course data are done in the parse function for time course;
                                        #  we then check the output from the parser in this function
  

 



  if(resp.type==samr.const.twoclass.unpaired.response | resp.type==samr.const.twoclass.unpaired.timecourse.response){
    if(sum(y==1)+sum(y==2) !=length(y)){
      stop("Error in input response data: values must be 1 or 2")
    }
  }

  if(resp.type==samr.const.twoclass.paired.response | resp.type==samr.const.twoclass.paired.timecourse.response ){
    if(sum(y)!=0){
      stop("Error in input response data: values must be -1, 1, -2, 2, etc")
    }
    if(sum(table(y[y>0])!=abs(table(y[y<0]))) ){
      stop("Error in input response data: values must be -1, 1, -2, 2, etc")
    }
  }



  if(resp.type==samr.const.oneclass.response | resp.type==samr.const.oneclass.timecourse.response){
    if(sum(y==1) !=length(y)){
      stop("Error in input response data: values must all be 1")
    }
  }
  

 

  if(resp.type==samr.const.multiclass.response){
    tt=table(y)
    nc=length(tt)
    if(sum(y<=nc & y> 0) <length(y)){
      stop("Error in input response data: values must be 1,2, ... number of classes")
    }

    for(k in 1:nc){
      if(sum(y==k)<2){
        stop("Error in input response data: there must be >1 sample per class")
      }}
  }


  if(resp.type==samr.const.quantitative.response){
    if(!is.numeric(y)){
      stop("Error in input response data: values must be numeric")
    }
  }

    if(resp.type==samr.const.survival.response){
      if(is.null(status)){
          stop("Error in input response data: error in censoring indicator")
        }
      if(!is.numeric(y) | sum(y<0) >0){
        stop("Error in input response data: survival times  must be numeric and nonnegative")
        if(sum(status==0) +sum(status==1) !=length(status)){
          stop("Error in input response data: censoring indicators must be 0 (censored) or 1 (failed)")

        }
      }
      if(sum(status==1)<1){
          stop("Error in input response data: there are no uncensored observations")
        }
    }

return()
}

