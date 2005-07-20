

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


  sd <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m2) + (n1-1) * varr(x[, y==1], meanx=m1) )*(1/n1+1/n2)/(n1+n2-2) )

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
  x <- x*matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=TRUE)
  m <- rowMeans(x)
  sd <- sqrt( varr(x, meanx=m)/n )
  dif.obs <- m/(sd + s0)
  return(list(tt=dif.obs, numer=m,sd=sd))

}

patterndiscovery.func=function(x,s0=0, eigengene.number=1){
  a=mysvd(x, n.components=eigengene.number)
v=a$v[,eigengene.number]

# here we try to guess the most interpretable orientation for the eigengene
om=abs(a$u[, eigengene.number]) >quantile(abs(a$u[, eigengene.number]),.95)
if(median(a$u[om,eigengene.number])<0){
 v=-1.0*v
}

aa=quantitative.func(x,v,s0=s0)
eigengene=cbind(1:nrow(a$v),v)
dimnames(eigengene)=list(NULL,c("sample number","value"))
  return(list(tt=aa$tt, numer=aa$numer, sd=aa$sd, eigengene=eigengene))


}




paired.ttest.func <- function(x,y,s0=0,useden=TRUE){
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



cox.func <- function(x,y,censoring.status,s0=0){
  scor <- coxscor(x,y, censoring.status)$scor
  sd <- sqrt(coxvar(x,y, censoring.status))
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
  scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=TRUE)*mm*mm,1,sum)))

  sd <- sqrt(rowSums(v)*(1/sum(nn-1))*sum(1/nn))
  tt <- scor/(sd+s0)
  mm.stand=t(scale(t(mm),center=FALSE,scale=sd))
  return(list(tt=tt, numer=scor, sd=sd,stand.contrasts=mm.stand))

}





#quantitative.func <- function(x,y,s0=0){
#  yy <- y-mean(y)
#  temp <- x%*%yy
#mx=rowMeans(x)
#sxx <-rowSums( (x-mx%*%t(rep(1,ncol(x))))^2 )
#
#  scor <- temp/sxx
#  b0hat <- mean(y)-scor*mx
#  yhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+x*matrix(scor,nrow=nrow(x),ncol=ncol(x))
#  ty <- matrix(y,nrow=nrow(yhat),ncol=ncol(yhat),byrow=TRUE)
#  sigma <- sqrt(rowSums((ty-yhat)^2)/(ncol(yhat)-2))
#  sd <- sigma/sqrt(sxx)
#  tt <- scor/(sd+s0)
#  return(list(tt=tt, numer=scor, sd=sd))
#
#}
quantitative.func  <-
function(x,y,s0=0){

# regression of x on y

my=mean(y)
  yy <- y-my
  temp <- x%*%yy
mx=rowMeans(x)
syy= sum(yy^2)

  scor <- temp/syy
  b0hat <- mx-scor*my
  xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+y*matrix(scor,nrow=nrow(x),ncol=ncol(x))
  sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
  sd <- sigma/sqrt(syy)
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



detec.slab <- function(samr.obj, del, min.foldchange) {

  ## find genes above and below the slab of half-width del

 # this calculation is tricky- for consistency,  the slab condition picks
# all genes that are beyond the first departure from the slab
# then the fold change condition is applied (if applicable)

  n <- length(samr.obj$tt)
  tt <- samr.obj$tt
  evo <- samr.obj$evo
  numer <- samr.obj$tt*(samr.obj$sd+samr.obj$s0)
  tag <- order(tt)

  pup <- NULL
  foldchange.cond.up=rep(T,length(evo))
  foldchange.cond.lo=rep(T,length(evo))

  

  if(!is.null(samr.obj$foldchange[1]) & (min.foldchange>0)){

    foldchange.cond.up= samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo= samr.obj$foldchange <= 1/min.foldchange

  }

 o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0  ]
  if(length(o1) > 0) {
    o1 <- o1[1]
    o11 <- o1:n
    o111 <- rep(F, n)
    o111[tag][o11] <- T
    pup <- (1:n)[o111 & foldchange.cond.up]
  }
  plow <- NULL
  o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0  ]
  if(length(o2) > 0) {
    o2 <- o2[length(o2)]
    o22 <- 1:o2
    o222 <- rep(F, n)
    o222[tag][o22] <- T
    plow <- (1:n)[o222 & foldchange.cond.lo]
  }


  

  return(list(plow=plow, pup=pup))

}



sumlengths <- function(aa) {

  length(aa$pl) + length(aa$pu)

}



samr.compute.delta.table <- function(samr.obj, min.foldchange=0, dels=NULL, nvals=50) {
                                        # computes delta table, starting with samr object "a", for nvals values of delta


 lmax=sqrt(max(abs(sort(samr.obj$tt)-samr.obj$evo)))


if(is.null(dels)){dels=(seq(0,lmax, length=nvals)^2)}

  col=matrix(1,nrow=length(samr.obj$evo),ncol=nvals)

  ttstar0 <- samr.obj$ttstar0
  tt <- samr.obj$tt
  n <- samr.obj$n
  evo <- samr.obj$evo
  nsim <- ncol(ttstar0)
  res1 <- NULL



    foldchange.cond.up=matrix(T,nrow=nrow(samr.obj$ttstar),ncol=ncol(samr.obj$ttstar))
    foldchange.cond.lo=matrix(T,nrow=nrow(samr.obj$ttstar),ncol=ncol(samr.obj$ttstar))

    if(!is.null(samr.obj$foldchange[1]) & (min.foldchange>0)){

      foldchange.cond.up= samr.obj$foldchange.star >= min.foldchange

      foldchange.cond.lo= samr.obj$foldchange.star <= 1/min.foldchange
    }

  

cutup=rep(NA, length(dels))
cutlow=rep(NA, length(dels))
g2=rep(NA, length(dels))
errup=matrix(NA,ncol=length(dels),nrow=ncol(samr.obj$ttstar0))
errlow=matrix(NA,ncol=length(dels),nrow=ncol(samr.obj$ttstar0))

  for(ii in 1:length(dels)) {
    
cat(ii,fill=TRUE)     
     ttt <- detec.slab(samr.obj,dels[ii], min.foldchange)
    cutup[ii] <- 10e9
    if(length(ttt$pup>0)){ cutup[ii] <- min(samr.obj$tt[ttt$pup])}
    cutlow[ii] <- -10e9
    if(length(ttt$plow)>0){cutlow[ii] <- max(samr.obj$tt[ttt$plow])}
   g2[ii]=sumlengths(ttt)

errup[,ii]=colSums(samr.obj$ttstar0> cutup[ii] & foldchange.cond.up)
errlow[,ii]=colSums(samr.obj$ttstar0< cutlow[ii] & foldchange.cond.lo)
}

  
    s <- sqrt(apply(errup,2,var)/nsim + apply(errlow,2,var)/nsim)
    gmed <- apply(errup + errlow,2,median)
    g90=apply(errup + errlow,2,quantile, .90)


    res1 <-  cbind(samr.obj$pi0*gmed, samr.obj$pi0*g90, g2, samr.obj$pi0*gmed/g2, samr.obj$pi0*g90/g2, cutlow,cutup)

       
  res1 <- cbind(dels, res1)

# remove rows with #called=0
#om=res1[,4]==0
#res1=res1[!om,,drop=F]

# remove duplicate rows with same # of genes called

#omm=!duplicated(res1[,4])
#res1=res1[omm,,drop=F]


  dimnames(res1) <- list(NULL, c("delta", "# med false pos",
                                 "90th perc false pos", "# called", "median FDR", "90th perc FDR","cutlo","cuthi"))

  return(res1)
}



detec.horiz <- function(samr.obj,cutlow,cutup, min.foldchange){

  ## find genes above or below horizontal cutpoints

  dobs <- samr.obj$tt
  n <- length(dobs)
  foldchange.cond.up=rep(T,n)
  foldchange.cond.lo=rep(T,n)

  

  if(!is.null(samr.obj$foldchange[1]) & (min.foldchange>0)){
    foldchange.cond.up= samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo= samr.obj$foldchange <= 1/min.foldchange
  }

  

  pup <- (1:n)[dobs> cutup & foldchange.cond.up]
  plow <- (1:n)[dobs< cutlow & foldchange.cond.lo]

  

  return(list(plow=plow,pup=pup))

  

}





samr.plot <- function(samr.obj, del, min.foldchange=0) {

## make observed-expected plot
## takes foldchange into account too



  LARGE=10e9

  b <-  detec.slab(samr.obj, del, min.foldchange)

  bb <- c(b$pup,b$plow)
  b1= LARGE
  b0=-LARGE

  

  if(!is.null(b$pup)){b1 <- min(samr.obj$tt[b$pup])}
  if(!is.null(b$plow)){b0 <- max(samr.obj$tt[b$plow])}

  c1 <- (1:samr.obj$n)[sort(samr.obj$tt)>=b1]
  c0 <- (1:samr.obj$n)[sort(samr.obj$tt)<=b0]
  c2 <- c(c0,c1)

  

  foldchange.cond.up=rep(T,length(samr.obj$evo))
  foldchange.cond.lo=rep(T,length(samr.obj$evo))

  

  if(!is.null(samr.obj$foldchange[1]) & (min.foldchange>0)){
    foldchange.cond.up= samr.obj$foldchange >= min.foldchange
     foldchange.cond.lo= samr.obj$foldchange <= 1/min.foldchange
  }

  

  col=rep(1,length(samr.obj$evo))
  col[b$plow]=3
  col[b$pup]=2

  if(!is.null(samr.obj$foldchange[1]) & (min.foldchange>0)){
    col[!foldchange.cond.lo & !foldchange.cond.up]=1
  }

  col.ordered=col[order(samr.obj$tt)]

    ylims <- range(samr.obj$tt)
    xlims <- range(samr.obj$evo)

    plot(samr.obj$evo,sort(samr.obj$tt),xlab="expected score", ylab="observed score",ylim=ylims, 
         xlim=xlims, type="n")

    points(samr.obj$evo,sort(samr.obj$tt),col=col.ordered)

       

    abline(0,1)
    abline(del,1,lty=2)
    abline(-del,1,lty=2)

}



localfdr <- function(samr.obj,  min.foldchange, perc=.01, df=10) {

  ## estimates local fdr at score "d", using SAM object "samr.obj"
  ## "d" can be a vector of d scores
  ## returns estimate of symmetric fdr  as a percentage

# this version uses a 1% symmetric window, and does not estimate fdr in
# windows  having fewer than 100 genes


  ## to use: first run SAM in Splus, and then pass the resulting fit object to
  ## localfdr
  ## NOTE: at most 20 of the perms are used to estimate the fdr (for speed sake)


# I try two window shapes: symmetric and an assymetric one
# currently I use the symetric window to estimate the  local fdr

ngenes=length(samr.obj$tt)
mingenes=50

# perc is increased, in order to get at least mingenes in a window
perc=max(perc,mingenes/length(samr.obj$tt))

 nperms.to.use=min(20,ncol(samr.obj$ttstar))
nperms=ncol(samr.obj$ttstar)

d=seq(sort(samr.obj$tt)[1], sort(samr.obj$tt)[ngenes], length=100)
  ndscore <- length(d)
  dvector <- rep(NA,ndscore)

  ind.foldchange=rep(T,length(samr.obj$tt))
  if(!is.null(samr.obj$foldchange[1]) & min.foldchange>0){
    ind.foldchange= (samr.obj$foldchange>=  min.foldchange) | (samr.obj$foldchange<=  min.foldchange)
  }


fdr.temp=function(temp, dlow, dup, pi0, ind.foldchange){
return(sum(pi0*(temp>=dlow & temp<=dup &ind.foldchange)))}


  for(i in 1:ndscore)

    {
      pi0<-samr.obj$pi0
      r <- sum(samr.obj$tt<d[i])
      r22 <-round(max(r-length(samr.obj$tt)*perc/2, 1))
      dlow.sym <- sort(samr.obj$tt)[r22]

      if(d[i]<0)
        {
          r2 <- max(r-length(samr.obj$tt)*perc/2, 1)
          r22= min(r+length(samr.obj$tt)*perc/2, length(samr.obj$tt))

          dlow <- sort(samr.obj$tt)[r2]
          dup=sort(samr.obj$tt)[r22]

        }

      r22 <- min(r+length(samr.obj$tt)*perc/2, length(samr.obj$tt))
      dup.sym <- sort(samr.obj$tt)[r22]

      if(d[i]>0)
        {
          r2 <- min(r+length(samr.obj$tt)*perc/2, length(samr.obj$tt))
          r22 <- max(r-length(samr.obj$tt)*perc/2, 1)
          dup <- sort(samr.obj$tt)[r2]
          dlow <- sort(samr.obj$tt)[r22]

        }
      o <- samr.obj$tt>=dlow & samr.obj$tt<= dup & ind.foldchange
      oo <- samr.obj$tt>=dlow.sym & samr.obj$tt<= dup.sym & ind.foldchange


      nsim <- ncol(samr.obj$ttstar)
      fdr <- rep(NA,nsim)
      fdr2 <- fdr

          if(!is.null(samr.obj$foldchange[1]) &  min.foldchange>0){
             temp=as.vector(samr.obj$foldchange.star[,1:nperms.to.use])
            ind.foldchange=(temp >=  min.foldchange) | (temp <=  min.foldchange)
          }

         temp=samr.obj$ttstar0[,sample(1:nperms, size=nperms.to.use)]

          fdr <-median(apply(temp,2,fdr.temp,dlow, dup, pi0, ind.foldchange))
          fdr.sym <-median(apply(temp,2,fdr.temp,dlow.sym, dup.sym, pi0, ind.foldchange))

      fdr <- 100*fdr/sum(o)
      fdr.sym <- 100*fdr.sym/sum(oo)



      dlow.sym <- dlow.sym
      dup.sym <-dup.sym
      dlow <- dlow
      dup <- dup

      dvector[i] <- fdr.sym

    }

om=!is.na(dvector) & (dvector!=Inf)
aa=smooth.spline(d[om], dvector[om], df=df)
  return(list(smooth.object=aa, perc=perc,df=df))

}



predictlocalfdr= function(smooth.object, d){
 yhat=predict(smooth.object, d)$y
yhat=pmin(yhat,100)
yhat=pmax(yhat,0)

  return(yhat)

}



samr.compute.siggenes.table=function(samr.obj,del, data, delta.table,  min.foldchange=0, all.genes=FALSE){



  ## computes significant genes table, starting with samr object "a" and "delta.table"

  ##  for a  **single** value del
  ## if all.genes is true, all genes are printed (and value of del is ignored)


  if(!all.genes){
    sig=detec.slab(samr.obj, del, min.foldchange)
  }

  if(all.genes){
    p=length(samr.obj$tt)
    pup=(1:p)[samr.obj$tt>=0]
    plo=(1:p)[samr.obj$tt<0]
    sig=list(pup=pup, plo=plo)
  }

  aa=localfdr(samr.obj, min.foldchange)

  if(length(sig$pup)>0){fdr.up=predictlocalfdr(aa$smooth.object,samr.obj$tt[sig$pup])}

  if(length(sig$plo)>0){ fdr.lo=predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$plo])}

qvalues=NULL
if(length(sig$pup)>0 | length(sig$plo)>0){
  qvalues=qvalue.func(samr.obj,sig, delta.table)
}

  
  res.up=NULL
  res.lo=NULL

  done=FALSE
     # two class unpaired or paired (folchange is reported) 

  if((samr.obj$resp.type==samr.const.twoclass.unpaired.response | samr.obj$resp.type==samr.const.twoclass.paired.response) ){

    if(!is.null(sig$pup)){
      res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],samr.obj$tt[sig$pup],samr.obj$numer[sig$pup],samr.obj$sd[sig$pup],samr.obj$foldchange[sig$pup],qvalues$qvalue.up, fdr.up)

      dimnames(res.up)=list(NULL,c("Row","Gene ID","Gene Name", "Score(d)", "Numerator(r)","Denominator(s+s0)", "Fold Change", "q-value(%)","localfdr(%)"))
    }

    if(!is.null(sig$plo)){
      res.lo=cbind(sig$plo+1,data$genenames[sig$plo],data$geneid[sig$plo],samr.obj$tt[sig$plo],samr.obj$numer[sig$plo],samr.obj$sd[sig$plo],   samr.obj$foldchange[sig$plo],qvalues$qvalue.lo, fdr.lo)

      dimnames(res.lo)=list(NULL,c("Row","Gene ID","Gene Name", "Score(d)", "Numerator(r)","Denominator(s+s0)", "Fold Change", "q-value(%)","localfdr(%)"))


    }

    done=TRUE
  }

# multiclass
  if(samr.obj$resp.type==samr.const.multiclass.response){

    if(!is.null(sig$pup)){
      res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],samr.obj$tt[sig$pup],samr.obj$numer[sig$pup],samr.obj$sd[sig$pup],samr.obj$stand.contrasts[sig$pup,],qvalues$qvalue.up, fdr.up)

      collabs.contrast=paste("contrast-",as.character(1:ncol(samr.obj$stand.contrasts)),sep="")
      dimnames(res.up)=list(NULL,c("Row","Gene ID","Gene Name", "Score(d)", "Numerator(r)","Denominator(s+s0)",collabs.contrast,  "q-value(%)","localfdr(%)"))
    }

    res.lo=NULL
    done=TRUE
  }





 #all other cases


  if(!done){
    if(!is.null(sig$pup)){
      res.up=cbind(sig$pup+1,data$genenames[sig$pup],data$geneid[sig$pup],samr.obj$tt[sig$pup],samr.obj$numer[sig$pup],samr.obj$sd[sig$pup], 
        samr.obj$foldchange[sig$pup],qvalues$qvalue.up, fdr.up)

      dimnames(res.up)=list(NULL,c("Row","Gene ID","Gene Name", "Score(d)", "Numerator(r)","Denominator(s+s0)","q-value(%)","localfdr(%)"))
    }

    if(!is.null(sig$plo)){
      res.lo=cbind(sig$plo+1,data$genenames[sig$plo],data$geneid[sig$plo],samr.obj$tt[sig$plo],samr.obj$numer[sig$plo],samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo],qvalues$qvalue.lo, fdr.lo)

      dimnames(res.lo)=list(NULL,c("Row","Gene ID","Gene Name", "Score(d)", "Numerator(r)","Denominator(s+s0)","q-value(%)","localfdr(%)"))
    }
    done=TRUE
  }


  if(!is.null(res.up)){ 
    o1=order(-samr.obj$tt[sig$pup])
    res.up=res.up[o1,,drop=F]
  }

  if(!is.null(res.lo)){
    o2=order(samr.obj$tt[sig$plo])
    res.lo=res.lo[o2,,drop=F]
  }


  
  color.ind.for.multi=NULL
  if(samr.obj$resp.type==samr.const.multiclass.response & !is.null(sig$pup)){
    color.ind.for.multi=1*(samr.obj$stand.contrasts[sig$pup,]>samr.obj$stand.contrasts.95[2]) + (-1)*(samr.obj$stand.contrasts[sig$pup,]<samr.obj$stand.contrasts.95[1])
  }

  ngenes.up=nrow(res.up)
  if(is.null(ngenes.up)){ngenes.up=0}
   ngenes.lo=nrow(res.lo)
  if(is.null(ngenes.lo)){ngenes.lo=0}
  
  return(list(genes.up=res.up,genes.lo=res.lo, color.ind.for.multi=color.ind.for.multi, ngenes.up=ngenes.up, ngenes.lo=ngenes.lo))
}



qvalue.func=function(samr.obj,sig, delta.table){

# returns q-value as a percentage (out of 100)
  

  LARGE=10e9
  qvalue.up=rep(NA,length(sig$pup))
  o1=sig$pup
  cutup=delta.table[,8]
  FDR=delta.table[,5]

  ii=0

  for(i in o1){

    o= abs(cutup-samr.obj$tt[i])
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
    o= abs(cutlo-samr.obj$tt[i])
    o[is.na(o)]=LARGE
    oo=(1:length(o))[o==min(o)]
    oo=oo[length(oo)]
    ii=ii+1
    qvalue.lo[ii]= FDR[oo]

  }

# any qvalues that are missing, are set to 1 (the highest value)

qvalue.lo[is.na(qvalue.lo)]=1
qvalue.up[is.na(qvalue.up)]=1

# ensure that each qvalue vector is monotone non-increasing

o1=order(samr.obj$tt[sig$plo])
qv1=qvalue.lo[o1]
qv11=qv1

if(length(qv1)>1){
for(i in 2:length(qv1)){
  if(qv11[i]<qv11[i-1]){qv11[i]=qv11[i-1]}
}
qv111=qv11
qv111[o1]=qv11
}
else{qv111=qv1}

o2=order(samr.obj$tt[sig$pup])
qv2=qvalue.up[o2]
qv22=qv2
if(length(qv2)>1){
for(i in 2:length(qv2)){
  if(qv22[i]>qv22[i-1]){qv22[i]=qv22[i-1]}
}
qv222=qv22
qv222[o2]=qv22
}
else{qv222=qv2}


  return(list(qvalue.lo=100*qv111,qvalue.up=100*qv222))

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

  

  nc <- ncol(x)/2
  o <- 1:nc
  o1 <- rep(0,ncol(x)/2);o2 <- o1
  for(j in 1:nc){o1[j] <- (1:ncol(x))[y==-o[j]]}
  for(j in 1:nc){o2[j] <- (1:ncol(x))[y==o[j]]}


  d <- x[,o2]/x[,o1]

  m <- rowMeans(d)

  return(m)

}



est.s0<-function(tt,sd,s0.perc=seq(0,1, by=.05)){


  ## estimate s0 (exchangeability) factor for denominator.
  ## returns the actual estimate s0 (not a percentile)

  a<-cut(sd,breaks=quantile(sd,seq(0,1,len=101)),labels=F)
  a[is.na(a)]<-1
  cv.sd<-rep(0,length(s0.perc))

  for(j in 1:length(s0.perc)){
    w<-quantile(sd,s0.perc[j])
    w[j==1]<-0
    tt2<-tt*sd/(sd+w)
    tt2[tt2==Inf]=NA
    sds<-rep(0,100)

    for(i in 1:100){
      sds[i]<-mad(tt2[a==i], na.rm=TRUE)
    }

    cv.sd[j]<-sqrt(var(sds))/mean(sds)
  }

  o=(1:length(s0.perc))[cv.sd==min(cv.sd)]

# we don;t allow taking s0.aht to be  0th percentile when min sd is 0
  s0.hat=quantile(sd[sd!=0],s0.perc[o])

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





samr.missrate <- function(samr.obj, del, delta.table, quant=NULL){

# returns miss rate as a percentage

if(is.null(quant)){
if(samr.obj$resp.type!=samr.const.multiclass.response){
  quant=c(0, .05,.10, .15, .20, .25, .75, .80, .85,.90,.95, 1.0 )
  }
if(samr.obj$resp.type==samr.const.multiclass.response){
  quant=c( .75, .80, .85,.90,.95, 1.0)
  }

}


  ## estimate miss rate from sam object "a"

 o= abs(delta.table[,1]-del)
 oo=(1:nrow(delta.table))[o==min(o)]


cut.lo=delta.table[oo,7]
cut.up=delta.table[oo,8]

ooo=samr.obj$tt> cut.lo & samr.obj$tt < cut.up
  cuts=quantile(samr.obj$tt[ooo],quant)

  ncuts <- length(cuts)


  ngenes <- rep(NA,ncuts)
  ngenes0 <- rep(NA,ncuts)
  ngenes2 <- rep(NA,ncuts)


  missrate <- rep(NA,ncuts)


  nperm=ncol(samr.obj$ttstar)


  for(j in 1:(ncuts-1)){
    ngenes2[j] <- sum(samr.obj$tt >cuts[j] & samr.obj$tt< cuts[j+1])
    ngenes0[j] <- sum(samr.obj$ttstar >cuts[j] & samr.obj$ttstar< cuts[j+1])/nperm
    missrate[j] <- (ngenes2[j]- samr.obj$pi0*ngenes0[j])/ngenes2[j]
    missrate[j] <- max(missrate[j],0)
  }
cuts=round(cuts,3)
res=matrix(NA,ncol=3,nrow=ncuts-1)
missrate=round(missrate,4)
for(i in 1:(ncuts-1)){
res[i,1]=paste(as.character(quant[i]),as.character(quant[i+1]), sep=" -> ")
res[i,2]=paste(as.character(cuts[i]),as.character(cuts[i+1]), sep=" -> ")
res[i,3]=100*missrate[i]
}

dimnames(res)=list(NULL,c("Quantiles","Cutpoints", "Miss Rate(%)"))

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
# generates all perms of the vector elem
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

sample.perms <-function(elem, nperms) {
# randomly generates  nperms of the vector elem
 res=permute.rows(matrix(elem,nrow=nperms,ncol=length(elem),byrow=T))
  return(res)
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
             block.perms[[j]]=sample.perms(y[blocky==j], nperms=nperms)
           }
         for(j in 1:nblocks){
             permsy=cbind(permsy,block.perms[[j]])
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

  
last5char=rep(NA,n)
last3char=rep(NA,n)
for( i in 1:n){
 last3char[i]=substring(y[i],nchar(y[i])-2, nchar(y[i]))
 last5char[i]=substring(y[i],nchar(y[i])-4, nchar(y[i]))
}

if(sum(last3char=="End")!= sum(last5char=="Start")){
 stop("Error in format of  time course data: a Start or End tag is missing")
}

  y.act=rep(NA,n)
  timey=rep(NA,n)
  person.id=rep(NA,n)
  k=1
  end.flag=FALSE
  person.id[1]=1
  if(substring(y[1],nchar(y[1])-4,nchar(y[1]))!="Start"){ 
    stop("Error in format of  time course data: first cell should have a Start tag")
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
   num=(1:nrow(tt))[apply(tt,1,junk)>1]
    stop(paste("Error in format of  time course data, timecourse #",as.character(num)) )
  }
  
  npeople=length(unique(person.id))

  newx=matrix(NA,nrow=nrow(x),ncol=npeople)

  for(j in 1:npeople){
    jj=person.id==j
    tim=timey[jj]
    xc=t(scale(t(x[,jj]),center=TRUE,scale=FALSE))
    if(time.summary.type=="slope"){
      newx[,j]=quantitative.func(xc,tim-mean(tim))$numer
    }
    if(time.summary.type=="signed.area"){
      newx[,j]=timearea.func(x[,jj],tim)$numer
    }
  }

  y.unique=y.act[!duplicated(person.id)]
  return(list(y=y.unique,  x=newx))
}


check.format=function(y, resp.type, censoring.status=NULL){

# here i do some format checks for the input data$y
# note that checks for time course data are done in the parse function for time course;
#  we then check the output from the parser in this function
  

 



  if(resp.type==samr.const.twoclass.unpaired.response | resp.type==samr.const.twoclass.unpaired.timecourse.response){
    if(sum(y==1)+sum(y==2) !=length(y)){
      stop(paste("Error in input response data: response type ",
resp.type, " specified; values must be 1 or 2"))
    }
  }

  if(resp.type==samr.const.twoclass.paired.response | resp.type==samr.const.twoclass.paired.timecourse.response ){
    if(sum(y)!=0){
      stop(paste("Error in input response data: response type ",
resp.type, " specified; values must be -1, 1, -2, 2, etc"))
    }
    if(sum(table(y[y>0])!=abs(table(y[y<0]))) ){
      stop(paste("Error in input response data:  response type ",
resp.type, " specified; values must be -1, 1, -2, 2, etc"))
    }
  }



  if(resp.type==samr.const.oneclass.response | resp.type==samr.const.oneclass.timecourse.response){
    if(sum(y==1) !=length(y)){
      stop(paste("Error in input response data: response type ",
resp.type, " specified;  values must all be 1"))
    }
  }
  

 

  if(resp.type==samr.const.multiclass.response){
    tt=table(y)
    nc=length(tt)
    if(sum(y<=nc & y> 0) <length(y)){
      stop(paste("Error in input response data: response type ",
resp.type, " specified; values must be 1,2, ... number of classes"))
    }

    for(k in 1:nc){
      if(sum(y==k)<2){
        stop(paste("Error in input response data: response type ",
resp.type, " specified; there must be >1 sample per class"))
      }}
  }


  if(resp.type==samr.const.quantitative.response){
    if(!is.numeric(y)){
      stop(paste("Error in input response data: response type", resp.type, " specified; values must be numeric"))
    }
  }

    if(resp.type==samr.const.survival.response){
      if(is.null(censoring.status)){
          stop(paste("Error in input response data: response type ",
resp.type, " specified; error in censoring indicator"))
        }
      if(!is.numeric(y) | sum(y<0) >0){
        stop(paste("Error in input response data:  response type ",resp.type, " specified; survival times  must be numeric and nonnegative"))
        if(sum(censoring.status==0) +sum(censoring.status==1) !=length(censoring.status)){
          stop(paste("Error in input response data: response type ",
resp.type, " specified; censoring indicators must be 0 (censored) or 1 (failed)"))

        }
      }
      if(sum(censoring.status==1)<1){
          stop(paste("Error in input response data:   response type ",resp.type, " specified; there are no uncensored observations"))
        }
    }

return()
}
mysvd<-function(x,  n.components=NULL){
# finds PCs of matrix x
  p<-nrow(x)
  n<-ncol(x)

# center the observations (rows)

 feature.means<-rowMeans(x)
x<- t(scale(t(x),center=feature.means,scale=F))


  if(is.null(n.components)){n.components=min(n,p)}
  if(p>n){
    a<-eigen(t(x)%*%x)
    v<-a$vec[,1:n.components,drop=FALSE]
    d<-sqrt(a$val[1: n.components,drop=FALSE])

      u<-scale(x%*%v,center=FALSE,scale=d)


    return(list(u=u,d=d,v=v))
  }
  else{

      junk<-svd(x,LINPACK=TRUE)
      nc=min(ncol(junk$u), n.components)
      return(list(u=junk$u[,1:nc],d=junk$d[1:nc],
                  v=junk$v[,1:nc]))
}
}
       
permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

