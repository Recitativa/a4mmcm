ppX <- function(x,alpha=.6,beta=sqrt((1-alpha)/alpha)) {
	c1 = sqrt(2/pi)*2
	y = 0
	for (i in 1:length(x)) {
		if(x[i]>=0)
			y[i] <- c1*exp(-x[i]^2/2)*(1-pnorm(x[i]*beta)) 
		else
			y[i] <- c1*exp(-x[i]^2/2)*(1-pnorm(-x[i]/beta))
		}
	return(y)
	}
	
FppX <- function(end=Inf,alpha=.6) {
	y <- 0
	pX <- function(x) {ppX(x,alpha)}
	for (i in 1:length(end)) {
		y[i] <- integrate(pX,-Inf,end[i])$value
		}
	return(y)
	}




readFile <- function(inFilename) {
  fin <- file(inFilename, open="rb")
  Rb <- readBin(fin,integer())
  Re <- readBin(fin,integer())-1
  nRQ <- 2**(Rb-1) +1
  RQuantiles <- (2**(Rb-1)):(2**Rb)/(2**Rb)

    
  rawdat <- vector() 
  Nrows <- 0; 
  repeat {
    n = nRQ*(Re-Rb+1+1) # length of each record
    v <- readBin(fin, double(), n=n)
    if(length(v)!=n) break;
    Nrows <- Nrows+1
    rawdat <- append(rawdat, v)
  }
  adat <- array(rawdat, dim=c(nRQ,Re-Rb+1+1,Nrows),
                dimnames = list(
                  paste("Q",
                        formatC(RQuantiles*100,
                                width=2, flag="0"),
                        sep = "_"),
                  c(paste("A", Rb:Re,sep = "_"),"Dense"),
                  NULL
                  )
                )
  adat <- aperm(adat,perm=c(3,2,1))
  ret <- list(adat = adat, Rb=Rb, Re=Re,
              nRQ=nRQ, RQuantiles =RQuantiles)
  close(fin)
  return(ret)
}



run <- function(inFilename = "nout_4_20.bin") {
  boutname <- sub(".bin$","", inFilename)
  ret <- readFile(inFilename)
  Re = ret$Re
  Rb = ret$Rb
  nRQ = ret$nRQ
  RQuantiles <- ret$RQuantiles
  adat <- ret$adat


  dErr <- adat[,1:(Re-Rb+1),] - adat[,rep("Dense",Re-Rb+1),]
  AdErr <- abs(dErr)
  mAdErr <- apply(AdErr, c(2,3), mean,trim=.05)
  fmAdErr <- as.data.frame(mAdErr)
  fmAdErr$P <- Rb:Re

  pdf(paste(boutname,".lines.pdf",sep=""))
  plot(c(Rb,Re),
       c(min(mAdErr),max(mAdErr)),
       type="n")
  cols <- rainbow(nRQ)
  for(i in 1:nRQ) {
    lines(Rb:Re, mAdErr[,i], type="b",col=cols[i])
  }
  dev.off()
  
  pdf(paste(boutname,".l.lines.pdf",sep=""))
  l2mAdErr <- log(mAdErr)/log(2)
  plot(c(Rb,Re),
       c(min(l2mAdErr),max(l2mAdErr)),
       type="n", xlab = "k", ylab="log(|Err|)/log(2)",
       main="log(|Err|) v.s. k, under different Quantiles")
  cols <- rainbow(nRQ)
  for(i in 1:nRQ) {
    lines(Rb:Re, l2mAdErr[,i], type="b",col=cols[i])
  }

  
  dev.off()
  alm <- function(x) {return(lm(log(Q)~P, data.frame(Q=x,P=Rb:Re)))}
  rLM <- apply(mAdErr, c(2), alm)
  getp <- function(x) {return(coef(x)["P"])}
  PrLM <- sapply(rLM, getp)
  pdf(paste(boutname,".rato.pdf",sep=""))
  plot(-as.vector(PrLM),RQuantiles, xlab="Quantile", ylab="Empirical rato",
       main="Empirical discrete Error Rato under different Quantiles")
  dev.off()


  PPP <- Re-Rb+2
  cols <- heat.colors(2*PPP)[1:PPP]
  # draw the distribution of quantiles
  tt <- c(seq(-6,6,length.out=300))
  for(i in 1:nRQ) {
    alpha <- RQuantiles[i]
    print(alpha)
    FX <- FppX(tt, alpha)
    pdf(paste(boutname,'_',alpha, '.quantile.pdf',sep=""))
    plot(tt,FX,type='l',col='black',
         main=sprintf('compare with Bb=%g Be=%g,alpha=%g',Rb,Re+1,alpha),
         ylab="Empirical Cumulative Distribution Function",xlab='')
    for(j in 1:PPP) {
      EeulerSim <- ecdf(adat[,j,i])
      lines(tt,EeulerSim(tt),type='l',col=cols[j])
    }
    dev.off()
    }
  dev.off()
  return(0)
}
