## Theoretical distribution of quantiles
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

## distribution of quantiles
## i.e. FppX(x,alpha) = P(X_alpha < x)
FppX <- function(end=Inf,alpha=.6) {
	y <- 0
	pX <- function(x) {ppX(x,alpha)}
	for (i in 1:length(end)) {
		y[i] <- integrate(pX,-Inf,end[i])$value
		}
	return(y)
	}



## read data from binary file
readFile <- function(inFilename) {
  fin <- file(inFilename, open="rb")
  Rb <- readBin(fin,integer())
  Re <- readBin(fin,integer())
  P2 <- readBin(fin,integer())
  Nseg <- readBin(fin,integer())
  nRQ <- 2**(Rb-1) + 1
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
  ## adat is a three dimentional array of data with [i,j,k]
  ## i -- i-th record
  ## j -- different k: Rb, Rb+1, Rb+2, ..., Re, Dense(see below)
  ## k -- Quantiles 2^(Rb-1)/2^Rb, 2^(Rb-1)+1/2^Rb, ...., 1
  adat <- aperm(adat,perm=c(3,2,1))
  ret <- list(adat = adat, Rb=Rb, Re=Re, P2 = P2, Nseg = Nseg, 
              nRQ=nRQ, ndata=n, RQuantiles =RQuantiles)
  close(fin)
  return(ret)
}


## main function to treat the data
run <- function(inFilename) {
  ## file name pattern: ?????.bin, then get ????? part. 
  boutname <- sub(".bin$","", inFilename)
  ## read data from file
  ret <- readFile(inFilename)
  attach(ret)
  ttt <- 4

  ## Formula: AdErr = X_k -  X_Dense
  dErr <- adat[,1:(Re-Rb+1),] - adat[,rep("Dense",Re-Rb+1),]
  AdErr <- dErr
  ## Formula: mAdErr = |mean(X_k -  X_Dense)|
  ## mAdErr is a 2-dim array, 1st-dimension is k, 2nd-dim is Quantile
  mAdErr <- abs(apply(AdErr, c(2,3), mean,trim=.05))[1:(Re-Rb+1-ttt)]

  ## plot the lines of mAdErr for different Quantile
  pdf(paste(boutname,".lines.pdf",sep=""),pointsize=8)
  plot(c(Rb,Re-ttt),
       c(min(mAdErr),max(mAdErr)),
       type="n", xlab = "k", ylab="|Err|",
       main="|Err| v.s. k, under different Quantiles")
  cols <- rainbow(nRQ)
  legend(x="topright", paste("",RQuantiles),
         col=cols, lty=1, ncol=3)
  for(i in 1:nRQ) {
    lines(Rb:(Re-ttt), mAdErr[1:(Re-Rb-ttt+1),i], type="b",col=cols[i])
  }
  dev.off()

  ## plot log(|Err|)/log(2)
  pdf(paste(boutname,".l.lines.pdf",sep=""),pointsize=8)
  l2mAdErr <- log(mAdErr)/log(2)
  plot(c(Rb,Re-ttt),
       c(min(l2mAdErr),max(l2mAdErr)),
       type="n", xlab = "k", ylab="log(|Err|)/log(2)",
       main="log(|Err|) v.s. k, under different Quantiles")
  cols <- rainbow(2*nRQ)[(nRQ+1):(2*nRQ)]
  legend(x = "topright", paste("",RQuantiles),
         col=cols, lty=1,ncol=3)
 
  for(i in 1:nRQ) {
    lines(Rb:(Re-ttt), l2mAdErr[1:(Re-Rb-ttt+1),i],  type="b",col=cols[i])
  }
  dev.off()
  
  ## regression for each quantile & plot the result.
  alm <- function(x) {return(lm(log(Q)~P, data.frame(Q=x[1:(Re-Rb-ttt+1)],P=Rb:(Re-ttt))))}
  rLM <- apply(mAdErr, c(2), alm)
  getp <- function(x) {return(coef(x)["P"])}
  PrLM <- sapply(rLM, getp)
  pdf(paste(boutname,".rato.pdf",sep=""))
  plot(RQuantiles, -as.vector(PrLM), ylab="Quantile", xlab="Empirical rato",
       main=sprintf(
         "Empirical discrete Error Rato under different Quantiles, sim=%g",
         ndata))
  dev.off()

  if(FALSE) {
    ## compare the empirical distribution and theoretical distribution.
    ## each alpha a pdf file, and for each k a line.
    ## The black line is theoretic and blue line is the "dense" one.
    PPP <- Re-Rb+2
    cols <- heat.colors(2*PPP)[1:PPP]
    cols[PPP] = "blue"
                                        # draw the distribution of quantiles
    tt <- c(seq(-6,6,length.out=300))
    for(i in 1:nRQ) {
      alpha <- RQuantiles[i]
      print(alpha)
      FX <- FppX(tt, alpha)
      pdf(paste(boutname,'_',alpha, '.quantile.pdf',sep=""))
      plot(tt,FX,type='l',col='black',
           main=sprintf('compare with Bb=%g Be=%g,alpha=%g, sim=%g',
             Rb,Re,alpha,ndata),
           ylab="Empirical Cumulative Distribution Function",xlab='X')
      for(j in 1:PPP) {
        EeulerSim <- ecdf(adat[,j,i])
        lines(tt,EeulerSim(tt),type='l',col=cols[j])
      }
      dev.off()
    }
    detach(ret)
    dev.off()
  }
  return(0)
}
