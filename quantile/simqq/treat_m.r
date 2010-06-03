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
readFile <- function(inFilenames,NN=200000) {
  #Rb =Rm=Re=nRQ=RQuantiles=adat=NULL;
  Nrows <- 0;
  for(i in 1:length(inFilenames)) {
    name <- inFilenames[i]
    fin <- file(name , open="rb")
    Rb <- readBin(fin,integer())
    Rm <- readBin(fin,integer())
    Re <- readBin(fin,integer())
    Nseg <- readBin(fin,integer())
    if(i==1) {
      nRQ <- 2**(Rb-1) + 1
      RQuantiles <- (2**(Rb-1)):(2**Rb)/(2**Rb)
      adat <- array(dim=c(nRQ,Rm-Rb+2,NN),
                    dimnames = list(
                      paste("Q",
                            formatC(RQuantiles*100,
                                    width=2, flag="0"),
                            sep = "_"),
                      c(paste("A", Rb:Rm,sep = "_"),"Dense"),
                      NULL
                      ))
    }
    cat("Nseg",Nseg,"\n")
    cat("Reading data file...",name,"\n");flush(stdout());
    repeat {
      n = nRQ*(Rm-Rb+2) # length of each record
      v <- readBin(fin, double(), n=n)
      if(length(v)!=n) break;
      if(Nrows+1 == NN) break;
      Nrows <- Nrows+1
      adat[,,Nrows] <- v
    }
    close(fin)
    cat("I have read the file...\n");flush(stdout())
  }
  ## adat is a three dimentional array of data with [i,j,k]
  ## i -- i-th record
  ## j -- different k: Rb, Rb+1, Rb+2, ..., Rm, Dense(see below)
  ## k -- Quantiles 2^(Rb-1)/2^Rb, 2^(Rb-1)+1/2^Rb, ...., 1
  adat <- aperm(adat[,,1:Nrows],perm=c(3,2,1))
  ret <- list(adat = adat, Rb=Rb, Rm=Rm, Re = Re, Nseg = Nseg, 
              nRQ=nRQ, ndata=n, Nrows=Nrows, RQuantiles =RQuantiles)
  return(ret)
}

## ns <- c(paste("nout_4_16_25_1_0_",10:24,".bin",sep=""))


## main function to treat the data
run <- function(inFilenames,ttt=8,SP=1,AB=TRUE,DQ=FALSE,ZEE=FALSE,NN=200000) {
  ## read data from file
  inFilenames <- c(inFilenames)
  ret <- readFile(c(inFilenames),NN=NN)
  attach(ret)
  ## file name pattern: ?????.bin, then get ????? part. 
  boutname <- sub(".bin$","", inFilenames[1])

  if(ZEE)
    adat[,"Dense","Q_50"] <- 0;
  
  ## Formula: (A)dErr = X_k -  X_Dense
  dErr <- adat[,1:(Rm-Rb+1),] - adat[,rep("Dense",Rm-Rb+1),]
  AdErr <- abs(dErr);
  ## Formula: m(A)dErr = |mean(X_k -  X_Dense)|
  ## mAdErr is a 2-dim array, 1st-dimension is k, 2nd-dim is Quantile
  mdErr <- apply(dErr, c(2,3), mean)[1:(Rm-Rb+1-ttt),]
  mdErr <- abs(mdErr);
  mAdErr <- apply(AdErr, c(2,3), mean)[1:(Rm-Rb+1-ttt),]

  PRQ <- c(seq(1,nRQ,by=SP),nRQ); PRQ <- unique(PRQ);

  plotlines <- function(fname,Dat,xlab="s",ylab="Error") {
    pdf(fname,pointsize=8)
    plot(c(Rb,Rm-ttt),
         c(min(Dat),max(Dat)),
         type="n", xlab = xlab, ylab=ylab)
    cols <- rainbow(nRQ)
    legend(x="topright", paste("",RQuantiles[PRQ]),
           col=cols[PRQ], lty=1, ncol=3)
    for(i in PRQ) {
      lines(Rb:(Rm-ttt), Dat[,i], type="b",col=cols[i])
    }
    dev.off()
  }
  
  ## plot the lines of mAdErr for different Quantile
  plotlines(paste(boutname, ".pdf",sep=""),mdErr,ylab="Error");
  plotlines(paste(boutname, "_abs.pdf",sep=""), mAdErr, ylab="|Error|");
  l2mdErr <- log(mdErr)/log(2);
  l2mAdErr <- log(mAdErr)/log(2);
  plotlines(paste(boutname, "_log.pdf",sep=""), l2mdErr, ylab="log(|Error|)/log(2)");
  plotlines(paste(boutname, "_abs_log.pdf",sep=""), l2mAdErr, ylab="log(|Error|)/log(2)");


  ## regression for each quantile & plot the result.
  plotrate <- function(fname, Dat,
                       xlab="Quantile",
                       ylab="Empirical convergence rate") {
    alm <- function(x) {return(lm(Q~P, data.frame(Q=x[1:(Rm-Rb-ttt+1)],P=Rb:(Rm-ttt))))}
    rLM <- apply(Dat, c(2), alm)
    getp <- function(x) {return(coef(x)["P"])}
    PrLM <- sapply(rLM, getp)
    PrLM <- -as.vector(PrLM)
    PSrLM <- sapply(PrLM,formatC,digits=3,format="f")
    pdf(fname)
    plot(RQuantiles[PRQ], PrLM[PRQ], xlab=xlab, ylab=ylab)
    text(RQuantiles[PRQ], PrLM[PRQ], labels=PSrLM[PRQ],pos=1,col="red",)
    dev.off()
  }
  plotrate(paste(boutname,"rato.pdf",sep=""), l2mdErr);
  plotrate(paste(boutname,"abs_rato.pdf",sep=""), l2mAdErr);

  if(DQ) {
    ## compare the empirical distribution and theoretical distribution.
    ## each alpha a pdf file, and for each k a line.
    ## The black line is theoretic and blue line is the "dense" one.
    PPP <- Rm-Rb+2
    cols <- heat.colors(2*PPP)[1:PPP]
    cols[PPP] = "blue"
                                        # draw the distribution of quantdiles
    tt <- c(seq(-6,6,length.out=300))
    for(i in 1:nRQ) {
      alpha <- RQuantiles[i]
      print(alpha)
      FX <- FppX(tt, alpha)
      pdf(paste(boutname,'_',alpha, '.quantile.pdf',sep=""))
      plot(tt,FX,type='l',col='black',
           main=sprintf('compare with Bb=%g Be=%g,alpha=%g, sim=%g',
             Rb,Rm,alpha,Nrows),
           ylab="Empirical Cumulative Distribution Function",xlab='X')
      for(j in 1:PPP) {
        EeulerSim <- ecdf(adat[,j,i])
        lines(tt,EeulerSim(tt),type='l',col=cols[j])
      }
      dev.off()
    }
    dev.off()
  }
  detach(ret)
  return(ret$Nrows)
}

gg <- function(x) {
  return(x*pnorm(x)+1/sqrt(2*pi)*exp(-x^2/2))
}

CC <- function(mu=0.4,sigma=1,T=1,alpha=0.5) {
  mu1 <- mu*sqrt(alpha*T)/sigma;
  mu2 <- -mu*sqrt((1-alpha)*T)/sigma;
  ret <- -sigma*sqrt(T)*( (2*gg(mu1)-mu1)*sqrt(alpha)
                         -(2*gg(mu2)-mu2)*sqrt(1-alpha))
  return(ret)
}

CCC <- function(mu=0, sigma=1,T=1,alpha=0.5,N=2^4) {
  return(log(abs(CC(mu,sigma,T,alpha)/4/N))/log(2))
}
