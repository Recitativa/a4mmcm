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
      #nRQ <- 2**(Rb-1) + 1
      RQuantiles <- (2**(Rb-1)):(2**Rb)/(2**Rb)
      nRQ <- length(RQuantiles)
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
      if(length(v)!=n) {
        if(length(v)!=0)
          cat("sprising endding!\n");
        break;}
      Nrows <- Nrows+1
      adat[,,Nrows] <- v
      if(Nrows == NN) break;
    }
    close(fin)
    cat("I have read the file...\n");flush(stdout())
    if(Nrows==NN) break; # break if the number of records reach NN 
  }
  ## adat is a three dimentional array of data with [i,j,k]
  ## i -- i-th record
  ## j -- different k: Rb, Rb+1, Rb+2, ..., Rm, Dense(see below)
  ## k -- Quantiles 2^(Rb-1)/2^Rb, 2^(Rb-1)+1/2^Rb, ...., 1
 
## adat <- aperm(adat[,,1:Nrows],perm=c(3,2,1))
  adat <- adat[,,1:Nrows]
  ret <- list(adat = adat, Rb=Rb, Rm=Rm, Re = Re, Nseg = Nseg, 
              nRQ=nRQ, ndata=n, Nrows=Nrows, RQuantiles =RQuantiles)
  return(ret)
}

## ns <- c(paste("nout_4_16_25_1_3_",100:114,".bin",sep=""))


## main function to treat the data
run <- function(inFilenames,
                bbb=0,ttt=0,SP=1,
                AB=TRUE,DQ=FALSE, ZEE=FALSE,
                NN=200000,PRQ=NA) {
  ## read data from file
  inFilenames <- c(inFilenames)
  ret <- readFile(inFilenames,NN=NN)
  attach(ret)
  ## file name pattern: ?????.bin, then get ????? part. 
  boutname <- sub(".bin$","", inFilenames[1])

  if(ZEE) {
    adat[1,"Dense",] <- 0;
  }
  
  ## Formula: (A)dErr =   X_Dense - X_n
  ## Formula: m(A)dErr = |mean(X_k -  X_Dense)|
  ## mAdErr is a 2-dim array, 1st-dim is Quantile, 2ed-dimension is k, 
  dErr <- adat[,rep("Dense",Rm-Rb+1),] -  adat[,1:(Rm-Rb+1),]
  mdErr <- apply(dErr, c(1,2), mean)[,1:(Rm-Rb+1-ttt)]
  dErr <- abs(dErr);
  mAdErr <- apply(dErr, c(1,2), mean)[,1:(Rm-Rb+1-ttt)]
  #sdErr <- apply(dErr, c(1,2),sd)[1:(Rm-Rb+1-ttt),]/sqrt(Nrows)
  #sdAErr <- apply(AdErr, c(2,3),sd)[1:(Rm-Rb+1-ttt),]/sqrt(Nrows)

  if(is.na(PRQ))
     PRQ <- c(seq(1,nRQ,by=SP),nRQ); PRQ <- unique(PRQ);

  plotlines <- function(fname,Dat,sdDat=NULL,xlab="s",ylab="Error") {
    pdf(fname,pointsize=8)
    plot(c(Rb,Rm-ttt),
         range(Dat),
         type="n", xlab = xlab, ylab=ylab)
    cols <- rainbow(nRQ)
    legend(x="topright", paste("",RQuantiles[PRQ]),
           col=cols[PRQ], lty=1, ncol=3)
    for(i in PRQ) {
      lines(Rb:(Rm-ttt), Dat[i,], type="b",col=cols[i])
    }
    dev.off()
  }
  
  ## plot the lines of mAdErr for different Quantile
  plotlines(paste(boutname, ".pdf",sep=""),mdErr,ylab="Error");
  plotlines(paste(boutname, "a.pdf",sep=""), mAdErr, ylab="|Error|");

  mdErr <- abs(mdErr);
  l2mdErr <- log(mdErr)/log(2);
  l2mAdErr <- log(mAdErr)/log(2);
  plotlines(paste(boutname, "log.pdf",sep=""), l2mdErr, ylab="log(Error)/log(2)");
  plotlines(paste(boutname, "alog.pdf",sep=""), l2mAdErr, ylab="log(|Error|)/log(2)");

  ## regression for each quantile & plot the result.
  plotrate <- function(fname, Dat,
                       xlab=expression(alpha),
                       ylab=expression(lambda)) {
    alm <- function(x) {return(lm(Q~P, data.frame(Q=x[(bbb+1):(Rm-Rb-ttt+1)],P=(Rb+bbb):(Rm-ttt))))}
    rLM <- apply(Dat, c(1), alm)
    getp <- function(x) {return(coef(x)["P"])}
    PrLM <- sapply(rLM, getp)
    PrLM <- -as.vector(PrLM)
    PSrLM <- sapply(PrLM,formatC,digits=3,format="f")
    pdf(fname)
    plot(RQuantiles[PRQ], PrLM[PRQ], xlab=xlab, ylab=ylab)
    text(RQuantiles[PRQ], PrLM[PRQ], labels=PSrLM[PRQ],pos=1,col="red",)
    dev.off()

    geti <- function(x) {return(coef(x)[1])}
    PrI <- sapply(rLM,geti)
    PrI <- as.vector(PrI)
    return(PrI)
  }
  PrI  <- plotrate(paste(boutname,"rato.pdf",sep=""), l2mdErr);
  PrAI <- plotrate(paste(boutname,"arato.pdf",sep=""), l2mAdErr);

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
      cat("treating ", alpha,"\n")
      FX <- FppX(tt, alpha)
      pdf(paste(boutname,'_',alpha, '.quantile.pdf',sep=""))
      plot(tt,FX,type='l',col='black',
           main=sprintf('compare with Bb=%g Be=%g,alpha=%g, sim=%g',
             Rb,Rm,alpha,Nrows),
           ylab="Empirical Cumulative Distribution Function",xlab='X')
      for(j in 1:PPP) {
        EeulerSim <- ecdf(adat[i,j,])
        lines(tt,EeulerSim(tt),type='l',col=cols[j])
      }
      dev.off()
    }
  }
  TrI <- CCC(alpha=RQuantiles[PRQ]);
  ftxt <- file(paste(boutname,".txt",sep=""),"wt")
  write(Nrows,file=ftxt)
  write.table(l2mdErr,file=ftxt)
  write.table(PrI,file=ftxt)
  write.table(TrI,file=ftxt)
  close(ftxt)
  detach(ret)
  return(list(Nrows=ret$Nrows, l2mdErr=l2mdErr,PrI=PrI, TrI=TrI))
}

gg <- function(x) {
  return(x*pnorm(x)+exp(-x^2/2)/sqrt(2*pi))
}

CC <- function(mu=0.4,sigma=1,T=1,alpha=0.5) {
  mu1 <- mu*sqrt(alpha*T)/sigma;
  mu2 <- -mu*sqrt((1-alpha)*T)/sigma;
  ret <- -sigma*sqrt(T)/4 *( (2*gg(mu1)-mu1)/sqrt(alpha)
                         -(2*gg(mu2)-mu2)/sqrt(1-alpha))
  return(ret)
}

CCC <- function(mu=0, sigma=1,T=1,alpha=0.5,N=1) {
  return(log(CC(mu,sigma,T,alpha)/N)/log(2))
}

plotI <- function(mu=1:10,n=500,ylim=c(-1.5,1.5)) {
  plot(0,0, xlim=c(0,1), ylim=ylim,
       xlab=expression(alpha),
       ylab="coefficient",
       type='n');
  cols <- rainbow(length(mu))
  alpha <- seq(0,1,length.out=n)
  for(i in 1:length(mu)) 
    lines(alpha,CC(mu=mu[i],alpha=alpha),type='l',col=cols[i])
  legend(x="topleft",
         as.expression(
             lapply(mu, function(x) bquote(mu == .(x))),
             col=cols,lty=1)
         )
}
