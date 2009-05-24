readFile <- function(inFilename) {
  fin <- file(inFilename, open="rb")
  P2 <- readBin(fin, integer())
  Rb <- readBin(fin,integer())
  Re <- readBin(fin,integer())
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
                  c("Dense", paste("A", Rb:Re,sep = "_")),
                  NULL
                  )
                )
  adat <- aperm(adat,perm=c(3,2,1))
  ret <- list(adat = adat, P2=P2, Rb=Rb, Re=Re,
              nRQ=nRQ, RQuantiles =RQuantiles)
  close(fin)
  return(ret)
}



run <- function(inFilename = "sout_22_ 4_17_22.bin") {
  boutname <- sub(".bin$","", inFilename)
  ret <- readFile(inFilename)
  Re = ret$Re
  Rb = ret$Rb
  nRQ = ret$nRQ
  dErr <- ret$adat[,2:(Re-Rb+2),] - ret$adat[,rep("Dense",Re-Rb+1),]
  AdErr <- abs(dErr)
  mAdErr <- apply(AdErr, c(2,3), mean)
  fmAdErr <- as.data.frame(mAdErr)
  fmAdErr$P <- Rb:Re

  pdf(paste(boutname,"lines.pdf"))
  plot(c(Rb,Re),
       c(min(mAdErr),max(mAdErr)),
       type="n")
  cols <- rainbow(nRQ)
  for(i in 1:nRQ) {
    lines(Rb:Re, mAdErr[,i], type="b",col=cols[i])
  }
  dev.off()
  
  pdf(paste(boutname,"l.lines.pdf"))
  l2mAdErr <- log(mAdErr)/log(2)
  plot(c(Rb,Re),
       c(min(l2mAdErr),max(l2mAdErr)),
       type="n")
  cols <- rainbow(nRQ)
  for(i in 1:nRQ) {
    lines(Rb:Re, l2mAdErr[,i], type="b",col=cols[i])
  }

  
  dev.off()
  alm <- function(x) {return(lm(log(Q)~P, data.frame(Q=x,P=Rb:Re)))}
  rLM <- apply(mAdErr, c(2), alm)
  getp <- function(x) {return(coef(x)["P"])}
  PrLM <- sapply(rLM, getp)
  pdf(paste(boutname,"rato.pdf"))
  plot(PrLM)
  dev.off()
}
