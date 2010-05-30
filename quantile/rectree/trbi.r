readtxt <- function(fn,name="S90") {
  tt <- read.table(fn,skip=1);
  dd[name] <<- NA; # ,<<- operator global variable 
  for(i in 1:length(tt$V1)) dd[tt$V1[i],name] <<- tt$V2[i];
}

extr <- function(name="S90",sep=2) {
  n <- length(dd$STEP)
  rname <- paste("R",name,sep="")
  dd[rname] <<- NA;
  for(i in 1:sep) {
    for(j in 0:(n/sep-1)) {
      id1 <- j*sep+i;
      id2 <- (j+1)*sep+i;
      if(id2<=n) {
        t <- dd$STEP[id2]/dd$STEP[id1]
        dd[id2,rname] <<- (t*dd[id2,name]-dd[id1,name])/(t-1);
      }
    }
  }
}

## ops <- c("S90","S95","S100","S105");
## ddb <- dd[10:38,c("STEP",ops)]

extra <- function(ddb, ops,tor=0.001) {
  lam <- seq(.2,2.1,length.out=20); #set{0.1,.2,.3...1}
  n <- length(ddb$STEP)
  seps <- seq(2,20,by=1);# seq gives a sequence
  dn <- list(ddb$STEP, ops, seps, lam); # give a jogged structure
  res <- array(dim=lapply(dn,length),dimnames=dn) #array gives a multidimitional array
  err <- array(dim=lapply(dn,length),dimnames=dn)
  for(name in ops) { 
    for(sep in seps) {
      for(i in 1:n) {
        id1 <- i;
        id2 <- i+sep;
        if(id2<=n) {
          t <- ddb$STEP[id2]/ddb$STEP[id1]
          res[paste(ddb$STEP[id2]), #paste change objects into string
              name,paste(sep),] <- 
                ((t^lam)*ddb[id2,name]-ddb[id1,name])/(t^lam-1)
        }
      }
    }
  }
  bbb <- array(c(1.62866,3.21382,5.65519,8.96151),dimnames=list(ops)) # dimnames must be a list 
  for(name in ops)
    err[,name,,] <- res[,name,,] - bbb[name] # if indicator not specified, operator acordingly
  err <- err^2; 
  merr <- apply(err,c(1,3,4),sum) #apply a function to matrix
  mm <- min(merr,na.rm=TRUE)
  rp <- list();
  for(sep in paste(seps))
    for(s in paste(ddb$STEP))
      for(ll in paste(lam)) {
        v <- merr[s,sep,ll]
        if (!is.na(v) && v< tor && as.integer(s)>35)
          rp <- append(rp,paste(sep, as.integer(s)-as.integer(sep),s,ll,paste(err[s,,sep,ll]),merr[s,sep,ll]))
      }
  return(list(res=res,err=err,merr=merr,rp=rp))
}

## treat res
## RESULT: sep=5, step=36, lam=0.5  

ptab <- function(data=df, ind=3:8) {
  tt <- ""
  for(i in ind) {
    cat(paste(data[i,]),sep=" & ")
    cat("\\\\\n") 
  }
  return(tt)
}

  
#dd = data.frame(STEP=1:40)
