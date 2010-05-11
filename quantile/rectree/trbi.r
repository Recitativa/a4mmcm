readtxt <- function(fn,name="S90") {
  tt <- read.table(fn,skip=1);
  dd[name] <<- NA;
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

ptab <- function(data=df, ind=3:8) {
  tt <- ""
  for(i in ind) {
    cat(paste(data[i,]),sep=" & ")
    cat("\\\\\n") 
  }
  return(tt)
}

  
dd = data.frame(STEP=1:40)
