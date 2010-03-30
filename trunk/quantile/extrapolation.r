richard <- function(d,p=.5) {
  n <- length(d$h)
  aa <- array(rep(0.0,n*n),dim=c(n,n))
  aa[,1] <- d$v;
  tt <- d$h[2:n]/d$h[1:n-1]
  if( length(p)==1) p<- 1:n*p
  for(i in 1:(n-1))
    aa[1:(n-i),i+1] <- (tt[1:(n-i)]^p[i]*aa[1:(n-i),i]-
                        aa[2:(n-i+1)]
                        )/(tt[1:(n-i)]^p[i]-1)
  d$aa <- aa
  return(d)
}

aitken <- function(v) {
  n <- length(v)
  return(v[1:(n-2)]-diff(v)[-n+1]^2/diff(v,differences=2))
}

ind <- seq(30,38,by=2)
h <- 1/ind
v <- c(1.49669, 1.5078, 1.51042,1.51703,1.51923)
v <- c(1.80059, 1.78771, 1.77631,1.76377,1.75524)
d <- data.frame(ind,h,v)

richard(d)
