require(graphics)

Pr <- function(t, x, sigma, mu) {
  A <- exp(2*mu*x/(sigma^2))
  B <- 1-pnorm((x+mu*t)/(sigma*sqrt(t)))
  C <- pnorm((x-mu*t)/(sigma*sqrt(t)))
  val <- -A*B + C
  return(val)
}



## main function to treat the data
run <- function(T = .25, alpha=.8, sigma=.2, r=.05, S=100, K=95, n=10000) {
  x <- c(seq(0,5,length.out=1000),seq(2,100,length.out=100))
  yy <- seq(0,1, length.out=1000);
  mu <- r-sigma^2/2
  PM1 <- Pr(alpha*T, x, sigma, mu)
  PM2 <- Pr((1-alpha)*T, x, sigma,-mu)
  
  invPM1 <- approxfun(PM1,x)
  invPM2 <- approxfun(PM2,x)
  plot(PM1,x,type='l',col='red');
  lines(PM2,x,type='l',col='blue');
  inM1 <- invPM1(yy)
  inM2 <- invPM2(yy)
  lines(yy,inM1,col='green');
  lines(yy,inM2,col='green');
  
  simM1 <- invPM1(runif(n))
  simM2 <- invPM2(runif(n))
  Sat <- S*exp(simM1-simM2)
  V <- rep(0,n);
  for(j in 1:n) {
    V[j] = max(Sat[j]-K,0)
  }
  V <- exp(-r*T)*V
  
  M <- mean(V)
  print(t.test(V))
  return(V)
}

