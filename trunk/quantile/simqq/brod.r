
spm <- function(m, T, sigma, mu) {
  # W = sigma(mu/sigma t + B_t)
  PP <- pnorm((m-mu*T)/sqrt(T)/sigma)-exp(2*mu*m/sigma^2)*pnorm((-m-mu*T)/sqrt(T)/sigma)
  return(PP)
}

rpm <- function(T,sigma,mu) {
  y = runif(1)
  #ff <- function(m) return(abs(spm(m,T,sigma,mu)-y));
  ff <- function(m) return(spm(m,T,sigma,mu)-y);
  xx = uniroot(ff,c(-10,10));
  x = xx$root
  return(x)
}

bpmn <- function(T,sigma, mu,n=1) {
  B <- rnorm(n,mu*T,sigma*sqrt(T))
  M <- (B+sqrt(B^2+rexp(n)*(2*sigma^2*T)))/2
  return(M)
}

bpm <-function(T,sigma, mu) {
  nu <- mu/sigma
  B <- rnorm(1,nu*T,sqrt(T))
  M <- B/2+sqrt(B*B-2*T*log(runif(1)))/2
  M <- M*sigma
  return(M)
}

#S=100;K=100;r=0.05;sigma=0.2;T=1;alpha=0.5;
#source("brod.r");gc();system.time(rr<-price(n=10000000));cat(rr$C,rr$SD,"\n")
 

price <- function(S=100,K=100,r=0.05,sigma=0.2,T=1,alpha=0.5,n=100000) {
  t1 <- alpha*T; t2 <- (1-alpha)*T;
  mu = r - sigma^2/2

  X1 <- bpmn(t1,sigma,mu,n)
  X2 <- bpmn(t2,sigma,-mu,n)
  QQ <- X1 - X2
  PP <- pmax(0,S*exp(QQ)-K)
  C <- mean(PP)*exp(-r*T)
  SD <- sd(PP)/sqrt(n-1)*exp(-r*T)
  return(list(QQ=QQ, PP=PP,C=C, SD=SD))
}


