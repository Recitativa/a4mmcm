
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

bpm <- function(T,sigma, mu) {
  nu <- mu/sigma
  B <- rnorm(1,nu*T,sqrt(T))
  M <- B/2+sqrt(B*B-2*T*log(runif(1)))/2
  M <- M*sigma
  return(M)
}

S=100;K=100;r=0.05;sigma=0.2;T=1;alpha=0.5;
price <- function(S=100,K=100,r=0.05,sigma=0.2,T=1,alpha=0.5,n=100000) {
  t1 <- alpha*T; t2 <- (1-alpha)*T;
  mu = r - sigma^2/2
  PP <- array(dim=n)
  QQ <- array(dim=n)
  for(i in 1:n) {
    X1 <- bpm(t1,sigma,mu)
    X2 <- bpm(t2,sigma,-mu)
    Q <- X1 - X2
    QQ[i] <- Q
    PP[i] <- max(0,S*exp(Q)-K) 
  }
  PP <- PP*exp(-r*T)
  C <- mean(PP)
  SD <- sd(PP)/sqrt(n-1)
  return(list(QQ=QQ, PP=PP,C=C, SD=SD))
}


