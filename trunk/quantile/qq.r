phi <- function(t, x, u) {
  iphi <- function(s) {
    return((sqrt(2/(pi*s))*exp(-(u^2)*s/2)
            -2*u*(1-pnorm(u*sqrt(s))))
           *(2*u+sqrt(2/(pi*(t-s)))*exp(-u^2/2*(t-s))
              - 2*u*(1-pnorm(u*sqrt(t-s)))))
  }
  return(0.5*integrate(iphi,0,t*x)$value)
}

G <- function(t, x, alpha, u) {
  z <- rep(0,length(x))
  for(j in 1:length(x)) {
    xx<- x[j]
    iG <- function(s) {
      y <- rep(0,length(s))
      for(i in 1:length(s)) {
        ss <- s[i]
        h <- abs(xx)/sqrt(2*pi*ss) * exp(((abs(xx)-u*ss)^2)/(2*ss))
        y[i] <- h*phi(t-ss,alpha-ss/t,u)
        if(y[i]>1000000) y[i] <- 1000000
      }
      #if(!is.finite(sum(y))) {print(xx); print(s); print(y)}
      return(y)
    }
    z[j] <- integrate(iG,0,t*alpha)$value
  }
  return(z)
}

pppi <- function(T=0.25, c=95, alpha=.8,r=.05, x0=100, sigma=.2) {
  bT <- exp(-r*T)
  ipppi <- function(y) {
    return(G(T,1/sigma*log(y/x0),alpha,r/sigma-sigma/2))
  }
  return(bT*integrate(ipppi,c,200)$value+ c*bT*G(T,1/sigma*log(c/x0),alpha,r/sigma-sigma/2))
}

pppi()
