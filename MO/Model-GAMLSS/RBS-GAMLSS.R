
RBS <- function (mu.link = "log" , sigma.link="log")
{
  mstats = checklink("mu.link", "BirnbaumSaunders", substitute(mu.link),
                     c("sqrt","log","identity"))
  dstats = checklink("sigma.link", "BirnbaumSaunders", substitute(sigma.link)
                     ,c("sqrt", "log", "identity"))
  structure(list(family = c("RBS","BirnbaumSaunders"),
                 parameters = list(mu=TRUE,sigma=TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 #the first derivative of the likelihood with respect to the location parameter mu
                 dldm = function(y,mu,sigma) #first derivate of log-density respect to mu
                 {
                   mustart = 1/(2*mu)
                   ystart =  ((sigma+1)*y)/(4*mu*mu) - (sigma^2)/(4*(sigma+1)*y) + sigma/((sigma*y) + y + (sigma*mu))
                   
                   dldm = ystart-mustart
                   dldm
                 },
                 #the expected second derivative of the likelihood with respect to the location parameter mu
                 d2ldm2 = function(mu,sigma) {        #expected of second derivate of log-density respect to mu
                   d2ldm2 =  - sigma/(2*mu*mu) - ((sigma/(sigma+1))^2)*Ims(mu,sigma)
                   d2ldm2
                 },
                 
                 #the first derivative of the likelihood with respect to the scale parameter sigma
                 dldd = function(y,mu,sigma) {      #first derivate log-density respect to sigma
                   sigmastart  = -(sigma)/(2*(sigma+1))
                   y2start   = (y+mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (mu*sigma*(sigma+2))/(4*y*((sigma+1)^2))
                   dldd  = y2start-sigmastart
                   dldd
                 },
                 #the expected second derivative of the likelihood with respect to the scale parameter sigma ok
                 d2ldd2 = function(mu,sigma) {      #expected of second derivate log-density respect to sigma
                   lss =  ((sigma^2) + (3*sigma) + 1)/(2*sigma*sigma*((sigma+1)^2))
                   d2ldd2 = -lss - ((mu^2)/((sigma+1)^4))*Ims(mu,sigma)
                   d2ldd2
                 },
                 #the expected cross derivative of the likelihood with respect to both the location mu and scale parameter sigma
                 d2ldmdd = function(mu,sigma) {   #expected of partial derivate of log-density respect to mu and sigma
                   lms = 1/(2*mu*(sigma+1))
                   d2ldmdd = - lms - ((mu*sigma)/((sigma+1)^3))*Ims(mu,sigma)
                   d2ldmdd
                 },
                 
                 G.dev.incr = function(y,mu,sigma,...) -2*dRBS(y,mu,sigma,log=TRUE),
                 rqres = expression(rqres(pfun = "pRBS", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- y + mean(y)/2 }),
                 sigma.initial = expression({sigma <- rep(1,length(y)) }),
                 mu.valid = function(mu) all(mu>0) ,
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0),
                 mean = function(mu,sigma) mu,
                 variance = function(mu, sigma) (mu*mu)*((2*sigma+5)/((sigma+1)^2)) ),
            class = c("gamlss.family","family"))
}




sigmatil=function(y)
{
  s = mean(y)
  r = 1/mean(1/y)
  alphatil = (2*( (s/r)^(1/2)  - 1))^(1/2)
  dest = 2/(alphatil^2 )
  return(dest)
}


Ims <- function(mu,sigma)
{
  esp = function(mu=1,sigma=1)
  {
    
    integral=function(aest)
    {
      fu=function(u)
      {
        w1 = (1 / ((1 +u^2)*(u^2)))
        w2 = (exp((-1 /(2*aest^2) )*((u - 1/u)^2)))
        (w1*w2)
      }
      return(integrate(fu,0,Inf)$value)
    }
    
    const = function(alpha,beta)
    {
      const = 1/(alpha*beta*beta*sqrt(2*pi))
      return(const)
    }
    
    alpha = sqrt(2/sigma)
    beta = (mu*sigma)/(sigma+1)
    e = const(alpha,beta)*integral(alpha)
    return(e)
  }
  
  res <- mapply(esp, mu,sigma)
  
  res
}




resrbs=function(y,mu,sigma)
{
  
  z = -1/(2*mu) - (sigma^2)/(4*(sigma+1)*y) + ((sigma+1)*y)/(4*mu*mu) + sigma/((sigma*y) + y + (sigma*mu))
  v = sigma/(2*mu*mu) + ((sigma*sigma)/((sigma+1)*(sigma+1)))*Ims(mu,sigma)
  res = z/sqrt(v)
  return(res)
}



dRBS<-function(x, mu=1, sigma=1, log=FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))  stop(paste("x must be positive", "\n", ""))
  log.lik =  0.5*(sigma - log(mu) + log(sigma+1) - log(16*pi)) - (3/2)*log(x) - ((sigma+1)/(4*mu))*x - ((mu*sigma*sigma)/(4*(sigma+1)))*(1/x)  + log(x + ((mu*sigma)/(sigma+1)))
  if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
  fy
}


pRBS <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{       if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  cdf1 <- pnorm((1/a)*(sqrt(q/b) - sqrt(b/q)))
  cdf <- cdf1
  
  ## the problem with this approximation is that it is not working with
  ## small sigmas and produce NA's. So here it is a solution
  if (any(is.na(cdf)))
  {
    index <- seq(along=q)[is.na(cdf)]
    for (i in index)
    {
      cdf[i] <- integrate(function(x)
        dcbs(x, alpha = a[i], beta = b[i], log=FALSE), 0.001, q[i] )$value
    }
  }
  
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
  cdf
}


qRBS = function (p, mu = 0.5, sigma = 1, lower.tail = TRUE,
                 log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive ", "\n", ""))
  if (any(sigma < 0))  #In this parametrization  sigma = phi
    stop(paste("sigma must be positive", "\n", ""))
  
  
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p <= 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  
  suppressWarnings(q <- qcbs(p ,alpha = a, beta = b, lower.tail = TRUE, log.p = FALSE))
  q
}



rRBS = function(n, mu=1, sigma=1)
{
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  a = sqrt(2/sigma)
  b = (mu*sigma)/(sigma+1)
  r = rcbs(n, alpha = a, beta = b)
  r
}






dcbs <- function(x, alpha = 1, beta = 1, log = FALSE){
  if(!is.numeric(x)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- x
  c   <-(1 / sqrt(2 * pi))
  u   <- (alpha^(-2)) * ((x / beta) + (beta / x) - 2)
  e   <- exp((-1 / 2 ) * u)
  du  <- (x^(-3 / 2) * (x + beta)) /
    (2 * alpha * sqrt (beta))
  pdf <- c * e * du
  if (log==TRUE){pdf <-log(pdf)}
  return(pdf)
}


#'@rdname BS
#'@importFrom stats pnorm
#'@export


pcbs <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if(!is.numeric(q)||!is.numeric(alpha)||!is.numeric(beta)){
    stop("non-numeric argument to mathematical function")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  x   <- q
  s   <- (x / beta)
  a   <- ((1 / alpha) * (s^(1 / 2) - s^(-1 / 2)))
  cdf <- pnorm(a, 0, 1)
  if (lower.tail == FALSE){cdf <-(1 - cdf)}
  if (log.p == TRUE){cdf <-log(cdf)}
  return(cdf)
}

#'@rdname BS
#'@importFrom stats qnorm
#'@export
#'

qcbs <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE){
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  if (log.p == TRUE){p  <- log(p)}
  if (lower.tail == FALSE){p <- (1 - p)}
  q   <- beta * (((alpha * qnorm(p, 0, 1) / 2) +
                    sqrt(((alpha * qnorm(p, 0, 1) / 2)^2) +
                           1)))^2
  return(q)
}

#'@rdname BS
#'@importFrom stats rnorm
#'@export
#'
rcbs <- function(n, alpha = 1, beta = 1)
{
  if (!is.numeric(n)||!is.numeric(alpha)||!is.numeric(beta))
  {stop("non-numeric argument to mathematical function")}
  if (n == 0){stop("value of n must be greater or equal then 0")}
  if (alpha <= 0){stop("alpha must be positive")}
  if (beta <= 0){stop("beta must be positive")}
  z   <- rnorm(n, 0, 1)
  t   <- (beta/4)*(alpha*z + sqrt((alpha*z)^2 + 4))^2
  return(t)
}




