require(numDeriv)
require(flexsurv)
require(gamlss)

MO.RBS <- function (mu.link = "log", sigma.link="log", nu.link = "logit")
{
  mstats <- checklink(   "mu.link", "Marshall Olkin Birnbaum-Saunders", substitute(mu.link), 
                         c("1/mu^2","sqrt", "log", "identity","own"))
  dstats <- checklink("sigma.link", "Marshall Olkin Birnbaum-Saunders", substitute(sigma.link), 
                      c("sqrt","log", "identity", "own"))
  vstats <- checklink(   "nu.link", "Marshall Olkin Birnbaum-Saunders", substitute(nu.link),    
                         c("sqrt", "log", "identity", "own"))
  structure(
    list(family = c("MO.RBS", "Marshall Olkin Reparameterized Birnbaum-Saunders"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         
         dldm = function(y,mu,sigma,nu){ #----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiMO.RBS(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiMO.RBS(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu){#----------------------------------------------------- ok  
           lpdf<-function(t,mu,x,nu){log(dauxiMO.RBS(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiMO.RBS(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,x){log(dauxiMO.RBS(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok 
           lpdf<-function(t,mu,sigma,x){log(dauxiMO.RBS(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
         d2ldmdd = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiMO.RBS(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,x,nu){log(dauxiMO.RBS(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiMO.RBS(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiMO.RBS(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         d2ldddv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiMO.RBS(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiMO.RBS(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,...) 
         { 
           -2*dMO.RBS(y,mu,sigma,nu,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pMO.RBS", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)) ,
         #mu.initial = expression(mu <- rep(mean(y), length(y))),
         #sigma.initial = expression(sigma <- rep(1.2, length(y))), #OK
         mu.initial = expression({mu <- y + mean(y)/2 }),
         sigma.initial = expression({sigma <- rep(1,length(y)) }),
         nu.initial = expression(nu <- rep(0.45,length(y))), #)k
         mu.valid = function(mu) all(mu > 0), 
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu)  all(nu >= 0.01) && all(nu <= .99),
         y.valid = function(y)  all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dMO.RBS<- function(x, mu = 2, sigma = 2, nu = 0.5, log = FALSE){
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  g <- dRBS(x,mu=mu,sigma=sigma)
  G <- pRBS(x,mu=mu,sigma=sigma)
  f <- (nu*g)/(1-(1-nu)*(1-G))^2
  if(log==FALSE) fy<-f else fy<-log(f)
  fy
}    
#-----------------------------------------------------------------  
pMO.RBS<- function(q, mu = 2, sigma = 2, nu = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", ""))  
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
  G <- pRBS(q,mu=mu,sigma=sigma)
  FF1 <- 1- (nu*(1-G))/(G+nu*(1-G))
  
  if(lower.tail==TRUE) cdf<-FF1 else cdf<- 1-FF1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#-----------------------------------------------------------------  
hMO.RBS <-function(x, mu = 2, sigma = 0.5, nu = 0.5){ 
  hazard=dMO.RBS(x,mu,sigma,nu)/(1-pMO.RBS(x,mu,sigma,nu))
  hazard
}
#-----------------------------------------------------------------  
qMO.RBS <-  function(p, mu=2, sigma=0.5, nu=0.5, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", ""))  
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  u1 <- 1-((1-p)/(1+(nu-1)*p))
  q <- qRBS(u1,mu=mu,sigma=sigma)
  q
}
#-----------------------------------------------------------------  
rMO.RBS <- function(n, mu=2, sigma=0.5, nu=0.5){
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
 # if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", ""))  
  uni<- runif(n = n,0,1)
  r <- qMO.RBS(uni,mu =mu, sigma =sigma, nu=nu)
  r
}
#-----------------------------------------------------------------  


dauxiMO.RBS <- function(t,mu,sigma,nu){ 
  g <-  dRBS(t,mu=mu,sigma=sigma)
  G <- pRBS(t,mu=mu,sigma=sigma)
  f <- (nu*g)/(1-(1-nu)*(1-G))^2
  f}


