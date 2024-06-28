
#Packages

require(gamlss)
require(gamlss.cens)
require(survival)

#Loading the proposed models
source("Model-GAMLSS/MORBS-GAMLSS.r")
source("Model-GAMLSS/MOWEI-GAMLSS.r")
source("Model-GAMLSS/RBS-GAMLSS.r")


dados <- read.table("efron.txt",header = T)



## Weibull model
fitW.Null <- gamlss(Surv(tempo,cens)~1,
                    sigma.formula =~1,family=cens(WEI3)
                    ,n.cyc=200,data=dados)
summary(fitW.Null,type = "qr", 
        robust=T,  hessian.fun = "R")
fitW.Media <- gamlss(Surv(tempo,cens)~x1,
                     sigma.formula =~1,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Media,type = "qr", 
        robust=T,  hessian.fun = "R")
fitW.Shape <- gamlss(Surv(tempo,cens)~1,
                     sigma.formula =~x1,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Shape,type = "qr", 
        robust=T,  hessian.fun = "R")
fitW.Full <- gamlss(Surv(tempo,cens)~x1,
                    sigma.formula =~x1,family=cens(WEI3),n.cyc=200,data=dados)
summary(fitW.Full,type = "qr", 
        robust=T,  hessian.fun = "R")


## RBS Model
fitRBS.Null <- gamlss(Surv(tempo,cens)~1,
                      sigma.formula =~1,family=cens(RBS)
                      ,n.cyc=200,data=dados)
summary(fitRBS.Null,type = "qr", 
        robust=T,  hessian.fun = "R")

fitRBS.Media <- gamlss(Surv(tempo,cens)~x1,
                       sigma.formula =~1,family=cens(RBS),n.cyc=200,data=dados)
summary(fitRBS.Media,type = "qr", 
        robust=T,  hessian.fun = "R")

fitRBS.Shape <- gamlss(Surv(tempo,cens)~1,
                       sigma.formula =~x1,family=cens(RBS),n.cyc=200,data=dados)
summary(fitRBS.Shape,type = "qr", 
        robust=T,  hessian.fun = "R")

fitRBS.Full <- gamlss(Surv(tempo,cens)~x1,
                      sigma.formula =~x1,family=cens(RBS),n.cyc=200,data=dados)
summary(fitRBS.Full,type = "qr", 
        robust=T,  hessian.fun = "R")


## MO-W

fitMOW.Null <- gamlss(Surv(tempo,cens)~1,
                      sigma.formula =~1,family=cens(MO.WEI)
                      ,n.cyc=2000,data=dados)
summary(fitMOW.Null,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMOW.Media <- gamlss(Surv(tempo,cens)~x1,
                       sigma.formula =~1,family=cens(MO.WEI),n.cyc=1250,data=dados)
summary(fitMOW.Media,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMOW.Shape <- gamlss(Surv(tempo,cens)~1,
                       sigma.formula =~x1,family=cens(MO.WEI),n.cyc=1350,data=dados)
summary(fitMOW.Shape,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMOW.Full <- gamlss(Surv(tempo,cens)~x1,
                      sigma.formula =~x1,family=cens(MO.WEI),n.cyc=1200,
                      data=dados)
summary(fitMOW.Full,type = "qr", 
        robust=T,  hessian.fun = "R")




## MO-RBS

fitMORBS.Null <- gamlss(Surv(tempo,cens)~1,
                        sigma.formula =~1,family=cens(MO.RBS)
                        ,n.cyc=1510,data=dados,
                        nu.start = 0.9,
                        mu.start = fitW.Null$mu.fv)
summary(fitMORBS.Null,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMORBS.Media <- gamlss(Surv(tempo,cens)~x1,
                         sigma.formula =~1,family=cens(MO.RBS),
                         n.cyc=50,data=dados)
summary(fitMORBS.Media,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMORBS.Shape <- gamlss(Surv(tempo,cens)~1,
                         sigma.formula =~x1,family=cens(MO.RBS),
                         n.cyc=1250,data=dados)
summary(fitMORBS.Shape,type = "qr", 
        robust=T,  hessian.fun = "R")

fitMORBS.Full <- gamlss(Surv(tempo,cens)~x1,
                        sigma.formula =~x1,family=cens(MO.RBS),
                        n.cyc=1500,data=dados)
summary(fitMORBS.Full,type = "qr", 
        robust=T,  hessian.fun = "R")



