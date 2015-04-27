## construct confidence intervals from Wald var-cov matrices:
set.seed(101)
dd <- data.frame(y=1:20+2*rnorm(20),x=1:20)
m1 <- lm(y~x,data=dd)
library("MASS")
## pick 1000 MVN values
pars <- mvrnorm(1000,mu=coef(m1),Sigma=vcov(m1))
pred1 <- apply(pars,1,function(x) x[1]+x[2]*dd$x)
quant1 <- t(apply(pred1,1,quantile,c(0.025,0.25,0.5,0.75,0.975)))
matplot(dd$x,quant1,type="l",
        col=c(2,3,1,3,2),
        lty=c(2,3,1,3,2))

matplot(dd$x,pred1,col=adjustcolor("black",alpha=0.02),
        type="l",lty=1)
matlines(dd$x,quant1,type="l",
        col=c(2,4,1,4,2),
        lty=c(2,4,1,4,2))

## do the same thing for MCMC results
