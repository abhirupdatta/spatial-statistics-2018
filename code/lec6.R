##Predictive process example
rm(list=ls())
library(MBA)
library(fields)
library(spBayes)

##Simulated some data
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 500
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))

##Set up spLM call
n.samples <- 2000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.05, "sigma.sq"=0.05, "tau.sq"=0.05)

priors <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
               "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.report <- 100
verbose <- TRUE

##Call full GP and predictive process GP mosels
burn.in <- floor(0.75*n.samples)

##Full GP
m.gp <- spLM(y~X-1, coords=coords, starting=starting,
             tuning=tuning, priors=priors, cov.model=cov.model,
             n.samples=n.samples, verbose=verbose, n.report=n.report)

m.gp <- spRecover(m.gp, start=burn.in, thin=2)

## PP GP with 36 knots
m.pp.gp.25 <- spLM(y~X-1, coords=coords, knots=c(5,5,0.1), starting=starting,
                   tuning=tuning, priors=priors, cov.model=cov.model,
                   n.samples=n.samples, verbose=verbose, n.report=n.report)

m.pp.gp.25 <- spRecover(m.pp.gp.25, start=burn.in, thin=2)

## PP GP with 64 knots
m.pp.gp.64 <- spLM(y~X-1, coords=coords, knots=c(8,8,0.1), starting=starting,
                   tuning=tuning, priors=priors, cov.model=cov.model,
                   n.samples=n.samples, verbose=verbose, n.report=n.report)

m.pp.gp.64 <- spRecover(m.pp.gp.64, start=burn.in, thin=2)

## Timing
m.gp$run.time
m.pp.gp.25$run.time
m.pp.gp.64$run.time

## Summary cov parameters
round(summary(m.gp$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.25$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.64$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

## DIC 
spDiag(m.gp)$DIC
spDiag(m.pp.gp.25)$DIC
spDiag(m.pp.gp.64)$DIC

## Summary random effects
m.gp.w.hat <- apply(m.gp$p.w.recover.samples, 1, median)

m.pp.25.w.hat <- apply(m.pp.gp.25$p.w.recover.samples, 1, median)

m.pp.64.w.hat <- apply(m.pp.gp.64$p.w.recover.samples, 1, median)


## Interpolate
surf.w <- mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.gp <- mba.surf(cbind(coords, m.gp.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.25 <- mba.surf(cbind(coords, m.pp.25.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.64 <- mba.surf(cbind(coords, m.pp.64.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est

dev.new()
par(mfrow=c(2,2))
image.plot(surf.w, main="True w")
image.plot(surf.gp, main="GP estimated w")
image.plot(surf.pp.25, main="PPGP knots 25 w"); points(m.pp.gp.25$knot.coords, pch=19)
image.plot(surf.pp.64, main="PPGP knots 64 w"); points(m.pp.gp.64$knot.coords, pch=19)
