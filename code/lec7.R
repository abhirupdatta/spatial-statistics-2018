##Predictive process example
rm(list=ls())
library(MBA)
library(fields)
library(spBayes)
library(spNNGP)

## first part same as last class (full GP and predictive process)
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


## NNGP (m=15) 
m.s.5 <- spNNGP(y~X-1, coords=coords, starting=starting, method="sequential", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, return.neighbors = TRUE, n.omp.threads=2)


m.s.15 <- spNNGP(y~X-1, coords=coords, starting=starting, method="sequential", n.neighbors=15,
                 tuning=tuning, priors=priors, cov.model=cov.model,
                 n.samples=n.samples, return.neighbors = TRUE, n.omp.threads=2)


## Timing
m.gp$run.time
m.pp.gp.25$run.time
m.pp.gp.64$run.time
m.s.5$run.time
m.s.15$run.time

## Summary cov parameters
round(summary(m.gp$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.25$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.pp.gp.64$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.s.5$p.theta.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.s.15$p.theta.samples)$quantiles[,c(3,1,5)],2)


## Summary random effects
m.gp.w.hat <- apply(m.gp$p.w.recover.samples, 1, median)

m.pp.25.w.hat <- apply(m.pp.gp.25$p.w.recover.samples, 1, median)

m.pp.64.w.hat <- apply(m.pp.gp.64$p.w.recover.samples, 1, median)

m.s.15.w.hat <- apply(m.s.15$p.w.samples, 1, median)

m.s.5.w.hat <- apply(m.s.5$p.w.samples, 1, median)

## Interpolate
surf.w <- mba.surf(cbind(coords, w), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.gp <- mba.surf(cbind(coords, m.gp.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.25 <- mba.surf(cbind(coords, m.pp.25.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.pp.64 <- mba.surf(cbind(coords, m.pp.64.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.nngp.15 <- mba.surf(cbind(coords, m.s.15.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
surf.nngp.5 <- mba.surf(cbind(coords, m.s.5.w.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est

dev.new()
par(mfrow=c(2,3))
image.plot(surf.w, main="True w")
image.plot(surf.gp, main="GP estimated w")
image.plot(surf.pp.25, main="PPGP knots 25 w"); points(m.pp.gp.25$knot.coords, pch=19)
image.plot(surf.pp.64, main="PPGP knots 64 w"); points(m.pp.gp.64$knot.coords, pch=19)
image.plot(surf.nngp.5, main="NNGP (m=5) estimated w")
image.plot(surf.nngp.15, main="NNGP (m=15) estimated w")

dev.new()
par(mfrow=c(2,3))
plot(w, m.gp.w.hat,  xlab="True w", ylab="Full GP Posterior median w")
plot(w, m.pp.25.w.hat,  xlab="True w", ylab="PPGP 36 Posterior median w")
plot(w, m.pp.64.w.hat,  xlab="True w", ylab="PPGP 64 Posterior median w")
plot(w, m.s.15.w.hat,  xlab="True w", ylab="NNGP 5 Posterior median w")
plot(w, m.s.15.w.hat,  xlab="True w", ylab="NNGP 15 Posterior median w")
plot(m.gp.w.hat, m.s.15.w.hat,  xlab="Full GP Posterior median w", ylab="NNGP 15 Posterior median w")

#' 
#' # Harvard Forest canopy height analysis using spNNGP and leaflet package
#' 
#' Here we consider forest canopy height (m) data 
#' measured using the NASA Goddard's LiDAR Hyperspectral 
#' and Thermal (G-LiHT) Airborne Imager over a subset of 
#' Harvard Forest Simes Tract, MA, collected in Summer 2012. 
#' This is a sampling LiDAR system that only records strips 
#' of canopy height across the landscape. We would like to 
#' use the Harvard Forest data to assess if the current density 
#' of LiDAR measurements can be reduced, which would allow for 
#' wider strips to be collected. Ultimately, interest is in 
#' creating wall-to-wall maps of forest canopy height with 
#' associated uncertainty.
#' 
#' Let's load the necessary packages and canopy height data which are part of the `spNNGP` package. Here too, we subset the data and divide it into a model and testing set.
## ---- message=FALSE------------------------------------------------------
library(geoR)
library(raster)
library(leaflet)

data(CHM)

CHM <- CHM[CHM[,3]>0,]

set.seed(1)
mod <- sample(1:nrow(CHM), 25000)
ho <- sample((1:nrow(CHM))[-mod], 10000)

CHM.mod <- CHM[mod,]
CHM.ho <- CHM[ho,]

#' 
#' Let's again start with a `leaflet` basemap then overlay the canopy height data. Recall, `leaflet` maps expect data to be in geographic coordinate system (i.e., longitude and latitude), so we first need reproject the CHM data (just for visualization purposes, we'll fit the model using the projected coordinates).
#' 
## ------------------------------------------------------------------------
chm.r <- rasterFromXYZ(CHM)
proj4string(chm.r) <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
chm.r.ll <- projectRaster(chm.r, crs="+proj=longlat +datum=WGS84")

pal <- colorNumeric(rev(terrain.colors(50)), domain = values(chm.r.ll), na.color = "transparent")

base.map <- leaflet(width="100%") %>%
  addProviderTiles("Esri.WorldImagery", group="Satellite") %>%
  addProviderTiles("Esri.WorldShadedRelief", group="Terrain")

base.map %>%
  addRasterImage(chm.r.ll, colors = pal, opacity = 1, group="Canopy height") %>%
  addLegend("bottomright", pal = pal, values = values(chm.r.ll), opacity = 1, title = "<center>Canopy height (m)</center>") %>%
  addLayersControl(
    baseGroup = c("Satellite", "Terrain"),
    overlayGroups = c("Canopy height"),
    options = layersControlOptions(collapsed = FALSE)
  )

 
#' Let's try and fit a variogram to the data to get a sense of the spatial structure. These `variog` function calculates the $n\times n$ Euclidean distance matrix to construct the empirical variogram. When $n$ is large this will you will likely run out of memory, so you might need to consider only a subset of your data.
#' 
## ---- fig.align="center"-------------------------------------------------
sub <- 1:10000

#note, max intesite distance is ~1.5km
v <- variog(coords=CHM.mod[sub,1:2], data=CHM.mod[sub,3], uvec=(seq(0, 500, length=30))) 

plot(v, xlab="Distance (m)")

#' 
#' Now let's fit some spatial regression models using NNGP random effects.
## ------------------------------------------------------------------------
n.samples <- 1000

starting <- list("phi"=3/50, "sigma.sq"=15, "tau.sq"=2.5)

tuning <- list("phi"=0.05, "sigma.sq"=0.01, "tau.sq"=0.01)

priors <- list("phi.Unif"=c(3/1000, 3/10), "sigma.sq.IG"=c(2, 10), "tau.sq.IG"=c(2, 5))

cov.model <- "exponential"

##Response model 
m.r <- spNNGP(CHM.mod[,3] ~ 1, coords=CHM.mod[,1:2], starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=100)

save(m.r,file="mr.RData")
#load("mr.Rdata")

round(summary(m.r$p.beta.samples)$quantiles[c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],3)

m.r$run.time

#' Now prediction for the holdout set. 
#' 
## ---- fig.align="center"-------------------------------------------------
burn.in <- floor(0.5*n.samples)

p.r <- spPredict(m.r, X.0 = as.matrix(rep(1,nrow(CHM.ho))), coords.0 = CHM.ho[,c("x","y")],
                 start=burn.in, thin=2, n.report=5000, n.omp.threads=2)

##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(2, 10)

cov.model <- "exponential"

g <- 5
theta.alpha <- as.matrix(expand.grid(seq(3/1000,3/10,length.out=g),
  seq(0.001,0.5,length.out=g)))

colnames(theta.alpha) <- c("phi", "alpha")

m.c <- spConjNNGP(CHM.mod[,3] ~ 1, coords=CHM.mod[,1:2], n.neighbors = 10,
                  X.0 = as.matrix(rep(1,nrow(CHM.ho))), coords.0 = CHM.ho[,1:2],
                  k.fold = 5, score.rule = "crps",
                  n.omp.threads = 2,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.c$sigma.sq.hat
m.c$theta.alpha
tau.sq <- m.c$theta.alpha[2]*m.c$sigma.sq.hat
tau.sq

round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],3)

dev.new()
crps.surf <- mba.surf(m.c$k.fold.scores[,c("phi","alpha","crps")], no.X=100, no.Y=100)$xyz.est
image.plot(crps.surf, xlab="phi", ylab="alpha=tau^2/sigma^2", main="CRPS (lower is better)")
points(m.c$theta.alpha, col="white", pch=2)

#### run times ###
m.r$run.time
m.c$run.time

#### true and predicted y scatterplot
dev.new()
par(mfrow=c(1,2))
plot(CHM.ho[,3], y.hat.r, main="Response NNGP model",
     xlab="True canopy height", ylab="Posterior predictive distribution mean")
plot(CHM.ho[,3],m.c$y.0.hat, main="Conjugate NNGP model",
     xlab="True canopy height", ylab="Posterior predictive distribution mean")

#### out of sample RMSPE and CP
y.hat.r <- apply(p.r$p.y.0, 1, mean)

### coverage probability ###
cp.r=mean((apply(p.r$p.y.0, 1, quantile, 0.025) < CHM.ho[,3]) & 
            (apply(p.r$p.y.0, 1, quantile, 0.975) > CHM.ho[,3]))

cp.r

### coverage probability ###
cp.c=mean((abs(CHM.ho[,3] - y.hat.c) <  qnorm(0.975)*sqrt(m.c$y.0.hat.var)))
cp.c

### RMSPE 
sqrt(mean((CHM.ho[,3]- y.hat.r)^2))

#### out of sample RMSPE and CP
y.hat.c <- m.c$y.0.hat

### RMSPE 
sqrt(mean((CHM.ho[,3]- y.hat.c)^2))
mean(CHM.ho[,3]) ## showing the scale of y to understand how the relative scale of RMSPE

##########################################
## Now something really BIG-N just for fun
set.seed(1)
mod <- sample(1:nrow(CHM), 500000)

CHM.mod <- CHM[mod,]

theta.alpha <- c(0.07, 0.13)
names(theta.alpha) <- c("phi", "alpha")

m.c.big <- spConjNNGP(CHM.mod[,3] ~ 1, coords=CHM.mod[,1:2], n.neighbors = 10,
                      theta.alpha = theta.alpha,
                      n.omp.threads = 2,
                      sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.c.big$run.time


