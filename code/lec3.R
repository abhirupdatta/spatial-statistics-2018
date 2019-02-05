library("ggplot2")
library("RColorBrewer")
library("akima")
library("fields")
library("ggmap")
library("maps")
library("spBayes")
library("MBA")
library("classInt")
library("plotrix")
library("geoR")
library("sp")
library("maptools")
library("rgdal")
library("classInt")
library("lattice")
library("raster")
library("sf")
library("mvtnorm")
library("MCMCpack")

################## part a ########################

###################################################
### understanding prediction variances of WEF data from Lecture 2
###################################################
WEF.dat=read.csv("../data/WEFsmall.csv")

WEF.dat$logDBH=log(WEF.dat$DBH)

set.seed(123)
ind=sample(1:nrow(WEF.dat),100,replace=FALSE)

### holdout data to assess RMSPE ###
WEF.out=WEF.dat[ind,]
WEF.in=WEF.dat[-ind,]
rm("WEF.dat")

### diameter at breast height for the trees
logDBH <- WEF.in$logDBH

coords <- as.matrix(WEF.in[,c("East_m","North_m")])

col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow",  "orange", "red"))
col.pal <- col.br(5)

### spatial mle ###
mle <- likfit(coords=coords, data=logDBH, trend = trend.spatial(~Species,WEF.in), ini.cov.pars=c(0.1,40),
              nugget = 0.25,cov.model="exponential",nospatial=TRUE)

mle


### kriged surface ##
## prediction locations ###
WEF.pred=read.csv("../data/WEFpred.csv")

krigsurf1=krige.conv(coords=coords, data=logDBH,
  locations=WEF.pred[,c("East_m","North_m")],krige=krige.control(type.krige="OK",obj.model=mle,
  trend.d=trend.spatial(~Species,WEF.in),trend.l=trend.spatial(~Species,WEF.pred)))

pred=krigsurf1$predict
predvar=krigsurf1$krige.var

predsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],pred), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

predvarsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predvar), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predvarsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=rev(terrain.colors(25)))


#### recovering w (setting trend.l to be zero) ####
krigsurf2=krige.conv(coords=coords, data=logDBH,
  locations=WEF.pred[,c("East_m","North_m")],krige=krige.control(type.krige="OK",obj.model=mle,
  trend.d=trend.spatial(~Species,WEF.in),trend.l=0*trend.spatial(~Species,WEF.pred))
  ,output=output.control(signal=TRUE))

predw=krigsurf2$predict
predwvar=krigsurf2$krige.var

predwsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predw), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predwsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

predwvarsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predwvar), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predwvarsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=rev(terrain.colors(25)))
points(coords,pch=16)

#### variance of the predicted linear regression mean ####
Xpred=as.matrix(unclass(trend.spatial(~Species,WEF.pred)))
predreg=Xpred%*%mle$beta
predregvar=diag(Xpred%*% mle$beta.var %*% t(Xpred))

predregsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predreg), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predregsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

predregvarsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predregvar), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predregvarsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=rev(terrain.colors(25)))

#### why variance of the linear regression component is high in some patches ####
WEF.pred$Species[which(predregvar>0.02)]
length(which(predregvar>0.02))
table(WEF.pred$Species)

table(WEF.in$Species)

### plot of Species type ###
spnum=as.numeric(WEF.in$Species)
col.pal2 <- col.br(length(unique(spnum)))

dev.new()
plot(coords, col=col.pal2[spnum], pch=19, cex=1.5, main="", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=col.pal2, 
       legend=levels(WEF.in$Species), bty="n")

predregvarsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predregvar), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predregvarsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=rev(terrain.colors(25)))
points(WEF.pred[which(WEF.pred$Species=="GF"),1:2],col="cyan",pch=16,cex=2)