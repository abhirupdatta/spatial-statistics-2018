library(classInt)
library(MBA)
library(fields)
library(geoR)
library(RColorBrewer)

###################################################
### dataset 3 ###
data3=read.csv("../data/dataset3.csv")

### function for plotting interpolated surface of a column of a data table
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

myplot=function(tab,colname){
  
  surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
  dev.new()
  image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
  
}

myplot(data3,"y")
myplot(data3,"x")


lmobj3=lm(y~x,data=data3)

data3$res=lmobj3$residuals

myplot(data3,"res")

### empirical variograms ###
max.dist <- 0.75*max(rdist(data3[,1:2]))
bins <- 20

dev.new()
vario3 <- variog(coords=data3[,1:2], data=data3$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario3,pch=16)

#### adding the coordinates into the regression
lmobj3s=lm(y~x+sx+sy,data=data3)
data3$res2=lmobj3s$residuals
myplot(data3,"res2")

dev.new()
vario3s <- variog(coords=data3[,1:2], data=data3$res2, uvec=(seq(0, max.dist, length=bins)))
plot(vario3s,pch=16)

### analysis using a Gaussian Process ###
### spatial mle ###
mle <- likfit(coords=data3[,1:2], data=data3[,4], trend = trend.spatial(~x,data3),
  ini.cov.pars=c(0.12,0.25),nugget = 0.02,cov.model="exponential",nospatial=TRUE)

mle

## model comparison ##
mle$AIC
mle$BIC

mle$nospatial$AIC
mle$nospatial$BIC

### in sample predictions ##
sp.pred.obj <- krige.conv(coords=data3[,1:2], data=data3[,4],
  locations=data3[,1:2],krige=krige.control(type.krige="OK",obj.model=mle,
  trend.d=trend.spatial(~x,data3),trend.l=trend.spatial(~x,data3)))

sp.pred=sp.pred.obj$predict
dev.new()
plot(data3$y,sp.pred)
max(abs(data3$y - sp.pred))

### in sample predictions (signal = TRUE) ##
sp.pred.obj <- krige.conv(coords=data3[,1:2], data=data3[,4],
  locations=data3[,1:2],krige=krige.control(type.krige="OK",obj.model=mle,
  trend.d=trend.spatial(~x,data3),trend.l=trend.spatial(~x,data3)),output=output.control(signal=TRUE))

sp.pred=sp.pred.obj$predict
data3$res3=data3$y-sp.pred
myplot(data3,"res3")

dev.new()
vario3s <- variog(coords=data3[,1:2], data=data3$res3, uvec=(seq(0, max.dist, length=bins)))
plot(vario3s,pch=16)



###################################################
### WEF data
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

### logDBH quantile based color coding of the locations
quant <- classIntervals(logDBH, n=5, style="quantile")
brks <- round(quant$brks, 2)
quant <- classIntervals(logDBH, n=5, style="fixed",
                        fixedBreaks = brks)

quant.col <- findColours(quant, col.pal)

dev.new()
plot(coords, col=quant.col, pch=19, cex=1.5, main="", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=attr(quant.col, "palette"), 
       legend=names(attr(quant.col, "table")), bty="n",cex=1.3)

### plot of interpolated surface using mba package ###
surf <- mba.surf(cbind(coords,logDBH), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

### plot of Species type ###
spnum=as.numeric(WEF.in$Species)
col.pal2 <- col.br(length(unique(spnum)))

dev.new()
plot(coords, col=col.pal2[spnum], pch=19, cex=1.5, main="", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=col.pal2, 
       legend=levels(WEF.in$Species), bty="n")

### Linear regression ###
lm.logDBH <- lm(logDBH~Species, data=WEF.in)
summary(lm.logDBH)
logDBH.resid <- resid(lm.logDBH)

surf <- mba.surf(cbind(coords,logDBH.resid), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

### variogram of raw data and residuals ###
max.dist=0.5*max(rdist(coords))
bins=20

dev.new()
vario.logDBH <- variog(coords=coords, data=logDBH, uvec=(seq(5, max.dist, length=bins)))
plot(vario.logDBH,pch=16)


dev.new()
vario.logDBH.resid <- variog(coords=coords, data=logDBH.resid, uvec=(seq(0, max.dist, length=bins)))
plot(vario.logDBH.resid,pch=16)


### spatial mle ###
mle <- likfit(coords=coords, data=logDBH, trend = trend.spatial(~Species,WEF.in), ini.cov.pars=c(0.1,100),
              nugget = 0.25,cov.model="exponential",nospatial=TRUE)

mle

### in sample predictions (signal = TRUE) ##
sp.pred.obj <- krige.conv(coords=coords, data=logDBH,
    locations=WEF.in[,c("East_m","North_m")],krige=krige.control(type.krige="OK",obj.model=mle,
    trend.d=trend.spatial(~Species,WEF.in),trend.l=trend.spatial(~Species,WEF.in)),output=output.control(signal=TRUE))

sp.pred=sp.pred.obj$predict
logDBH.sp.resid=logDBH-sp.pred

surf <- mba.surf(cbind(coords,logDBH.sp.resid), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

dev.new()
vario.logDBH.sp.resid <- variog(coords=coords, data=logDBH.sp.resid, uvec=(seq(0, max.dist, length=bins)))
plot(vario.logDBH.sp.resid,pch=16)


## model comparison ##
mle$AIC
mle$BIC

mle$nospatial$AIC
mle$nospatial$BIC

### out of sample predictions ###
### RMSPE ###
krig_mlefit=krige.conv(coords=coords, data=logDBH,
    locations=WEF.out[,c("East_m","North_m")],krige=krige.control(type.krige="OK",obj.model=mle,
    trend.d=trend.spatial(~Species,WEF.in),trend.l=trend.spatial(~Species,WEF.out)))

pred_spatial=krig_mlefit$predict
rmspe_spatial=sqrt(mean((pred_spatial-WEF.out$logDBH)^2))

pred_lm=as.vector(as.matrix(trend.spatial(~Species,WEF.out))%*%lm.logDBH$coefficients)
rmspe_lm=sqrt(mean((pred_lm-WEF.out$logDBH)^2))

rmspe_spatial
rmspe_lm


### CP ###
CI_spatial=pred_spatial+1.96*sqrt(krig_mlefit$krige.var)%*%t(c(-1,1))  ## confidence interval ##
CP_spatial=mean(CI_spatial[,1]<WEF.out$logDBH & CI_spatial[,2]>WEF.out$logDBH) ## coverage probability ##
CIW_spatial=mean(CI_spatial[,2]-CI_spatial[,1]) ## confidence interval width ##

CP_spatial
CIW_spatial

N=nrow(WEF.out)
CI_lm=pred_lm+1.96*summary(lm.logDBH)$sigma*cbind(-rep(1,N),rep(1,N))
CP_lm=mean(CI_lm[,1]<WEF.out$logDBH & CI_lm[,2]>WEF.out$logDBH)
CIW_lm=mean(CI_lm[,2]-CI_lm[,1])

CP_lm
CIW_lm

### kriged surface ##
## prediction locations ###
WEF.pred=read.csv("../data/WEFpred.csv")

krigsurf_mlefit=krige.conv(coords=coords, data=logDBH,
    locations=WEF.pred[,c("East_m","North_m")],krige=krige.control(type.krige="OK",obj.model=mle,
    trend.d=trend.spatial(~Species,WEF.in),trend.l=trend.spatial(~Species,WEF.pred)))

pred=krigsurf_mlefit$predict
predsd=sqrt(krigsurf_mlefit$krige.var)

predsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],pred), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

predsdsurf <- mba.surf(cbind(WEF.pred[,c("East_m","North_m")],predsd), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(predsdsurf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=rev(terrain.colors(25)))

