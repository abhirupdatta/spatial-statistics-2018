##### Code for Lecture 1 #####
library(MBA)
library(fields)
library(classInt)
library(geoR)
library(RColorBrewer)

### dataset 1 ###
data1=read.csv("../data/dataset1.csv")

col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

### plot raw data pounts ###
dev.new()
plot(data1$sx,data1$sy,pch=16,xlab="Easting (m)", ylab="Northing (m)",
     col=findColours(classIntervals(data1$y, n=10, style="equal"), col.pal))
dev.new()
plot(data1$sx,data1$sy,pch=16,xlab="Easting (m)", ylab="Northing (m)",
     col=findColours(classIntervals(data1$x, n=10, style="equal"), col.pal))

### function for plotting interpolated surface of a column of a data table
myplot=function(tab,colname){
  
  surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
  dev.new()
  image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
  
}

## plotting surfaces ##
dev.new()
myplot(data1,"y")
myplot(data1,"x")

### linear regression on x
lmobj1=lm(y~x,data=data1)
data1$res=lmobj1$residuals

myplot(data1,"res")


### datasets 2 and 3
data2=read.csv("../data/dataset2.csv")
data3=read.csv("../data/dataset3.csv")

myplot(data2,"y")
myplot(data3,"y")

lmobj2=lm(y~x,data=data2)
lmobj3=lm(y~x,data=data3)

data2$res=lmobj2$residuals
data3$res=lmobj3$residuals

myplot(data2,"res")
myplot(data3,"res")

### empirical variograms ###
max.dist <- 0.75*max(rdist(data1[,1:2]))
bins <- 20

dev.new()
vario1raw <- variog(coords=data1[,1:2], data=data1$y, uvec=(seq(0, max.dist, length=bins)))
plot(vario1raw,pch=16)

dev.new()
vario1 <- variog(coords=data1[,1:2], data=data1$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario1,pch=16)

dev.new()
vario2 <- variog(coords=data2[,1:2], data=data2$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario2,pch=16)

dev.new()
vario3 <- variog(coords=data3[,1:2], data=data3$res, uvec=(seq(0, max.dist, length=bins)))
plot(vario3,pch=16)

#### adding the coordinates into the regression
lmobj2s=lm(y~x+sx+sy,data=data2)
data2$res2=lmobj2s$residuals
myplot(data2,"res2")

dev.new()
vario2s <- variog(coords=data2[,1:2], data=data2$res2, uvec=(seq(0, max.dist, length=bins)))
plot(vario2s,pch=16)

lmobj3s=lm(y~x+sx+sy,data=data3)
data3$res2=lmobj3s$residuals
myplot(data3,"res2")

dev.new()
vario3s <- variog(coords=data3[,1:2], data=data3$res2, uvec=(seq(0, max.dist, length=bins)))
plot(vario3s,pch=16)