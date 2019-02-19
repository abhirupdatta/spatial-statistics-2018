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
library("coda")

### function for plotting interpolated surface of a column of a data table
myplot=function(tab,colname){
  
  surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
  dev.new()
  image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
  
}
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)


### Gibbs sampler simulated dataset 3 from lecture 1 ####
data3=read.csv("../data/dataset3.csv")

### data and covariate surface plots
n=nrow(data3)

X=cbind(1,data3$x)
y=data3$y

S=data3[,c("sx","sy")]
dmat=as.matrix(dist(S))


#### Chain length and burn in ####
N=10000
Nb=5001

#### true values ####
truesigs=0.25
truetaus=0.01
truebeta=c(0.2,-0.3)
truephi=2

####### Gibbs sampler (fixing phi at the variogram produced estimate) ##########

### linear regression and empirical variogram to choose phi  ###
ols=lm(y~x,data=data3)
data3$res=ols$residuals

max.dist <- 0.75*max(iDist(data3[,1:2]))
bins <- 20
vario3 <- variog(coords=data3[,1:2], data=data3$res, 
                 uvec=(seq(0, max.dist, length=bins)))
vfit3 <-variofit(vario3, ini.cov.pars=c(0.1,1), ##sigma^2 and 1/phi 
                 cov.model="exponential", minimisation.function="optim",
                 nugget=0.01, weights="equal")
phi=1/vfit3$cov.pars[2]

R=exp(-phi*dmat)
Rinv=chol2inv(chol(R))
gram=chol2inv(chol(t(X)%*%X))
M=gram%*%t(X)

### creating matrices and vectors to store the samples ###
betamat=matrix(0,N,2)
wmat=matrix(0,N,n)
sigsvec=tausvec=rep(0,N)

### initial values ###
beta=rep(0,2)
sigs=taus=1

### hyper-params ###
asig=1
atau=0.5
bsig=btau=2
res=y-X%*%beta

### Gibbs sampler ###
### set.seed(1)
set.seed(12345) ### (picking a different seed to demonstrate convergence)
for(i in 1:N){
    Vw=chol2inv(chol(Rinv/sigs+diag(n)/taus))
    w=as.vector(rmvnorm(1,Vw%*%res/taus,Vw))
    beta=as.vector(rmvnorm(1,M%*%(y-w),taus*gram))
    res=y-X%*%beta
    sigs=rinvgamma(1,asig+n/2,bsig+t(w)%*%Rinv%*%w/2)
    taus=rinvgamma(1,atau+n/2,btau+sum((res-y)^2)/2)
    betamat[i,]=beta
    wmat[i,]=w
    sigsvec[i]=sigs
    tausvec[i]=taus
    if(i %% 10 == 0) print(i)
}

MCMCmat=cbind(betamat,wmat,sigsvec,tausvec)
colnames(MCMCmat)= c(paste0("beta",1:length(beta)),paste0("w",1:length(w)),"sigs","taus")
#write.csv(MCMCmat,"data3MCMCmat.csv",row.names=F,quote=F)
write.csv(MCMCmat,"data3MCMCmat_seed12345.csv",row.names=F,quote=F)
samples12345=MCMCmat

### loading the pre-saved samples as the MCMC takes a while
samples=read.csv("data3MCMCmat.csv")  ### loading the results for seed 1
samples12345=read.csv("data3MCMCmat_seed12345.csv") ### loading the results for seed 12345

##### convergence diagnostics #####
##### trace plots ####
dev.new()
par(mfrow=c(2,3))
plot(samples[,1],lwd=2,main="beta0",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,1],samples12345[,1]))
lines(samples12345[,1],lwd=2,col="blue")

plot(samples[,2],lwd=2,main="beta1",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,2],samples12345[,2]))
lines(samples12345[,2],lwd=2,col="blue")

plot(samples[,n+3],lwd=2,main="sigs",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,n+3],samples12345[,n+3]))
lines(samples12345[,n+3],lwd=2,col="blue")

plot(samples[,n+4],lwd=2,main="taus",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,n+4],samples12345[,n+4]))
lines(samples12345[,n+4],lwd=2,col="blue")

plot(samples[,2+77],lwd=2,main="w77",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,2+77],samples12345[,2+77]))
lines(samples12345[,2+77],lwd=2,col="blue")

plot(samples[,2+498],lwd=2,main="w498",xlab="",ylab="",col="red",
     type="l",ylim=range(samples[,2+498],samples12345[,2+498]))
lines(samples12345[,2+498],lwd=2,col="blue")

##### Gelman-Rubin diagnostics #####
mcmc1=as.mcmc(samples)
mcmc2=as.mcmc(samples12345)
mcmclist=mcmc.list(mcmc1,mcmc2)
gr=gelman.diag(mcmclist)

dev.new()
plot(gr$psrf[-c(1:2,503:504),1],ylab="Gelamn-Rubin shrink factor",xlab="i",ylim=range(gr$psrf[,1]))
points(c(1:2,499:500),gr$psrf[c(1:2,503:504),1],col=c("orange","red","cyan","blue"),cex=1.5,pch=16)
legend("topleft",c("beta0","beta1","sigs","taus","w(s_i)"),bty='n',
       col=c("orange","red","cyan","blue","black"),pch=c(rep(16,4),1),pt.cex=c(rep(1.5,4),1))

gr$mpsrf

### THIS TAKES A LONG TIME ###
gelman.plot(mcmclist) ### this will create plot for all 504 variables ####


#### density plots ####
dev.new()
par(mfrow=c(2,3))
plot(density(samples[-(1:Nb),1]),lwd=2,main="beta0",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),1]),lwd=2,col="blue")

plot(density(samples[-(1:Nb),2]),lwd=2,main="beta1",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),2]),lwd=2,col="blue")

plot(density(samples[-(1:Nb),n+3]),lwd=2,main="sigs",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),n+3]),lwd=2,col="blue")

plot(density(samples[-(1:Nb),n+4]),lwd=2,main="taus",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),n+4]),lwd=2,col="blue")

plot(density(samples[-(1:Nb),2+77]),lwd=2,main="w77",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),2+77]),lwd=2,col="blue")

plot(density(samples[-(1:Nb),2+498]),lwd=2,main="w498",xlab="",ylab="",col="red",type="l")
lines(density(samples12345[-(1:Nb),2+498]),lwd=2,col="blue")


### model comparison using DIC for spatial LM ###
dicgen=function(samples,Nmcmc,Nburn,y){
    #wmat=post$w[Nburn:N,]
    post=samples[Nburn:Nmcmc,]
    n1=Nmcmc-Nburn+1
    n=ncol(samples)-4
    wlist=paste0('w',1:n)
    
    meanmat=t(t(rep(1,n)))%*%post[,"beta1"]+X[,2,drop=F]%*%post[,"beta2"]+t(post[,wlist])  ### Xbeta + w
    res=y-meanmat ### y - Xbeta - w
    qf=colSums(res*res) ##  (y - Xbeta - w)'(y - Xbeta - w))
    
    avgdev=mean(qf/taus)+n*mean(log(post[,"taus"])) ## avg deviance
    
    mumean=mean(post[,"beta1"])+X[,2,drop=F]%*%mean(post[,"beta2"])+colMeans(post[,wlist]) ## average mean
    tausmean=mean(post[,"taus"]) ## average tauusq
    devavg=sum((y-mumean)*(y-mumean)/tausmean)+n*log(tausmean) ## deviance at avg mean and variance 
    
    list(dic=2*avgdev-devavg,pd=avgdev-devavg) ## dic and pd
}


dicgen(samples,N,Nb,y)

### recovering w ####
pbsample <- samples[Nb:N,] ## post burn-in samples

### quantiles of each paraneters ###
qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
qtls[,c(1:2,n+3,n+4)]

what=qtls[2,3:(n+2)]
data3$what=what
myplot(data3,"what")

### plotting kriged surface ###
set.seed(1)
subsample=pbsample[sample(1:nrow(pbsample),100),]

xo=yo=seq(0,1,0.02)
so=expand.grid(xo,yo)

Do=rdist(so,S)

### we can calculate this as phi is fixed here 
c=exp(-phi*Do)
weights=c%*%Rinv

### we will use only 100 posterior samples 
subsample=as.matrix(subsample)
wpredmean=subsample[,3:(n+2)]%*%t(weights)

Xo=cbind(1,0.5*sin(10*so[,1]*so[,2])+1*(0.5-so[,1])^2)

### kriging using composition sampling wo|w,params,y then yo|wo,w,params,y
set.seed(1)
predmat=sapply(1:nrow(so),function(i,c,weights,wpredmean,subsample,Xo){
    wpredvar=pbsample[,n+3]*(1-sum(c[i,]*weights[i,]))
    wo=rnorm(100,wpredmean[,i],sqrt(wpredvar))
    yo=as.vector(subsample[,1:2]%*%Xo[i,]) + wo + rnorm(100,rep(0,100),sqrt(subsample[,n+4]))
    yo
},c,weights,wpredmean,subsample,Xo)

surface_krig_tab=cbind(so,Xo[,2],apply(predmat,2,median),apply(predmat,2,var))
colnames(surface_krig_tab)=c("sx","sy","x","yhat","vyhat")

myplot(surface_krig_tab,"yhat")
myplot(surface_krig_tab,"vyhat")

############ Example: Metropolis Algorithm #############
n=100
theta=3
sigmasq=1
set.seed(1)
Y=rnorm(n,theta,sqrt(sigmasq))

tausq=5
mu=0
posterior_loglikelihood=function(theta,Y,n,sigmasq,mu,tausq){
    -0.5*n*(mean(Y)-theta)^2/sigmasq-0.5*(theta-mu)^2/tausq
}
N=10000
Nb=5001
lambda=0.1
theta_sample_metropolis=flag=rep(0,N)
thetainit=0
for(i in 1:N){
    thetastar=rnorm(1,thetainit,sqrt(lambda))
    logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
    logu=log(runif(1,0,1))
    if(logu <= logr){
        theta_sample_metropolis[i]=thetastar
        flag[i]=1
    }	else	
    {
        theta_sample_metropolis[i]=thetainit
    }
    thetainit=theta_sample_metropolis[i]
}

acceptance_ratio=sum(flag)/length(flag)
acceptance_ratio

### second chain ###
theta_sample_metropolis2=flag2=rep(0,N)
thetainit2=10
for(i in 1:N){
    thetastar=rnorm(1,thetainit,sqrt(lambda))
    logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit2,Y,n,sigmasq,mu,tausq)
    logu=log(runif(1,0,1))
    if(logu <= logr){
        theta_sample_metropolis2[i]=thetastar
        flag2[i]=1
    }	else	
    {
        theta_sample_metropolis2[i]=thetainit2
    }
    thetainit2=theta_sample_metropolis2[i]
}

dev.new()
plot(theta_sample_metropolis, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="Trace plots")
lines(theta_sample_metropolis2, col="blue")

### direct sample ###
theta_sample=rnorm(N-Nb+1,(n*mean(Y)/sigmasq+mu/tausq)/(n/sigmasq+1/tausq),sqrt(1/(n/sigmasq+1/tausq)))

dev.new()
plot(density(theta_sample),xlab="theta",ylab="p(theta)",main="Posterior density")
lines(density(theta_sample_metropolis[Nb:N]), col="red")
lines(density(theta_sample_metropolis2[Nb:N]), col="blue")
legend("topright",legend=c("Direct","MA 1","MA 2"), col=c("black","red","blue"), lwd=rep(1,3))

####### Tuning ########
lambda=10
theta_sample_metropolis_low=flaglow=rep(0,N)
thetainit=0
for(i in 1:N){
    thetastar=rnorm(1,thetainit,sqrt(lambda))
    logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
    logu=log(runif(1,0,1))
    if(logu <= logr){
        theta_sample_metropolis_low[i]=thetastar
        flaglow[i]=1
    }	else	
    {
        theta_sample_metropolis_low[i]=thetainit
    }
    thetainit=theta_sample_metropolis_low[i]
}

acceptance_ratio_low=sum(flaglow)/length(flaglow)
acceptance_ratio_low

plot(theta_sample_metropolis_low, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="lambda=10")

lambda=0.0001
theta_sample_metropolis_high=flaghigh=rep(0,N)
thetainit=0
for(i in 1:N){
    thetastar=rnorm(1,thetainit,sqrt(lambda))
    logr=posterior_loglikelihood(thetastar,Y,n,sigmasq,mu,tausq)-posterior_loglikelihood(thetainit,Y,n,sigmasq,mu,tausq)
    logu=log(runif(1,0,1))
    if(logu <= logr){
        theta_sample_metropolis_high[i]=thetastar
        flaghigh[i]=1
    }	else	
    {
        theta_sample_metropolis_high[i]=thetainit
    }
    thetainit=theta_sample_metropolis_high[i]
}

acceptance_ratio_high=sum(flaghigh)/length(flaghigh)
acceptance_ratio_high

plot(theta_sample_metropolis_high, type="l", col="red", xlab="Iteration i", ylab="theta_i", main="lambda=0.0001")

######  Example 2: Jacobian ######
library(MCMCpack)
n=100
sigmasq=4
set.seed(1)
Y=rnorm(n,0,sqrt(sigmasq))

alpha=2
beta=1
posterior_loglikelihood=function(sigmasq,Y,n,alpha,beta){
    -(alpha+n/2+1)*log(sigmasq)-(sum(Y^2)/2+beta)/sigmasq
}
N=10000
Nb=5001
lambda=0.1
sigmasq_sample_metropolis=flag=rep(0,N)
sigmasqinit=1
for(i in 1:N){
    sigmasqstar=exp(rnorm(1,log(sigmasqinit),sqrt(lambda)))
    logr=posterior_loglikelihood(sigmasqstar,Y,n,alpha,beta)+
        log(sigmasqstar)-posterior_loglikelihood(sigmasqinit,Y,n,alpha,beta)-log(sigmasqinit)
    logu=log(runif(1,0,1))
    if(logu <= logr){
        sigmasq_sample_metropolis[i]=sigmasqstar
        flag[i]=1
    }	else	
    {
        sigmasq_sample_metropolis[i]=sigmasqinit
    }
    sigmasqinit=sigmasq_sample_metropolis[i]
}

acceptance_ratio=sum(flag)/length(flag)
acceptance_ratio

sigmasq_sample=rinvgamma(N-Nb+1,alpha+n/2,beta+sum(Y^2)/2)

dev.new()
plot(density(sigmasq_sample),xlab="theta",ylab="p(theta)",main="Posterior density")
lines(density(sigmasq_sample_metropolis[Nb:N]), col="red")
legend("topright",legend=c("Direct","MA"), col=c("black","red"), lwd=rep(1,2))

############ dataset3 MH for the marginalized model using Nimble ##############
library(nimble)
data3=read.csv("../data/dataset3.csv")

n=nrow(data3)

X=cbind(1,data3$x)
y=data3$y

S=data3[,c("sx","sy")]
dmat=as.matrix(dist(S))

gpcov <- nimbleFunction(
    run = function(dmat=double(2),phi=double(0),n=double(0)){
        returnType(double(2))
        M=matrix(0,n,n)
        for(i in 1:n) 
            for(j in 1:n)
                M[i,j]=exp( -phi * dmat[i,j])
        #M=exp( -ph * dmat)
        return(M)
    })


bayesSpatialLM <- nimbleCode({
    sigs <- 1/invsigmasq
    taus <- 1/invtausq
    G[,] <- gpcov(dmat[,],phi,n)
    V[,] <- sigs*G[,]+ taus * I[,]
    C[,] <- chol(V[,])
    mu[] <- X[,]%*%beta[]
    y[] ~ dmnorm(mu[], cholesky=C[,], prec_param=0)  ## marginalized model
    
    ### parameter priors ###
    beta[] ~ dmnorm(mub[], cholesky=Cb[,], prec_param=0)
    invsigmasq ~ dgamma(2, rate=1)
    invtausq ~ dgamma(0.01, 0.01)
    phi ~  dunif(0,10)
})

constants <- list(n = n, dmat=dmat,
                  ones = rep(1,n), I=diag(n), Cb=1000*diag(2), mub=rep(0,2)
)

data <- list(y = y, X=X)

dimensions = list(G = c(n, n),
                  V = c(n, n),
                  C = c(n, n),
                  Cb = c(2, 2),
                  mub = c(2),
                  mu = c(n),
                  ones = c(n),
                  I = c(n,n),
                  beta = c(2),
                  y = c(n),
                  X= c(n,2),
                  G = c(n,n))

model <- nimbleModel(code=bayesSpatialLM, constants=constants, data=data,
                     dimensions = dimensions,check = FALSE)

Cmodel <- compileNimble(model)  ## compiling model 
modelconf <- configureMCMC(Cmodel,print=TRUE) 
modelconf$addMonitors(c('beta','sigs','taus','phi'))

Rmcmc <- buildMCMC(modelconf)
nimbleOptions(showCompilerOutput = FALSE) 
Cmcmc <- compileNimble(Rmcmc, project = model) ## compiling MCMC
set.seed(1)
Cmcmc$run(100) ## running MCMC

samples <- as.matrix(Cmcmc$mvSamples)
#write.csv(samples,"data3nimble.csv",quote=F,row.names=F)

### loading the pre-saved samples with 10000 iterations
### samples = read.csv("data3nimble.csv") 

N=10000
Nb=5001
pbsample <- samples[Nb:N,] ## post burn-in samples

qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
round(qtls,3)

#### true values ####
truesigs=0.25
truetaus=0.01
truebeta=c(0.2,-0.3)
truephi=2

#### true values and posterior distribution and quantiles ####
dev.new()
plot(density(pbsample[,"beta.1."]),lwd=2,main="beta0",xlab="",ylab="")
abline(v=qtls[,"beta.1."],lwd=2)
abline(v=truebeta[1],col="red",lwd=2)

dev.new()
plot(density(pbsample[,"beta.2."]),lwd=2,main="beta1",xlab="",ylab="")
abline(v=qtls[,"beta.2."],lwd=2)
abline(v=truebeta[2],col="red",lwd=2)

dev.new()
plot(density(pbsample[,"sigs"]),lwd=2,main="sigs",xlab="",ylab="")
abline(v=qtls[,"sigs"],lwd=2)
abline(v=truesigs,col="red",lwd=2)

dev.new()
plot(density(pbsample[,"taus"]),lwd=2,main="taus",xlab="",ylab="")
abline(v=qtls[,"taus"],lwd=2)
abline(v=truetaus,col="red",lwd=2)

dev.new()
plot(density(pbsample[,"phi"]),lwd=2,main="phi",xlab="",ylab="")
abline(v=qtls[,"phi"],lwd=2)
abline(v=truephi,col="red",lwd=2)

#### Posterior recovery of w ####
set.seed(1)
subsample=pbsample[sample(1:nrow(pbsample),100),]
wsamples=apply(subsample,1,function(params,y,X,dmat){
    beta=params[1:2]
    taus=params[7]
    sigs=params[6]
    phi=params[5]
    Rinv=chol2inv(chol(exp(-phi*dmat)))
    Vw=chol2inv(chol(Rinv/sigs+diag(n)/taus))
    w=as.vector(rmvnorm(1,Vw%*%(y-X%*%beta)/taus,Vw))
},y,X,dmat  )

what=apply(wsamples,1,median)
data3$what=what
myplot(data3,"what")


##### spBayes package ######
## BEF data ###
## Data preliminaries
data(BEF.dat)
BEF.dat <- BEF.dat[BEF.dat$ALLBIO02_KGH>0,]
bio <- BEF.dat$ALLBIO02_KGH*0.001;
log.bio <- log(bio)
## Extract the coordinates
coords <- as.matrix(BEF.dat[,c("XUTM","YUTM")])

## Make a surface plot
x.res <- 100; y.res <- 100

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)")
points(coords)

### variogram on the residuals ###
BEF.dat$log.bio=log.bio
ols=lm(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,data=BEF.dat)
BEF.dat$res=ols$residuals

max.dist <- 0.5*max(iDist(BEF.dat[,2:3]))
bins <- 20

vario <- variog(coords=BEF.dat[,2:3], data=BEF.dat$res, 
                 uvec=(seq(0, max.dist, length=bins)))

dev.new()
plot(vario,pch=16)

vfit <-variofit(vario, ini.cov.pars=c(0.04,3/6000), ##sigma^2 and 1/phi 
                cov.model="exponential", minimisation.function="optim",
                nugget=0.08, weights="equal")
vfit

p <- 6 ## This is the number of columns in the design matrix
## Set the prior mean and precision for the regression
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

## For use with bayesGeostatExact, do the following
phi <- 1/vfit$cov.pars[2] ## Set the spatial range (from the variogram)
alpha <- vfit$nugget/vfit$cov.pars[1] ## Set the nugget/partial-sill ratio
sigma.sq.prior.shape <- 2.0 ## Set IG shape for sigma.sq (partial sill)
sigma.sq.prior.rate <- 0.08 ## Set IG scale for sigma.sq (partial sill)

## Run bayesGeostatExact to deliver exact posterior samples
set.seed(1)
sp.exact <- bayesGeostatExact(
    log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
    data=BEF.dat, coords=coords, n.samples=1000,
    beta.prior.mean=beta.prior.mean,
    beta.prior.precision=beta.prior.precision,
    cov.model="exponential",
    phi=phi, alpha=alpha,
    sigma.sq.prior.shape=sigma.sq.prior.shape,
    sigma.sq.prior.rate=sigma.sq.prior.rate,
    sp.effects=FALSE)

##Produce the posterior summaries
round(summary(sp.exact$p.samples)$quantiles,3)

## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
set.seed(1)
bef.sp <- spLM(log.bio~ELEV+SLOPE+SUM_02_TC1+SUM_02_TC2+SUM_02_TC3,
               data=BEF.dat, coords=coords, starting=list("phi"=3/200,"sigma.sq"=0.08,
                                                          "tau.sq"=0.02), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
               priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),
                           "tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)

round(summary(mcmc(bef.sp$p.theta.samples))$quantiles,3)

## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

### DIC using spDiag 
spDiag(bef.sp)$DIC

## Obtain trace plots for regression coefficients
dev.new()
par(mfrow=c(3,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

round(summary(mcmc(bef.sp$p.beta.recover.samples))$quantiles,3)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)

## Plot the spatial residual mean surface and a map of sd's
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords, ols$residuals), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="LM residuals")
surf <- mba.surf(cbind(coords, w.hat.mu), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean spatial effects (w(s))")

## loading shp files for predictions
BEF.shp <- readOGR("../data/BEF-data/BEF_bound.shp")
shp2poly <- BEF.shp@polygons[[1]]@Polygons[[1]]@coords
BEF.poly <- as.matrix(shp2poly)
BEF.grids <- readGDAL("../data/BEF-data/dem_slope_lolosptc_clip_60.img")

## Construct the prediction design matrix for the entire grid extent.
pred.covars <- cbind(BEF.grids[["band1"]], BEF.grids[["band2"]], BEF.grids[["band3"]], BEF.grids[["band4"]], BEF.grids[["band5"]])
pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)


## Extract the coordinates of the BEF bounding polygon vertices and use the pointsInPoly (spBayes) function to obtain the desired subset of the prediction design matrix and associated prediction coordinates (i.e., pixel centroids).
pred.coords <- SpatialPoints(BEF.grids)@coords
pointsInPolyOut <- pointsInPoly(BEF.poly, pred.coords)
pred.covars <- pred.covars[pointsInPolyOut,]
pred.coords <- pred.coords[pointsInPolyOut,]

bef.bio.pred <- spPredict(bef.sp, start=burn.in, thin=2, pred.coords=pred.coords, pred.covars=pred.covars)

## Mapping the predicted values
bef.bio.pred.mu = apply(bef.bio.pred$p.y.predictive.samples,1,mean)
bef.bio.pred.sd = apply(bef.bio.pred$p.y.predictive.samples,1,sd)

surf <- mba.surf(cbind(coords, log.bio), no.X=x.res, no.Y=x.res, extend=TRUE, sp=TRUE)$xyz.est
#surf <- surf [!is.na(over(surf, BEF.shp)),]
surf <- surf [!is.na((over(surf, BEF.shp)))[,1],]
surf <- as.image.SpatialGridDataFrame(surf)
z.lim <- range(surf[["z"]], na.rm=TRUE)

pred.grid <- as.data.frame(list(pred.coords, pred.mu=bef.bio.pred.mu, pred.sd=bef.bio.pred.sd))
coordinates(pred.grid) = c("x", "y")
gridded(pred.grid) <- TRUE
pred.mu.image <- as.image.SpatialGridDataFrame(pred.grid["pred.mu"])

par(mfrow=c(1,2))
image.plot(surf, axes=TRUE, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Log metric tons of biomass")
plot(BEF.shp, add=TRUE)
image.plot(pred.mu.image, zlim=z.lim, col=tim.colors(25), xaxs = "r", yaxs = "r", main="Mean predicted log metric tons of biomass")
plot(BEF.shp, add=TRUE)
