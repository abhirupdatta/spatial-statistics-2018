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

#### Monte Carlo integration ####
n=1000
set.seed(1)
x.samples=rnorm(n)

dev.new()
plot(density(x.samples))
lines(seq(-4,4,length=200),dnorm(seq(-4,4,length=200)),col="red")

mean(exp(x.samples))

dev.new()
plot(density(sin(x.samples)))
quantile(sin(x.samples),c(0.025,0.975))


#### Composition-sampling univariate Normal-Inverse Gamma ####
theta=1
sigsq=1
n=1000

set.seed(123)
y=rnorm(n,theta,sqrt(sigsq))

N=10000
sigsq.samples=rinvgamma(N,(n-1)/2,sum((y-mean(y))^2)/2)
theta.samples=rnorm(N,mean(y),sqrt(sigsq.samples/n))

dev.new()
plot(density(sigsq.samples))
abline(v=sigsq)
mean(sigsq.samples)
quantile(sigsq.samples,c(0.025,0.5,0.975))

plot(density(theta.samples))
abline(v=theta)
mean(theta.samples)
quantile(theta.samples,c(0.025,0.5,0.975))

########## Example: Gibbs' sampling ###########
library(mvtnorm)
n=100
set.seed(1)
X=rnorm(n,0,1)
Y=2+5*X+rnorm(n,0,1)
gamma0=100
gamma1=10

N=10000
Nb=5001

beta0_sample_gibbs=beta1_sample_gibbs=rep(0,N)
beta1init=10
for(i in 1:N){
  beta0_sample_gibbs[i]=rnorm(1,sum(Y-X*beta1init)/(n+1/gamma0),sqrt(1/(n+1/gamma0)))
  beta0init=beta0_sample_gibbs[i]
  beta1_sample_gibbs[i]=rnorm(1,sum(X*(Y-beta0init))/(sum(X^2)+1/gamma1),sqrt(1/(sum(X^2)+1/gamma1)))
  beta1init=beta1_sample_gibbs[i]
}

Xmat=cbind(rep(1,n),X)
betavar=solve(t(Xmat)%*%Xmat+diag(1/c(gamma0,gamma1)))
beta_sample=rmvnorm(N-Nb+1, betavar%*%t(Xmat)%*%Y, betavar)

dev.new()
plot(density(beta_sample[,1]),xlab="theta",ylab="p(theta)",main="Posterior density beta_0")
lines(density(beta0_sample_gibbs[Nb:N]), col="red")
legend("topright",legend=c("Direct","Gibbs'"), col=c("black","red"), lwd=rep(1,2))

dev.new()
plot(density(beta_sample[,2]),xlab="theta",ylab="p(theta)",main="Posterior density beta_1")
lines(density(beta1_sample_gibbs[Nb:N]), col="red")
legend("topright",legend=c("Direct","Gibbs'"), col=c("black","red"), lwd=rep(1,2))


### Gibbs sampler simulated dataset 3 from lecture 1 ####
data3=read.csv("../data/dataset3.csv")

### function for plotting interpolated surface of a column of a data table
myplot=function(tab,colname){
  
  surf <- mba.surf(tab[,c("sx","sy",colname)], no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
  dev.new()
  image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))
  
}

### data and covariate surface plots
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

n=nrow(data3)

X=cbind(1,data3$x)
y=data3$y

S=data3[,c("sx","sy")]
dmat=as.matrix(dist(S))


#### Chain length and burn in ####
N=10000
Nb=5001

#### true values and posterior distribution and quantiles ####
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
set.seed(1)
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
samples=MCMCmat

### loading the pre-saved samples as the MCMC takes a while
### samples=read.csv("data3MCMCmat.csv")

pbsample <- samples[Nb:N,] ## post burn-in samples

qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
qtls[,c(1:2,n+3,n+4)]

dev.new()
plot(density(pbsample[,1]),lwd=2,main="beta0",xlab="",ylab="")
abline(v=qtls[,1],lwd=2)
abline(v=truebeta[1],col="red",lwd=2)

dev.new()
plot(density(pbsample[,2]),lwd=2,main="beta1",xlab="",ylab="")
abline(v=qtls[,2],lwd=2)
abline(v=truebeta[2],col="red",lwd=2)

dev.new()
plot(density(pbsample[,n+3]),lwd=2,main="sigs",xlab="",ylab="")
abline(v=qtls[,n+3],lwd=2)
abline(v=truesigs,col="red",lwd=2)

dev.new()
plot(density(pbsample[,n+4]),lwd=2,main="taus",xlab="",ylab="")
abline(v=qtls[,n+4],lwd=2)
abline(v=truetaus,col="red",lwd=2)

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
  logr=posterior_loglikelihood(sigmasqstar,Y,n,alpha,beta)+log(sigmasqstar)-posterior_loglikelihood(sigmasqinit,Y,n,alpha,beta)-log(sigmasqinit)
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

### loading the pre-saved samples
### samples = read.csv("data3nimble.csv")

pbsample <- samples[Nb:N,] ## post burn-in samples

qtls <- apply(pbsample,2,quantile,c(0.025,0.5,0.95))
qtls

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
