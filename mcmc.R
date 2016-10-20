verbose=TRUE
#setwd('~/r/flu/yilin/May_version/')
library(Rcpp)
library(RcppArmadillo)
if(verbose)cat('Compiling C++ now\n')
Rcpp::sourceCpp("code/C_functions.cpp")
if(verbose)cat('Finished compiling C++ already\n')
if(verbose)cat('Reading R functions\n')
source("code/R_functions.r")


if(verbose)cat('Reading data and inits from file\n')
dataset=import_data()
theta = import_inits(guided=FALSE)
rm(MATlow,MATtime,MATupp,cohort,mohphat,inits)
oneiterationplots(theta,dataset,'initialplots.pdf')
#logLikelihood(theta$sigma,theta$phi, theta$B, theta$L, theta$M, theta$S, theta$I, theta$T, data$low, data$upp, data$time, data$phat, data$cohort)

BURNINs=2000
MCMCits=10000 # Number of iterations
n=dim(dataset$time)[1] # Number of individuals
# Where output will be saved:
z=rep(0,MCMCits);y=rep(0,n);x=matrix(0,n,max(dataset$time));a=matrix(0,MCMCits,range(dataset$time)[2]);b=matrix(0,MCMCits,878)
storage=list(sigma=z,phi=z,muB=z,nuB=z,muL=z,nuL=z,muM=z,nuM=z,muS=z,nuS=z,logp=z,logl=z,
             tinf=y,tinf2=y,inf=y,
             titres=x,titres2=x,ltitres=x,ltitres2=x,
             CAR=a,GMT=a,PSP=a,PSC=a, B=b,S=b,T=b,M=b,I=b,L=b)
rm(x,y,z,a,b)


if(verbose)cat('Starting burnin loop at',date(),'\n')
## Main MCMC loop
for(iteration in 1:BURNINs)
{
  ## Track progress
  if(iteration%%10==0)cat('BURN==> ',iteration,' ',sum(theta$I[dataset$cohort==1]),' ',round(100*mean(theta$I[dataset$cohort==1])),'% [LL = ',round(theta$logl),']\n',sep='')
  if(iteration%%200==0)oneiterationplots(theta,dataset,'latestplots.pdf')
  ## Update parameters
  theta = propose_hyper(1,theta)
  theta = propose_hyper(2,theta)
  theta = propose_hyper(3,theta)
  theta = propose_hyper(4,theta)
  for(individual in 1:n)theta = propose_indiv(individual,theta,dataset)
  theta = propose_sigma(theta,dataset)
  theta = propose_phi(theta,dataset)
}

if(verbose)cat('Starting MCMC loop at',date(),'\n')
for(iteration in 1:MCMCits)
{
  ## Track progress
  if(iteration%%10==0)cat('MCMC==> ',iteration,' ',sum(theta$I[dataset$cohort==1]),' ',round(100*mean(theta$I[dataset$cohort==1])),'% [LL = ',round(theta$logl),']\n',sep='')
  if(iteration%%200==0)oneiterationplots(theta,dataset,'latestplots.pdf')
  
  ## Update parameters
  theta = propose_hyper(1,theta)
  theta = propose_hyper(2,theta)
  theta = propose_hyper(3,theta)
  theta = propose_hyper(4,theta)
  for(individual in 1:n)theta = propose_indiv(individual,theta,dataset)
  theta = propose_sigma(theta,dataset)
  theta = propose_phi(theta,dataset)
  
  # Store stuff
  storage = store(theta,storage,iteration,MCMCits,dataset)
}

if(verbose)cat('Finished MCMC loop at',date(),'\n')

plot(storage$logl+storage$logp)

dumptofile(theta,storage)

#rm(MCMCits,n,storage,iteration)


