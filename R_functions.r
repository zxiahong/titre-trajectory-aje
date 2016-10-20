import_data=function()
{
  ## Read data: MOH data
  source('mesoput/mohphat.r')
  source('mesoput/MATtime.r')
  source('mesoput/MATlow.r')
  source('mesoput/MATupp.r')
  source('mesoput/cohort.r')
  
  ## Create dataset for inference
  dataset=list(time=MATtime,low=MATlow,upp=MATupp,cohort=cohort,phat=mohphat)
  return(dataset)
}

whichfile=function(a,b){if(file.exists(a))return(a);return(b)}

import_inits=function(guided=TRUE)
{
  if(guided)source('mesoput/smart_inits.r')
  if(!guided)source(whichfile('mesoput/estimated_inits.r','mesoput/smart_inits.r'))
  theta=inits
  theta$T=theta$T+(theta$T==0) # ie make 0s --> 1s
  theta$logp=theta$logl=0 # placeholder
  return(theta)
}



logf=function(indi=1,theta,data) ## log likelihood contribution
{
  output = logLikelihood1(indi-1,theta$sigma,theta$phi, theta$B, theta$L, theta$M, theta$S, theta$I, theta$T, data$low, data$upp, data$time, data$phat, data$cohort)
  #output = logLikelihood(theta$sigma,theta$phi, theta$B, theta$L, theta$M, theta$S, theta$I, theta$T, data$low, data$upp, data$time, data$phat, data$cohort)
  return(output)
}

logprior=function(theta)
{
  output = logPrior(theta$sigma,theta$phi, theta$I,theta$B, theta$L, theta$M, theta$S, theta$muB, theta$nuB, theta$muL, theta$nuL, theta$muM, theta$nuM, theta$muS, theta$nuS)
  return(output)
}

loglikelihooder=function(theta,data)
{
  output = logLikelihood(theta$sigma,theta$phi, theta$B, theta$L, theta$M, theta$S, theta$I, theta$T, data$low, data$upp, data$time, data$phat, data$cohort)
  return(output)
}


calculate_titre=function(diffs, eB, eL, eM, eS, Ii) # diffs can be a vector
{
  output = length(diffs)
  for(di in 1:length(diffs))
    output[di] = titre(diffs[di], eB, eL, eM, eS, Ii)
  return(output)
}

propose_sigma=function(theta,data)
{
  reject=FALSE
  theta$logl = loglikelihooder(theta,data)
  theta$logp = logprior(theta)
  old=theta
  theta$sigma = rnorm(1,theta$sigma,0.05)
  if(!reject)
  {
    theta$logl = loglikelihooder(theta,data)
    theta$logp = logprior(theta)
    logap = theta$logl + theta$logp - old$logl - old$logp
    if(-rexp(1)>logap) reject = TRUE
  }
  if(reject)return(old)
  return(theta)
}

propose_phi=function(theta,data)
{
  reject=FALSE
  theta$logl = loglikelihooder(theta,data)
  theta$logp = logprior(theta)
  old=theta
  theta$phi = rnorm(1,theta$phi,0.1)
  if(!reject)
  {
    theta$logl = loglikelihooder(theta,data)
    theta$logp = logprior(theta)
    logap = theta$logl + theta$logp - old$logl - old$logp
    if(-rexp(1)>logap) reject = TRUE
  }
  if(reject)return(old)
  return(theta)
}

gibb=function(a,nu)
{
  n=length(a)
  tausq=1/nu
  mu0=0;kappa0=1;
  v0=1;tau0=1
  abar=mean(a)
  
  mu=rnorm(1,(mu0*kappa0^2+abar*n*tausq)/(kappa0^2+n*tausq),1/sqrt(kappa0^2+n*tausq))
  tausq=rgamma(1,(v0+n)*0.5,v0/(2*tau0^2)+0.5*(n-1)*var(a)+n*0.5*(abar-mu)^2)
  nu=1/tausq
  return(c(mu,nu))
}



propose_hyper=function(parid,theta)
{
  if(parid==1)
  {
    temp=gibb(theta$B,theta$nuB)
    theta$muB = temp[1]
    theta$nuB = temp[2]
  }
  if(parid==2)
  {
    temp=gibb(theta$L,theta$nuL)
    theta$muL = temp[1]
    theta$nuL = temp[2]
  }
  if(parid==3)
  {
    temp=gibb(theta$M,theta$nuM)
    theta$muM = temp[1]
    theta$nuM = temp[2]
  }
  if(parid==4)
  {
    temp=gibb(theta$S,theta$nuS)
    theta$muS = temp[1]
    theta$nuS = temp[2]
  }
  return(theta)
}

propose_indiv=function(ind,theta,data)
{
  #if(!data$cohort[ind]) ##### TODO
  #{
  reject=FALSE
  old=theta
  
  
  
  if(data$cohort[ind]==1)
  {
    theta$I[ind] = rbinom(1,1,0.5)
    
  }
   # muBi=mean(storage$B[,ind]); #if(iteration==1) muBi=theta$B[ind]
   # muLi=mean(storage$L[,ind]); #if(iteration==1) muLi=theta$L[ind]
   # muMi=mean(storage$M[,ind]); #if(iteration==1) muMi=theta$M[ind]
   # muSi=mean(storage$S[,ind]); #if(iteration==1) muSi=theta$S[ind]
  
  muBi=theta$B[ind]
  muLi=theta$L[ind]
  muMi=theta$M[ind]
  muSi=theta$S[ind]
  #
  
  muT=theta$T[ind]
  theta$B[ind] = rnorm(1,muBi,0.01)
  
  if(theta$I[ind]==1){           
    theta$L[ind] = rnorm(1,muLi,0.1)
    theta$M[ind] = rnorm(1,muMi,0.1)
    theta$S[ind] = rnorm(1,muSi,0.1)
    
    
    if(data$cohort[ind]==1){
      theta$T[ind] = round(rnorm(1,muT,1))
      if(theta$T[ind]<1)reject=TRUE
      if(theta$T[ind]>max(data$time))reject=TRUE
    }
    
    
  }
  
  
  if(!reject)
  {
    if((old$I[ind]-theta$I[ind])==0){
      logap=logf(ind,theta,data)+ logprior(theta) -logf(ind,old,data)-logprior(old)      
    }
    
    if(((old$I[ind]-theta$I[ind])==-1) && (data$cohort[ind]==1)){  ## birth
      logap=logf(ind,theta,data)+ logprior(theta) -logf(ind,old,data)-logprior(old)-log(dnorm(theta$L[ind],muLi,0.1))-log(dnorm(theta$M[ind],muMi,0.1))-log(dnorm(theta$S[ind],muSi,0.1))-log(dnorm(theta$T[ind],muT,1))
    }
    if(((old$I[ind]-theta$I[ind])==1) && (data$cohort[ind]==1)){  ## death
      logap=logf(ind,theta,data)+ logprior(theta) -logf(ind,old,data)-logprior(old)+log(dnorm(theta$L[ind],muLi,0.1))+log(dnorm(theta$M[ind],muMi,0.1))+log(dnorm(theta$S[ind],muSi,0.1))+log(dnorm(theta$T[ind],muT,1))
    }
    
    if(((old$I[ind]-theta$I[ind])==-1) && (data$cohort[ind]==0)){  ## birth
      logap=logf(ind,theta,data)+ logprior(theta) -logf(ind,old,data)-logprior(old)-log(dnorm(theta$L[ind],muLi,0.1))-log(dnorm(theta$M[ind],muMi,0.1))-log(dnorm(theta$S[ind],muSi,0.1))
    }
    if(((old$I[ind]-theta$I[ind])==1) && (data$cohort[ind]==0)){  ## death
      logap=logf(ind,theta,data)+ logprior(theta) -logf(ind,old,data)-logprior(old)+log(dnorm(theta$L[ind],muLi,0.1))+log(dnorm(theta$M[ind],muMi,0.1))+log(dnorm(theta$S[ind],muSi,0.1))
    }
    
    logu = -rexp(1)
    if(logu>logap) reject = TRUE
  }
  #print(reject)
  if(reject)return(old)
  
  #}
  return(theta)
}


store=function(theta,storage,iteration,MCMCits,data)
{
  output=storage
  output$sigma[iteration] = theta$sigma
  output$phi[iteration] = theta$phi
  output$muB[iteration]=theta$muB
  output$nuB[iteration]=theta$nuB
  output$muL[iteration]=theta$muL
  output$nuL[iteration]=theta$nuL
  output$muM[iteration]=theta$muM
  output$nuM[iteration]=theta$nuM
  output$muS[iteration]=theta$muS
  output$nuS[iteration]=theta$nuS
  output$logp[iteration]=theta$logp
  output$logl[iteration]=theta$logl
  output$B[iteration,]=theta$B
  output$S[iteration,]=theta$S
  output$M[iteration,]=theta$M
  output$L[iteration,]=theta$L
  output$I[iteration,]=theta$I
  output$T[iteration,]=theta$T

  output$tinf=storage$tinf+theta$I*theta$T/MCMCits
  output$tinf2=storage$tinf2+theta$I*theta$T^2/MCMCits
  output$inf=storage$inf+theta$I/MCMCits
  temp=titerwrapper(theta,data)
  ltemp = log(temp)
  output$titres=storage$titres+temp/MCMCits
  output$titres2=storage$titres2+temp^2/MCMCits
  output$ltitres=storage$ltitres+ltemp/MCMCits
  output$ltitres2=storage$ltitres2+ltemp^2/MCMCits
  
  co=which(data$cohort==1)
  output$GMT[iteration,]=colMeans(temp[co,])
  output$PSP[iteration,]=colMeans(temp[co,]>=3)
  temp2=temp[,1]
  delta=temp-temp2
  vari=2*exp(theta$sigma)^2
  epsilon=matrix(rnorm(length(delta),0,sqrt(vari)),dim(delta)[1],dim(delta)[2])
  diff=delta+epsilon
  output$PSC[iteration,]=colMeans(diff[co,]>2)
  
  temp=theta$T*theta$I + (1-theta$I)*(max(data$time)+1)
  temp2=tabulate(temp[co])
  output$CAR[iteration,] = cumsum(temp2[-length(temp2)])/sum(temp2)
  return(output)
}

dumptofile=function(theta,storage)
{
  inits = theta
  dump('inits','estimated_inits.r')
  
  temp=data.frame(sigma=storage$sigma,phi=storage$phi,B=storage$B,M=storage$M,S=storage$S,L=storage$L, T=storage$T,I=storage$I,
                  muB=storage$muB,nuB=storage$nuB,muL=storage$muL,nuL=storage$nuL,muM=storage$muM,
                  nuM=storage$nuM,muS=storage$muS,nuS=storage$nuS,logp=storage$logp,logl=storage$logl)
  write.table(temp,'CODA_main.csv',row.names=FALSE,col.names=TRUE,sep=',')
  
  temp=data.frame(tinf=storage$tinf,tinf2=storage$tinf2,inf=storage$inf)
  write.table(temp,'CODA_infection_status.csv',row.names=FALSE,col.names=TRUE,sep=',')
  
  write.table(storage$titres,'CODA_titres.csv',row.names=FALSE,col.names=TRUE,sep=',')
  write.table(storage$titres2,'CODA_titres2.csv',row.names=FALSE,col.names=TRUE,sep=',')
  write.table(storage$ltitres,'CODA_ltitres.csv',row.names=FALSE,col.names=TRUE,sep=',')
  write.table(storage$ltitres2,'CODA_ltitres2.csv',row.names=FALSE,col.names=TRUE,sep=',')
  
  write.table(storage$CAR,'CODA_CAR.csv',row.names=FALSE,col.names=TRUE,sep=',')
  write.table(storage$GMT,'CODA_GMT.csv',row.names=FALSE,col.names=TRUE,sep=',')
  write.table(storage$PSP,'CODA_PSP.csv',row.names=FALSE,col.names=TRUE,sep=',')
}

titerwrapper=function(theta,data)
{
  return(titers(theta$B, theta$L, theta$M, theta$S, theta$I, theta$T, max(data$time)))
}

plotfit=function(i,storage,data)
{
  temp=storage$ltitres[i,];temp2=storage$ltitres2[i,];temp3=sqrt(temp2-temp*temp)
  plot(exp(temp),type='l',ylim=c(0,10),xlab='Time',ylab='Titre');lines(exp(temp+1.96*temp3),col=8);lines(exp(temp-1.96*temp3),col=8)
  for(j in 1:7)
  {
    if(data$low[i,j]>-99)
    {
      lines(rep(data$time[i,j],2),data$low[i,j]+c(0,1),col=2)
    }
  }
}

oneiterationplots=function(theta,data,filename)
{
  nt=max(data$time)
  cm=1/2.54
  pdf(filename,height=10*cm,width=10*cm)
  for(i in 1:dim(data$low)[1])
  {
    eB=exp(theta$B[i]);
    eL=exp(theta$L[i]);
    eM=exp(theta$M[i]);
    eS=exp(theta$S[i]);
    kappai = eM*eM/(eS*eS); # shape
    thetai = eM/(eS*eS);    # rate
    mean=rep(eB,nt);
    if(theta$I[i]==1)
    {
      diff = (1:nt)-theta$T[i]+0.00001;
      mean=calculate_titre(diff, eB, eL, eM, eS, 1)
      #mean=mean+dgamma(diff+0.00001,kappai,thetai)*eL*(diff>0)
    }
    par(mai=c(2,2,0.5,0.5)*cm)
    plot(mean,type='l',ylim=c(0,10),xlab='Time',ylab='Titre')
    text(36,9,paste('individual',i))
    for(j in 1:7)
    {
      if(data$low[i,j]>-99)
      {
        lines(rep(data$time[i,j],2),data$low[i,j]+c(0,1),col=2)
      }
    }
  }
  dev.off()
}
