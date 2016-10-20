// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;




// [[Rcpp::export]]
double logPrior(double sigma, double phi, vec I,vec B, vec L, vec M, vec S, double muB, double nuB, double muL, double nuL, double muM, double nuM, double muS, double nuS)
{
  double output = 0.0;
  int n = B.n_elem;
//  int nt = titres.n_cols;

  // Log prior from sigma and phi: constant so can ignore
  output += 0.0;

  // Log prior from BLMS parameters, univariate normals:
  for(int i=1; i<n; i++)
  {
    output += R::dnorm(B(i),muB,sqrt(nuB),1);
    
   if(I(i)==1){
      output += R::dnorm(L(i),muL,sqrt(nuL),1);
      output += R::dnorm(M(i),muM,sqrt(nuM),1);
      output += R::dnorm(S(i),muS,sqrt(nuS),1);
    }
    
  }
  
  // Log prior from hyper parameters: make constant so can ignore
  output += 0.0;

  return output;
}

// [[Rcpp::export]]
double titre(double diff, double eB, double eL, double eM, double eS, int Ii)
{
  //titre(diff,eB,eL,eM,eS,Ii)
  double mode = eM+2.0;
  double var  = eS*eS;
  double scale = (-1.0*mode +sqrt(mode*mode+4.0*var))/2.0;
  double shape = var/(scale*scale);
  if(scale<0.0)return 99.9;
  if(shape<0.0)return 99.9;
  double output = eB;

  if(Ii==1)
  {
    if(diff>0.0)
    {
      output+=R::dgamma(diff,shape,scale,0)*eL;
    }
  }
  return output;
}

// [[Rcpp::export]]
mat titers(vec B, vec L, vec M, vec S, vec I, vec T, int nt)
{ 
  int n = B.n_elem;
  mat output(n,nt);

  for(int i=0;i<n;++i)
  {
    double eB=exp(B(i));
    double eL=exp(L(i));
    double eM=exp(M(i));
    double eS=exp(S(i));    
    for(int j=0;j<nt;++j)
    {
      double diff = j+1-T(i);
      double mean=titre(diff,eB,eL,eM,eS,I(i));
      output(i,j)=mean;
    }
  }
  
  return(output);
}


// [[Rcpp::export]]
double logLikelihood(double sigma, double phi, vec B, vec L, vec M, vec S, vec I, vec T, mat low, mat upp, mat time, vec phat, vec cohort)
{ 
  double output = 0.0;
  double MAXtitre=11.0;
  int n = B.n_elem;
  int nobs = low.n_cols;
  
  // Contribution from titres
  for(int i=0;i<n;++i)
  {
    double eB=exp(B(i));
    double eL=exp(L(i));
    double eM=exp(M(i));
    double eS=exp(S(i));
    double esigma = exp(sigma);
    for(int j=0;j<nobs;++j)
    {
      if(low(i,j)>-98) // placeholder for no observation
      {
        double diff = time(i,j)-T(i);
        double mean=titre(diff,eB,eL,eM,eS,I(i));
        double diff2=R::pnorm(upp(i,j),mean,esigma,1,0)-R::pnorm(low(i,j),mean,esigma,1,0);
        if(diff2<0.0000001)diff2=0.0000001;
        output += log(diff2);
      }
    }
    double maxtitre=titre(eM,eB,eL,eM,eS,I(i));
    if(maxtitre>MAXtitre)output -= 999999.9;
  }
  
  
  // Contribution from infection status and time
  int nt = phat.n_elem;
  vec h(nt);
  vec logh(nt);
  vec H(nt);
  double ephi=exp(phi);
  for(int j=0; j<nt; j++){h(j) = ephi;}
  H(0)= h(0);
  for(int j=1; j<nt; j++){H(j) = H(j-1)+h(j);}
  for(int j=0; j<nt; j++){logh(j) = log(h(j));}
 // double temp1 = 0.0;
 // double temp2 = 0.0;
  for(int i=0; i<n; i++)
  {
    if(cohort(i)==1)
    {
      if(I(i)==1)output+= logh(T(i)-1)-H(T(i)-1);
      if(I(i)==0)output+= -1.0*H(nt-1);
      //temp1+= logh(T(i)-1)-H(T(i)-1);
      //temp2+= -1.0*H(nt-1);
    }
  }
  //cout << temp1 << " / " << temp2 << std::endl;
  
  return(output);
}



// [[Rcpp::export]]
double logLikelihood1(int ind, double sigma, double phi, vec B, vec L, vec M, vec S, vec I, vec T, mat low, mat upp, mat time, vec phat, vec cohort)
{ 
  double MAXtitre=11.0;
  double output = 0.0;
//  int n = B.n_elem;
  int nobs = low.n_cols;
  
  // Contribution from titres
  int i = ind;
  double eB=exp(B(i));
  double eL=exp(L(i));
  double eM=exp(M(i));
  double eS=exp(S(i));
  double esigma = exp(sigma);
  for(int j=0;j<nobs;++j)
  {
    if(low(i,j)>-98) // placeholder for no observation
    {
      double diff = time(i,j)-T(i);
      double mean=titre(diff,eB,eL,eM,eS,I(i));
      double diff2=R::pnorm(upp(i,j),mean,esigma,1,0)-R::pnorm(low(i,j),mean,esigma,1,0);
      if(diff2<0.0000001)diff2=0.0000001;
      
      output += log(diff2);
      //if(mean>MAXtitre)output -= 999999.9;
    }
  }
  double maxtitre=titre(eM,eB,eL,eM,eS,I(i));
  if(maxtitre>MAXtitre)output -= 999999.9;
  
  // Contribution from infection status and time
  int nt = phat.n_elem;
  vec h(nt);
  vec logh(nt);
  vec H(nt);
  double ephi=exp(phi);
  for(int j=0; j<nt; j++){h(j) = ephi;}
  H(0)= h(0);
  for(int j=1; j<nt; j++){H(j) = H(j-1)+h(j);}
  for(int j=0; j<nt; j++){logh(j) = log(h(j));}

  if(cohort(i)==1)
  {
    if(I(i)==1)output+= logh(T(i)-1)-H(T(i)-1);
    if(I(i)==0)output+= -1.0*H(nt-1);
  }
  
  
  return(output);
}


/* o */
