/* User functions for the Rb9 example in section 5.5.3 of Ph.D. thesis.
   Data is for tumour counts is taken from 
   Haigis and Dove, 2003 (see thesis for full reference). */

#include<stdio.h>
#include<math.h>

#define tpi 6.283185307179586477

extern double loggamma(double x);

extern double sdrand();

/* Hyperparameters */
double alpha1=2.0,alpha2=1.0,beta1=0.1,beta2=2.0;

/* Internal functions used in required user functions */
double boxm();
 
/* Function to return number of models */
void getkmax(int *kmax){
  *kmax=10;
  return;
}

/* Function to return the dimension of each model */ 
void getnk(int kmax,int *nk){
  int k;
  for(k=0;k<6;k++){
    nk[k]=4;
  }
  for(k=6;k<kmax;k++){
    nk[k]=5;
  }
  return;
}

/* Function to return initial conditions for RWM runs */
void getic(int k, int nkk, double *rwm){
  int j,ql,qk;
  int nlambda[10],nkappa[10];
  double u;

  for(j=0;j<10;j++){
    nlambda[j]=3;
    nkappa[j]=1;
  }
  nlambda[7]=4;
  nlambda[8]=4;
  nlambda[9]=4;
  nkappa[6]=2;

  ql=nlambda[k];
  qk=nkappa[k];
    
  for(j=0;j<ql;j++){
    u=exp(boxm());
    rwm[j]=43.87879*u;
  }
  for(j=ql;j<ql+qk;j++){
    u=exp(boxm());
    rwm[j]=2.152937*u;
  }
}

/* Function to return log of posterior up to additive const at (k,theta)
   likelihood returned in llh1 - prior settings as in Thesis chapter 3 */

double lpost(int k,int nkk,double *theta, double *llh1){

  int i,j,ql,qk,cumnobs;
  int nobs[4],nlambda[10],nkappa[10],pindic[4];
  double lp,llh,lambda[4],kappa[4],dlgdummy1,dlgdummy2,kappamin1;
  int X[66]={121,169,112,199,80,121,194,140,131,199,262,
		121,140,166,150,103,5,15,13,9,15,13,13,9,18,
		12,8,7,16,11,12,8,14,12,20,12,8,11,10,10,10,
		7,8,7,8,10,11,7,4,6,9,7,5,7,3,7,4,11,15,10,
		6,10,6,12,6,11};

  for(i=0;i<nkk;i++){
    if(theta[i]<0){
      lp=-1000000.0;
      return lp;
    }
  }
  nobs[0]=16;
  nobs[1]=17;
  nobs[2]=15;
  nobs[3]=18;

  for(j=0;j<10;j++){
    nlambda[j]=3;
    nkappa[j]=1;
  }
  nlambda[7]=4;
  nlambda[8]=4;
  nlambda[9]=4;
  nkappa[6]=2;

  ql=nlambda[k];
  qk=nkappa[k];

  pindic[0]=1;
  pindic[1]=0;
  pindic[2]=0;
  pindic[3]=1;
  if(k==3||k==9){
    pindic[1]=1;
  }
  if(k==2||k==9){
    pindic[2]=1;
  }
  if(k==0||k==4||k==7){
    pindic[3]=0;
  }

  lambda[0]=theta[0];
  lambda[1]=theta[1];
  if(k<4||k==6){
    lambda[2]=theta[1];
  }
  else{
    lambda[2]=theta[2];
  }
  if(k<7){
    lambda[3]=theta[2];
  }
  else{
    lambda[3]=theta[3];
  }
  if(k<7){
    kappa[0]=theta[3];
  }
  else{
    kappa[0]=theta[4];
  }

  /* know next 3 lines aren't always true but for cases which 
     it is wrong don't use the incorrect values of kappa */
  kappa[1]=kappa[0];
  kappa[2]=kappa[0];
  kappa[3]=kappa[0];
  if(k==6){
    kappa[3]=theta[4];
  }

  /* prior */
  lp=0.0;
  for(i=0;i<ql;i++){
    lp+=alpha1*log(beta1)+(alpha1-1.0)*log(theta[i])
      -beta1*theta[i]-loggamma(alpha1);
  }
  for(i=ql;i<ql+qk;i++){
    lp+=alpha2*log(beta2)+(alpha2-1.0)*log(theta[i])-
      beta2*theta[i]-loggamma(alpha2);
  }

  /* likelihood */
  llh=0.0;
  cumnobs=0;
  
  for(i=0;i<4;i++){
    if(pindic[i]==0){
      /*Poisson*/ 
      for(j=0;j<nobs[i];j++){
	dlgdummy1=X[cumnobs]+1;
	llh+=(-lambda[i]+X[cumnobs]*log(lambda[i])-loggamma(dlgdummy1));
        cumnobs++;
      }
    }
    else{
      /* Negative Binomial */
      kappamin1=1.0/kappa[i];
      for(j=0;j<nobs[i];j++){
	dlgdummy1=X[cumnobs]+1;
	dlgdummy2=X[cumnobs]+kappamin1;
	llh+=(X[cumnobs]*log(lambda[i])+loggamma(dlgdummy2)-
	      loggamma(dlgdummy1)+kappamin1*log(kappamin1)-
	      loggamma(kappamin1)-(X[cumnobs]+kappamin1)*
	      log(lambda[i]+kappamin1));
	cumnobs++;
      }
    }   
  }

  *llh1=llh;
  lp+=llh;
  return lp;
      
} 

double boxm(){

  /*Function for returning single N(0,1) variable */

  double u1,u2,out;

  u1=sdrand();
  u2=sdrand();
  u1=tpi*u1;
  u2=sqrt(-2.0*log(u2));
  out=u2*cos(u1);

  return out;

}
