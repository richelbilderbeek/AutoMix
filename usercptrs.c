/* User functions for the rescaled change point example in section 5.5.2 
   of Ph.D. thesis.
   Original data is for times of coal mining disasters and is taken from 
   Green, 1995 (see thesis for full reference). */

#include<stdio.h>
#include<math.h>

#define tpi 6.283185307179586477

/* Function loggamma is in file gammafns.c bundled with the AutoMix software */

extern double loggamma(double x);

/* Hyperparameters */
/* beta is now 200/1459 (this means mean is 1459 times orig and var is 1459^2
   times orig - in line with transformation of times*/
double alpha=1.0,beta=0.137,lambda=3.0;

/* Function to return number of models */
void getkmax(int *kmax){
  *kmax=6;
  return;
}

/* Function to return the dimension of each model */ 
void getnk(int kmax,int *nk){
  int k;
  for(k=0;k<kmax;k++){
    nk[k]=2*k+3;
  }
  return;
}

/* Function to return initial conditions for RWM runs */
void getic(int k, int nkk, double *rwm){
  int j;

  for(j=0;j<k+2;j++){
    rwm[j]=alpha/beta;
  }
  for(j=1;j<k+2;j++){
    rwm[k+1+j]=(28.04*j)/(k+2);
  }
}

/* Function to return log of posterior up to additive const at (k,theta)
   likelihood returned in llh1 */

double lpost(int k,int nkk,double *theta,double *llh1){

  int i,j,nsteps=(k+1),nsofar,nj;
  double lp,llh,logl,abcon,top;
  double h[nsteps+1],s[nsteps+2],ds[nsteps+1];

  /* Rescaled data */

  double y[191]={0.05,0.16,0.24,0.24,0.33,0.34,0.34,0.35,0.49,0.55,
		 0.56,0.58,0.63,0.78,1.35,1.38,1.39,1.41,1.54,1.6,1.65,
		 1.78,1.79,1.85,1.99,2.28,2.3,2.4,2.47,2.48,2.5,2.55,
		 2.69,2.71,2.79,2.98,2.99,3.05,3.2,3.24,3.25,3.62,3.75,
		 3.77,3.84,3.87,3.96,3.99,3.99,4.16,4.22,4.22,4.44,4.48,
		 4.5,4.57,4.62,4.64,4.71,4.73,4.79,4.89,4.9,4.92,5.01,
		 5.04,5.05,5.19,5.21,5.29,5.32,5.45,5.54,5.83,5.89,5.98,
		 6,6.09,6.24,6.24,6.24,6.5,6.53,6.54,6.56,6.7,6.71,6.81,
		 6.81,6.82,6.87,6.93,7.02,7.05,7.14,7.27,7.39,7.43,7.5,
		 7.54,7.75,7.79,7.83,7.84,7.84,7.97,8.21,8.22,8.28,8.48,
		 8.55,8.63,8.76,8.92,8.94,8.95,8.99,9.05,9.11,9.34,9.53,
		 9.56,9.71,9.79,9.81,10.08,10.18,10.43,10.64,10.88,11.09,
		 11.28,11.34,11.35,12.17,12.62,12.94,13.53,13.57,13.65,
		 13.96,14.3,14.34,14.43,14.55,14.73,14.86,15.01,15.4,
		 15.72,15.87,16.43,16.78,17.91,17.94,18.17,19.07,19.31,
		 19.82,19.96,20.05,20.24,20.25,20.29,20.49,20.5,20.75,
		 20.96,21.19,21.2,21.43,21.65,21.87,22.24,22.34,22.39,
		 22.64,22.66,22.67,22.78,22.81,22.9,24.02,24.04,24.19,
		 24.19,24.2,25.14,26.76,27.41,27.84};
  
  logl=log(28.04);
  abcon=alpha*log(beta)-loggamma(alpha);
  h[0]=theta[0];
  s[0]=0.0;
  /* work out heights (rates) h and steps (change points) s and dist between 
     steps ds */
  for(i=1;i<=nsteps;i++){
    h[i]=theta[i];
    s[i]=theta[nsteps+i];
    ds[i-1]=s[i]-s[i-1];
    
  }
  ds[nsteps]=28.04-s[nsteps];
  s[nsteps+1]=28.04;

  for(i=0;i<=nsteps;i++){
    if((h[i]<=0.0)||ds[i]<=0.0){
      lp=-100000.0;
      return lp;
    }
  }

  /* Prior (equivalent to settings detailed in Green, 1995) */ 
  lp=-lambda+nsteps*log(lambda)-loggamma((double)(nsteps+1));

  for(i=0;i<=nsteps;i++){
    lp+=(abcon+(alpha-1.0)*log(h[i])-beta*h[i]);
    lp+=log(ds[i]);
  }

  lp+=(loggamma(2.0*(nsteps+1))-(2.0*nsteps+1.0)*logl);

  nsofar=0;
  top=s[1];
  j=0;

  /* Likelihood */
  llh=0.0; 
  for(i=0;i<191;i++){
    if(y[i]>top){
      nj=i-nsofar;
      nsofar=i;
      llh+=(nj*log(h[j])-h[j]*ds[j]);
      j++;
      if(j>nsteps){
	return lp;
      }
      top=s[j+1];
    }
  }
  nj=191-nsofar;
  llh+=nj*log(h[j])-h[j]*ds[j];
  lp+=llh;
  *llh1=llh;
  return lp;
      
} 
