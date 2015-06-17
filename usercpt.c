/* User functions for the change point example in section 5.5.2 of Ph.D. 
   thesis.
   Data is for times of coal mining disasters and is taken from 
   Green, 1995 (see thesis for full reference). */

#include<stdio.h>
#include<math.h>

#define tpi 6.283185307179586477

/* Function loggamma is in file gammafns.c bundled with the AutoMix software */

extern double loggamma(double x);

/* Hyperparameters */

double alpha=1.0,beta=200.0,lambda=3.0;

/* Function to return number of models */
void getkmax(int *kmax){
  *kmax=6;
  return;
}

/* Function to return the dimension of each model */ 
void getnk(int kmax,int *nk){
  int k;
  for(k=0;k<kmax;k++){
    /* model k has k+1 changepoints and k+2 rates */
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
    rwm[k+1+j]=(40907.0*j)/(k+2);
  }
}

/* Function to return log of posterior up to additive const at (k,theta)
   likelihood returned in llh1 */

double lpost(int k,int nkk,double *theta,double *llh1){

  int i,j,nsteps=(k+1),nsofar,nj;
  double lp,llh,logl,abcon,top;
  double h[nsteps+1],s[nsteps+2],ds[nsteps+1];

  /* Data kindly provided by PJG */

  double y[191]={74,231,354,356,480,492,496,506,722,802,814,847,
		 913,1145,1971,2011,2023,2052,2242,2339,2404,2590,
		 2613,2705,2902,3333,3349,3503,3598,3623,3642,3720,
		 3922,3958,4068,4344,4360,4448,4673,4726,4743,5281,
		 5468,5502,5603,5644,5783,5825,5826,6076,6156,6159,
		 6483,6539,6570,6666,6736,6777,6870,6894,6985,7128,
		 7144,7171,7315,7360,7366,7574,7603,7715,7758,7951,
		 8085,8505,8600,8725,8759,8886,9104,9106,9106,9484,
		 9520,9535,9566,9781,9792,9929,9933,9948,10020,
		 10116,10240,10290,10410,10613,10789,10844,10937,
		 10996,11311,11370,11431,11432,11445,11634,11979,
		 11999,12080,12366,12480,12588,12776,13009,13037,
		 13059,13120,13198,13297,13623,13898,13952,14169,
		 14282,14314,14702,14853,15214,15526,15880,16187,
		 16462,16540,16557,17762,18406,18873,19744,19792,
		 19915,20371,20869,20918,21049,21231,21486,21680,
		 21904,22470,22932,23160,23966,24483,26126,26180,
		 26506,27818,28166,28911,29128,29248,29523,29543,
		 29609,29901,29905,30273,30580,30916,30935,31264,
		 31594,31906,32442,32587,32662,33026,33063,33082,
		 33238,33285,33414,35044,35073,35290,35297,35315,
		 36673,39039,39991,40623};

  
  logl=log(40907.0);
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
  ds[nsteps]=40907.0-s[nsteps];
  s[nsteps+1]=40907.0;

  for(i=0;i<=nsteps;i++){
    if((h[i]<=0.0)||ds[i]<=0.0){
      lp=-10000.0;
      return lp;
    }
  }

  /* Prior (settings detailed in Green, 1995) */ 

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
