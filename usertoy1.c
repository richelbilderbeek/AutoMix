/* User functions for the toy example in section 5.5.1 of Ph.D. thesis. */

#include<stdio.h>
#include<math.h>

#define tpi 6.283185307179586477

extern double sdrand();

/* Internal functions used in required user functions */
double boxm();

/* Function to return number of models */
void getkmax(int *kmax){
  *kmax=2;
  return;
}

/* Function to return the dimension of each model */ 
void getnk(int kmax,int *nk){
  int k;
  for(k=0;k<kmax;k++){
    nk[k]=k+1;
  }
  return;
}

/* Function to return initial conditions for RWM runs */
void getic(int k, int nkk, double *theta){
  int j;
  for(j=0;j<nkk;j++){
    theta[j]=boxm();
  }
}

/* Function to return log target distribution up to additive const at (k,theta)
   value also returned in llh */

double lpost(int k,int nkk,double *theta, double *llh){
  int j1,j2,l1;
  double work[nkk],mu[nkk+1][nkk];
  double sig[nkk+1][nkk][nkk],detsig[nkk+1],w[nkk+1];
  double lp,lptemp;

  if(k==0){
    mu[0][0]=-3;
    mu[1][0]=2;
    sig[0][0][0]=2.0;
    sig[1][0][0]=1.0;
    w[0]=0.2;
    w[1]=0.8;
  }
  if(k==1){
    mu[0][0]=0;
    mu[0][1]=3;
    mu[1][0]=-4;
    mu[1][1]=1;
    mu[2][0]=4;
    mu[2][1]=1;
    sig[0][0][0]=2;
    sig[0][1][0]=0;
    sig[0][1][1]=0.7071068;
    sig[1][0][0]=1.414214;
    sig[1][1][0]=1.060660;
    sig[1][1][1]=0.9354143;
    sig[2][0][0]=1.414214;
    sig[2][1][0]=-1.060660;
    sig[2][1][1]=0.9354143;
    w[0]=1.0/3.0;
    w[1]=1.0/3.0;
    w[2]=1.0/3.0;
  }

  for(l1=0;l1<(nkk+1);l1++){
    detsig[l1]=1.0;
    for(j1=0;j1<nkk;j1++){
      (detsig[l1])*=sig[l1][j1][j1];
    }
  }

  lp=0.0;
  for(l1=0;l1<(nkk+1);l1++){  
    for(j1=0;j1<nkk;j1++){
      work[j1]=theta[j1]-mu[l1][j1];
    }
    for(j1=0;j1<nkk;j1++){
      for(j2=0;j2<j1;j2++){
        (work[j1])-=sig[l1][j1][j2]*work[j2];
      }
      (work[j1])/=sig[l1][j1][j1]; 
    }
    lptemp=0.0;
    for(j1=0;j1<nkk;j1++){
      lptemp+=(work[j1]*work[j1]);
    }
    lptemp=w[l1]*exp(-0.5*lptemp)*pow(tpi,-nkk/2.0)/(detsig[l1]);    
    lp+=lptemp;
  }
  if(k==0){
    lp=log(0.3*lp);
  }
  else{
    lp=log(0.7*lp);
  }

  *llh=lp;
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
