/* User functions for the AIDS clinical trial example in section 5.5.4 of 
   Ph.D. thesis.
   Example is taken from Han & Carlin, 2001 (see thesis for full reference). */

#include<stdio.h>
#include<math.h>

/* Data provided in user library - bundled with AutoMix software */
#include "ddidata.h"

#define tpi 6.283185307179586477
#define pi 3.141592653589793238
#define max(A,B) ((A)>(B)?(A):(B))

extern double loggamma(double x);

extern double sdrand();

/* hyperparameters */
double a=3.0,b=0.005;
int rho=24;
double c0[9]={10.0,0.0,0.0,0.0,0.0,0.0,-3.0,0.0,0.0};
double c1[6]={10.0,0.0,0.0,0.0,-3.0,0.0};
double D0min1[9]={0.25,1.0,1.0,100.0,1.0,1.0,1.0,1.0,1.0};
double D1min1[6]={0.25,1.0,100.0,1.0,1.0,1.0};
double R0[3][3]={{4.0,0.0,0.0},{0.0,1/16.0,0.0},{0.0,0.0,1.0/16.0}};
double R1[2][2]={{4.0,0.0},{0.0,1/16.0}};

/* Internal functions used in required user functions */
double boxm2();

void chol2(int nkk, double **B,int *posdef);

double det2(int nkk,double **B);

double lnormprob2(int nkk, double **B, double *mu, double *datai);

/* Function to return number of models */
void getkmax(int *kmax){
  *kmax=2;
  return;
}

/* Function to return the dimension of each model */ 
void getnk(int kmax,int *nk){
  
  nk[0]=16;
  nk[1]=10;
  return;

}

/* Function to return initial conditions for RWM runs */

void getic(int k, int nkk, double *rwm){
  int i,j,j1,j2,posdef;
  double u,alpha[9],gamma[6],Vmin1[3][3],Umin1[2][2],V[3][3],U[2][2];
  double sigmasq,tausq;
  double **Ctest,**temptest;

  /* Note boxm2() is used to return random initial conditions */

 BEGIN_IC:
  if(k==0){

    for(j=0;j<9;j++){
      alpha[j]=c0[j]+sqrt(1.0/D0min1[j])*boxm2();
    }
    for(j=0;j<3;j++){
      for(j1=0;j1<3;j1++){
	Vmin1[j][j1]=0;
        V[j][j1]=0;
      }
      u=exp(0.01*boxm2());
      Vmin1[j][j]=u;
      V[j][j]=1.0/u;
    }
    sigmasq=100.0*exp(0.1*boxm2());
    for(j=0;j<9;j++){
      rwm[j]=alpha[j];
    }
    rwm[9]=Vmin1[0][0];
    rwm[10]=Vmin1[1][0];
    rwm[11]=Vmin1[1][1];
    rwm[12]=Vmin1[2][0];
    rwm[13]=Vmin1[2][1];
    rwm[14]=Vmin1[2][2];
    rwm[15]=sigmasq;

    temptest=(double**)malloc(3*sizeof(double));
    Ctest=(double**)malloc(5*sizeof(double));
    for(j=0;j<5;j++){
      Ctest[j]=(double*)malloc(5*sizeof(double));
      if(j<3){
	temptest[j]=(double*)malloc(5*sizeof(double));
      }
    }

    for(i=0;i<467;i++){
      for(j=0;j<3;j++){
        for(j1=0;j1<5;j1++){
	  temptest[j][j1]=0.0;
	}
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<3;j2++){
	    temptest[j][j1]+=V[j][j2]*W[i][j1][j2];
	  }
	}
      }
      for(j=0;j<5;j++){
	for(j1=0;j1<5;j1++){
	  Ctest[j][j1]=0.0;
	}
      }
      for(j=0;j<S[i];j++){
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<3;j2++){
	    Ctest[j][j1]+=W[i][j][j2]*temptest[j2][j1];
	  }
	}
      }
      for(j=0;j<S[i];j++){
	Ctest[j][j]+=sigmasq;
      }
      posdef=1;
      chol2(S[i],Ctest,&posdef);
      if(posdef==0){
	goto BEGIN_IC;
      }    
    }
  }
  else{

    for(j=0;j<6;j++){
      gamma[j]=c1[j]+sqrt(1.0/D1min1[j])*boxm2();
    }
    for(j=0;j<2;j++){
      for(j1=0;j1<2;j1++){
	Umin1[j][j1]=0;
        U[j][j1]=0;
      }
      u=exp(0.01*boxm2());
      Umin1[j][j]=u;
      U[j][j]=1.0/u;
    }
    tausq=100.0*exp(0.1*boxm2());
    for(j=0;j<6;j++){
      rwm[j]=gamma[j];
    }
    rwm[6]=Umin1[0][0];
    rwm[7]=Umin1[1][0];
    rwm[8]=Umin1[1][1];
    rwm[9]=tausq;

    temptest=(double**)malloc(2*sizeof(double));
    Ctest=(double**)malloc(5*sizeof(double));
    for(j=0;j<5;j++){
      Ctest[j]=(double*)malloc(5*sizeof(double));
      if(j<2){
	temptest[j]=(double*)malloc(5*sizeof(double));
      }
    }

    for(i=0;i<467;i++){
      for(j=0;j<2;j++){
        for(j1=0;j1<5;j1++){
	  temptest[j][j1]=0.0;
	}
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<2;j2++){
	    temptest[j][j1]+=U[j][j2]*Q[i][j1][j2];
	  }
	}
      }
      for(j=0;j<5;j++){
	for(j1=0;j1<5;j1++){
	  Ctest[j][j1]=0.0;
	}
      }
      for(j=0;j<S[i];j++){
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<2;j2++){
	    Ctest[j][j1]+=Q[i][j][j2]*temptest[j2][j1];
	  }
	}
      }
      for(j=0;j<S[i];j++){
	Ctest[j][j]+=tausq;
      }
      posdef=1;
      chol2(S[i],Ctest,&posdef);
      if(posdef==0){
	goto BEGIN_IC;
      }
    }
  }
}

/* Function to return log of posterior up to additive const at (k,theta)
   likelihood returned in llh1 - prior settings as in Han and Carlin, 2001*/

double lpost(int k,int nkk,double *theta,double *llh1){

  int i,j,j1,j2,posdef;
  double lp,llh,alpha[9],gamma[6],**V,**U,sigmasq,tausq,**Vmin1,**Umin1;
  double detVmin1,detUmin1;
  double temp[3][5],**C,detC,mu[5],datai[5],**tempmat,**tempmat2;
  double **R0min1,**R1min1;
  double Y[467][5];

  C=(double**)malloc(5*sizeof(double));
  if(k==0){
    V=(double**)malloc(3*sizeof(double));
    Vmin1=(double**)malloc(3*sizeof(double));
    R0min1=(double**)malloc(3*sizeof(double));
    tempmat=(double**)malloc(3*sizeof(double));
    tempmat2=(double**)malloc(3*sizeof(double));
    for(j=0;j<5;j++){
      C[j]=(double*)malloc(5*sizeof(double));
      if(j<3){
	V[j]=(double*)malloc(3*sizeof(double));
	Vmin1[j]=(double*)malloc(3*sizeof(double));
	tempmat[j]=(double*)malloc(3*sizeof(double)); 
	tempmat2[j]=(double*)malloc(3*sizeof(double)); 
	R0min1[j]=(double*)malloc(3*sizeof(double));
      }
    }
    for(j=0;j<3;j++){
      for(j1=0;j1<3;j1++){
	R0min1[j][j1]=0.0;
      }
    }
    R0min1[0][0]=0.25;
    R0min1[1][1]=16.0;
    R0min1[2][2]=16.0;
  }
  else{
    U=(double**)malloc(2*sizeof(double));
    Umin1=(double**)malloc(2*sizeof(double));
    R1min1=(double**)malloc(2*sizeof(double));
    tempmat=(double**)malloc(2*sizeof(double));
    tempmat2=(double**)malloc(2*sizeof(double));
    for(j=0;j<5;j++){
      C[j]=(double*)malloc(5*sizeof(double));
      if(j<2){
	U[j]=(double*)malloc(2*sizeof(double));
	Umin1[j]=(double*)malloc(2*sizeof(double));
	tempmat[j]=(double*)malloc(2*sizeof(double)); 
	tempmat2[j]=(double*)malloc(2*sizeof(double)); 
	R1min1[j]=(double*)malloc(2*sizeof(double));
      } 
    }
    for(j=0;j<2;j++){
      for(j1=0;j1<2;j1++){
	R1min1[j][j1]=0.0;
      }
    }
    R1min1[0][0]=0.25;
    R1min1[1][1]=16.0;
  }

  for(i=0;i<467;i++){
    for(j=0;j<5;j++){
      Y[i][j]=sqrt(counts[i][j]);
    }
  }

  if(k==0){
    for(j1=0;j1<9;j1++){
      alpha[j1]=theta[j1];
    }
    for(j1=0;j1<3;j1++){
      for(j2=0;j2<=j1;j2++){
	Vmin1[j1][j2]=theta[9+(j1*j1+j1)/2+j2];
        Vmin1[j2][j1]=Vmin1[j1][j2];
      }
    }
    sigmasq=theta[15];
    if(sigmasq<0){
      for(j=0;j<5;j++){
	free(C[j]);
	if(j<3){
	  free(R0min1[j]);
	  free(tempmat[j]);
	  free(tempmat2[j]);
	  free(V[j]);
	  free(Vmin1[j]);
	}
      }
      free(C);
      free(R0min1);
      free(tempmat);
      free(tempmat2);
      free(V);
      free(Vmin1);
      lp=-10000000.0;
      *llh1=lp;
      return lp;
    }
     
  }
  else if(k==1){
    for(j1=0;j1<6;j1++){
      gamma[j1]=theta[j1];
    }
    for(j1=0;j1<2;j1++){
      for(j2=0;j2<=j1;j2++){
	Umin1[j1][j2]=theta[6+(j1*j1+j1)/2+j2];
        Umin1[j2][j1]=Umin1[j1][j2];
      }
    }
    tausq=theta[9];
    if(tausq<0){
      for(j=0;j<5;j++){
	free(C[j]);
	if(j<2){
	  free(R1min1[j]);
	  free(tempmat[j]);
	  free(tempmat2[j]);
	  free(U[j]);
	  free(Umin1[j]);
	}
      }
      free(C);
      free(R1min1);
      free(tempmat);
      free(tempmat2);
      free(U);
      free(Umin1);
      lp=-10000000.0;
      *llh1=lp;
      return lp;
    }
     
  }
  else{
    for(j=0;j<5;j++){
      free(C[j]);
      if(j<2){
	free(R1min1[j]);
	free(tempmat[j]);
	free(tempmat2[j]);
	free(U[j]);
	free(Umin1[j]);
      }
    }
    free(C);
    free(R1min1);
    free(tempmat);
    free(tempmat2);
    free(U);
    free(Umin1);
    lp=-10000000.0;
    *llh1=lp;
    return lp;
  }

  
  /* first check V^(-1) or U^(-1) is positive definite  */
  if(k==0){
    for(j=0;j<3;j++){
      for(j1=0;j1<=j;j1++){
	tempmat[j][j1]=Vmin1[j][j1];
      }
    }
    posdef=1;
    chol2(3,tempmat,&posdef);
    if(posdef==0){
      for(j=0;j<5;j++){
	free(C[j]);
	if(j<3){
	  free(R0min1[j]);
	  free(tempmat[j]);
	  free(tempmat2[j]);
	  free(V[j]);
	  free(Vmin1[j]);
	}
      }
      free(C);
      free(R0min1);
      free(tempmat);
      free(tempmat2);
      free(V);
      free(Vmin1);
      lp=-10000000.0;
      *llh1=lp;
      return lp;
    }
  }
  else{
    for(j=0;j<2;j++){
      for(j1=0;j1<=j;j1++){
	tempmat[j][j1]=Umin1[j][j1];
      }
    }
    posdef=1;
    chol2(2,tempmat,&posdef);
    if(posdef==0){
      for(j=0;j<5;j++){
	free(C[j]);
	if(j<2){
	  free(R1min1[j]);
	  free(tempmat[j]);
	  free(tempmat2[j]);
	  free(U[j]);
	  free(Umin1[j]);
	}
      }
      free(C);
      free(R1min1);
      free(tempmat);
      free(tempmat2);
      free(U);
      free(Umin1);
      lp=-10000000.0;
      *llh1=lp;
      return lp;
    }
  }

  /* work out V or U (will be used in likelihood) */ 

  if(k==0){

    /* tempmat is matrix sqrt of V^(-1). tempmat2 will be inverse of this
       so that ((tempmat2)^T)(tempmat2)=V 
       and (tempmat)((tempmat)^T)=V^(-1) (same with U below) */ 

    for(j=0;j<3;j++){
      for(j1=0;j1<3;j1++){
	tempmat2[j][j1]=0.0;
      }
      tempmat2[j][j]=1.0;
    }
    for(j=0;j<3;j++){
      for(j1=0;j1<=j;j1++){
	for(j2=0;j2<j;j2++){
	  tempmat2[j][j1]-=tempmat[j][j2]*tempmat2[j2][j1];
	}
	tempmat2[j][j1]/=tempmat[j][j];
      }
    }
    for(j=0;j<3;j++){
      for(j1=0;j1<3;j1++){
        V[j][j1]=0.0;
        for(j2=0;j2<3;j2++){
	  V[j][j1]+=tempmat2[j2][j]*tempmat2[j2][j1];
	}
      }
    }
  }
  else{
    for(j=0;j<2;j++){
      for(j1=0;j1<2;j1++){
	tempmat2[j][j1]=0.0;
      }
      tempmat2[j][j]=1.0;
    }
    for(j=0;j<2;j++){
      for(j1=0;j1<=j;j1++){
	for(j2=0;j2<j;j2++){
	  tempmat2[j][j1]-=tempmat[j][j2]*tempmat2[j2][j1];
	}
	tempmat2[j][j1]/=tempmat[j][j];
      }
    }
    for(j=0;j<2;j++){
      for(j1=0;j1<2;j1++){
        U[j][j1]=0.0;
        for(j2=0;j2<2;j2++){
	  U[j][j1]+=tempmat2[j2][j]*tempmat2[j2][j1];
	}
      }
    }
  }



  /* Start with prior contribution */
  
  if(k==0){
    lp=1.0; 
    /* Normal for alpha */
    for(j=0;j<9;j++){
      lp*=D0min1[j];
    }
    lp=0.5*log(lp);
    lp-=4.5*log(tpi); 
    for(j=0;j<9;j++){
      lp-=0.5*(alpha[j]-c0[j])*(alpha[j]-c0[j])*D0min1[j];      
    }

    /* Wishart for V^(-1) */

    detVmin1=pow(det2(3,tempmat),2.0);

    lp+=((rho-3.0-1.0)/2.0)*log(detVmin1);
    for(j=0;j<3;j++){
      lp-=0.5*rho*R0[j][j]*Vmin1[j][j];
    }
    lp-=(rho/2.0)*log(pow(rho,-3.0)*det2(3,R0min1));
    lp-=(rho*3.0/2.0)*log(2.0);
    lp-=(3.0*2.0/4.0)*log(pi); 
    for(j=0;j<3;j++){
      lp-=loggamma((double)(rho-j)/2.0);
    }
    /* Inverse gamma for sigmasq */
    lp+=(-(a+1.0)*log(sigmasq)-a*log(b)-1.0/(b*sigmasq)-loggamma(a));
  }
  else{
    lp=1.0; 
    /* Normal for gamma */
    for(j=0;j<6;j++){
      lp*=D1min1[j];
    }
    lp=0.5*log(lp);
    lp-=3.0*log(tpi); 
    for(j=0;j<6;j++){
      lp-=0.5*(gamma[j]-c1[j])*(gamma[j]-c1[j])*D1min1[j];      
    }

    /* Wishart for U^(-1) */
    
    detUmin1=pow(det2(2,tempmat),2.0);
  
    lp+=((rho-2.0-1.0)/2.0)*log(detUmin1);
    for(j=0;j<2;j++){
      lp-=0.5*rho*R1[j][j]*Umin1[j][j];
    }
    lp-=(rho/2.0)*log(pow(rho,-2.0)*det2(2,R1min1));
    lp-=(rho*2.0/2.0)*log(2.0);
    lp-=(2.0*1.0/4.0)*log(pi); 
    for(j=0;j<2;j++){
      lp-=loggamma((double)(rho-j)/2.0);
    }
    /* Inverse gamma for tausq */
    lp+=(-(a+1.0)*log(tausq)-a*log(b)-1.0/(b*tausq)-loggamma(a));
  }


  /* Look at likelihood contribution */
  llh=0.0;
  if(k==0){
    for(i=0;i<467;i++){
      for(j=0;j<3;j++){
        for(j1=0;j1<5;j1++){
	  temp[j][j1]=0.0;
	}
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<3;j2++){
	    temp[j][j1]+=V[j][j2]*W[i][j1][j2];
	  }
	}
      }
      for(j=0;j<5;j++){
	for(j1=0;j1<5;j1++){
	  C[j][j1]=0.0;
	}
      }
      for(j=0;j<S[i];j++){
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<3;j2++){
	    C[j][j1]+=W[i][j][j2]*temp[j2][j1];
	  }
	}
      }
      for(j=0;j<S[i];j++){
	C[j][j]+=sigmasq;
      }
      posdef=1;
      chol2(S[i],C,&posdef);
      if(posdef==0){
	for(j=0;j<5;j++){
	  free(C[j]);
	  if(j<3){
	    free(R0min1[j]);
	    free(tempmat[j]);
	    free(tempmat2[j]);
	    free(V[j]);
	    free(Vmin1[j]);
	  }
	}
	free(C);
	free(R0min1);
	free(tempmat);
	free(tempmat2);
	free(V);
	free(Vmin1);	
	lp=-10000000.0;
	*llh1=lp;
	return lp;
      }
  
      for(j=0;j<S[i];j++){
	mu[j]=0;
        for(j1=0;j1<9;j1++){
	  mu[j]+=X[i][j][j1]*alpha[j1];
	}
      }

      j1=0;
      for(j=0;j<5;j++){
	if(Y[i][j]<90){
	  datai[j1]=Y[i][j];
	  j1++;
	}
      }

      llh+=lnormprob2(S[i],C,mu,datai);     
    }
  }
  else{
    for(i=0;i<467;i++){
      for(j=0;j<2;j++){
        for(j1=0;j1<5;j1++){
	  temp[j][j1]=0.0;
	}
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<2;j2++){
	    temp[j][j1]+=U[j][j2]*Q[i][j1][j2];
	  }
	}
      }
      for(j=0;j<5;j++){
	for(j1=0;j1<5;j1++){
	  C[j][j1]=0.0;
	}
      }
      for(j=0;j<S[i];j++){
	for(j1=0;j1<S[i];j1++){
	  for(j2=0;j2<2;j2++){
	    C[j][j1]+=Q[i][j][j2]*temp[j2][j1];
	  }
	}
      }
      for(j=0;j<S[i];j++){
	C[j][j]+=tausq;
      }
      posdef=1;
      chol2(S[i],C,&posdef);
      if(posdef==0){
	for(j=0;j<5;j++){
	  free(C[j]);
	  if(j<2){
	    free(R1min1[j]);
	    free(tempmat[j]);
	    free(tempmat2[j]);
	    free(U[j]);
	    free(Umin1[j]);
	  }
	}
	free(C);
	free(R1min1);
	free(tempmat);
	free(tempmat2);
	free(U);
	free(Umin1);
	lp=-10000000.0;
	*llh1=lp;
	return lp;
      }

      for(j=0;j<S[i];j++){
	mu[j]=0;
        for(j1=0;j1<6;j1++){
	  mu[j]+=P[i][j][j1]*gamma[j1];
	}
      }
      j1=0;
      for(j=0;j<5;j++){
	if(Y[i][j]<90){
	  datai[j1]=Y[i][j];
	  j1++;
	}
      }

      llh+=lnormprob2(S[i],C,mu,datai);     
    }
  }
  *llh1=llh;
  lp+=llh;

  if(k==0){
    for(j=0;j<5;j++){
      free(C[j]);
      if(j<3){
	free(R0min1[j]);
	free(tempmat[j]);
	free(tempmat2[j]);
	free(V[j]);
        free(Vmin1[j]);
      }
    }
    free(C);
    free(R0min1);
    free(tempmat);
    free(tempmat2);
    free(V);
    free(Vmin1);
  }
  else{
    for(j=0;j<5;j++){
      free(C[j]);
      if(j<2){
	free(R1min1[j]);
	free(tempmat[j]);
	free(tempmat2[j]);
	free(U[j]);
	free(Umin1[j]);
      }
    }
    free(C);
    free(R1min1);
    free(tempmat);
    free(tempmat2);
    free(U);
    free(Umin1);
  } 
	  
  return lp;
      
} 

double boxm2(){

  /*Function for returning single N(0,1) variable */

  double u1,u2,out;

  u1=sdrand();
  u2=sdrand();
  u1=tpi*u1;
  u2=sqrt(-2.0*log(u2));
  out=u2*cos(u1);

  return out;

}

void chol2(int nkk, double **A,int *posdef){

  /* Function for performing cholesky decompositon of B and returns result 
   in the same matrix - adapted from PJG Fortran function*/  
    
  int j1,j2,j3;
  double sum;

  for(j1=0;j1<nkk;j1++){
    sum=A[j1][j1];
    for(j2=0;j2<j1;j2++){
      sum-=pow(A[j1][j2],2);
    }
    if(sum<0){
      *(posdef)=0;
	return;
    }
    A[j1][j1]=sqrt(sum);
  
    for(j2=j1+1;j2<nkk;j2++){
      sum=A[j2][j1];
      for(j3=0;j3<j1;j3++){
	sum-=A[j2][j3]*A[j1][j3];
      }
      A[j2][j1]=sum/A[j1][j1];
    }
  }
}

double det2(int nkk,double **B){

  /* Function for evaluating determinant of an nkk by nkk matrix B */

  int j1;
  double out;
  out=1.0;
  for(j1=0;j1<nkk;j1++){
    out*=B[j1][j1];
  }
  return out;
}

double lnormprob2(int nkk, double **B, double *mu, double *datai){

  /* Function for evaluating the log of p.d.f. of Multivariate normal
   with mean mu and sqrt of cov matrix B, at pt datai */

  int j1,j2;
  double work[nkk];
  double out;
  
  for(j1=0;j1<nkk;j1++){
    work[j1]=datai[j1]-mu[j1];
  }
  for(j1=0;j1<nkk;j1++){
    for(j2=0;j2<j1;j2++){
      (work[j1])-=B[j1][j2]*work[j2];
    }
    (work[j1])/=B[j1][j1]; 
  }
  out=0.0;
  for(j1=0;j1<nkk;j1++){
    out+=(work[j1]*work[j1]);
  }
  out=-0.5*out;
  out-=(nkk/2.0)*log(tpi);
  out-=log(det2(nkk,B));    
  return out;

}
