/* The AutoMix program. 

   Last edited 25/11/04.
   Developed by David Hastie, Department of Mathematics,
   University of Bristol, UK as a part of a submission for the 
   degree of Ph.D. This Ph.D. was supervised by Prof. Peter Green (PJG), 
   University of Bristol. Special thanks also to Dr. Christophe Andrieu CA), 
   University of Bristol, for advice on adaptive schemes and mixture fitting.

   The AutoMix sampler is free for personal and academic use, but must 
   reference the sampler as instructed below.  For commercial
   use please permission must be sought from the author. To seek permission
   for such use please send an e-mail to d_hastie@hotmail.com 
   outlining the desired usage.  

   Use of the AutoMix sampler is entirely at the user's own risk. It is the
   responsibility of the user to ensure that any conclusions made through the 
   use of the AutoMix sampler are valid. The author accepts no responsibility 
   whatsoever for any loss, financial or otherwise, that may arise in 
   connection with the use of the sampler.   

   The AutoMix sampler is available from http://www.davidhastie.me.uk/AutoMix
   Although the sampler may be modified and redistributed, the author
   encourages users to register at the above site so that updates of the 
   software can be received.

   Before use, please read the README file bundled with this software.   
   
   Users should reference the sampler as instructed on the AutoMix website
   (see above). Initially this is likely to be the Ph.D. thesis that 
   introduces the AutoMix sampler. However, this will hopefully change to 
   be a published paper in the not too distant future.  */

/* Standard library files */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>

#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))

/* Global constants (please feel free to change as required)

   nkmaxmax = maximum dimension of any one model under consideration
   kmaxmax = maximum number of models
   Lkmaxmax = initial number of mixture components fitted in stage 2 of 
              AutoMix algorithm */

#define nkmaxmax 20
#define kmaxmax 15
#define Lkmaxmax 30
#define tpi 6.283185307179586477
#define pi 3.141592653589793238
#define logrtpi 0.5*log(tpi)

/* --- Internal functions (described below) ----------------- */

void gauss(double *z,int nkk);

void rt(double *z,int nkk,int dof);

void chol(int nkk,double **B);

void perm(double *work,int nkk);

double ltprob(int dof,double z,double *constt);

double lnormprob(int k,int nkk,int l,double ***mu, 
		double ****B, double *datai);

double det(int k, int nkk, int l, double ****B);

/* --- External functions -----------------*/
  
/* Random number functions sdrni and sdrand can be found in sd.c
   bundled with this software. References for these functions can be 
   found within sd.c */

extern void sdrni(unsigned long *seed); 

extern double sdrand();

/* Functions rgamma and loggamma can be found in gammafns.c bundled with 
   with this software. */

extern double rgamma(double s);

extern double loggamma(double s);

/* Function sokal is found in file sokal.c bundled with this software. */

extern void sokal(int n, double *xreal, double *var, double *tau, int *m);

/* --- User supplied functions ------------ */

/* Functions must be supplied in user***.c file (see e.g. usertoy1.c, 
   usercpt.c etc bundled with this software). 

   Descriptions: 
   1. lpost(&k,theta,&llh) 
   This should be a c function written by the user that evaluates 
   log posterior (up to an additive constant) at (k,theta). The function 
   can also return the likelihood at this point in llh.
   2. getkmax(&kmax)
   This should be a c function written by the user that returns the 
   number of models kmax.
   3. getnk(kmax, nk)
   This should be a c function written by the user that returns the dimensions
   nk for model k=1,...,kmax.
   4. getic(k,nkk,rwm)
   This should be a c function written by the user that returns the 
   possibly random starting point for the rwm in stage 1 of the AutoMix 
   sampler */

extern double lpost(int k, int nkk, double *theta, double *llh);

extern void getkmax(int *kmax);
     
extern void getnk(int kmax,int nk[kmax]);

extern void getic(int k, int nkk, double *rwm);

/* ---main program-------------------------- */

int main(int argc,char *argv[]){

  /*---Section 1 Declare Variables -------------------------*/
  
  /* ---clock variables ---------------------- */
  clock_t starttime,endtime;
  double timesecs;
  
  /* ---indexing variables ------------------- */
  int t1,t2,i1,i2,j1,j2,k1,l1,l2,sweep,remain; 
  
  /* ---counting variables ------------------- */
  int nsweep,count,nsweep2,naccrwmb,naccrwms,nacctd,ntryrwmb,ntryrwms,ntrytd;
  int nburn,nsokal,nkeep,keep,nsweepr,*nacc,*ntry,*ksummary;
  
  /* ---command line reading parameters ------ */
  int numargs,sametest; 
  char word[20],selector[3],iparam[18];    
  
  /* ---filename variables ------------------- */
  int check;
  FILE *fpk,*fpl,*fpt[kmaxmax],*fpcf,*fpmix,*fplp,*fpp,*fpac,*fpad;
  char fname[18],fname1[18],kno[6];
  
  /* ---random no. variables ----------------- */ 
  unsigned long seed;
  double u,constt;
  int dof;

  /* ---logical variables -------------------- */
  int doperm;
 
  /* ---State parameters and variables ------- */
  int k,kn,kmax,nkk,nkkn,nkmax,lendata;
  int *nk;
  double *theta,*thetan;
  double **data,*propk,*pk;

  /* ---Mixture parameters --------------------*/
  int l,ln,Lkk,Lkkn,Lkkmin,nparams,ldel,Lkmax;
  int *Lk;
  double ***mu,****B,****BBT,**detB;
  double ***mumin,****Bmin,****BBTmin;
  double **lambda,**lambdamin,**w,*logw,**lpdatagivenl,*palloc,*pallocn;
  double tol=0.00001,costfn,costfnnew,costfnmin,minlambda;

  /* ---RWM parameters ------------------------*/
  double *rwm,*rwmn,Z[1],*Znkk;
  double **sig,gamma,accept,alphastar=0.25;

  /* ---Probabilities ------------------------ */
  double lp,lpn,logratio,llh,llhn;

  /* ---working arrays and variables --------- */
  int indic,stop,natann,forceann,mode;
  int *init;
  double sum,sigma,wnew,wnewl1,*sumw,sumwnew,sumlambda,*work,thresh;
  double *datamean,**M1;

  /* ---autocorrelation variables ------------ */
  double *xr,var,tau;
  int m;

  /* ---adaptation parameters ---------------- */
  int adapt,reinit,nreinit;
  double pkllim;


  /* --- Section 2 - Read in Comand Line Variables ----------------- */

  starttime=clock(); 
 
  /* Definition of command line variables and explanation
  
     Prog variable ~ Command line variable ~ Explanation 

     mode ~ m ~ 0 if mixture fitting, 1 if user supplied mixture params, 
                2 if AutoRJ
     nsweep ~ N ~ no. of reversible jump sweeps in stage 3
     nsweep2 ~ n ~ max(n,10000*nk,100000) sweeps in within-model RWM in stage 1
     adapt ~ a ~ 1 if RJ adaptation done in stage 3, 0 otherwise
     doperm ~ p ~ 1 if random permutation done in stage 3 RJ move, 0 otherwise 
     seed ~ s ~ random no. seed, 0 uses clock seed
     dof ~ t ~ 0 if Normal random variables used in RWM and RJ moves, otherwise
               specify integer degrees of freedom of student t variables
     fname ~ f ~ filename base */


  /* Default values */

  nsweep=100000; 
  nsweep2=100000;
  numargs=argc-1;
  strcpy(fname,"output");
  doperm=1;
  seed=0;
  mode=0; 
  adapt=1;  
  dof=0;

  /* Override defaults if user supplies command line options */

  if(numargs>0){
    for(t1=1;t1<=numargs;t1++){
      
      strcpy(word,argv[t1]);
      for(t2=0;t2<2;t2++){
	selector[t2]=word[t2];
      }
      selector[2]='\0';
      for(t2=0;t2<17;t2++){
        iparam[t2]=word[t2+2];
      }
      iparam[17]='\0';
      
      
      sametest=strcmp(selector,"-f");  
      if(sametest==0){ 
	strcpy(fname,iparam);
	continue;
      }
      sametest=strcmp(selector,"-N");  
      if(sametest==0){ 
	nsweep=atoi(iparam);
        continue;
      }

      sametest=strcmp(selector,"-n");  
      if(sametest==0){ 
	nsweep2=max(atoi(iparam),100000);
        continue;
      }

      sametest=strcmp(selector,"-s");
      if(sametest==0){
        seed=atoi(iparam);
        continue;
      }
      sametest=strcmp(selector,"-p");
      if(sametest==0){
        doperm=atoi(iparam);
        continue;
      }
      sametest=strcmp(selector,"-m");
      if(sametest==0){
        mode=atoi(iparam);
        continue;
      }

      sametest=strcmp(selector,"-a");
      if(sametest==0){
        adapt=atoi(iparam);
        continue;
      }

      sametest=strcmp(selector,"-t");
      if(sametest==0){
        dof=atoi(iparam);
        continue;
      }

    }
  }

  sdrni(&seed);

    
  /* --- Section 3 - Initial File handling ---------------------  */
    
  sprintf(fname1,fname);
  strcat(fname1,"_log.data");  
  fpl = fopen(fname1,"w");
  sprintf(fname1,fname);
  strcat(fname1,"_pk.data");  
  fpp = fopen(fname1,"w");
  sprintf(fname1,fname);
  strcat(fname1,"_ac.data");  
  fpac = fopen(fname1,"w");
  sprintf(fname1,fname);
  strcat(fname1,"_adapt.data");  
  fpad = fopen(fname1,"w");
  sprintf(fname1,fname);
  strcat(fname1,"_cf.data");  
  fpcf = fopen(fname1,"w");

  /* Print user options to log file */

  fprintf(fpl,"seed: %ld\n",seed);
  fprintf(fpl,"m: %d\n",mode);
  fprintf(fpl,"a: %d\n",adapt);
  fprintf(fpl,"p: %d\n",doperm); 
  fprintf(fpl,"n: %d\n",nsweep2);
  fprintf(fpl,"N: %d\n",nsweep);

  /* Check user has supplied mixture parameters if trying to use mode 1.
     If not default back to mode 0 */

  if(mode==1){
    sprintf(fname1,fname);
    strcat(fname1,"_mix.data");  
    if((fpmix = fopen(fname1,"r"))==NULL){
      printf("\nMixture file doesn't exist:");
      printf("\nContinuing using RWM to estimate parameters");
      mode=0;
    }
  }  
  
  /* --- Section 4.0 - Read in key variables from user functions -*/ 

  getkmax(&kmax);
  if(kmax>kmaxmax){
    printf("\nError:kmax too large \n");
    return 0;
  }
  else if(kmax<0){
    printf("\nError:negative kmax \n");
    return 0;
  }

  nk=(int*)malloc(kmax*sizeof(int));
  Lk=(int*)malloc(kmax*sizeof(int));
  ksummary=(int*)malloc(kmax*sizeof(int));

  getnk(kmax,nk);  
  nkmax=nk[0];
  ksummary[0]=0;
  for(k1=1;k1<kmax;k1++){
    nkmax=max(nk[k1],nkmax);
    ksummary[k1]=0;
  }

  lambda=(double**)malloc(kmax*sizeof(double));
  lambdamin=(double**)malloc(kmax*sizeof(double));
  mu=(double***)malloc(kmax*sizeof(double));
  mumin=(double***)malloc(kmax*sizeof(double));
  BBT=(double****)malloc(kmax*sizeof(double));
  BBTmin=(double****)malloc(kmax*sizeof(double));
  B=(double****)malloc(kmax*sizeof(double));
  Bmin=(double****)malloc(kmax*sizeof(double));
  detB=(double**)malloc(kmax*sizeof(double));
  sig=(double**)malloc(kmax*sizeof(double));
  for(k1=0;k1<kmax;k1++){
    nkk=nk[k1];
    lambda[k1]=(double*)malloc(Lkmaxmax*sizeof(double));
    lambdamin[k1]=(double*)malloc(Lkmaxmax*sizeof(double));
    mu[k1]=(double**)malloc(Lkmaxmax*sizeof(double));
    mumin[k1]=(double**)malloc(Lkmaxmax*sizeof(double));
    BBT[k1]=(double***)malloc(Lkmaxmax*sizeof(double));
    BBTmin[k1]=(double***)malloc(Lkmaxmax*sizeof(double));
    B[k1]=(double***)malloc(Lkmaxmax*sizeof(double));
    Bmin[k1]=(double***)malloc(Lkmaxmax*sizeof(double));
    detB[k1]=(double*)malloc(Lkmaxmax*sizeof(double));
    sig[k1]=(double*)malloc(nkk*sizeof(double));
    for(l1=0;l1<Lkmaxmax;l1++){
      mu[k1][l1]=(double*)malloc(nkk*sizeof(double));
      mumin[k1][l1]=(double*)malloc(nkk*sizeof(double));
      BBT[k1][l1]=(double**)malloc(nkk*sizeof(double));
      BBTmin[k1][l1]=(double**)malloc(nkk*sizeof(double));
      B[k1][l1]=(double**)malloc(nkk*sizeof(double));
      Bmin[k1][l1]=(double**)malloc(nkk*sizeof(double));
      for(j1=0;j1<nkk;j1++){
	BBT[k1][l1][j1]=(double*)malloc(nkk*sizeof(double));
	BBTmin[k1][l1][j1]=(double*)malloc(nkk*sizeof(double));
	B[k1][l1][j1]=(double*)malloc(nkk*sizeof(double));
	Bmin[k1][l1][j1]=(double*)malloc(nkk*sizeof(double));
      }
    }
  } 
        
  
  /* --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) --- */

  /* These parameters are used if mode 1 (m=1) of the AutoMix sampler is 
     used. Note that if the parameters are unavailable or inconsistent with 
     the user supplied functions or unavailable go straight to section 5.2 
     where initial within-model RWM are performed */

  if(mode>2){
    printf("\nInvalid mode entered. Mode must be 0,1,2");
    return -100;
  }
  else if(mode==1){
    if((check=fscanf(fpmix,"%d",&k1))==EOF){
      printf("\nEnd of file encountered before parameters read:"); 
      printf("\nContinuing using RWM to estimate parameters");
      mode=0;     
      goto RWMSTART;  
    }
    if(k1!=kmax){     
      printf("\nFile kmax contradicts getkmax function:");
      printf("\nContinuing using RWM to estimate parameters");
      mode=0;     
      goto RWMSTART;     
    }
    for(k1=0;k1<kmax;k1++){
      if((check=fscanf(fpmix,"%d",&j1))==EOF){
	printf("\nEnd of file encountered before parameters read:"); 
	printf("\nContinuing using RWM to estimate parameters");
	mode=0;     
	goto RWMSTART;
      }
      if(j1!=nk[k1]){
	printf("\nFile kmax contradicts getnk function:");
	printf("\nContinuing using RWM to estimate parameters");
	mode=0;     
	goto RWMSTART;     
      }
    }
    for(k1=0;k1<kmax;k1++){
      nkk=nk[k1];
      for(j1=0;j1<nkk;j1++){
	if((check=fscanf(fpmix,"%lf",&(sig[k1][j1])))==EOF){
	  printf("\nEnd of file encountered before parameters read:"); 
	  printf("\nContinuing using RWM to estimate parameters");
	  mode=0;     
	  goto RWMSTART;
	}      
      }
      if((check=fscanf(fpmix,"%d",&(Lk[k1])))==EOF){
	printf("\nEnd of file encountered before parameters read:"); 
	printf("\nContinuing using RWM to estimate parameters");
	mode=0;     
	goto RWMSTART;
      }
      Lkk=Lk[k1];
      for(l1=0;l1<Lkk;l1++){
	if((check=fscanf(fpmix,"%lf",&(lambda[k1][l1])))==EOF){
	  printf("\nEnd of file encountered before parameters read:"); 
	  printf("\nContinuing using RWM to estimate parameters");
	  mode=0;     
	  goto RWMSTART;
	}
	for(j1=0;j1<nkk;j1++){
	  if((check=fscanf(fpmix,"%lf",&(mu[k1][l1][j1])))==EOF){
	    printf("\nEnd of file encountered before parameters read:"); 
	    printf("\nContinuing using RWM to estimate parameters");
	    mode=0;     
	    goto RWMSTART;
	  }
	}
	for(j1=0;j1<nkk;j1++){
	  for(j2=0;j2<=j1;j2++){
	    if((check=fscanf(fpmix,"%lf",&(B[k1][l1][j1][j2])))==EOF){
	      printf("\nEnd of file encountered before parameters read:"); 
	      printf("\nContinuing using RWM to estimate parameters");
	      mode=0;     
	      goto RWMSTART;
	    }
	  }
	}
      }
      sumlambda=0.0;
      for(l1=0;l1<Lkk;l1++){
	sumlambda+=lambda[k1][l1];
      }
      if(sumlambda<0.99999||sumlambda>1.00001){
	printf("\nComponents weights read do not sum to one for k=%d:",k1); 
	printf("\nContinuing using RWM to estimate parameters");
	mode=0;     
	goto RWMSTART;
      }
      if(sumlambda<1.0||sumlambda>1.0){
	for(l1=0;l1<Lkk;l1++){
	  lambda[k1][l1]/=sumlambda;
	}
      }
    }
    if(!(fpmix==NULL)){
      fclose(fpmix);
    }

  }
  else if(mode==0||mode==2){
  
  /* --- Section 5.2 - Within-model runs if mixture parameters unavailable -*/ 

  RWMSTART:

    for(k1=0;k1<kmax;k1++){
      /* --- Section 5.2.1 - RWM Within Model (Stage 1) -------*/
              
      nkk=nk[k1];
      nparams=nkk+(nkk*(nkk+1))/2;
      lendata=1000*nkk; 
      nsweepr=max(nsweep2,10000*nkk);
      nburn=nsweepr/10;     
      nsweepr+=nburn;
      data=(double**) malloc(lendata*sizeof(double));
      for(i1=0;i1<lendata;i1++){
	data[i1]=(double*)malloc(nkk*sizeof(double));
      }
      rwm=(double*)malloc(nkk*sizeof(double));   
      rwmn=(double*)malloc(nkk*sizeof(double));
      nacc=(int*)malloc(nkk*sizeof(int));
      ntry=(int*)malloc(nkk*sizeof(int));
      Znkk=(double*)malloc(nkk*sizeof(double));
      init=(int*)malloc(nkk*sizeof(int));
      
      printf("\nRWM for Model %d",k1+1);
      fprintf(fpl,"\nRWM for Model %d",k1+1);
      fprintf(fpcf,"RWM for Model %d\n",k1+1);
      fprintf(fpad,"RWM for Model %d\n",k1+1);
      fflush(NULL);
      getic(k1,nkk,rwm);
      for(j1=0;j1<nkk;j1++){
	rwmn[j1]=rwm[j1];
        sig[k1][j1]=10.0;       
	nacc[j1]=0;
	ntry[j1]=0;
      } 
      lp=lpost(k1,nkk,rwm,&llh);

      i2=0;
      remain=nsweepr;
      for(sweep=1;sweep<=nsweepr;sweep++){
        remain--;
	if((sweep>=nburn)&&(fmod((sweep-nburn),((nsweepr-nburn)/10))<tol)){
	  printf("\nNo. of iterations remaining: %d",remain);
	  fflush(NULL);
	}
        u=sdrand();
	if(sweep>nburn&&u<0.1){     	    
	  if(dof>0){
	    rt(Znkk,nkk,dof);
	  }
	  else{
	    gauss(Znkk,nkk);
	  }
	  for(j1=0;j1<nkk;j1++){
	    rwmn[j1]=rwm[j1]+sig[k1][j1]*Znkk[j1];
	  }
	  lpn=lpost(k1,nkk,rwmn,&llhn);
	  if(sdrand()<exp(max(-30.0,min(0.0,lpn-lp)))){
	    for(j1=0;j1<nkk;j1++){
	      rwm[j1]=rwmn[j1];
	    }
	    lp=lpn;
            llh=llhn;
	  }
	}
	else{
	  gamma=10.0*pow(1.0/(sweep+1),2.0/3.0);
	  for(j1=0;j1<nkk;j1++){
	    rwmn[j1]=rwm[j1];
	  }
	  for(j1=0;j1<nkk;j1++){
	    if(dof>0){
	      rt(Z,1,dof);
	    }
	    else{
	      gauss(Z,1);
	    }
	    rwmn[j1]=rwm[j1]+sig[k1][j1]*Z[0];	
	    lpn=lpost(k1,nkk,rwmn,&llhn);
	    accept=min(1,exp(max(-30.0,min(0.0,lpn-lp))));
	    if(sdrand()<accept){
              (nacc[j1])++;
              (ntry[j1])++; 
	      rwm[j1]=rwmn[j1];
	      lp=lpn;
              llh=llhn;         
	      sig[k1][j1]=max(0,sig[k1][j1]-gamma*(alphastar-1));
	    }
	    else{
              (ntry[j1])++;
	      rwmn[j1]=rwm[j1];
	      sig[k1][j1]=max(0,sig[k1][j1]-gamma*(alphastar));
	    }	   
	  }
	}
	if(remain<(10000*nkk)&&fmod(remain,10.0)<0.05){
	  for(j1=0;j1<nkk;j1++){
	    data[i2][j1]=rwm[j1];
	  }
	  i2++;
	}
        if(fmod(sweep,100.0)<0.05){
          for(j1=0;j1<nkk;j1++){
	    fprintf(fpad,"%lf %lf ",sig[k1][j1],
		    (double)nacc[j1]/(double)ntry[j1]);
          }
          fprintf(fpad,"\n"); 
	}
	
      }
      free(init);
      free(rwm);
      free(rwmn);
      free(nacc);
      free(ntry);
      free(Znkk);
      
      /* --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)- */ 
      /* Note only done if mode 0 (m=0) if mode m=2, go to section 5.2.3*/
      /* Mixture fitting done component wise EM algorithm described in 
	 Figueiredo and Jain, 2002 (see thesis for full reference) */ 

      printf("\nMixture Fitting: Model %d",k1+1);      
      if(mode==0){
	Lkk=Lkmaxmax;
	init=(int*)malloc(Lkk*sizeof(int));
	l1=0;
	while(l1<Lkk){
	  indic=0;
	  u=sdrand(); 
	  init[l1]=(int)floor(lendata*u);
	  if(l1>0){
	    for(l2=0;l2<l1;l2++){
	      if(init[l2]==init[l1]){
		indic=1;
		break;
	      }
	    }
	  }
	  if(indic==0){
	    l1++;
	  }
	}     
	
	datamean=(double*)malloc(nkk*sizeof(double));
	M1=(double**)malloc(nkk*sizeof(double));
	for(j1=0;j1<nkk;j1++){
	  M1[j1]=(double*)malloc(nkk*sizeof(double));
	}
	for(j1=0;j1<nkk;j1++){
	  datamean[j1]=0.0;
	  for(i1=0;i1<lendata;i1++){
	    datamean[j1]+=data[i1][j1];
	  }
	  datamean[j1]/=((double)lendata);
	}
	for(j1=0;j1<nkk;j1++){
	  for(j2=0;j2<nkk;j2++){
	    M1[j1][j2]=0;
	    for(i1=0;i1<lendata;i1++){
	      M1[j1][j2]+=(data[i1][j1]-datamean[j1])*
		(data[i1][j2]-datamean[j2]);
	    }
	    M1[j1][j2]/=((double)lendata);
	  }
	}
	sigma=0.0;
	for(j1=0;j1<nkk;j1++){
	  sigma+=M1[j1][j1];
	}
	sigma/=(10.0*nkk);
	
	for(l1=0;l1<Lkk;l1++){
	  for(j1=0;j1<nkk;j1++){
	    mu[k1][l1][j1]=data[init[l1]][j1];
	    BBT[k1][l1][j1][j1]=sigma;
	    B[k1][l1][j1][j1]=BBT[k1][l1][j1][j1];
	    for(j2=0;j2<j1;j2++){
	      BBT[k1][l1][j1][j2]=0.0;  
	      B[k1][l1][j1][j2]=BBT[k1][l1][j1][j2];
	    }
	  }
	  chol(nkk,B[k1][l1]);
	  lambda[k1][l1]=1.0/Lkk;
	}
	
	w=(double**)malloc(lendata*sizeof(double));
	logw=(double*)malloc(Lkk*sizeof(double));
	lpdatagivenl=(double**)malloc(lendata*sizeof(double));
	for(i1=0;i1<lendata;i1++){
	  w[i1]=(double*)malloc(Lkk*sizeof(double));
	  lpdatagivenl[i1]=(double*)malloc(Lkk*sizeof(double));
	}
	
	for(i1=0;i1<lendata;i1++){
	  sum=0.0;
	  for(l1=0;l1<Lkk;l1++){
	    lpdatagivenl[i1][l1]=lnormprob(k1,nkk,l1,mu,B,data[i1]);
	    logw[l1]=log(lambda[k1][l1])+lpdatagivenl[i1][l1];
	    w[i1][l1]=exp(logw[l1]);
	    sum+=w[i1][l1];
	  }
	  for(l1=0;l1<Lkk;l1++){
	    w[i1][l1]/=sum;
	  }
	}
	
	sumw=(double*)malloc(Lkk*sizeof(double));
	
	stop=0;
	count=0;
	
	while(!stop){
	  count++;
	  l1=0;
	  natann=0;
	  forceann=0;
	  while(l1<Lkk){
	    sumwnew=0.0;
	    for(l2=0;l2<Lkk;l2++){
	      sumw[l2]=0.0;
	      for(i1=0;i1<lendata;i1++){
		sumw[l2]+=w[i1][l2];
	      }
	      wnew=max(0.0,(sumw[l2]-nparams/2.0));
	      if(l2==l1){
		wnewl1=wnew;
	      }
	      sumwnew+=wnew;
	    }
	    lambda[k1][l1]=wnewl1/sumwnew;
	    sumlambda=0.0;
	    for(l2=0;l2<Lkk;l2++){
	      sumlambda+=lambda[k1][l2];
	    }
	    for(l2=0;l2<Lkk;l2++){
	      lambda[k1][l2]/=sumlambda;
	    }
	    
	    if(lambda[k1][l1]>0.005){ 
	      /*changed to 0.005 from 0.0 -renormalise else */
	      for(j1=0;j1<nkk;j1++){
		mu[k1][l1][j1]=0.0;
		for(i1=0;i1<lendata;i1++){
		  mu[k1][l1][j1]+=data[i1][j1]*w[i1][l1];
		}
		mu[k1][l1][j1]/=sumw[l1];
		
		for(j2=0;j2<=j1;j2++){
		  BBT[k1][l1][j1][j2]=0.0;
		  for(i1=0;i1<lendata;i1++){
		    BBT[k1][l1][j1][j2]+=(data[i1][j1]-mu[k1][l1][j1])*
		      (data[i1][j2]-mu[k1][l1][j2])*w[i1][l1];
		  }
		  BBT[k1][l1][j1][j2]/=sumw[l1];
		  B[k1][l1][j1][j2]=BBT[k1][l1][j1][j2];
		}
	      }
	      
	      chol(nkk,B[k1][l1]);
	      
	      for(i1=0;i1<lendata;i1++){
		lpdatagivenl[i1][l1]=lnormprob(k1,nkk,l1,mu,B,data[i1]);
	      }
	      l1++;
	      
	    }
	    else{
              if(fmod(Lkk,5)<0.05){
		printf("\n");
	      }
	      printf("%d(%d-n) ",Lkk,count);
	      natann=1;
	      if(l1<(Lkk-1)){
		for(l2=l1;l2<(Lkk-1);l2++){
		  lambda[k1][l2]=lambda[k1][l2+1];
		  for(j1=0;j1<nkk;j1++){
		    mu[k1][l2][j1]=mu[k1][l2+1][j1];
		    for(j2=0;j2<=j1;j2++){
		      BBT[k1][l2][j1][j2]=BBT[k1][l2+1][j1][j2];
		      B[k1][l2][j1][j2]=B[k1][l2+1][j1][j2];
		    }
		  }
		  for(i1=0;i1<lendata;i1++){
		    lpdatagivenl[i1][l2]=lpdatagivenl[i1][l2+1];
		  }
		}
	      }
	      Lkk--;
	      sumlambda=0.0;
	      for(l2=0;l2<Lkk;l2++){
		sumlambda+=lambda[k1][l2];
	      }
	      for(l2=0;l2<Lkk;l2++){
		lambda[k1][l2]/=sumlambda;
	      }
	      
	    }
	    
	    lpn=0.0;
	    for(i1=0;i1<lendata;i1++){
	      sum=0.0;
	      for(l2=0;l2<Lkk;l2++){
		logw[l2]=log(lambda[k1][l2])+lpdatagivenl[i1][l2];
		w[i1][l2]=exp(logw[l2]);
		sum+=w[i1][l2];
	      }
	      if(sum>0){
		for(l2=0;l2<Lkk;l2++){
		  w[i1][l2]/=sum;
		}
		lpn+=log(sum);
	      }
	      else{
		/* if no component fits point well make equally likely */
		for(l2=0;l2<Lkk;l2++){
		  w[i1][l2]=1.0/Lkk;
		}
		lpn+=(-500.0);
	      }
	    }
	  }
	
	  sum=0.0;
	  for(l1=0;l1<Lkk;l1++){
	    sum+=log(lendata*lambda[k1][l1]/12.0);
	  }
	  costfnnew=(nparams/2.0)*sum+(Lkk/2.0)*log(lendata/12.0)+
	    Lkk*(nparams+1)/2.0-lpn;
	  
	  if(count==1){
	    costfn=costfnnew;
	  }
	  if(count==1||costfnnew<costfnmin){
	    Lkkmin=Lkk;
	    costfnmin=costfnnew;
	    for(l1=0;l1<Lkk;l1++){
	      lambdamin[k1][l1]=lambda[k1][l1];
	      for(j1=0;j1<nkk;j1++){
		mumin[k1][l1][j1]=mu[k1][l1][j1];
		for(j2=0;j2<=j1;j2++){
		  Bmin[k1][l1][j1][j2]=B[k1][l1][j1][j2];
		}
	      }
	    }
	  }
	  if((fabs(costfn-costfnnew)<min(tol*fabs(costfn),0.01))&&(count>1)){
	    if(Lkk==1){
	      stop=1;
	    }
	    else{
              if(fmod(Lkk,5)<0.05){
		printf("\n");
	      }
	      printf("%d(%d-f) ",Lkk,count);
	      forceann=2;
	      minlambda=lambda[k1][0];
	      ldel=0;
	      for(l1=1;l1<Lkk;l1++){
		if(minlambda>lambda[k1][l1]){
		  minlambda=lambda[k1][l1];
		  ldel=l1;
		}
	      }
	      if(ldel<(Lkk-1)){
		for(l1=ldel;l1<(Lkk-1);l1++){
		  lambda[k1][l1]=lambda[k1][l1+1];
		  for(j1=0;j1<nkk;j1++){
		    mu[k1][l1][j1]=mu[k1][l1+1][j1];
		    for(j2=0;j2<=j1;j2++){
		      BBT[k1][l1][j1][j2]=BBT[k1][l1+1][j1][j2];
		      B[k1][l1][j1][j2]=B[k1][l1+1][j1][j2];
		    }
		  }
		  for(i1=0;i1<lendata;i1++){
		    lpdatagivenl[i1][l1]=lpdatagivenl[i1][l1+1];
		  }
		}
	      }
	      Lkk--;
	      sumlambda=0.0;
	      for(l1=0;l1<Lkk;l1++){
		sumlambda+=lambda[k1][l1];
	      }
	      for(l1=0;l1<Lkk;l1++){
		lambda[k1][l1]/=sumlambda;
	      }
	      
	      lpn=0.0;
	      for(i1=0;i1<lendata;i1++){
		sum=0.0;
		for(l2=0;l2<Lkk;l2++){
		  logw[l2]=log(lambda[k1][l2])+lpdatagivenl[i1][l2];
		  w[i1][l2]=exp(logw[l2]);
		  sum+=w[i1][l2];
		}
		if(sum>0){
		  for(l2=0;l2<Lkk;l2++){
		    w[i1][l2]/=sum;
		  }
		  lpn+=log(sum);
		}
		else{
		  /* if no component fits point well make equally likely */
		  for(l2=0;l2<Lkk;l2++){
		    w[i1][l2]=1.0/Lkk;
		  }
		  lpn+=(-500.0);
		}
		
	      }
	      
	      sum=0.0;
	      for(l1=0;l1<Lkk;l1++){
		sum+=log(lendata*lambda[k1][l1]/12.0);
	      }
	      costfnnew=(nparams/2.0)*sum+(Lkk/2.0)*log(lendata/12.0)+
		Lkk*(nparams+1)/2.0-lpn;
	    }
	  }
	  if(count>5000){
	    stop=1;
	  }
	  costfn=costfnnew;
	  fprintf(fpcf,"%d %lf %lf %d\n",Lkk,lpn,costfnnew,(natann+forceann)); 
	  fflush(NULL);
	}
      
	for(j1=0;j1<nkk;j1++){
	  free(M1[j1]);
	}
	free(M1);
	free(datamean);
	for(i1=0;i1<lendata;i1++){
	  free(data[i1]);
	  free(w[i1]);
	  free(lpdatagivenl[i1]);
	}
	free(w);
	free(lpdatagivenl);
	free(data);
	free(logw);
	free(sumw);
	free(init);
	Lk[k1]=Lkkmin;
	for(l1=0;l1<Lkkmin;l1++){
	  lambda[k1][l1]=lambdamin[k1][l1];
	  for(j1=0;j1<nkk;j1++){
	    mu[k1][l1][j1]=mumin[k1][l1][j1];
	  }
	  for(j1=0;j1<nkk;j1++){
	    for(j2=0;j2<=j1;j2++){
	      B[k1][l1][j1][j2]=Bmin[k1][l1][j1][j2];
	    }
	  }
	}
      }
      else if(mode==2){
	/* --- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --*/
	/* Note only done if mode 2 (m=2).*/
        Lk[k1]=1;
        lambda[k1][0]=1.0;
        for(j1=0;j1<nkk;j1++){
          mu[k1][0][j1]=0.0; 
	  for(i1=0;i1<lendata;i1++){
	    mu[k1][0][j1]+=data[i1][j1];
	  }
	  mu[k1][0][j1]/=((double)lendata);
	}
        for(j1=0;j1<nkk;j1++){
	  for(j2=0;j2<=j1;j2++){
	    B[k1][0][j1][j2]=0.0;
            for(i1=0;i1<lendata;i1++){
	      B[k1][0][j1][j2]+=(data[i1][j1]-mu[k1][0][j1])*
		(data[i1][j2]-mu[k1][0][j2]);
	    }
	    B[k1][0][j1][j2]/=((double)(lendata-1));
	  }
	}
	chol(nkk,B[k1][0]);

	for(i1=0;i1<lendata;i1++){
	  free(data[i1]);
	}
	free(data);
      }
    }
  }

  /* Print mixture parameters to file (log and mix files) for reference 
     and use in future runs. */

  sprintf(fname1,fname);
  strcat(fname1,"_mix.data");  
  fpmix = fopen(fname1,"w");
  fprintf(fpmix,"%d\n",kmax);
  for(k1=0;k1<kmax;k1++){
    fprintf(fpmix,"%d\n",nk[k1]);
  }

  for(k1=0;k1<kmax;k1++){
    fprintf(fpl,"\nModel:%d\n",k1+1);
    Lkk=Lk[k1];
    nkk=nk[k1];
    fprintf(fpl,"\nARW params:\n");
    for(j1=0;j1<nkk;j1++){
      fprintf(fpl,"%lf ",sig[k1][j1]);
    }
    fprintf(fpl,"\n");
    fprintf(fpl,"\nLkk:%d\n",Lkk);
    for(j1=0;j1<nkk;j1++){
      fprintf(fpmix,"%lf\n",sig[k1][j1]);
    }
    fprintf(fpmix,"%d\n",Lkk);
    for(l1=0;l1<Lkk;l1++){
      fprintf(fpl,"\nComponent:%d\n",l1+1);
      fprintf(fpl,"lambda:%lf\n",lambda[k1][l1]);
      fprintf(fpmix,"%lf\n",lambda[k1][l1]);
      fprintf(fpl,"mu:\n");	
      for(j1=0;j1<nkk;j1++){
	fprintf(fpl,"%lf ",mu[k1][l1][j1]);
	fprintf(fpmix,"%lf\n",mu[k1][l1][j1]);
      }
      fprintf(fpl,"\nB:\n");
      for(j1=0;j1<nkk;j1++){
	for(j2=0;j2<=j1;j2++){
	  fprintf(fpl,"%lf ",B[k1][l1][j1][j2]);
	  fprintf(fpmix,"%lf\n",B[k1][l1][j1][j2]);
	}
	fprintf(fpl,"\n");
      }
      detB[k1][l1]=det(k1,nkk,l1,B);
    }
  }
  fflush(NULL); 

  /* --Section 6 - Secondary file handling -------------*/

  sprintf(fname1,fname);
  strcat(fname1,"_k");
  strcat(fname1,".data");
  fpk=fopen(fname1,"w");
  sprintf(fname1,fname);
  strcat(fname1,"_lp");
  strcat(fname1,".data");
  fplp=fopen(fname1,"w");
  
  for(k1=0;k1<kmax;k1++){
    sprintf(fname1,fname);
    sprintf(kno,"%d",k1+1);
    strcat(fname1,"_theta");  
    strcat(fname1,kno);
    strcat(fname1,".data");
    fpt[k1]=fopen(fname1,"w");   
  }

  /* --Section 7 - Final initialisation of variables ----*/
    
  naccrwmb=0;
  ntryrwmb=0;
  naccrwms=0;
  ntryrwms=0;
  nacctd=0;
  ntrytd=0;

  constt=100000.0;
  Lkmax=Lk[0];
  for(k1=1;k1<kmax;k1++){
    Lkmax=max(Lkmax,Lk[k1]);
  }
  k=(int)floor(kmax*sdrand());
  nkk=nk[k];
  Lkk=Lk[k];
  
  theta=(double*)malloc(nkmax*sizeof(double));
  thetan=(double*)malloc(nkmax*sizeof(double));
  work=(double*)malloc(nkmax*sizeof(double));
  palloc=(double*)malloc(Lkmax*sizeof(double));
  pallocn=(double*)malloc(Lkmax*sizeof(double));
  Znkk=(double*)malloc(nkmax*sizeof(double));
  propk=(double*)malloc(kmax*sizeof(double));
  pk=(double*)malloc(kmax*sizeof(double));

  getic(k,nkk,theta);

  lp=lpost(k,nkk,theta,&llh);
 
  for(k1=0;k1<kmax;k1++){
    pk[k1]=1.0/kmax;
    if(k1==k){
      propk[k1]=1.0;
    }
    else{
      propk[k1]=0.0;
    }
  }
  nreinit=1;
  reinit=0;
  pkllim=1.0/10.0;
    
  tol=0.5/nsweep;
  nburn=max(10000,(int)(nsweep/10));

  nsokal=1;
  nkeep=nsweep/(2*nsokal);
  nkeep=(int)pow(2.0,min(15,(int)(log(nkeep)/log(2.0)+0.001)));
  keep=nburn+(nsweep-nkeep*nsokal);
  xr=(double*)malloc(nkeep*sizeof(double));
 

 
 /* -----Start of main loop ----------------*/
  for(sweep=1;sweep<=(nburn+nsweep);sweep++){

    /* --Section 8 - RWM within-model moves ---*/


    /* Every 10 sweeps to block RWM */
    if(fmod(sweep,10)<0.05){
      ntryrwmb++;
      if(dof>0){
	rt(Znkk,nkk,dof);
      }
      else{
	gauss(Znkk,nkk);
      }
      for(j1=0;j1<nkk;j1++){
	thetan[j1]=theta[j1]+sig[k][j1]*Znkk[j1];	
      }
      lpn=lpost(k,nkk,thetan,&llhn);
      if(sdrand()<exp(max(-30.0,min(0.0,lpn-lp)))){
	naccrwmb++;
	for(j1=0;j1<nkk;j1++){
	  theta[j1]=thetan[j1];
	}
	lp=lpn;
	llh=llhn;
      }	
    }
    else{
      /* else do component-wise RWM */
      for(j1=0;j1<nkk;j1++){
	thetan[j1]=theta[j1];
      }
      for(j1=0;j1<nkk;j1++){
	ntryrwms++;
        if(dof>0){
	  rt(Z,1,dof);
	}
	else{
	  gauss(Z,1);
	}
	thetan[j1]=theta[j1]+sig[k][j1]*Z[0];	
	lpn=lpost(k,nkk,thetan,&llhn);
	if(sdrand()<exp(max(-30.0,min(0.0,lpn-lp)))){
	  naccrwms++;
	  theta[j1]=thetan[j1];
	  lp=lpn;
          llh=llhn;
	}
	else{
	  thetan[j1]=theta[j1];
	}	   
      }
    }

    /* --- Section 9 Reversible Jump Moves -------------------*/

    /* --Section 9.1 - Allocate current position to a component --*/

    ntrytd++;
    if(Lkk>1){
      sum=0.0;
      for(l1=0;l1<Lkk;l1++){
	palloc[l1]=log(lambda[k][l1])+lnormprob(k,nkk,l1,mu,B,theta);
	palloc[l1]=exp(palloc[l1]);
	sum+=palloc[l1];
      }
      if(sum>0){
	for(l1=0;l1<Lkk;l1++){
	  palloc[l1]/=sum;       
	}
      }
      else{
	for(l1=0;l1<Lkk;l1++){
	  palloc[l1]=1.0/Lkk;
	}
      }
      u=sdrand();
      thresh=0.0;
      for(l1=0;l1<Lkk;l1++){
        thresh+=palloc[l1];
	if(u<thresh){
	  l=l1;
	  break;
	}
      }
    }
    else{
      l=0;
      palloc[l]=1.0;	
    }
       

    /* --Section 9.2 - Standardise state variable --------------- */

    for(j1=0;j1<nkk;j1++){
      work[j1]=theta[j1]-mu[k][l][j1];	
    }
    for(j1=0;j1<nkk;j1++){
      for(j2=0;j2<j1;j2++){
	work[j1]=work[j1]-B[k][l][j1][j2]*work[j2];
      }      	
      work[j1]=work[j1]/B[k][l][j1][j1];
    }

    /* --Section 9.3 - Choose proposed new model and component ----*/

    if(kmax==1){
      kn=k;
      logratio=0.0;
    }
    else{
      gamma=pow(1.0/(sweep+1),(2.0/3.0));
      u=sdrand();
      thresh=0.0;
      for(k1=0;k1<kmax;k1++){
	thresh+=pk[k1];
	if(u<thresh){
	  kn=k1;
          break;
	}
      }
      logratio=log(pk[k])-log(pk[kn]);
    }

    nkkn=nk[kn];
    Lkkn=Lk[kn];

    u=sdrand();
    thresh=0.0;
    for(l1=0;l1<Lkkn;l1++){
      thresh+=lambda[kn][l1];
      if(u<thresh){
	ln=l1;
	break;
      }
    }

    /* --Section 9.4 Propose new state ----------------*/ 

    if(nkk<nkkn){
      if(dof>0){
	rt(&(work[nkk]),nkkn-nkk,dof);
	for(j1=nkk;j1<nkkn;j1++){
	  logratio-=ltprob(dof,work[j1],&constt);
	}
      }
      else{
	gauss(&(work[nkk]),nkkn-nkk);
	for(j1=nkk;j1<nkkn;j1++){
	  logratio+=0.5*pow(work[j1],2.0)+logrtpi;
	}
      }
      if(doperm){
	perm(work,nkkn);
      }
    }
    else if(nkk==nkkn){
      if(doperm){
	perm(work,nkk);
      }
    }
    else{
      if(doperm){
	perm(work,nkk);
      }
      if(dof>0){
	for(j1=nkkn;j1<nkk;j1++){
	  logratio+=ltprob(dof,work[j1],&constt);
	}
      }
      else{
	for(j1=nkkn;j1<nkk;j1++){
	  logratio-=(0.5*pow(work[j1],2.0)+logrtpi);
	}
      }
    }
       
    for(j1=0;j1<nkkn;j1++){
      thetan[j1]=mu[kn][ln][j1];
      for(j2=0;j2<=j1;j2++){
	thetan[j1]+=B[kn][ln][j1][j2]*work[j2];
      }
    }

    /* --Section 9.5 - Work out probability of allocating to component
       for acceptance ratio (reverse move) ---------*/

    if(Lkkn>1){
      sum=0.0;
      for(l1=0;l1<Lkkn;l1++){
	pallocn[l1]=log(lambda[kn][l1])+lnormprob(kn,nkkn,l1,mu,B,thetan);
	pallocn[l1]=exp(pallocn[l1]);
	sum+=pallocn[l1];
      }
      if(sum>0){
	for(l1=0;l1<Lkkn;l1++){
	  pallocn[l1]/=sum;       
	}
      }
      else{
	for(l1=0;l1<Lkkn;l1++){
	  pallocn[l1]=1.0/Lkkn;
	}
      }
    }
    else{
      pallocn[ln]=1.0;	
    }

    /* --Section 9.6 - Work out acceptance probability  and new state --*/

    lpn=lpost(kn,nkkn,thetan,&llhn);

    logratio+=(lpn-lp);
    logratio+=(log(pallocn[ln])-log(palloc[l]));
    logratio+=(log(lambda[k][l])-log(lambda[kn][ln]));
    logratio+=(log(detB[kn][ln])-log(detB[k][l])); 

    if(sdrand()<exp(max(-30.0,min(0.0,logratio)))){
      for(j1=0;j1<nkkn;j1++){
	theta[j1]=thetan[j1];
      }
      lp=lpn;
      llh=llhn;
      k=kn;
      nkk=nkkn;
      Lkk=Lkkn;
      nacctd++;
    }

    if(adapt==1){
      if(sweep>nburn){
	for(k1=0;k1<kmax;k1++){
	  if(k1==k){
	    propk[k1]=1.0;
	  }
	  else{
	    propk[k1]=0.0;
	  }
	  pk[k1]+=(gamma*(propk[k1]-pk[k1]));
	}
	for(k1=0;k1<kmax;k1++){
	  if(pk[k1]<pkllim){
	    reinit=1;
	  }
	}
	if(reinit==1){
	  reinit=0;
	  nreinit++;
	  pkllim=1.0/(10.0*nreinit);   
	  for(k1=0;k1<kmax;k1++){
	    pk[k1]=1.0/kmax;
	  }
	}
      }
    }

    if(sweep>nburn){
      (ksummary[k])++;
    }
    /* --- Section 10 - Write variables to files --------- */

    if(sweep>nburn){

      fprintf(fpk,"%d\n",k+1);
      fprintf(fplp,"%lf %lf\n",lp,llh);
      for(k1=0;k1<kmax;k1++){
        fprintf(fpp,"%lf ",pk[k1]);
      }
      fprintf(fpp,"\n");
      for(j1=0;j1<nkk;j1++){
        fprintf(fpt[k],"%lf ",theta[j1]);
      }
      fprintf(fpt[k],"\n");
    }
      
    if(sweep>keep&&fmod(sweep-keep,nsokal)<0.005){
      xr[((sweep-keep)/nsokal)-1]=k;
    }
    
    if(sweep==1){
      printf("\nBurning in");
    }
    if((sweep<=nburn)&&(fmod(sweep,(nburn/10))<tol)){
      printf(" .");
      fflush(NULL);
    }
    
    if(sweep==(nburn+1)){
      printf("\nStart of main sample:");
    }
    
    if((sweep>nburn)&&(fmod(sweep-nburn,(nsweep/10))<tol)){
      printf("\nNo. of iterations remaining: %d",nsweep+nburn-sweep);
    }
    fflush(NULL);
  }
  printf("\n");


  /* --- Section 11 - Write log file ----------------------*/

  sokal(nkeep,xr,&var,&tau,&m);
  fprintf(fpl,"\nAutocorrelation Time:\n");
  fprintf(fpl,"nkeep:%d, nsokal:%d, var:%lf, tau:%lf\n",nkeep,nsokal,
	  var,tau);
  for(i1=0;i1<m;i1++){
    fprintf(fpac,"%lf\n",xr[i1]);
  }

  fprintf(fpl,"\nPosterior Model Probabilities:\n");
  for(k1=0;k1<kmax;k1++){
    fprintf(fpl,"Model %d: %lf\n",k1+1,(double)ksummary[k1]/(double)nsweep);
  }

  fprintf(fpl,"\nAcceptance Rates:\n");
  fprintf(fpl,"Block RWM: %lf\n",(double)naccrwmb/(double)ntryrwmb);
  fprintf(fpl,"Single RWM: %lf\n",(double)naccrwms/(double)ntryrwms);
  fprintf(fpl,"Auto RJ: %lf\n",(double)nacctd/(double)ntrytd);
 
  endtime=clock();
  timesecs=(endtime-starttime)/((double)CLOCKS_PER_SEC);
  fprintf(fpl,"\nRun time:\n");
  fprintf(fpl,"Time: %lf\n",timesecs);
     
  return 0;
  
} 
  
 
void gauss(double *z,int n){

  /* Uses Box mueller method to simulate n N(0,1) variables and stores them 
     in z */
  
  int i,n1;
  double u,v;

  n1=n-1;
  for(i=0;i<n1;i+=2){
    u=sqrt(-2.0*log(sdrand()));
    v=tpi*sdrand();
    z[i]=u*sin(v);
    z[i+1]=u*cos(v);
  }
  if(fmod(n,2)<0.5){
    return;
  }
  else{
    u=sqrt(-2.0*log(sdrand()));
    v=tpi*sdrand();
    z[n-1]=u*sin(v);
  }
  return;
}

void rt(double *z,int n,int dof){

  /* Simulates n random t variable with dof degrees of freedom
     by simulating standard normals and chi-squared random variables. 
     Chi-squared rvs simulated by rgamma function that simulates random gamma
     variables (see gammafns file for details)*/

  int j1;
  double s;

  gauss(z,n);
  s=0.5*dof;
  for(j1=0;j1<n;j1++){
    z[j1]/=sqrt(rgamma(s)/s);
  }
  return;
}

void chol(int nkk, double **A){

  /* Performs cholesky decompositon of A and returns result in the 
     same matrix - adapted from PJG Fortran function*/
    
  int j1,j2,j3;
  double sum;

  for(j1=0;j1<nkk;j1++){
    sum=A[j1][j1];
    for(j2=0;j2<j1;j2++){
      sum-=pow(A[j1][j2],2);
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
  
void perm(double *work, int nkk){

  /* Randomly permutes the nkk-vector work */ 
  int j1,j2;
  double temp;

  for(j1=0;j1<(nkk-1);j1++){
    j2=j1+(int)((nkk-j1)*sdrand());  
    if(j2!=j1){
      temp=work[j2];
      work[j2]=work[j1];
      work[j1]=temp;
    }
  }
  return;
}
    
double ltprob(int dof,double z,double *constt){

  /* Evaluates the log of p.d.f. of a t variable with dof degrees of freedom
     at point z */

  double out;

  /* only calculate const of proportionality once */
  if((*constt)>10000.0){
    *constt=loggamma(0.5*(dof+1))-loggamma(0.5*dof)-0.5*log(dof*pi);
  }
  out=(*constt)-0.5*(dof+1)*log(1.0+pow(z,2.0)/dof);
  return out;
}
  

double lnormprob(int k,int nkk,int l,double ***mu, 
		double ****B, double *datai){

  /* Evaluates log of p.d.f. for a multivariate normal for model
     k, of dimension nkk, component l. The summary of means and 
     sqrt of cov matrices (for all models and all component)
     are supplied in mu and B */ 

  int j1,j2;
  double work[nkk];
  double out;
  
  for(j1=0;j1<nkk;j1++){
    work[j1]=datai[j1]-mu[k][l][j1];
  }
  for(j1=0;j1<nkk;j1++){
    for(j2=0;j2<j1;j2++){
      (work[j1])-=B[k][l][j1][j2]*work[j2];
    }
    (work[j1])/=B[k][l][j1][j1]; 
  }
  out=0.0;
  for(j1=0;j1<nkk;j1++){
    out+=(work[j1]*work[j1]);
  }
  out=-0.5*out-(nkk/2.0)*log(tpi)-log(det(k,nkk,l,B));    
  return out;

}

double det(int k, int nkk, int l, double ****B){

  /* Evaluates the determinant of a matrix in B corresponding to model k, 
     component l. */ 
  int j1;
  double out;
  out=1.0;
  for(j1=0;j1<nkk;j1++){
    out*=B[k][l][j1][j1];
  }
  return out;
}










