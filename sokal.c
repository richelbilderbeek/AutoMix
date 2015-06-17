/* Functions to estimates integrated autocorrelation time using 
   method of Sokal. Taken from PJG function sokal.f
   Note that the definition is the sum from minus infinity to 
   infinity of the autocorrelation function, hence twice Sokal's 
   definition. */

#include <math.h>
#include <stdio.h>

#define pi 3.141592653589793

void fastfr(int nin, double *xreal, double *ximag);

void sokal(int n, double *xreal, double *var, double *tau, int *m){

  int i1;
  double ximag[n],c,sum;

  if(n>pow(2.0,20.0)){
    printf("\nAuto-correlation length exceeded");
    return;
  }

  for(i1=0;i1<n;i1++){
    ximag[i1]=0.0;
  }

  /* Use FFT to compute autocorrelations (in xr)  and variance (in var). */

  fastfr(n,xreal,ximag);

  for(i1=0;i1<n;i1++){
    xreal[i1]=pow(xreal[i1],2.0)+pow(ximag[i1],2.0);
    ximag[i1]=0.0;
  }

  xreal[0]=0.0;
  
  fastfr(n,xreal,ximag);
  *var=xreal[0]/((double)n*(n-1));
  c=1.0/xreal[0];

  for(i1=0;i1<n;i1++){
    xreal[i1]=xreal[i1]*c;
  }

  /* Use Sokal's adaptive truncated periodogram method to estimate 
     integrated autocorrelation time (in tau). */
  
  sum=-(1.0/3.0);
  for(i1=0;i1<n;i1++){
    sum+=xreal[i1]-(1.0/6.0);
    if(sum<0.0){
      goto CALC_TAU;
    }
  }
  
 CALC_TAU:
  *tau=2.0*(sum+i1/6.0);
  *m=i1+1;

  return;
}

void fastfr(int nin, double *xreal, double *ximag){

  /* radix 4 complex discrete fast Fourier transform
     without usual normalisation
     Eric Renshaw -> PJG 20 March 1987 -> DIH, March 2004

     xreal = array which on input contains real part of data
     for transformation and on output gives real part of the
     result,type real, dimension iabs(isize) or greater.
     ximag = array which on input contains imaginary part
     for transformation and on output gives imaginary part of
     the result, type real, dimension iabs(isize) or greater.
     isize = integer variable or constant specifying length and
     type of transform. The length is iabs(isize) which must be
     a power of 2, minimum length 4, maximum 2**20. Illegal length
     leads to warning and return. If isize is positive the forward 
     transform is calculated. If negative the inverse transform is found.

     The transform is defined by
     out(r) = sum(j = 0 to n-1) in(j)*exp(-2*pi*i*r*j/n)
			if isize = n > 0,
     out(r) = sum(j = 0 to n-1) in(j)*exp(2*pi*i*r*j/n)
			if isize = -n < 0,
     for r = 0,1,2,...,(n-1),
       where i = sqrt(-1), and both in(j) and out(j)
       are stored in xreal(j+1)+i*ximag(j+1)              */

  
  double z,bcos,bsin,temp,cw1,cw2,cw3,sw1,sw2,sw3;
  double xs0,xs1,xs2,xs3,ys0,ys1,ys2,ys3,x1,x2,x3,y1,y2,y3;

  int i0,i1,i2,i3,ul[20],n,counta,countb,time,indic,test;
  int j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,j19;

  n = nin<0 ? -nin : nin;

  if(n<4){
    printf("\nFast Fourier length too short");
    return;
  }
  
  /* if this is to be an inverse transform, conjugate the data */

  if(nin<0){
    for(i1=0;i1<n;i1++){
      ximag[i1]=-ximag[i1];
    }
  }

  test=n;
  while(test>1){
    if(test&1){
      printf("\nFast Fourier length must be power of 2");
      return;
    }
    test>>=1;
  }
    
  counta=n/4;
  time=0;
  
  while(counta>0){
    countb=counta*4;
    time+=2;

    /* do the transforms required by this stage */

    z=pi/countb;
    bcos=-2.0*pow(sin(z),2.0);
    bsin=sin(2.0*z);
    
    cw1=1.0;
    sw1=0.0;

    /* this is the main calculation of radix 4 transforms */

    for(j1=0;j1<counta;j1++){
      for(j2=j1;j2<n;j2+=countb){
	i0=j2;
	i1=i0+counta;
	i2=i1+counta;
	i3=i2+counta;
	xs0=xreal[i0]+xreal[i2];
	xs1=xreal[i0]-xreal[i2];
	ys0=ximag[i0]+ximag[i2];
	ys1=ximag[i0]-ximag[i2];
	xs2=xreal[i1]+xreal[i3];
	xs3=xreal[i1]-xreal[i3];
	ys2=ximag[i1]+ximag[i3];
	ys3=ximag[i1]-ximag[i3];

	xreal[i0]=xs0+xs2;
	ximag[i0]=ys0+ys2;
	
	x1=xs1+ys3;
	y1=ys1-xs3;
	x2=xs0-xs2;
	y2=ys0-ys2;
	x3=xs1-ys3;
	y3=ys1+xs3;

	if(j1==0){
	  xreal[i2]=x1;
	  ximag[i2]=y1;
	  xreal[i1]=x2;
	  ximag[i1]=y2;
	  xreal[i3]=x3;
	  ximag[i3]=y3;
	}
	else{

	  /* multiply by twiddle factors if required */

	  xreal[i2]=x1*cw1+y1*sw1;
	  ximag[i2]=y1*cw1-x1*sw1;
	  xreal[i1]=x2*cw2+y2*sw2;
	  ximag[i1]=y2*cw2-x2*sw2;
	  xreal[i3]=x3*cw3+y3*sw3;
	  ximag[i3]=y3*cw3-x3*sw3;
	}
      }
      if(j1<(counta-1)){

	/* calculate a new set of twiddle factors */
	
	z=cw1*bcos-sw1*bsin+cw1;
	sw1=bcos*sw1+bsin*cw1+sw1;
	temp=1.5-0.5*(z*z+sw1*sw1);
	cw1=z*temp;
	sw1=sw1*temp;
	cw2=cw1*cw1-sw1*sw1;
	sw2=2.0*cw1*sw1;
	cw3=cw1*cw2-sw1*sw2;
	sw3=cw1*sw2+cw2*sw1;
      }
    }

    indic=0;
    if(counta>1){
      /* set up the transform split for the next stage */
      counta/=4;
      indic=1;
    }
    else{
      counta=0;
    }
  }
  if(indic){
   
    /* this is the calculation of a radix two stage */
    for(j1=0;j1<n;j1+=2){
      temp=xreal[j1]+xreal[j1+1];
      xreal[j1+1]=xreal[j1]-xreal[j1+1];
      xreal[j1]=temp;
      temp=ximag[j1]+ximag[j1+1];
      ximag[j1+1]=ximag[j1]-ximag[j1+1];
      ximag[j1]=temp;
    }
    time++;
  }
  
  /* if this was an inverse transform, conjugate the result */

  if(nin<0){
    for(j1=0;j1<n;j1++){
      ximag[j1]=-ximag[j1];
    }
  }

  /* unscramble the result */

  if(time>20){
    printf("\nFast Fourier length too long");
    return;
  }

  i1=20-time;
  for(j1=0;j1<i1;j1++){
    ul[j1]=1;
  }
  if(i1==0){
    ul[0]=2;
    i1++;
  }
  for(j1=i1;j1<20;j1++){
    ul[j1]=2*ul[j1-1];
  }

  i0=0;
  for(j0=0;j0<ul[0];j0++){
    for(j1=j0;j1<ul[1];j1+=ul[0]){
      for(j2=j1;j2<ul[2];j2+=ul[1]){
	for(j3=j2;j3<ul[3];j3+=ul[2]){
	  for(j4=j3;j4<ul[4];j4+=ul[3]){
	    for(j5=j4;j5<ul[5];j5+=ul[4]){
	      for(j6=j5;j6<ul[6];j6+=ul[5]){
		for(j7=j6;j7<ul[7];j7+=ul[6]){
		  for(j8=j7;j8<ul[8];j8+=ul[7]){
		    for(j9=j8;j9<ul[9];j9+=ul[8]){
		      for(j10=j9;j10<ul[10];j10+=ul[9]){
			for(j11=j10;j11<ul[11];j11+=ul[10]){
			  for(j12=j11;j12<ul[12];j12+=ul[11]){
			    for(j13=j12;j13<ul[13];j13+=ul[12]){
			      for(j14=j13;j14<ul[14];j14+=ul[13]){
				for(j15=j14;j15<ul[15];j15+=ul[14]){
				  for(j16=j15;j16<ul[16];j16+=ul[15]){
				    for(j17=j16;j17<ul[17];j17+=ul[16]){
				      for(j18=j17;j18<ul[18];j18+=ul[17]){
					for(j19=j18;j19<ul[19];j19+=ul[18]){
					  if((i0-j19)<0){
					    temp=xreal[i0];
					    xreal[i0]=xreal[j19];
					    xreal[j19]=temp;
					    temp=ximag[i0];
					    ximag[i0]=ximag[j19];
					    ximag[j19]=temp;
					  }
					  i0++;
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return;
}
