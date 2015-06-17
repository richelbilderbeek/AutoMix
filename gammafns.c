/* Routines calculating quantities related to the gamma function */

#include<math.h>

#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))

extern double sdrand();

/* Taken from algama.f (PJG) - converted to C using appropriate machine
   constants by DIH 04/11/03 */

double loggamma(double X){

  /* This routine calculates the LOG GAMMA function for a positive real
     argument X.  Computation is based on an algorithm outlined in
     references 1 and 2.  The program uses rational functions that
     theoretically approximate LOG GAMMA to at least 18 significant
     decimal digits.  The approximation for X > 12 is from reference
     3, while approximations for X < 12.0 are similar to those in
     reference 1, but are unpublished.  The accuracy achieved depends
     on the arithmetic system, the compiler, the intrinsic functions,
     and proper selection of the machine-dependent constants.   
    
     Values taken from float.h and algama.f (for XBIG)

     ---------------------------------------------------------------

     Explanation of machine-dependent constants

  beta   - radix for the floating-point representation
  maxexp - the smallest positive power of beta that overflows
  XBIG   - largest argument for which LN(GAMMA(X)) is representable
           in the machine, i.e., the solution to the equation
	           LN(GAMMA(XBIG)) = beta**maxexp
  XINF   - largest machine representable floating-point number;
           approximately beta**maxexp.
  EPS    - The smallest positive floating-point number such that
           1.0+EPS > 1.0
  FRTBIG - Rough estimate of the fourth root of XBIG  

     ---------------------------------------------------------------

     Error returns

     The program returns the value XINF for X <= 0.0 or when
     overflow would occur.  The computation is believed to 
     be free of underflow and overflow.

     ---------------------------------------------------------------
     References:

     1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
     1967, pp. 198-203.
     
     2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
     1969.
     
     3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
     York, 1968.
     
  -----------------------------------------------------------------
  Start of code 
  -----------------------------------------------------------------*/

  int I;
  double CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT68;
  double RES,SQRTPI,THRHAL,TWELVE,TWO,XBIG,XDEN,XINF;
  double XM1,XM2,XM4,XNUM,Y,YSQ,ZERO;
  double C[7],P1[8],P2[8],P4[8],Q1[8],Q2[8],Q4[8];
 
  /*---------------------------------------------------------------------
  Mathematical constants
  ---------------------------------------------------------------------*/

  ONE=1.0E0;
  HALF=0.5E0;
  TWELVE=12.0E0;
  ZERO=0.0E0;
  FOUR=4.0E0;
  THRHAL=1.5E0;
  TWO=2.0E0;
  PNT68=0.6796875E0;
  SQRTPI=0.9189385332046727417803297E0;  /* eh? */

  /*---------------------------------------------------------------------
    Machine dependent parameters
    -------------------------------------------------------------------*/

  XBIG=2.55E305;
  XINF=1.79E308;
  EPS=2.22E-16;
  FRTBIG=2.25E76;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    approximation over (EPS,1.5).
    -------------------------------------------------------------------*/
  D1=-5.772156649015328605195174E-1;
  P1[0]=4.945235359296727046734888E0;
  P1[1]=2.018112620856775083915565E2;
  P1[2]=2.290838373831346393026739E3;
  P1[3]=1.131967205903380828685045E4;
  P1[4]=2.855724635671635335736389E4;
  P1[5]=3.848496228443793359990269E4;
  P1[6]=2.637748787624195437963534E4;
  P1[7]=7.225813979700288197698961E3;
  Q1[0]=6.748212550303777196073036E1;
  Q1[1]=1.113332393857199323513008E3;
  Q1[2]=7.738757056935398733233834E3;
  Q1[3]=2.763987074403340708898585E4;
  Q1[4]=5.499310206226157329794414E4;
  Q1[5]=6.161122180066002127833352E4;
  Q1[6]=3.635127591501940507276287E4;
  Q1[7]=8.785536302431013170870835E3;
 
 /*---------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    Approximation over (1.5,4.0).
    ------------------------------------------------------------------*/
  
  D2=4.227843350984671393993777E-1;
  P2[0]=4.974607845568932035012064E0;
  P2[1]=5.424138599891070494101986E2;
  P2[2]=1.550693864978364947665077E4;
  P2[3]=1.847932904445632425417223E5;
  P2[4]=1.088204769468828767498470E6;
  P2[5]=3.338152967987029735917223E6;
  P2[6]=5.106661678927352456275255E6;
  P2[7]=3.074109054850539556250927E6;
  Q2[0]=1.830328399370592604055942E2;
  Q2[1]=7.765049321445005871323047E3;
  Q2[2]=1.331903827966074194402448E5;
  Q2[3]=1.136705821321969608938755E6;
  Q2[4]=5.267964117437946917577538E6;
  Q2[5]=1.346701454311101692290052E7;
  Q2[6]=1.782736530353274213975932E7;
  Q2[7]=9.533095591844353613395747E6;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    Approximation over (4.0,12.0).
    -------------------------------------------------------------------*/

  D4=1.791759469228055000094023E0;
  P4[0]=1.474502166059939948905062E4;
  P4[1]=2.426813369486704502836312E6;
  P4[2]=1.214755574045093227939592E8;
  P4[3]=2.663432449630976949898078E9;
  P4[4]=2.940378956634553899906876E10;
  P4[5]=1.702665737765398868392998E11;
  P4[6]=4.926125793377430887588120E11;
  P4[7]=5.606251856223951465078242E11;
  Q4[0]=2.690530175870899333379843E3;
  Q4[1]=6.393885654300092398984238E5;
  Q4[2]=4.135599930241388052042842E7;
  Q4[3]=1.120872109616147941376570E9;
  Q4[4]=1.488613728678813811542398E10;
  Q4[5]=1.016803586272438228077304E11;
  Q4[6]=3.417476345507377132798597E11;
  Q4[7]=4.463158187419713286462081E11;
  
  /*---------------------------------------------------------------------
    Coefficients for minimax approximation over (12, INF).
    -------------------------------------------------------------------*/
  C[0]=-1.910444077728E-03;
  C[1]=8.4171387781295E-04;
  C[2]=-5.952379913043012E-04;
  C[3]=7.93650793500350248E-04;
  C[4]=-2.777777777777681622553E-03;
  C[5]=8.333333333333333331554247E-02;
  C[6]=5.7083835261E-03;

  /*----------------------------------------------------------------------
    0 < X <= EPS 
    --------------------------------------------------------------------*/
  Y=X;
  if((Y>0)&&(Y<=XBIG)){
    if(Y<=EPS){
      RES=-log(Y);
    }
    else if(Y<=THRHAL){
 /*-----------------------------------------------------------------------
   EPS < X <= 1.5
   ---------------------------------------------------------------------*/
      if(Y<PNT68){
	CORR=-log(Y);
	XM1=Y;
      }
      else{
	CORR=ZERO;
	XM1=(Y-HALF)-HALF;
      }

      if((Y<=HALF)||(Y>=PNT68)){
	XDEN=ONE;
	XNUM=ZERO;
	for(I=0;I<8;I++){
	  XNUM=XNUM*XM1+P1[I];
	  XDEN=XDEN*XM1+Q1[I];
	}
	RES=CORR+(XM1*(D1+XM1*(XNUM/XDEN)));
      }
      else{
	XM2=(Y-HALF)-HALF; /*Is XM2 symbol used to agree with other 2 symbols*/
	XDEN=ONE;
	XNUM=ZERO;
	for(I=0;I<8;I++){
	  XNUM=XNUM*XM2+P2[I];
	  XDEN=XDEN*XM2+Q2[I];
	}
	RES=CORR+XM2*(D2+XM2*(XNUM/XDEN));
      }
    }
    else if(Y<=FOUR){

  /*---------------------------------------------------------------------
    1.5 < X <= 4.0
    -------------------------------------------------------------------*/
      XM2=Y-TWO;
      XDEN=ONE;
      XNUM=ZERO;
      for(I=0;I<8;I++){
	XNUM=XNUM*XM2+P2[I];
	XDEN=XDEN*XM2+Q2[I];
      }
      RES=XM2*(D2+XM2*(XNUM/XDEN));
    }
    else if(Y<=TWELVE){

  /*----------------------------------------------------------------------
    4.0 < X <= 12.0
    --------------------------------------------------------------------*/
      XM4=Y-FOUR;
      XDEN=-ONE;
      XNUM=ZERO;
      for(I=0;I<8;I++){
	XNUM=XNUM*XM4+P4[I];
	XDEN=XDEN*XM4+Q4[I];
      }
      RES=D4+XM4*(XNUM/XDEN);
    }
    else{
              
  /*---------------------------------------------------------------------
    X > 12.0,
    -------------------------------------------------------------------*/
      RES=ZERO;
      if(Y<=FRTBIG){
	RES=C[6];
	YSQ=Y*Y;
	for(I=0;I<6;I++){
	  RES=RES/YSQ+C[I];
	}
      }
      RES=RES/Y;
      CORR=log(Y);
      RES=RES+SQRTPI-HALF*CORR;
      RES=RES+Y*(CORR-ONE);
    }
  }
  else{

  /*----------------------------------------------------------------------
    Return for bad arguments
    --------------------------------------------------------------------*/
    RES=XINF;
  }

  /*----------------------------------------------------------------------
    Final return
    --------------------------------------------------------------------*/
  return(RES);

}


/* Function for generating a Gamma(s,1) random variable 
   Converted to C from PJG rgamma FORTRAN function
   Note: (1/t)*GAMMA(s,1) is GAMMA(s,t) */

double rgamma(double s){

  double exp1,b,c1,c2,c3,c4,c5,u1,u2,w,bu,out;

  exp1=exp(1.0);

  if(s<1){
    b=(s+exp1)/exp1;
    c1 = 1.0/s;
  LAB_1:
    bu=b*sdrand();
    if(bu<=1.0){
      out=exp(max(-30.0,c1*log(bu)));
      if(sdrand()>=exp(-out)){
	goto LAB_1;
      }
    }
    else{
      out=-log((b-bu)/s);
      if(sdrand()>=pow(out,(s-1.0))){
	goto LAB_1;
      }
    }
  }
  else if(s==1.0){
    out=-log(sdrand());
  }
  else{
    c1=s-1.0;
    c2=(s-1.0/(6.0*s))/c1;
    c3=2.0/c1;
    c4=c3+2.0;
    c5=1.0/sqrt(s);
  LAB_2:
    u1=sdrand();
    u2=sdrand();
    if(s>2.5){
      u1=u2+c5*(1.0-1.86*u1);
    }
    if(u1<=0.0||u1>=1.0){
      goto LAB_2;
    }
    w=c2*u2/u1;
    if((c3*u1+w+1.0/w)<=c4){
      goto LAB_3;
    }
    if((c3*log(u1)-log(w)+w)>=1.0){
      goto LAB_2;
    }
  LAB_3:	
    out=c1*w;
  }

  return out;

}     
