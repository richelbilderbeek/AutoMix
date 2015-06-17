##########################################################################
# README file - The AutoMix sampler
#  
# Last edited 26/09/06.
# Sampler developed by David Hastie, University of Bristol, UK as a part of a 
# submission for the degree of Ph.D. This Ph.D. was supervised by 
# Prof. Peter Green (PJG), University of Bristol. Special thanks also 
# to Dr. Christophe Andrieu (CA), University of Bristol, for advice on 
# adaptive schemes and mixture fitting. 
#
# The AutoMix sampler is free for personal and academic use, but users must 
# reference the sampler as instructed below.  For commercial
# use permission must be sought from the author. To seek permission
# for such use please send an e-mail to d_hastie@hotmail.com 
# outlining the desired usage.  
#
# Use of the AutoMix sampler is entirely at the user's own risk. It is the
# responsibility of the user to ensure that any conclusions made through the 
# use of the AutoMix sampler are valid. The author accepts no responsibility 
# whatsoever for any loss, financial or otherwise, that may arise in 
# connection with the use of the sampler.   
#
# The AutoMix sampler is available from http://www.davidhastie.me.uk/AutoMix
# The sampler may be modified and redistributed as desired but the author
# encourages users to register at the above site so that notice can be 
# received of updates of the software.
#
# Before use, please read this README file bundled with this software.   
#
# 26/09/06 UPDATE
# Bug fixed in change point user files. Comment had been added 
# without closing the comment, meaning function loggamma was not defined.
############################################################################

##### INSTALLATION of the AutoMix sampler #####

The AutoMix package is a C program for Unix-like systems, implementing 
the automatic reversible jump MCMC sampler of the same name
described in Chapters 4, 5, and 6 of David Hastie's Ph.D. thesis

The package distributed consists of:
      main C source file	 automix.c
      auxiliary routines in C	 gammafns.c sd.c sokal.c
      makefile			 Makefile
      example user files	 usertoy1.c usertoy2.c usercpt.c 
				 usercptrs.c userrb9.c userddi.c
      header file ddi example    ddidata.h 
      readme file		 README.txt

The package is provided as both a gzipped tarfile and as a compressed 
windows archive. For a free Unix-like shell that runs under Windows, 
we suggest Cygwin (http:\\www.cygwin.com).

The program has been run successfully:
      using GNU compilers and Cygwin under Windows on Intel PCs
      using GNU compilers on Sun Workstations

We believe that with the appropriate changes to the Makefile the 
program will run as required on other platforms. 

##### USING the AutoMix sampler #####

The AutoMix package is a C program for Unix-like systems, implementing 
the automatic reversible jump MCMC sampler of the same name
described in Chapters 4, 5, and 6 of David Hastie's Ph.D. thesis
These notes assume familiarity with this material.

In particular, potential users should carefully understand the 
limitations of using the AutoMix sampler. The reliability of results 
from this sampler depends on many factors, including the scaling of 
the parameters and the degree of multimodality of the within-model
conditionals of the target distribution. 

To run the sampler for a particular problem the program must be compiled 
with a user-provided file containing four C functions that define the 
problem in question. Examples of such files for problems considered 
within the aforementioned thesis are bundled with the archived AutoMix 
software. More details of the functions are provided below.

The output of the sampler is in the form of several text files 
summarising the MCMC sampler. We discuss these files below. 
Subsequent analysis may be performed on the output with the use of 
statistical packages. We recommend R (http://www.stats.bris.ac.uk/R/) 
as a good free choice of such a package.  

##### WRITING USER FUNCTIONS for the AutoMix sampler #####

For each example to which the user wishes to apply the AutoMix sampler
the user must supply a file (which is linked at compile time, see Makefile 
for examples) containing four user functions written in C. A familiarity 
with the C programming language is assumed.

These functions must have names and arguments detailed below and return 
the following information:

1.  A function to get the number of models:

    void getkmax(int *kmax)
    
    On exit from the function the (integer) number of models under 
    consideration should be contained in *kmax.

2.  A function to get the dimension of each model:

    void getnk(int kmax,int *nk){

    Given the number of models kmax, on exit from the function the 
    vector nk should contain the dimensions of the models. In particular
    the dimension of model k should be in n[k-1], for k=1,...,kmax (Note
    that although models run from 1,...,kmax, the vectors run from
    the element 0 in C ) 

3.  A function to get initial conditions (possibly random) for
    the RWM for a given model:

    void getic(int k, int nkk, double *rwm){

    Given the model index k, and the dimension nkk of that model, 
    on exit from the function the vector rwm should contain the initial 
    conditions. In particular, rwm[j-1] should contain 
    the (possibly random) intial state for component j of the parameter
    vector associated with model k, (j=1,...,nkk).
    
    If random initial states are used, the user must also declare the 
    random number functions in the file (see the example files).   

4.  A function to return the log of the target function pi evaluated at 
    a given point in the state space, up to an additive constant.

    double lpost(int k,int nkk,double *theta,double *llh1){

    Given the model index k, and parameter vector theta (of dimension nkk),
    the function must return the log of the target function (up to 
    an additive constant) evaluated at this point. If pi is a posterior
    distribution, the double *llh1 should contain the likelihood evaluated
    at this point (although this is only necessary for returning the 
    likelihood to output file, and can contain any other value if preferred).

Any other functions used by these four functions should also be declared 
in the user file. They should also be defined if they are not defined in 
other files that are linked at compile time.

The examples provided, with comments, show typical examples of these user 
files for the problems under consideration. 

##### COMPILING the AutoMix sampler #####

To compile the AutoMix sampler for the examples included within
the AutoMix package, the user should edit the Makefile (supplied) so that
the appropriate C compiler is used (currently set to be the GNU sampler 
gcc). Other aspects of the Makefile may be edited if required (for example
the compile flags, or libraries) 

Typing
    make all

in the shell, in the AutoMix folder where the AutoMix distribution was 
unzipped to, will compile all the programs that were distributed with this 
package. By default, for each example two executable programs are created 
the first using optimisation (for faster run times) and the second
using the debugging flags to enable the user to debug the programs. The 
supplied programs have all been debugged but we acknowledge that use
of a debugger can often help to understand how the program works. 

Any of the programs can also be made individually, with the appropriate
make command (for example "make amtoy1" , see the Makefile for further 
examples).

The executables have form 
    amNAME     (for the optimised version)
    amNAMEd    (for the version that can be debugged)

where NAME is the name of the example.

To remove the executables and object (.o) files, type  
    make clean   

To compile the sampler for a new example, the Makefile can be edited to 
contain the appropriate commands (remembering to include command for
compiling any dependencies) for the new program.

##### RUNNING the AutoMix sampler #####

The sampler is run by typing the name of the executable, followed
by run-time flags separated by spaces. The run-time flags control
a number of options within the sampler. If flags are not supplied
default values are used. The flags can be used in any order.

The flags can be summarised as follows (I is assumed to be a positive
integer):

  -mD     controls the mode of the sampler. D=0 is mixture fitting;
	  D=1 skips stage 1 and 2 if a file containing the mixture
	  parameters is supplied; D=2 fits AutoMix version of 
	  AutoRJ sampler (see Green, 2003 - full reference in thesis).
	  (Default uses D=0).  
  -nI     run the sampler for max(I,nkk*10000,100000) iterations
	  in the stage 1 RWM for each model k. (Default uses I=100000)
  -NI	  run the sampler for I Reversible jump iterations in stage 3.
	  (Default uses I=100000).
  -sI     initialises the random number generator with seed I. 
	  (Default uses clock as seed).
  -aA     controls whether or not adaptation is done in stage 3 RJ. If
	  A=0 no adaptation is done, if A=1 adaptation is done. (Default
	  has A=1).
  -pP     controls whether or not random permutation is done in stage 3 
	  RJ. If P=0 no permutation is done, if P=1 permutation is done.
	  (Default has P=0).
  -tI	  Controls whether standard Normal or t distributed variables are 
	  used in RWM and in RJ moves. If I=0 Normal variables are used,
	  otherwise t-distributed variables with I degrees of freedom are 
	  used. (Default I=0). 
  -fF	  Uses the string F as the bases for filenames (e.g. if F=output,
	  filenames are output_log.data, output_mix.data etc). (Default
	  is F=output)

As an example, typing
   amtoy1 -m0 -N1000000 -p1 -ftoy1

runs the optimised mixture fitting version of the toy1 problem 
(see thesis, section 5.5.1) with 1 million RJ sweeps, enabling
permutation and storing the output in files of the type toy1_***.data 	  

Running the sampler produces a summary of how the run is progressing.

For each of the models:
In stage 1 a countdown of the number of iterations remaining is
printed to screen;
In stage 2 a summary of the mixture fitting is printed to screen. 
This summary consists of a countdown of the number of components in 
the current mixture, with the iteration number that the last component 
was removed and an indicator n if the component was annihilated 
naturally, and f if the annihilation was forced.

In the RJ stage 3 a countdown of the number of iterations remaining
is printed to screen. 

No summary statistics are printed to screen. Instead all output from
the sampler is written to files.

##### OUTPUT from the AutoMix sampler #####

The following files are returned by the AutoMix sampler (assuming the
filestem is "output")

  output_ac.data   

  A file containing the estimated autocorrelation coefficient for the 
  chain. The coefficients go up to the order that was used to 
  compute the autocorrelation time using Sokal's method 
  (see Green and Han, 1992)
  
  
  output_adapt.data		   

  A file summarising the AAP adaptation process for each model - adapting
  the scale parameter of the single component RWM. For a model with nk 
  components, there are 2*nk columns, with the odd columns showing the 
  evolution of the scale parameters and the even columns showing the 
  evolution of the average RWM  acceptance rate for the corresponding 
  scale parameter. 


  output_cf.data

  A file summarising the evolution of the mixture fitting for each 
  model. In particular, for each model the first column records
  the number of mixture components at the current iteration, 
  the third column summarises the cost function (to be minimised
  to find the best mixture) and the final column is an indicator
  function which takes the value 1 if a component is annihilated 
  naturally at that iteration, 2 if a component is removed
  by forced annihiliation. (The second column is a term that 
  contibutes to the cost function, see Figueiredo and Jain, 2002). 
  

  output_k.data

  A file recording the model index at each sweep of the 
  RJ sampler (after burn-in).

 
  output_log.data  

  A log file, containing important run statistics, including 
  run-time options, run seed, mixture and RWM parameters for each 
  model, and acceptance rates and autocorrelation statistics for the 
  RJ stage 3. Model probabilities are also returned in this file.
 

  output_lp.data

  A file recording the log posterior (1st column) and log likelihood
  (2nd column) at each sweep of the RJ sampler (after burn-in).


  output_mix.data
 
  A file containing the fitted mixture parameters for the each model. 
  The file is for use by the AutoMix program and is not easily readable
  by the user. The file contains the number of models and the dimension 
  of each model. Then for each model in turn, the file records the 
  adapted random walk scaling parameters, the number of components 
  of the mixture fitted to that model, and for each component in turn 
  the weight, the mean and the lower triangle of the B matrix. It is 
  this file that is required by the AutoMix sampler if it is run in 
  mode 1.


  output_pk.data

  A file containing the evolution of the model jumping proposal 
  parameters for the RJ stage 3. Column k is the probability of
  jumping to model k.


and for each model k...

  output_thetak.data

  A file recording the parameters conditional on model k at each sweep.  
 

##### REFERENCING the AutoMix sampler #####

The following reference should be used:

Hastie, D. I. (2004) Towards Automatic Reversible Jump Markov
Chain Monte Carlo, Ph.D. Thesis, Statistics Group, University 
of Bristol.

It is intended that this reference will change to an academic paper 
 (currently in preparation) in the not too distant future. 
