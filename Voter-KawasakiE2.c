#include <math.h>
#include <stdio.h>
#define NRANSI
#include <time.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
main()
{ 
 FILE *file3; 
 file3=fopen("VK100-case1.d","w");
 int Ln=100; // Side of the square lattice
 int sig[101][101]; // 2D  array denoting lattice state. At any site (i,j) variable value 0 is susceptible, 1 is infected, and 2 is recovered. 
int r1,r2,r3; // Number of susceptible, infected and recovered sites. 
 int pSI,pSR,pIR; // Number of nearest neighbor pairs corresponding to susceptible-infected (pSI), susceptible-recovered (pSR),and infected-recovered (pIR)
 double delta=1.0, gamma=10, lambda=0.5; // Probability rates of transition between states I ->R, R->I, I ->R respectively.
 double alpha=0.5; // Rate of exchange of ANY nearest-neighbor site.
 int tlen =24000; //Number of time steps of the code. 
 long double tstep =0.01; // Duration of a time step. 
 int i,j,i1,k1,k; // Iteration dummy variables. 
 double ran2,ran3,ran4; // Random variables.
 int nbr,m,p1, p2,iu,id,ku,kd;
 int w1,w2,w3,w4; 
gsl_rng *ran=gsl_rng_alloc(gsl_rng_mt19937); //Pointer for generation of random numbers. 
gsl_rng_set(ran,time(NULL));
srand((unsigned)time(NULL));
/* Initial state is assigned where every site is either suscpetible (0) or infected (1) */
for (i=1; i <=Ln;i++)   
{
	for (j=1;j<=Ln; j++)
 {
  ran2=  gsl_rng_uniform(ran);
   if (ran2<0.5)                          
    sig[i][j]=0;
   else 
    sig[i][j]=1;  
  }
 }
 fprintf(file3, "Apha = %lf Gamma= %lf  Lambda = %lf \n\n",alpha, gamma, lambda);
 fprintf (file3,"%s \n", " P(S)   P(I)   P(R)    P(S,I)   P(S,R)   P(I,R)");
/* Main system evolution Loop */
for (j=1;j<tlen -1;j++)
 { 
 for ( i1=1;i1<=Ln;i1++)
  {
  	for (k1=1;k1<=Ln;k1++)
    {
      /*Random selection of a site (i,k) */
      ran2 = gsl_rng_uniform(ran); 
      i =(Ln)*ran2 +1.0;
      ran2 = gsl_rng_uniform(ran);
      k =ran2*(Ln) +1.0;
      /*Periodic Boundary Condition */
      /* iu, id,ku,kd are variables used to determine nearest neighbors of (i,k) */
      if (i==Ln)
       {
     	 iu=1;
         id=Ln-1;
       }
     else if (i==1)
     {
     	id=Ln;
      iu=2;
     }   
     else
     {
     	iu=i+1;
     id=i-1;
     } 
     
     if (k==Ln)
     {
     ku=1;
    kd=Ln-1;
    }
    else  if (k==1)
     {
     	kd=Ln;
        ku=2;
      }
      else
       {
      	 ku=k+1;
         kd=k -1;
       }  
        /*Infected case */
        if (sig[i][k]==1)
	     {
	     ran2 = gsl_rng_uniform(ran); 
            if (ran2 < delta*tstep)  //Determining probability of transition.
             sig[i][k]=2;
       	 }	
       	 /* Recovered case */
	     else if (sig[i][k]==2)
	       { 
	   	     ran2 = gsl_rng_uniform(ran); 
		  if (ran2<gamma*tstep)
		       sig[i][k]=0;  
		  }
		  /*Susceptible case */
         else
         {
		   w1= 2*sig[id][k]-sig[id][k]*sig[id][k];       // this function is 1 if site (id,k) is infected; if not, it is 0.  
  			 w2= 2*sig[iu][k]-sig[iu][k]*sig[iu][k];
		   w4= 2*sig[i][ku]-sig[i][ku]*sig[i][ku];
			    w3= 2*sig[i][kd]-sig[i][kd]*sig[i][kd]; 
		   		  nbr=w1+w2+w4+w3;                   // number of nearest neighbor sites that are infected.
					 ran2 =gsl_rng_uniform(ran);
			  if (ran2 <lambda*nbr*tstep)
			  sig[i][k]=1;
			}
			
                 ran2 = gsl_rng_uniform(ran);
                   i =Ln*ran2+ 1.0;
                ran2 = gsl_rng_uniform(ran);
                 k =ran2*Ln +1.0;
				ran3 = gsl_rng_uniform(ran); 
              /* Random Exchange */
              if (ran3< (alpha)*tstep)
               {
               	p1=sig[i][k];
			  sig[i][k]=sig[iu][k];
			  sig[iu][k]=p1;
			   }  
			ran4 = gsl_rng_uniform(ran); 
			  	if (ran4<alpha*tstep)
               {
               	p1=sig[i][kd];
			  sig[i][kd]=sig[i][k];
			  sig[i][k]=p1;
			   } 
			
		 }
      } 
     
/* Computation of average across the lattice and nearest-neighbor (two point) correlations */
  if(j%2000==0 && j>(tlen/2))  // System is allowed to run half the time (to acheive steady state) before averages are taken 
   {
        k= floor(j*4/tlen);
        r1=r2=r3=0;
        pSI=pSR=pIR=0;
        for (i=1;i<=Ln;i++)
                {
                    for (k=1;k<=Ln;k++)
                        {
                            
                            if (sig[i][k]==1)
                                r2=r2+1;
                            else if (sig[i][k]==0)
                                r1=r1+1;
                            else
                                r3=r3+1;
                            if (k!=Ln)
                              {
                                   if (sig[i][k]==0 && sig[i][k+1]==1)
                                        pSI=pSI+1;
                                   if (sig[i][k]==0 && sig[i][k+1]==2)
                                        pSR =pSR+1;
                                   if (sig[i][k]==1 && sig[i][k+1]==2)
                                        pIR=pIR+1;
                               }
                            else
                               {
                                   if (sig[i][k]==0 && sig[i][1]==1)
                                      pSI=pSI+1;
                                   if (sig[i][k]==0 && sig[i][1]==2)
                                      pSR =pSR+1;
                                   if (sig[i][k]==1 && sig[i][1]==2)
                                      pIR=pIR+1;
                                }
                            }
                   }
                fprintf(file3,"%lf %lf %lf %lf %lf %lf\n", (r1*1.0)/(Ln*Ln),(r2*1.0)/(Ln*Ln),(r3*1.0)/(Ln*Ln),(pSI*1.0)/(Ln*Ln),(pSR*1.0)/(Ln*Ln),(pIR*1.0)/(Ln*Ln));
     }

    }
  fclose(file3);
 }



   
      




