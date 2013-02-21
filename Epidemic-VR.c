#include <math.h>
#include <stdio.h>
#define NRANSI
#include <time.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
float invsq(float l_av,double r);
main()
{ 
 FILE *file3; 
 file3=fopen("Epidem-VR-1","w");// File to save the data.
 int Ln=100; // Side of the square lattice
 int sig[101][101]; // 2D  array denoting lattice state. At any site (i,j) variable value 0 is susceptible, 1 is infected, and 2 is recovered. 
 int r1,r2,r3; // Number of susceptible, infected and recovered sites. 
 int pSI,pSR,pIR; // Number of nearest neighbor pairs corresponding to susceptible-infected (pSI), susceptible-recovered (pSR),and infected-recovered (pIR)
 double sum=0;	
 int count=0;
 double lam,delta=1.0, gamma=0.5, lam_av=0.6; // Probability rates of transition between states I ->R, R->I, I ->R respectively.
 double alpha=0.0; // Rate of exchange of ANY nearest-neighbor site.
 int tlen =40000; //Number of time steps of the code. 
 long double tstep =0.1; // Duration of a time step. 
 int i,j,i1,k1,k; // Iteration dummy variables. 
 gsl_rng *ran=gsl_rng_alloc(gsl_rng_mt19937); //Pointer for generation of random numbers. 
 gsl_rng_set(ran,time(NULL));
 srand((unsigned)time(NULL));
 double ran2,ran3,ran4; // Random variables.
 int nbr,m,p1, p2,iu,id,ku,kd;
 int w1,w2,w3,w4; 
 int i_comb,k_comb;
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
 fprintf(file3, "Apha = %lf Gamma= %lf  Lambda = %lf \n\n",alpha, gamma, lam_av);
 fprintf (file3,"%s \n", " P(S)   P(I)   P(R)    P(S,I)   P(S,R)   P(I,R)");
/* Main system evolution Loop */
for (j=1;j<tlen -1;j++)
 { 
 for ( i1=1;i1<=Ln;i1++)
  {
  	for (k1=1;k1<=Ln;k1++)
    {
      /*Random selection of a site (i_comb,k_comb) on the lattice (2*Ln, Ln)*/
      ran2 = gsl_rng_uniform(ran); 
      i_comb =(2*Ln)*ran2 +1.0; 
      ran2 = gsl_rng_uniform(ran);
      k_comb =ran2*(Ln) +1.0;
      if (i_comb>Ln) 
        {
          if (sig[i_comb -Ln][k_comb]==1)
             {
              i=i_comb -Ln;
              k=k_comb;
             }
          else
            break;
        }

        
        
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
             if (ran2<0.5)
              {                
	        ran2 = gsl_rng_uniform(ran); 
                if (ran2 < delta*tstep)  //Determining probability of transition.
                sig[i][k]=2;
       	      }
             else
              {
	           ran3 = gsl_rng_uniform(ran); 
	           lam=invsq(lam_av,ran3); 
                   //lam=lam_av;
                   //printf("%lf \n",lam);
                   count+=1;
                   sum+=lam;
	           ran3 = gsl_rng_uniform(ran); 
                   if (ran3<lam*tstep)
                   {
                        if (sig[id][k]==0)
                          sig[id][k]=1;    
                        if (sig[iu][k]==0)
                          sig[iu][k]=1;    
                        if (sig[i][ku]==0)
                          sig[i][ku]=1;    
                        if (sig[i][kd]==0)
                          sig[i][kd]=1;    
                  }
              }
 
         }	
       	 /* Recovered case */
	     else if (sig[i][k]==2)
	       { 
	   	     ran2 = gsl_rng_uniform(ran); 
	             if (ran2<gamma*tstep)
		       sig[i][k]=0;  
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
 printf("%lf \n",sum/count*1.0);
 }


float invsq(float l_av,double r)
{
	double norm_lam_min=0.1, norm_lam_max=4;
	double lam_min=l_av*norm_lam_min;
	double lam_max=l_av*norm_lam_max; 
	double lam_const=1.0/(1.0/lam_min - 1.0/lam_max);
        double renorm = (1.0/norm_lam_min -1.0/norm_lam_max)/(log(norm_lam_max/norm_lam_min));
        double a = lam_const/lam_max;
        double y = a+r;
        return (renorm*lam_const)/y ;
}
   
      




