#include <math.h>
#include <stdio.h>
#define NRANSI
#include <time.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
double invsq(double l_av,double r, double ren);
 gsl_rng *ran;   
const double norm_lam_min=0.1, norm_lam_max=3;
main()
{  
 ran = gsl_rng_alloc(gsl_rng_mt19937); //Pointer for generation of random numbers.
 gsl_rng_set(ran,time(NULL));
const double renorm = (1.0/norm_lam_min -1.0/norm_lam_max)/(log(norm_lam_max/norm_lam_min));
printf("%lf \n",renorm);
int nmc=5;
 double lam,delta=1.0, gamma=0.5, lam_av=0.6; // Probability rates of transition between states I ->R, R->I, I ->R respectively.
// srand((unsigned)time(NULL));
 FILE *file3; 
 file3=fopen("Epidem-VR","w");// File to save the data.
 fprintf(file3, "The distribution of the infection rate lambda is chosen from an inverse square distribution between %lf and %lf with mean %lf \n",norm_lam_min*renorm,norm_lam_max*renorm, lam_av*renorm);
 int Ln=100; // Side of the square lattice
 int sig[101][101]; // 2D  array denoting lattice state. At any site (i,j) variable value 0 is susceptible, 1 is infected, and 2 is recovered. 
 int r1[11],r2[11],r3[11]; // Number of susceptible, infected and recovered sites. 
 int pSI[11],pSR[11],pIR[11]; // Number of nearest neighbor pairs corresponding to susceptible-infected (pSI), susceptible-recovered (pSR),and infected-recovered (pIR)
 double sum=0;	
 int count=0;
 double alpha=0.0; // Rate of exchange of ANY nearest-neighbor site.
 int tlen =80000; //Number of time steps of the code. 
 int tlen_sam=tlen/4;
 int sam_int=tlen_sam/10;
 int sam_ind=0;
 long double tstep =0.1; // Duration of a time step. 
 int rep,i,j,i1,k1,k; // Iteration dummy variables. 
 double ran2,ran3,ran4; // Random variables.
 int nbr,m,p1, p2,iu,id,ku,kd;
 int w1,w2,w3,w4; 
 int i_comb,k_comb;
/* Initial state is assigned where every site is either suscpetible (0) or infected (1) */
for (rep=1;rep<=nmc;rep++)
{
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
count=0;
sum=0;
 for (i=0;i<len(r1);i++)
   {
        r1[i]=r2[i]=r3[i]=0;
        pSI[i]=pSR[i]=pIR[i]=0;
    }
 fprintf(file3, "\n \nApha = %lf Gamma= %lf  Lambda = %lf \n\n",alpha, gamma, lam_av);
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

      else
        {
              i=i_comb -Ln;
              k=k_comb;
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
	           lam=invsq(lam_av,ran3,renorm); 
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
  if(j%sam_int==0 && j>(t_len -tlen_sam))  // System is allowed to run half the time (to acheive steady state) before averages are taken 
   {
       sam_ind=floor((j-t_len)/sam_int)
        for (i=1;i<=Ln;i++)
                {
                    for (k=1;k<=Ln;k++)
                        {
                            
                            if (sig[i][k]==1)
                                r2[sam_ind]=r2[sam_ind]+1;
                            else if (sig[i][k]==0)
                                r1[sam_ind]=r1[sam_ind]+1;
                            else
                                r3=r3+1;
                            if (k!=Ln)
                              {
                                   if (sig[i][k]==0 && sig[i][k+1]==1)
                                        pSI[sam_ind]=pSI[sam_ind]+1;
                                   if (sig[i][k]==0 && sig[i][k+1]==2)
                                        pSR[sam_ind] =pSR[sam_ind]+1;
                                   if (sig[i][k]==1 && sig[i][k+1]==2)
                                        pIR[sam_ind]=pIR[sam_ind]+1;
                               }
                            else
                               {
                                   if (sig[i][k]==0 && sig[i][1]==1)
                                      pSI[sam_ind]=pSI[sam_ind]+1;
                                   if (sig[i][k]==0 && sig[i][1]==2)
                                      pSR[sam_ind] =pSR[sam_ind]+1;
                                   if (sig[i][k]==1 && sig[i][1]==2)
                                      pIR[sam_ind]=pIR[sam_ind]+1;
                                }
                            }
                   }
       }
 
  }
      for (i=1;i<11;i++)
      {
         r1_av=r1_av+r1[i];
         r2_av=r2_av+r2[i];
         r3_av=r3_av+r3[i];
         pSI_av=pSI_av + pSI[i];
         pSR_av=pSRav + pSR[i];
         pIR_av=pIR_av + pIR[i];
      } 
        fprintf(file3,"%lf %lf %lf %lf %lf %lf\n", (r1[rep]*1.0)/(Ln*Ln),(r2*1.0)/(Ln*Ln),(r3*1.0)/(Ln*Ln),(pSI*1.0)/(Ln*Ln),(pSR*1.0)/(Ln*Ln),(pIR*1.0)/(Ln*Ln));

    
}
  fclose(file3);
 printf("%lf \n",sum/count*1.0);
 }


double invsq(double l_av,double r,double ren)
{
	double lam_min=l_av*ren*norm_lam_min;
	double lam_max=l_av*ren*norm_lam_max; 
	double lam_const=1.0/(1.0/lam_min - 1.0/lam_max);
        double a = lam_const/lam_max;
        double y = a+r;
        return (lam_const)/y ;
}
   
      




