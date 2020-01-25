#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>
#include"nbody.h"
#include<omp.h>

/*  GLOBAL VARIABLES  */

// cosmological parameters read from input file  "input.nbody_com"

float  vhh, // Hubble parameter in units of 100 km/s/Mpc
  vomegam, // Omega_matter; total matter density (baryons+CDM) parameter
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon
  sigma_8_present ,//  Last updated value of sigma_8 (Presently WMAP)
  vnn; // Spectral index of primordial Power spectrum

long N1,N2,N3;// box dimension (grid) 
int NF, // Fill every NF grid point 
  Nbin; // Number of bins to calculate final P(k) (output)

float   LL; // grid spacing in Mpc

long    MM; // Number of particles
// global variables  (calculated )

int zel_flag=1, // memory allocation for zel is 3 times that for nbody
  fourier_flag; //for fourier transfrom

float  DM_m, // Darm matter mass of simulation particle in 10^10 M_sun h^-1 unit 
  norm, // normalize Pk
  pi=M_PI;

io_header    header1;

// arrays for storing data
float ***ro; // for density/potential
fftwf_plan p_ro; // for FFT
fftwf_plan q_ro; // for FFT


//end of declaration of global variables 

void main()
{
  FILE  *inp, *outpp;
  int ii,jj;
  int Nflag, Noutput;
  int oflag,  // desired output format
    pk_flag; //for power spectrum calculation 
  long seed,*no; // seed for random phases
  float **rra, **vva; // particle position and velocity
  float  aa,delta_aa,delta_aap,afin,aa_i,vpk;  
  float DD,vaa,*nz; 
  char fname[20];
  
  char file[100],num[8];
  
  double t,T=omp_get_wtime(),  // for timing 
    *power, *powerk, *kmode;
  
  float rho_c=2.7755*1.e11; //rho_c in units of h^2 M_sun/Mpc^3
  
  /*---------------------------------------------------------------------------*/
  /* Read input parameters for the simulation from the file "input.nbody_comp" */
  /*---------------------------------------------------------------------------*/
  inp=fopen("input.nbody_comp","r");
  fscanf(inp,"%ld%d",&seed,&Nbin);
  fscanf(inp,"%f%f%f%f",&vhh,&vomegam,&vomegalam,&vnn);
  fscanf(inp,"%f%f",&vomegab,&sigma_8_present);
  fscanf(inp,"%ld%ld%ld%d%f",&N1,&N2,&N3,&NF,&LL);
  fscanf(inp,"%d%d",&oflag,&pk_flag);
  fscanf(inp,"%f%f",&vaa,&delta_aa);  /* time step, final scale factor*/
  fscanf(inp,"%d",&Noutput);
  
  nz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput 
  
  for(ii=0;ii<Noutput;ii++)
    fscanf(inp,"%f",&nz[ii]);
  
  fclose(inp);
  MM=(N1*N2*N3)/pow(NF,3);
  DD=Df(vaa); // growing mode 
  
  DM_m=vomegam*rho_c*pow(vhh,3.)*(1.0*N1)*(1.0*N2)*(1.0*N3)*powf(LL,3.)/(MM*1.e10); //mass per particle in 10^10 M_sun h^-1 unit
  printf("DM_m= %e 10^10 M_sun \t L_box=%e Mpc\n", DM_m/vhh, N1*LL);
  
  /*---------------------------done inputting parameters-----------------*/
  
  /*---------------------------------------------------------------------*/
  /*                           initialize power spectrum                 */
  /*---------------------------------------------------------------------*/
  
  float tcmb=2.728; // CMBR temperature 
  TFset_parameters(vomegam*vhh*vhh,vomegab/vomegam,tcmb);
  
  /*---------------------------------------------------------------------*/
  /*                           done intitializing power spectrum         */
  /*---------------------------------------------------------------------*/
  /*              normalizing  the power spectrum  using sigma_8         */
  /*---------------------------------------------------------------------*/
  norm=1.;
  norm=simp(sigma_func,0.00001,3.5,100000);
  norm=pow(sigma_8_present,2.)/norm;       //normalization factor for Pk(k)
  /*---------------------------------------------------------------------*/
  /*              Normalization of powerspectrum done                    */
  /*---------------------------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  /*                   allocate memory  for particle positions           */
  /*---------------------------------------------------------------------*/
  
  rra= allocate_float_2d(MM,3);
  vva= allocate_float_2d(MM,3);
  
  /*---------------------------------------------------------------------*/
  /*----------allocate memory for power spectrum and k modes-------------*/
  power = calloc((size_t)Nbin,sizeof(double));
  powerk = calloc((size_t)Nbin,sizeof(double));
  kmode = calloc((size_t)Nbin,sizeof(double));
  no = calloc((size_t)Nbin,sizeof(long));  
  /*----------------------------------------------------------------*/
  
  Setting_Up_Memory_For_Ro(vaa); 
  
  /*----------------------------------------------------------------*/
  
  /* get values of delta in fourier Space */
  
  t=omp_get_wtime();
  
  delta_fill(&seed);  //random phases for Fourier modes
  
  printf("ok delta_fill time = %e\n",omp_get_wtime()-t);  
  
  //*------------------------------------------------------------------------*//
  
  if(pk_flag==1)
    {
      t=omp_get_wtime();
      
      outpp=fopen("pk.inp","w");  
      
      calpow(0, Nbin, power, powerk, kmode, no); // calculates power spectrum of (Delta(k))
      
      for(ii=0;ii<Nbin;++ii)
	fprintf(outpp,"%e %e %e %ld\n",kmode[ii],power[ii],DD*DD*2.*pi*pi*powerk[ii],no[ii]);
      
      fclose(outpp);
      
      printf("ok cal_pow time = %e\n",omp_get_wtime()-t);  
    }
  
  //---------------------------------------------------------------------
  
  t=omp_get_wtime();
  
  Zel_move_gradphi(vaa,rra,vva);                // move particles using ZA
  
  printf("ok Zel_move time = %e\n",omp_get_wtime()-t);
  
  /*--------------------------ZA complete--------------------------------*/
  
  if(pk_flag==1)
    {
      t=omp_get_wtime();
      cic(rra);
      printf("ok cic = %e\n",omp_get_wtime()-t);
      
      outpp=fopen("pk.zel","w");
      calpow(1, Nbin, power, powerk, kmode, no); // calculates power spectrum of (Delta rho)
      
      for(ii=0;ii<Nbin;++ii)
	fprintf(outpp,"%e %e %e %ld\n", kmode[ii], power[ii], DD*DD*2.*pi*pi*powerk[ii], no[ii]);
      
      fclose(outpp);
    }
  /*----------------------------------------------------------------------*/
  
  aa=vaa;  //initial scale factor
  
  Update_v(aa, -1.*delta_aa, rra, vva); // Initialize v after ZA
  
  Nflag=0;
  
  //-------------------------Starting Nbody Block---------------------------//
  
  for(jj=0;jj<Noutput;jj++)
    {
      afin=1.0/(nz[jj]+1.0);
      
      t=omp_get_wtime();
      
      /*---------------------------first step ---------------------------------*/
      
      while((afin-aa)<delta_aa)
	{
	  delta_aa=delta_aa*0.5;
	}
      
      Update_v(aa, 0.5*delta_aa, rra, vva);  // half step for v
      Update_x(aa+0.5*delta_aa, delta_aa, rra, vva); // full step for x
      
      /*----------------------first step done---------------------------------*/
      
      aa=aa+delta_aa;
      
      while (aa + delta_aa <=afin)
	{
	  Update_v( aa, delta_aa, rra, vva);
	  Update_x( aa+0.5*delta_aa, delta_aa, rra, vva);
	  aa=aa+delta_aa;
	  Nflag++;
	} /* NBody loop end */
      
      /*-----------------------last step--------------------------------------*/
      
      delta_aap=afin-aa;
      printf("delta_aap=%e\n",delta_aap);
      
      if(delta_aap>0.)
	{
	  Update_v(aa, 0.5*delta_aa, rra, vva);
	  Update_x(aa, delta_aap, rra, vva);
	  
	  aa=aa+delta_aap;
	  Update_v(aa, delta_aap, rra, vva);
	}
      
      printf("(z=%3.1f) time for  Nbody loop = %e\t",nz[jj], omp_get_wtime()-t);
      
      /*------------------- done last step--------------------------------------*/
      
      /*----------------------print final power spectrum -----------------------------------*/
      
      if(pk_flag==1)
	{
	  strcpy(file,"pk.nbody");
	  sprintf(num,"%3.1f",nz[jj]);
	  strcat(file,num);
	  outpp=fopen(file,"w");
	  
	  cic(rra);
	  calpow(1, Nbin, power, powerk, kmode, no);
	  DD=Df(aa); // growing mode
	  
	  for(ii=0;ii<Nbin;++ii)
	    fprintf(outpp,"%e %e %e %ld\n", kmode[ii], power[ii], DD*DD*2.*pi*pi*powerk[ii], no[ii]);
	  
	  fclose(outpp);
	}
      /*----------------------print nbody output -----------------------------------*/
      
      t=omp_get_wtime();
      
      strcpy(file,"output.nbody_");
      sprintf(num,"%3.3f",nz[jj]);
      strcat(file,num);
      
      write_multiout(file,seed,oflag,rra,vva,aa);
      printf("time for write output = %e\n",omp_get_wtime()-t);
 
      aa=afin;
    }
  
  printf("number of step = %d\n",Nflag);
  
  free(rra);
  free(vva);
  free(ro);

  printf("done. Total time taken = %dhr %dmin %dsec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);  
}
