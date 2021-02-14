#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>
#include"nbody.h"
#include<omp.h>
/*---------------------GLOBAL VARIABLES declaration---------------------*/
// Cosmological parameters read from input file  "input.zel"

extern float  vhh, // Hubble parameter 
  vomegam, // Omega_matter; total matter density (baryons+CDM) parameter
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon
  tcmb, //  CMB temperature
  sigma_8_present ,//  Last updated value of sigma_8 (Presently WMAP)
  vnn; // Spectral index of primordial Power spectrum

extern long N1,N2,N3;// box dimension (grid) 
extern int NF, // Fill every NF grid point 
  Nbin; // Number of bins to calculate final P(k) (output)

extern long    MM; // Number of particles 
extern float   LL; // grid spacing in Mpc

// global variables  (calculated )

extern int zel_flag, // memory allocation for zel is 3 times that for nbody
  fourier_flag; //for fourier transfrom

extern float   DM_m, // Darm matter mass of simulation particle (calculated) 
  norm, // normalize Pk
  pi;


extern   io_header  header1; // structure for header 
extern float ***ro;  // for density or potential
extern fftwf_plan p_ro; // for FFT
extern fftwf_plan q_ro; // for FFT


// arrays for storing data
static float ***va; // for grad phi
static float  norml, CC; // normalization constants
static float rho_b_inv, // (Mean number of particles/grid cell)^{-1}
  Cx,Cy,Cz, // conversion from i,j,k to kx,ky,kz
  Lcube,vol,  // L^3, Volume=L^3*N1*N2*N3
  tpibyL;// 2 pi /LL
 
static fftwf_plan q_va; // for FFT
static omp_lock_t *ro_lck;

/*********************************************************************************************************/
void Setting_Up_Memory_For_Ro(float av)
{  
  // for multiple  threads
  fftwf_init_threads();
  
  fftwf_plan_with_nthreads(omp_get_max_threads());
  // fftwf_plan_with_nthreads(1);
  
  omp_set_num_threads(omp_get_max_threads());
  //  omp_set_num_threads(1);
  
  printf("No of threads = %d\n",omp_get_max_threads());
  // done multi thread 
  
  ro = allocate_fftwf_3d(N1,N2,N3+2);
  va = allocate_fftwf_3d(N1,N2,N3+2); /* The last dimension gets padded because of using REAL FFT */
  
  ro_lck=(omp_lock_t*)&(va[0][0][0]); 
  /****************************************************/	
  /* Creating the plans for forward and reverse FFT's */
  
  p_ro = fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro[0][0][0]), (fftwf_complex*)&(ro[0][0][0]), FFTW_ESTIMATE);  
  q_ro = fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro[0][0][0]), &(ro[0][0][0]), FFTW_ESTIMATE);
  q_va = fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(va[0][0][0]), &(va[0][0][0]), FFTW_ESTIMATE);
  
  /****************************************************/
  /* Fwd FFT followed immediately by Reverse FFT gives results reduced by a product N1*N2*N3 */
  /* The FFT is unnormalized.. So we introduce a normalization factor */
  
  norml=1./((1.0*N1)*(1.0*N2)*(1.0*N3));
  Cx=2.*pi/(N1*LL);  Cy=2.*pi/(N2*LL);   Cz=2.*pi/(N3*LL);
  Lcube=powf(LL,3.);
  vol=Lcube/norml;
  CC=pi*pi*vol*powf(Df(av),2.); //  (2*pi^2 Volume DD^2)/2  DD-> initial growing mode 
  rho_b_inv=1.0/(norml*MM); // inverse of density in units of particles per unit cell
  tpibyL=2.0*pi/LL;// 2 pi /LL 
  
  printf("norml=%e\tCC=%e\trho_b_inv=%e\n",norml,CC,rho_b_inv);
}

/*******************************************************************************************/

void delta_fill(long* seed)
{
  long i,j,k,m,jj1,jj2,kk;
  long index1,index2;
  float amp,phasev,DD;
  fftwf_complex *A;
  
  float ran1(long*);
  float gasdev(long*);
  
  float p(long,long,long); // returns amplitude of delta(k)
  
  
  A=(fftwf_complex*)&(ro[0][0][0]);
  
  /********************FILLING POINTS***********************/
  
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      for(k=0;k<2;++k)
	{
	  if (i==0 && j==0 && k==0)
	    {
	      index1=0;
	      A[index1][0]=0.0;
	      A[index1][1]=0.0;
	    }
	  else
	    {
	      index1 = i*(N1/2)*N2*(N3/2+1) + j*(N2/2)*(N3/2+1) + k*(N3/2);
	      
	      amp=p(i*N1/2,j*N2/2,k*N3/2);
	      A[index1][0]=amp*gasdev(seed);
	      A[index1][1]=0.0;
	    }
	}
  
  /********************FILLING LINES***********************/
  
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	for(k=1;k<N1/2;k++)
	  {
	    amp=p(k,i*N2/2,j*N3/2);
	    
	    index1 =k*N2*(N3/2+1) +i*(N2/2)*(N3/2+1) + j*(N3/2);
	    
	    A[index1][0]=amp*gasdev(seed);
	    A[index1][1]=amp*gasdev(seed);
	    
	    index2 =(N1-k)*N2*(N3/2+1) +i*(N2/2)*(N3/2+1) + j*(N3/2);
	    
	    A[index2][1]=A[index1][0];
	    A[index2][1]=-A[index1][1];
	  }
	
	for(k=1;k<N2/2;k++)
	  {
	    amp=p(i*N1/2,k,j*N3/2);
	    
	    index1 =i*(N1/2)*N2*(N3/2+1)+ k*(N3/2+1) + j*(N3/2);
	    index2 =i*(N1/2)*N2*(N3/2+1)+(N2-k)*(N3/2+1) + j*(N3/2);
	    
	    A[index1][0]=amp*gasdev(seed);
	    A[index1][1]=amp*gasdev(seed);
	    
	    A[index2][0]=A[index1][0];
	    A[index2][1]=-A[index1][1];
	  }
	
	for(k=1;k<N3/2;k++)
	  {
	    amp=p(i*N1/2,j*N2/2,k);
	    
	    index1 =i*(N1/2)*N2*(N3/2+1)+ j*(N2/2)*(N3/2+1) + k;
	    
	    A[index1][0]=amp*gasdev(seed);
	    A[index1][1]=amp*gasdev(seed);
	  }
      }
  
  
  /**********************FILLING PLANES***********************/
  for(i=0;i<2;i++)
    for(m=0;m<2;m++)
      {
	for(j=1;j<N2/2;j++)
	  for(k=1;k<N3/2;k++)
	    {
	      amp=p(i*N1/2,j,k);
	      
	      jj1=(1-m)*j + m * (N2-j);	    
	      index1 =i*(N1/2)*N2*(N3/2+1)+ jj1*(N3/2+1) + k;
	      
	      A[index1][0]=amp*gasdev(seed);
	      A[index1][1]=amp*gasdev(seed);	      
	    }
	
	for(j=1;j<N1/2;j++)
	  for(k=1;k<N3/2;k++)
	    {
	      amp=p(j,i*N2/2,k);
	      
	      jj1=(1-m)*j + m * (N1-j);
	      index1 =jj1*N2*(N3/2+1)+ i*(N2/2)*(N3/2+1) + k;
	      
	      A[index1][0]=amp*gasdev(seed);
	      A[index1][1]=amp*gasdev(seed);
	    }	      
	
	for(j=1;j<N1/2;j++)
	  for(k=1;k<N2/2;k++)
	    {
	      amp=p(j,k,i*N3/2);
	      
	      jj1=(1-m)*j + m * (N1-j);
	      kk=k;
	      index1 =jj1*N2*(N3/2+1)+ kk*(N3/2+1)+i*N3/2;
	      
	      jj2=m*j + (1-m) * (N1-j);
	      kk=N2-k;
	      index2 =jj2*N2*(N3/2+1)+ kk*(N3/2+1)+i*N3/2;
	      
	      A[index1][0]=amp*gasdev(seed);
	      A[index1][1]=amp*gasdev(seed);
	      
	      A[index2][0]=A[index1][0];
	      A[index2][1]=-A[index1][1];
	    }
      }
  /***************************FILLING CUBES***********************/ 
  
  for(i=1;i<N1/2;i++)
    for(j=1;j<N2/2;j++)
      for(k=1;k<N3/2;k++)
	{
	  amp=p(i,j,k);
	  
	  index1 =i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  A[index1][0]=amp*gasdev(seed);
	  A[index1][1]=amp*gasdev(seed);
	  
	  index1 =(N1-i)*N2*(N3/2+1) + j*(N3/2+1) + k;
	  A[index1][0]=amp*gasdev(seed);
	  A[index1][1]=amp*gasdev(seed);
	  
	  index1 =i*N2*(N3/2+1) + (N2-j)*(N3/2+1) + k;
	  A[index1][0]=amp*gasdev(seed);
	  A[index1][1]=amp*gasdev(seed);
	  
	  index1 = (N1-i)*N2*(N3/2+1) + (N2-j)*(N3/2+1) + k;
	  A[index1][0]=amp*gasdev(seed);
	  A[index1][1]=amp*gasdev(seed);
	}
} 
/**************************************************************************************************/

float p(long i,long j,long k)
// returns amplitude of delta(k)
{
  float kk,val,Pk(float),kx,ky,kz;
  
  kx=Cx*i;
  ky=Cy*j;
  kz=Cz*k;
  
  kk = (float)sqrt((double)(kx*kx+ ky*ky + kz*kz)); 
  
  val = (float)sqrt(CC*Pk(kk));
  return(val);
}


/*******************************************************************************************************/
/*******************************************************************************************************/
/******************** TO FIND POWER SPECTRUM **************************/

void calpow(int f_flag, int Nbin, double* power, double* powerk, double* kmode, long *no)
{
  long i, j, k, a, b, c;
  int d;
  fftwf_complex *comp_ro;
  float fac1,fac2,fac3, fac, m, scale;
  long index,index1,index2,c1=0;
  
  // f-flag==0 -> ro contains delta(k)
  // f-flag==1 -> ro contains delta(x)
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  // do forward  transform to get delta(k)
  
  if(f_flag==1)
    {
      for(i=0;i<N1;i++)
	for(j=0;j<N2;j++)
	  for(k=0;k<N3;k++)
	    ro[i][j][k]=ro[i][j][k]*Lcube;
      
      fftwf_execute(p_ro);  // ro mow contains delta(k)
    }
  
  comp_ro = (fftwf_complex *)&(ro[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  
  for(i=0;i<Nbin;++i)
    {
      power[i]=0.0;
      powerk[i]=0.0;
      kmode[i]=0.0;
      no[i]=0;
    }
  
  // dor logarithmic bins
  
  scale=log10(0.5*N1)/Nbin;
  
  // ----------------------
  
  
  /**************** BINNING POWER SPECTRA **********************/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  a=i;
	  b=j;
	  c=k;
	  
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	  
	  d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
	      powerk[d] += Pk((float)(tpibyL*m))*1.0;
	      kmode[d] += m*1.;
	      no[d]=no[d]+1;
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  powerk[d] += Pk((float)(tpibyL*m))*1.0; 
		  kmode[d] += m*1.;
		  no[d]=no[d]+1;
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      /* m*(2 * pi/LL) is |k| */
	      /* m=1/2 corresponds to kmode[Nbin-1] i.e. Nyquits */
	      
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if (d>=0 && d<Nbin)
		{
		  power[d] += ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  powerk[d] += Pk((float)(tpibyL*m))*1.0;
		  kmode[d] += m*1.;
                  no[d]=no[d]+1;
		}
	    } /* end k for loop */
	  
	}
    }
  
  
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0)
	{
	  power[i] =  power[i]/(1.*no[i]*vol);
	  powerk[i] = powerk[i]/(1.*no[i]);
	  kmode[i]=tpibyL*kmode[i]/(1.*no[i]);
	}
    }
  
  // do back transform to get delta(x)
  if(f_flag==1)
    {
      for(i=0;i<N1;i++)
	{
	  index1 = i*N2*(N3/2+1) ;	  
	  for(j=0;j<N2;j++)
	    {
	      index2=index1 + j*(N3/2+1) ;
	      for(k=0;k<N3/2;k++)
		{
		  index=index2 + k;
		  
		  comp_ro[index][0]=comp_ro[index][0]/vol;
		  comp_ro[index][1]=comp_ro[index][1]/vol;
		}
	    }
	}
      /*  now convert the array back to real space */
      fftwf_execute(q_ro);
    }
}

/**********************************************************************/
/******************** TO FIND POWER SPECTRUM **************************/

void calpow_k(int f_flag, float kmin, float kmax, int Nbin,double* power, double* powerk, double* kmode, long *no)
{
  long i, j, k, a, b, c;
  int d;
  fftwf_complex *comp_ro;
  float fac1,fac2,fac3, fac, m, scale;
  long index,index1,index2;
  
  // f-flag==0 -> ro contains delta(k)
  
  // f-flag==1 -> ro contains delta(x)
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  
  // do forward  transform to get delta(k)
  
  if(f_flag==1)
    {
      for(i=0;i<N1;i++)
	for(j=0;j<N2;j++)
	  for(k=0;k<N3;k++)
	    ro[i][j][k]=ro[i][j][k]*Lcube;
      
      fftwf_execute(p_ro);  // ro mow contains delta(k)
    }
  
  comp_ro = (fftwf_complex *)&(ro[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  
  for(i=0;i<Nbin;++i)
    {
      power[i]=0.0;
      powerk[i]=0.0;
      kmode[i]=0.0;
      no[i]=0;
    }
  
  // dor logarithmic bins
  
  scale=(log10(kmax)-log10(kmin))/Nbin;
  
  // ----------------------
  
  /**************** BINNING POWER SPECTRA **********************/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  a=i;
	  b=j;
	  c=k;
	  
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	  
	  d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
	      powerk[d] += Pk((float)(m))*1.0;
	      kmode[d] += m*1.;
	      no[d]=no[d]+1;
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      
	      index = index2 + k;
	      
	      m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	      
	      d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  powerk[d] += Pk((float)(m))*1.0;
		  kmode[d] += m*1.;
		  no[d]=no[d]+1;
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      index = index2 + k;
	      
	      m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	      
	      d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	      
	      if (d>=0 && d<Nbin)
		{
		  power[d] += ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  powerk[d] += Pk((float)(m))*1.0;
		  kmode[d] += m*1.;
                  no[d]=no[d]+1;
		}
	    } /* end k for loop */
	  
	}
    }
  
  
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0)
	{
	  power[i] =  power[i]/(1.*no[i]*vol);
	  powerk[i] = powerk[i]/(1.*no[i]);
	  kmode[i] = kmode[i]/(1.*no[i]);
	}
    }
  
  // do back transform to get delta(x)
  if(f_flag==1)
    {
      for(i=0;i<N1;i++)
	{
	  index1 = i*N2*(N3/2+1) ;	  
	  for(j=0;j<N2;j++)
	    {
	      index2=index1 + j*(N3/2+1) ;
	      for(k=0;k<N3/2;k++)
		{
		  index=index2 + k;
		  
		  comp_ro[index][0]=comp_ro[index][0]/vol;
		  comp_ro[index][1]=comp_ro[index][1]/vol;
		}
	    }
	}
      /*  now convert the array back to real space */
      fftwf_execute(q_ro);
    }
}


/*******************************************************************************************************/

/*******************************************************************************************************/

void grad_phi(int ix)
{
  //  at start of this function: ro contains delta(k)
  //  at end of function va  contains  the ix component of Grad[Laplacian^-1[delta(x)]]
  // this is used in ZA
  
  long ii,jj,kk;
  long index;
  float AA,a0,a1,a2;
  fftwf_complex *vc,*roc;
  
  vc=(fftwf_complex*)&(va[0][0][0]);
  roc=(fftwf_complex*)&(ro[0][0][0]);
  
#pragma omp parallel for private(jj,kk,a0,a1,a2,AA,index)
  for(ii=0;ii<N1;ii++)
    {   
      a0 = (ii>N1/2)? Cx*(ii-N1) : Cx*ii; // kx
      
      for(jj=0;jj<N2;jj++)
	{
	  a1 = (jj>N2/2)? Cy*(jj-N2) : Cy*jj; // ky
	  
	  kk=0;
	  a2=0.;
	  index = (ii*N2+jj)*(N3/2+1) + kk;	    	    
	  
	  
	  AA = pow(a0,2.)+pow(a1,2.)+pow(a2,2.);  // AA=K^2
	  
	  AA = (fabs(AA)>1.e-10)? (((ix-1)*(ix-2)/2)*a0 - ix*(ix-2)*a1 + (ix*(ix-1)/2)*a2)/AA : 0.0;  // AA=k_ix/K^2 or 0 
	  
	  vc[index][0]=-1.0*AA*roc[index][1]/vol;
	  vc[index][1]=AA*roc[index][0]/vol;	    	
	  
	  for(kk=1;kk<N3/2+1;kk++)
	    {	    
	      a2=kk*Cz; // kz
	      
	      index = (ii*N2+jj)*(N3/2+1) + kk;
	      
	      AA= pow(a0,2.)+pow(a1,2.)+pow(a2,2.);	// AA=K^2
	      AA=(((ix-1)*(ix-2)/2)*a0 - ix*(ix-2)*a1 + (ix*(ix-1)/2)*a2)/AA;  // AA=k_ix/K^2
	      
	      vc[index][0] = -1.0*AA*roc[index][1]/vol;
	      vc[index][1] = AA*roc[index][0]/vol;
	      
	    }
	}
    }
  
  fftwf_execute(q_va); /* Take fourier transform */
}

//############################################################################################
//############################################################################################

void Zel_move_gradphi(float av,float **rra,float **vva)
{
  float  DD, // growing mode
    vfac; // converts to peculiar velocity
  
  int ii;
  long jj, kk, ll, N, pos;
  long pin;
  
  
  vfac = av*av*Hf(av)*ff(av)/LL;
  
  for(ii=0;ii<3;++ii)
    {
      grad_phi(ii); // Calculate ix component of -Grad[Laplacian^{_1}[delta]]  in real space
      
#pragma omp parallel for private(jj,kk,ll,pos,pin,N)
      for(jj=0;jj<N1/NF;jj++)
	for(kk=0;kk<N2/NF;kk++)
	  for(ll=0;ll<N3/NF;ll++)
	    {
	      pin = jj*(N2/NF)*(N3/NF) + kk*(N3/NF) + ll;
	      
	      pos= (ii-1)*(ii-2)*jj/2 - ii*(ii-2)*kk + ii*(ii-1)*ll/2;

	      rra[pin][ii]=(float)(NF*pos)+ va[NF*jj][NF*kk][NF*ll]/LL;
	      vva[pin][ii]=vfac*va[NF*jj][NF*kk][NF*ll];
	      
	      // periodic boundary condition 
	      
	      N = (ii-1)*(ii-2)*N1/2 - ii*(ii-2)*N2 + ii*(ii-1)*N3/2;
	      
	      /* to ensure the values are not negative */
	      rra[pin][ii]=rra[pin][ii]+N;
	      rra[pin][ii]= rra[pin][ii]-1.0*N*(int)(floor(rra[pin][ii])/(1.*N));
	    }
      /*------------------- done ---------------------*/
    }
  
} /* end of  Zel_move_gradphi */

//######################################################################################

void cic(float **rra)
// //For a  given particle distribution 
// This uses Cloud in Cell(CIC)  to calculate (1 + delta) on the  Grid 
{  
  long i, j, k, ix, jy, kz;
  int ii, jj, kk; 
  long  pin,index;
  float wx,wy,wz;

  /* Clear out the array ro. ******/
#pragma omp parallel for private(i,j,k,index)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
	{
	  index = (i*N2+j)*N3 + k;
	  ro[i][j][k] = 0.0;
	  omp_init_lock(&ro_lck[index]);
	}

   
  /********************************/
#pragma omp parallel for private(i, j, k, ii, jj, kk, wx, ix, wy, jy, wz, kz, index)
  for(pin=0;pin<MM;pin++)
    { /* begin particle index loop */
      /* (a/b/c)[0] or (a/b/c)[1] can never be greater than (N1/N2/N3) */
      // left bottom corner of cell containing the particle
      
      i = (long)floor(rra[pin][0]);
      j = (long)floor(rra[pin][1]);
      k = (long)floor(rra[pin][2]);
      
      /* for each of the 8 corner points */
      for(ii=0;ii<=1;ii++)
	{
	  wx=fabs(1.-rra[pin][0]+i-ii)*rho_b_inv; // divide by mean density    
	  ix=(i+ii)%N1;
	  
	  for(jj=0;jj<=1;jj++)
	    {
	      wy=fabs(1.-rra[pin][1]+j-jj);
	      jy=(j+jj)%N2;
	      
	      for(kk=0;kk<=1;kk++)
		{ 
		  wz=fabs(1.-rra[pin][2]+k-kk);       
		  kz=(k+kk)%N3;

		  index=(ix*N2+jy)*N3 + kz;

		  omp_set_lock(&ro_lck[index]);
		  ro[ix][jy][kz]+=wx*wy*wz;
		  omp_unset_lock(&ro_lck[index]);
		}
	    }
	} /* end of 8 grid corners loop>	*/
    } /* end of each particle loop */

} /* end of function cic */

//############################################################################################
//############################################################################################

void cic_sampled(float **rra, int *s_indx)
// //For a  given particle distribution 
// This uses Cloud in Cell(CIC)  to calculate (1 + delta) on the  Grid 
{  
  long i, j, k, ix, jy, kz;
  int ii, jj, kk; 
  long  pin,index;
  float wx,wy,wz;

  /* Clear out the array ro. ******/
#pragma omp parallel for private(i,j,k,index)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
	{
	  index = (i*N2+j)*N3 + k;
	  ro[i][j][k] = 0.0;
	  omp_init_lock(&ro_lck[index]);
	}

  /********************************/
#pragma omp parallel for private(i, j, k, ii, jj, kk, wx, ix, wy, jy, wz, kz, index)
  for(pin=0;pin<MM;pin++)
   if(s_indx[pin]==-1)
    { /* begin particle index loop */
      /* (a/b/c)[0] or (a/b/c)[1] can never be greater than (N1/N2/N3) */
      // left bottom corner of cell containing the particle
      
      i = (long)floor(rra[pin][0]);
      j = (long)floor(rra[pin][1]);
      k = (long)floor(rra[pin][2]);
      
      /* for each of the 8 corner points */
      for(ii=0;ii<=1;ii++)
	{
	  wx=fabs(1.-rra[pin][0]+i-ii)*rho_b_inv; // divide by mean density    
	  ix=(i+ii)%N1;
	  
	  for(jj=0;jj<=1;jj++)
	    {
	      wy=fabs(1.-rra[pin][1]+j-jj);
	      jy=(j+jj)%N2;
	      
	      for(kk=0;kk<=1;kk++)
		{ 
		  wz=fabs(1.-rra[pin][2]+k-kk);       
		  kz=(k+kk)%N3;

		  index=(ix*N2+jy)*N3 + kz;

		  omp_set_lock(&ro_lck[index]);
		  ro[ix][jy][kz]+=wx*wy*wz;
		  omp_unset_lock(&ro_lck[index]);
		}
	    }
	} /* end of 8 grid corners loop>	*/
    } /* end of each particle loop */

} /* end of function cic_sampled */

//##########################################################################################

void Get_phi(int f_flag)  // calculates Laplacian^{-1}[ ro]]
{
  long ii,jj,kk;
  long index;
  float AA,a0,a1,a2;
  fftwf_complex *roc;
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  // do forward  transform to get ro(k)
  if(f_flag==1)
    {
      for(ii=0;ii<N1;ii++)
	for(jj=0;jj<N2;jj++)
	  for(kk=0;kk<N3;kk++)
	    ro[ii][jj][kk]*=Lcube;  
      
      fftwf_execute(p_ro);  // ro mow contains delta(k)
    }
  
  //------------------------------------------
  
  roc=(fftwf_complex*)&(ro[0][0][0]);
  
#pragma omp parallel for private(jj,kk,a0,a1,a2,AA,index)
  for(ii=0;ii<N1;ii++)
    {   
      a0=ii*pi/N1; // kx *LL/2
      a0=pow(2.*sin(a0)/LL,2.);
      
      for(jj=0;jj<N2;jj++)
	{
	  a1=jj*pi/N2; // ky *LL/2
	  a1=pow(2.*sin(a1)/LL,2.);
	  
	  kk=0;	  
	  index = (ii*N2+jj)*(N3/2+1) + kk;
	    	    
	  AA=a0+a1;
	  AA=(fabs(AA)>1.e-10)? 1/AA:0.;  //  AA=k_ix/K^2 or 0 
	  
	  roc[index][0] = -1.*AA*roc[index][0]/vol;
	  roc[index][1] = -1.*AA*roc[index][1]/vol;	    	
	  
	  for(kk=1;kk<N3/2+1;kk++)
	    {	    
	      a2=kk*pi/N3; // kz *LL/2
	      a2=pow(2.*sin(a2)/LL,2.);
	      
	      index = (ii*N2+jj)*(N3/2+1) + kk;	    	    
	      AA=1./(a0+a1+a2);	       
	      
	      roc[index][0] = -1.*AA*roc[index][0]/vol;
	      roc[index][1] = -1.*AA*roc[index][1]/vol;	    		   
	      
	    }
	}
    }
  
  fftwf_execute(q_ro); /* Take fourier transform */
}

//#############################################################################################

void Update_v(float av,float delta_aa,float **rra,float **vva)  // udates v if (delta a >0) else intializes v
{
  long pin;
  int ii,jj,kk;
  long xp,yp,zp,xn,yn,zn; /* (x/y/z)(p/n) = (x/y/z)(previous/next) */  
  long a,b,c,ix,jy,kz;

  float wx,wy,wz,g0,g1,g2;
  float xx,yy,zz;
  float coeff,Hf(float aa),vflag;

  cic(rra); // calculate delta(x) from particle distribution
  Get_phi(1);// Calculate  [Laplacian^{_1}[delta]]  in real space

 
  if(delta_aa>0.)
    {
      coeff=delta_aa*1.5* vomegam/(LL*av*av*Hf(av)); // delta_aa scaled 
      vflag=1.0;
    }
  else
    {
      coeff = av*av*Hf(av)*ff(av)/LL;
      vflag=0.0;
    }
  
  coeff/=(2.*LL); // divide by 2 LL for Grad on grid 
  
#pragma omp parallel for private(a,b,c,xx,yy,zz,g0,g1,g2,ii,jj,kk,ix,wx,xp,xn,jy,wy,yn,yp,kz,wz,zp,zn)	   
  for(pin=0;pin<MM;pin++)
    { 
      /* left most corner of the cube enclosing the particle */
      
      a = (long)floor(rra[pin][0]);
      b = (long)floor(rra[pin][1]);
      c = (long)floor(rra[pin][2]);
      
      /* particle co-ordinates itself */
      xx = rra[pin][0];
      yy = rra[pin][1];
      zz = rra[pin][2];
      
      g0=0;
      g1=0;
      g2=0;
      
      /* begin 8 corners loop */
      for(ii=0;ii<=1;ii++)
	{
	  ix = (a+ii)%N1;
	  wx=fabs(1.- xx +a-ii);    
	  xp = ((ix-1)==-1) ? N1-1: ix-1;
	  xn = ((ix+1)==N1) ? 0   : ix+1;
	  
	  for(jj=0;jj<=1;jj++)
	    {
              jy = (b+jj)%N2;
	      wy=fabs(1.- yy +b-jj);       
              yp = ((jy-1)==-1) ? N2-1: jy-1;
              yn = ((jy+1)==N2) ? 0   : jy+1;
	      
	      for(kk=0;kk<=1;kk++)
		{ 		  
		  kz = (c+kk)%N3;
		  wz=fabs(1.- zz +c-kk);      
		  
		  zp = ((kz-1)==-1) ? N3-1: kz-1;
		  zn = ((kz+1)==N3) ? 0   : kz+1;
		  
		  /* ix,jy,kz are the current co-ordinates of the cube vertex point */
		  /* calculating the difference from the respective corner */
		  
		  g0 += wx*wy*wz*(ro[xp][jy][kz] - ro[xn][jy][kz]);
		  g1 += wx*wy*wz*(ro[ix][yp][kz] - ro[ix][yn][kz]);
		  g2 += wx*wy*wz*(ro[ix][jy][zp] - ro[ix][jy][zn]);
		} 
	    }
	} /* End of 8 corners loop */
      //for (ii=0;ii<3;++ii)
      vva[pin][0] = vflag*vva[pin][0]+ g0*coeff;
      vva[pin][1] = vflag*vva[pin][1]+ g1*coeff;
      vva[pin][2] = vflag*vva[pin][2]+ g2*coeff;
    } /* end particle index  */        
  
} /* end update_v_FDA */

//####################################################################################

void Update_x(float aa,float delta_aa,float **rra,float **vva)
// updates x
{
  float coeff;
  int  ii;
  long N;
  long pin;
  float Hf(float aa);
  
  coeff=delta_aa/(aa*aa*aa*Hf(aa)); // delta_aa scaled 
  
#pragma omp parallel for private(ii,N)
  for(pin=0;pin<MM;pin++)  /* begin particle index loop */
    {
      for (ii=0;ii<3;++ii)   /* begin co-ordinate x/y/z loop */
	{
	  
	  rra[pin][ii]=rra[pin][ii]+vva[pin][ii]*coeff;	    	    
	  
	  /* imposing periodic boundary condition */
	  N = (ii-1)*(ii-2)*N1/2 - ii*(ii-2)*N2 + ii*(ii-1)*N3/2;
	  
	  /* to ensure the values are not negative */
	  rra[pin][ii] = rra[pin][ii]+N;
	  rra[pin][ii] = rra[pin][ii]-1.0*N*(int)(floor(rra[pin][ii])/(1.*N));

	} /* end of x-y-z co-od loop */
    } /* end of particle loop */
} /* end function */

//####################################################################################

/***************************************************************************/
/*                           IO  FUNCTIONS                                 */
/***************************************************************************/

int write_multiout(char *fname,long int seed,int output_flag,float **rra,float **vva,float vaa)
{
  FILE *fp1; 
  long   ii; 
  fp1=fopen(fname,"w");   
  int dummy, kk; 
  // set header structure values 
  
  for(kk=0;kk<6;++kk) 
    { 
      header1.npart[kk]=0; 
      header1.mass[kk]=0.; 
      header1.npartTotal[kk]=0; 
    } 
  header1.npart[1]=(long)MM; //DM particles in this snapshot file 
  header1.mass[1]=DM_m; //DM particle mass in units of 10^10 M_sun/h 
  header1.npartTotal[1]=(long)MM; //Total DM particles in simulation 
  header1.time=vaa;     //scale factor of  nbody output     
  header1.redshift=1./vaa-1.; 
  header1.flag_sfr=0; 
  header1.flag_cooling=0; 
  header1.flag_feedback=0; 
  header1.num_files=1;  //no. of files in each snapshot 
  header1.BoxSize=N1*LL*1000*vhh;//simulation box size in kpc/h 
  header1.Omega0=vomegam;    //Omega_m   
  header1.OmegaLambda=vomegalam;     //OmegaLambda     
  header1.HubbleParam=vhh;     //HubbleParam 
  header1.Omegab=vomegab;    //Omega_b      
  header1.sigma_8_present=sigma_8_present;  
  header1.Nx=N1;
  header1.Ny=N2;
  header1.Nz=N3; 
  header1.LL=LL; 
  header1.output_flag=output_flag; 
  header1.in_flag=zel_flag; 
  header1.seed=seed; 
  // done setting header  
  if(output_flag!=1) 
    { 
      //final paticle positions stored in kpc/h unit and velocity written in km/sec unit  
      for (ii=0;ii<MM;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]*LL*1000.*vhh;//coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]*LL*1000.*vhh; 
 	  rra[ii][2]=rra[ii][2]*LL*1000.*vhh; 
	  
 	  vva[ii][0]=vva[ii][0]*LL*vhh*100./vaa ;//peculiar velocities in km/sec 
 	  vva[ii][1]=vva[ii][1]*LL*vhh*100./vaa; 
 	  vva[ii][2]=vva[ii][2]*LL*vhh*100./vaa; 
 	} 
    }
  
  // writing header
  
  fwrite(&dummy,sizeof(dummy),1,fp1);     
  fwrite(&header1,sizeof(io_header),1,fp1);     
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  // writing data  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++) 
    fwrite(&rra[ii][0],sizeof(float),3,fp1); 
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++) 
    fwrite(&vva[ii][0],sizeof(float),3,fp1); 
  fwrite(&dummy,sizeof(dummy),1,fp1); 

  fclose(fp1);


  if(output_flag!=1) 
    { 
      //paticle positions restored in grid units  
      for (ii=0;ii<MM;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]/(LL*1000.*vhh);  //coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]/(LL*1000.*vhh); 
 	  rra[ii][2]=rra[ii][2]/(LL*1000.*vhh); 
	  
 	  vva[ii][0]=(vva[ii][0]*vaa)/(LL*vhh*100.);  //peculiar velocities in km/sec 
 	  vva[ii][1]=(vva[ii][1]*vaa)/(LL*vhh*100.);
 	  vva[ii][2]=(vva[ii][2]*vaa)/(LL*vhh*100.);
 	} 
    }

} 
/*-------------------------------------------------------------------------------------------------------*/
int write_output(char *fname,long int seed,int output_flag,float **rra,float **vva,float vaa)
{
  FILE *fp1; 
  long   ii; 
  fp1=fopen(fname,"w");   
  int dummy, kk; 
  // set header structure values 
  
  for(kk=0;kk<6;++kk) 
    { 
      header1.npart[kk]=0; 
      header1.mass[kk]=0.; 
      header1.npartTotal[kk]=0; 
    } 
  header1.npart[1]=(long)MM; //DM particles in this snapshot file 
  header1.mass[1]=DM_m; //DM particle mass in units of 10^10 M_sun/h 
  header1.npartTotal[1]=(long)MM; //Total DM particles in simulation 
  header1.time=vaa;     //scale factor of  nbody output     
  header1.redshift=1./vaa-1.; 
  header1.flag_sfr=0; 
  header1.flag_cooling=0; 
  header1.flag_feedback=0; 
  header1.num_files=1;  //no. of files in each snapshot 
  header1.BoxSize=N1*LL*1000*vhh;//simulation box size in kpc/h 
  header1.Omega0=vomegam;    //Omega_m   
  header1.OmegaLambda=vomegalam;     //OmegaLambda     
  header1.HubbleParam=vhh;     //HubbleParam 
  header1.Omegab=vomegab;    //Omega_b      
  header1.sigma_8_present=sigma_8_present;  
  header1.Nx=N1;
  header1.Ny=N2;
  header1.Nz=N3; 
  header1.LL=LL; 
  header1.output_flag=output_flag; 
  header1.in_flag=zel_flag; 
  header1.seed=seed; 
  // done setting header  
  if(output_flag!=1) 
    { 
      //final paticle positions stored in kpc/h unit and velocity written in km/sec unit  
      for (ii=0;ii<MM;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]*LL*1000.*vhh;//coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]*LL*1000.*vhh; 
 	  rra[ii][2]=rra[ii][2]*LL*1000.*vhh; 
	  
 	  vva[ii][0]=vva[ii][0]*LL*vhh*100./vaa ;//peculiar velocities in km/sec 
 	  vva[ii][1]=vva[ii][1]*LL*vhh*100./vaa; 
 	  vva[ii][2]=vva[ii][2]*LL*vhh*100./vaa; 
 	} 
    } 
  fwrite(&dummy,sizeof(dummy),1,fp1);     
  fwrite(&header1,sizeof(io_header),1,fp1);     
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  // header  written 
  // writing data  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++) 
    fwrite(&rra[ii][0],sizeof(float),3,fp1); 
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++) 
    fwrite(&vva[ii][0],sizeof(float),3,fp1); 
  fwrite(&dummy,sizeof(dummy),1,fp1); 

  fclose(fp1);
} 
/*-------------------------------------------------------------------------------------------------------*/

int write_sampled(char *fname,int *ss_indx,long fact,long int seed,int output_flag,float **rra,float **vva,float vaa)
{
  FILE *fp1; 
  long   ii; 
  fp1=fopen(fname,"w");   
  int dummy, kk; 
  // set header structure values 
  
  for(kk=0;kk<6;++kk) 
    { 
      header1.npart[kk]=0; 
      header1.mass[kk]=0.; 
      header1.npartTotal[kk]=0; 
    } 
  header1.npart[1]=(long)MM/fact; //DM particles in this snapshot file 
  header1.mass[1]=DM_m*fact*1.; //DM particle mass in units of 10^10 M_sun/h 
  header1.npartTotal[1]=(long)MM/fact; //Total DM particles in simulation 
  header1.time=vaa;     //scale factor of  nbody output     
  header1.redshift=1./vaa-1.; 
  header1.flag_sfr=0; 
  header1.flag_cooling=0; 
  header1.flag_feedback=0; 
  header1.num_files=1;  //no. of files in each snapshot 
  header1.BoxSize=N1*LL*1000*vhh;//simulation box size in kpc/h 
  header1.Omega0=vomegam;    //Omega_m   
  header1.OmegaLambda=vomegalam;     //OmegaLambda     
  header1.HubbleParam=vhh;     //HubbleParam 
  header1.Omegab=vomegab;    //Omega_b      
  header1.sigma_8_present=sigma_8_present;  
  header1.Nx=N1;
  header1.Ny=N2;
  header1.Nz=N3; 
  header1.LL=LL; 
  header1.output_flag=output_flag; 
  header1.in_flag=zel_flag; 
  header1.seed=seed; 
  // done setting header  
  printf("MM= %ld\n",MM);
   printf("DM = %e \n",DM_m);
  if(output_flag!=1) 
    { 
      //final paticle positions stored in kpc/h unit and velocity written in km/sec unit  
      for (ii=0;ii<MM;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]*LL*1000.*vhh;//coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]*LL*1000.*vhh; 
 	  rra[ii][2]=rra[ii][2]*LL*1000.*vhh; 
	  
 	  vva[ii][0]=vva[ii][0]*LL*vhh*100./vaa ;//peculiar velocities in km/sec 
 	  vva[ii][1]=vva[ii][1]*LL*vhh*100./vaa; 
 	  vva[ii][2]=vva[ii][2]*LL*vhh*100./vaa; 
 	} 
    } 
  fwrite(&dummy,sizeof(dummy),1,fp1);     
  fwrite(&header1,sizeof(io_header),1,fp1);     
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  // header  written 
  // writing data  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++)
  {
   if(ss_indx[ii]!=0) 
       fwrite(&rra[ii][0],sizeof(float),3,fp1); 
  }
   fwrite(&dummy,sizeof(dummy),1,fp1); 
  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++)
  { 
   if(ss_indx[ii]!=0) 
     fwrite(&vva[ii][0],sizeof(float),3,fp1); 
  }
   fwrite(&dummy,sizeof(dummy),1,fp1); 

  fclose(fp1);
} 
/*-------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/
int read_output(char *fname, int read_flag,long int *seed,int *output_flag,int *in_flag,float **rra,float **vva,float *aa)
{
  FILE *fp1;
  long ii;
  float vaa;
  fp1=fopen(fname,"r");
  
  int kk,dummy;
  
  // header reading
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&header1,sizeof(io_header),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  
  
  vaa=(float)header1.time;     //scale factor of  nbody output
  *aa=vaa;
  
  MM=(long)header1.npartTotal[1]; //Total DM particles in this simulation
  DM_m=(float)header1.mass[1]; //DM particle mass in units of 10^10 M_sun/h
  
  vomegam=(float)header1.Omega0;    //Omega_
  vomegalam=(float)header1.OmegaLambda;     //OmegaLambda
  vhh=(float)header1.HubbleParam;     //HubbleParam
  vomegab=(float)header1.Omegab;    //Omega_b
  sigma_8_present=(float)header1.sigma_8_present;
  N1=(long)header1.Nx;
  N2=(long)header1.Ny;
  N3=(long)header1.Nz;
  LL=(float)header1.LL;
  *output_flag=header1.output_flag; // input units ? (!=1) => kp/h; else grid
  *in_flag=header1.in_flag; //  input file generated by ? 1 => zel  else nbody
  *seed=header1.seed;
  
  if(read_flag!=1)
    {
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(&rra[ii][0],sizeof(float),3,fp1);
	  
	  if(*output_flag!=1)
	    {
	      rra[ii][0]=rra[ii][0]/(LL*1000.*vhh);//coordinates in grid
	      rra[ii][1]=rra[ii][1]/(LL*1000.*vhh);
	      rra[ii][2]=rra[ii][2]/(LL*1000.*vhh);
	      
	      rra[ii][0] = rra[ii][0]-1.0*N1*(long)(floor(rra[ii][0])/(1.*N1));
	      rra[ii][1] = rra[ii][1]-1.0*N2*(long)(floor(rra[ii][1])/(1.*N2));  // imposing periodic boundary condition
	      rra[ii][2] = rra[ii][2]-1.0*N3*(long)(floor(rra[ii][2])/(1.*N3));
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
      
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(&vva[ii][0],sizeof(float),3,fp1);
	  
	  if (*output_flag!=1)
	    {
	      vva[ii][0]=vva[ii][0]/(LL*vhh*100./vaa);//velocities in 
	      vva[ii][1]=vva[ii][1]/(LL*vhh*100./vaa);
	      vva[ii][2]=vva[ii][2]/(LL*vhh*100./vaa);
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
    }
  
  fclose(fp1);
}
	  
//*************************************************************************
//              END of   IO  FUNCTIONS
//*************************************************************************


//***********************************************************************************
//                read nbody header and sorted halo catalogue
//***********************************************************************************

void read_fof(char *fname, int read_flag,int *output_flag, long *totcluster, float **halo, float *aa)
{
  FILE *fp1;
  long ii,t;
  float vaa;
  fp1=fopen(fname,"r");
  
  int kk,dummy;
  
  //***********************************************************************************
  //                               header reading
  //***********************************************************************************
  
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&header1,sizeof(io_header),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  

  
  vaa=(float)header1.time;     //scale factor of  nbody output
  *aa=vaa;
  
  MM=(long)header1.npartTotal[1]; //Total DM particles in this simulation
  DM_m=(float)header1.mass[1]; //DM particle mass in units of 10^10 M_sun/h
  
  vomegam=(float)header1.Omega0;    //Omega_
  vomegalam=(float)header1.OmegaLambda;     //OmegaLambda
  vhh=(float)header1.HubbleParam;     //HubbleParam
  vomegab=(float)header1.Omegab;    //Omega_b
  sigma_8_present=(float)header1.sigma_8_present;
  N1=(long)header1.Nx;
  N2=(long)header1.Ny;
  N3=(long)header1.Nz;
  LL=(float)header1.LL;
  *output_flag=header1.output_flag; // input units ? (!=1) => kp/h; else grid
  //*in_flag=header1.in_flag; //  input file generated by ? 1 => zel  else nbody
  
  //***********************************************************************************
  //                        reading total number of haloes
  //***********************************************************************************
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&t,sizeof(long),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  *totcluster=t;
  
  //***********************************************************************************
  //                        reading halo position and velocity
  //***********************************************************************************
  
  if(read_flag!=1)
    {
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0; ii<t; ii++)
	fread(&halo[ii][0],sizeof(float),7,fp1); 
      fread(&dummy,sizeof(dummy),1,fp1);
      
      if(*output_flag!=1)
	for(ii=0; ii<t; ii++)
	  {
	    halo[ii][0]=halo[ii][0]/DM_m;

	    halo[ii][1]=halo[ii][1]/(LL*1000.*vhh);   //coordinates in grid
	    halo[ii][2]=halo[ii][2]/(LL*1000.*vhh);
	    halo[ii][3]=halo[ii][3]/(LL*1000.*vhh);
	    halo[ii][1] = halo[ii][1]-1.0*N1*(long)(floor(halo[ii][1])/(1.*N1));
	    halo[ii][2] = halo[ii][2]-1.0*N2*(long)(floor(halo[ii][2])/(1.*N2));   // imposing periodic boundary condition
	    halo[ii][3] = halo[ii][3]-1.0*N3*(long)(floor(halo[ii][3])/(1.*N3));
	    
	    halo[ii][4]=halo[ii][4]/(LL*vhh*100./vaa);  //velocities
	    halo[ii][5]=halo[ii][5]/(LL*vhh*100./vaa);
	    halo[ii][6]=halo[ii][6]/(LL*vhh*100./vaa);
	  } 
    }
  
  fclose(fp1);
}
//***********************************************************************************
//                  done read nbody header and sorted halo catalogue
//***********************************************************************************
