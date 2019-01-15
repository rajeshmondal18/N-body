#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include "power1.h"

/*
  PK(k) is the power spectrum divided by (2 \pi^2) : Pk = A k^n T(k) /(2 pi^2)
  kk is wavenumber (in Mpc^{-1}) and power spectrum is in units (Mpc^3). 
  Notice there are no 'h's in  the definitions. 
*/
extern float norm,vhh;

float Pk(float kk)
{
  float  y,ffb,ffcdm;

  /* Modified on 15.02.2008, using WMAP 7 year data if required please change the value of sigma_8 in the line below. It is presently using sigma_8 = 0.809 */

  //norm=pow(0.809,2.)/(1.424307e-06);
  
  if (kk>1.e-6)
    {
      y=norm*kk*pow(TFfit_onek(kk,&ffb,&ffcdm),2.);
      
    }
  else
    {
      y=0.;
    }

  return(y);
}

float sigma_func(float kk)
{
  float y,R;
  
  R=8./vhh;
  y=3.*(sin(R*kk)-(R*kk)*cos(R*kk))*pow(R*kk,-3);
  y=(kk*kk*Pk(kk)*y*y);
  return(y);
}


float simp(float (*func)(float),float a,float b,int N)
{
	float x,sum1=0.0,sum2=0.0;
	int i;
	float h=(float)((b-a)/N),z;
	
	
	for(x=a+fabs(h);x<=b;x=x+(2*fabs(h)))
	{
		sum1=sum1+(4.0*(*func)(x));
	}
	for(x=a+(2*fabs(h));x<=b;x=x+(2*fabs(h)))
	{
		sum2=sum2+(2.0*(*func)(x));
	}
	z=((sum1+sum2+(*func)(a)+(*func)(b))*fabs(h)/3);
	return(z);
}
