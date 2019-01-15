/*    ----------------------------------------------------
This has various functions - needed in cosmology
 
Df(aa) - calculates normalized D(a) (growing mode) at  scale factor 
aa=1/(1+z) D(1) = 1
ff(aa) - calculates f(a) (dlnD/dlna) at  scale factor aa=1/(1+z)

ONLY works for FRW universe with ordinary matter (vomegam), cosmological 
constant(vomegalam) and curvature(1.0-vomegam-vomegalam).

Also calculates a random number (routine copied from NRC) ran1(long *)
  
/*---------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#define FUNC(x) ((*func)(x))
#define EPS1 1.0e-4 /* Fractional accuracy */
#define JMAX 20

extern float vhh,vomegam,vomegalam;
/* Hubble parameter,Omega-Matter,Cosmological Const */


/*---------------------------------------------------------------------------*/
/*             Does integration using trapezoidal rule                      */
/*---------------------------------------------------------------------------*/

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}


float qtrap(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	int j;
	float s,olds;
	
	float dummy;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS1*fabs(olds)) return s;
		olds=s;
	}
	printf("\nToo many steps in routine qtrap...exiting routine...\n");
	return 0.0; /* Never gets here normally */
}
/*---------------------------------------------------------------------------*/
/*Trapezoidal Rule done */
/*---------------------------------------------------------------------------*/


float Hf(float aa)
{
  return(sqrt(vomegam*pow(aa,-3.) + (1.0-vomegam-vomegalam)*pow(aa,-2.)+ vomegalam));
}

float qf(float aa)
{
  return(0.5*(vomegam*pow(aa,-3.0) - 2.0*vomegalam)/(Hf(aa)*Hf(aa)));
}

float Integrandf(float aa)
{
return(pow(aa*Hf(aa),-3.0));
}

float Integralf(float aa)
{
  /* The Integrand blows up at 0.0 so we take a small step 0.0001 */
return(  qtrap(Integrandf,0.00001,aa)                    );
}

/*---------------------------------------------------------------------------*/
/*Returns the growing mode of denisty perturbation normalized to unity at present */
/*---------------------------------------------------------------------------*/

float Df(float aa)
{
return( Hf(aa)*Integralf(aa)/( Hf(1.0)*Integralf(1.0) ) );
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*Returns the logarithmic derivative of the growing mode of denisty perturbatio */
/*---------------------------------------------------------------------------*/

float ff(float aa)
{
return( (1.0/( aa*aa*pow(Hf(aa),3.0)*Integralf(aa) ) ) - (1.0 + qf(aa) )  );
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*                   function   ran1 copied from Numerical Recipes*/
/*---------------------------------------------------------------------------*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef EPS1
#undef FUNC
#undef JMAX
/*---------------------------------------------------------------------------*/
/*                            ran1 done*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*                   function gasdev copied from Numerical Recipes*/
/*---------------------------------------------------------------------------*/
float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/*---------------------------------------------------------------------------*/
/*                            gasdev done*/
/*---------------------------------------------------------------------------*/
