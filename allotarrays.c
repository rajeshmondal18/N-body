#include<stdlib.h>
#include<fftw3.h>

float  ***allocate_fftwf_3d(long N1,long N2,long N3)
{
  long ii,jj;
  long asize,index;
  float ***phia, *phi;

  phia=(float ***) fftwf_malloc (N1 *  sizeof(float **));


  for(ii=0;ii<N1;++ii)
      phia[ii]=(float **) fftwf_malloc (N2 *  sizeof(float *));

  asize = N1*N2;
  asize = asize*N3;

  if(!(phi = (float *) calloc(asize,sizeof(float))))
    {
      printf("error in allocate_fftwf_3d");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    for(jj=0;jj<N2;++jj)
      {
	index = N2*N3;
	index = index*ii + N3*jj;
	phia[ii][jj]=phi+ index;
      }
  return(phia);
}

float **allocate_float_2d(long N1,int N2) /* Typically N1 is number of particles and N2 is number of dimensions namely 3 */
{
  float **xxa, *xx;
  long ii;

  xxa=(float**)malloc(N1 *  sizeof(float*));
  if(!(xx = (float *) calloc((size_t)(N1*N2),sizeof(float))))
    {
      printf("error in allocate_float_2d");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    xxa[ii]=xx + N2*ii ;

return(xxa);
}
