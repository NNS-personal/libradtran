/************************************************************************
 * $Id: bandec.h 3445 2018-12-12 13:53:59Z bernhard.mayer $
 ************************************************************************/

/*--------------------------------------------------------------------
 * ask Robert Buras
 * ulrike 16.06.2010
 *--------------------------------------------------------------------*/

#ifndef __bandec_h
#define __bandec_h

#if defined (__cplusplus)
extern "C" {
#endif

void bandec5d (double **a, long n, double **al,
	       long indx[]);

void banbks5d (double **a, long n, double **al,
	       long indx[], double b[]);

void bandec(float **a, unsigned long n, int m1, int m2, float **al,
	    unsigned long indx[], float *d);
  
void banbks (float **a, unsigned long n, int m1, int m2, float **al,
	     unsigned long indx[], float b[]);

#endif
