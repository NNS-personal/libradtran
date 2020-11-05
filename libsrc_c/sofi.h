/************************************************************************
 * $Id: sofi.h 3432 2018-10-09 14:18:27Z Claudia.Emde $
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

#ifndef __sofi_h
#define __sofi_h
#include "mystic.h"

int sample_photons_sofi(float lambda, double dx, double ratio, double sdist,
			int Nx,int limb, double *pd);

void generate_photon_sofi ( atmosphere_struct* atmos, photon_struct* p, 
                            double sza, double phi, double* pd );

#endif
