/************************************************************************
 * $Id: mystic_3d.h 3351 2018-02-21 15:00:33Z tobias.koelling $
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

#ifndef _MYSTIC_3D_H
#define _MYSTIC_3D_H 1

#include "mystic.h"

#if defined (__cplusplus)
extern "C" {
#endif

int step3D ( photon_struct     *p,
	     atmosphere_struct *atmos,
	     double             step,
	     int                bcond,
	     int                photonpath,
	     int                spherical3D,
	     int                visualize );

int intersection3D ( photon_struct     *p,
		     atmosphere_struct *atmos,
		     double             tau,
		     double             tausca,
		     double            *step);

#if defined (__cplusplus)
}
#endif

#endif /* _MYSTIC_3D_H */

