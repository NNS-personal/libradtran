/************************************************************************
 * $Id: alis.c 3509 2019-11-18 14:15:06Z Claudia.Emde $
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

#include <math.h>
#include <stdlib.h>

#include "alis.h"
#include "ascii.h"

/***********************************************************************************/
/* Function: spectral_is_weight                                           @62_30i@ */
/* Description: Calculate total absorption weight for spectral importance          */
/*              sampling (corresponing to Eq. 14 of ALIS paper (Emde et al., 2011).*/
/*                                                                                 */
/*                                                                                 */
/* Parameters:  Output:                                                            */
/*              totweight_spectral     spectral weight                             */
/*                                                                                 */ 
/*              Input:                                                             */
/*              p                      photon structure of local estimate photon   */
/*              photon                 photon structure of "real" photon           */
/*              atmos                  atmosphere structure                        */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/* Date: 2011-07-20                                                                */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int spectral_is_weight ( double           **totweight_spectral,
			 photon_struct     *p,
			 photon_struct     *photon,
			 atmosphere_struct *atmos )
{
  double static *tau_spectral=NULL;
  int iv=0, kc=0; 
  
  if (tau_spectral==NULL)
    tau_spectral = calloc((size_t) atmos->nlambda_abs, sizeof(double));

  for (iv=0; iv<atmos->nlambda_abs; iv++){
    tau_spectral[iv]=0.0;
    
    for (kc=0; kc<atmos->Nz; kc++)
      /* Spectral absorption optical thickness */
      tau_spectral[iv]+= p->pathlength_per_layer[kc] *
        ( atmos->kabs_spectral[MCCAOTH_TOT][iv][kc] - 
          (atmos->kabs->prof [MCCAOTH_TOT])[kc]
	  + atmos->ksca_spectral[MCCAOTH_TOT][iv][kc] - 
          (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[kc]);
  }
  
  /* XXX RPB: this has nothing to do with NEW_REFLECT, but the correct thing
     here is p instead of photon. However, I did not dare to commit without
     more testing,*/ 
  /* ???? CE: I think only for albedo weight p is needed */
  for (iv=0; iv<atmos->nlambda_abs; iv++)
    (*totweight_spectral)[iv] = exp(-tau_spectral[iv])*photon->q_spectral[iv]*
      photon->q2_spectral[iv] 
      * p->q_albedo_spectral[iv];
  
  return 0; 
}

/***********************************************************************************/
/* Function: concentration_is_weight                                      @62_30i@ */
/* Description: Calculate total absorption weight for spectral importance          */
/*              sampling (corresponing to Eq. 14 of ALIS paper (Emde et al., 2011).*/
/*                                                                                 */
/*                                                                                 */
/* Parameters:  Output:                                                            */
/*              totweight_concentration  concentration weight                      */
/*                                                                                 */ 
/*              Input:                                                             */
/*              p                      photon structure of local estimate photon   */
/*              photon                 photon structure of "real" photon           */
/*              atmos                  atmosphere structure                        */
/* Known bugs:                                                                     */
/* Author: Claudia Emde and Marius Schmidl                                         */
/* Date: 2012-03-06                                                                */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int concentration_is_weight ( double           **totweight_concentration,
			      photon_struct     *p,
			      photon_struct     *photon, 
			      atmosphere_struct *atmos )
{
  double tauabs_concentration=0.0, tausca_concentration=0.0; 
  int ic=0, kc=0; 

  for (ic=0; ic<atmos->Nc; ic++){
    tauabs_concentration=0.0;
    tausca_concentration=0.0;
    
    for (kc=0; kc<atmos->Nz; kc++){
      
      /* Concentration absorption optical thickness */
      tauabs_concentration+= p->pathlength_per_layer[kc] *
        ( - (atmos->kabs->prof [MCCAOTH_AER])[kc] + atmos->kabs_scaled[ic][kc]);
      
      /* Concentration optical thickness, only applied for aerosol scattering */
      tausca_concentration+= p->pathlength_per_layer[kc] *
        (  -(atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL] ->prof [MCCAOTH_AER])[kc]
           + atmos->ksca_scaled[ic][kc]);
    }
    
    
    /* Total weight for concentration importance sampling: */
    
   /*  fprintf(stderr, "ic %d tausca_concentration %g tauabs_concentration %g *photon->q_concentration[ic] %g *photon->q2_concentration[ic] %g\n", ic, tausca_concentration, tauabs_concentration, photon->q_concentration[ic], p->q2_concentration[ic]); */
    
    (*totweight_concentration)[ic] = 
      exp(-tauabs_concentration-tausca_concentration)*photon->q_concentration[ic]*p->q2_concentration[ic];
    
  }
  return 0; 
}

void write_jacobian_result(sample_struct* sample,
                           atmosphere_struct* atmos, 
                           jacobian_result_field* jacobian,
                           FILE* backjacfile,
                           int is,
                           int js)
{
  int kc=0, isp=0;

  for (kc=0; kc<atmos->Nz; kc++) {
    fprintf (backjacfile, "%.8e %.8e %.6e",
             ((double) is + 0.5) * sample->delX,
             ((double) js + 0.5) * sample->delY,
             atmos->Z[kc]);

    for (isp=1; isp<=jacobian->n_caoth; isp++)
      
      fprintf (backjacfile, " %.6e %.6e",
               jacobian->jacobian_t[is][js][isp][0][kc],
               jacobian->jacobian_t[is][js][isp][1][kc]);
              
    fprintf (backjacfile, "\n");
    
  }
  
}

/***********************************************************************************/
/* Function: calloc_jacobian_result_field                                 @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct jacobian_result_field and initialize structure.     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

jacobian_result_field *calloc_jacobian_result_field (int std, int Nx, int Ny, int Nz, int n_caoth)
{
  int ijac=0, isp=0, is=0, js=0;

  jacobian_result_field *temp = calloc(1, sizeof(jacobian_result_field));

  if (temp==NULL) return NULL;

  temp->Nx=Nx;
  temp->Ny=Ny;
  temp->Nz_jac=Nz;
  temp->n_caoth=n_caoth;

  if ((
       temp->jacobian_t = calloc ((size_t) Nx, sizeof(double ****))
       )==NULL ) return NULL;
  

  for (is=0; is<Nx; is++) {
    if ((
	 temp->jacobian_t[is] = calloc ((size_t) Ny, sizeof(double ***))
	 )==NULL ) return NULL;

    for (js=0; js<Ny; js++) {
      if ((
	   temp->jacobian_t[is][js] = calloc ((size_t) n_caoth+1, sizeof(double **))
	   )==NULL ) return NULL;

      for (isp=1; isp<=n_caoth; isp++) {
      if ((
	   temp->jacobian_t[is][js][isp] = calloc ((size_t) 2, sizeof(double *))
	   )==NULL ) return NULL;

	for (ijac=0; ijac<2; ijac++) {
	  if ((
	       temp->jacobian_t[is][js][isp][ijac] = calloc ((size_t) Nz, sizeof(double))
	       )==NULL ) return NULL;
	}
      }
    }
  }

  if (std) {
    if ((
	 temp->jacobian_t2 = calloc ((size_t) Nx, sizeof(double ****))
	 )==NULL ) return NULL;
  

    for (is=0; is<Nx; is++) {
      if ((
	   temp->jacobian_t2[is] = calloc ((size_t) Ny, sizeof(double ***))
	   )==NULL ) return NULL;

      for (js=0; js<Ny; js++) {
	if ((
	     temp->jacobian_t2[is][js] = calloc ((size_t) n_caoth+1, sizeof(double **))
	     )==NULL ) return NULL;

	for (isp=1; isp<=n_caoth; isp++) {
	  if ((
	       temp->jacobian_t2[is][js][isp] = calloc ((size_t) 2, sizeof(double *))
	       )==NULL ) return NULL;

	  for (ijac=0; ijac<2; ijac++) {
	    if ((
		 temp->jacobian_t2[is][js][isp][ijac] = calloc ((size_t) Nz, sizeof(double))
		 )==NULL ) return NULL;
	  }
	}
      }
    }

    if ((
	 temp->jacobian_tact = calloc ((size_t) n_caoth+1, sizeof(double **))
	 )==NULL ) return NULL;

    for (isp=1; isp<=n_caoth; isp++) {
      if ((
	   temp->jacobian_tact[isp] = calloc ((size_t) 2, sizeof(double *))
	   )==NULL ) return NULL;

      for (ijac=0; ijac<2; ijac++) {
	if ((
	     temp->jacobian_tact[isp][ijac] = calloc ((size_t) Nz, sizeof(double))
	     )==NULL ) return NULL;
      }
    }
  }

  return temp;
}

/***********************************************************************************/
/* Function: free_jacobian_result_field                                   @62_30i@ */
/* Description:                                                                    */
/*  Free memory for struct jacobian_result_field.                                  */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void free_jacobian_result_field (jacobian_result_field *res)
{
  int isp=0, ijac=0, is=0, js=0;

  if (res->jacobian_t != NULL) {
    for (is=0; is<res->Nx; is++) {
      for (js=0; js<res->Ny; js++) {
	for (isp=1; isp<=res->n_caoth; isp++) {
	  for (ijac=0; ijac<2; ijac++)
	    free(res->jacobian_t[is][js][isp][ijac]);
	  free(res->jacobian_t[is][js][isp]);
	}
	free(res->jacobian_t[is][js]);
      }
      free(res->jacobian_t[is]);
    }
    free(res->jacobian_t);
  }

  if (res->jacobian_t2 != NULL) {
    for (is=0; is<res->Nx; is++) {
      for (js=0; js<res->Ny; js++) {
	for (isp=1; isp<=res->n_caoth; isp++) {
	  for (ijac=0; ijac<2; ijac++)
	    free(res->jacobian_t2[is][js][isp][ijac]);
	  free(res->jacobian_t2[is][js][isp]);
	}
	free(res->jacobian_t2[is][js]);
      }
      free(res->jacobian_t2[is]);
    }
    free(res->jacobian_t2);

    for (isp=1; isp<=res->n_caoth; isp++) {
      for (ijac=0; ijac<2; ijac++)
	free(res->jacobian_tact[isp][ijac]);
      free(res->jacobian_tact[isp]);
    }
    free(res->jacobian_tact);
  }

  free(res);
}

/***********************************************************************************/
/* Function: calloc_photon_jacobian_part                                     @62_30i@ */
/* Description: Allocate photon, the part for Jacobians or lidar                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void calloc_photon_jacobian_part ( photon_struct *p,
				   sample_struct *sample,
				   int            Nz,
				   int            source,
				   int            n_caoth )
{
  int j=0, isp=0;

  /* allocate memory for jacobian matrix evaluation */
  /* this is quite expensive!!! only needed in generate_photon, can we do this better? */
  /* increases computational time relative to pure MC by 10% */

  /* BCA: introduce logical FIRST, only allocate last dimension again and again;
     or simlpy set to zero in generate_photon */

  if (sample->LLE_jacobian || sample->abs_jacobian) {

    p->Nz_jac = Nz;
    p->n_caoth = n_caoth;

    p->q_jacobian = calloc ((size_t) 2, sizeof(double **));

    for (j=0; j<2; j++)
      p->q_jacobian[j] = calloc ((size_t) n_caoth, sizeof(double *));

    for (j=0; j<2; j++)
      for (isp=0; isp<n_caoth; isp++)
	p->q_jacobian[j][isp] = calloc ((size_t) Nz, sizeof(double));

    /* CE: Check whether q_jacobian could also be used instead of new variable */
    ASCII_calloc_double(&p->q_jacobian_sca, n_caoth, Nz);
    
    p->r_jacobian = calloc ((size_t) Nz, sizeof(double));

    p->lest.pdir_sct = calloc((size_t) n_caoth+1, sizeof(double *));
    for (isp=0; isp<=n_caoth; isp++)
      p->lest.pdir_sct[isp] = calloc((size_t) sample->nstokes, sizeof(double));
  }
  else {
    p->q_jacobian = NULL;
    p->r_jacobian = NULL;
    p->lest.pdir_sct = NULL;
  }
}
