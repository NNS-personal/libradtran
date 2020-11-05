/*--------------------------------------------------------------------
 * $Id: molecular3d.c 3451 2019-01-28 00:28:05Z bernhard.mayer $ 
 *
 * This file is part of libRadtran.
 * Copyright (c) 1997-2017 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>                                                          
#include <stdlib.h>
#include <string.h>
#include "uvspec.h"
#include "ascii.h"
#include "redistribute.h"
#include "ancillary.h"
#include "allocnd.h"


#if HAVE_LIBNETCDF
#include <netcdf.h>
#define RUN_NC(command, msg) do { \
    int status = (command); \
    if (status != NC_NOERR) { \
      fprintf(stderr, "netCDF error %d: \"%s\" %s from %s\n", status, nc_strerror(status), (msg), filename); \
      return status; \
    } \
} while(0)
#endif


#ifndef MAX_LENGTH_OF_LINE
#define MAX_LENGTH_OF_LINE     65536
#endif

/************************************/
/* prototypes of internal functions */
/************************************/
static int read_atmosphere_3d (char *filename, atm_out_struct *out,
                               float* mol_mass, int interpol_method_temper, int *interpol_method_gas,
                               int quiet, int verbose);

static int write_dummy_profile_file(char *filename, 
				    int quiet, int verbose);

static int redistribute_3D_temp (float ****data, int nx, int ny,
                                 float *z_old, int nlyr_old, 
                                 float *z_new, int nlyr_new, double *temper_1D); 


int setup_molecular3d(input_struct input, output_struct *output)
{
  int status=0; 
  
  
  if(!input.atmosphere3d)
    return status;
  
  else{
    /* Setup 1D calculations for trace gases that are not 3D */
    /* must be set to 0, because first 1D atmospheric profiles are set up */
    output->molecular3d=0;
    
    /* set dens column of included species to 0 for 1D calculations, should be specified in input file, first consider only H2O,
     use 1e-8 because 0.0 produces nan in Fu band 18.*/
    if (input.ck_scheme == CK_FU)
      input.atm.column[MOL_H2O]=1.e-8;
    else
      input.atm.column[MOL_H2O]=0.0;
    
    /* write dummy profile file */
    status = write_dummy_profile_file(input.atmosphere3d_filename,
				      input.quiet, input.verbose); 
    
    if (status!=0) {
      fprintf (stderr, "Error %d creating dummy profile file from %s\n",
	       status, input.atmosphere3d_filename);
      return status;
    }
  }
  
  return status; 
  }


int setup_optprop_molecular3d(input_struct input, output_struct *output){

  int status = 0;
  char function_name[]="setup_optprop_molecular3d";
  char file_name[]="molecular3d.c";
  int ix=0, iy=0, lc=0, iv=0, iq=0;
  double ***tmp_molabs=NULL, ***tmp_rayleigh=NULL;
  float *tmp_zd_common=NULL, *tmp_zd=NULL; 
  int tmp_nlev_common=0, tmp_nlev=0; 
  double *tmp_temper=NULL;   
  
  if(!input.atmosphere3d)
    return status;
  
  else{
    output->molecular3d=1;

    /* We want to use the same functions for the calculation of
       molecular absorption for 1D and 3d, therefore we need to save
       the variables output->atm.optprop.tau_molabs_r,
       output->atm.optprop.tau_rayleigh_r . These are reserved for 1D
       optical properties and temporarily also used for 3D optical
       properties. */
    
    /* copy 1D molecular optical properties to temporary arrays */
    ASCII_calloc_double_3D_arylen_restricted (&tmp_molabs, output->atm.nlev_common-1, output->wl.nlambda_r,output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r); 
    ASCII_calloc_double_3D_arylen_restricted (&tmp_rayleigh, output->atm.nlev_common-1, output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r);
    
    for (lc=0; lc< output->atm.nlev_common-1; lc++){
      for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++){
        for(iq=0; iq<output->atm.nq_r[iv]; iq++){
          tmp_molabs[lc][iv][iq]=output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq];
          tmp_rayleigh[lc][iv][iq]=output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq];
	  //	  fprintf (stdout, "R1 %d %e\n", lc, tmp_rayleigh[lc][iv][iq]);
	}
      }
    }

    /* We assume that z-grid of the 3d molecular atmosphere is the
       same as that of the "profile_files", therefore the common
       z-grid should not be changed after including 3D molecular
       properties. Variables are changed in functions below, so we
       need to store common z-grid temporarily. */

    tmp_nlev=output->atm.nlev;
    tmp_zd = calloc ( output->atm.nlev, sizeof (float));
    for (lc=0; lc< output->atm.nlev; lc++){
       tmp_zd[lc]=output->atm.zd[lc];
    }
    
    tmp_nlev_common=output->atm.nlev_common;
    tmp_zd_common = calloc ( output->atm.nlev_common, sizeof (float));
    tmp_temper = calloc ( output->atm.nlev_common, sizeof (double));
    for (lc=0; lc< output->atm.nlev_common; lc++)
      tmp_zd_common[lc]=output->atm.zd_common[lc];
    
    for (lc=0; lc< output->atm.nlev_common; lc++)
      tmp_temper[lc]=output->atm.microphys.temper[0][0][lc];

    /* Calculate 3D molecular optical properties */
    status = read_atmosphere_3d (input.atmosphere3d_filename,
                                 &(output->atm), 
                                 input.atm.mol_mass, input.atm.interpol_method_temper, input.atm.interpol_method_gas,
                                 input.quiet, input.verbose);
    if (status!=0) {
      fprintf (stderr, "Error %d reading 3D atmosphere from %s\n",
               status, input.atmosphere3d_filename);
      return status;
    }

    /*This is needed if we have 3D O2 distributions (can probably be removed.) */
    /* BM: If you allow horizontal pressure and temperature variations, you also need to allow horizontal O2 variations :-) */
    /* for (ix=0; ix<output->atm.Nxatm; ix++) { */
    /*   for (iy=0; iy<output->atm.Nyatm; iy++) { */
    /*  for (lc=0; lc<output->atm.nlev; lc++) { */
    /*    fact = output->atm.microphys.dens[MOL_O2][ix][iy][lc]* 1e-23; /\* O4 cross section is scaled by 1e46 due *\/ */
    /*    /\* to float precision limitations         *\/ */
    /*    /\* account for that here. ak 20110404     *\/ */
    /*    output->atm.microphys.dens[MOL_O4][ix][iy][lc] = fact*fact; */
    /*  } */
    /*   } */
    /* } */
    


    /**** Read cross section files, like O3, BrO, OclO, NO2, etc. ****/
    pmesg (" ... calling setup_crs(), reading cross section files\n", input.verbose);
    status = setup_crs (input, output); /* in molecular.c */
    if (status!=0) {
      fprintf (stderr, "Error %d setting up absorption cross sections for 3D atmosphere in %s (%s)\n", status, function_name, file_name);
      return status;
    }
    


    /**** Calculate Rayleigh scattering cross section ****/
    pmesg (" ... calling setup_rayleigh(), calculating Rayleigh scattering for 3D atmosphere \n", input.verbose);
    status = setup_rayleigh (input, output); /* in molecular.c */ 
    if (status!=0) {
      fprintf (stderr, "Error %d setting up Rayleigh cross sections for 3D atmosphere in %s (%s)\n", status, function_name, file_name);
      return status;
    }

    /****  Optical properties of trace gases ****/
    pmesg (" ... calling setup_gases(), generating optical properties of trace gases for 3D atmosphere \n", input.verbose);
    status = setup_gases (input, output);  /* in atmosphere.c */
    if (status!=0) {
      fprintf (stderr, "Error %d setting up optical properties of trace gases for 3D atmosphere in %s (%s)\n", status, function_name, file_name);
      return status;
    }


    /**** Redistribute 3D molecular profiles to required vertical resolution ****/
    pmesg (" ... redistributing 3D molecular optical properties\n", input.verbose);

    /* common z_grid */
    free(output->atm.zd_common);
    /* ASCII_free_float_3D (output->atm.microphys.temper, output->atm.Nxatm, output->atm.Nyatm); */
    /* ASCII_calloc_float_3D (&output->atm.microphys.temper, 1, 1, tmp_nlev_common); */
    output->atm.zd_common=calloc(tmp_nlev_common, sizeof(float));
    output->atm.nlev_common=tmp_nlev_common;
    output->atm.nlyr_common=tmp_nlev_common-1;
    
    for (lc=0; lc< output->atm.nlev_common; lc++)
      output->atm.zd_common[lc]=tmp_zd_common[lc];
        
    /* for (lc=0; lc< output->atm.nlev_common; lc++) */
    /*   output->atm.microphys.temper[0][0][lc]=tmp_temper[lc]; */

    redistribute_molecular(&(output->atm.optprop), output->wl.nlambda_r, 
                           output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
                           output->atm.nq_r,
                           output->nipa,
                           output->atm.zd, output->atm.nlev-1,
                           output->atm.Nxatm, output->atm.Nyatm,
                           output->atm.zd_common, output->atm.nlev_common-1);
    
    
    /* redistribute temperature */
    redistribute_3D_temp(&(output->atm.microphys.temper), output->atm.Nxatm, output->atm.Nyatm, output->atm.zd, output->atm.nlev, output->atm.zd_common, output->atm.nlev_common, tmp_temper); 
    /* for (ix=0; ix<output->atm.Nxatm; ix++) */
    /*   for (iy=0; iy<output->atm.Nyatm; iy++) */
    /*  for (lc=0; lc< output->atm.nlev_common; lc++) */
    /*    fprintf(stderr, "T after redist ix %d iy %d lc %d z %.2f T %.2f \n", ix, iy, lc, output->atm.zd_common[lc], output->atm.microphys.temper[ix][iy][lc]); */

    
    ASCII_free_float_3D (output->atm.microphys.temper_avg, output->atm.Nxatm, output->atm.Nyatm);
    ASCII_calloc_float_3D (&output->atm.microphys.temper_avg, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common-1);
    
    for (ix=0; ix<output->atm.Nxatm; ix++)
      for (iy=0; iy<output->atm.Nyatm; iy++)
        for (lc=0; lc< output->atm.nlev_common-1; lc++)
          output->atm.microphys.temper_avg[ix][iy][lc]=0.5*(output->atm.microphys.temper[ix][iy][lc]+
                                                         output->atm.microphys.temper[ix][iy][lc+1]);
    
    

    /* copy optical properties of 3D species to mc structure. */
    /*
    ASCII_calloc_float_5D_arylen_restricted (&output->mc.kabs3D, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common-1, output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r);
    ASCII_calloc_float_5D_arylen_restricted (&output->mc.ksca3D, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common-1, output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r);
    */

    if(!(output->mc.kabs3D = calloc_float_5D_restricted ((size_t) output->atm.Nxatm, (size_t) output->atm.Nyatm, (size_t) output->atm.nlev_common-1,
							 (size_t) output->wl.nlambda_rte_lower, (size_t) output->wl.nlambda_rte_upper, output->atm.nq_r, "output->mc.kabs3D"))) return -1;

    if(!(output->mc.ksca3D = calloc_float_5D_restricted ((size_t) output->atm.Nxatm, (size_t) output->atm.Nyatm, (size_t) output->atm.nlev_common-1,
							 (size_t) output->wl.nlambda_rte_lower, (size_t) output->wl.nlambda_rte_upper, output->atm.nq_r, "output->mc.ksca3D"))) return -1;




    for (ix=0; ix<output->atm.Nxatm; ix++){
      for (iy=0; iy<output->atm.Nyatm; iy++){
        for (lc=0; lc< output->atm.nlev_common-1; lc++){ 
          for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++){
            for(iq=0; iq<output->atm.nq_r[iv]; iq++){
	      
	      /* BM: memory is doubled here (only temporarily, but anyway): output->mc.kabs3D is exactly the same as output->atm.optprop.tau_molabs_r */
	      /* does that make sense, considering that this might be memory critical?                                                                */

	      /* BM: careful! kabs3D and ksca3D are actually absorption and scattering coefficients while tau_molabs and tau_rayleigh */
	      /* are actually optical thicknesses. Need to divide by delta_z to convert!                                              */

	      output->mc.kabs3D[ix][iy][lc][iv][iq]=output->atm.optprop.tau_molabs_r[ix][iy][lc][iv][iq]/(output->atm.zd_common[lc] - output->atm.zd_common[lc+1]);

	      output->mc.ksca3D[ix][iy][lc][iv][iq]=output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq]/(output->atm.zd_common[lc] - output->atm.zd_common[lc+1]);

	      
	      /* fprintf(stderr, "ix %d iy %d lc %d iv %d iq %d  3doptprop z %.2f  tauabs %g tausca %g kabs %g ksca %g dz %g \n", */
	      /* 	      ix, iy, lc, iv, iq, output->atm.zd_common[lc], output->atm.optprop.tau_molabs_r[ix][iy][lc][iv][iq], output->atm.optprop.tau_rayleigh_r[ix][iy][lc][iv][iq],  output->mc.kabs3D[ix][iy][lc][iv][iq], output->mc.ksca3D[ix][iy][lc][iv][iq],output->atm.zd[lc] - output->atm.zd[lc+1]  ); */
	      
	      /* BM: optionally set 3D Rayleigh to zero */
	      // output->mc.ksca3D[ix][iy][lc][iv][iq]=1e-15;   // Rayleigh 1
            }
          }
        }
      }
    }
    
    /* cross sections not needed anymore */
    /*
    ASCII_free_float_5D_restricted (output->atm.optprop.tau_molabs_r, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common-1,
				    output->wl.nlambda_r,output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper); 
    ASCII_free_float_5D_restricted (output->atm.optprop.tau_rayleigh_r, output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev_common-1,
				    output->wl.nlambda_r,output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper);
    */

    free_float_5D_restricted (output->atm.optprop.tau_molabs_r, output->wl.nlambda_rte_lower);
    free_float_5D_restricted (output->atm.optprop.tau_rayleigh_r, output->wl.nlambda_rte_lower);
    

    switch(input.ck_scheme) {
    /* ??? this is highly experimental: BM 7.9. 2018; only element (0,0,0,MOL_AIR,iv) is used below, hence we don't need Nx x Ny array ??? */
    case CK_FU:
    case CK_FILE:          /* ??? need to check if that is appropriate here ??? */
    case CK_AVHRR_KRATZ:   /* ??? need to check if this is appropriate here ??? */
    case CK_LOWTRAN:       /* ??? need to check if this is appropriate here ??? */
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
      break;
      
    /* ??? all other spectral schemes unchanged, as introduced by CE ??? */
    default: 
      
      ASCII_free_float_5D(output->crs.crs, 
			  output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev, MOL_NN);
      break;
    }


    
    /* copy back 1D molecular optical properties. */
    /*
    ASCII_calloc_float_5D_arylen_restricted (&output->atm.optprop.tau_molabs_r, 1, 1, output->atm.nlev_common-1, output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,output->atm.nq_r); 
    ASCII_calloc_float_5D_arylen_restricted (&output->atm.optprop.tau_rayleigh_r, 1, 1, output->atm.nlev_common-1, output->wl.nlambda_r, output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r);
    */
    
    if(!(output->atm.optprop.tau_molabs_r = calloc_float_5D_restricted (1, 1, (size_t) output->atm.nlev_common-1,
									(size_t) output->wl.nlambda_rte_lower, (size_t) output->wl.nlambda_rte_upper, output->atm.nq_r, "output->atm.optprop.tau_molabs_r"))) return -1;

    if(!(output->atm.optprop.tau_rayleigh_r = calloc_float_5D_restricted (1, 1, (size_t) output->atm.nlev_common-1,
									  (size_t) output->wl.nlambda_rte_lower, (size_t) output->wl.nlambda_rte_upper, output->atm.nq_r, "output->atm.optprop.tau_rayleigh_r"))) return -1;

    

    //fprintf(stderr, "dim alloc nlyr %d nlam %d nlamlow %d nlamup %d nq[0] %d\n", output->atm.nlev_common-1, output->wl.nlambda_r,output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper, output->atm.nq_r[0]);
    
    free(output->atm.zd);
    output->atm.zd=calloc(tmp_nlev, sizeof(float));
    output->atm.nlev=tmp_nlev;
    output->atm.nlyr=tmp_nlev-1;
    for (lc=0; lc< output->atm.nlev; lc++)
      output->atm.zd[lc]=tmp_zd[lc];
    
    for (lc=0; lc< output->atm.nlev_common-1; lc++){
      for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++){
        for(iq=0; iq<output->atm.nq_r[iv]; iq++){
          /* fprintf(stderr, "lc %d iv %d iq %d z %.2f abs %g sca %g \n", lc, iv, iq, output->atm.zd_common[lc], tmp_molabs[lc][iv][iq], tmp_rayleigh[lc][iv][iq]);  */
          output->atm.optprop.tau_molabs_r[0][0][lc][iv][iq]=tmp_molabs[lc][iv][iq];
          output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq]=tmp_rayleigh[lc][iv][iq];

	  /* ??? BM: quick and dirty hack: set 1D Rayleigh scattering to 0; need to change that for cases where 1D and 3D atmospheres are combined ??? */
	  output->atm.optprop.tau_rayleigh_r[0][0][lc][iv][iq]=0;  // Rayleigh 2
        }
      }
    }

    /* free temporary variables */
    
    ASCII_free_double_3D (tmp_molabs, output->atm.nlev_common-1, output->wl.nlambda_r);
    ASCII_free_double_3D (tmp_rayleigh, output->atm.nlev_common-1, output->wl.nlambda_r);
    
    free(tmp_zd_common);
    free(tmp_zd);
    free(tmp_temper);

    return status; 
  }
}


/***********************************************************************************/
/* Function: optical_properties_atmosphere3D                                       */
/* Description:                                                                    */
/*  Calculates 3D Rayleigh scattering coefficient and 3D gas absorption            */
/*  coefficient. Only molecules are included in this function, other atmospheric   */
/*  constituents should be included with options ic_file 3D, wc_file 3D,           */
/*  and profile_file 3D.                                                           */
/*                                                                                 */
/* Author: Claudia Emde                                                            */
/* Date:   May 2017                                                                */
/*                                                                                 */
/***********************************************************************************/
int optical_properties_molecular3d (input_struct   input,
                                    output_struct *output,
                                    caoth3d_out_struct *caoth3d,
                                    int            iv,
                                    int            iq)
{
  int ix=0, iy=0, lc=0, start_lc=0;
  int isp = 0, isp_mol=0; 
  int nlev = output->atm.nlev_common;
  int nlyr = nlev-1;
  int Nx = output->atm.Nxatm;
  int Ny = output->atm.Nyatm;  
  //double ***tmp_kabs3D;
  //double ***tmp_ksca3D; 
  
  
  //ASCII_calloc_double_3D (&tmp_kabs3D, Nx, Ny, nlyr);
  //ASCII_calloc_double_3D (&tmp_ksca3D, Nx, Ny, nlyr);
  
  //output->mc.z       = (float *) calloc (nlev, sizeof(float));
  ASCII_calloc_float_3D (&output->mc.temper, Nx, Ny, nlyr+1);
  
  /* ??? refind should also be 3D */
  //output->mc.refind  = (float *) calloc (nlev, sizeof(float));

  //fprintf(stderr, "3DAbs optical_properties_atmosphere3D nlyr %d\n", nlyr); 
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      /* for (lc=0; lc<nlyr; lc++) { */
        
      /*   /\* cross section depends on the subband                  *\/ */
      /*   switch(input.ck_scheme) { */
      /*   case CK_FU: */
      /*     /\* only in the Fu and Liou case, the Rayleigh scattering *\/ */
      /*     /\* cross section depends on the subband                  *\/ */
      /*     tmp_ksca3D[ix][iy][lc] = output->mc.ksca3D[ix][iy][lc][iv][iq]; */
      /*     break; */
          
      /*   case CK_KATO: */
      /*   case CK_KATO2: */
      /*   case CK_KATO2_96: */
      /*   case CK_KATO2ANDWANDJI: */
      /*   case CK_AVHRR_KRATZ: */
      /*   case CK_FILE: */
      /*   case CK_LOWTRAN: */
      /*   case CK_CRS: */
      /*   case CK_REPTRAN: */
      /*   case CK_REPTRAN_CHANNEL: */
      /*     tmp_ksca3D[ix][iy][lc] = output->mc.ksca3D[ix][iy][lc][iv][0]; */
      /*     break; */
          
      /*   default: */
      /*     fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme); */
      /*     return -1; */
          
      /*     break; */
      /*   } */
        
      /*   switch (output->atm.molabs) { */
      /*   case MOLABS_CALC: */
      /*   case MOLABS_LOOKUP: */
      /*   case MOLABS_FILE_MONO: */
      /*   case MOLABS_FILE_SPEC: */
      /*     tmp_kabs3D [ix][iy][lc] = output->mc.kabs3D[ix][iy][lc][iv][iq]; */
      /*     break; */
          
      /*   case MOLABS_NONE: */
      /*     tmp_kabs3D [ix][iy][lc] = 0; */
      /*     break; */
          
      /*   default: */
      /*     fprintf (stderr, "Error, unknown molecular absorption option %d\n",  */
      /*              output->atm.molabs); */
      /*     return -1; */
      /*   } */
      /* } */
      for (lc=0; lc<nlyr+1; lc++) {
	if (output->atm.molabs ==MOLABS_NONE )
	  output->mc.kabs3D[ix][iy][lc][iv][iq]=0.0;
	
        output->mc.temper[ix][iy][nlyr-lc] = output->atm.microphys.temper[ix][iy][lc];
      }
    }
  }
  
  /* for (lc=0; lc<=nlyr; lc++) { */
  /*   output->mc.z     [nlyr-lc] = output->atm.zd[lc]; */
  /* } */
  
  // ext and ssa in caoth3d are only allocated for 3D layers, thus we need different indices here 
  
  for (isp=0; isp<input.n_caoth; isp++)
    if (strcmp(output->caoth3d[isp].name, "molecular_3d") ==0)
      isp_mol=isp; 

  //fprintf(stderr, "name %s input.n_caoth %d isp_mol %d\n", output->caoth3d[isp_mol].name, input.n_caoth, isp_mol);
   
  start_lc = nlyr - output->caoth3d[isp_mol].nthreed;  

  //fprintf(stderr, "3DAbs start_lc %d nthreed %d \n", start_lc, output->caoth3d[isp_mol].nthreed);
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      for (lc=start_lc; lc<nlyr; lc++){
        /* fprintf(stderr, "ix %d iy %d lc %d z %.2f ksca %.3g  kabs %.3g \n",ix, iy, lc, output->atm.zd_common[lc], tmp_ksca3D[ix][iy][lc], tmp_kabs3D[ix][iy][lc]); */

	switch(input.ck_scheme) {
        case CK_FU:
          /* only in the Fu and Liou case, the Rayleigh scattering */
          /* cross section depends on the subband                  */
	  output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy]=0.001*(output->mc.ksca3D[ix][iy][lc][iv][iq]+output->mc.kabs3D[ix][iy][lc][iv][iq]);
	  if (output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy] > 0)
	    output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0.001*output->mc.ksca3D[ix][iy][lc][iv][iq]/output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy];
	  else   /* avoid dividing 0/0 in case of no extinction */
	    output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0;
	  break;
	case CK_KATO:
        case CK_KATO2:
        case CK_KATO2_96:
        case CK_KATO2ANDWANDJI:
        case CK_AVHRR_KRATZ:
        case CK_FILE:
        case CK_LOWTRAN:
        case CK_CRS:
        case CK_REPTRAN:
        case CK_REPTRAN_CHANNEL:
	  output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy]=0.001*(output->mc.ksca3D[ix][iy][lc][iv][0]+output->mc.kabs3D[ix][iy][lc][iv][iq]);
	  if (output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy] > 0)
	    output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0.001*output->mc.ksca3D[ix][iy][lc][iv][0]/output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy];
	  else   /* avoid dividing 0/0 in case of no extinction */
	    output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0;
	  break; 
	default:
          fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", input.ck_scheme);
          return -1;
          
          break;
        }
      
       
	/* output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy]=0.001*(tmp_ksca3D[ix][iy][lc]+tmp_kabs3D[ix][iy][lc]); */
	/* if (output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy] > 0) */
	/*   output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0.001*tmp_ksca3D[ix][iy][lc]/output->caoth3d[isp_mol].ext[nlyr-lc-1][ix][iy]; */
	/* else   /\* avoid dividing 0/0 in case of no extinction *\/ */
	/*   output->caoth3d[isp_mol].ssa[nlyr-lc-1][ix][iy]=0; */
      }
    }
  }

  /* free temporary variables */
  //ASCII_free_double_3D (tmp_kabs3D, Nx, Ny);
  //ASCII_free_double_3D (tmp_ksca3D, Nx, Ny);

  return 0; 
}

 

/*********************************************************/
/* Function: read_atmosphere_3d                          */
/* Description:                                          */
/* Read 3D atmosphere file and write to atm_out_struct.  */
/*                                                       */ 
/* Autor:                                                */
/* March 2017    C. Emde   Created                       */
/*********************************************************/
static int read_atmosphere_3d (char *filename, atm_out_struct *out,
                               float* mol_mass, int interpol_method_temper, int *interpol_method_gas,
                               int quiet, int verbose)
{
  int i=0, number=0, status=0;
  int ix=0, iy=0, iz=0;

  int rows=0, min_columns=0, max_columns=0, max_length=0;

  char line[MAX_LENGTH_OF_LINE+1]="";
  char *string=NULL;

  FILE *file=NULL;;
  char **array=NULL;
  char *dummy=NULL;
  
  if (!quiet) 
    fprintf (stderr, " ... reading 3D atmosphere file %s\n", filename);

#if HAVE_LIBNETCDF
  int ixstart=0, iystart=0;
  size_t N;

  size_t n=0;
  int read_netcdf=0;

  int    ncid   =0, idd_nx =0, idd_ny =0, idd_nz=0;
  int    id_z=0, id_press=0, id_temp=0, id_qH2O=0;
  
  float ***tmpdata=NULL;

#endif


  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {
    read_netcdf=1;

    if (!quiet)
      fprintf (stderr, " ... reading 3D atmosphere data from netCDF file %s\n", filename);

    RUN_NC(nc_inq_dimid(ncid, "nx", &idd_nx), "reading nx");
    RUN_NC(nc_inq_dimlen(ncid, idd_nx, &n), "reading nx");
    out->Nxatm=n;

    RUN_NC(nc_inq_dimid(ncid, "ny", &idd_ny), "reading ny");
    RUN_NC(nc_inq_dimlen(ncid, idd_ny, &n), "reading ny");
    out->Nyatm=n;

    RUN_NC(nc_inq_dimid(ncid, "nz_lev", &idd_nz), "reading nz");
    RUN_NC(nc_inq_dimlen(ncid, idd_nz, &n), "reading nz");

    out->nlev = n;
    out->nlyr = out->nlev-1;
    out->nlev_common = n;
    out->nlyr_common = out->nlev_common-1;
    
    RUN_NC(nc_get_att_double(ncid, NC_GLOBAL, "dx", &(out->dxatm)), "reading dx");  
    RUN_NC(nc_get_att_double(ncid, NC_GLOBAL, "dy", &(out->dyatm)), "reading dy");
  }
  else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", 
	       filename);
#endif
  
    /* check input file */
    status =  ASCII_checkfile (filename, 
			       &rows,
			       &min_columns,
			       &max_columns,
			       &max_length);
  
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n",
	       status, filename);
      return status;
    }
    
    if (rows<2) {
      fprintf (stderr, "Error: found less than two rows in %s\n", filename);
      return -1;
    }
    
    if (min_columns<3) {
      fprintf (stderr, "Error: found less than three columns in %s\n", filename);
      return -1;
    }
  
    /* reset string to the beginning of the memory block */
    string = line;
  
    /* open file */
    if ( (file = fopen(filename, "r")) == NULL)  
      return ASCIIFILE_NOT_FOUND;
  
  
    /* first line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);
    
      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }
  
    if (number<3) {
      fprintf (stderr, "Error: expected at least 3 values in 1st line, found only %d\n", number);
      return -1;
    }
  
    /* first 3 values are number of cells */
    out->Nxatm = strtol (array[0], &dummy, 0);
    out->Nyatm = strtol (array[1], NULL, 0);
    //out->Nzatm = strtol (array[2], NULL, 0) -1;

    out->nlev = strtol (array[2], NULL, 0);
    out->nlyr = out->nlev-1;
    out->nlev_common = strtol (array[2], NULL, 0);
    out->nlyr_common = out->nlev_common-1;
    /* all functions that are used for 1D and 3D should use Nzatm instead of nlev, nlyr, because those refer only to 1D profiles */
  
    //fprintf(stderr, "3DAbs read_atmosphere3D nlyr %d \n", out->nlyr); 
  
    /* unit of watervapor could be specified with 4th number, to be implemented later */
    /* if (number>3) */
    /*   *wv_unit = strtol (array[3], NULL, 0); */
  
    /* 2nd line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);
    
      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }
  
    if (number!= out->nlev+2) {
      fprintf (stderr, "Error: found %d z levels, expected %d\n",
	       number-2, out->nlev);
      return -1;
    }
  
    /* horizontal grid sizes, in km */
    out->dxatm = strtod (array[0], NULL);
    out->dyatm = strtod (array[1], NULL);
  
    if (!quiet) {
      fprintf (stderr, " ... reading %d x %d x %d data points from %s\n",
	       out->Nxatm, out->Nyatm, out->nlyr, filename);
      /* fprintf (stderr, " ... water vapor unit = %d\n", *wv_unit); */
      fprintf (stderr, " ... atmosphere grid size in x direction: %g km\n", out->dxatm);
      fprintf (stderr, " ... atmosphere grid size in y direction: %g km\n", out->dyatm);
    }

#if HAVE_LIBNETCDF
  }
#endif


  
  /* free and reallocate memory for altitude levels */
  free(out->zd);
  free(out->zd_common);
  out->zd = calloc(out->nlev, sizeof(float)); 
  out->zd_common = calloc(out->nlev, sizeof(float)); 
  
  /* 3DAbs for now we assume the same 3D grid for the full atmosphere, so we do not need */
  /* the variables z_3d_layer and threed */
  /* out->z_3d_layer = calloc ((size_t) out->Nlyr_3d, sizeof(float)); */
  /* out->threed = calloc ((size_t) out->Nlyr_3d, sizeof(int)); */
 
  /* free and reallocate memory for pressure, temperature, water vapor ... */
  ASCII_free_float_3D (out->microphys.press, 1, 1);
  ASCII_calloc_float_3D (&out->microphys.press, out->Nxatm, out->Nyatm, out->nlev);
  ASCII_free_float_3D (out->microphys.temper, 1, 1);
  ASCII_calloc_float_3D (&out->microphys.temper, out->Nxatm, out->Nyatm, out->nlev);
  ASCII_free_float_4D (out->microphys.dens, MOL_NN, 1, 1);
  ASCII_calloc_float_4D (&out->microphys.dens, MOL_NN, out->Nxatm, out->Nyatm, out->nlev);
  
  //3DAbs avg may be allocated later 
  ASCII_calloc_float_3D (&out->microphys.temper_avg, out->Nxatm, out->Nyatm, out->nlev-1);
  ASCII_calloc_float_4D (&out->microphys.dens_avg, MOL_NN, out->Nxatm, out->Nyatm, out->nlev-1);

#if HAVE_LIBNETCDF
  if (read_netcdf==1) {
    /* read from netcdf */
    if (!quiet)
      fprintf (stderr, " ... reading data from netcdf\n");
  }
  else {
#endif
    /* read from ASCII */
    if (!quiet)
      fprintf (stderr, " ... reading data from ASCII\n");
#if HAVE_LIBNETCDF
  }
#endif



#if HAVE_LIBNETCDF
  if (read_netcdf==1) {
    size_t zstart = 0;
    size_t nz_lev = out->nlev;
    RUN_NC(nc_inq_varid(ncid, "z", &id_z), "reading z");
    RUN_NC(nc_get_vara_float(ncid, id_z, &zstart, &nz_lev, out->zd), "reading z");

    /* reverse profile and copy to out->zd_common[] */
    for (i=0;i<out->nlev;i++) 
      out->zd_common[i] = out->zd[out->nlev-1-i];

    for (i=0;i<out->nlev;i++) 
      out->zd[i] = out->zd_common[i];

    if(!(tmpdata = calloc_float_3D(out->nlev, out->Nxatm, out->Nyatm, "tmpdata"))) return -1;
    
    size_t tstart[3] = {iystart, ixstart, 0};
    size_t tcount[3] = {out->Nyatm, out->Nxatm, out->nlev};
    ptrdiff_t tstride[3] = {1, 1, 1};
    ptrdiff_t timap[3] = {1, out->Nyatm, out->Nyatm * out->Nxatm};

    /* the reading procedure using the temporary array is a bit awkward but necessary because */
    /* the file format of the 3D water vapor file is the same as the 3D cloud file, but       */
    /* the target arrays are structured differently, wc[iz][ix][iy] vs. press[ix][iy][iz]     */

    /* read press */
    RUN_NC(nc_inq_varid(ncid, "press", &id_press), "reading press");
    RUN_NC(nc_get_varm_float(ncid, id_press, tstart, tcount, tstride, timap, **tmpdata), "reading press"); 
    /* ??? General question: "threed" is not set here. This doesn't cause any trouble because      */
    /* ??? later in optical_properties_atmosphere3D() it's assumed anyway that all layers are 3D   */
    /* ??? This should be changed eventually, because there is potential for saving memory         */
    /* ??? and computational time if 1D and 3D absorption cross sections can be combined, and      */
    /* ??? and tracing would be done differently in 1D and 3D, as it is done for all other caoth's */                                
    //check_threed(out->microphys.press, out->nlev, out->Nxatm * out->Nyatm, 0.0, threed);
    
    /* copy tmpdata to target array and invert z (top to bottom) */
    for (ix=0;ix<out->Nxatm; ix++)
      for (iy=0;iy<out->Nyatm; iy++)
	for (iz=0;iz<out->nlev; iz++)
	  out->microphys.press[ix][iy][iz] = tmpdata[out->nlev-1-iz][ix][iy];


    /* read temperature */
    RUN_NC(nc_inq_varid(ncid, "temp", &id_temp), "reading temperature");
    RUN_NC(nc_get_varm_float(ncid, id_temp, tstart, tcount, tstride, timap, **tmpdata), "reading temperature"); 
    //check_threed(out->microphys.temper, out->nlev, out->Nxatm * out->Nyatm, 0.0, threed);
    
    /* copy tmpdata to target array and invert z (top to bottom) */
    for (ix=0;ix<out->Nxatm; ix++)
      for (iy=0;iy<out->Nyatm; iy++)
	for (iz=0;iz<out->nlev; iz++)
	  out->microphys.temper[ix][iy][iz] = tmpdata[out->nlev-1-iz][ix][iy];

    
    /* read specific humidity */ 
    RUN_NC(nc_inq_varid(ncid, "qH2O", &id_qH2O), "reading specific humidity");
    RUN_NC(nc_get_varm_float(ncid, id_qH2O, tstart, tcount, tstride, timap, **tmpdata), "reading specific humidity"); 
    //check_threed(out->microphys.dens[MOL_H2O], out->nlev, out->Nxatm * out->Nyatm, 0.0, threed);
    
    /* copy tmpdata to target array and invert z (top to bottom) */
    for (ix=0;ix<out->Nxatm; ix++)
      for (iy=0;iy<out->Nyatm; iy++)
	for (iz=0;iz<out->nlev; iz++)
	  out->microphys.dens[MOL_H2O][ix][iy][iz] = tmpdata[out->nlev-1-iz][ix][iy];

    /* free memory of temporary data */
    free_float_3D(tmpdata);

    /* close netcdf file */
    nc_close (ncid);
  }
  else {
#endif
    /* rest of line 2 */
    /* set altitude levels */
    for (i=2; i<out->nlev +2; i++){
      /*altitudes in out->zd are sorted in descending order,*/
      /* 3D input files in ascending order */
      out->zd[i-2]        = (float) strtod(array[out->nlev +3-i], NULL); 
      out->zd_common[i-2] = (float) strtod(array[out->nlev +3-i], NULL);
    }
    free(array);


    for (i=2; i<rows; i++) {
      number=0;
      while (number==0) {
	dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
	status = ASCII_parsestring (string, &array, &number);
	
	if (status!=0) {
	  fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	  return status;
	}
      }
      
      /* cell indices */
      ix = strtol(array[0], NULL, 0) - 1;
      iy = strtol(array[1], NULL, 0) - 1;
      iz = strtol(array[2], NULL, 0) - 1;
      
      if (ix<0 || ix>=out->Nxatm) {
	fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix+1, i+1);
	return -1;
      }
      
      if (iy<0 || iy>=out->Nyatm) {
	fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy+1, i+1);
	return -1;
      }
      
      if (iz<0 || iz>=out->nlev) {
	fprintf (stderr, "Error, iz = %d out of bounds in line %d\n", iz+1, i+1);
	return -1;
      }
      
      out->microphys.press [ix][iy][out->nlyr -iz] = strtod (array[3], NULL);
      out->microphys.temper[ix][iy][out->nlyr -iz]  = strtod (array[4], NULL);
      out->microphys.dens[MOL_H2O] [ix][iy][out->nlyr -iz]  = strtod (array[5], NULL);
      
      /*other trace gases can follow in further columns, to be implemented. */
      
      free(array);
    }

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", 
	       rows-2, filename);

    /* close input file */
    fclose (file);
    
#if HAVE_LIBNETCDF
  }
#endif


  // print 3D p,T,qH2O profile
  //  for (ix=0; ix<out->Nxatm; ix++)
  //    for (iy=0; iy<out->Nyatm; iy++)
  //      for (iz=0;iz<out->nlev;iz++)
  //	fprintf (stderr, "%d %d %d %f %f %f %f\n", ix, iy, iz, out->zd[iz], out->microphys.press[ix][iy][iz], out->microphys.temper[ix][iy][iz], out->microphys.dens[MOL_H2O][ix][iy][iz]);
    

  /* calculate air number density from pressure and temperature, as well as layer average quantities */
  for (ix=0; ix<out->Nxatm; ix++){
    for (iy=0; iy<out->Nyatm; iy++){
      for (iz=0; iz<out->nlev; iz++){
        out->microphys.dens[MOL_AIR][ix][iy][iz]=out->microphys.press[ix][iy][iz] * 1E-04 / (BOLTZMANN*out->microphys.temper[ix][iy][iz]);
        
        /* convert from kg/kg (specific humidity) to concentration 1/cm^3  */
        out->microphys.dens[MOL_H2O] [ix][iy][iz] = mol_mass[MOL_AIR]/mol_mass[MOL_H2O]*out->microphys.dens[MOL_H2O][ix][iy][iz]*out->microphys.dens[MOL_AIR][ix][iy][iz];
      }
      
      /* layer average temperature */
      status = average_dens (out->microphys.temper[ix][iy],
			     out->microphys.dens[MOL_AIR][ix][iy], 
			     out->zd_common, out->nlyr+1,
			     interpol_method_temper, 
			     &(out->microphys.temper_avg[ix][iy]), NO);

      if (status!=0) {
	fprintf (stderr, "Error %d calculating average temperature (line %d, function '%s' in '%s')\n", 
		 status, __LINE__, __func__, __FILE__ );
	return status;
      }
      
      /* layer average air density  */
      status = average_dens (out->microphys.dens[MOL_AIR][ix][iy],
			     out->microphys.dens[MOL_AIR][ix][iy],
			     out->zd_common, out->nlyr+1,
			     interpol_method_gas[MOL_H2O],
			     &(out->microphys.dens_avg[MOL_AIR][ix][iy]), NO); 

      if (status!=0) {
	fprintf (stderr, "Error %d calculating average air density (line %d, function '%s' in '%s')\n", 
		 status, __LINE__, __func__, __FILE__ );
	return status;
      }

      /* layer average water vapor */
      status = average_dens (out->microphys.dens[MOL_H2O][ix][iy],
			     out->microphys.dens[MOL_AIR][ix][iy],
			     out->zd_common, out->nlyr+1,
			     interpol_method_gas[MOL_H2O],
			     &(out->microphys.dens_avg[MOL_H2O][ix][iy]), NO); 

      if (status!=0) {
	fprintf (stderr, "Error %d calculating average water vapor density (line %d, function '%s' in '%s')\n", 
		 status, __LINE__, __func__, __FILE__ );
	return status;
      }
    }
  }
  
  

  return 0; 
}
  

/*********************************************************/
/* Function: write_dummy_profile_file                    */
/* Description:                                          */
/* Generate dummy profile_file from 3d molecular file    */
/*                                                       */ 
/* Author:                                                */
/* August 2017    C. Emde   Created                      */
/*********************************************************/
static int write_dummy_profile_file(char *filename, 
				    int quiet, int verbose){

  int i=0, number=0, status=0, ix=0, iy=0, iz=0, lc=0;
  
  int rows=0, min_columns=0, max_columns=0, max_length=0;
  
  char line[MAX_LENGTH_OF_LINE+1]="";
  char *string=NULL;

  FILE *file=NULL;;
  char **array=NULL;
  char *dummy=NULL;
  float *zd=NULL, *tmpzd=NULL; 
  FILE *file_dummy=NULL;
  int Nxatm=0, Nyatm=0, nlev=0;
  double dxatm=0.0, dyatm=0.0; 
  int *threed=NULL;

#if HAVE_LIBNETCDF
  int ixstart=0, iystart=0;
  size_t N;

  size_t n=0;
  int read_netcdf=0;

  int    ncid   =0, idd_nx =0, idd_ny =0, idd_nz=0;
  int    id_z=0, id_press=0, id_temp=0, id_qH2O=0;
  
  float ***tmpdata=NULL;
#endif

  if (!quiet) 
    fprintf (stderr, " ... generate dummy profile file from 3d atmosphere file %s\n", filename);


  #if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {
    read_netcdf=1;

    if (!quiet)
      fprintf (stderr, " ... reading Cloud data from netCDF file %s\n", filename);

    RUN_NC(nc_inq_dimid(ncid, "nx", &idd_nx), "reading nx");
    RUN_NC(nc_inq_dimlen(ncid, idd_nx, &n), "reading nx");
    Nxatm=n;

    RUN_NC(nc_inq_dimid(ncid, "ny", &idd_ny), "reading ny");
    RUN_NC(nc_inq_dimlen(ncid, idd_ny, &n), "reading ny");
    Nyatm=n;

    RUN_NC(nc_inq_dimid(ncid, "nz_lev", &idd_nz), "reading nz");
    RUN_NC(nc_inq_dimlen(ncid, idd_nz, &n), "reading nz");

    nlev = n;
    
    RUN_NC(nc_get_att_double(ncid, NC_GLOBAL, "dx", &dxatm), "reading dx");  
    RUN_NC(nc_get_att_double(ncid, NC_GLOBAL, "dy", &dyatm), "reading dy");

    /* read z levels */
    tmpzd = calloc(nlev, sizeof(float)); 
    zd = calloc(nlev, sizeof(float)); 

    size_t zstart = 0;
    size_t nz_lev = nlev;
    RUN_NC(nc_inq_varid(ncid, "z", &id_z), "reading z");
    RUN_NC(nc_get_vara_float(ncid, id_z, &zstart, &nz_lev, tmpzd), "reading z");

    /* reverse profile */
    for (i=0;i<nlev;i++) 
      zd[i] = tmpzd[nlev-1-i];

    free (tmpzd);
  }
  else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", 
	       filename);
#endif

    /* check input file */
    status =  ASCII_checkfile (filename, 
			       &rows,
			       &min_columns,
			       &max_columns,
			       &max_length);
  
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n",
	       status, filename);
      return status;
    }
  
    if (rows<2) {
      fprintf (stderr, "Error: found less than two rows in %s\n", filename);
      return -1;
    }
  
    if (min_columns<3) {
      fprintf (stderr, "Error: found less than three columns in %s\n", filename);
      return -1;
    }
  
    /* reset string to the beginning of the memory block */
    string = line;
  
    /* open file */
    if ( (file = fopen(filename, "r")) == NULL)  
      return ASCIIFILE_NOT_FOUND;
  
  
    /* first line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);
    
      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }
  
    if (number<3) {
      fprintf (stderr, "Error: expected at least 3 values in 1st line, found only %d\n", number);
      return -1;
    }
  
    /* first 3 values are number of cells */
    Nxatm = strtol (array[0], &dummy, 0);
    Nyatm = strtol (array[1], NULL, 0);
    nlev = strtol (array[2], NULL, 0);

    if (nlev<2) {
      fprintf (stderr, "Error, less than 2 levels in %s\n", filename);
      return -1;
    }

  
    /* 2nd line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);
    
      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }
  
    if (number!= nlev+2) {
      fprintf (stderr, "Error: found %d z levels, expected %d\n",
	       number-2, nlev);
      return -1;
    }
  
    /* horizontal grid sizes, in km */
    dxatm = strtod (array[0], NULL);
    dyatm = strtod (array[1], NULL);
    zd = calloc(nlev, sizeof(float)); 

    /* rest of line 2 */
    /* set altitude levels */
    for (i=2; i<nlev+2; i++){
      /* altitudes in out->zd are sorted in descending order,*/
      /* 3D input files in ascending order */
      zd[i-2] = (float) strtod(array[nlev+3-i], NULL); 
    }
    free(array);
    
#if HAVE_LIBNETCDF
  }
#endif
  


  
  /* BM, read rest of file - we need to know which layers have 3D absorption */
  threed = calloc (nlev-1, sizeof(int));

#if HAVE_LIBNETCDF
  if (read_netcdf==1) {

    // ??? commented out 12.9.2018, BM. At present, all layers are assumed 3D anyway, hence we simply set
    // ??? the whole vector threed[] to 1. Need to reactivate later!
    // ???
    // ??? BEGIN checkthreed
    //    if(!(tmpdata = calloc_float_3D(nlev, Nxatm, Nyatm, "tmpdata"))) return -1;
    //    
    //    size_t tstart[3] = {iystart, ixstart, 0};
    //    size_t tcount[3] = {Nyatm, Nxatm, nlev};
    //    ptrdiff_t tstride[3] = {1, 1, 1};
    //    ptrdiff_t timap[3] = {1, Nyatm, Nyatm * Nxatm};
    //
    //    /* the reading procedure using the temporary array is a bit awkward but necessary because */
    //    /* the file format of the 3D water vapor file is the same as the 3D cloud file, but       */
    //    /* the target arrays are structured differently, wc[iz][ix][iy] vs. press[ix][iy][iz]     */
    //
    //    /* read press */
    //    RUN_NC(nc_inq_varid(ncid, "press", &id_press), "reading press");
    //    RUN_NC(nc_get_varm_float(ncid, id_press, tstart, tcount, tstride, timap, **tmpdata), "reading press"); 
    //
    //    /* if any grid box in a level has a value > 0, set threed (in inverse order!) to 1 */
    //    for (ix=0;ix<Nxatm; ix++)
    //      for (iy=0;iy<Nyatm; iy++)
    //	for (iz=0;iz<nlev-1; iz++)
    //	  if (tmpdata[nlev-1-iz][ix][iy]>=0)
    //	    threed[iz] = 1;
    //
    //
    //    
    //    /* read temperature */
    //    RUN_NC(nc_inq_varid(ncid, "temp", &id_temp), "reading temperature");
    //    RUN_NC(nc_get_varm_float(ncid, id_temp, tstart, tcount, tstride, timap, **tmpdata), "reading temperature"); 
    //    
    //    /* if any grid box in a level has a value > 0, set threed (in inverse order!) to 1 */
    //    for (ix=0;ix<Nxatm; ix++)
    //      for (iy=0;iy<Nyatm; iy++)
    //	for (iz=0;iz<nlev-1; iz++)
    //	  if (tmpdata[nlev-1-iz][ix][iy]>=0)
    //	    threed[iz] = 1;
    //
    //    
    //    /* read specific humidity */ 
    //    RUN_NC(nc_inq_varid(ncid, "qH2O", &id_qH2O), "reading specific humidity");
    //    RUN_NC(nc_get_varm_float(ncid, id_qH2O, tstart, tcount, tstride, timap, **tmpdata), "reading specific humidity"); 
    //    
    //    /* if any grid box in a level has a value > 0, set threed (in inverse order!) to 1 */
    //    for (ix=0;ix<Nxatm; ix++)
    //      for (iy=0;iy<Nyatm; iy++)
    //	for (iz=0;iz<nlev-1; iz++)   
    //	  if (tmpdata[nlev-1-iz][ix][iy]>=0)
    //	    threed[iz] = 1;
    //
    //    /* BM: This is a bit tricky! Since water vapor is specified at levels, we only want a layer to be 3D */
    //    /* if lower and upper level are 3D - to be fixed when 3D and 1D can actually be combined ???         */
    //
    //    /* free memory of temporary data */
    //    free_float_3D(tmpdata);
    // ??? END checkthreed

    fprintf (stderr, "Attention! It is not checked if there are any non-3D layers in the 3D atmosphere file\n");
    fprintf (stderr, "since this is not possible at the moment anyway; however, once this is possible,\n");
    fprintf (stderr, "the check should be activated again in molecular3d.c\n");
    
    /* ??? instead, assume all layers are 3D */
    for (iz=0;iz<nlev-1; iz++)
      threed[iz] = 1;
    
    /* close netcdf file */
    nc_close (ncid);
  }
  else {
#endif
  
    for (i=2; i<rows; i++) {
      number=0;
      while (number==0) {
	dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
	status = ASCII_parsestring (string, &array, &number);
      
	if (status!=0) {
	  fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	  return status;
	}
      }
      
      if (number!=6) {
	fprintf (stderr, "Error: found %d columns in row %d, expected 6\n", number, i);
	return -1;
      }
      
      /* BM: This is a bit tricky! Since water vapor is specified at levels, we only want a layer to be 3D */
      /* if lower and upper level are 3D - to be fixed when 3D and 1D can actually be combined ???         */
      iz = strtol(array[2], NULL, 0) - 1;
      if (iz<nlev-1)     /* BM: prevent adressing "layer" nlev-1 */
	threed[iz] = 1;
      
      free(array);
    }

    /* close input file */
    fclose (file);
    
#if HAVE_LIBNETCDF
  }
#endif
  
    
  for (lc=0; lc<nlev-1; lc++)
    if (threed[lc] && !quiet)
      fprintf (stderr, "FOUND absorber in layer %d\n", lc);
  
  /* write dummy profile file */
  if ( (file_dummy = fopen("profile_mol3d_dummy.dat", "w")) == NULL)  
    return ASCIIFILE_NOT_FOUND;

  fprintf(file_dummy, "%d %d %d 3\n", Nxatm, Nyatm, nlev-1);
  fprintf(file_dummy, "%f %f ",dxatm, dyatm);
  for (i=0; i<nlev; i++)
    fprintf(file_dummy, "%f ",  zd[nlev-i-1]);
  fprintf(file_dummy, "\n");

  /* BM, write all layers with 3D absorption */
  for (lc=0; lc<nlev-1; lc++)
    if (threed[lc])
      fprintf(file_dummy, "1 1 %d 0.0 10 \n", lc+1);

  fclose(file_dummy);
  
  return 0;
}

  
/****************************************************************/
/* Regrid a 3D float array [ix, iy, iz] from a given input grid */
/* (z_old) to a given output grid.                              */
/****************************************************************/

static int redistribute_3D_temp (float ****data, int nx, int ny,
                                 float *z_old, int nlev_old, 
                                 float *z_new, int nlev_new, double *temper_1D) 
{
  int lc=0, status=0;
  int ix=0, iy=0;
  int nlev_tmp=0, start_lc=0; 
  float ***tmp=NULL, *z_tmp=NULL, *data_tmp=NULL; 

  ASCII_calloc_float_3D (&tmp, nx, ny, nlev_old); 
  /* copy original stuff to temporary array */

  /* common grid range, assume that 3D grid is defined inside lower part of 1D standard atmosphere */
  lc=0;
  while (z_new[lc]>z_old[0])
    lc++;

  start_lc=lc;
  nlev_tmp=nlev_new-start_lc;
  z_tmp=calloc(nlev_tmp, sizeof(float));

  for (lc=0; lc<nlev_new-start_lc; lc++)
    z_tmp[lc]=z_new[lc+start_lc];
  
  for (ix=0; ix<nx; ix++)
    for (iy=0; iy<ny; iy++)
      for (lc=0; lc<nlev_old; lc++){
        /* fprintf(stderr, "before ix %d iy %d lc %d  T %g \n", ix, iy, lc,  (*data)[ix][iy][lc]); */
        tmp[ix][iy][lc] = (*data)[ix][iy][lc];
      }
  
 
  data_tmp = calloc(nlev_tmp, sizeof(float));
  ASCII_free_float_3D (*data, nx, ny);
  ASCII_calloc_float_3D (data, nx, ny, nlev_new);
  
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      
      status = arb_wvn (nlev_old, z_old, tmp[ix][iy],
                        nlev_tmp, z_tmp, data_tmp,
                        1,1);
      
      if (status != 0) {
        fprintf(stderr, "Error %d returned by arb_wvn() in redistribute_temp (molecular_3d.c)\n", status);
        return status;
      }
      
      for (lc=0; lc<start_lc; lc++){
        (*data)[ix][iy][lc] = temper_1D[lc];
        /* fprintf(stderr, "after ix %d iy %d lc %d T %g \n", ix, iy, lc, (*data)[ix][iy][lc]); */
      }
      for (lc=start_lc; lc<nlev_new; lc++){
        (*data)[ix][iy][lc] = data_tmp[lc-start_lc];
        /* fprintf(stderr, "after ix %d iy %d lc %d T %g \n", ix, iy, lc, (*data)[ix][iy][lc]); */
      }
    }
  }
  
  
  ASCII_free_float_3D (tmp, nx, ny);
  free(data_tmp);
  free(z_tmp);
  
  return status;
}

