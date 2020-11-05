/************************************************************************
 * $Id: complex_surface.c 3524 2019-12-23 01:49:09Z bernhard.mayer $
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 * Authors of complex srfc: Marc Schwaerzel, Claudia Emde and Fabian Jakub
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if HAVE_LIBNETCDF
#include <netcdf.h>

#define RUN_NC(command, msg) do {               \
    int status = (command);                     \
    if (status != NC_NOERR) {                                           \
      fprintf(stderr, "netCDF error %d: \"%s\" %s from %s\n", status, nc_strerror(status), (msg), filename); \
      return status;                                                    \
    }                                                                   \
  } while(0)

#endif

#include "complex_surface.h"
#include "allocnd.h"


#ifdef HAVE_STAR_ENGINE
#include "star/s3d.h"
#include "star/s3daw.h"
#include "rsys/mem_allocator.h"
#include "math.h"
#endif

#ifdef HAVE_STAR_ENGINE
/* Initialize Star-Engine objects used for raytracing */
void surf_get_indices(const unsigned primitive_id, unsigned indices[3], void* data){  

  /* The following array lists the indices toward the 3D vertices of each
   * triangle.

   *       ,2----,3
   *     ,'0 \1,'
   *    0----1'
   * e.g. id=1 -> indices={1,3,2}
   */
  triangular_surface_struct* srfc = (triangular_surface_struct*)data;
  indices[0] = srfc->triangles[primitive_id][0];
  indices[1] = srfc->triangles[primitive_id][1];
  indices[2] = srfc->triangles[primitive_id][2];
}


void surf_get_position(const unsigned vertex_id, float position[3], void* data)  
{
  triangular_surface_struct* srfc = (triangular_surface_struct*)data;
  /* Memory layout of the vertex positions
   *   2.....3
   *  /     /
   * 0-----1
  */
  position[0] = srfc->vertices[vertex_id][0];
  position[1] = srfc->vertices[vertex_id][1];
  position[2] = srfc->vertices[vertex_id][2];
}

int s4vs_discard_self_hit
  (const struct s3d_hit* hit,
   const float ray_org[3],
   const float ray_dir[3],
   void* ray_data,
   void* filter_data)
{
  const struct s3d_primitive* prim_from = ray_data;

  /* Avoid unused variable warn */
  (void)ray_org, (void)ray_dir, (void)filter_data;
  return prim_from ? S3D_PRIMITIVE_EQ(prim_from, &hit->prim) : 0;
}

void set_attach_surf(struct s3d_device * s3d, struct s3d_scene * scene, triangular_surface_struct* srfc){  

  struct s3d_shape* surf;
  struct s3d_vertex_data vertex_attribs[2];

  vertex_attribs[0].usage = S3D_POSITION;
  vertex_attribs[0].type = S3D_FLOAT3;
  vertex_attribs[0].get = surf_get_position;
  vertex_attribs[1] = S3D_VERTEX_DATA_NULL;

  s3d_shape_create_mesh(s3d, &surf);
  s3d_mesh_setup_indexed_vertices(surf, srfc->N_triangles/*#triangles*/, surf_get_indices, srfc->N_vertices/*#vertices*/, vertex_attribs, 1, (void*)srfc);
  s3d_mesh_set_hit_filter_function(surf, s4vs_discard_self_hit,(void*)srfc);/* NULL);*/
  s3d_scene_attach_shape(scene, surf);
  s3d_shape_ref_put(surf) ;
}

struct t_star_engine* init_star_engine(
        triangular_surface_struct *triangular_surface,
        int quiet)
{
    if (!quiet) fprintf(stderr, "Initializing Star-Engine-Tracer!\n");
    struct t_star_engine* star_engine = malloc(sizeof( struct t_star_engine ));

    struct s3d_device *dev;
    struct s3d_scene_view *view;
    struct s3d_scene  *scene;
    s3d_device_create(NULL, NULL, 0, &dev);                      /*initialize device*/
    s3d_scene_create(dev, &scene);                  /*initialize scene*/
    set_attach_surf(dev, scene, triangular_surface);
    s3d_scene_view_create(scene, S3D_TRACE, &view); /*initialize geometry*/
    if (!quiet) fprintf(stderr, "Initializing Star-Engine-Tracer ... done!\n");

    star_engine->dev  = dev;
    star_engine->view = view;
    star_engine->scene= scene;
    S3D(device_ref_put(dev));
    S3D(scene_ref_put(scene));
    printf("Initializing Star-Engine-Tracer done sucessfully!\n");
    return star_engine;
}
#endif


/*******************************************************************/
/* Read surface data and setup triangular surface.                 */
/*******************************************************************/

int setup_triangular_surface( char                      *filename,
                              triangular_surface_struct *triangular_surface,
                              int                       quiet)
{


#if HAVE_LIBNETCDF

  int status=0;

  int ncid=0;
  int idd_Nvert=0, idd_Ncoord=0, idd_Ntriangl=0, idd_Nvert_ix, idd_Nsurface_props=0;
  int id_vertices=0, id_triangles=0, id_surface_props_list=0,
    id_surface_props_ind=0;

  size_t n=0;


  /* open as netCDF file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading triangular surface data from netCDF file %s\n", filename);

    RUN_NC(nc_inq_dimid(ncid, "Nvert", &idd_Nvert), "reading Nvert");
    RUN_NC(nc_inq_dimlen(ncid, idd_Nvert, &n), "reading Nvert");
    triangular_surface->N_vertices = n;

    RUN_NC(nc_inq_dimid(ncid, "Ncoord", &idd_Ncoord), "reading Ncoord");
    RUN_NC(nc_inq_dimlen(ncid, idd_Ncoord, &n), "reading Ncoord");
    triangular_surface->N_coords = n;

    RUN_NC(nc_inq_dimid(ncid, "Ntriangl", &idd_Ntriangl), "reading Ntriangl");
    RUN_NC(nc_inq_dimlen(ncid, idd_Ntriangl, &n), "reading Ntriangl");
    triangular_surface->N_triangles = n;

    RUN_NC(nc_inq_dimid(ncid, "Nvert_ix", &idd_Nvert_ix), "reading Nvert_ix");
    RUN_NC(nc_inq_dimlen(ncid, idd_Nvert_ix, &n), "reading Nvert_ix");
    triangular_surface->N_vert_ind = n;

    RUN_NC(nc_inq_dimid(ncid, "Nsurface_props", &idd_Nsurface_props), "reading Nsurface_props");
    RUN_NC(nc_inq_dimlen(ncid, idd_Nsurface_props, &n), "reading Nsurface_props");
    triangular_surface->N_surface_prop = n;

    /* allocate memory for vertices */
    if(!(triangular_surface->vertices = calloc_double_2D(triangular_surface->N_vertices, triangular_surface->N_coords, "vertices" ))) return -1;

    RUN_NC(nc_inq_varid(ncid, "vertices", &id_vertices), "reading vertices");
    RUN_NC(nc_get_var_double(ncid, id_vertices, *triangular_surface->vertices), "reading vertices");

    if(!(triangular_surface->triangles = calloc_int_2D(triangular_surface->N_triangles, triangular_surface->N_vert_ind, "triangles" ))) return -1;

    RUN_NC(nc_inq_varid(ncid, "triangles", &id_triangles), "reading triangles");
    RUN_NC(nc_get_var_int(ncid, id_triangles, *triangular_surface->triangles), "reading triangles");

    triangular_surface->surface_props_list = calloc ( triangular_surface->N_surface_prop, sizeof(double));

    RUN_NC(nc_inq_varid(ncid, "surface_properties_list", &id_surface_props_list), "reading surface_properties_list");
    RUN_NC(nc_get_var_double(ncid, id_surface_props_list, triangular_surface->surface_props_list), "reading surface_properties_list");

    triangular_surface->surface_props_ind = calloc( triangular_surface->N_triangles, sizeof(int));

    RUN_NC(nc_inq_varid(ncid, "surface_properties_index", &id_surface_props_ind), "reading surface_properties_index");
    RUN_NC(nc_get_var_int(ncid, id_surface_props_ind, triangular_surface->surface_props_ind), "reading surface_properties_index");

    nc_close (ncid);
  }
  else {
    fprintf (stderr, "Error: unable to read triangular surface properties file %s. It is expected to be in netcdf format. \n", filename);
    return -1;
  }

#endif

  triangular_surface->star_engine_tracer = NULL;
#ifdef HAVE_STAR_ENGINE
  triangular_surface->star_engine_tracer = init_star_engine(triangular_surface, quiet);
#endif

  return status;
}

int cross_triangular_surface(
    const float origin[3],
    const float dir[3],
    triangular_surface_struct* srfc,  
    double *distance,
    int *primitive_id)
{

#ifdef HAVE_STAR_ENGINE
    if(srfc->star_engine_tracer) {
        struct s3d_hit hit = S3D_HIT_NULL;
        float rng[] = {0,HUGE_VAL};      /*big range for how far we look*/
        
        s3d_scene_view_trace_ray(srfc->star_engine_tracer->view, origin, dir, rng, NULL, &hit);
        if(S3D_HIT_NONE(&hit)) {
            *distance = -1;
        } else {
            *distance = hit.distance;
	        *primitive_id = hit.prim.prim_id;
        }

        return 0;
    }
#endif
      *distance=-1;
    return 1; // did not use any intersection method?
}

int get_triangle_albedo(
    int primitive_id,
    double *triangle_albedo,
    triangular_surface_struct* srfc){  
 
    *triangle_albedo = srfc->surface_props_list[srfc->surface_props_ind[primitive_id]];
    return 1;
}

void crossProduct(double a[3], double b[3], double c[3])
{
	double d[] = {a[2]*b[2] - a[2]*b[1], a[2]*b[0] - 
			a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
	for(int i=0; i<3; i++){
		c[i] = d[i];
	}
}

void normalize(double u[3], double *U){
    double norm = sqrt (u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );   
    for (int i=0;i<3;i++){U[i]=u[i]/norm;}
}

int get_triangle_normal(
    int primitive_id,
    double *triangle_normal,
    triangular_surface_struct* srfc){
    double u[3] = {
         srfc->vertices[srfc->triangles[primitive_id][1]][0]-srfc->vertices[srfc->triangles[primitive_id][0]][0],
         srfc->vertices[srfc->triangles[primitive_id][1]][1]-srfc->vertices[srfc->triangles[primitive_id][0]][1],
         srfc->vertices[srfc->triangles[primitive_id][1]][2]-srfc->vertices[srfc->triangles[primitive_id][0]][2]
    };
    /* normalize u */
    double U[3];
    normalize(u, &U);
    double v[3] = {
         srfc->vertices[srfc->triangles[primitive_id][2]][0]-srfc->vertices[srfc->triangles[primitive_id][0]][0],
         srfc->vertices[srfc->triangles[primitive_id][2]][1]-srfc->vertices[srfc->triangles[primitive_id][0]][1],
         srfc->vertices[srfc->triangles[primitive_id][2]][2]-srfc->vertices[srfc->triangles[primitive_id][0]][2]
    };
    /* normalize V */
    double V[3];
    normalize(v, &V);
    double normal_vector[3];
    crossProduct(U, V, normal_vector); 
    for(int i=0;i<3;i++){triangle_normal[i] = normal_vector[i];}
    return 1;
}

void get_triangle_zmax(
    double *triangle_zmax,
    triangular_surface_struct *srfc){
    double zmax = 0.0;
    for(int i=0;i<srfc->N_vertices;i++){
      if(srfc->vertices[i][2]>zmax){
        zmax = srfc->vertices[i][2];
      }
    }
    *triangle_zmax = zmax;
}

