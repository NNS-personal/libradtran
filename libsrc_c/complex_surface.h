/************************************************************************
 * $Id: complex_surface.h 3503 2019-10-04 13:07:25Z Claudia.Emde $
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
typedef struct t_star_engine {
#ifdef HAVE_STAR_ENGINE
    struct s3d_device     *dev;
    struct s3d_scene_view *view;
    struct s3d_scene      *scene;
    //struct s3d_hit *hit;/*TODO*/
#endif
} t_star_engine;

typedef struct {
    size_t N_vertices;
    size_t N_triangles;
    size_t N_coords;
    size_t N_vert_ind;
    size_t N_surface_prop;
    double **vertices;
    int **triangles;
    double *surface_props_list;
    int *surface_props_ind;
    struct t_star_engine *star_engine_tracer;
} triangular_surface_struct;

int setup_triangular_surface( char                   *filename,
                              triangular_surface_struct *triangular_surface,
                              int                    quiet);

int cross_triangular_surface(
    const float origin[3],
    const float dir[3],
    triangular_surface_struct* srfc, 
    double *distance,
    int *primitive_id);

int get_triangle_albedo(
    int primitive_id,
    double *triangle_albedo,
    triangular_surface_struct* srfc);

int get_triangle_normal(
    int primitive_id,
    double *triangle_normal,
    triangular_surface_struct* srfc);

void get_triangle_zmax(
    double *triangle_zmax,
    triangular_surface_struct *srfc);

