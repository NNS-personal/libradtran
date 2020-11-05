#include "allocnd.h"

#include <stdio.h>

#define CALLOC_1D(type, name) \
type * calloc_##name##_1D(size_t sx, const char* varname) { \
    type * data = calloc(sx, sizeof(type)); \
    if(!data) { \
        fprintf (stderr, "Error allocating memory for %s\n", varname); \
    } \
    return data; \
} \
 \
void free_##name##_1D(type * data) { \
    if(data) { \
        free(data); \
    } \
}

#define CALLOC_2D(type, name) \
type ** calloc_##name##_2D(size_t sx, size_t sy, const char* varname) { \
    type * data = calloc_##name##_1D(sx * sy, varname); \
    if(!data) return NULL; \
    type ** p1 = malloc(sx * sizeof(type*)); \
    if(!p1) { \
        free_##name##_1D(data); \
        fprintf (stderr, "Error allocating memory for 1st pointer level of %s\n", varname); \
        return NULL; \
    } \
   for(size_t i = 0; i < sx; ++i) {		\
        p1[i] = data + i * sy; \
    } \
    return p1;  \
} \
 \
void free_##name##_2D(type ** data) { \
    if(data) { \
        free_##name##_1D(*data); \
        free(data); \
    } \
}

#define CALLOC_3D(type, name) \
type *** calloc_##name##_3D(size_t sx, size_t sy, size_t sz, const char* varname) { \
    type ** p1 = calloc_##name##_2D(sx * sy, sz, varname); \
    if(!p1) return NULL; \
    type *** p2 = malloc(sx * sizeof(type**)); \
    if(!p2) { \
        free_##name##_2D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < sx; ++i) { \
        p2[i] = p1 + i * sy; \
    } \
    return p2; \
} \
 \
void free_##name##_3D(type *** data) { \
    if(data) { \
        free_##name##_2D(*data); \
        free(data); \
    } \
}

#define CALLOC_4D(type, name) \
  type **** calloc_##name##_4D(size_t sx, size_t sy, size_t sz, size_t dim4, const char* varname) { \
    type *** p1 = calloc_##name##_3D(sx * sy, sz, dim4, varname);	\
    if(!p1) return NULL; \
    type **** p2 = malloc(sx * sizeof(type***)); \
    if(!p2) { \
        free_##name##_3D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < sx; ++i) { \
        p2[i] = p1 + i * sy; \
    } \
    return p2; \
} \
 \
void free_##name##_4D(type **** data) { \
    if(data) { \
      free_##name##_3D(*data);			\
      free(data); \
    } \
}

#define CALLOC_5D(type, name) \
  type ***** calloc_##name##_5D(size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, const char* varname) { \
    type **** p1 = calloc_##name##_4D(sx * sy, sz, dim4, dim5, varname);	\
    if(!p1) return NULL; \
    type ***** p2 = malloc(sx * sizeof(type****)); \
    if(!p2) { \
        free_##name##_4D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < sx; ++i) { \
        p2[i] = p1 + i * sy; \
    } \
    return p2; \
} \
 \
void free_##name##_5D(type ***** data) { \
    if(data) { \
        free_##name##_4D(*data); \
        free(data); \
    } \
}




#define CALLOC_1D_RESTRICTED(type, name) \
type * calloc_##name##_1D_restricted(size_t s1, const char* varname) { \
    type * data = calloc(s1, sizeof(type)); \
    if(!data) { \
        fprintf (stderr, "Error allocating memory for %s\n", varname); \
    } \
  return data;					\
} 

#define CALLOC_2D_RESTRICTED(type, name) \
  type ** calloc_##name##_2D_restricted(size_t s1, size_t s2_lower, size_t s2_upper, int *s3, const char* varname) { \
  int sum=s3[s2_lower];	      					\
  size_t *cum=calloc(s2_upper+1, sizeof(size_t));				\
  for(size_t i = s2_lower+1; i <= s2_upper; ++i) {					\
      cum[i] = cum[i-1]+s3[i-1];\
      sum+=s3[i];					\
    } \
  									\
    type * data = calloc_##name##_1D_restricted(s1*sum, varname); \
    if(!data) return NULL; \
    type ** p1 = malloc(s1 * (s2_upper - s2_lower + 1) * sizeof(type*)); \
    if(!p1) { \
        free_##name##_1D(data); \
        fprintf (stderr, "Error allocating memory for 1st pointer level of %s\n", varname); \
        return NULL; \
    } \
    size_t counter=0;				\
    for(size_t i = 0; i < s1; ++i) {				\
      for(size_t j = s2_lower; j <= s2_upper; ++j) {		\
        p1[counter++] = data + sum*i + cum[j];	\
      }		\
    }		\
    return p1;  \
} 

#define CALLOC_3D_RESTRICTED(type, name) \
  type *** calloc_##name##_3D_restricted(size_t s1, size_t s2_lower, size_t s2_upper, int *s3, const char* varname) { \
    type ** p1 = calloc_##name##_2D_restricted(s1, s2_lower, s2_upper, s3, varname);	\
    if(!p1) return NULL; \
    type *** p2 = malloc(s1 * sizeof(type**)); \
    if(!p2) { \
        free_##name##_2D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < s1; ++i) { \
     /* subtract s2_lower because wavelength index is addressed from lower ... upper rather than 0 ... upper - lower */ \
     p2[i] = p1 + i * (s2_upper - s2_lower + 1) - s2_lower;	\
    }	       \
    return p2;					\
} \
void free_##name##_3D_restricted(type *** data, size_t s4_lower) {			\
    if(data) { \
      /* undo "subtract s2_lower ..." from CALLOC_3D_RESTRICTED(type, name) */ \
      free_##name##_2D(*data+s4_lower);		\
      free(data); \
    } \
}

#define CALLOC_4D_RESTRICTED(type, name) \
  type **** calloc_##name##_4D_restricted(size_t s1, size_t s2, size_t s3_lower, size_t s3_upper, int *s4, const char* varname) { \
    type *** p1 = calloc_##name##_3D_restricted(s1 * s2, s3_lower, s3_upper, s4, varname); \
    if(!p1) return NULL; \
    type **** p2 = malloc(s1 * sizeof(type***)); \
    if(!p2) { \
        free_##name##_3D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < s1; ++i) { \
        p2[i] = p1 + i * s2; \
    } \
    return p2; \
} \
  void free_##name##_4D_restricted(type **** data, size_t s4_lower) {			\
    if(data) { \
      free_##name##_3D_restricted(*data, s4_lower);	\
      free(data); \
    } \
}



#define CALLOC_5D_RESTRICTED(type, name) \
  type ***** calloc_##name##_5D_restricted(size_t s1, size_t s2, size_t s3, size_t s4_lower, size_t s4_upper, int *s5, const char* varname) { \
    type **** p1 = calloc_##name##_4D_restricted(s1 * s2, s3, s4_lower, s4_upper, s5, varname); \
    if(!p1) return NULL; \
    type ***** p2 = malloc(s1 * sizeof(type****)); \
    if(!p2) { \
        free_##name##_4D(p1); \
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname); \
        return NULL; \
    } \
    for(size_t i = 0; i < s1; ++i) { \
        p2[i] = p1 + i * s2; \
    } \
    return p2; \
} \
void free_##name##_5D_restricted(type ***** data, size_t s4_lower) { \
    if(data) { \
      free_##name##_4D_restricted(*data, s4_lower);	\
        free(data); \
    } \
}







#define CALLOC_ND(type, name)			\
CALLOC_1D(type, name) \
CALLOC_2D(type, name) \
CALLOC_3D(type, name) \
CALLOC_4D(type, name) \
CALLOC_5D(type, name) \
CALLOC_1D_RESTRICTED(type, name) \
CALLOC_2D_RESTRICTED(type, name) \
CALLOC_3D_RESTRICTED(type, name) \
CALLOC_4D_RESTRICTED(type, name) \
CALLOC_5D_RESTRICTED(type, name)

CALLOC_ND(char, char);
CALLOC_ND(unsigned char, uchar);
CALLOC_ND(int, int);
CALLOC_ND(float, float);
CALLOC_ND(double, double);

float *** calloc_float_3D_sparse(size_t nz, size_t nx, size_t ny, int *allocz, const char * varname) {
    size_t sz = 0;
    for(size_t i = 0; i < nz; ++i) {
        if(allocz[i]) ++sz;
    }
    float ** p1 = calloc_float_2D(sz * nx, ny, varname);
    if(!p1) return NULL;
    float *** p2 = malloc(nz * sizeof(float**));
    if(!p2) {
        free_float_2D(p1);
        fprintf (stderr, "Error allocating memory for 2st pointer level of %s\n", varname);
        return NULL;
    }
    size_t i_sparse = 0;
    for(size_t i = 0; i < nz; ++i) {
        if(allocz[i]) {
            p2[i] = p1 + i_sparse * nx;
            ++i_sparse;
        } else {
            p2[i] = NULL;
        }
    }
    return p2;
}
