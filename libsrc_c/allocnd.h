#ifndef ALLOCND_H
#define ALLOCND_H

#include <stdlib.h>

#define CALLOC_ND(type, name) \
type * calloc_##name##_1D(size_t sx, const char* varname); \
void free_##name##_1D(type * data); \
\
type ** calloc_##name##_2D(size_t sx, size_t sy, const char* varname); \
void free_##name##_2D(type ** data); \
\
type *** calloc_##name##_3D(size_t sx, size_t sy, size_t sz, const char* varname); \
void free_##name##_3D(type *** data); \
\
type **** calloc_##name##_4D(size_t sx, size_t sy, size_t sz, size_t dim4, const char* varname); \
void free_##name##_4D(type **** data); \
\
type ***** calloc_##name##_5D(size_t sx, size_t sy, size_t sz, size_t dim4, size_t dim5, const char* varname); \
void free_##name##_5D(type ***** data); \
\
type * calloc_##name##_1D_restricted(size_t sx, const char* varname); \
\
 type ** calloc_##name##_2D_restricted(size_t s1, size_t s2_lower, size_t s2_upper, int *s3, const char* varname); \
\
 type *** calloc_##name##_3D_restricted(size_t s1, size_t s2_lower, size_t s2_upper, int *s3, const char* varname); \
 void free_##name##_3D_restricted(type *** data, size_t s4_lower);  \
\
 type **** calloc_##name##_4D_restricted(size_t s1, size_t s2, size_t s3_lower, size_t s3_upper, int *s4, const char* varname); \
 void free_##name##_4D_restricted(type **** data, size_t s4_lower);  \
\
 type ***** calloc_##name##_5D_restricted(size_t s1, size_t s2, size_t s3, size_t s4_lower, size_t s4_upper, int *s5, const char* varname); \
 void free_##name##_5D_restricted(type ***** data, size_t s4_lower); 

CALLOC_ND(char, char);
CALLOC_ND(unsigned char, uchar);
CALLOC_ND(int, int);
CALLOC_ND(float, float);
CALLOC_ND(double, double);

float *** calloc_float_3D_sparse(size_t nz, size_t nx, size_t ny, int *allocz, const char * varname);

#undef CALLOC_ND

#endif
