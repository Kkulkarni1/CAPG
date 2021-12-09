/**
 * @file myfun.h
 * @author Yudi Zhang
 * @date 7/25/18
 *
 * PRINT AND INITIALIZE MATRIX VECTOR AND 3D ARRAY
 *
 * [KSD] See io.[ch] for such functionality.  Also, you should avoid requiring
 * files, such as haplotype_data.h, that limit the general functionality of
 * this code, which should ideally work for all vectors and matrices.
 *
 */

#ifndef MYFUN_H
#define MYFUN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "haplotype_data.h"    /* [KSD] You should avoid this. */

#define PRINT_FORMAT ".0lf"

#define PRINT_VECTOR(a,n) do {                                                 \
if ((a) != NULL) {                                                             \
size_t ARRAY_H1RESERVED;                                               \
for (ARRAY_H1RESERVED = 0; ARRAY_H1RESERVED < (n); ++ARRAY_H1RESERVED) \
printf("%"PRINT_FORMAT"\t", (a)[ARRAY_H1RESERVED] + 0.);       \
printf("\n");                                                  \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "             \
"out of memory!\n",  __FILE__, __func__, __LINE__);            \
}                                                                              \
} while (0)


#define PRINT_MATRIX(a,m,n) do {                                               \
if ((a) != NULL) {                                                             \
size_t ARRAY_H2RESERVED;                                               \
for (ARRAY_H2RESERVED = 0; ARRAY_H2RESERVED < (m); ++ARRAY_H2RESERVED) \
PRINT_VECTOR((a)[ARRAY_H2RESERVED], (n));                      \
printf("\n");                                                  \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "             \
"out of memory!\n",  __FILE__, __func__, __LINE__);            \
}                                                                              \
} while (0)

#define PRINT_3DARRAY(a,k,m,n) do {                                               \
if ((a) != NULL) {                                                             \
size_t ARRAY_H3RESERVED;                                               \
for (ARRAY_H3RESERVED = 0; ARRAY_H3RESERVED < (k); ++ARRAY_H3RESERVED) \
PRINT_MATRIX((a)[ARRAY_H3RESERVED],(m),(n));                      \
printf("\n");                                                  \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "             \
"out of memory!\n",  __FILE__, __func__, __LINE__);            \
}                                                                              \
} while (0)


#define SETZERO_VECTOR(a,n) do {                                               \
if ((a) != NULL) {                                                             \
size_t ARRAY_H1RESERVED;                                               \
for (ARRAY_H1RESERVED = 0; ARRAY_H1RESERVED < (n); ++ARRAY_H1RESERVED) \
(a)[ARRAY_H1RESERVED] = 0;                                     \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "             \
"out of memory!\n",  __FILE__, __func__, __LINE__);            \
}                                                                              \
} while (0)


#define SETZERO_MATRIX(a,m,n) do {                                             \
if ((a) != NULL) {                                                             \
size_t ARRAY_H2RESERVED;                                               \
for (ARRAY_H2RESERVED = 0; ARRAY_H2RESERVED < (m); ++ARRAY_H2RESERVED) \
SETZERO_VECTOR((a)[ARRAY_H2RESERVED], (n));                    \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "             \
"out of memory!\n",  __FILE__, __func__, __LINE__);            \
}                                                                              \
} while (0)


#define SETZERO_3DARRAY(a,k,m,n) do {                                          \
if ((a) != NULL) {                                                             \
size_t ARRAY_H3RESERVED;                                                \
for (ARRAY_H3RESERVED = 0; ARRAY_H3RESERVED < (k); ++ARRAY_H3RESERVED)  \
SETZERO_MATRIX((a)[ARRAY_H3RESERVED],(m),(n));                  \
} else {                                                                       \
fprintf(stderr, "*** in file %s, function %s(), line %d: "              \
"out of memory!\n",  __FILE__, __func__, __LINE__);             \
}                                                                               \
} while (0)

#endif /* MYFUN_H */

