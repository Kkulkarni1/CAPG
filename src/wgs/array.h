/**
 * @file array.h
 * @author Rouben Rostamian, rostamian@umbc.edu
 * @author Ranjan Maitra, maitra@iastate.edu
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Thu Dec  6 08:29:14 CST 2012
 *
 * o This file defines the following macros:
 *
 *   - MAKE_1ARRAY(a,n)              make 1D array of length "n"
 *   - CMAKE_1ARRAY(a,n)             make 1D array of length "n", using calloc
 *   - REALLOC_1ARRAY(a,n,n1)        reallocate 1D array
 *   - MAKE_2ARRAY(a,m,n)            make 2D array of dims "m x n"
 *   - CMAKE_2ARRAY(a,m,n)           make 2D array of dims "m x n"
 *   - REALLOC_2ARRAY(a,m,n,m1,n1)   reallocate 2D array of dims "m1 x n1"
 *   - MAKE_3ARRAY(a,l,m,n)          make 3D array of dims "l x m x n"
 *   - CMAKE_3ARRAY(a,l,m,n)         make 3D array of dims "l x m x n"
 *   - REALLOC_3ARRAY(a,l,m,n,l1,m1,n1)
 *   - MAKE_4ARRAY(a,k,l,m,n)        make 4D array of dims "k x l x m x n"
 *   - CMAKE_4ARRAY(a,k,l,m,n)       make 4D array of dims "k x l x m x n"
 *   - MAKE_2JAGGED_ARRAY(a,m,v)     make 2D jagged array; v is array of len m
 *   - CMAKE_2JAGGED_ARRAY(a,m,v)    make 2D jagged array; v is array of len m
 *   - MAKE_3JAGGED_ARRAY(a,m,n,v)   make 3D jagged matrix; v is array of len n
 *   - CMAKE_3JAGGED_ARRAY(a,m,n,v)  make 3D jagged matrix; v is array of len n
 *   - COPY_1ARRAY(a,b,n)            copy 1D array of length "n"
 *   - COPY_2ARRAY(a,b,n,m)          copy 2D array of length "m x n"
 *   - COPY_3ARRAY(a,b,l,m,n)        make 3D array of dimensions "l x m x n"
 *   - COPY_4ARRAY(a,b,k,l,m,n)      make 4D array of dimensions "k x l x m x n"
 *   - COPY_2JAGGED_ARRAY(a,b,v)     make 2D jagged array; v is array of len m
 *   - COPY_3JAGGED_ARRAY(a,b,v)     make 3D jagged matrix; v is array of len n
 *
 *   - FREE_1ARRAY(a)         free memory allocated by MAKE_1ARRAY()
 *   - FREE_2ARRAY(a)         free memory allocated by MAKE_2ARRAY()
 *   - FREE_3ARRAY(a)         free memory allocated by MAKE_3ARRAY()
 *   - FREE_4ARRAY(a)         free memory allocated by MAKE_4ARRAY()
 *
 * o Additionally, it defines the following convenience macros as synonyms:
 *
 *   - MAKE_VECTOR(a,n)       same as MAKE_1ARRAY(a,n)
 *   - CMAKE_VECTOR(a,n)      same as CMAKE_1ARRAY(a,n)
 *   - REALLOC_VECTOR(a,n,n1) same as REALLOC_1ARRAY(a,n,n1)
 *   - FREE_VECTOR(a)         same as FREE_1ARRAY(a)
 *   - MAKE_MATRIX(a,m,n)     same as MAKE_2ARRAY(a,m,n)
 *   - CMAKE_MATRIX(a,m,n)    same as CMAKE_2ARRAY(a,m,n)
 *   - FREE_MATRIX(a)         same as FREE_2ARRAY(a)
 *   - REALLOC_MATRIX(a,m,n,m1,n1)  same as REALLOC_2ARRAY(a,m,n,m1,na)
 *
 * o Additionally, it declares and uses the identifiers:
 *
 *   - ARRAY_H2RESERVED
 *   - ARRAY_H3RESERVED
 *   - ARRAY_H4RESERVED
 *
 *   within local blocks within the macros.  THESE IDENTIFIERS
 *   ARE RESERVED.  The user should not use identifiers with this
 *   names within the scope of these macros.
 *
 * o vector/matrix/array elements can be of _any_ type.
 *
 * o If malloc() fails during the execution of any of the MAKE_* macros:
 *     - An "out of memory" message is printed to stderr.
 *     - The macro's first argument is set to NULL.
 *
 * o After a call to any of the FREE_*() macros, the macro's argument
 *   is set to NULL.
 *
 * o The FREE_*() macros can be applied to previously FREE_*()-ed
 *   arrays with no ill effect.
 *
 * o Note that macro arguments are evaluated more than once.
 *
 * o The current versions lead to fragmented memory.  TODO write other versions
 *   that distribute memory differently.
 *
 * o Customization
 *
 *    When malloc() returns NULL, MAKE_1ARRAY prints a message on stderr
 *    and the program continues.  The user can alter this behavior by
 *    redefining the MAKE_1ARRAY macro.  For instance, to call exit()
 *    whenever malloc() fails, the user can do:
 *
 *    #undef MAKE_1ARRAY
 *    #define MAKE_1ARRAY(a,n) do {                                          \
 *        (a) = malloc((n) * sizeof *(a));                                   \
 *        if ((a)==NULL) {                                                   \
 *            fprintf(stderr, "*** in file %s, function %s(), line %d: "     \
 *                    "out of memory!\n",  __FILE__, __func__, __LINE__);    \
 *            exit(EXIT_FAILURE);                                            \
 *        }                                                                  \
 *    } while (0)
 *
 *    Since only MAKE_1ARRAY calls malloc() explicitly, this change affects
 *    the behavior of not only MAKE_1ARRAY but all other MAKE_* macros as well.
 *
 *
 * ---- SAMPLE USAGE -------------
 *#include "array.h"
 *int main(void)
 *{
 *    float ***a;            // can use any other type instead of "float"
 *    size_t p=3, q=4, r=5;  // will make a 3x4x5 3-D array of float
 *
 *    MAKE_3ARRAY(a, p, q, r);
 *    if (a==NULL)
 *        return EXIT_FAILURE;
 *
 *    a[2][0][1] = 3.14;
 *    printf("%g \n", a[2][0][1]);
 *    FREE_3ARRAY(a);
 *    return EXIT_SUCCESS;
 *}
 *---- END OF SAMPLE USAGE -------
 *
 *
 *   - June 2000
 *   - Revised Jul 2000
 *   - Revised Sep 2002
 *   - Revised Jun 2003
 *     . macros don't call exit anymore
 *     . changed macro names to all-caps
 *   - Revised Aug 2003
 *     . changed reserved names to all caps
 *   - Revised Dec 2003
 *     . changed array index types to size_t (were int)
 *     . added an "out of memory" message, printed to stderr
 *     . most FREE* macros now come before MAKE* macros,
 *       for possible improved efficiency in preprocessing
 *   - Revised Dec 2012
 *     . added CMAKE_* versions calling calloc() instead of malloc()
 *     . added *MAKE_*JAGGED_ARRAY versions for creating jagged arrays where
 *       last dimension is variable, depending on index of penultimate
 *       dimension
 *     . added *COPY*ARRAY macros for copying arays using efficient memcpy()
 *       calls
 */


#ifndef H_ARRAY_
#define H_ARRAY_

#include <stdio.h>
#include <stdlib.h>

/* pointer to desired MAKE_1ARRAY variant */
#define MAKE_1ARRAY MAKE_1ARRAY_DEFAULT
#define CMAKE_1ARRAY CMAKE_1ARRAY_DEFAULT
#define REALLOC_1ARRAY REALLOC_1ARRAY_DEFAULT

/* ---------- 1D arrays ---------------------- */

#define MAKE_1ARRAY_DEFAULT(a,n) do {                                        \
    (a) = malloc((n) * sizeof *(a));                                         \
    if ((a)==NULL)                                                           \
        fprintf(stderr, "*** in file %s, function %s(), line %d: "           \
                "out of memory!\n",  __FILE__, __func__, __LINE__);          \
} while (0)

#define CMAKE_1ARRAY_DEFAULT(a,n) do {                                       \
    (a) = calloc((n), sizeof *(a));                                          \
    if ((a)==NULL)                                                           \
        fprintf(stderr, "*** in file %s, function %s(), line %d: "           \
                "out of memory!\n",  __FILE__, __func__, __LINE__);          \
} while (0)

#define REALLOC_1ARRAY_DEFAULT(a,n,n1) do {                                  \
	void *RESERVED_PTR;                                                  \
	if ((n) == (n1))                                                     \
		break;                                                       \
	RESERVED_PTR = realloc((a), (n1) * sizeof *(a));                     \
	if (RESERVED_PTR == NULL)                                            \
		fprintf(stderr, "*** in file %s, function %s(), line %d: "   \
			"out of memory!\n",  __FILE__, __func__, __LINE__);  \
	else                                                                 \
		(a) = RESERVED_PTR;                                          \
} while (0)

#define COPY_1ARRAY(a,b,n) do {                                                \
	(a) = memcpy((a), (b), (n) * sizeof *(a));                             \
	if ((a) == NULL)                                                       \
		fprintf(stderr, "*** in file %s, function %s(), line %d: "     \
			"memcpy failed!\n", __FILE__, __func__, __LINE__);     \
} while (0)

#define FREE_1ARRAY(a)  do {                                                 \
    free(a);                                                                 \
    a = NULL;                                                                \
} while (0)

/* ---------- 2D arrays ---------------------- */

/* note: parenthesize first arg because it may be given as `*a' */
#define MAKE_2ARRAY(a,m,n) do {                                              \
    size_t ARRAY_H2RESERVED;                                                 \
    MAKE_1ARRAY(a,(m)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[m] = NULL;                                                           \
    for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(size_t)(m);                   \
	ARRAY_H2RESERVED++) {                                                \
        MAKE_1ARRAY((a)[ARRAY_H2RESERVED],(n));                              \
        if ((a)[ARRAY_H2RESERVED]==NULL) {                                   \
            FREE_2ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define CMAKE_2ARRAY(a,m,n) do {                                             \
    size_t ARRAY_H2RESERVED;                                                 \
    MAKE_1ARRAY((a), (m)+1);                                                 \
    if ((a) == NULL)                                                         \
        break;                                                               \
    (a)[m] = NULL;                                                           \
    for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(size_t)(m);                   \
	ARRAY_H2RESERVED++) {                                                \
        CMAKE_1ARRAY((a)[ARRAY_H2RESERVED],(n));                             \
        if ((a)[ARRAY_H2RESERVED]==NULL) {                                   \
            FREE_2ARRAY((a));                                                \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define REALLOC_2ARRAY(a,m,n,m1,n1) do {                                     \
	size_t ARRAY_H2RESERVED;                                             \
	if ((m) != (m1)) {                                                   \
		REALLOC_1ARRAY((a), (m) + 1, (m1) + 1);                      \
		if ((a) == NULL)                                             \
			break;                                               \
		(a)[m1] = NULL;                                              \
	}                                                                    \
	if ((n) != (n1)) {                                                   \
		for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(size_t)(m1);      \
			ARRAY_H2RESERVED++) {                                \
			REALLOC_1ARRAY((a)[ARRAY_H2RESERVED], n, n1);        \
			if ((a)[ARRAY_H2RESERVED] == NULL) {                 \
				FREE_2ARRAY((a));                            \
				break;                                       \
			}                                                    \
		}                                                            \
	}                                                                    \
} while (0)

#define MAKE_2JAGGED_ARRAY(a,m,v) do {                                         \
	size_t ARRAY_H2RESERVED;                                               \
	MAKE_1ARRAY(a, (m)+1);                                                 \
	if ((a) == NULL)                                                       \
		break;                                                         \
	(a)[m] = NULL;                                                         \
	for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(size_t)(m);                 \
		ARRAY_H2RESERVED++) {                                          \
		MAKE_1ARRAY((a)[ARRAY_H2RESERVED], (v)[ARRAY_H2RESERVED]);     \
		if ((a)[ARRAY_H2RESERVED] == NULL) {                           \
			FREE_2ARRAY(a);                                        \
			break;                                                 \
		}                                                              \
	}                                                                      \
} while (0)

#define CMAKE_2JAGGED_ARRAY(a,m,v) do {                                        \
	size_t ARRAY_H2RESERVED;                                               \
	MAKE_1ARRAY(a, (m)+1);                                                 \
	if ((a) == NULL)                                                       \
		break;                                                         \
	(a)[m] = NULL;                                                         \
	for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(size_t)(m);                 \
		ARRAY_H2RESERVED++) {                                          \
		CMAKE_1ARRAY((a)[ARRAY_H2RESERVED], (v)[ARRAY_H2RESERVED]);    \
		if ((a)[ARRAY_H2RESERVED] == NULL) {                           \
			FREE_2ARRAY(a);                                        \
			break;                                                 \
		}                                                              \
	}                                                                      \
} while (0)

#define COPY_2ARRAY(a,b,n) do {                                                \
	size_t ARRAY_H2RESERVED;                                               \
	for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL;                  \
		ARRAY_H2RESERVED++) {                                          \
		COPY_1ARRAY((a)[ARRAY_H2RESERVED], (b)[ARRAY_H2RESERVED], (n));\
		if ((a)[ARRAY_H2RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define COPY_2JAGGED_ARRAY(a,b,v) do {                                         \
	size_t ARRAY_H2RESERVED;                                               \
	for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL;                  \
		ARRAY_H2RESERVED++) {                                          \
		COPY_1ARRAY((a)[ARRAY_H2RESERVED], (b)[ARRAY_H2RESERVED],      \
			(v)[ARRAY_H2RESERVED]);                                \
		if ((a)[ARRAY_H2RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define FREE_2ARRAY(a) do {                                                  \
    size_t ARRAY_H2RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL; ARRAY_H2RESERVED++)\
        FREE_1ARRAY((a)[ARRAY_H2RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)


/* ---------- 3D arrays ---------------------- */

#define MAKE_3ARRAY(a,p,q,r) do {                                            \
    size_t ARRAY_H3RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(size_t)(p);                   \
	ARRAY_H3RESERVED++) {                                                \
        MAKE_2ARRAY((a)[ARRAY_H3RESERVED],(q),(r));                          \
        if ((a)[ARRAY_H3RESERVED]==NULL) {                                   \
            FREE_3ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define CMAKE_3ARRAY(a,p,q,r) do {                                           \
    size_t ARRAY_H3RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(size_t)(p);                   \
	ARRAY_H3RESERVED++) {                                                \
        CMAKE_2ARRAY((a)[ARRAY_H3RESERVED],(q),(r));                         \
        if ((a)[ARRAY_H3RESERVED]==NULL) {                                   \
            FREE_3ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define REALLOC_3ARRAY(a,l,m,n,l1,m1,n1) do {                                \
	size_t ARRAY_H3RESERVED;                                             \
	if ((l) != (l1)) {                                                   \
		REALLOC_1ARRAY((a), (l) + 1, (l1) + 1);                      \
		if ((a) == NULL)                                             \
			break;                                               \
		(a)[l1] = NULL;                                              \
	}                                                                    \
	if ((m) != (m1) || (n) != (n1)) {                                    \
		for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(size_t)(l1);      \
			ARRAY_H3RESERVED++) {                                \
			REALLOC_2ARRAY((a)[ARRAY_H3RESERVED], n, m, n1, m1); \
			if ((a)[ARRAY_H3RESERVED] == NULL) {                 \
				FREE_2ARRAY((a));                            \
				break;                                       \
			}                                                    \
		}                                                            \
	}                                                                    \
} while (0)

#define MAKE_3JAGGED_ARRAY(a,p,q,v) do {                                       \
	size_t ARRAY_H3RESERVED;                                               \
	MAKE_1ARRAY((a), (p)+1);                                               \
	if ((a) == NULL)                                                       \
		break;                                                         \
	(a)[p] = NULL;                                                         \
	for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(size_t)(p);                 \
		ARRAY_H3RESERVED++) {                                          \
		MAKE_2JAGGED_ARRAY((a)[ARRAY_H3RESERVED], (q), (v));           \
		if ((a)[ARRAY_H3RESERVED] == NULL) {                           \
			FREE_3ARRAY(a);                                        \
			break;                                                 \
		}                                                              \
	}                                                                      \
} while (0)

#define CMAKE_3JAGGED_ARRAY(a,p,q,v) do {                                      \
	size_t ARRAY_H3RESERVED;                                               \
	MAKE_1ARRAY((a), (p)+1);                                               \
	if ((a) == NULL)                                                       \
		break;                                                         \
	(a)[p] = NULL;                                                         \
	for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(size_t)(p);                 \
		ARRAY_H3RESERVED++) {                                          \
		CMAKE_2JAGGED_ARRAY((a)[ARRAY_H3RESERVED], (q), (v));          \
		if ((a)[ARRAY_H3RESERVED] == NULL) {                           \
			FREE_3ARRAY(a);                                        \
			break;                                                 \
		}                                                              \
	}                                                                      \
} while (0)

#define COPY_3ARRAY(a,b,r) do {                                                \
	size_t ARRAY_H3RESERVED;                                               \
	for (ARRAY_H3RESERVED=0; (a)[ARRAY_H3RESERVED]!=NULL;                  \
		ARRAY_H3RESERVED++) {                                          \
		COPY_2ARRAY((a)[ARRAY_H3RESERVED], (b)[ARRAY_H3RESERVED], (r));\
		if ((a)[ARRAY_H3RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define COPY_3JAGGED_ARRAY(a,b,v) do {                                         \
	size_t ARRAY_H3RESERVED;                                               \
	for (ARRAY_H3RESERVED=0; (a)[ARRAY_H3RESERVED]!=NULL;                  \
		ARRAY_H3RESERVED++) {                                          \
		COPY_2JAGGED_ARRAY((a)[ARRAY_H3RESERVED],                      \
			(b)[ARRAY_H3RESERVED], (v));                           \
		if ((a)[ARRAY_H3RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define FREE_3ARRAY(a) do {                                                  \
    size_t ARRAY_H3RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H3RESERVED=0; (a)[ARRAY_H3RESERVED]!=NULL; ARRAY_H3RESERVED++)\
        FREE_2ARRAY((a)[ARRAY_H3RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

/* ---------- 4D arrays ---------------------- */

#define MAKE_4ARRAY(a,p,q,r,s) do {                                          \
    size_t ARRAY_H4RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(size_t)(p);                   \
	ARRAY_H4RESERVED++) {                                                \
        MAKE_3ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(s));                      \
        if ((a)[ARRAY_H4RESERVED]==NULL) {                                   \
            FREE_4ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define CMAKE_4ARRAY(a,p,q,r,s) do {                                         \
	size_t ARRAY_H4RESERVED;                                             \
	MAKE_1ARRAY(a,(p)+1);                                                \
	if (a==NULL)                                                         \
		break;                                                       \
	(a)[p] = NULL;                                                       \
	for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(size_t)(p);               \
		ARRAY_H4RESERVED++) {                                        \
		CMAKE_3ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(s));             \
		if ((a)[ARRAY_H4RESERVED]==NULL) {                           \
			FREE_4ARRAY(a);                                      \
			break;                                               \
		}                                                            \
	}                                                                    \
} while (0)

#define MAKE_4JAGGED_ARRAY(a,p,q,r,v) do {                                   \
    size_t ARRAY_H4RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(size_t)(p);                   \
	ARRAY_H4RESERVED++) {                                                \
        MAKE_3JAGGED_ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(v));               \
        if ((a)[ARRAY_H4RESERVED]==NULL) {                                   \
            FREE_4ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define CMAKE_4JAGGED_ARRAY(a,p,q,r,v) do {                                  \
    size_t ARRAY_H4RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(size_t)(p);                   \
	ARRAY_H4RESERVED++) {                                                \
        CMAKE_3JAGGED_ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(v));              \
        if ((a)[ARRAY_H4RESERVED]==NULL) {                                   \
            FREE_4ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define COPY_4ARRAY(a,b,s) do {                                                \
	size_t ARRAY_H4RESERVED;                                               \
	for (ARRAY_H4RESERVED=0; (a)[ARRAY_H4RESERVED]!=NULL;                  \
		ARRAY_H4RESERVED++) {                                          \
		COPY_3ARRAY((a)[ARRAY_H4RESERVED], (b)[ARRAY_H4RESERVED], (s));\
		if ((a)[ARRAY_H4RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define COPY_4JAGGED_ARRAY(a,b,v) do {                                         \
	size_t ARRAY_H4RESERVED;                                               \
	for (ARRAY_H4RESERVED=0; (a)[ARRAY_H4RESERVED]!=NULL;                  \
		ARRAY_H4RESERVED++) {                                          \
		COPY_3JAGGED_ARRAY((a)[ARRAY_H4RESERVED],                      \
			(b)[ARRAY_H4RESERVED], (v));                           \
		if ((a)[ARRAY_H4RESERVED] == NULL)                             \
			break;                                                 \
	}                                                                      \
} while (0)

#define FREE_4ARRAY(a) do {                                                  \
    size_t ARRAY_H4RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H4RESERVED=0; (a)[ARRAY_H4RESERVED]!=NULL; ARRAY_H4RESERVED++)\
        FREE_3ARRAY((a)[ARRAY_H4RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

/* ---------- synonyms ---------------------- */

#define MAKE_VECTOR(a,n)		MAKE_1ARRAY(a,n)
#define CMAKE_VECTOR(a,n)		CMAKE_1ARRAY(a,n)
#define REALLOC_VECTOR(a,n,n1)		REALLOC_1ARRAY(a,n,n1)
#define MAKE_MATRIX(a,m,n)		MAKE_2ARRAY(a,m,n)
#define CMAKE_MATRIX(a,m,n)		CMAKE_2ARRAY(a,m,n)
#define REALLOC_MATRIX(a,m,n,m1,n1)	REALLOC_2ARRAY(a,m,n,m1,n1)

#define FREE_VECTOR(a)      FREE_1ARRAY(a)
#define FREE_MATRIX(a)      FREE_2ARRAY(a)

#endif /* H_ARRAY_ */
