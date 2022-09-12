/**
 * \file order.h
 * Functions to permute indices to order|sort arrays in as|decending order, like
 * <a href="http://www.r-project.org/">R</a>'s order() function.
 * For more information, see order.c.
 *
 * @author David Faden, dfaden@gmail.com
 * @author Karin S. Dorman kdorman@iastate.edu
 * @date August 21, 2005
 * @date Sun Sep 11 23:58:31 CDT 2011
 */

#ifndef __ORDER_H__
#define __ORDER_H__

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "constants.h"
#include "math.h"
#include "array.h"

size_t *qorder_int(const int* base, size_t numElements);
size_t *qorder_size_t(const size_t* base, size_t numElements);
size_t *qorder_double(const double* base, size_t numElements);
size_t *qorder_string(char* const* base, size_t numElements);

/**
 * Comparison function for two objects.
 * Examples of this type of function include #int_compare(), #double_compare(),
 * and #string_compare().
 * @param v1 void pointer to first object
 * @param v2 void pointer to second object
 * @return
 * 	-1 if v1 is "less than" v2,
 * 	0 if v1 "equals" v2, though see RETURN_RAND_CMP
 * 	1 if v1 is "greater than" v2.
 */
typedef int (*ComparisonFunc)(const void *v1, const void *v2);

/* some valid comparison functions */
int int_compare(const void *, const void *);
int size_t_compare(const void *, const void *);
int uint_compare(const void *, const void *);
int double_compare(const void *, const void *);
int string_compare(const void *, const void *);

/**
 * Comparison macro for numeric data types.
 * Returns -1 if a < b, 1 if a > b, and 0 otherwise from a function.
 * Useful to return from #ComparisonFunc and #CompareVectorElts functions.
 */
#define RETURN_CMP(a,b) if ((a) < (b)) {  return -1; } \
	else if ((a) == (b)) { return 0; }            \
	else { return 1; }

#define CMP(a,b) (((a) < (b)) ? -1 : ((a) == (b)) ? 0 : 1)

#define REV_RETURN_CMP(a,b) if ((a) < (b)) { return 1; }  \
	else if ((a) == (b)) { return 0; }                \
	else { return -1; }


/**
 * Random comparison macro for numeric data types.
 * Returns -1 if a < b, 1 if a > b, and -1 or 1 if a == b.
 * Useful to return from #ComparisonFunc and #CompareVectorElts functions.
 */
#define RETURN_RAND_CMP(a, b) if ((a) < (b)) { return -1; } \
	else if ((a) > (b)) { return 1;}                    \
	else if (rand()/RAND_MAX < 0.5) {return -1; }       \
	else { return 1; }


size_t* order(const void* base, size_t numElements, size_t size,
	      ComparisonFunc compare);


/**
 * Macro to define an #order() wrapper function.
 *
 * This macro takes a TYPE argument and defines a wrapper function named
 * order_TYPE.  The function takes two arguments, a const TYPE * pointer for an
 * array of TYPE elements, and a size_t, indicating the length of the array.
 * The function body defines a #ComparisonFunc-compliant function called compare_TYPE and
 * then calls the #order() function for arrays of type TYPE.  For example, the
 * call MAKE_ORDER_FUNC(int) would define a function order_int(const int *,
 * size_t), which is identical to #orderInt().

 *
 * @param TYPE data type for which you want an order function
 */
#define MAKE_ORDER_FUNC(TYPE) size_t* order_ ## TYPE (const TYPE* a, size_t len) \
{                                                                                \
	int compare_ ## TYPE (const void* v1, const void* v2) {                  \
		TYPE t1 = *((const TYPE*)v1);                                    \
		TYPE t2 = *((const TYPE*)v2);                                    \
		if (t1 < t2) { return -1; }                                      \
		else if (t1 == t2) { return 0; }                                 \
		else { return 1; }                                               \
	}                                                                        \
	return order(a, len, sizeof( TYPE ), compare_ ## TYPE );                 \
}

/**
 * Comparison function for comparing two elements in a vector.  Examples of
 * this type of function include #int_elt_compare, #size_t_elt_compare,
 * #int_elt_rev_compare, #size_t_elt_rev_compare, #double_elt_compare, and
 * #string_elt_compare.
 *
 * @param v vector of elements which are being sorted or ordered
 * @param i index of first element to compare
 * @param j index of second element to compare
 * @return -1, 0, 1 to indicate elt i is <, =, or > elt j.
 */
typedef int (*CompareEltsFunc)(const void *v, size_t i, size_t j);

/**
 * Swapper function for swapping two elements in a vector.  Examples of this
 * type of function include #int_swap, #size_t_swap, #double_swap, and
 * #string_swap.
 *
 * @param v vector of elements, two of which are to be swapped
 * @param i index of first element to be swapped
 * @param j index of second element to be swapped
 */
typedef void (*SwapperFunc)(void *v, size_t i, size_t j);

void quicksort(void *, CompareEltsFunc, SwapperFunc, size_t, size_t);
void quicksort_double_in_situ(double *, size_t, size_t);
void quickselect_double_in_situ(double *, size_t, size_t, size_t);

/**
 * But you do not need our implementation of quicksort.  Following
 * function is a demo of using C's qsort() to sort doubles.  You can call
 * these functions on your own: see definition in order.c.  All you need
 * is a #ComparisonFunc(), such as #double_compare
 */
void qsort_double_in_situ(double const * vec, size_t numElements);


/* some valid CompareEltsFunc functions */
int int_elt_compare(const void *, size_t, size_t);
int size_t_elt_compare(const void *, size_t, size_t);
int int_elt_rev_compare(const void *, size_t, size_t);		/* reverse order */
int size_t_elt_rev_compare(const void *, size_t, size_t);	/* reverse order */
int double_elt_compare(const void *, size_t, size_t);
int string_elt_compare(const void *, size_t, size_t);

/* numeric comparison macro */
#define COMPARE(a,b) ((a) < (b) ? -1 : (a) == (b) ? 0 : 1)
#define REV_COMPARE(a,b) ((a) < (b) ? 1 : (a) == (b) ? 0 : -1)

/* some valid SwapperFunc functions */
void int_swap(void *, size_t, size_t);
void size_t_swap(void *, size_t, size_t);
void double_swap(void *, size_t, size_t);
void string_swap(void *, size_t, size_t);

/* some swap macros */
#define SWAP_double(a,b) {register double temp=(a);(a)=(b);(b)=temp;}
#define SWAP_int(a,b) {register int temp=(a);(a)=(b);(b)=temp;}
#define SWAP_size_t(a,b) {register size_t temp=(a);(a)=(b);(b)=temp;}
#define SWAP_SIZE_T(a,b) {register SIZE_T temp=(a);(a)=(b);(b)=temp;}
#define SWAP_string(a,b) {register char *temp=(a);(a)=(b);(b)=temp;}

/* sorting that uses less memory, but requires C99 or greater standard */
#if (defined C99 || defined C11)

/**
 * Comparison function for elements within an array of simple data types.
 * Examples of this type of function include #compare_int_elts(),
 * #compare_double_elts(), and #compare_string_elts().
 *
 * @param v vector
 * @param index ordered indices (these are ordered rather than elts in v)
 * @param left left index
 * @param right right index
 * @return -1 if v[index[left]] < v[index[right]],
 *   1 if v[index[left]] > v[index[right]], 0 otherwise
 */
typedef int (*CompareVectorElts)(const void *v, size_t *index, size_t left, size_t right, va_list);

/* some order functions */
void with_index_quicksort(void *, size_t *, CompareVectorElts, size_t, size_t, ...);
void with_index_quickselect(void *, size_t *, CompareVectorElts, size_t, size_t, size_t, ...);
size_t *index_quicksort(void *, CompareVectorElts, size_t, size_t, ...);
size_t *order_int_quicksort(int *, size_t);
size_t *order_size_t_quicksort(size_t *, size_t);
size_t *order_double_quicksort(double *, size_t);
void order_double_quicksort_with_index(double *, size_t *, size_t);
void order_size_t_quicksort_with_index(size_t *, size_t *, size_t);
void select_int_with_index(int *, size_t *, size_t, size_t);
void select_double_with_index(double *, size_t *, size_t, size_t);
void reverse_select_double_with_index(double *, size_t *, size_t, size_t);
SIZE_T select_index(double *, SIZE_T, SIZE_T, SIZE_T);


/* some valid comparison functions */
int compare_int_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_size_t_elts(const void *, size_t *, size_t, size_t, va_list);
int reverse_compare_size_t_elts(const void *, size_t *, size_t, size_t, va_list);
int reverse_compare_uint_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_double_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_ulong_elts(const void *, size_t *, size_t, size_t, va_list);
int reverse_compare_double_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_string_elts(const void *, size_t *, size_t, size_t, va_list);

#endif

/* heap sort stuff */
size_t *heap_order_double(double *array, size_t size, size_t first_k);
int heap_order_with_index_double(double *array, size_t size, size_t first_k, size_t *sorted_index);
void heap_sort_double_in_situ(double *array, size_t size, size_t first_k);
size_t *heap_order_by_coord_double(double **array, size_t size, size_t first_k, size_t coord);

size_t *heap_order_int(int *array, size_t size, size_t first_k);
size_t *heap_order_intrev(int *array, size_t size, size_t first_k);
int heap_order_with_index_int(int *array, size_t size, size_t first_k, size_t *sorted_index);
void heap_sort_int_in_situ(int *array, size_t size, size_t first_k);

size_t *heap_order_size_t(size_t *array, size_t size, size_t first_k);
size_t *heap_order_size_trev(size_t *array, size_t size, size_t first_k);
int heap_order_with_index_size_t(size_t *array, size_t size, size_t first_k, size_t *sorted_index);
void heap_sort_size_t_in_situ(size_t *array, size_t size, size_t first_k);

#endif
