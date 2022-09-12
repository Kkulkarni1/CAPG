/**
 * \file order.c
 *
 * Routines that return a permutation of indices to put an array in sorted
 * order, like <a href="http://www.r-project.org">R</a>'s order function.
 *
 * For arrays of simple data types, see #int_order_simple() for integers and
 * #index_quicksort() for any array of simple data types, provided a
 * #CompareVectorElts function is provided.  You can write your own
 * #CompareVectorElts function, but see #compare_int_elts,
 * #compare_double_elts, and #compare_string_elts for comparing
 * elements of an integer, double, or const string array.
 *
 * See #order() for a generic function that can order indices for any type of
 * object for which the user can provide a #ComparisonFunc.  You can write
 * your own #ComparisonFunc function, but see #int_compare, #double_compare, and
 * #string_compare for comparing int, double, or strings.  The functions
 * #qorder_int(), #qorder_double(), and #qorder_string() are wrapper functions
 * for #order() using these #ComparisonFunc functions, respectively.
 *
 * @author David Faden, dfaden@gmail.com
 * @author Karin S. Dorman, kdorman@iastate.edu
 * @date August 21, 2005
 * @date Sun Sep 11 14:43:36 CDT 2011
 */

#include "order.h"

/*** BEGIN: Sorting with C library quicksort qsort() ***/

/**
 * Pair struct for agnostic sorting of objects.
 */
typedef struct {
	size_t index;		/*!< index of object in array */
	const void* datum;	/*!< object in array */
	ComparisonFunc compare;	/*!< pointer to comparison function */
} Pair;

/**
 * Compare two #Pair objects using user-provided comparison function.
 *
 * Generic comparison function for two objects.  The arguments are of type
 * #Pair, which provides a comparison function for doing the actual comparison.
 * If the comparison functions of the two objects match, then use it to compare
 * the two objects.
 *
 * @param v1 pointer to first object
 * @param v2 pointer to second object
 * @return -1 if v1 < v2, 0 if v1 == v2, and 1 otherwise
 */
int comparePairs(const void* v1, const void* v2)
{
	const Pair* p1 = v1;
	const Pair* p2 = v2;
	assert(p1->compare == p2->compare);
	return p1->compare(p1->datum, p2->datum);
} /* comparePairs */

/**
 * Function returning index order of an arbitrary array of objects.
 *
 * This function uses the C library function qsort to find index order
 * that would sort arbitrary objects in increasing order according to a
 * user-defined comparison.  This function is very flexible with regards
 * to the types of objects it can sort, but requires considerable memory.
 * In particular it needs to allocate two pointers and a size_t index
 * for each element in the array.  See #order_int_quicksort() for
 * index ordering of simple integer arrays or #index_quicksort() for
 * index ordering of arrays of an arbitrary, simple data type.
 *
 * The caller is responsible for freeing the returned object.
 *
 * @param base array of objects
 * @param numElements length of array
 * @param size size in bytes of elements in array
 * @param compare a function that is capable of comparing objects in the array
 * @return indices of array in order to produce objects in increasing order
 */
size_t *order(const void* base, size_t numElements, size_t size,
				ComparisonFunc compare)
{
	Pair *pairs = malloc(sizeof(Pair) * numElements);
	size_t *indices = malloc(sizeof(size_t) * numElements);
	size_t i;
	const char *p = base;

	for (i = 0; i < numElements; ++i, p += size) {
		pairs[i].index = i;
		pairs[i].datum = p;
		pairs[i].compare = compare;
	}

	qsort(pairs, numElements, sizeof(Pair), comparePairs);

	for (i = 0; i < numElements; ++i)
		indices[i] = pairs[i].index;

	free(pairs);

	return indices;
} /* order */


/**
 * A #ComparisonFunc-compliant function for two integers.
 * @param v1 pointer to first integer
 * @param v2 pointer to second integer
 * @return comparison value, see #RETURN_CMP.
 */
int int_compare(const void* v1, const void* v2)
{
	int i1 = *((const int*)v1);
	int i2 = *((const int*)v2);
	RETURN_CMP(i1, i2);
} /* int_compare */

/**
 * A #ComparisonFunc-compliant function for two unsigned int integers.
 * @param v1 pointer to first unsigned int
 * @param v2 pointer to second unsigned int
 * @return comparison value, see #RETURN_CMP.
 */
int uint_compare(const void* v1, const void* v2)
{
	unsigned int i1 = *((unsigned int const *)v1);
	unsigned int i2 = *((unsigned int const *)v2);
	RETURN_CMP(i1, i2);
} /* uint_compare */

/**
 * A #ComparisonFunc-compliant function for two unsigned int integers.
 * @param v1 pointer to first unsigned int
 * @param v2 pointer to second unsigned int
 * @return comparison value, see #RETURN_CMP.
 */
int rev_uint_compare(const void* v1, const void* v2)
{
	unsigned int i1 = *((unsigned int const *)v1);
	unsigned int i2 = *((unsigned int const *)v2);
	REV_RETURN_CMP(i1, i2);
} /* rev_uint_compare */

/**
 * A #ComparisonFunc-compliant function for two size_t integers.
 * @param v1 pointer to first size_t
 * @param v2 pointer to second size_t
 * @return comparison value, see #RETURN_CMP.
 */
int size_t_compare(const void* v1, const void* v2)
{
	size_t i1 = *((const size_t*)v1);
	size_t i2 = *((const size_t*)v2);
	RETURN_CMP(i1, i2);
} /* size_t_compare */

/**
 * A #ComparisonFunc-compliant function for two size_t integers.
 * @param v1 pointer to first unsigned int
 * @param v2 pointer to second unsigned int
 * @return comparison value, see #RETURN_CMP.
 */
int rev_size_t_compare(const void* v1, const void* v2)
{
        size_t i1 = *((unsigned int const *)v1);
        size_t i2 = *((unsigned int const *)v2);
        REV_RETURN_CMP(i1, i2);
} /* rev_size_t_compare */

/**
 * A #ComparisonFunc-compliant function for two doubles.
 * @param v1 pointer to first double
 * @param v2 pointer to second double
 * @return comparison value, see #RETURN_CMP.
 */
int double_compare(const void* v1, const void* v2)
{
	double d1 = *((const double*)v1);
	double d2 = *((const double*)v2);
	RETURN_CMP(d1, d2);
} /* double_compare */

/**
 * A #ComparisonFunc-compliant function for two strings.
 * @param v1 pointer to first string
 * @param v2 pointer to second string
 * @return comparison value, see #RETURN_CMP.
 */
int string_compare(const void* v1, const void* v2)
{
	return strcmp(*(const char**)v1, *(const char**)v2);
} /* string_compare */

/**
 * Wrapper for #order() to produce index order of integer array.
 * @param base array of integers
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* qorder_int(const int* base, size_t numElements)
{
	return order(base, numElements, sizeof(int), int_compare);
} /* qorder_int */

/**
 * Wrapper for #order() to produce index order of size_t array.
 * @param base array of size_t
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* qorder_size_t(const size_t* base, size_t numElements)
{
	return order(base, numElements, sizeof(size_t), size_t_compare);
} /* qorder_size_t */

/**
 * Wrapper for #order() to produce index order of double array.
 * @param base array of doubles
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* qorder_double(const double* base, size_t numElements)
{
	return order(base, numElements, sizeof(double), double_compare);
} /* qorder_double */

/**
 * Wrapper for #order() to produce index order of string array.
 * @param base array of strings
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* qorder_string(char* const* base, size_t numElements)
{
	return order(base, numElements, sizeof(const char*), string_compare);
} /* qorder_string */

/**
 * Example wrapper for in situ sort of vector using C library qsort().
 *
 * @param vec vector to sort
 * @param size length of vector
 */
void qsort_double_in_situ(double const *vec, size_t size)
{
	qsort((void *) vec, size, sizeof(double), double_compare);
} /* qsort_double_in_situ */

/*** END: Sorting with C library quicksort qsort() ***/



/*** BEGIN: our first implementation of quicksort ***/

#define SWAP(a, b, c) do { \
	(c) = (a);         \
	(a) = (b);         \
	(b) = (c);         \
} while (0);

/* [TODO] should inline these */
void int_swap(void *in_v, size_t i, size_t j)
{
	int tmp;
	int *v = (int *) in_v;
	SWAP(v[i], v[j], tmp);
} /* int_swap */

void size_t_swap(void *in_v, size_t i, size_t j)
{
	size_t tmp;
	size_t *v = (size_t *) in_v;
	SWAP(v[i], v[j], tmp);
} /* size_t_swap */

void double_swap(void *in_v, size_t i, size_t j)
{
	double tmp;
	double *v = (double *) in_v;
	SWAP(v[i], v[j], tmp);
} /* double_swap */

void string_swap(void *in_v, size_t i, size_t j)
{
	char *tmp;
	char **v = (char **) in_v;
	SWAP(v[i], v[j], tmp);
} /* string_swap */

int int_elt_compare(const void *in_v, size_t i, size_t j)
{
	const int *v = (const int *) in_v;
	RETURN_CMP(v[i], v[j]);
} /* int_elt_compare */

int size_t_elt_compare(const void *in_v, size_t i, size_t j)
{
	const size_t *v = (const size_t *) in_v;
	RETURN_CMP(v[i], v[j]);
} /* size_t_elt_compare */

int int_elt_rev_compare(const void *in_v, size_t i, size_t j)
{
	const int *v = (const int *) in_v;
	REV_RETURN_CMP(v[i], v[j]);
} /* int_elt_compare */

int size_t_elt_rev_compare(const void *in_v, size_t i, size_t j)
{
	const size_t *v = (const size_t *) in_v;
	REV_RETURN_CMP(v[i], v[j]);
} /* size_t_elt_compare */

int double_elt_compare(const void *in_v, size_t i, size_t j)
{
	const double *v = (const double *) in_v;
	RETURN_CMP(v[i], v[j]);
} /* double_elt_compare */

int string_elt_compare(const void *in_v, size_t i, size_t j)
{
	const char **v = (const char **) in_v;
	return strcmp(v[i], v[j]);
} /* string_elt_compare */

size_t vpartition(void *vec, CompareEltsFunc compare, SwapperFunc swapper, size_t left, size_t right, size_t pivot);

void
quicksort(void *vec, CompareEltsFunc compare, SwapperFunc swapper,
	size_t left, size_t right)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;

	if (left < right) {
		middle = (right - left)/2 + left;
		compare1 = compare(vec, left, middle);	/* !!! */
		compare2 = compare(vec, middle, right);
		compare3 = compare(vec, left, right);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		pivot = vpartition(vec, compare, swapper, left, right, pivot);
		if (pivot)
			quicksort(vec, compare, swapper, left, pivot - 1);
		quicksort(vec, compare, swapper, pivot + 1, right);
	}
} /* quicksort */

size_t
vpartition(void *vec, CompareEltsFunc compare, SwapperFunc swapper,
	size_t left, size_t right, size_t pivot)
{
	size_t i, new_pivot;

	/* move pivot to end of array */
	swapper(vec, right, pivot);

	/* sort this partition */
	new_pivot = left;
	for (i = left; i < right; i++) {
		if (compare(vec, i, right) < 0) {
			swapper(vec, new_pivot, i);
			new_pivot++;
		}
	}

	swapper(vec, new_pivot, right);

	return new_pivot;
} /* vpartition */

size_t vpartition_double_in_situ(double *vec, size_t left, size_t right, size_t pivot);

/**
 * Our implementation of quick sort for double arrays.  This can also be done
 * with the generic #quicksort(), but we wrote this to test the cost of using
 * function pointers.
 *
 * @param vec vector of doubles to sort
 * @param left left index of portion to sort
 * @param right right index of portion to sort
 */
void
quicksort_double_in_situ(double *vec, size_t left, size_t right)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;

	if (left < right) {
		middle = (right - left)/2 + left;
		compare1 = COMPARE(vec[left], vec[middle]);	/* !!! */
		compare2 = COMPARE(vec[middle], vec[right]);
		compare3 = COMPARE(vec[left], vec[right]);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		pivot = vpartition_double_in_situ(vec, left, right, pivot);
		if (pivot)
			quicksort_double_in_situ(vec, left, pivot - 1);
		quicksort_double_in_situ(vec, pivot + 1, right);
	}
} /* quicksort_double_in_situ */

/**
 * Quick select for double arrays.  Find the k smallest doubles in a vector.
 *
 * @param vec array of doubles to partially sort
 * @param left leftmost index (inclusive) to sort
 * @param right rightmost index (inclusive) to sort
 * @param k sort up to kth smallest element
 */
void
quickselect_double_in_situ(double *vec, size_t left, size_t right, size_t k)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;

	if (left < right) {
		middle = (right - left)/2 + left;
		compare1 = COMPARE(vec[left], vec[middle]);	/* !!! */
		compare2 = COMPARE(vec[middle], vec[right]);
		compare3 = COMPARE(vec[left], vec[right]);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		pivot = vpartition_double_in_situ(vec, left, right, pivot);
		if (k < pivot)
			quickselect_double_in_situ(vec, left, pivot - 1, k);
		else if (k > pivot)
			quickselect_double_in_situ(vec, pivot + 1, right, k);
	}
} /* quickselect_double_in_situ */

size_t
vpartition_double_in_situ(double *vec, size_t left, size_t right, size_t pivot)
{
	size_t i, new_pivot;

	/* move pivot to end of array */
	SWAP_double(vec[right], vec[pivot]);

	/* sort this partition */
	new_pivot = left;
	for (i = left; i < right; i++) {
		if (COMPARE(vec[i], vec[right]) < 0) {
			SWAP_double(vec[new_pivot], vec[i]);
			new_pivot++;
		}
	}

	SWAP_double(vec[new_pivot], vec[right]);

	return new_pivot;
} /* vpartition_double_in_situ */


/*** END: our own implementation of quicksort ***/


/*** BEGIN: heap sort algorithms ***/

/* Some testing reveals that sorting a copy of the array in heap sort is faster
 * than just sorting the index.  In addition, our implementation of quick sort
 * is slower than heap sort on randomly ordered arrays, I guess because of the
 * recursive structure of the quick sort.  Because we may someday want heap sort
 * that does not make a copy of the original data, even though it may be
 * slower, there is commented code below indicating how one can adjust to
 * avoid ordering on a parallel copy of the data. [KSD]
 */

/* These functions are in groups that work together:
 * - insert_elt_*, delete_min_*, heap_order_*_internal, and heap_order_* are
 *   called by functions like heap_order_double, heap_order_int, etc.  These
 *   functions make a copy of the array (the heap), an index, and sort both.
 * - insert_elt_with_index_*, delete_min_with_index_*,
 *   heap_order_with_index_*_internal, and heap_order_with_index_* are called by
 *   functions like heap_order_with_index_double, heap_order_with_index_size_t,
 *   etc.  These functions use an incoming index, which may select a portion of
 *   the original array.  They create a copy of the selected (and maybe already
 *   ordered portion of the array), an index, and then order.  The changed order
 *   is transferred to the original index array.
 * - insert_*_in_situ, delete_min_*_in_situ, heap_sort_*_in_situ are called by
 *   heap_sort_double_in_situ.  These functions sort the incoming array
 *   directly without making and sorting a copy or index.
 */

/* Ugh, this is where C++ templating could help! */

/**
 * Insert element into heap.
 *
 * @param heap heap
 * @param heap_size length of heap
 * @param element selected element of vector to sort
 * @param index indices of vector elements in heap
 * @param iter element of index
 */
#define MAKE_INSERT_ELT(TYPE, rev, CMP) \
void insert_elt_ ## TYPE ## rev (TYPE *heap, size_t *heap_size,\
	TYPE element, size_t *index, size_t iter)\
{\
	size_t now;\
\
	/* insert element at iter into last position of heap */\
	heap[*heap_size] = element;\
	index[*heap_size] = iter;\
	now = ++(*heap_size);\
\
        /* adjust its position ?? */\
        while (((now >> 1) >= 1) &&\
		(CMP(heap[(now >> 1) - 1], element) > 0)) {\
		heap[now - 1] = heap[(now >> 1) - 1];\
		index[now - 1] = index[(now >> 1) - 1];\
                now = now >> 1;\
        }\
	index[now-1] = iter;\
	heap[--now] = element;\
} /* MAKE_INSERT_ELT() */

MAKE_INSERT_ELT(double, , COMPARE)
MAKE_INSERT_ELT(int, , COMPARE)
MAKE_INSERT_ELT(size_t, , COMPARE)
MAKE_INSERT_ELT(int, rev, REV_COMPARE)
MAKE_INSERT_ELT(size_t, rev, REV_COMPARE)

/**
 * Part of heap sort: delete minimum element and replace next smallest element
 *
 * - heap[1] is the minimum element.  Remove heap[1].  Size of the heap is
 *   decreased.
 * - Put last element in heap[1] and see if it fits.
 * - If it does not fit, take minimum element among both its children and
 *   replace parent with it.
 * - Again see if the last element fits in that place.
 * - This function returns the index of the first element of the heap: this is
 *   the index of the next smallest element.
 *
 * @param heap heap
 * @param index
 * @param heap_size current size of heap
 * @return next smallest element
 */
#define MAKE_DELETE_MIN(TYPE, rev, CMP) \
size_t delete_min_ ## TYPE ## rev (TYPE *heap, size_t *index, size_t *heap_size)\
{\
\
	TYPE last_element;\
	size_t child, now, lastindex, index_min = index[0];\
\
	/* get last element and its index */\
	(*heap_size)--;\
	last_element = heap[*heap_size];\
	lastindex = index[*heap_size];\
\
        /* now refers to the index at which we are now */\
        for (now = 1; now*2 <= *heap_size ; now = child) {\
        	/* child is the index of the element which is minimum among\
		 * both the children with indexes i*2 and i*2 + 1*/\
                child = now*2;\
\
		/* child!=heap_size beacuse heap[heap_size+1] does not exist,\
		 * which means it has only one child */\
		if (child != (*heap_size)\
			&& CMP(heap[child], heap[child - 1]) < 0)\
			child++;\
\
		/* To check if the last element fits or not it suffices to\
		 * check if the last element is less than the minimum element\
		 * among both the children */\
		if (CMP(last_element, heap[child - 1]) > 0) {\
			heap[now - 1] = heap[child - 1];\
			index[now - 1] = index[child - 1];\
		} else /* it fits there */\
			break;\
	}\
	heap[now - 1] = last_element;\
	index[now - 1] = lastindex;\
        return index_min;\
} /* MAKE_DELETE_MIN */

MAKE_DELETE_MIN(double, , COMPARE)
MAKE_DELETE_MIN(int, , COMPARE)
MAKE_DELETE_MIN(size_t, , COMPARE)
MAKE_DELETE_MIN(int, rev, REV_COMPARE)
MAKE_DELETE_MIN(size_t, rev, REV_COMPARE)

/**
 * Heap sort of a vector.  Actually this is a macro that defines a function
 * to perform heap sorts on different types of vectors.
 *
 * @param array the vector to sort
 * @param size the length of the vector
 * @param sorted_index indices of vector; will be in sort order
 * @param first_k provide indices of first k smallest values
 * @return error status
 */
#define MAKE_HEAP_ORDER(TYPE, rev) \
size_t *heap_order_ ## TYPE ## rev (TYPE *array, size_t size, size_t first_k)\
{\
	TYPE *heap = malloc(size * sizeof *heap);\
	size_t *index = malloc(size * sizeof *index);\
	size_t iter = 0, heap_size = 1;\
	size_t last_index =  MIN(first_k, size >> 1);\
\
	if (!heap)\
		return NULL;\
	if (!index) {\
		free(heap);\
		return NULL;\
	}\
\
	/* initialize heap */\
	heap[0] = array[iter];\
	index[0] = iter;\
\
        for (iter = 1; iter < size; iter++)\
		insert_elt_ ## TYPE ## rev(heap, &heap_size, array[iter],\
			index, iter);\
\
	for(iter = 0; iter < first_k; iter++)\
		index[size - iter - 1] = delete_min_ ## TYPE ## rev(heap,\
			index, &heap_size);\
\
	for (iter = 0; iter < last_index; iter++)\
		SWAP_size_t(index[iter], index[size - iter - 1]);\
\
	free(heap);\
	return index;\
} /* MAKE_HEAP_ORDER() */

MAKE_HEAP_ORDER(double, )
MAKE_HEAP_ORDER(int, )
MAKE_HEAP_ORDER(size_t, )
MAKE_HEAP_ORDER(int, rev)
MAKE_HEAP_ORDER(size_t, rev)

/**
 * Heap sort of a column in a matrix.  Actually this is a macro that defines a
 * function to perform heap sorts on different types of vectors by one of their
 * coordinates.
 *
 * @param array the vector to sort
 * @param size the length of the vector
 * @param sorted_index indices of vector; will be in sort order
 * @param first_k provide indices of first k smallest values
 * @param coord coordinate to sort on
 * @return error status
 */
#define MAKE_HEAP_ORDER_BY_COORD(TYPE, rev) \
size_t *heap_order_by_coord_ ## TYPE ## rev (TYPE **array, size_t size,\
	size_t first_k, size_t coord)\
{\
	TYPE *heap = malloc(size * sizeof *heap);\
	size_t *index = malloc(size * sizeof *index);\
	size_t iter = 0, heap_size = 1;\
	size_t last_index =  MIN(first_k, size >> 1);\
\
	if (!heap)\
		return NULL;\
	if (!index) {\
		free(heap);\
		return NULL;\
	}\
\
	/* initialize heap */\
	heap[0] = array[iter][coord];\
	index[0] = iter;\
\
        for (iter = 1; iter < size; iter++)\
		insert_elt_ ## TYPE ## rev(heap, &heap_size,\
			array[iter][coord], index, iter);\
\
	for(iter = 0; iter < first_k; iter++)\
		index[size - iter - 1] = delete_min_ ## TYPE ## rev(heap,\
			index, &heap_size);\
\
	for (iter = 0; iter < last_index; iter++)\
		SWAP_size_t(index[iter], index[size - iter - 1]);\
\
	free(heap);\
	return index;\
} /* MAKE_HEAP_ORDER_BY_COORD() */

MAKE_HEAP_ORDER_BY_COORD(double, )


/**
 * Heap sort of a double vector given an index vector.  The index vector may
 * already partially sort or subselect a larger array.  If so, heap sort will
 * start with the input order and sort only the subselected entries of the
 * input vector.
 *
 * @param array the vector to sort
 * @param size the length of the vector
 * @param first_k provide indices of first k smallest values
 * @param caller_index indices of vector; will be ordered to sort array
 * @return error status
 */
#define MAKE_HEAP_ORDER_WITH_INDEX(TYPE) \
int heap_order_with_index_ ## TYPE(TYPE *array, size_t size, size_t first_k,\
	size_t *caller_index)\
{\
	TYPE *heap = NULL;		/* local, re-ordered copy of array */\
	size_t *index = NULL;		/* local, re-ordered index of heap */\
	size_t *new_index = NULL;\
	size_t iter = 0, heap_size = 1;\
	size_t last_index =  MIN(first_k, size >> 1);\
\
	heap = malloc(size * sizeof *heap);\
	index = malloc(size * sizeof *index);\
\
	if (!heap || !index) {\
		if (heap)\
			free(heap);\
		return 1;\
	}\
\
	/* initialize heap */\
	heap[0] = array[caller_index[0]];	/* extract elts from array */\
	index[0] = iter;\
/*fprintf(stderr, "caller_index[%zu]: %zu\n", (size_t)0, caller_index[0]);*/\
        for (iter = 1; iter < size; iter++)\
{\
/*fprintf(stderr, "caller_index[%zu]: %zu\n", iter, caller_index[iter]);*/\
		insert_elt_ ## TYPE(heap, &heap_size,\
			array[caller_index[iter]], index, iter);\
}\
\
	for(iter = 0; iter < first_k; iter++)\
		index[size - iter - 1] = delete_min_ ## TYPE(heap, index,\
			&heap_size);\
\
	for (iter = 0; iter < last_index; iter++)\
		SWAP_size_t(index[iter], index[size - iter - 1]);\
\
	free(heap);\
\
	/* retrieve and store the new ordering in caller_index */\
	new_index = malloc(size * sizeof *new_index);\
	if (!new_index) {\
		free(index);\
		return 1;\
	}\
\
	for (iter = 0; iter < size; iter++)\
{\
/*fprintf(stderr, "index[%zu]: %zu %zu\n", iter, index[iter], caller_index[index[iter]]);*/\
		new_index[iter] = caller_index[index[iter]];\
}\
\
	memcpy(caller_index, new_index, size * sizeof *new_index);\
	free(new_index);\
\
	return 0;\
} /* MAKE_HEAP_ORDER_WITH_INDEX() */

MAKE_HEAP_ORDER_WITH_INDEX(double)
MAKE_HEAP_ORDER_WITH_INDEX(int)
MAKE_HEAP_ORDER_WITH_INDEX(size_t)

/**
 * See comments for non-in-situ version.
 */
#define MAKE_INSERT_IN_SITU(TYPE) \
void insert_ ## TYPE ## _in_situ(TYPE *heap, size_t *heap_size, TYPE element)\
{\
        size_t now;\
        heap[*heap_size] = element;\
	now = ++(*heap_size);\
        while (((now >> 1) >= 1) && (heap[(now >> 1) - 1] > element)) {\
		heap[now - 1] = heap[(now >> 1) - 1];\
                now = now >> 1;\
        }\
        heap[now - 1] = element;\
} /* MAKE_INSERT_IN_SITU */

MAKE_INSERT_IN_SITU(double)
MAKE_INSERT_IN_SITU(int)
MAKE_INSERT_IN_SITU(size_t)

#define MAKE_DELETE_MIN_IN_SITU(TYPE) \
TYPE delete_min_ ## TYPE ## _in_situ(TYPE *heap, size_t *heap_size)\
{\
	TYPE last_element, min_element = heap[0];\
	size_t child, now;\
\
	last_element = heap[--(*heap_size)];\
\
	for (now = 1; now*2 <= *heap_size ; now = child) {\
		child = now*2;\
		if (child != (*heap_size) && heap[child] < heap[child - 1])\
			child++;\
		if (last_element > heap[child - 1])\
			heap[now - 1] = heap[child - 1];\
		else\
			break;\
	}\
	heap[now - 1] = last_element;\
	return min_element;\
} /* MAKE_DELETE_MIN_IN_SITU() */

MAKE_DELETE_MIN_IN_SITU(double)
MAKE_DELETE_MIN_IN_SITU(int)
MAKE_DELETE_MIN_IN_SITU(size_t)

#define MAKE_HEAP_SORT_IN_SITU(TYPE) \
void heap_sort_ ## TYPE ## _in_situ(TYPE *array, size_t size, size_t first_k) \
{\
	size_t iter, heap_size = 1;\
	size_t last_index = MIN(first_k, size >> 1);\
\
	for (iter = 1; iter < size; iter++)\
		insert_ ## TYPE ## _in_situ(array, &heap_size, array[iter]);\
        for (iter = 0; iter < first_k; iter++)\
		array[size - iter - 1] = delete_min_ ## TYPE ## _in_situ(array,\
			&heap_size);\
	for (iter = 0; iter < last_index; iter++)\
		SWAP_size_t(array[iter], array[size - iter - 1]);\
}/* MAKE_HEAP_SORT_IN_SITU() */

MAKE_HEAP_SORT_IN_SITU(double)
MAKE_HEAP_SORT_IN_SITU(int)
MAKE_HEAP_SORT_IN_SITU(size_t)

/*** END: heap sort algorithms ***/


/*** BEGIN: our second implementation of quicksort and quickselect ***/

/* sorting that uses less memory, but requires C99 or greater standard */
#if (defined C99 || defined C11)
size_t index_partition(void *, size_t *, CompareVectorElts, size_t, size_t, size_t, va_list);
void vindex_quicksort(void *, size_t *, CompareVectorElts, size_t, size_t, va_list);
void vindex_quickselect(void *, size_t *, CompareVectorElts, size_t, size_t, size_t, va_list);

/**
 * A #index_quicksort() wrapper for index sorting of integer arrays.
 *
 * This function allocates an index array, and then uses the quick sort
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec integer vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_int_quicksort(int *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_int_elts, 0, n-1);
} /* order_int_quicksort */

/**
 * A #index_quicksort() wrapper for index sorting of size_t positive integer
 * arrays.
 *
 * This function allocates an index array, and then uses the quick sort
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec size_t vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_size_t_quicksort(size_t *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_size_t_elts, 0, n-1);
} /* order_size_t_quicksort */

/**
 * A #index_quicksort() wrapper for index sorting of double arrays.
 *
 * This function allocates an index array, and then uses the quick sort
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort() directly with the appropriate comparisonFunc.
 *
 * @param vec double vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_double_quicksort(double *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_double_elts, 0, n-1);
} /* order_double_quicksort */

void
order_double_quicksort_with_index(double *vec, size_t *idx, size_t n)
{
	with_index_quicksort((void *)vec, idx, compare_double_elts, 0, n-1);
} /* order_double_quicksort_with_index */

void
order_size_t_quicksort_with_index(size_t *vec, size_t *idx, size_t n)
{
	with_index_quicksort((void *)vec, idx, compare_size_t_elts, 0, n-1);
} /* order_size_t_quicksort_with_index */

void
order_uint_quicksort_with_index(unsigned int *vec, size_t *idx, size_t n)
{
	with_index_quicksort((void *)vec, idx, compare_uint_elts, 0, n-1);
} /* order_uint_quicksort_with_index */

void
reverse_order_uint_quicksort_with_index(unsigned int *vec, size_t *idx, size_t n)
{
	with_index_quicksort((void *)vec, idx, reverse_compare_uint_elts, 0, n-1);
} /* reverse_order_uint_quicksort_with_index */

void
reverse_order_size_t_quicksort(size_t *vec, size_t *idx, size_t n)
{
	with_index_quicksort((void *)vec, idx, reverse_compare_size_t_elts, 0, n - 1);
} /* reverse_order_size_t_quicksort_with_index */

/**
 * A #with_index_quickselect() wrapper for index sorting of integer arrays.
 *
 * This function uses the quick select algorithm to partially sort the index in
 * increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec integer vector to index sort
 * @param n length of vector
 * @param k sort up to kth smallest element
 * @return indices partially sorted in order of increasing array value
 */
void
select_int_with_index(int *vec, size_t* idx, size_t n, size_t k)
{
	with_index_quickselect((void *)vec, idx, compare_int_elts, 0, n-1, k);
} /* select_int_with_index */

/**
 * A #with_index_quickselect() wrapper for index sorting of double arrays.
 *
 * This function uses the quick select algorithm to partially sort the index in
 * increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec double vector to index sort
 * @param n length of vector
 * @param k sort up to kth smallest element
 * @return indices partially sorted in order of increasing array value
 */
void
select_double_with_index(double *vec, size_t* idx, size_t n, size_t k)
{
	with_index_quickselect((void *)vec, idx, compare_double_elts, 0, n-1, k);
} /* select_double_with_index */

void
reverse_select_double_with_index(double *vec, size_t* idx, size_t n, size_t k)
{
	with_index_quickselect((void *)vec, idx, reverse_compare_double_elts, 0, n-1, k);
} /* select_double_with_index */

/**
 * Quick sort by index for arrays of simple data types.
 * The permuted index it returns would put the array in first argument vec into
 * ascending order.  Order is determined by the user-provided last argument
 * compare.
 *
 * @param vec array of simple data types to order
 * @param left left most index (inclusive)
 * @param right right most index (inclusive)
 * @param compare #CompareVectorElts for comparing elements in array vec
 * @return ordered index array
 */
size_t *
index_quicksort(void *vec, CompareVectorElts compare,
	size_t left, size_t right, ...)
{
	size_t i, last;
	va_list vl;
	size_t *index;

	/* we will actuall only sort indices */
	last = right - left + 1;
	MAKE_VECTOR(index, last);

	if (!index)
		return NULL;

	for (i = 0; i < last; i++)
		index[i] = i;

	va_start(vl, right);

	vindex_quicksort(vec, index, compare, left, right, vl);

	va_end(vl);

	return index;
} /* index_quicksort */

/**
 * Quick sort with an index provided by the user (avoids
 * allocation/deallocation over multiple calls).
 *
 * @param vec array of simple data types to order
 * @param left leftmost index (inclusive)
 * @param right rightmost index (inclusive)
 * @param compare #CompareVectorElts for comparing elements in array vec
 */
void
with_index_quicksort(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, ...)
{
	va_list vl;
	va_start(vl, right);
	vindex_quicksort(vec, index, compare, left, right, vl);
	va_end(vl);
} /* with_index_quicksort */

/**
 * Quick select with index array provided by user.
 *
 * @param vec array of simple data types to partially order
 * @param left leftmost index (inclusive)
 * @param right rightmost index (inclusive)
 * @param compare #CompareVectorElts for comparing elements in array vec
 */
void
with_index_quickselect(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, size_t k, ...)
{
	va_list vl;
	va_start(vl, k);
	vindex_quickselect(vec, index, compare, left, right, k, vl);
	va_end(vl);
} /* with_index_quickselect */

/**
 * va_list version of index_quicksort
 */
void
vindex_quicksort(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, va_list vl)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;
	va_list vl2;

	if (left < right) {
		/* choose median of first, middle, last elements (Wikipedia) */
		middle = (right - left)/2 + left;
		va_copy(vl2, vl);
		compare1 = compare(vec, index, left, middle, vl2);	/* !!! */
		va_end(vl2);
		va_copy(vl2, vl);
		compare2 = compare(vec, index, middle, right, vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		compare3 = compare(vec, index, left, right, vl2);
		va_end(vl2);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		va_copy(vl2, vl);
		pivot = index_partition(vec, index, compare, left, right, pivot,
			vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		if (pivot) vindex_quicksort(vec, index, compare, left, pivot - 1, vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		vindex_quicksort(vec, index, compare, pivot + 1, right, vl2);
		va_end(vl2);
	}
} /* vindex_quicksort */

/**
 * Quick select.
 *
 * @param vec array of elements
 * @param index indices of elements
 * @param compare function used to compare elements
 * @param left leftmost index inclusive
 * @param right rightmost index inclusive
 * @param k seek kth smallest element
 */
void
vindex_quickselect(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, size_t k, va_list vl)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;
	va_list vl2;

	if (left < right) {
		/* choose median of first, middle, last elements (Wikipedia) */
		middle = (right - left)/2 + left;
		va_copy(vl2, vl);
		compare1 = compare(vec, index, left, middle, vl2);	/* !!! */
		va_end(vl2);
		va_copy(vl2, vl);
		compare2 = compare(vec, index, middle, right, vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		compare3 = compare(vec, index, left, right, vl2);
		va_end(vl2);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		va_copy(vl2, vl);
		pivot = index_partition(vec, index, compare, left, right, pivot,
			vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		if (k < pivot) vindex_quickselect(vec, index, compare, left, pivot - 1, k, vl2);
		else vindex_quickselect(vec, index, compare, pivot + 1, right, k, vl2);
		va_end(vl2);
	}
} /* vindex_quickselect */

/**
 * Quick sort partition algorithm to work with #index_quicksort().  See
 * <a href="http://en.wikipedia.org/wiki/Quicksort">Wikipedia</a> for more
 * information.
 *
 * @param vec integer vector whose ordering is sought
 * @param index order vector
 * @param left left index (inclusive)
 * @param right right index (inclusive)
 * @param pivot input pivot index
 * @param compare #ComparisonFunc
 * @return new pivot
 */
size_t
index_partition(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, size_t pivot, va_list vl)
{
//const double **dvec = (const double **)vec;
	size_t i, new_pivot, swapper;
	va_list vl2;

	/* move pivot to end of array */
	swapper = index[right];
	index[right] = index[pivot];
	index[pivot] = swapper;

	/* sort this partition */
	new_pivot = left;
	for (i = left; i < right; i++) {
		va_copy(vl2, vl);
		if (compare(vec, index, i, right, vl2) < 0) {
			swapper = index[new_pivot];
			index[new_pivot] = index[i];
			index[i] = swapper;
			new_pivot++;
		}
		va_end(vl2);
	}

	/* move pivot to appropriate position */
	swapper = index[new_pivot];
	index[new_pivot] = index[right];
	index[right] = swapper;

	return new_pivot;
} /* index_partition */

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an integer array.
 * @param v vector of integers
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_int_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const int *iv = (int *)v;
	RETURN_CMP(iv[index[left]], iv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an size_t array.
 * @param v vector of size_t
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_size_t_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const size_t *iv = (size_t *)v;
	RETURN_CMP(iv[index[left]], iv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an unsigned int array.
 * @param v vector of unsigned ints
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_uint_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const unsigned int *dv = (unsigned int *)v;
	RETURN_CMP(dv[index[left]], dv[index[right]]);
}


/**
 * A #CompareVectorElts-compliant function for comparing two elements in a double array.
 * @param v vector of doubles
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_double_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const double *dv = (double *)v;
	RETURN_CMP(dv[index[left]], dv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in a unsigned long array.
 * @param v vector of unsigned long
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_ulong_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const unsigned long *ulv = (unsigned long *)v;
	RETURN_CMP(ulv[index[left]], ulv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in a string array.
 * @param v vector of const strings
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_string_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const char **sv = (const char **)v;
	return strcmp(sv[index[left]], sv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an size_t array.
 * @param v vector of size_t
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int reverse_compare_size_t_elts(const void *v, size_t *index, size_t left,
	size_t right, va_list vl)
{
	UNUSED(vl);
	const size_t *iv = (size_t *)v;
	REV_RETURN_CMP(iv[index[left]], iv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an unsigned int array.
 * @param v vector of unsigned ints
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int reverse_compare_uint_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const unsigned int *dv = (unsigned int *)v;
	REV_RETURN_CMP(dv[index[left]], dv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in a double array.
 * @param v vector of doubles
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int reverse_compare_double_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	UNUSED(vl);
	const double *dv = (double *)v;
	REV_RETURN_CMP(dv[index[left]], dv[index[right]]);
}

/*** END: our second implementation of quicksort and quickselect ***/



/*** failed attempt at median algorithm ***/

SIZE_T select_index(double *x, SIZE_T left, SIZE_T right, SIZE_T k);
SIZE_T partition(double *x, SIZE_T left, SIZE_T right, SIZE_T idx);

SIZE_T
median_of_medians(double *x, SIZE_T left, SIZE_T right)
{
	SIZE_T i;
	SIZE_T k = 2;
	SIZE_T nmedians = (right - left) / 5;
	SIZE_T sleft, sright;
	SIZE_T med_idx;
	double d;

/*
fprintf(stderr, "median(%u, %u)\n", (unsigned)left, (unsigned)right);
for (i = left; i <= right; i++)
fprintf(stderr, ", %f", x[i]);
fprintf(stderr, "\n");
*/
	for (i = 0; i <= nmedians; i++) {
		sleft = left + i*5;
		sright = sleft + 5 - 1;
		if (sright > right) {
			sright = right;
			k = (sright - sleft) / 2;
		}
		med_idx = select_index(x, sleft, sright, k);
//fprintf(stderr, "Selected %f as possible median.\n", x[med_idx]);
		d = x[left + i];
		x[left + i] = x[med_idx];
		x[med_idx] = d;
	}
	return select_index(x, left, left + nmedians, nmedians/2);
} /* median_of_medians */

SIZE_T
select_index(double *vec, SIZE_T left, SIZE_T right, SIZE_T k)
{
	SIZE_T pvt_idx, new_pvt_idx, pd;

/*
fprintf(stderr, "select_index(%u, %u, %u): ", (unsigned)left, (unsigned)right, (unsigned)k);
for (SIZE_T i = left; i <= right; i++)
fprintf(stderr, ", %f", vec[i]);
fprintf(stderr, "\n");
*/

	if (right - left <= 9) {
		SIZE_T *order = order_double_quicksort(&vec[left], right - left + 1);
		pd = left + order[k];
		free(order);
//fprintf(stderr, "Returning %u\n", (unsigned)pd);
		return pd;
	}

	do {
		pvt_idx = median_of_medians(vec, left, right);

		new_pvt_idx = partition(vec, left, right, pvt_idx);
//fprintf(stderr, "new_pvt_idx = %u (%f)\n", (unsigned) new_pvt_idx, vec[new_pvt_idx]);
		pd = new_pvt_idx - left; /* + 1; */
		if (pd == k) {
//fprintf(stderr, "select_index() returning: %f\n", vec[new_pvt_idx]);
			return new_pvt_idx;
}
		else if (k < pd)
			right = new_pvt_idx - 1;
		else {
			k = k - pd - 1;/* + 1; KSD*/
			left = new_pvt_idx + 1;
		}
//fprintf(stderr, "%u %u %u\n", (unsigned)left, (unsigned)right, (unsigned)k);
	} while (1);
} /* select_index */

SIZE_T
partition(double *x, SIZE_T left, SIZE_T right, SIZE_T idx)
{
	SIZE_T sidx, i;
	double v = x[idx], d;

//fprintf(stderr, "partition(%u, %u, %u) on %f\n", (unsigned)left, (unsigned)right, (unsigned)idx, v);

	x[idx] = x[right];
	x[right] = v;
	sidx = left;
	for (i = left; i < right; i++) {
		if (x[i] <= v) {
			d = x[i];
			x[i] = x[sidx];
			x[sidx] = d;
			sidx++;
		}
	}
	d = x[sidx];
	x[sidx] = x[right];
	x[right] = d;
/*
for (i = left; i <= right; i++)
fprintf(stderr, " %f", x[i]);
fprintf(stderr, "\n");
*/
	return sidx;
} /* partition */

#endif
