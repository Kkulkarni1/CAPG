#include <stdlib.h>
#include <math.h>

#include "math.h"
#include "brent.h"
#include "error.h"
#include "order.h"

/**
 * Compute squared Euclidean distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return squared distance between two vectors
 */
double sq_euc_dis(double *x, double *y, size_t p)
{
	double sum = 0;
	size_t i;

	 for (i = 0; i < p; i++)
		  sum += SQ(x[i] - y[i]);
	 return sum;
} /* sq_euc_dis */

/**
 * Compute Euclidean distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return distance between two vectors
 */
double euc_dis(double *x, double *y, size_t p)
{
	return sqrt(sq_euc_dis(x, y, p));
} /* euc_dis */

/**
 * Compute Frobenius norm of difference between two vectorized matrices.
 *
 * @param x	vectorized square matrix of length p*p
 * @param y	vectorized square matrix of length p*p
 * @param p	dimension of one side
 * @param diag	include diagonal in distance
 * @return	distance between two vectors
 */
double frobenius_norm(double *x, double *y, size_t p, int diag)
{
	double fn = 0;
	for (size_t i = 0; i < p; ++i)
		for (size_t j = 0; j < p; ++j)
			if (diag || i != j)
				fn += SQ(x[i*p + j] - y[i*p + j]);

	return sqrt(fn);
} /* frobenius_norm */

/**
 * Compute Hamming distance between two p-vectors.
 * [KSD] duplicate of kmodes.c:hd()
 *
 * @param x	character vector of length p
 * @param y	character vector of length p
 * @param p	length of vector
 * @return	Hamming distance between two vectors
 */
size_t hamming_char_dis(char *x, char *y, size_t p)
{
	size_t hd = 0;
	for (size_t i = 0; i < p; ++i)
		hd += x[i] != y[i];
	return hd;
} /* hamming_char_dis */

/**
 * Determine if given element is in a vector.
 * 
 * @param array	pointer to vector of numbers
 * @param num	element to find in vector
 * @param l	length of vector
 * @return	1/0 for true/false
 **/
int find_uint(unsigned int *array, unsigned int num, size_t l)
{
	for (size_t i = 0 ; i < l; i++)
		if (num == array[i])
			return 1;
	return 0;
}/* find_uint */


/**
 * Find medium from sorted vector x[start:end_inclusive].
 *
 * @param x		vector of numeric values
 * @param start		first index to include in dataset
 * @param end_inclusive	last index to include in dataset
 * @return		median
 */
double median(double *x, int start, int end_inclusive)
{
	int size = end_inclusive - start + 1;
	int m;

	if (size <= 0)
		return NAN;

	m = start + size / 2;
	if (size % 2)
		return x[m];
	return (x[m - 1] + x[m]) / 2.0;
}/* median */


/**
 * Find mimimum, lower quartile, median, upper quartile and maximum of
 * a set x_len of numeric values.
 *
 * @param[in]	x	vector of numeric values
 * @param[out]	result	length-5 vector of result statistics
 * @param[in]	x_len	length of vector x
 * @return		error status
 */
inline int five_number_summary(double *x, double *result, int x_len)
{
	int i, m, lower_end;

	for (i = 0; i < x_len; i++)
		if (x[i] != x[i])
			return 1;

	qsort(x, x_len, sizeof(double), double_compare);
	result[0] = x[0];
	result[2] = median(x, 0, x_len - 1);
	result[4] = x[x_len - 1];
	m = x_len / 2;
	lower_end = (x_len % 2) ? m : m - 1;
	result[1] = median(x, 0, lower_end);
	result[3] = median(x, m, x_len - 1);

	return NO_ERROR;
}/* five_number_summary */


/**
 * Brent's method to find root of a univariate function.
 *
 * @param a	lower limit
 * @param b	upper limit
 * @param f	function to zero
 * @param param	other arguments to f()
 * @return	max_{a<=x<=b} f(x)
 */
double brent(double a, double b, double (f)(double, void *), void *param)
{
	int status = 0;
	double value = a;

	do {
		value = f(value, param);
		value = local_min_rc(&a, &b, &status, value);
	} while (status);

	return value;
} /* brent */

/**
 * Brent's method for fast univariate optimization.
 *
 * @param a	capture 0
 * @param b	capture 0
 * @param eps	absolute tolerance
 * @param f	function to optimize
 * @param param	parameters needed by function f()
 * @return	optimization value
 */
double brent_wikipedia(double a, double b, double eps,
			double (f)(double, void *), void *param)
{
	double delta = eps;
	double fa = f(a, param);
	double fb = f(b, param);
	double c, s, d, fc, fs, fsb, fbc, fcd, tmp;
	char mflag = 1;

	if (fa * fb >= 0)
		return NAN;

	if (fabs(fa) < fabs(fb)) {
		tmp = a;
		a = b;
		b = tmp;
	}

	c = a;
	s = 0;
	d = 0;
	fc = fa;
	fs = fb;

	while (fs > 0 && fabs(b - a) > eps) {

		if (fa != fc && fb != fc)
			s = a*fb*fc/((fa-fb) * (fa - fc)) + b*fa*fc/((fb - fa)*(fb - fc))
				+ c*fa*fb/((fc - fa) * (fc - fb));
		else
			s = b - fb * (b - a) / (fb - fa);

		fsb = fabs(s - b);
		fbc = fabs(b - c);
		fcd = fabs(c - d);

		if (s < (3*a + b)/4 || s > b || (mflag && fsb >= fbc/2)
			|| (!mflag && fsb >= fcd/2)
			|| (mflag && fbc < delta)
			|| (!mflag && fcd < delta)) {
			s = (a + b) / 2;
			mflag = 1;
		} else {
			mflag = 0;
		}

		fs = f(s, param);
		d = c;
		c = b;

		if (fa * fb < 0)
			b = s;
		else
			a = s;

		/* swap */
		if (fabs(fa) < fabs(fb)) {
			tmp = a;
			a = b;
			b = a;
		}
	}
	return s;
}/* brent_wikipedia */
