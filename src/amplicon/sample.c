#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sample.h"
#include "error.h"

void sample(unsigned int N, unsigned int n, unsigned int *idx)
{
	unsigned int t = 0, m = 0;
	double u;

	while (m < n) {
		u = (double) rand() / RAND_MAX;

		if ( (N - t)*u >= n - m )
			++t;
		else
			idx[m++] = t++;
	}
} /* sample */

/**
 * Sample n from N without replacement.  See Vitter, J. S. (1987) ``An Efficient
 * Algorithm for Sequential Random Sampling'' ACM Transactions on Mathematical
 * Software. 13(1): 58--67.
 *
 */
void sample_vitter(unsigned int N, unsigned int n, unsigned int *idx)
{
	if (n >= N)
		return;

	double dn = (double) n;
	double dN = (double) N;
	double ninv = 1.0 / n;
	double vprime = exp(log((double) rand() / RAND_MAX) * ninv);
	unsigned int qu1 = 1 + N - n;
	double dqu1 = (double) qu1;
	int nainv = -13.;	/* recommended by Vitter */
	unsigned int thres = -nainv * n;
	unsigned int S=0, limit;
	double nmin1inv, U, X, dnS, y1, y2, top, bottom, V, q;
	unsigned int k = 0;

	while (n > 1 && thres < N) {
		nmin1inv = 1./(-1. + dn);
		do {
			do {
				X = dN * (-vprime + 1.);
				fprintf(stderr, "X = %f\n", X);
				S = (int) X;
				if (S < qu1)
					break;
				vprime = exp(log((double) rand() / RAND_MAX) * ninv);
			} while (1);
			U = (double) rand() / RAND_MAX;
			dnS = (double) -S;
			y1 = exp(log(U * dN/dqu1) * nmin1inv);
			vprime = y1 * (-X/dN + 1.) * (dqu1 / (dnS + dqu1));
			if (vprime <= 1.)
				break;
			y2 = 1.;
			top = -1. + dN;
			if (S < n - 1) {
				bottom = -dn + dN;
				limit = N - S;
			} else {
				bottom = dN + dnS - 1.;
				limit = qu1;
			}
			for (unsigned int t = N - 1; t >= limit; --t) {
				y2 = (y2 * top) / bottom;
				top = top - 1.;
				bottom = bottom - 1.;
			}
			if (dN / (dN - X) >= y1 * exp(log(y2)*nmin1inv)) {
				vprime = exp(log((double) rand() / RAND_MAX) * nmin1inv);
				break;
			}
			vprime = exp(log((double) rand() / RAND_MAX) * ninv);
		} while (1);
		idx[k++] = S;
		N = (N - 1) - S;
		dN = (dN - 1.) + dnS;
		n = n - 1;
		dn = dn - 1.;
		ninv = nmin1inv;
		qu1 = -S + qu1;
		dqu1 = dnS + dqu1;
		thres = thres + nainv;
	}
	if (n > 1) {
		top = N - n;
		while (n >= 2) {
			V = (double) rand() / RAND_MAX;
			q = top / dN;
			while (q > V) {
				S = S + 1;
				top = top - 1.;
				dN = dN - 1.;
				q = (q * top) / dN;
			}
			idx[k++] = S;
			dN = dN - 1.;
			n = n - 1;
		}
	} else {
		S = (int) (N * vprime);
		idx[k++] = S;
	}
} /* sample */

/**
 * Adjust sample_vitter() for ampliclust 
 * Sample n from N without replacement.
 *
 * @return		error status
 */
int random_sample(size_t N, size_t n,  size_t *D_idx, size_t *s_idx)
{
	if (n > N)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"invalid random sample");

	size_t t = 0, m = 0;
	double u;

	while (m < n) {
		u = (double) rand() / RAND_MAX;

		if ( (N - t)*u >= n - m )
			++t;
		else
			s_idx[m++] = D_idx[t++];
	}
	return NO_ERROR;
} /* random_sample */
