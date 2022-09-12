/**
 * @file cluster.c
 * @author Karin Dorman
 *
 * General utility functions for clustered data.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "cluster.h"
#include "math.h"
#include "error.h"

static unsigned int *__pt = NULL;
static unsigned int *__pa = NULL;
static unsigned int **__pat = NULL;

static inline int allocate_statics(unsigned int k_assigned, unsigned int k_true)
{
	__pt = calloc(k_true, sizeof *__pt);
	__pa = calloc(k_assigned, sizeof *__pa);
	unsigned int *tmp = calloc(k_assigned * k_true, sizeof **__pat);
	__pat = malloc(k_assigned * sizeof *__pat);

	if (!__pt || !__pa || !tmp || !__pat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, NULL);

	for (unsigned int k = 0; k < k_assigned; ++k) {
		__pat[k] = tmp;
		tmp += k_true;
	}

	return NO_ERROR;
} /* allocate_statics */

static inline void reset_statics(unsigned int k_assigned, unsigned int k_true)
{
	memset(__pt, 0, k_true * sizeof *__pt);
	memset(__pa, 0, k_assigned * sizeof *__pa);
	memset(*__pat, 0, k_assigned * k_true * sizeof **__pat);
} /* reset_statics */

void free_cluster_statics()
{
	if (__pt) {
		free(__pt);
		__pt = NULL;
	}
	if (__pa) {
		free(__pa);
		__pa = NULL;
	}
	if (__pat) {
		free(__pat[0]);
		free(__pat);
		__pat = NULL;
	}
} /* free_cluster_statics */


/**
 * Compute (normalized) mutual information.
 *
 * @param iclass_assigned inferred class assignment in {0,...,k} (n x 1)
 * @param iclass_true true class assignment in {0,...,k} (n x 1)
 * @param n dimension of first two arguments
 * @param k_assigned number of assigned classes
 * @param k_true number of true classes
 * @param mi_variant what type of mutual information
 * @param normalization what type of normalization
 * @return (normalized) mutual information
 */
double
mutual_information(unsigned int *iclass_assigned, unsigned int *iclass_true,
	 unsigned int n, unsigned int k_assigned, unsigned int k_true,
	 int mi_variant, int normalization)
{
	unsigned int i, j;
	double ht = 0, ha = 0, hab = 0, mi = 0, idx;

	if (!__pt && allocate_statics(k_assigned, k_true))
		return 0.;
	else
		reset_statics(k_assigned, k_true);

	for (i = 0; i < n; i++) {
		__pt[iclass_true[i]]++;
		__pa[iclass_assigned[i]]++;
		__pat[iclass_assigned[i]][iclass_true[i]]++;
	}

	if (mi_variant) {
		for (j = 0; j < k_true; j++)
			ht += __pt[j] * (log(n) - log(__pt[j])) / n;

		for (j = 0; j < k_assigned; j++)
			ha += __pa[j] * (log(n) - log(__pa[j])) / n;
	}

	for (i = 0; i < k_true; i++)
		for (j = 0; j < k_assigned; j++) {
			if (__pat[j][i])
				mi += __pat[j][i] * (log(n) + log(__pat[j][i]) - log(__pt[i]) - log(__pa[j])) / n;
			if (__pat[j][i] && normalization == JOINT_ENTROPY_SCALING)
				hab += __pat[j][i] * (log(n) - log(__pat[j][i])) / n;
		}

	if (mi_variant == MUTUAL_INFORMATION)
		return (mi);
	else if (mi_variant == VARIATION_OF_INFORMATION) {
		idx = ht + ha - 2*mi;
		return (idx >= 0 ? idx : 0);
	} else if (mi_variant == NORMALIZED_VARIATION_OF_INFORMATION) {
		idx = (ht + ha - 2*mi)/(ht + ha - mi);
		return (idx >= 0 ? idx : 0);
	} else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == MIN_SCALING)
		return (mi / MIN(ht, ha));
	else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == GEOMETRIC_MEAN_SCALING)
		return (mi / sqrt(ht*ha));
	else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == MEAN_SCALING)
		return (2 * mi / (ht + ha));
	else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == MAX_SCALING)
		return (mi / MAX(ht, ha));
	else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == JOINT_ENTROPY_SCALING)
		return (mi / hab);
	else if (mi_variant == NORMALIZED_MUTUAL_INFORMATION
		&& normalization == NO_MI_SCALING)
		mmessage(ERROR_MSG, INTERNAL_ERROR, "Normalized mutual "
			"information requires scaling: see option -N.");
	else if (mi_variant == ADJUSTED_MUTUAL_INFORMATION)
		mmessage(ERROR_MSG, INTERNAL_ERROR, "Adjusted mutual "
			"information not implemented.");
	else if (mi_variant == STANDARDIZED_MUTUAL_INFORMATION)
		mmessage(ERROR_MSG, INTERNAL_ERROR, "Standarized mutual "
			"information not implemented.");

	return -1;
} /* mutual_information */

/**
 * Compute biological homogeneity index (bhi) [Datta & Datta (2006)].
 * Notice, this index is NAN if any cluster has 0 or 1 members.  Also,
 * computing this index is SLOW for large samples: o(Kn^2).
 *
 * @param iclass inferred class assignments in {0, ..., K-1} (n x 1)
 * @param iclass_true true class assignments in {0, ..., K-1} (n x 1)
 * @param n dimension of first two arguments
 * @param k_inferred number of classes
 * @return bhi
 */
double
biological_homogeneity_index(unsigned int *iclass_inferred,
	unsigned int *iclass_true, unsigned int n, unsigned int k_inferred)
{
	unsigned int k;
	unsigned int i, j, ncnt;
	unsigned int nm1 = n - 1;
	unsigned int nmatch;
	double bhi = 0;

	if (!iclass_true || !iclass_inferred)
		return -1;

	for (k = 0; k < k_inferred; k++) {
		nmatch = 0;
		ncnt = (iclass_inferred[n-1] == k);
		for (i = 0; i < nm1; i++) {
			if (iclass_inferred[i] == k)
				ncnt++;
			else
				continue;
			for (j = i + 1; j < n; j++)
				if (iclass_inferred[i] == iclass_inferred[j]
					&& iclass_true[i] == iclass_true[j])
					nmatch++;
		}
		if (ncnt > 1)
			bhi += 2.0 * nmatch / ncnt / (ncnt-1);
	}
	return(bhi / k_inferred);
} /* biological_homogeneity_index */

/**
 * Compute one of adjusted rand, rand, or e index.
 *
 * @param iclass_inferred estimated class assignments (1 x n)
 * @param iclass_true true class assignments (1 x n)
 * @param n number of observations
 * @param k_inferred number of classes in inferred solution
 * @param k_true number of classes in true solution
 * @param rand_type type of index to compute
 * @return index
 */
double
cluster_index(unsigned int *iclass_inferred, unsigned int *iclass_true,
	unsigned int n, unsigned int k_inferred, unsigned int k_true,
	int rand_type)
{
	unsigned int k, j;
	unsigned int i;
	unsigned int ninf_sq = 0, ntru_sq = 0, ninf_sq_x_ntru_sq = 0;
	double tmpa, tmpb;
	double sumsq = 0, sumprs = 0, sumnrpr = 0, sumncpr = 0;
	double rtn = NAN;

	if (!__pt && allocate_statics(k_inferred, k_true))
		return 0;
	else
		reset_statics(k_inferred, k_true);

	for (i = 0; i < n; i++) {
		__pt[iclass_true[i]]++;
		__pa[iclass_inferred[i]]++;
		__pat[iclass_inferred[i]][iclass_true[i]]++;
	}

	for (k = 0; k < k_inferred; k++) {
		tmpa = __pa[k] * __pa[k];
		ninf_sq += tmpa;
		for (j = 0; j < k_true; j++) {
			tmpb = __pt[j] * __pt[j];
			if (!k)
				ntru_sq += tmpb;
			ninf_sq_x_ntru_sq += tmpa * tmpb;
		}
	}


	if (rand_type == E_INDEX) {
		rtn = ninf_sq_x_ntru_sq / (n * (n-1) + n * n / (n-1)) - (ninf_sq + ntru_sq) / (n-1);
		rtn /= n * (n-1) / 2.0;
		return rtn;
	}

	for (k = 0; k < k_inferred; k++) {
		sumnrpr += __pa[k] * (__pa[k] - 1) / 2.0;
		for (j = 0; j < k_true; j++) {
			sumsq += __pat[k][j] * __pat[k][j];
			sumprs += __pat[k][j] * (__pat[k][j] - 1) / 2.0;
			if (!k)
				sumncpr += __pt[j] * (__pt[j] - 1) / 2.0;
		}
	}

	if (rand_type == RAND_INDEX) {
		rtn = 1.0 + (2.0 * sumsq - ninf_sq - ntru_sq) / (n * (n-1));
	} else if (rand_type == ADJUSTED_RAND_INDEX) {
		tmpa = 2.0 * sumnrpr * sumncpr / (n * (n-1));
		rtn = (sumprs - tmpa) / ((sumnrpr + sumncpr)/2.0 - tmpa);
	}

	return rtn;
} /* cluster_index */
