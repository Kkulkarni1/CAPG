/**
 * @file kmodes.c
 *
 * Implements various k-modes algorithms and initializations.
 *
 * WARNING: You must call reset_k() BEFORE calling kmodes_*() function with new
 * K and free_kmodes() BEFORE calling kmodes_*() function with new data.
 *
 * Implement Hartigan and Wong's Algorithm AS 136 published in Applied
 * Statistics (1979) 28(1):100 for optimizing the objective function of
 * k-modes for categorical data, introduced in Huang's A Fast Clustering
 * Algorithm to Cluster Very Large Categorical Data Sets in Data Mining
 * published in In Research Issues on Data Mining and Knowledge Discovery
 * (1997), 1--8.
 *
 * This code also implements Huang's algorithm and several initialization methods.
 *
 * Code based on Applied Statistics algorithms (C) Royal Statistical Society
 * 1979. Adapted for C by Ranjan Maitra, Baltimore, 07/12/02
 * identical to kmns.c except for indx which is a number here.
 * Further adjusted by K. Dorman, Ames, IA, 9/12.
 * Adapted to kmodes.c by K. Dorman, Ames, IA, 6/17.
 */

#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "kmodes.h"
#include "array.h"
#include "order.h"
#include "error.h"
#include "math.h"
#include "io.h"
#include "io_kmodes.h"
#include "sample.h"

/* data allocation */
size_t allocate_and_compute_category_counts(data_t **x, unsigned int n, unsigned int p, unsigned int K, int wgt);
size_t allocate_and_compute_nj(data_t **x, unsigned int n, unsigned int p);
int allocate_and_compute_njc(data_t **x, unsigned int n, unsigned int p, size_t ncat);
int allocate_nkjc(unsigned int K, unsigned int p, size_t ncat);
int allocate_hw_memory(unsigned int n, unsigned int p, unsigned int K);
int set_nkjc(unsigned int K, unsigned int p, size_t l);
int reset_nkjc(size_t ***nkjc, unsigned int K, unsigned int p);
void compute_nkjc(size_t ***nkjc, data_t **x, unsigned int *ic, unsigned int n, unsigned int p);

/* Hartigan and Wong algorithm */
void optra(data_t **a, unsigned int m, unsigned int n, data_t **c, unsigned int k, unsigned int *ic1,
	unsigned int *ic2, unsigned int *nc, unsigned int *live, unsigned int *cd, double *d, unsigned int *indx,
	unsigned int *nj, size_t **nt, size_t ***nkt, data_t **c2, int wgt, int use_qtran);

void qtran(data_t **a, unsigned int m, unsigned int n, data_t **c, unsigned int k, unsigned int *ic1,
	unsigned int *ic2, unsigned int *nc, unsigned int *live, double *d, unsigned int *indx,
	unsigned int *nj, size_t **nt, size_t ***nkt, data_t **c2, int wgt);

static inline double cost_of_membership(size_t **nkt, data_t *a, data_t *c, data_t *c2, unsigned int p);
static inline double weighted_cost_of_membership(size_t **nkt, data_t *a, data_t *c, data_t *c2, unsigned int p, size_t **nt, int fxn_debug);
static inline double cost_to_join(size_t **nkt, data_t *a, data_t *c, unsigned int p, double bcost);
static inline double weighted_cost_to_join(size_t **nkt, data_t *a, data_t *c, unsigned int p, double bcost, size_t **nt, int fxn_debug);
static inline void update_modes(size_t ***nkt, data_t *a, data_t **c, data_t **c2, unsigned int l1, unsigned int l2, unsigned int p, unsigned int *nj);
static inline int compute_modes(size_t ***nkjc, data_t **modes, data_t **mmodes, unsigned int K, unsigned int p);
static inline double compute_cost(size_t ***nkjc, double *cost, data_t **modes, unsigned int K, unsigned int p);
static inline double compute_weighted_cost(size_t ***nkjc, double *cost, data_t **modes, unsigned int K, unsigned int p);
double compute_criterion(data_t **x, data_t **seeds, unsigned int *ic, double *criterion, unsigned int K, unsigned int n, unsigned int p, int wght);

/* distance functions */
static inline double hd(data_t *x, data_t *y, unsigned int p, int weight);
static inline double hd_min(data_t *x, data_t *y, unsigned int p, int wgt, double dprev);
static inline double hd_coord(data_t *x, data_t *y, unsigned int j, int weight);

/* k-modes initialization routines */
int kmodes_init_random_seeds(data_t **x, unsigned int n, unsigned int p, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
int kmodes_init_h97(data_t **x, unsigned int n, unsigned int p, unsigned int K, unsigned int k1, int wgt, int rdm, data_t **seeds, unsigned int *sd_idx);
int kmodes_init_hd17(data_t **x, unsigned int n, unsigned int p, unsigned int K, unsigned int k1, int wgt, data_t **seeds, unsigned int *sd_idx);
int kmodes_init_av07(data_t **X, unsigned int n, unsigned int p, unsigned int K, unsigned int k1, int weight, data_t **seeds, unsigned int *sd_idx, int greedy);
int compare_data(data_t **x, unsigned int i, unsigned int j, unsigned int n);
int compare_data_to_seed(data_t **x, unsigned int i, data_t *seed, unsigned int p);

//#define __KMODES_DEBUGGING__
#ifdef __KMODES_DEBUGGING__

/* debugging functions */
int allocate_nkjc_debug(unsigned int K, unsigned int p, size_t ncat);
double debug_cost_to_move(data_t **x, unsigned int n, unsigned int p, unsigned int K, unsigned int *ic, unsigned int i, unsigned int l1, unsigned int l2, int wgt);
static inline double compute_cost_from_nkjc(size_t ***nkjc, unsigned int K, unsigned int p, int wgt);
void compare_nkjc(size_t ***nkjc1, size_t ***nkjc2, unsigned int K, unsigned int p);
void compute_and_compare_modes(size_t ***nkjc, data_t **modes, data_t **mmodes, unsigned int K, unsigned int p);
void summarize_clusters(char const *str, size_t ***nkjc, data_t **c, data_t **c2, unsigned int K, unsigned int p);

static size_t ***__nkjc_debug = NULL;
#endif


static unsigned int *__nj = NULL;	/* no. categories at each coordinate */
static size_t **__njc = NULL;		/* cat. counts in each coordinate */
static size_t ***__nkjc = NULL;		/* cat. counts w/in coords in clusters */
static double *__dens = NULL;
static double *__dis = NULL;
static unsigned int *__ic2 = NULL;
static unsigned int *__live = NULL;
static unsigned int *__cd = NULL;
static data_t **__c2 = NULL;		/* minor modes */





/**
 * Hartigan and Wong's algorithm.
 *
 * @param x		mXn data of m observations, n dimensions, all
 *			categorical, each 0...n_j, where n_j is number of
 *			distinct categories at site j, j=0...n-1.
 * @param n		number of observations
 * @param p		dimension of data
 * @param c		matrix of kxn modes (to be updated)
 * @param K		number of clusters
 * @param ic1		cluster assignment of each observation
 * @param nc		number of members in each cluster (to be filled)
 * @param max_iter	maximum number of iterations
 * @param csd		k-modes minimized criterion (to be filled)
 * @param ifault	error code (to be filled)
 * @param iter		number of iterations used
 * @param opt		pointer to structure with run options
 * 	weight		1/0 indicate if to use weights
 * 	init_update	update mode after each addition
 *	use_qtran	use quick-transfer step
 * @return		cost
 */
double
kmodes_hw(data_t **x, unsigned int n, unsigned int p, data_t **c,
	unsigned int K, unsigned int *ic1, unsigned int *nc,
	unsigned int max_iter, double *csd, int *ifault, unsigned int *iter,
	kmodes_options *opt)
//       	int weight, int init_update, int use_qtran)
{
	unsigned int k;		/* indices */
	unsigned int i, j, l = 0;
	unsigned int indx = 0;	/* no. observations evaluated without swap */
	double da, db, dc;	/* to calculate distances */


	*ifault = KMODES_NO_ERROR;

	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		return INFINITY;
	}

	/* pre-compute category counts */
	if (!(l = allocate_and_compute_category_counts(x, n, p, K,
						opt->weighted))) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

	/* trivial input: user requests 1 cluster */
	if (K == 1) {
		nc[0] = n;
		/* assign all observations to cluster 1 [PARALLEL] */
		for (i = 0; i < n; i++)
			ic1[i] = 0;
		/* compute centers [PARALLEL] */
		for (j = 0; j < p; j++) {
			unsigned int max = 0;
			for (l = 0; l < __nj[j]; l++) {
				__njc[j][l] = 0;
				for (i = 0; i < n; i++)
					__njc[j][(int) x[i][j]]++;
				if (max < __njc[j][l]) {
					max = __njc[j][l];
					c[0][j] = l;
				}
			}
		}
		/* compute sum-of-distances [PARALLEL] */
		csd[0] = 0;
		return compute_criterion(x, c, ic1, csd, K, n, p,
							opt->weighted);

	/* trivial input: user requests m clusters */
	} else if (K == n) {
		for (i = 0; i < n; i++) {
			nc[i] = 1;	/* 1 member per cluster */
			memcpy(c[i], x[i], p * sizeof **c);
			ic1[i] = i;
			csd[i] = 0;
		}

		return 0.;
	}

	/* now the interesting case of 1 < k < m */

	if (!__live && allocate_hw_memory(n, p, K)) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}


	for (k = 0; k < K; k++) {
		nc[k] = 0;	/* initialize counts per cluster */

		/* initialize minor mode to be NOT initial modes */
		for (j = 0; j < p; j++)
			__c2[k][j] = c[k][j] == 0 ? 1 : 0;
	}

	/* space for counts of categories per cluster per site */
	/* l = total number categories */
	if (set_nkjc(K, p, l)) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

#ifdef __KMODES_DEBUGGING__
	allocate_nkjc_debug(K, p, l);
#endif

	/* assign observations to clusters and find two closest centres,
	 * ic1[i] and ic2[i] [PARALLEL] */
	for (i = 0; i < n; i++) {

		/* initially guess first and second center */
		ic1[i] = 0;
#ifndef __KMODES_NO_QTRANS___
		__ic2[i] = 1;
#endif

		/* compute distances to first and second */
		da = hd(x[i], c[0], p, opt->weighted);
		db = hd(x[i], c[1], p, opt->weighted);

		/* swap first and second, if second closer */
		if (da > db) {
			ic1[i] = 1;
#ifndef __KMODES_NO_QTRANS__
			__ic2[i] = 0;
#endif
			dc = da;
			da = db;
			db = dc;
		}

		/* for remaining modes */
		for (k = 2; k < K; ++k) {

			/* compute distance */
			dc = 0.;
			for (j = 0; j < p; ++j) {
				dc += hd_coord(c[k], x[i], j, opt->weighted);
				if (dc > db)	/* losing already */
					break;
			}

			/* it is closer than second closest */
#ifndef __KMODES_NO_QTRANS__
			if (dc < db) {

				/* second closest so far */
				if (dc >= da) {
					db = dc;
					__ic2[i] = k;

				/* closest so far */
				} else {
					db = da;
					__ic2[i] = ic1[i];
					da = dc;
					ic1[i] = k;
				}
			}
#else
			if (dc < da) {
				db = da;
				da = dc;
				ic1[i] = k;
			}
#endif
		}

		nc[ic1[i]]++;
		for (j = 0; j < p; j++)
			__nkjc[ic1[i]][j][(int) x[i][j]]++;

		/* update modes (should be default behavior) */
		if (opt->init_update)
			for (j = 0; j < p; ++j) {
				if (__nkjc[ic1[i]][j][(int) x[i][j]]
					> __nkjc[ic1[i]][j][(int) c[ic1[i]][j]])
					c[ic1[i]][j] = x[i][j];
				else if (__nkjc[ic1[i]][j][(int) x[i][j]]
					== __nkjc[ic1[i]][j][(int) c[ic1[i]][j]]
					&& x[i][j] < c[ic1[i]][j])
					c[ic1[i]][j] = x[i][j];
			}
	} /* find closest centers */

	/* initialize pre-computed values needed by algorithm */
	for (k = 0; k < K; ++k) {

		/* error if initialization routine produces empty cluster */
		if (nc[k] == 0) {
			*ifault = KMODES_NULL_CLUSTER_ERROR;
			return INFINITY;
		}

		/* live[k] stores obs. last transferred to/from cluster k */
		if (max_iter) {
			__live[k] = n + 1;	/* all clusters start live */
			__cd[k] = n + 1;	/* all distances need computing */
		}
	}

	/* if not iterating, simply return cost after initialization */
	if (!max_iter)
		return opt->weighted
			? compute_weighted_cost(__nkjc, csd, c, K, p)
			: compute_cost(__nkjc, csd, c, K, p);
	//fprintf(stderr, "Initial cost: %f\n", compute_cost(__nkjc, csd, c, K, p));

	/* compute minor modes */
	if (compute_modes(__nkjc, opt->init_update ? NULL : c, __c2, K, p)) {
		*ifault = KMODES_INTERNAL_ERROR;
		return INFINITY;
	}

#ifdef __KMODES_DEBUGGING__
	fprintf(stderr, "in kmodes_hw()\n");
	reset_nkjc(__nkjc_debug, K, p);
	compute_nkjc(__nkjc_debug, x, ic1, n, p);
	compare_nkjc(__nkjc, __nkjc_debug, K, p);
	compute_and_compare_modes(__nkjc_debug, c, __c2, K, p);
#endif

	/* iterate optimal and quick transfer stages until max. iterations */
	for (i = 0; i < max_iter; i++) {

		/* optimal transfer stage: each point reallocated,
		 * if appropriate to cluster which best reduces csd */
		optra(x, n, p, c, K, ic1, __ic2, nc, __live, __cd, __dis, &indx,
						__nj, __njc, __nkjc, __c2,
						opt->weighted, opt->use_qtran);
		/* m observations processed and no change, we're done */
		if (indx == n)
			break;

		/* quick transfer stage: consider moving each point i to ic2[i]
		 * and keep looping until no more change.  */
#ifndef __KMODES_NO_QTRANS__
		if (opt->use_qtran)
			qtran(x, n, p, c, K, ic1, __ic2, nc, __live, __dis,
				&indx, __nj, __njc, __nkjc, __c2, opt->weighted);

		/* two clusters: no need for additional optimal transfer */
		if (opt->use_qtran && K == 2)
			break;
#endif
	}

	*iter = i;

	/* maximum iterations exceeded */
	if ((indx != n) && (K != 2))
		*ifault = KMODES_EXCEED_ITER_WARNING;

	/* error if initialization routine produces empty cluster */
	for (k = 0; k < K; ++k)
		if (nc[k] == 0) {
			*ifault = KMODES_NULL_CLUSTER_ERROR;
			break;
		}

	/* summing d[i] is not the same thing */
	return opt->weighted
		? compute_weighted_cost(__nkjc, csd, c, K, p)
		: compute_cost(__nkjc, csd, c, K, p);

} /* kmodes_hw */

int allocate_hw_memory(unsigned int n, unsigned int p, unsigned int K)
{
	/* allocate second best information */
#ifndef __KMODES_NO_QTRANS__
	__ic2 = malloc(n * sizeof *__ic2);
	if (!__ic2)
		return KMODES_MEMORY_ERROR;
#endif
	__c2 = malloc(K * sizeof *__c2);
	/* allocate memory for minimum distance, live and cost sets */
	__dis = malloc(n * sizeof *__dis);
	__live = malloc(K * sizeof *__live);
	__cd = malloc(K * sizeof *__cd);

	/* allocate c2 memory as a single block */
	data_t *tmp = malloc(K * p * sizeof **__c2);

	if (!__c2 || !tmp || !__dis || !__live || !__cd)
		return KMODES_MEMORY_ERROR;

	for (unsigned int k = 0; k < K; k++) {
		__c2[k] = tmp;
		tmp += p;
	}

	return NO_ERROR;
} /* allocate_hw_memory */


/**
 * Optimal-transfer stage of K-modes algorithm.
 *
 * This is the optimal transfer stage.  Each point is re-allocated, if
 * possible, to the cluster that will induce a maximum reduction in the total
 * within-cluster sum of distances.  If a cluster was updated in the last
 * quick-transfer stage, it belongs to the live set throughout this stage.
 * Otherwise, at each step, it is not in the live set if it has not been
 * updated in the last $m$ optimal transfer steps.
 *
 * @param a	data observations (m x n)
 * @param m	number of observations
 * @param n	number of dimensions
 * @param c	centers (k x n)
 * @param k	number of clusters
 * @param ic1	closest and (1 x n)
 * @param ic2	next closest center (1 x n)
 * @param nc	number in each cluster
 * @param live	live[k] = i+1 if obs. i transferred to/from cluster k (1 x k)
 * @param cd	indicates if dis to closest center NOT needed (1 x n)
 * @param d	criterion decrease if remove i from its current cluster (1 x n)
 * @param indx	first 0, number of consecutive observations not transferred
 * @param nj	number of categories at each site (1 x n)
 * @param nt	count of each category in full dataset (n x nj)
 * @param nkt	count of each category in each cluster (k x n x nj)
 * @param c2	minor modes of each cluster (k x n)
 * @param wgt	use weighted distances
 * @param qt	use quick-transfer step
 */
void optra(data_t **a, unsigned int m, unsigned int n, data_t **c, unsigned int k,
	unsigned int *ic1, unsigned int *ic2, unsigned int *nc, unsigned int *live,
	unsigned int *cd, double *d, unsigned int *indx, unsigned int *nj, size_t **nt,
	size_t ***nkt, data_t **c2, int wgt, int qt)
{
	unsigned int i, i1;		/* observation indices */
#ifdef __KMODES_NO_QTRANS__
	UNUSED(ic2);
#else
	unsigned int ll;
#endif
	unsigned int l, l1, l2 = 0;		/* closest indices */
	double r2, de;

#ifdef __KMODES_DEBUGGING__
	double min_diff;
	unsigned int best_l;
#endif

	/* one loop through all observations */
	for (i = 0, i1=1; i < m; ++i, ++i1) {
		(*indx)++;	/* one more observation to be processed */
		l1 = ic1[i];	/* current closest center */

		/* no transfer if observation i is only member of its cluster */
		if (nc[l1] != 1) {

#ifndef __KMODES_NO_QTRANS__
			l2 = ic2[i];	/* current 2nd closest center */
			ll = l2;	/* current contender (this is not l1) */
#endif

			/* centroid c[l1] updated since this observation was
			 * last considered: re-compute membership cost
			 * note: qtran unsyncs cd and live by keeping membership
			 * cost updated and leaving clusters live, so this
			 * cost does not need recalculation after qtran() */
			if (i < cd[l1])
				d[i] = wgt ? weighted_cost_of_membership(
					nkt[l1], a[i], c[l1], c2[l1], n, nt,
					ABSOLUTE_SILENCE)
					: cost_of_membership(nkt[l1], a[i],
					c[l1], c2[l1], n);

			/* cluster with min. join cost is likely to be
			 * current ll = l2 */
			r2 = INFINITY;
#ifndef __KMODES_NO_QTRANS__
			if (i1 < live[l1] || i1 < live[ll])
				r2 = wgt ? weighted_cost_to_join(nkt[ll], a[i],
					c[ll], n, INFINITY, nt,
					ABSOLUTE_SILENCE) : cost_to_join(
					nkt[ll], a[i], c[ll], n, INFINITY);
#endif

#ifdef __KMODES_DEBUGGING__
			min_diff = INFINITY;
			best_l = l1;
#endif

			/* set l2 to the new closest mode */
			for (l = 0; l < k; ++l) {

#ifdef __KMODES_DEBUGGING__
				if (l != l1) {
					/* compute in an independent way */
					de = debug_cost_to_move(a, m, n, k, ic1,
								 i, l1, l, wgt);
					debug_msg(1, 0, "Cost of moving member "
						"%zu from %u to %u is %.0f\n",
								i, l1, l, de);
					if ((i1 < live[l1] || i1 < live[l])
						&& (de < min_diff ||
						(de == min_diff && l < l2))) {
						min_diff = de;
						best_l = l;
					}
				}
#endif

				/* l1 and ll (orig l2) already computed */
#ifndef __KMODES_NO_QTRANS__
				if (l == l1 || l == ll)
#else
				if (l == l1)
#endif
					continue;

				/* two clusters to compare are not live for
				 * this observation */
				if (i1 >= live[l1] && i1 >= live[l])
					continue;

				/* to save time, abort calculation as soon as
				 * cost exceeds r2 */
				de = wgt ? weighted_cost_to_join(nkt[l], a[i],
					c[l], n, r2, nt, ABSOLUTE_SILENCE)
					: cost_to_join(nkt[l], a[i], c[l], n,
					r2);

				/* less costly cluster */
				if (de < r2 || (de == r2 && l < l2)) {
					/* [TODO] create l3 and maintain record of second closest */
					r2 = de;	/* cost of joining l */
					l2 = l;
				}
			}
#ifdef __KMODES_DEBUGGING__
			/* */
			if (!isinf(min_diff) && (
				min_diff != r2 - d[i] || best_l != l2)) {
				fprintf(stderr, "Cost to be in %u\n", l1);
				wgt ? weighted_cost_of_membership(nkt[l1], a[i],
					c[l1], c2[l1], n, nt, DEBUG_I)
					: cost_of_membership(nkt[l1], a[i],
					c[l1], c2[l1], n);
				fprintf(stderr, "Cost to join %u\n", l2);
				wgt ? weighted_cost_to_join(nkt[l2], a[i], c[l2],
					n, INFINITY, nt, DEBUG_I)
					: cost_to_join(nkt[l2], a[i], c[l2], n,
					INFINITY);
				fprintf(stderr, "Cost to join %u\n", best_l);
				wgt ? weighted_cost_to_join(nkt[best_l], a[i],
					c[best_l], n, INFINITY, nt, DEBUG_I)
					: cost_to_join(nkt[best_l], a[i],
					c[best_l], n, INFINITY);
				summarize_clusters("optra0:", nkt, c, c2, k, n);
				fprintf(stderr, "Observation %u: ", i);
				fprint_data_ts(stderr, a[i], n, 1, 1);
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
					"[min_diff] %.0f != %.0f - %.0f = %.0f "
					"for move %u -> %u (%u)\n", min_diff,
					r2, d[i], r2 - d[i], ic1[i], l2,
					best_l));
			}
#endif

			/* no transfer improves criterion, l2 is new ic2[i] */
			if (r2 > d[i] || (r2 == d[i] && l1 < l2))
#ifndef __KMODES_NO_QTRANS__
				ic2[i] = l2;
#endif
			/* transfer i from l1 to l2, l1 is new ic2[i] */
			else {
				*indx = 0;	/* obs.s since last transfer */
				live[l1] = m + i1;
				live[l2] = m + i1;
				cd[l1] = live[l1];	/* now cd and live */
				cd[l2] = live[l2];	/* have same info */

				/* update cluster counts & assignments */
				nc[l1]--;
				nc[l2]++;
				ic1[i] = l2;
#ifndef __KMODES_NO_QTRANS__
				ic2[i] = l1;
#endif

				update_modes(nkt, a[i], c, c2, l1, l2, n, nj);

			}


#ifdef __KMODES_DEBUGGING__
			if (__nkjc_debug) {
				reset_nkjc(__nkjc_debug, k, n);
				compute_nkjc(__nkjc_debug, a, ic1, m, n);
				compare_nkjc(__nkjc, __nkjc_debug, k, n);
				compute_and_compare_modes(__nkjc_debug, c, c2, k, n);
			}
#endif

		}

		/* no transfers: algorithm has converged */
		if (*indx == m)
			return;
	}

	/* there was >= 1 transfer; adjust last obs. transferred
	 * to/from l */
	for (l = 0; l < k; ++l) {
		live[l] -= m;
		if (qt)
			cd[l] = 0;	/* qtran keeps distances up-to-date */
		else
			cd[l] = live[l];
	}

	return;
} /* optra */


/**
 * Quick transfer stage of Hartigan and Wong's K-means algorithm.
 *
 * This is the quick transfer stage. ic1[i] is the cluster observation i
 * belongs to.  ic2[i] is the next closest cluster, the one i is most likely to
 * transfer to. For each observation i, ic1[i] & ic2[i] are switched, if
 * necessary, to reduce within-cluster sum of squares.  The cluster centres are
 * updated after each step. In the optimal transfer stage, live[l] indicates
 * which observation was last transferred in/out of cluster l.  Later
 * observations have compared against the new center; earlier observations have
 * not yet.
 *
 * @param a data observations
 * @param m number of observations
 * @param n number of dimensions
 * @param c centers
 * @param k number of clusters
 * @param ic1 closest and
 * @param ic2 next closest center
 * @param nc number in each cluster
 * @param live set to i if i last observation transferred to/from cluster
 * @param d distance to assigned centroid for each point
 * @param indx first 0, number of consecutive observations not transferred
 * @param nj number of categories at each site
 * @param nt count of each category in full dataset
 * @param nkt count of each category in each cluster
 * @param c2 second maximal category count in each cluster
 * @param wgt use weighted distances
 */
void qtran(data_t **a, unsigned int m, unsigned int n, data_t **c, unsigned int k, unsigned int *ic1,
	unsigned int *ic2, unsigned int *nc, unsigned int *live, double *d, unsigned int *indx,
	unsigned int *nj, size_t **nt, size_t ***nkt, data_t **c2, int wgt)
{
	unsigned int l1, l2, j;	/* cluster indices */
	double da;
	unsigned int mp1 = m + 1;
	unsigned int i;		/* current observation */
	unsigned int istep = 0;	/* total observations processed this call */
	unsigned int loop_cnt = 1;	/* observations processed since last transfer */

	for (i = 0; i < mp1; i++, istep++, loop_cnt++) {
		if (i == m) i = 0;	/* restart loop */

		/* no transfer during last loop; we're done */
		if (loop_cnt > m) {

			/* there was at least one quick transfer, so
			 * add affect clusters to live set */
			if (!(*indx))
				for (j = 0; j < k ;j++)
					/* affected clusters revived */
					if (live[j] > m)
						live[j] = m + 1;
			return;
		}

		/* no transfer if obs. i is sole member of cluster */
		if (nc[ic1[i]] == 1)
			continue;

		l1 = ic1[i];
		l2 = ic2[i];

		/* l1 changed; recompute membership cost */
		if (istep <= live[l1])
			d[i] = wgt ? weighted_cost_of_membership(nkt[l1], a[i],
				c[l1], c2[l1], n, nt, ABSOLUTE_SILENCE)
				: cost_of_membership(nkt[l1], a[i],
				c[l1], c2[l1], n);

		/* closest centers changed; transfer possible */
		if (istep <= live[l1] || istep <= live[l2]) {

			/* compute join cost to next closest */
			da = wgt ? weighted_cost_to_join(nkt[l2], a[i], c[l2],
				n, d[i], nt, ABSOLUTE_SILENCE)
				: cost_to_join(nkt[l2], a[i], c[l2],
				n, d[i]);

#ifdef __KMODES_DEBUGGING__
			/* verify cost of change and exit if incorrect */
			double de = debug_cost_to_move(a, m, n, k, ic1, i, l1,
								l2, wgt);
			debug_msg(1, 0, "Cost of moving member %zu from %u to "
						"%u is %.0f\n", i, l1, l2, de);
			if (de < 0 && da >= d[i]) {
				fprintf(stderr, "Cost to be in %u\n", l1);
				wgt ? weighted_cost_of_membership(nkt[l1], a[i],
					c[l1], c2[l1], n, nt, DEBUG_I)
					: cost_of_membership(nkt[l1], a[i],
					c[l1], c2[l1], n);
				fprintf(stderr, "Cost to join %u\n", l2);
				wgt ? weighted_cost_to_join(nkt[l2], a[i],
					c[l2], n, INFINITY, nt, DEBUG_I) :
					cost_to_join(nkt[l2], a[i], c[l2], n,
					INFINITY);
				summarize_clusters("optra0:", nkt, c, c2, k, n);
				fprintf(stderr, "Observation %u: ", i);
				fprint_data_ts(stderr, a[i], n, 1, 1);
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
					"from %u to %u (%.0f vs. %.0f)\n",
					l1, l2, de, da - d[i]));
			}
#endif

			/* less costly cluster */
			if (da < d[i] || (da == d[i] && l2 < l1)) {

				loop_cnt = 0;
				*indx = 0;
				live[l1] = istep + m;
				live[l2] = istep + m;

				/* update cluster counts & assignments */
				nc[l1]--;
				nc[l2]++;
				ic1[i] = l2;
				ic2[i] = l1;

				update_modes(nkt, a[i], c, c2, l1, l2, n, nj);
			}
		}


#ifdef __KMODES_DEBUGGING__
		/* verify modes and exit if incorrect */
		if (__nkjc_debug) {
			reset_nkjc(__nkjc_debug, k, n);
			compute_nkjc(__nkjc_debug, a, ic1, m, n);
			compare_nkjc(__nkjc, __nkjc_debug, k, n);
			compute_and_compare_modes(__nkjc_debug, c, c2, k, n);
		}
#endif

	}
} /* qtran */


/**
 * Compute cost for an observation to join new cluster.  There are four relevant
 * cases shown in the table below (see update_modes for others).  Note that if
 * the mode changes, then the weighted cost function needs to be recomputed in
 * its entirety.
 *
 *	Case	Mode	Tie	Minor	Tie	HD	Description
 *	----	----	----	-----	---	--	-----------
 *      0	old	-	old	-	1	add non-mode
 *	1	old	-	old	-	0	add to mode
 *	2 	new	break	new	-	0	add lo rank, untie mode
 *	3	old	form	-	-	1	add lo rank, tie mode
 *	 a			old	-		add to mmode
 *	 b			new	break		add lo rank, untie mmode
 *	 c			new	form		add hi rank, tie mmode
 *	4	new	form	new	-	1	add hi rank, tie mode
 *	5	old	none	new	form	0	add hi rank, tie mmode
 * 	-----------------------------------------------------------
 *	Mode: new indicate mode has changed, ow old
 *	Tie: break means a tie was broken, form means tie was formed, ow none
 *	HD: the unweighted cost to join at this site
 *	mmode = minor mode
 *
 * @param nkt	counts of each character at each position in proposed cluster
 * @param a	current observation
 * @param c	mode of proposed cluster
 * @param p	number of coordinates of data
 * @param bcost best achieved cost so far (branch-and-bound algorithm)
 * @return	correct cost if lower than bcost, otherwise some number > bcost
 */
static inline double cost_to_join(size_t **nkt, data_t *a, data_t *c, unsigned int p,
	double bcost)
{
	double de = 0.;

	for (unsigned int j = 0; j < p; ++j) {

		/* mode unchanged: add character less frequent than mode (0) */
		if (nkt[j][(int) c[j]] > nkt[j][(int) a[j]] + 1) {
			de += 1.;

		/* mode unchanged unless add hi rank:
		 *	unchanged: it costs to add new member (3)
		 *	  changed: it costs to keep old mode (4)
		 */
		} else if (nkt[j][(int) c[j]] == nkt[j][(int) a[j]] + 1) {
// [TODO] almost certainly superfluous?			&& c[j] != a[j])
			if (c[j] == a[j])	/* [TODO] delete if not used */
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
						"should not be possible!"));
			de += 1.;
		}

		/* higher cost than some previously known cost: abort */
		if (de > bcost)
			break;
	}
	return de;
} /* cost_to_join */


/**
 * As above, weighted version.
 */
static inline double weighted_cost_to_join(size_t **nkt, data_t *a, data_t *c,
	unsigned int p, double bcost, size_t **nt, int fxn_debug)
{
	double de = 0.;

	if (!nt)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "not possible"));

	for (unsigned int j = 0; j < p; ++j) {
		if (fxn_debug) {
			debug_msg(fxn_debug, 0, "%zu (%u):", j, a[j]);
			for (unsigned int l = 0; l < __nj[j]; ++l)
				fprintf(stderr, " %zu", nkt[j][l]);
			fprintf(stderr, " ");
		}

		/* mode unchanged: add character less frequent than mode (5) */
		if (nkt[j][(int) c[j]] > nkt[j][(int) a[j]] + 1) {
			de += (double) (nt[j][(int) c[j]] + nt[j][(int) a[j]])
				/ (nt[j][(int) c[j]] * nt[j][(int) a[j]]);
			debug_msg(QUIET <= fxn_debug, fxn_debug, "[1] %.0f\n",
									de);

		/* (w/ weights): add char that breaks tie (3-4) */
		} else if (nkt[j][(int) c[j]] == nkt[j][(int) a[j]] + 1) {
			/* add lo rank character that tie mode (3) */
			if (c[j] < a[j]) {
				de += (double) (nt[j][(int) c[j]] + nt[j][(int) a[j]])
					/ (nt[j][(int) c[j]] * nt[j][(int) a[j]]);
				debug_msg(QUIET <= fxn_debug, fxn_debug, "[2] "
								"%.0f\n", de);
			/* (w/ weights): add hi rank char. that ties mode (4) */
			} else if (a[j] < c[j]) {
				for (unsigned int l = 0; l < __nj[j]; ++l)
					if (l == (unsigned int) c[j])	/* previous mode */
						de += (double) nkt[j][l] *
							((nt[j][(int) a[j]] + nt[j][l])
							/(nt[j][(int) a[j]] * nt[j][l]));
					else if (l == (unsigned int) a[j])	/* new mode */
						de -= (double) nkt[j][l] *
							((nt[j][(int) c[j]] + nt[j][l])
							/(nt[j][(int) c[j]] * nt[j][l]));
					else
						de += (double) nkt[j][l] *
							((nt[j][(int) a[j]] + nt[j][l])
							/(nt[j][(int) a[j]] * nt[j][l])
							- (nt[j][(int) c[j]] + nt[j][l])
							/(nt[j][(int) c[j]] * nt[j][l]));
				debug_msg(QUIET >= fxn_debug, fxn_debug, "[3] "
								"%.0f\n", de);
			}

		/* (w/ weights): add low rank char that breaks tie (2) */
		} else if (nkt[j][(int) c[j]] == nkt[j][(int) a[j]]
			&& c[j] < a[j]) {
			for (unsigned int l = 0; l < __nj[j]; ++l)
				if (l == (unsigned int) c[j])	/* previous mode */
					de += (double) nkt[j][l] *
						((nt[j][(int) a[j]] + nt[j][l])
						/(nt[j][(int) a[j]] * nt[j][l]));
				else if (l == (unsigned int) a[j])	/* new mode */
					de -= (double) nkt[j][l] *
						((nt[j][(int) c[j]] + nt[j][l])
						/(nt[j][(int) c[j]] * nt[j][l]));
				else
					de += (double) nkt[j][l] *
						((nt[j][(int) a[j]] + nt[j][l])
						/(nt[j][(int) a[j]] * nt[j][l])
						- (nt[j][(int) c[j]] + nt[j][l])
						/(nt[j][(int) c[j]] * nt[j][l]));
			debug_msg(QUIET >= fxn_debug, fxn_debug, "[5] %.0f\n", de);
		}

		/* higher cost than some previously known cost: abort */
		if (de > bcost)
			break;
	}
	return de;
} /* weighted_cost_to_join */


/**
 * Compute cost of observation belonging to given cluster.
 *
 * @param nkt	count of each character in current cluster
 * @param a	observation
 * @param c	cluster mode
 * @param c2	cluster minor mode
 * @param p	number of coordinates of data
 * @return	cost of membership
 */
static inline double cost_of_membership(size_t **nkt, data_t *a, data_t *c,
	data_t *c2, unsigned int p)
{
	double de = 0.;
	for (unsigned int j = 0; j < p; ++j) {

		/* removing non-modal character in ith observation */
		if (a[j] != c[j])
			de += 1.;

		/* removing modal character and departure will change mode
		 * third condition required for sites with one category
		 */
		else if (a[j] == c[j] && nkt[j][(int) c[j]]
			== nkt[j][(int) c2[j]] && c[j] < c2[j])
			de += 1.;
	}
	return de;
} /* cost_of_membership */


/**
 * Compute cost of observation belonging to given cluster.
 *
 * @param nkt	count of each character in current cluster
 * @param a	observation
 * @param c	cluster mode
 * @param c2	cluster minor mode
 * @param p	number of coordinates of data
 * @param nt	total counts of characters at each position (optional)
 * @param wgt	weighted distances
 * @param fxn_debug
 * @return	cost of membership
 */
static inline double weighted_cost_of_membership(size_t **nkt, data_t *a,
	data_t *c, data_t *c2, unsigned int p, size_t **nt, int fxn_debug)
{
	double de = 0.;
	for (unsigned int j = 0; j < p; ++j) {
		if (fxn_debug) {
			fprintf(stderr, "%u (%u):", j, a[j]);
			for (unsigned int l = 0; l < __nj[j]; ++l)
				fprintf(stderr, " %zu", nkt[j][l]);
			fprintf(stderr, " ");
		}

		/* removing non-modal character in ith observation */
		if (a[j] != c[j]) {
			de += (double) (nt[j][(int) c[j]] + nt[j][(int) a[j]])
				/ (nt[j][(int) c[j]] * nt[j][(int) a[j]]);
			debug_msg(QUIET >= fxn_debug, fxn_debug, "remove "
				"non-mode %u (%.0f)\n", a[j], de);

		/* removing modal character and departure will change mode */
		} else if (a[j] == c[j] && nkt[j][(int) c[j]]
			== nkt[j][(int) c2[j]] && c[j] < c2[j]) {
			de += (double) (nt[j][(int) c[j]] + nt[j][(int) c2[j]])
				/ (nt[j][(int) c[j]] * nt[j][(int) c2[j]]);
			debug_msg(QUIET >= fxn_debug, fxn_debug, "remove mode "
				"%u, new mode %u (%.0f)\n", a[j], c2[j], de);
		}
	}
	return de;
} /* weighted_cost_of_membership */


/**
 * Compute cost criterion. [PARALLEL]
 *
 * @param x		data (n x p)
 * @param seeds		current seeds (K x p)
 * @param ic		cluster assignments (n x 1)
 * @param criterion	per cluster criterion to compute (K x 1) [optional]
 * @param K		number of clusters
 * @param n		number of observations
 * @param p		number of coordinates
 * @param wgt		weighted HD
 * @return		return sum of cluster costs
 */
static double *__criterion = NULL;	/* used if memory not provided */
double compute_criterion(data_t **x, data_t **seeds, unsigned int *ic,
	double *criterion, unsigned int K, unsigned int n, unsigned int p, int wgt)
{
	double *lcriterion;
	if (!criterion) {
		if (!__criterion) {
			__criterion = malloc(K * sizeof *__criterion);
			if (!__criterion)
				exit(mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"__criterion"));
		}
		lcriterion = __criterion;
	} else
		lcriterion = criterion;

	/* reset per-cluster cost */
	for (unsigned int k = 0; k < K; ++k)
		lcriterion[k] = 0.;

	/* compute per-cluster cost */
	for (unsigned int i = 0; i < n; ++i)	/* [PARALLEL] */
		lcriterion[ic[i]] += hd(x[i], seeds[ic[i]], p, wgt);

	/* compute total cost */
	double sum = 0;
	for (unsigned int k = 0; k < K; ++k)
		sum += lcriterion[k];

	return sum;
} /* compute_criterion */


/**
 * Compute the cost from sufficient statistics.
 *
 * @param nkjc	category counts per cluster per site (k x p x __nj)
 * @param cost	cost (k x 1)
 * @param c	modes (k x p)
 * @param K	number of clusters
 * @param p	number of coordinates
 * @return	cost
 */
static inline double compute_cost(size_t ***nkjc, double *cost, data_t **c,
	unsigned int K, unsigned int p)
{
	double tcost = 0;

	if (!__nj)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj = %p; __njc = "
			"%p!\n", (void *)__nj, (void *)__njc));

	for (unsigned int k = 0; k < K; ++k)
		cost[k] = 0.;
	for (unsigned int k = 0; k < K; ++k) {
		for (unsigned int j = 0; j < p; ++j) {

			for (unsigned int l = 0; l < c[k][j]; ++l)
				cost[k] += nkjc[k][j][l];
			for (unsigned int l = c[k][j] + 1; l < __nj[j]; ++l)
				cost[k] += nkjc[k][j][l];
		}
		tcost += cost[k];
	}

	return tcost;
} /* compute_cost */


/**
 * Compute the weighted cost from sufficient statistics.
 *
 * @param nkjc	category counts per cluster per site (k x p x __nj)
 * @param cost	cost (k x 1)
 * @param c	modes (k x p)
 * @param K	number of clusters
 * @param p	number of coordinates
 * @return	cost
 */
static inline double compute_weighted_cost(size_t ***nkjc, double *cost,
	data_t **c, unsigned int K, unsigned int p)
{
	double tcost = 0;

	if (!__nj || !__njc)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj = %p; __njc = "
			"%p!\n", (void *)__nj, (void *)__njc));

	for (unsigned int k = 0; k < K; ++k) {
		for (unsigned int j = 0; j < p; ++j) {

			for (unsigned int l = 0; l < c[k][j]; ++l)
				cost[k] += nkjc[k][j][l] * ((double)
					(__njc[j][(int) c[k][j]] + __njc[j][(int) l])
					/ (__njc[j][(int) c[k][j]] * __njc[j][(int) l]));
			for (unsigned int l = c[k][j] + 1; l < c[k][j]; ++l)
				cost[k] += nkjc[k][j][l] * ((double)
					(__njc[j][(int) c[k][j]] + __njc[j][(int) l])
					/ (__njc[j][(int) c[k][j]] * __njc[j][(int) l]));
		}
		tcost += cost[k];
	}

	return tcost;
} /* compute_weighted_cost */


/**
 * Update modes after a transfer.  Suppose c is the character moving from
 * source to target.  Process the rules in the table in order, such that earlier
 * cases are eliminated before considering later cases.
 *
 *			Source
 *	Case	Mode	Tie	Minor	Tie	Description
 *	----	----	----	-----	---	-----------
 *	a	new	make	-	-	remove lo rank, tie w/ mmode
 *	 1			old	make	old mode becomes mmode
 *	 2			new	make	tie w/ new lo rank
 *	b	new	break	?	?	old mode becomes mmode
 * 	---------------------------------------------------
 *
 *			Target
 *	Case	Mode	Tie	Minor	Tie	Description
 *	----	----	----	-----	---	-----------
 *	1	old	-	old	-	add to mode
 *	2 	new	break	new	-	add lo rank, untie mode
 *	3	old	make	-	-	add lo rank, tie mode
 *	3a			old	-	add minor mode
 *	3b			new	break	add lo rank, untie mmode
 *	3c			new	make	add hi rank, tie mmode
 *	4	new	make	new	-	add hi rank, tie mode
 *	5	old	none	new	make	add hi rank, tie mmode
 * 	---------------------------------------------------
 *
 * @param nkt	counts of categories in each cluster
 * @param a	transferred observation
 * @param c	current modes
 * @param c2	current minor modes
 * @param l1	cluster transferred from
 * @param l2	cluster transferred to
 * @param p	number of coordinates
 * @param nj	number of categories at each coordinate
 */
static inline void update_modes(size_t ***nkt, data_t *a, data_t **c, data_t **c2,
	unsigned int l1, unsigned int l2, unsigned int p, unsigned int *nj)
{
	data_t ll, l;

	/* update centers [PARALLEL] */
	for (unsigned int j = 0; j < p; j++) {
		/* update category counts */
		nkt[l1][j][(int) a[j]]--;
		nkt[l2][j][(int) a[j]]++;

		/* update 1st and 2nd modes of affected clusters */

		/* cluster source */
		/* new mode if mode demoted to equal minor mode of lower rank */
		if (a[j] == c[l1][j] && nkt[l1][j][(int) c[l1][j]] ==
			nkt[l1][j][(int) c2[l1][j]] && c2[l1][j] < c[l1][j]) {
			ll = c[l1][j];
			c[l1][j] = c2[l1][j];

			/* minor mode is old mode unless there is earlier tie */
			c2[l1][j] = ll;
			for (l = c[l1][j] + 1; l < ll; ++l)
				if (nkt[l1][j][l] == nkt[l1][j][ll]) {
					c2[l1][j] = l;
					break;
				}
		/* new mode if mode demoted below minor mode */
		} else if (a[j] == c[l1][j] && nkt[l1][j][(int) c[l1][j]]
			== nkt[l1][j][(int) c2[l1][j]] - 1) {
			ll = c[l1][j];
			c[l1][j] = c2[l1][j];

			/* look for minor mode
			 * with same frequency as last minor mode */
			for (l = c2[l1][j] + 1; l < nj[j]; ++l)
				if (nkt[l1][j][l]
					== nkt[l1][j][(int) c2[l1][j]]) {
					c2[l1][j] = l;
					break;
				}

			/* or a lower rank minor mode w/ same frequency as a[j] */
			if (l == nj[j]) {
				c2[l1][j] = ll;
				for (l = 0; l < ll; ++l)
					if (nkt[l1][j][l] == nkt[l1][j][ll]) {
						c2[l1][j] = l;
						break;
					}
			}

		/* new 2nd mode if 2nd mode demoted below a 3rd mode */
		} else if (a[j] == c2[l1][j]) {
			for (l = 0; l < nj[j]; ++l) {

				/* higher mode that is not 1st mode: it was
				 * tied with c2, take the first */
				if (nkt[l1][j][l] > nkt[l1][j][(int) c2[l1][j]]
					&& l != c[l1][j]) {
					c2[l1][j] = l;
					break;

				/* equal mode that is earlier in sequence: by
				 * convention take the first, but keep looking
				 * for a higher mode */
				} else if (nkt[l1][j][l]
					== nkt[l1][j][(int) c2[l1][j]]
					&& l < c2[l1][j])
					c2[l1][j] = l;
			}
		}

		/* cluster target */
		/* new 1st mode or new mode ties with 1st and has lower rank */
		if (nkt[l2][j][(int) a[j]] > nkt[l2][j][(int) c[l2][j]]
			|| (nkt[l2][j][(int) a[j]] == nkt[l2][j][(int) c[l2][j]]
			&& a[j] < c[l2][j])) {
			c2[l2][j] = c[l2][j];
			c[l2][j] = a[j];
		}

		/* new 2nd mode */
		else if (a[j] != c[l2][j] && (nkt[l2][j][(int) a[j]]	/* CORRECT */
		//else if ((nkt[l2][j][(int) a[j]]	/* BUGGY */
			> nkt[l2][j][(int) c2[l2][j]] ||
			(nkt[l2][j][(int) a[j]] == nkt[l2][j][(int) c2[l2][j]]
			&& a[j] < c2[l2][j]))) {
			c2[l2][j] = a[j];
		}
	}
} /* update_modes */


/**
 * Huang 1997 k-modes algorithm modified to allow variations on initialization.
 * By default, the algorithm takes K seeds and assigns observations in input
 * order to seeds, updating the seeds with every addition.  However, this does
 * not match the original k-means algorithm, which postpones mean updates until
 * all observations have been assigned.  To mimic the latter initialization, set
 * init_update to 0.
 *
 * @param x		data (n X p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param seeds		initial seeds (may not be set)
 * @param K		number of clusters
 * @param ic1		initial partition (may not be set)
 * @param nclass	size of each cluster (to set)
 * @param max_iter	maximum allowed iterations
 * @param cost		per-cluster cost
 * @param ifault	error status (to set)
 * @param iter		number of iterations (to set)
 * @param opt		pointer to structure with run options
 *	weight		use weights
 *	init_update	update modes during initialization
 */
double
kmodes_huang(data_t **x, unsigned int n, unsigned int p, data_t **seeds,
	unsigned int K, unsigned int *ic1, unsigned int *nclass,
	unsigned int max_iter, double *cost, int *ifault,
	unsigned int *iter, kmodes_options *opt)
// int weight, int init_update)
{

	unsigned int i, j, l=0, m=0;	/* data indices */
	unsigned int k;		/* cluster index */
	size_t max_nkjc;	/* new mode count */
	data_t max;		/* new mode */
	double dmin, d;
	char keep_going = 1;
	int use_minor_mode = 0;	/* use update_modes */

	if (opt->use_hartigan)
		use_minor_mode = 1;

	*ifault = KMODES_NO_ERROR;

	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		return INFINITY;
	}

	/* make sure category counts summary statistics are allocated and set */
	if (!(l = allocate_and_compute_category_counts(x, n, p, K,
							opt->weighted))) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

	/* trivial input: user requests 1 cluster */
	if (K == 1) {
		nclass[0] = n;
		/* assign all observations to cluster 1 [PARALLEL] */
		for (i = 0; i < n; i++)
			ic1[i] = 0;
		/* compute centers [PARALLEL] */
		for (j = 0; j < p; j++) {
			unsigned int max = 0;
			for (l = 0; l < __nj[j]; l++) {
				__njc[j][l] = 0;
				for (i = 0; i < n; i++)
					__njc[j][(int) x[i][j]]++;
				if (max < __njc[j][l]) {
					max = __njc[j][l];
					seeds[0][j] = l;
				}
			}
		}
		/* compute sum-of-distances [PARALLEL] */
		cost[0] = 0;
		return compute_criterion(x, seeds, ic1, cost, K, n, p,
							opt->weighted);

	/* trivial input: user requests m clusters */
	} else if (K == n) {
		for (i = 0; i < n; i++) {
			nclass[i] = 1;	/* 1 member per cluster */
			memcpy(seeds[i], x[i], p * sizeof **seeds);
			ic1[i] = i;
			cost[i] = 0;
		}

		return 0.;
	}

	/* now the interesting case of 1 < k < m */

#ifdef __KMODES_DEBUGGING__
	allocate_nkjc_debug(K, p, l);
#endif

	/* allocate or reset category counts per site per cluster */
	/* l = total number categories */
	if (set_nkjc(K, p, l)) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

	if (use_minor_mode && !__c2) {
		__c2 = malloc(K * sizeof *__c2);
		/* allocate c2 memory as a single block */
		data_t *tmp = malloc(K * p * sizeof **__c2);

		if (!__c2 || !tmp) {
			*ifault = KMODES_MEMORY_ERROR;
			return INFINITY;
		}

		for (k = 0; k < K; ++k) {
			__c2[k] = tmp;
			tmp += p;

			/* initialize minor mode to be NOT initial modes */
			for (j = 0; j < p; j++)
				__c2[k][j] = __nj[j] > 1 && seeds[k][j] == 0 ? 1 : 0;
		}
	} else if (use_minor_mode) {
		for (k = 0; k < K; ++k)
			for (j = 0; j < p; j++)
				__c2[k][j] = __nj[j] > 1 && seeds[k][j] == 0 ? 1 : 0;
	}

	for (k = 0; k < K; ++k) {
		nclass[k] = 0;
		cost[k] = 0.;
	}

	/* assign to initial clusters */
	for (i = 0; i < n; i++) {		/* for each observation */

		/* must assign to a partition */
		dmin = INFINITY;
		for (l = 0; l < K; ++l) {	/* consider each mode */
			d = hd(x[i], seeds[l], p, opt->weighted);

			/* tie: assign to lower index cluster (klaR) */
			if (d < dmin) {
				k = l;
				dmin = d;
			}
		}

		ic1[i] = k;		/* assign to cluster k */
		cost[k] += dmin;	/* increment cost */
		nclass[k]++;		/* increment cluster count */

		/* update type counts in cluster k & affected mode */
		for (j = 0; j < p; ++j) {
			__nkjc[ic1[i]][j][(int) x[i][j]]++;

			/* update affected mode: as per Huang1997 */
			if (opt->init_update) {
				if (use_minor_mode) {
					if (__nkjc[k][j][(int) x[i][j]] >
						__nkjc[k][j][(int) seeds[k][j]]) {
						__c2[k][j] = seeds[k][j];
						seeds[k][j] = x[i][j];
					} else if (__nkjc[k][j][(int) x[i][j]]
						== __nkjc[k][j][(int) seeds[k][j]]
						&& x[i][j] < seeds[k][j]) {
						__c2[k][j] = seeds[k][j];
						seeds[k][j] = x[i][j];
					} else if (x[i][j] != seeds[k][j]
						&& __nkjc[k][j][(int) x[i][j]]
						> __nkjc[k][j][(int) __c2[k][j]]) {
						__c2[k][j] = x[i][j];
					} else if (x[i][j] != seeds[k][j]
						&& __nkjc[k][j][(int) x[i][j]]
						== __nkjc[k][j][(int) __c2[k][j]]
						&& x[i][j] < __c2[k][j]) {
						__c2[k][j] = x[i][j];
					}
				} else {
					if (__nkjc[k][j][(int) x[i][j]]
						> __nkjc[k][j][(int) seeds[k][j]])
						seeds[k][j] = x[i][j];
					else if (__nkjc[k][j][(int) x[i][j]]
						== __nkjc[k][j][(int) seeds[k][j]]
						&& x[i][j] < seeds[k][j])
						seeds[k][j] = x[i][j];
				}
			}
		}
	}

	/* update or compute modes, if not already done */
	if (!opt->init_update && compute_modes(__nkjc, seeds, use_minor_mode
		? __c2 : NULL, K, p)) {
		*ifault = KMODES_INTERNAL_ERROR;
		return INFINITY;
	}

	for (k = 0; k < K; ++k)
		if (nclass[k] == 0) {
			*ifault = KMODES_NULL_CLUSTER_ERROR;
			return INFINITY;
		}
	
	if (!max_iter)
		return opt->weighted
			? compute_weighted_cost(__nkjc, cost, seeds, K, p)
			: compute_cost(__nkjc, cost, seeds, K, p);

	/* continue iterating until no updates possible */
	while (m++ < max_iter && keep_going) {

		keep_going = 0;
		for (i = 0; i < n; i++) {

			/* decide if going to move to another cluster */

			/* hartigan: verify decrease in criterion */
			if (opt->use_hartigan) {
				dmin = opt->weighted
					? weighted_cost_of_membership(
						__nkjc[ic1[i]], x[i],
						seeds[ic1[i]], __c2[ic1[i]],
						p, __njc, ABSOLUTE_SILENCE)
					: cost_of_membership(__nkjc[ic1[i]],
						x[i], seeds[ic1[i]],
						__c2[ic1[i]], p);
				k = ic1[i];
				for (l = 0; l < K; ++l) {
					if (l == k)
						continue;
					d = opt->weighted
						? weighted_cost_to_join(
							__nkjc[l], x[i],
							seeds[l], p, dmin,
							__njc, ABSOLUTE_SILENCE)
						: cost_to_join(__nkjc[l], x[i],
							seeds[l], p, dmin);
					if (d < dmin) {
						k = l;
						dmin = d;
					}
				}

			/* huang (MacQueen): just check distance to clusters */
			} else {
				dmin = INFINITY;
				for (l = 0; l < K; l++) {
					d = hd(x[i], seeds[l], p,
							opt->weighted);
					if (d < dmin) {
						k = l;
						dmin = d;
					}
				}
			}

			/* closest mode has changed */
			/* note: reassign to lower index cluster on tie (klaR) */
			if (k != ic1[i]) {

				keep_going = 1;

				/* update cluster counts */
				nclass[ic1[i]]--;
				nclass[k]++;

/* appears to be slower, but we need this code for hartigan updates
*/
				if (use_minor_mode) {
					update_modes(__nkjc, x[i], seeds, __c2,
						ic1[i], k, p, __nj);

				} else {

				/* update affected modes & coordinate counts */
				for (j = 0; j < p; ++j) {

					/* coordinate counts */
					__nkjc[k][j][(int) x[i][j]]++;
					__nkjc[ic1[i]][j][(int) x[i][j]]--;

					/* new mode */
					max = 0;
					max_nkjc = __nkjc[k][j][(int) max];
					for (l = 1; l < __nj[j]; ++l)
						if (__nkjc[k][j][l] > max_nkjc) {
							max = l;
							max_nkjc = __nkjc[k][j][
								(int) max];
						}
					seeds[k][j] = max;

					/* previous mode */
					max = 0;
					max_nkjc = __nkjc[ic1[i]][j][(int) max];
					for (l = 1; l < __nj[j]; ++l)
						if (__nkjc[ic1[i]][j][l]
							> max_nkjc) {
							max = l;
							max_nkjc = __nkjc[ic1[i]]
								[j][(int) max];
						}
					seeds[ic1[i]][j] = max;
				}
				}

				/* record new cluster id */
				ic1[i] = k;

#ifdef __KMODES_DEBUGGING__
				if (__nkjc_debug) {
					reset_nkjc(__nkjc_debug, K, p);
					compute_nkjc(__nkjc_debug, x, ic1, n, p);
					compare_nkjc(__nkjc, __nkjc_debug, K, p);
					compute_and_compare_modes(__nkjc_debug,
						seeds, use_minor_mode ? __c2
						: NULL, K, p);
				}
#endif

			}
		}
	}
	m--;

	if (m > max_iter && max_iter)
		*ifault = KMODES_EXCEED_ITER_WARNING;

	/* error if initialization routine produces empty cluster */
	for (k = 0; k < K; ++k)
		if (nclass[k] == 0) {
			*ifault = KMODES_NULL_CLUSTER_ERROR;
			return INFINITY;
		}

	/* summing d[i] is not the same thing */
	*iter = m;

	/* recompute criterion */
	return opt->weighted
		? compute_weighted_cost(__nkjc, cost, seeds, K, p)
		: compute_cost(__nkjc, cost, seeds, K, p);

} /* kmodes_huang */


/**
 * Lloyd's algorithm.  Update modes after all observations processed.
 *
 * @param x		data (n x p)
 * @param c		initial modes (k x p)
 * @param nclass	count in clsuters (1 x k)
 * @param ic1		cluster assignments (1 x n)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param max_iter	maximum allowed iterations
 * @param cost		cost per cluster
 * @param ifault	pointer to error status
 * @param iter		pointer to number of iterations
 * @param weight	used weighted distances
 */
double kmodes_lloyd(data_t **x, data_t **c, unsigned int *nclass, unsigned int *ic1,
	unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter,
	double *cost, int *ifault, unsigned int *iter, int weight)
{
	unsigned int i, j, k, l=0, m=0;
	size_t max_nkjc;
	double dmin, d;
	char keep_going = 1;
	*ifault = KMODES_NO_ERROR;

	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		return INFINITY;
	}

	/* make sure category counts summary statistics are allocated and set */
	if (!(l = allocate_and_compute_category_counts(x, n, p, K, weight))) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

	/* allocate category counts per site per cluster */
	/* l = total number categories */
	if (set_nkjc(K, p, l)) {
		*ifault = KMODES_MEMORY_ERROR;
		return INFINITY;
	}

#ifdef __KMODES_DEBUGGING__
	allocate_nkjc_debug(K, p, l);
#endif
    
	/* initialize cluster counts */
	for (k = 0; k < K; ++k)
		nclass[k] = 0;

	/* assign to initial clusters and count things */
	for (i = 0; i < n; i++) {		/* for each observation */
		dmin = INFINITY;
		for (l = 0; l < K; ++l) {	/* compare to each mode */
			d = hd(x[i], c[l], p, weight);
			if (d < dmin) {
				k = l;
				dmin = d;
			}
		}
		nclass[k]++;
		ic1[i] = k;
		for (j = 0; j < p; ++j)
			__nkjc[k][j][(int) x[i][j]]++;
	}

	/* update modes */
	for (k = 0; k < K; ++k)			/* for each cluster */
		for (j = 0; j < p; ++j) {	/* for each coordinate */
			c[k][j] = 0;
			max_nkjc = __nkjc[k][j][(int) c[k][j]];
			for (l = 1; l < __nj[j]; ++l)
				if (__nkjc[k][j][l] > max_nkjc) {
					c[k][j] = l;
					max_nkjc =
						__nkjc[k][j][(int) c[k][j]];
				}
		}

	while (m++ < max_iter && keep_going) {
		keep_going = 0;
		for (i = 0; i < n; i++) {
			/* identify closest mode */
			dmin = INFINITY;
			for (l = 0; l < K; l++) {
				d = hd(x[i], c[l], p, weight);
				if (d < dmin) {
					k = l;
					dmin = d;
				}
			}
			if (ic1[i] != k) {
				keep_going = 1;
				nclass[ic1[i]]--;
				nclass[k]++;
				for (j = 0; j < p; ++j) {
					__nkjc[k][j][(int) x[i][j]]++;
					__nkjc[ic1[i]][j][(int) x[i][j]]--;
				}
				ic1[i] = k;
			}
		}

		if (!keep_going)
			break;

		/* update all modes */
		for (k = 0; k < K; ++k)
			for (j = 0; j < p; ++j) {
				c[k][j] = 0;
				max_nkjc = __nkjc[k][j][(int) c[k][j]];
				for (l = 1; l < __nj[j]; ++l)
					if (__nkjc[k][j][l] > max_nkjc) {
						c[k][j] = l;
						max_nkjc = __nkjc[k][j][
							(int) c[k][j]];
					}
			}
	}

	if (m > max_iter)
		*ifault = KMODES_EXCEED_ITER_WARNING;

	*iter = m;

	return compute_criterion(x, c, ic1, cost, K, n, p, weight);

} /* kmodes_lloyd */


/**
 * Compute modes given a hard clustering.  Calculation is done based on
 * precomputed category counts in variable nkjc.  Pass NULL if you want
 * just modes or mmodes.
 *
 * @param nkjc		counts of categories at each site in each cluster
 * @param modes		modes to compute (K x p)
 * @param mmodes	minor modes to compute (K x p)
 * @param K		number of clusters
 * @param p		number of coordinates
 * @return		error status
 */
static inline int compute_modes(size_t ***nkjc, data_t **modes, data_t **mmodes,
	unsigned int K, unsigned int p)
{
	if (!nkjc || !__nj)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "invalid call");

	for (unsigned int k = 0; k < K; ++k)
		for (unsigned int j = 0; j < p; ++j) {
			size_t max_nkjc = 0, max2_nkjc = 0;
			unsigned int max_l = 0, max2_l = 0;
			for (unsigned int l = 0; l < __nj[j]; ++l) {
				if (nkjc[k][j][l] > max_nkjc) {
					max2_nkjc = max_nkjc;
					max2_l = max_l;
					max_nkjc = nkjc[k][j][l];
					max_l = l;
				} else if (nkjc[k][j][l] > max2_nkjc) {
					max2_nkjc = nkjc[k][j][l];
					max2_l = l;
				}
			}
			if (modes)
				modes[k][j] = max_l;
			if (mmodes)
				mmodes[k][j] = max2_l;
		}
	return NO_ERROR;
} /* compute_modes */


/**
 * Hamming distance function.
 *
 * @param x	observation (1 x p)
 * @param y	observation (1 x p)
 * @param p	number of coordinates
 * @param wgt	use weights?
 * @return	(weighted) Hamming distance
 */
static inline double hd(data_t *x, data_t *y, unsigned int p, int wgt)
{
	double d = 0;

	for (unsigned int j = 0; j < p; ++j)
		d += hd_coord(x, y, j, wgt);	//x[j] != y[j];

	return(d);
} /* hd */


/**
 * Hamming distance branch and bound, aborting as soon as distance exceeds
 * previously computed one.
 *
 * @param x	observation (1 x p)
 * @param y	observation (1 x p)
 * @param p	number of coordinates
 * @param wgt	use weights?
 * @param dprev	previous distance
 * @return	(weighted) Hamming distance
 */
static inline double hd_min(data_t *x, data_t *y, unsigned int p, int wgt,
	double dprev)
{
	double d = 0;

	for (unsigned int j = 0; j < p; ++j) {
		d += hd_coord(x, y, j, wgt);
		if (d >= dprev)
			return dprev;
	}

	return d;
} /* hd_min */


/**
 * Hamming distance contribution from single coordinate.
 *
 * @param x	observation (1 x p)
 * @param y	bservation (1 x p)
 * @param j	selected coordinate
 * @param wgt	use weights?
 * @return	(weighted) Hamming distance contribution from coordinate j.
 */
static inline double hd_coord(data_t *x, data_t *y, unsigned int j, int wgt)
{
	if (wgt && !__njc) {
		mmessage(ERROR_MSG, INTERNAL_ERROR, "function called without "
			"setup of __njc");
		exit(EXIT_FAILURE);
	}

	return( x[j] == y[j] ? 0. : wgt
		? (double) (__njc[j][(int) x[j]] + __njc[j][(int) y[j]])
			/ (__njc[j][(int) x[j]] * __njc[j][(int) y[j]])
		: 1.);
} /* hd_coord */


/**
 * Create kmodes options object and initialize all choices to default.
 *
 * @param opt	kmodes_options pointer
 * @return	error status
 */
int make_kmodes_options(kmodes_options **opt)
{
	*opt = malloc(sizeof **opt);
	if (!*opt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"kmodes_options");

	(*opt)->weighted = 0;
	(*opt)->init_update = 0;
	(*opt)->use_qtran = 0;
	(*opt)->use_hartigan = 0;

	return NO_ERROR;
} /* make_kmodes_options */


/**
 * Initialize k-modes from partition.
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param k		number of clusters
 * @param wgt		use weights
 * @param seeds		initial modes to be chosen by this function (k x p)
 * @param ic		cluster assignments (1 x n)
 * @return		error status
 */
int kmodes_init_from_partition(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	int wgt, data_t **seeds, unsigned int *ic)
{
	unsigned int l;

	if (!(l = allocate_and_compute_category_counts(x, n, p, K, wgt)))
		return KMODES_MEMORY_ERROR;

	if (set_nkjc(K, p, l))
		return KMODES_MEMORY_ERROR;

	compute_nkjc(__nkjc, x, ic, n, p);

	return compute_modes(__nkjc, seeds, NULL, K, p);

} /* kmodes_init_from_partition */


/**
 * Initialize k-modes.
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param k		number of clusters
 * @param k1		number of seeds deterministically selected
 * @param seeds		initial modes to be chosen by this function
 * @param sidx		seed indices to be set (if initializer's seeds are observations)
 * @param method 	method of initialization (see kmodes.h)
 * @param weight	use weighted Hamming distances in cost function
 * @return		error status
 */
int kmodes_init(data_t **x, unsigned int n, unsigned int p, unsigned int k,
	unsigned int k1, data_t **seeds, unsigned int *sidx, int method,
	int weight)
{
	if (method == KMODES_INIT_RANDOM_SEEDS) {
		return kmodes_init_random_seeds(x, n, p, k, k1, seeds, sidx);
	} else if (method == KMODES_INIT_H97) {
		return kmodes_init_h97(x, n, p, k, k1, weight, 0, seeds, sidx);
	} else if (method == KMODES_INIT_H97_RANDOM) {
		return kmodes_init_h97(x, n, p, k, k1, weight, 1, seeds, sidx);
	} else if (method == KMODES_INIT_HD17) {
		return kmodes_init_hd17(x, n, p, k, k1, weight, seeds, sidx);
	} else if (method == KMODES_INIT_CLB09) {
		return kmodes_init_clb09(x, n, p, k, k1, weight, 0, seeds, sidx,UNSINGED_INT);
	} else if (method == KMODES_INIT_CLB09_RANDOM) {
		return kmodes_init_clb09(x, n, p, k, k1, weight, 1, seeds, sidx,UNSINGED_INT);
	} else if (method == KMODES_INIT_AV07) {
		return kmodes_init_av07(x, n, p, k, k1, weight, seeds, sidx, 0);
	} else if (method == KMODES_INIT_AV07_GREEDY) {
		return kmodes_init_av07(x, n, p, k, k1, weight, seeds, sidx, 1);
	} else if (method == KMODES_INIT_USER_SEEDS) {
		return NO_ERROR;	/* already seeded */
	} else
		return mmessage(ERROR_MSG, KMODES_INVALID_INITIALIZATION_METHOD,
			kmodes_error(KMODES_INVALID_INITIALIZATION_METHOD));
	return NO_ERROR;
} /* kmodes_init */

/**
 * Randomly choose seeds from a known partition.  If options::K is larger than
 * the number of partitions, uniformly sample partitions, then distinct seeds
 * from within
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of seeds to select
 * @param seeds		place to put seeds (K x p)
 * @param sd_idx	place to store selected seed indices (1 x K)
 * @param id		partition (1 x n)
 * @return		error status
 */
int kmodes_init_random_from_partition(data_t **x, unsigned int n, unsigned int p,
	 unsigned int K, data_t **seeds, unsigned int *sd_idx, unsigned int *id)
{
	unsigned int true_k = 0;
	unsigned int *nc = NULL;
	unsigned int *sidx = NULL;
	int same;

	if (n < K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Requesting %u clusters with only %u reads.\n", K, n);

	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);

		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"seed index");
	} else
		sidx = sd_idx;

	/* determine the true number of clusters */
	for (unsigned int i = 0; i < n; ++i)
		if (id[i] > true_k)
			true_k = id[i];
	++true_k;

	/* count the observations in each cluster */
	nc = calloc(true_k, sizeof *nc);
	if (!nc)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nc");

	for (unsigned int i = 0; i < n; ++i)
		++nc[id[i]];

	/* more than K clusters */
	if (K >= true_k) {

		/* pick one seed from each cluster */
		for (unsigned int k = 0; k < true_k; ++k) {
			unsigned int j = (unsigned int)((double) rand()
				/ RAND_MAX * nc[k]);
			sidx[k] = j;
			memcpy(seeds[k], x[j], p * sizeof **x);
		}

		/* remainding seeds */
		for (unsigned int k = true_k; k < K; ++k) {
			/* select cluster with replacement */
			unsigned int m = (unsigned int)((double) rand()
				/ RAND_MAX * K);
			if (m == K) --m;

			/* select distinct observation */
			do {
				same = 0;
				sidx[k] = (unsigned int) ((double) rand()
					/ RAND_MAX * nc[k]);
				if (sidx[k] == nc[k]) --sidx[k];
				for (unsigned l = 0; l < sidx[k]; ++l)
					if (sidx[l] == sidx[k] ||
						!compare_data(x, sidx[l],
						sidx[k], p)) {
						same = 1;
						break;
					}
				memcpy(seeds[k], x[sidx[k]], p * sizeof **x);
			} while (same);
		}
	} else if (K < true_k) {
		/* select clusters without replacement */
		unsigned int *sel_k = malloc(K * sizeof *sel_k);
		sample(true_k, K, sel_k);
		for (unsigned int k = 0; k < K; ++k) {
			unsigned int j = (unsigned int)((double) rand()
				/ RAND_MAX * nc[sel_k[k]]);
			if (sidx) sidx[k] = j;
			memcpy(seeds[k], x[j], p * sizeof **x);
		}
		free(sel_k);
	}

	free(nc);

	return NO_ERROR;
} /* kmodes_init_random_from_partition */


/**
 * Initialize by randomly selecting seeds from a set.
 *
 * @param K		number of seeds to select
 * @param p		number of coordinates
 * @param n_ss		number of seeds in seedset
 * @param seeds		seeds to set (K x p)
 * @param seedset	choices to select from (n_ss x p)
 * @return		error status
 */
int kmodes_init_random_from_set(unsigned int K, unsigned int p,
	unsigned int n_ss, data_t **seeds, data_t **seedset)
{

	if (K >= n_ss)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Need seed set "
			"size (%u) to exceed K=%u\n", n_ss, K);

	/* if triggered check p*sizeof(data_t) can store uint & use pointers */
	if (n_ss > pow(2, 8*sizeof(data_t)))
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "");

	unsigned int k = 0, t = 0;

	while (k < K) {
		double u = rand() / (RAND_MAX + 1.);

		if ( (n_ss - t) * u >= K - k )
			++t;
		else	/* [TODO] assumes K fits in data_t */
			seeds[k++][0] = (data_t) t++;
	}

	for (k = 0; k < K; ++k)
		memcpy(seeds[k], seedset[(size_t) seeds[k][0]], p * sizeof **seeds);

	return NO_ERROR;
} /* kmodes_init_random_from_set */


/**
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param k1		number of seeds already selected
 * @param seeds		seeds (to be set)
 * @param sd_idx	indices of chosen observations
 * @return		error status
 */
int kmodes_init_random_seeds(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int k1, data_t **seeds, unsigned int *sd_idx)
{
	unsigned int *sidx = NULL;
	int same;

	if (n < K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Requesting %u clusters with only %u reads.\n", K, n);

	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);

		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"seed index");
	} else
		sidx = sd_idx;
	for (unsigned j = 0; j < k1; ++j) {
		sidx[j] = n;
		for (unsigned int i = 0; i < n; ++i)
			if (!compare_data_to_seed(x, i, seeds[j], p)) {
				sidx[j] = i;
				break;
			}
	}

	for (unsigned int i = k1; i < K; ++i) {

		do {
			same = 0;
			sidx[i] = (unsigned int) ((double) rand() / RAND_MAX * n);
                        for (unsigned int j = 0; j < i; ++j)
                                /* sample without replacement */
                                /* check for same or equal seeds */
                                if (sidx[i] == sidx[j] || (sidx[j] < n &&
                                        !compare_data(x, sidx[i], sidx[j], p))) {
                                        same = 1;
                                        break;
                                }
                } while (same);
                memcpy(seeds[i], x[sidx[i]], p * sizeof **x);
        }

	if (!sd_idx)
		free(sidx);

	return NO_ERROR;
} /* kmodes_init_random_seeds */


/**
 * Initialization of Huang97.
 *
 * @param x	nxp data
 * @param n	number of observations
 * @param p	number of coordinates
 * @param K	number of clusters
 * @param k1	seeds already selected
 * @param wgt	use weighted HD
 * @param rdm	randomize Huang97
 * @param seeds	store chosen seeds here
 * @return	error code
 */
int kmodes_init_h97(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int k1, int wgt, int rdm, data_t **seeds, unsigned int *sd_idx)
{
	unsigned int *sidx = NULL;
	unsigned int r, ncat = 0;
	int same;
	double d, dmin;

	if (k1 > 0)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Initialization "
			"method h97 cannot be used with set seeds.\n");

	/* allocate category counts: last arg force __njc */
	if (!(ncat = allocate_and_compute_category_counts(x, n, p, K, 1)))
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__nj");

	/* allocate space for seed indices */
	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);

		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"seed index");
	} else
		sidx = sd_idx;

	size_t *idx = malloc(ncat * sizeof *idx);
	if (!idx) {
		if (!sd_idx)
			free(sidx);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "idx");
	}

	/* do the sorts [TODO: only necessary if non-random] */
	for (unsigned int j = 0, m = 0; j < p; ++j) {
		for (unsigned int l = 0; l < __nj[j]; ++l)
			idx[m + l] = l;
		if (!rdm && __nj[j] > K)
			with_index_quickselect((void *) __njc[j], &idx[m],
				reverse_compare_uint_elts, 0, __nj[j] - 1, K);
		m += __nj[j];
	}


	/* generate each seed */
	for (unsigned int k = 0; k < K; ++k) {
		do {

			same = 0;
			r = k;

			/* select methodically from most frequent categories */
			for (unsigned int j = 0, m = 0; j < p; ++j) {
				if (rdm)
					r = ((double) rand() / RAND_MAX * __nj[j]);
				seeds[k][j] = idx[m + r];
				m += __nj[j];
				if (!rdm)
					r = r + 2 == K ? 0 : r + 1;
			}

			/* find closest observation to avoid null clusters */
			dmin = INFINITY;
			for (unsigned int i = 0; i < n; ++i) {
				d = hd(x[i], seeds[k], p, wgt);
				if (d < dmin) {
					sidx[k] = i;
					dmin = d;
				}
			}

			/* insure the newly seed is not same as previous seed */
                        for (unsigned int j = 0; j < k; ++j)
                                if (sidx[k] == sidx[j] ||
                                        !compare_data(x, sidx[k], sidx[j], p)) {
                                        same = 1;
                                        break;
                                }
		} while (same);
		memcpy(seeds[k], x[sidx[k]], p * sizeof **x);
	}

	if (!sd_idx && sidx)
		free(sidx);
	free(idx);

	return NO_ERROR;
} /* kmodes_init_h97 */


/**
 * Initialization implemented by accident in Python version of k-modes, an
 * interpretation of Huang97.
 *
 * @param x	nxp data
 * @param n	number of observations
 * @param p	number of coordinates
 * @param K	number of clusters
 * @param k1	number of seed already chosen
 * @param wgt	use weighted Hamming distance
 * @param seeds	store chosen seeds here
 * @return	error code
 */
int kmodes_init_hd17(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int k1, int wgt, data_t **seeds, unsigned int *sd_idx)
{
	unsigned int *sidx = NULL;
	unsigned int r, cfreq, ncat = 0;
	int same;
	double d, dmin;

	if (k1 > 0)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Initialization "
			"method hd17 cannot be used with set seeds.\n");

	/* allocate category counts: last arg force __njc */
	if (!(ncat = allocate_and_compute_category_counts(x, n, p, K, 1)))
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__nj");

	/* allocate space for seed indices */
	if (!sd_idx) {
		sidx = malloc(K * sizeof *sidx);

		if (!sidx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"seed index");
	} else
		sidx = sd_idx;

	/* generate each seed */
	for (unsigned int k = 0; k < K; ++k) {
		do {

			same = 0;

			/* select category by observed abundance */
			for (unsigned int j = 0; j < p; ++j) {
				data_t cat = 0;
				r = (unsigned int) ((double) rand() / RAND_MAX * n);
				cfreq = __njc[j][(int) cat];
				while (r > cfreq) {
					cfreq += __njc[j][(int) ++cat];
				}
				seeds[k][j] = cat;
			}

			/* find closest observation to avoid null clusters */
			dmin = INFINITY;
			for (unsigned int i = 0; i < n; ++i) {
				d = hd(x[i], seeds[k], p, wgt);
				if (d < dmin) {
					sidx[k] = i;
					dmin = d;
				}
			}

			/* insure the new seed is not same as previous seed */
                        for (unsigned int j = 0; j < k; ++j)
                                if (sidx[k] == sidx[j] ||
                                        !compare_data(x, sidx[k], sidx[j], p)) {
                                        same = 1;
                                        break;
                                }
		} while (same);
		memcpy(seeds[k], x[sidx[k]], p * sizeof **x);
	}

	if (!sd_idx)
		free(sidx);

	return NO_ERROR;
} /* kmodes_init_hd17 */


/**
 * Allocate space for density calculations (used by CLB09).  This and
 * allocate_distance are intended to speed up reinitialization by avoiding
 * reallocation of memory.
 *
 * @param x	data (n x p)
 * @param n	number of observations
 * @param p	number of coordinates
 * @return	error status
 */
int allocate_density(data_t **x, unsigned int n, unsigned int p)
{
	MAKE_1ARRAY(__dens, n);
	if (!__dens)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__dens");

	for (unsigned int i = 0; i < n; ++i) {
		__dens[i] = 0;
		for (unsigned int j = 0; j < p; ++j)
			__dens[i] += __njc[j][(int) x[i][j]];
		__dens[i] /= (double) n * p;
	}

	return NO_ERROR;
} /* allocate_density */


/**
 * Allocate space for minimum distance to seed set (used by CLB09).  This and
 * allocate_density are intended to speed up reinitialization by avoiding
 * reallocation of this memory.
 *
 * @param n	number of coordinates
 * @return	error status
 */
int allocate_distance(unsigned int n)
{
	MAKE_1ARRAY(__dis, n);
	return !__dis ? MEMORY_ALLOCATION : NO_ERROR;
} /* allocate_distance */


/**
 * K-modes initialization by Cao et al., 2009.
 *
 * @param x		data (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param k1		number of seeds already chosen
 * @param wgt		weighted Hamming distance?
 * @param rdm		make method random
 * @param seeds		copy of chosen observations (to be determined)
 * @param sd_idx	index of chosen observations (to be determined)
 * @return		error status
 */
int kmodes_init_clb09(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int k1, int wgt, int rdm, data_t **seeds, void *sd_idx, int type)
{
	int err = NO_ERROR;
	size_t ncat = 0;
	unsigned int idx = 0;
	double max = 0, dtmp, sum;
	size_t *sd_idx_s = NULL;
	unsigned int * sd_idx_int = NULL;

	if(type==_SIZE_T_P)
		sd_idx_s = sd_idx;
	else
		sd_idx_int = sd_idx;

	/* allocate category counts: last arg force __njc */
	if (!(ncat = allocate_and_compute_category_counts(x, n, p, K, 1)))
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__nj");

	/* allocate memory for density and minimum distance to seeds */
	if (!__dens && (err = allocate_density(x, n, p)))
		return err;

	if (!__dis && (err = allocate_distance(n)))
		return err;

	for (unsigned int i = 0; i < n; ++i) {
//		fprintf(stderr, "%u: %f\n", i, __dens[i]);
		if (__dens[i] > max) {
			max = __dens[i];
			idx = i;
		}
	}

	if (rdm == 2)
		idx = (unsigned int) ((double) rand() / RAND_MAX * n);

	/* first seed is most dense (or random) */
	if (!k1){
		memcpy(seeds[k1++], x[idx], p * sizeof **x);
		if (sd_idx && (type == _SIZE_T_P))
			sd_idx_s[0] = idx;
		else if(sd_idx && type == UNSINGED_INT)
			sd_idx_int[0] = idx;
	}

	/* compute distance to seed set */
	sum = 0;
	max = 0;
	for (unsigned int i = 0; i < n; ++i) {
		__dis[i] = hd(x[i], seeds[0], p, wgt);
		for (unsigned int j = 1; j < k1; ++j) {
			dtmp = hd(x[i], seeds[j], p, wgt);
			if (dtmp < __dis[i])
				__dis[i] = dtmp;
		}
		if (rdm) {
			sum += (dtmp = __dis[i] * __dens[i]);
			if (max < dtmp) {
				max = dtmp;
				idx = i;
			}
		}
	}

	for (unsigned int k = k1; k < K; ++k) {

		/* find most dense/distant observation... */
		if (rdm) {
			dtmp = ((double) rand() / RAND_MAX);
			idx = 0;
			max = __dens[idx] * __dis[idx];
			while (dtmp > max / sum) {
				++idx;
				max += __dens[idx] * __dis[idx];
			}
		}

		/* ...for next seed */
		memcpy(seeds[k], x[idx], p * sizeof **x);
		if (sd_idx && (type == _SIZE_T_P))
			sd_idx_s[k] = idx;
		else if(sd_idx && type == UNSINGED_INT)
			sd_idx_int[k] = idx;

		/* update distance to seed set */
		sum = 0;
		max = 0;
		for (unsigned int i = 0; i < n; ++i) {
			dtmp = hd(x[i], seeds[k], p, wgt);
			__dis[i] = MIN(__dis[i], dtmp);
			if (rdm) {
				sum += (dtmp = __dens[i] * __dis[i]);
				if (max < dtmp) {
					max = dtmp;
					idx = i;
				}
			}
		}
	}

	return NO_ERROR;

} /* kmodes_init_clb09 */

/**
 * Implement initialization method of Arthur and Vassilvitskii (2007) for
 * k-modes.
 *
 * @param x		data matrix (n x p)
 * @param n		number of observations
 * @param p		number of coordinates
 * @param K		number of clusters
 * @param k1		number of seeds already chosen
 * @param seeds		the chosen observations
 * @param greedy	greedy version
 * @param weight	use weighted Hamming distance
 * @return		error status
 */
int
kmodes_init_av07(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int k1, int weight, data_t **seeds, unsigned int *sd_idx,
	int greedy)
{
	double d, r;
	unsigned int i, j, l, m, s = 0;
	unsigned int log_k = MAX(1, log(K));
	double dmin;
	double dsum = 0;        /* sum of minimum distances to seed set */
	double *W = NULL;       /* min. distances to seed set */

	if (!greedy)
		log_k = 1;

        MAKE_VECTOR(W, n);

        if (!W)
		return KMODES_MEMORY_ERROR;

	if (k1)
		for (l = 0; l < n; ++l) {
			W[l] =  p;
			for (j = 0; j < k1; ++j) {
				d = hd(x[l], seeds[j], p, weight);
				if (d < W[l])
					W[l] = d;
			}
			dsum += W[l];
		}

	for (j = k1; j < K; ++j) {
		/* choose the first seed at random from the data */
		if (j == 0) {
			dmin = INFINITY;
			m = 0;
			do {
                                i = (unsigned int) ((double) rand() / RAND_MAX * n);	/* TODO: allows multiple selections of same seed */
				dsum = 0;
				for (l = 0; l < n; ++l) {
					W[l] = hd(x[l], x[i], p, weight);
					dsum += W[l];
				}
				if (dsum < dmin) {
					dmin = dsum;
					s = i;
				} else
					s = i;
				++m;
			} while (m < log_k);

			if (sd_idx)
				sd_idx[0] = s;
			memcpy(seeds[0], x[s], p * sizeof **seeds);
			dsum = dmin;
		} else {
			/* W/sum W is the cdf for choosing next seed */
			dmin = INFINITY;
			m = 0;
			do {
				r = dsum * ((double) rand() / RAND_MAX);
				for (i = 0, dsum = 0;
					(i < n) && (dsum < r);
					dsum += W[i++]);
				if (i)	/* rare case: r == 0 */
					i--;
				if (greedy) {
					dsum = 0;       /* wss */
					for (l = 0; l < n; l++)
						dsum += hd_min(x[l], x[i], p,
							weight, W[l]);
					if (dsum < dmin) {
						dmin = dsum;
						s = i;
					}
				} else
					s = i;
				m++;
			} while (m < log_k);

			if (sd_idx)
				sd_idx[j] = s;

			memcpy(seeds[j], x[s], p * sizeof **seeds);

			/* reset minimum distance, maybe */
			dsum = 0;
			for (i = 0; i < n; i++) {
				W[i] = hd_min(x[i], seeds[j], p, weight, W[i]);
				dsum += W[i];
			}
		}
	}

	FREE_VECTOR(W);

	return NO_ERROR;
} /* kmodes_init_av07 */


/**
 * Compare two observations to allow ordering and testing for equality.
 *
 * @param x	data (? x n)
 * @param i	index of first observation
 * @param j	index of second observation
 * @param n	number of coordinates
 * @return	-1, 0, 1 order indicator
 */
int compare_data(data_t **x, unsigned int i, unsigned int j, unsigned int n) {
	for (unsigned int l = 0; l < n; ++l) {
		if (x[i][l] > x[j][l])
			return 1;
		else if (x[i][l] < x[j][l])
			return -1;
	}
	return 0;
} /* compare_data */

int compare_data_to_seed(data_t **x, unsigned int i, data_t *seed, unsigned int p) {
	for (unsigned int j = 0; j < p; ++j) {
		if (x[i][j] > seed[j])
			return 1;
		else if (x[i][j] < seed[j])
			return -1;
	}
	return 0;
} /* compare_data_to_seed */


/**
 * Utility function allocates and computes category counts.
 *
 * @param x	data (n x p)
 * @param n	number of observations
 * @param p	number of coordinates
 * @param K	number of clusters
 * @param wgt	weighted?
 * @return	total number of categories across coordinates
 */
size_t allocate_and_compute_category_counts(data_t **x, unsigned int n, unsigned int p,
	unsigned int K, int wgt)
{
	size_t l = 0;
	unsigned int j;

	/* pre-compute number of categories per column */
	if (!__nj) {
		l = allocate_and_compute_nj(x, n, p);

		/* allocate space for counts of every category */
		if (l && (wgt || K == 1) && !__njc) {
			allocate_and_compute_njc(x, n, p, l);
			return l;
		}
	/* will need total number of categories, summed over coordinates */
	} else
		for (j = 0; j < p; ++j)
			l += __nj[j];

	return l;
} /* allocate_and_compute_category_counts */


/**
 * Allocate and compute number of categories at each coordinate.
 *
 * @param x	data (n x p)
 * @param n	number of observations
 * @param p	number of coordinates
 * @return	total number of categories summed across coordinates.
 */
size_t allocate_and_compute_nj(data_t **x, unsigned int n, unsigned int p) {
	size_t ncat = 0;
	__nj = calloc(p, sizeof *__nj);	/* MAKE_VECTOR(__nj, p); */
	if (!__nj)
		return ncat;
	/* find no. categories in each column */
	for (unsigned int j = 0; j < p; ++j) {
		for (unsigned int i = 0; i < n; ++i)
			if (x[i][j] > __nj[j]) __nj[j] = x[i][j];
		++__nj[j];	/* assumes 0 is one of the values */
		ncat += __nj[j];
	}
	return ncat;
} /* allocate_and_compute_nj */


/**
 * Allocate and compute number of each category per coordinate.  These data
 * are needed for weighted Hamming distance (see hd()).
 *
 * @param x	data (n x p)
 * @param n	number of observations
 * @param p	number of coordinates
 * @param ncat	total number of categories across coordinates
 * @return	error status
 */
int allocate_and_compute_njc(data_t **x, unsigned int n, unsigned int p, size_t ncat) {
	if (!__nj)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj needs to be "
			"allocated before allocating __njc");
	size_t *tmp = malloc(ncat * sizeof **__njc);

	if (!tmp)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "tmp");

	__njc = malloc(p * sizeof *__njc);

	if (!__njc) {
		free(tmp);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__njc");
	}

	for (unsigned int j = 0; j < p; ++j) {
		__njc[j] = tmp;
		tmp += __nj[j];
		for (unsigned int l = 0; l < __nj[j]; ++l)
			__njc[j][l] = 0;
		for (unsigned int i = 0; i < n; ++i)
			++__njc[j][(int) x[i][j]];
		for (unsigned int l = 0; l < __nj[j]; ++l) {
		}
	}

	return NO_ERROR;
} /* allocate_and_compute_njc */


/**
 * Set up or reset __nkjc.
 *
 * @param K	number of clusters
 * @param p	number of coordinates
 * @param l	total number of categories
 * @return	error status
 */
int set_nkjc(unsigned int K, unsigned int p, size_t l)
{
	if (!__nkjc && allocate_nkjc(K, p, l))
		return KMODES_MEMORY_ERROR;
	else	/* QUITE DANGEROUS: user must call reset_k()
		 * if change K or free_kmodes() if change data
		 * between calls to kmodes_*() algorithms. */
		return reset_nkjc(__nkjc, K, p);
} /* set_nkjc */


/**
 * Allocate memory for __nkjc.
 *
 * @param K	number of clusters
 * @param p	number of coordinates
 * @param ncat	total number of categories, summed across coordinates
 * @return	error status
 */
int allocate_nkjc(unsigned int K, unsigned int p, size_t ncat)
{
	if (!__nj)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj needs to be "
			"allocated before allocating __njc");

	size_t *tmp = calloc(ncat * K, sizeof ***__nkjc);

	if (!tmp)
		goto ALLOCATE_NKJC_ERROR;

	CMAKE_2ARRAY(__nkjc, K, p);

	if (!__nkjc)
		goto ALLOCATE_NKJC_ERROR;

	for (unsigned int k = 0; k < K; ++k)
		for (unsigned int j = 0; j < p; ++j) {
			__nkjc[k][j] = tmp;
			tmp += __nj[j];
		}

	return NO_ERROR;

ALLOCATE_NKJC_ERROR:
	if (tmp) free(tmp);
	FREE_2ARRAY(*__nkjc);

	return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__nkjc");

} /* allocate_nkjc */


/**
 * Reset category counts per cluster per coordinate to 0.
 *
 * @param nkjc	count of categories at each coordinate in each cluster
 * @param K	number of clusters
 * @param p	number of coordinates
 * @return	error status
 */
int reset_nkjc(size_t ***nkjc, unsigned int K, unsigned int p)
{
	if (!nkjc || !__nj)
		return KMODES_INTERNAL_ERROR;
	/* initialize category counts */
	for (unsigned int k = 0; k < K; k++)
		for (unsigned int j = 0; j < p; j++)
			for (unsigned int l = 0; l < __nj[j]; l++)
				nkjc[k][j][l] = 0;
	return NO_ERROR;
} /* reset_nkjc */


/**
 * Compute category counts per cluster per site.  This is mostly used for debugging
 * but it is also used to initialize from a partition.
 *
 * @param nkjc	counts of categories per cluster per site (K x p x __nj)
 * @param x	data (n x p)
 * @param ic	cluster assignments (1 x n)
 * @param n	number of observations
 * @param p	number of coordinates
 */
void compute_nkjc(size_t ***nkjc, data_t **x, unsigned int *ic, unsigned int n, unsigned int p)
{
	if (!nkjc)
		return;
	for (unsigned int i = 0; i < n; ++i)
		for (unsigned int j = 0; j < p; ++j)
			nkjc[ic[i]][j][(int) x[i][j]]++;
} /* compute_nkjc */


/**
 * Returns human-friendly description of a k-means error.
 *
 * @param err integer representation of error; something returned by
 *            #kmeans() or some k-means initialiation routine.
 * @return descriptive string explaining error or NULL if invalid error
 */
const char *kmodes_error(int err) {
	if (err == KMODES_EXCEED_ITER_WARNING)
		return "Exceeded maximum iterations";
	else if (err == KMODES_NULL_CLUSTER_ERROR)
		return "Inferred null cluster";
	else if (err == KMODES_CALLER_INPUT_ERROR)
		return "Requested K=0 or K>n, more clusters than observations";
	else if (err == KMODES_MEMORY_ERROR)
		return "Memory allocation error";
	else if (err == KMODES_INVALID_INITIALIZATION_METHOD)
		return "Unknown initialization method requested";
	else if (err == KMODES_INTERNAL_ERROR)
		return "Internal error";
	else
		return NULL;
} /* kmodes_error */

/**
 * Returns human-friendly description of k-modes initializatiom method.
 *
 * @param method	method
 * @return		string constant
 */
const char *kmodes_init_method(int method)
{
	if (method == KMODES_INIT_USER_SEEDS)
		return "user seeds";
	else if (method == KMODES_INIT_RANDOM_SEEDS)
		return "random seeds";
	else if (method == KMODES_INIT_H97)
		return "Huang, 1997";
	else if (method == KMODES_INIT_H97_RANDOM)
		return "randomization of Huang, 1997";
	else if (method == KMODES_INIT_HD17)
		return "Huang, 1997 as interpretted by Python author de Vos";
	else if (method == KMODES_INIT_CLB09)
		return "Cao et al., 2009";
	else if (method == KMODES_INIT_CLB09_RANDOM)
		return "randomization of Cao et al., 2009";
	else if (method == KMODES_INIT_AV07)
		return "k-means++ adapted for k-modes";
	else if (method == KMODES_INIT_AV07_GREEDY)
		return "greedy k-means++ adapted for k-modes";
	else if (method == KMODES_INIT_RANDOM_FROM_PARTITION)
		return "random initialization from true partition";
	else if (method == KMODES_INIT_RANDOM_FROM_SET)
		return "random initialization from set of seeds";
	else
		return NULL;
}/* kmodes_init_method */

const char *kmodes_algorithm(int algorithm)
{
	if (algorithm == KMODES_HUANG)
		return "Huang 1997";
	else if (algorithm == KMODES_LLOYD)
		return "Lloyd 1982";
	else if (algorithm == KMODES_HARTIGAN_WONG)
		return "Hartigan Wong 1979";
	else
		return NULL;
} /* kmodes_algorithm */


/**
 * Free static memory if changing data.
 */
void free_kmodes()
{
	if (__nj) {
		free(__nj);
		__nj = NULL;
	}
	if (__njc) {
		free(__njc[0]);
		free(__njc);
		__njc = NULL;
	}
	if (__nkjc) {
		free(__nkjc[0][0]);
		FREE_2ARRAY(__nkjc);
	}
	if (__dens) {
		free(__dens);
		__dens = NULL;
	}
	if (__dis) {
		free(__dis);
		__dis = NULL;
	}
	if (__criterion) {
		free(__criterion);
		__criterion = NULL;
	}

	if (__ic2) {
		free(__ic2);
		__ic2 = NULL;
	}

	if (__live) {
		free(__live);
		__live = NULL;
	}
	if (__cd) {
		free(__cd);
		__cd = NULL;
	}
	if (__c2) {
		free(__c2[0]);
		free(__c2);
		__c2 = NULL;
	}

#ifdef __KMODES_DEBUGGING__
	if (__nkjc_debug) {
		free(__nkjc_debug[0][0]);
		FREE_2ARRAY(__nkjc_debug);
	}
#endif
} /* free_kmodes */


/**
 * Free static memory if changing number of clusters.
 */
void reset_k()
{
	if (__nkjc) {
		free(__nkjc[0][0]);
		FREE_2ARRAY(__nkjc);
	}

	if (__criterion) {
		free(__criterion);
		__criterion = NULL;
	}
	if (__live) {
		free(__live);
		__live = NULL;
	}
	if (__cd) {
		free(__cd);
		__cd = NULL;
	}
	if (__c2) {
		free(__c2[0]);
		free(__c2);
		__c2 = NULL;
	}


#ifdef __KMODES_DEBUGGING__
	if (__nkjc_debug) {
		free(__nkjc_debug[0][0]);
		FREE_2ARRAY(__nkjc_debug);
	}
#endif
} /* reset_k */



#ifdef __KMODES_DEBUGGING__


/**
 * Allocate memory for __nkjc_debug used for debugging.
 *
 * @param K	number of clusters
 * @param p	number of coordinates
 * @param ncat	total number of categories, summed across coordinates
 * @return	error status
 */
int allocate_nkjc_debug(unsigned int K, unsigned int p, size_t ncat)
{
	if (!__nj)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj needs to be "
			"allocated before allocating __njc_copy");

	size_t *tmp = malloc(ncat * K * sizeof ***__nkjc_debug);

	if (!tmp)
		goto ALLOCATE_NKJC_ERROR;

	CMAKE_2ARRAY(__nkjc_debug, K, p);

	if (!__nkjc_debug)
		goto ALLOCATE_NKJC_ERROR;

	for (unsigned int k = 0; k < K; ++k)
		for (unsigned int j = 0; j < p; ++j) {
			__nkjc_debug[k][j] = tmp;
			tmp += __nj[j];
		}

	return NO_ERROR;

ALLOCATE_NKJC_ERROR:
	if (tmp) free(tmp);
	FREE_2ARRAY(__nkjc_debug);

	return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "__nkjc_debug");

} /* allocate_nkjc_debug */


/**
 * A debugging function to detect when there is a possible move.
 *
 * @param x	data (n x p)
 * @param n	number of observations
 * @param p	number of coordinates
 * @param K	number of centers
 * @param ic	cluster assignments
 * @param i	index of observation to move
 * @param l1	current cluster of observation i
 * @param l2	proposed cluster for observation i
 * @param wgt	use weighted Hamming distance?
 * @return	cost difference: new cost - old cost
 */
double debug_cost_to_move(data_t **x, unsigned int n, unsigned int p, unsigned int K,
	unsigned int *ic, unsigned int i, unsigned int l1, unsigned int l2, int wgt)
{
	double curr_cost, new_cost;

	reset_nkjc(__nkjc_debug, K, p);
	compute_nkjc(__nkjc_debug, x, ic, n, p);
	curr_cost = compute_cost_from_nkjc(__nkjc_debug, K, p, wgt);
	for (unsigned int j = 0; j < p; ++j) {
		__nkjc_debug[l1][j][(int) x[i][j]]--;
		__nkjc_debug[l2][j][(int) x[i][j]]++;
	}
	new_cost = compute_cost_from_nkjc(__nkjc_debug, K, p, wgt);

	return new_cost - curr_cost;
} /* debug_cost_to_move */


/**
 * Debugging function for computing the cost from sufficient statistics.
 *
 * @param nkjc	category counts per cluster per site (k x p x __nj)
 * @param K	number of clusters
 * @param p	number of coordinates
 * @param wgt	use weighted Hamming distance in costs
 * @return	cost
 */
static inline double compute_cost_from_nkjc(size_t ***nkjc, unsigned int K,
	unsigned int p, int wgt)
{
	double max;
	double cost = 0;
	unsigned int max_l;

	if (!__nj || (wgt && !__njc))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "__nj = %p; __njc = "
			"%p!\n", (void *)__nj, (void *)__njc));

	for (unsigned int k = 0; k < K; ++k) {
		for (unsigned int j = 0; j < p; ++j) {

			/* find mode */
			max = 0;
			max_l = 0;
			for (unsigned int l = 0; l < __nj[j]; ++l)
				if (max < nkjc[k][j][l]) {
					max = nkjc[k][j][l];
					max_l = l;
				}


			/* compute cost */
			for (unsigned int l = 0; l < __nj[j]; ++l) {
				if (l == max_l)
					continue;
				cost += nkjc[k][j][l] * (wgt ? (double)
					(__njc[j][(int) max_l] + __njc[j][(int) l])
					/ (__njc[j][(int) max_l] * __njc[j][(int) l])
					: 1.);
			}
		}
	}

	return cost;
} /* compute_cost_from_nkjc */


/**
 * Debugging function to verify modes and minor modes.  Simply aborts and exists
 * if there is a problem.
 *
 * @param nkjc		category counts per cluster per site (K x p x __nj)
 * @param modes 	modes to verify (K x p)
 * @param mmodes	minor modes to verify (K x p)
 * @param K		number of clusters
 * @param p		number of coordinates
 */
void compute_and_compare_modes(size_t ***nkjc, data_t **modes, data_t **mmodes,
	unsigned int K, unsigned int p)
{
	if (!nkjc || !__nj)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Not here!\n"));

	for (unsigned int k = 0; k < K; ++k)
		for (unsigned int j = 0; j < p; ++j) {
			size_t max_nkjc = 0, max2_nkjc = 0;
			unsigned int max_l = 0, max2_l = 0;
			for (unsigned int l = 0; l < __nj[j]; ++l)
				if (nkjc[k][j][l] > max_nkjc) {
					max2_nkjc = max_nkjc;
					max_nkjc = nkjc[k][j][l];
					max2_l = max_l;
					max_l = l;
				} else if (mmodes &&
					nkjc[k][j][l] > max2_nkjc) {
					max2_nkjc = nkjc[k][j][l];
					max2_l = l;
				}
			if (max_l != modes[k][j]) {
				summarize_clusters("compute_and_compare_modes",
					nkjc, modes, mmodes, K, p);
				exit(mmessage(DEBUG_MSG, INTERNAL_ERROR,
					"Cluster %u, coordinate %zu: mode is %u"
					" and should be %u\n", k, j, modes[k][j],
					max_l));
			}
			if (mmodes && max2_l != mmodes[k][j]) {
				summarize_clusters("compute_and_compare_modes",
					nkjc, modes, mmodes, K, p);
				exit(mmessage(DEBUG_MSG, INTERNAL_ERROR,
					"Cluster %u, cooordinate %zu: minor "
					"mode is %u and should be %u\n", k, j,
					mmodes[k][j], max2_l));
			}
		}
} /* compute_and_compare_modes */


/**
 * Debugging function compares to sets of category counts.  Simply exists if
 * there is no match.
 *
 * @param nkjc1	category counts per cluster per site (K x p x __nj)
 * @param nkjc2	category counts per cluster per site (K x p x __nj)
 * @param K	number of clusters
 * @param p	number of coordinates
 */
void compare_nkjc(size_t ***nkjc1, size_t ***nkjc2, unsigned int K, unsigned int p)
{
	if (!nkjc1 || !nkjc2 || !__nj)
		return;

	for (unsigned int k = 0; k < K; k++)
		for (unsigned int j = 0; j < p; j++)
			for (unsigned int l = 0; l < __nj[j]; l++)
				if (nkjc1[k][j][l] != nkjc2[k][j][l])
					exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
						"nkjc[%u][%zu][%zu]: %zu != "
						"%zu\n", k, j, l, nkjc1[k][j][l],
						nkjc2[k][j][l]));
} /* compare_nkjc */


/**
 * A debugging function that displays category counts per cluster per site.
 * It also displays the mode and minor mode for each cluster.
 *
 * @param str	user-friendly string to display
 * @param nkjc	counts of characters per cluster per position
 * @param c	modes
 * @param c2	minor modes
 * @param K	number of modes
 * @param p	number of coordinates in observations
 */
void summarize_clusters(char const *str, size_t ***nkjc, data_t **c,
	data_t **c2, unsigned int K, unsigned int p)
{
	if (!nkjc)
		return;

	debug_msg(0, DEBUG_I, "%s\n", str);
	for (unsigned int k = 0; k < K; ++k) {
		fprintf(stderr, "Cluster %u:", k);
		if (c)
			fprint_data_ts(stderr, c[k], p, 1, 0);
		if (c2) {
			fprintf(stderr, "; ");
			fprint_data_ts(stderr, c2[k], p, 1, 1);
		} else {
			fprintf(stderr, "\n");
		}
		for (unsigned int j = 0; j < p; ++j) {
			fprintf(stderr, "Coordinate %u:", j);
			for (unsigned int l = 0; l < __nj[j]; ++l)
				fprintf(stderr, " %zu", nkjc[k][j][l]);
			fprintf(stderr, "\n");
		}
	}
} /* summarize_clusters */

#endif
