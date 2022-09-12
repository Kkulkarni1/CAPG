/**
 * @file ampliclust.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 * Cluster amplicon reads.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#include "ampliclust.h"
#include "initialize.h"
#include "aecm.h"
#include "kmodes.h"
#include "cluster.h"
#include "statistics.h"
#include "math.h"
#include "io.h"
#include "hash.h"
#include "simulate.h"
#include "options.h"
#include "initialize_options.h"
#include "error.h"


/**
 * Cluster amplicon reads.
 *
 * @param dat	data object
 * @param mod	model object
 * @param ini	initializer object
 * @param mopt	model options
 * @param iopt	initializer options
 * @param cpi	function to call after each iteration
 * @param init	function used to initialize
 * @param obj	additional content needed by cpi and init
 * @return	error status
 */
int cluster_amplicons(data *dat, model *mod, initializer *ini, 
	model_options *mopt, initialize_options *iopt,
	call_per_iterate_func cpi, initialize_func init, void *obj)
{
	int err = NO_ERROR;	/* error status */

	if (!iopt->n_init)
		if ((err = init(dat, mod, ini, obj)))
			return err;

	/* repeated initializations of AECM */
	for (ini->n_inits = 0; ini->n_inits < iopt->n_init; ++ini->n_inits) {

		/* fast forward for debugging purposes
		if (ini->n_inits && ini->n_inits < 19)
			opt->inio->n_kmodes_iter = 0;
		else
			opt->inio->n_kmodes_iter = 10000;
		*/

		/* initialize parameters: AECM will start with E step */
		if ((err = init(dat, mod, ini, obj)))
			return err;

		/* fast forward for debugging
		if (ini->n_inits && ini->n_inits < 19) {
			cc_msg_cont(global_wp, DEBUG_I, DEBUG_I, "\n");
			continue;
		}
		*/

		/* AECM iterations until convergence */
		aecm(dat, mod, mopt);

		cpi(dat, mod, ini, obj);

	}/* initialization */

	return NO_ERROR;
} /* cluster_amplicons */


/**
 * Hard assignment of reads to clusters based on posterior probabilities.
 *
 * @param eik		matrix E (contains posterior probabilities )
 * @param K 		number of clusters (K or n_mix)
 * @param sample_size	number of reads
 * @param cs		cluster size (\par mod.n_mix x 1)
 * @param ci		cluster index (\par dat.sample_size x 1)
 * @param by_K		matrix E in K*n (0) or in n*K (1)
 * @return		error status
 */
int assign_clusters(double *eik, unsigned int K, size_t sample_size,
				unsigned int *cs, unsigned int *ci, int by_K)
{
	double max, tmp;
	unsigned int i, k, l = 0;

	for (k = 0; k < K; ++k)
		cs[k] = 0;

	for (i = 0; i < sample_size; ++i) {
		max = -INFINITY;
		for (k = 0; k < K; ++k) {
			tmp = by_K ? eik[sample_size*k + i] : eik[K*i + k];
//if (i == 24769) fprintf(stderr, " %f", tmp);
			if (tmp > max) {
				max = tmp;
				l = k;
			}
		}
//if (i == 24769) fprintf(stderr, " <= sequence %u\n", i);
		/* three places to store result */
		ci[i] = l;
		++cs[l];
	}

	return NO_ERROR;
} /* assign_clusters */
