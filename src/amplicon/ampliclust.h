/**
 * @file ampliclust.h
 * @author Karin S. Dorman
 *
 * Header file for ampliclust.
 */

#ifndef __H_AMPLICLUST__
#define __H_AMPLICLUST__

#include <stdlib.h>
#include <curses.h>

#include "data.h"
#include "model.h"
#include "initialize.h"
#include "model_options.h"
#include "initialize_options.h"


/**
 * Where to cache results.  Is this the current solution, the best local
 * solution, the optimal global solution, or the truth?
 */
enum {
	DEFAULT_VALUES,	/*<! cluster_id, criterion, etc. */
	BEST_VALUES,	/*<! best_cluster_id, best_criterion, etc. */
	OPTIMAL_VALUES,	/*<! optimal_cluster_id, etc. */
	TRUE_VALUES,	/*<! true parameter values */
	FROM_BEST,	/*<! set default from best parameters */
	FROM_TRUE	/*<! set default from true parameters */
};

/**
 * Function pointer.  For calling during each iteration of AECM.
 */
typedef void (*call_per_iterate_func)(data *, model *, initializer *, void *);
typedef int (*initialize_func)(data *, model *, initializer *, void *);


int cluster_amplicons(data *dat, model *mod, initializer *ini, model_options *mopt, initialize_options *iopt, call_per_iterate_func cpi, initialize_func init, void *obj);
int assign_clusters(double *eik, unsigned int K, size_t sample_size, unsigned int *cs, unsigned int *ci,int by_K);

#endif
