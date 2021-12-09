/**
 * @file run_ampliclust.h
 * @author Karin S. Dorman
 *
 * Header file for running vanilla ampliclust.
 */

#ifndef __H_RUN_AMPLICLUST__
#define __H_RUN_AMPLICLUST__

#include <stdio.h>
#ifdef USE_CURSES
#include <curses.h>
#endif

#include "simulate.h"
#include "options.h"

#ifdef USE_CURSES
extern WINDOW *wp;
#endif
extern FILE *active_fp;

typedef struct _run_info run_info;

/**
 * Store information about the best solution across multiple runs (varying 
 * initialization or K or whatever).
 */
struct _run_info {
	size_t *optimal_seed_idx;		/*<! seed yielding best soln */
	unsigned int *optimal_cluster_id;	/*<! optimal hard clustering */
	unsigned int *optimal_cluster_size;	/*<! optimal cluster sizes */

	simulator *sim;				/*<! simulator object if simulated */
	options *opt;				/*<! ampliclust options object */
};

typedef struct _ampliclust_obj ampliclust_obj;

int do_initialize(data *dat, model *mod, initializer *ini, void *obj);
void do_per_iterate(data *dat, model *mod, initializer *ini, void *obj);
int ampliclustK(data *dat, model *best_mod, options *opt, initializer *ini, run_info *ri);

int realloc_run_info(run_info *ri, size_t sample_size, unsigned int K);
void free_run_info(run_info *ri);


#endif
