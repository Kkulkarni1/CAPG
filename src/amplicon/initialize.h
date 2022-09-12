/**
 * @file initialize.h
 * @author Karin S. Dorman
 * @author Xiyu Peng
 *
 * Header file for initialization structs, defines and extern functions
 */

#ifndef __H_INITIALIZE__
#define __H_INITIALIZE__

#include "model.h"
#include "data.h"
#include "initialize_options.h"

/**
 * Initialization methods.  See kmodes.h for k-modes initialization methods.
 */
enum {
	INIT_PARTITION,			/*<! init. w/ partition */
	INIT_HAPLOTYPES,		/*<! init. w/ haplotypes */
	INIT_HAPLOTYPES_AND_PARTITION,	/*<! init. w/ both */
	INIT_SEEDS,			/*<! init. w/ selected observations */
	INIT_PARAMETERS,		/*<! init. w/ parameters */
	INIT_RANDOM_SEEDS,		/*<! init. by random seeds */
	INIT_KMODES,			/*<! init. by k-modes */
	INIT_HW_KMODES,			/*<! init. by Hartigan & Wong k-modes */
	INIT_CAO,               	/*<! init. by Cao's method */
	INIT_DEBLUR,			/*<! init. by deblur */
	INIT_IMDEBLUR,		          /*<!  init. by imdeblur (given K) */
	INIT_TRUE_VALUES,		/*<! init. with true parameter values */
	INIT_TRUE_PARTITION,		/*<! init. with true partition */
	NUM_INIT_METHODS
};

/**
 * Run methods.
 */
enum {
	EXHAUSTIVE_EM,	/*<! run EM to convergence on each initialization */
	/* the rest do multiple random initializations and run EM on best */
	INI_EM,		/*<! use criterion-based initializer */
	RND_EM,		/*<! use ll-evaluated random initializations */
	EM_EM,		/*<! use ll after n EM iterations */
	NUM_RUN_METHODS
};

/**
 * impro_deplur errors.
 */
enum {
	DEBLUR_NO_ERROR,		/*<! no errors */
	DEBLUR_EXCEED_MAX_ITERATIONS,	/*<! over maximum allowed iterations */
	DEBLUR_ASCENT_VIOLATION,	/*<! true abundance increase or under 0 */
	NUM_DEBLUR_ERRORS		/*<! number of deblur errors */
};

typedef struct _initializer initializer;

struct _initializer {
	int synced;			/*<! synced with data */
	unsigned int n_inits;		/*<! current initialization */
	unsigned int K;			/*<! number of clusters */

	/* current initialization */
	size_t *seed_idx;		/*<! seed indices */
	data_t **seeds;			/*<! access to seed data */
	unsigned int *seed_lengths;	/*<! seed lengths */
	unsigned int *cluster_id;	/*<! cluster assignments (k-modes) */
	double *criterion;		/*<! criterion per cluster (k-modes) */
	unsigned int *cluster_size;	/*<! cluster sizes (k-modes) */

	/* best initialization for rnd-em & ilk */
	double best_total;		/*<! current best criterion */
	size_t *best_seed_idx;		/*<! seed indices */
	data_t **best_modes;		/*<! estimated modes (k-modes) */
	double *best_criterion;		/*<! criterion per cluster (k-modes) */
	unsigned int *best_cluster_id;	/*<! cluster assignments (k-modes) */
	unsigned int *best_cluster_size;/*<! cluster sizes (k-modes) */

	/* sequence table with abundance (sorted) for Deblur and imdeblur */
	size_t *uniq_seq_idx;		/*<! unique sequences index */
	unsigned int *uniq_seq_count;	/*<! observed abundance of uniq seqs */
	double *abun_true;		/*<! estimated true abundance */
	double *p;			/*<! probability of being chosen */
	unsigned int *H;		/*<! idx of haplotypes in uniq seq */
	double *e_trans;		/*<! (K+1)*n expected number misreads */

	/* optimal initialization criterion across repeated initialization */
	double optimal_total;		/*<! optimal criterion across inits */

}; /* initializer */


int make_initializer(initializer **ini, data *dat, model *mod, initialize_options *opt, int tbd);
int sync_initializer(initializer *ini, data *dat, initialize_options *opt);
int realloc_initializer(initializer *ini, data *dat, initialize_options *opt);
int initialize(model *mod, data *dat, initializer *ini, initialize_options *opt);
int simple_initialize(data *dat, model *mod, initializer *ini, void *obj);
int read_initialization_file(char const * const filename, data *dat, initialize_options *opt, initializer *ini);
int pull_bootstraps(data *dat, initialize_options *opt, initializer *ini);
int initialize_model_parameters(model *mod, data *dat, initialize_options *opt, initializer *ini, int best);
void free_initializer(initializer *ini, initialize_options *opt);

/* [KSD, TODO] Move these to deblur.{hc} */
int deblur(data *dat, size_t *idx, unsigned int *count, double error_rate, unsigned int K,data_t **seeds,size_t *seed_idx, unsigned int* seed_lengths, double *p,double *abun_true,unsigned int *H);
int impro_deblur(initialize_options *opt,data *dat, initializer *ini, model *mod);
int cal_e_trans(data *dat, initialize_options *opt, double *e_trans, unsigned int select, double *error_profile, unsigned char *seq, int self);
int detect_false_positive(data *dat, model *mod, initializer *ini, double *pre_bic, double *pre_aic,unsigned int select,int use_aic, int error_profile, int final);
//int impro_deblur(initialize_options *opt,data *dat, size_t *idx, unsigned int *count, unsigned int K,
//	data_t **seeds,size_t *seed_idx, unsigned int* seed_lengths,int stoch, int conve);
int e_TrueAbun(initialize_options *opt, data *dat, double *e_trans, double *abun_true, unsigned int *H, size_t *idx, unsigned int count_i,unsigned int select, unsigned int K,unsigned int i, int conve);

#endif
