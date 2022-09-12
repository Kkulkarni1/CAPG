#ifndef __INITIALIZE_OPTIONS_H__
#define __INITIALIZE_OPTIONS_H__

#include "model_options.h"
#include "kmodes.h"

typedef struct _initialize_options initialize_options;

struct _initialize_options {
	/* since the initializer is meant to initialize model parameters, it
	 * needs to have access to the model object it needs to initialize
	 */
	model_options *modo;
	int estimate_K;			/*<! initializer can estimate K */

	int initialization_method;	/*<! see methods in initialize.h */
	unsigned int n_init;		/*<! number of random initializations */
	int run_method;			/*<! exhaustive, rnd-em or em-em */
	unsigned int n_inner_init;	/*<! number of initializer initializations */
	int kmodes_initialization_method;	/*<! see kmodes.h */
	kmodes_options *kmodes_opt;	/*<! kmodes options object */
	unsigned int n_kmodes_iter;	/*<! max. number of k-modes iterations */
	int kmodes_huang97;		/*<! use Huang97 first round in k-modes */
	char *initialization_file; 	/*<! name of file for initialization */
	char const *partition_file;	/*<! name of partition file */
	int assignment_method;		/*<! way to assign reads to partitions */
	double beta_epsilon;		/*<! less stringent for initialization */

	/* initialization: deblur and imdeblur */
	int evaluate_initialization;	/*<! compare true seeds vs hash table (output) */
	double error_rate;		/*<! mean read error rate per nucleotide */
	int stochastic_deplur;		/*<! stochastic version */
	int convergence_deplur;         /*<! use fixed point iteration */
	double error_adjust;		/*<! adjust the rate of error */
	int check_false_positive;	/*<! discard unsupported haplotypes */
	int remove_gamma;		/*<! remove gamma (no longer prob. model) */
	int use_error_profile;		/*<! use provided error profile */
	double low_bound;		/*<! abundance threshold */
};

int make_initialize_options(initialize_options **, model_options *);
int process_initialization_option(int argc, const char **argv, int i, void *opt, initialize_options *inio);

#endif
