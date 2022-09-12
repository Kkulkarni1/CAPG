/**
 * @file haplotype_option.h
 * @author Yudi Zhang
 */

#ifndef HAPLOTYPE_OPTION_H
#define HAPLOTYPE_OPTION_H

#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>

#include "constants.h"
#include "kmodes.h"

/**
 * Varieties of Lloyd's algorithm for fastq data.
 */
enum {
	FASTQ_LLOYDS,		/*<! naive Lloyd's algorithm */
	FASTQ_LLOYDS_EFFICIENT,	/*<! efficient Lloyd's algorithm */
	LLOYDS_OLD,		/*<! my Lloyd's?? */
	FASTQ_MACQUEEN,		/*<! MacQueen's algorithm */
	FASTQ_HW		/*<! Hartigan & Wong algorithm */
};

/**
 * Run options for k-haplotype.
 */
typedef struct _options options;

struct _options {
	/* model */
	unsigned int K;			/*<! number of clusters */
	
	/* data */
	int run_with_quals;		/*<! input file format: fastq */
	char const *datafile;		/*<! name of datafile */
	double n_effective_coordinates;	/*<! effective no. indep. coordinates */
	int subtract_one;		/*<! subtract one (to make 0-base) */

	/* run conditions */
	unsigned long seed;		/*<! random number seed [srand()] */
	kmodes_options *kopt;		/*<! kmodes run options passed to algorithms */
	unsigned int n_init;		/*<! no. random initializations */
	unsigned int n_inner_init;	/*<! no. random inner loop initializations */
	unsigned int n_max_iter;	/*<! max. number of iterations */
	int kmodes_algorithm;		/*<! run Lloyd's, Huang's or HW k-modes */
	int use_hartigan;		/*<! use Hartigan updates (Lloyd's | Huang's) */
	int use_qtran;			/*<! use quick-transfer phase in HW */
	int update_modes;		/*<! update modes */
	int weight;			/*<! weighting method (not implemented) */
	double seconds;			/*<! number of seconds to initialize */
	int shuffle;			/*<! shuffle input order */
	int continue_run;		/*<! continue previous run */
	double target;			/*<! target optimum: best from prev. run */
	unsigned int n_bootstrap;

	/* model */
	int estimate_k;			/*<! estimate k */
	unsigned int min_k;		/*<! minimum K */
	unsigned int max_k;		/*<! maximum K */
	unsigned int n_k;		/*<! number of K [min_k, max_k] */

	/* results */
	char const ***result_files;
	unsigned int *n_result_files;
	
	/* initialization */
	int init_method;		/*<! initialization method to use */
	unsigned int *seed_idx;		/*<! user-provided seed indices */
	char const *pfile;		/*<! partition file */
	char const *sfile;		/*<! seed file */
	char const *mfile;		/*<! mode file */
	char const *mfile_out;		/*<! name of mode file to output */
	unsigned int n_sd_idx;		/*<! tmp: length of seed_idx */
	data_t **seed_set;		/*<! available seeds */
	unsigned int n_seed_set;	/*<! number of seeds in seed set */
	
	/* output */
	char const *data_outfile;	/*<! name of datafile to write */
	char const *ini_file;		/*<! initialization data outfile */
	char const *soln_file;		/*<! solution file */
	int info;			/*<! level of information to output */
	int quiet;			/*<! be quiet */
	
	/* simulation: this stuff is in model (ampliclust) */
	int simulate;			/*<! request to simulate data */
	unsigned int sim_K;		/*<! number of simulated clusters */
	unsigned int true_column;	/*<! supervised data: column of truth */
	unsigned int *true_cluster;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */
	unsigned int true_K;		/*<! true number of clusters */
	data_t **true_modes;		/*<! true modes */
	/* options:true_K <= options:sim_K */
	int require_sim_K;		/*<! require true_K == sim_K */
	double *sim_alpha;		/*<! dirichlet parameters [optional] */
	double *sim_pi;			/*<! mixing proportions */
	double sim_between_t;		/*<! between variation */
	double sim_within_t;		/*<! within variation */
	double sim_within_prob, sim_between_prob;
	unsigned int sim_n_observations;/*<! number of observations */
	unsigned int sim_n_coordinates;	/*<! number of coordinates */
	data_t sim_n_categories;	/*<! number of categories */
	data_t **sim_modes;		/*<! simulation (or true) modes */
	unsigned int *sim_cluster;	/*<! simulated cluster assignments */
}; /* options */

int make_options(options **opt);
void free_options(options *opt);
int parse_options(options *opt, int argc, const char **argv);
int process_arg_p(int argc, char const **argv, int *i, int j, options *opt);

#endif /* HAPLOTYPE_OPTION_H */
