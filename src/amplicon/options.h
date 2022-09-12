/**
 * @file options.h
 * @author Karin S. Dorman
 *
 * Header file for options struct.
 */

#ifndef __H_OPTIONS__
#define __H_OPTIONS__

#include <stdio.h>
#ifdef USE_CURSES
#include <curses.h>
#endif

#include "kmodes.h"
#include "model_options.h"
#include "initialize_options.h"
#include "simulate_options.h"

/**
 * Methods to assign reads to provided centers.
 */
enum {
	HAMMING,	/*<! Hamming distance */
	QUALITY		/*<! likelihood with quality scores taken literally */
};

/* every top-level options object should declare its own global_wp and active_fp
 * if it ones to link with aecm.c and initialize.c.
 */
#ifdef USE_CURSES
extern WINDOW *global_wp;							// aecm.c, ampliclust.c, initialize.c			[TODO] move to options_model (change intent?)
#else
extern FILE *global_wp;
#endif
extern FILE *active_fp;								// aecm.c, ampliclust.c					[TODO] move to options_model (change intent?)

typedef struct _options options;

/**
 * Run options.
 */
struct _options {
	/* model */
	model_options *modo;
	model_options *sim_modo;
	initialize_options *inio;
	initialize_options *sim_inio;
	simulate_options *simo;

// [TODO] initialize_options already references model_options, simulate.c references initialization_options as per above
// [TODO] options->initialization_options->model_options, options->simulation_options->model_options [DONE]
// [TODO] since simulator uses initialization code to estimate itself from data, I think we need options->simulation_options->initialization_options [DONE]

	int do_simulation;	/*<! whether to do simulation */
	int do_estimation;	/*<! whether to estimate model */
	int estimate_K;		/*<! whether algorithm should estimate K */		// => options (ampliclust.c, initialize.c, stages.c, run_ampliclust.c)
											// [TODO] options_initialize should get a copy, stages.c can use global_options
	unsigned long seed;	/*<! random number seed [srand()] */			// => options (ampliclust.c)

	/* multi-stages */								// => options_stages
	int multi_stage_method;	/*<! use multi-stage method */
	size_t sample_size;	/*<! per-stage sample size */
	size_t simu_size;	/*<! number of simulation reads for LR test */
	double proportion;	/*<! proportion of the whole dataset to get final estimates */
	double epsilon_s;	/*<! epsilon for identifying null cluster */

	/* input */
	char const *fastq_file;		/*<! name of fastq input file */		// needed by run_ampliclust.c
	char const *reference_file;	/*<! name of reference file */			// needed by run_ampliclust.c
	char const *seed_file;		/*<! name of file with seeds */			// needed by run_ampliclust.c
	char const *offset_file;	/*<! name of offset file */			// needed by data.c

	/* output */
	char const *outfile;		/*<! ... */					// needed by ampiclust.c
	char *outfile_k;		/* output file for different K */		// not used => delete					[TODO] delete
	int use_curses;									// needed by ampliclust.c
}; /* options */

/**
 * Different quality score models (not really used ever).
 */
enum {
	DEFAULT_Q_MODEL,
	FIVE_BY_SIX,
};

int make_options(options **opt);
int parse_options(options *opt, int argc, const char **argv);
void free_options(options *opt);

#endif
