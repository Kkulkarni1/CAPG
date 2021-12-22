#ifndef __SIMULATE_OPTIONS_H__
#define __SIMULATE_OPTIONS_H__

#include "initialize_options.h"


typedef struct _simulate_options simulate_options;

struct _simulate_options {
	
	model_options *modo;
	initialize_options *inio;
	
	/* simulation */
	unsigned int true_K;		/*<! number of true clusters */
	unsigned int sim_n_reads;	/*<! number of reads */
	unsigned int sim_read_length;	/*<! length of reads */
	int background_model;		/*<! simulate with background model */
	int simulation_random;		/*<! parts of data to simulate */

	/* simulate based on data */
	int simulation_from_data;	/*<! parts of data from real data */
	int simulation_initialization;	/*<! when simulating from real data */
	int simulation_control;		/*<! alterations to data */
	int quality_change;		/*<! prescribed adjustment to quality scores */

	/* simulation: haplotypes */
	unsigned char *ancestor;	/*<! ancestor sequence of */
	double cluster_spread;		/*<! prescribed within-cluster spread */
	double haplotype_spread;	/*<! expected mutations b/w haplotypes */
	char const *hap_control_file;	/*<! file showing mutable positions */

	/* simulation: mixing proportions */
	double sim_pi_alpha;		/*<! Dirichlet parameter for pi */

	/* simulation: substitution probabilities */
	double sim_gamma_alpha;		/*<! Dirichlet parameter for gamma */
	/* simulation parameters are in model::true_* */
	
	/* input */
	char *simulation_infile;	/*<! simulation infile: initialization */
	
	/* output */
	char const *fastq_outfile;	/*<! name of fastq outfile */
	char const *simulation_outfile;	/*<! name of simulation outfile */
	
};

int make_simulate_options(simulate_options **, initialize_options *, model_options *);

#endif /* simulate_option_h */
