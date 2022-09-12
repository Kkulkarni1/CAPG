/**
 * @file simulate.h
 * @author Karin S. Dorman
 *
 * Header file for ampliclust structs, extern functions and defines.
 */

#ifndef __H_SIMULATE__
#define __H_SIMULATE__

#include "model.h"
#include "data.h"
#include "options.h"
#include "initialize.h"
#include "simulate_options.h"

/**
 * Simulation controls of clustering difficulty.
 */
#define CONTROL_NOTHING		0	/*<! do not alter inferred values */
#define SIMULATE_HAPLOTYPES	1	/*<! simulate haplotype separation */
#define ADJUST_ERROR		2	/*<! control within-cluster spread */
#define ADJUST_QUALITY		4	/*<! additively adjust quality scores */
#define CONTROL_ABUNDANCE	8	/*<! control the abundances of clusters */
#define EXTRA_HAPLOTYPE		16	/*<! introduce extra file as haplotypes */

/**
 * Point estimation parameters at simulation parameters in order to estimate
 * parameters from data or reverse.
 */
#define POINT_TRUTH	0
#define UNPOINT_TRUTH	1

/**
 * Parts of the data.  Used to set \ref options:simulation_random (what to 
 * randomize in simulated data) and \ref options:simulation_from_data (what
 * parameters to estimate from data for simulation).  Here is what the various
 * settings mean:
 *
 * simulation_random & QUALITIES && simulation_from_data & QUALITIES
 *	estimate model:true_lambda from data and simulate qualities
 * !(simulation_random & QUALITIES) && simulation_from_data & QUALITIES
 *	qualities taken from data
 * simulation_random & QUALITIES && !(simulation_from_data & QUALITIES)
 *	simulate model:true_lambda and simulate qualities
 *	[WARNING] ampliclust: not possible because lambda0|1 cannot be simulated
 * !(simulation_random & QUALITIES) && !(simulation_from_data & QUALITIES)
 *	[WARNING] non-sensical
 *
 * simulation_random & NUCLEOTIDES && simulation_from_data & NUCLEOTIDES
 *	mlogit: estimate model:true_beta from data and simulate nucleotides
 *	[WARNING] dada2: not possible because error_profile cannot be estimated
 * !(simulation_random & NUCLEOTIDES) && simulation_from_data & NUCLEOTIDES
 *	nucleotides taken from data (seems weird to do this)
 * simulation_random & NUCLEOTIDES && !(simulation_from_data & NUCLEOTIDES)
 *	mlogit: there is a hard-coded model:true_beta and simulate nucleotides
 *	[WARNING] ampliclust: not possible because delta cannot be simulated
 * !(simulation_random & NUCLEOTIDES) && !(simulation_from_data & NUCLEOTIDES)
 *	[WARNING] non-sensical
 *
 * simulation_random & MIXING && simulation_from_data & MIXING
 *	estimate model:true_pi from data and simulate read origins
 * !(simulation_random & MIXING) && simulation_from_data & MIXING
 *	take and assume partition from data
 * simulation_random & MIXING && !(simulation_from_data & MIXING)
 *	simulate model:true_pi and simulate read origins
 * !(simulation_random & MIXING) && !(simulation_from_data & MIXING)
 *	[WARNING] non-sensical
 *
 * simulation_random & HAPLOTYPES && simulation_from_data & HAPLOTYPES
 *	estimate model:true_haplotypes from data and simulate similar spread (not implemented)//
 * !(simulation_random & HAPLOTYPES) && simulation_from_data & HAPLOTYPES
 *	take and assume haplotypes from data
 * simulation_random & HAPLOTYPES && !(simulation_from_data & HAPLOTYPES)
 *	simulate model:true_haplotypes
 * !(simulation_random & HAPLOTYPES) && !(simulation_from_data & HAPLOTYPES)
 *	[WARNING] non-sensical
 */
#define NOTHING		0	/*<! nothing */
#define QUALITIES	1	/*<! read quality scores */
#define NUCLEOTIDES	2	/*<! read nucleotides */
#define MIXING		4	/*<! partition into clusters */
#define HAPLOTYPES	8	/*<! true haplotypes */

/**
 * Simulation methods.
 */
enum {
	SIMULATE_FROM_DATA,	/*<! simulate from data: see
				 * initialize.h:initialization_methods
				 */
	SIMULATE_FROM_MODEL_WITH_DATA,
				/*<! simulate from a model with partial
				 * information from data as per
				 * options::simulation_from_data
				 */
	SIMULATE_FROM_MODEL,	/*<! simulate from model: see models below */
	NUM_SIMULATION_METHODS
};


/**
 * Simulation models for haplotypes
 * */
enum{
	JC69,
	POSITION_CONTROL
};

typedef struct _simulator simulator;

struct _simulator {
	/* truth when data simulated */
	unsigned int *true_cluster_id;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */

	model *mod;			/*<! model used for simulation */
}; /* simulator */


int make_simulator(simulator **sim_in, data *dat, model *mod, simulate_options *opt);
int do_simulation(simulator *sim, model *mod, data *dat, initializer *ini, simulate_options *opt);
int simulate_read(data *dat, model *mod, unsigned int len, unsigned char *simu_read,unsigned char *simu_qual,unsigned int off);
int simu_hap(unsigned int length, unsigned int K, unsigned int * exp_mut, unsigned int * mut_position, double hap_spread, unsigned char * anc, unsigned char * hap,unsigned int type);
int initialize_from_true_partition(model *mod, data *dat, initializer *ini, simulator *sim, simulate_options *sopt, initialize_options *iopt);
int initialize_from_true_parameters(model *mod, model *sim_mod, data *dat, model_options *mopt, simulate_options *sopt);


#endif
