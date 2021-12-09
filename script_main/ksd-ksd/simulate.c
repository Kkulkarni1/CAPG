/**
 * @file simulate.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 *
 * Simulate amplicon reads.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "simulate.h"
#include "data.h"
#include "model.h"
#include "options.h"
#include "ampliclust.h"
#include "aecm.h"
#include "initialize.h"
#include "math.h"
#include "util.h"
#include "matrix_exponential.h"
#include "io.h"
#include "simulate_options.h"

int set_simulation_parameters(data *dat, model *mod, initializer *ini, simulate_options *opt);
void simulate_pi(model *mod, simulate_options *opt);
void simulate_gamma(model *mod, simulate_options *opt);
int simulate_lambda(model *mod, data *dat, simulate_options *opt);
int simulate_data(simulator *sim, data *dat, model *mod, simulate_options *opt);
int simulate_ancestor(unsigned char **seq, unsigned int len);
int simulate_haplotypes(unsigned int length, unsigned int K, unsigned int *mut_position, double hap_spread, unsigned char *anc, unsigned char *hap);
void simulate_beta(model *mod, simulate_options *sopt);
void simulate_position(simulator *sim, data *dat, model *mod, simulate_options *opt, unsigned int i, unsigned int j);
void point_params_at_simulation(model *mod, data *dat, simulate_options *opt, int point);
void estimate_lambda(model *mod, data *dat);

void simple_update(data *dat, model *mod, initializer *ini, void *obj);
int sim_simple_initialize(data *dat, model *mod, initializer *ini, void *obj);



/**
 * Make simulator object.
 *
 * @param sim_in	simulator object ot allocate
 * @param dat		pointer to data object
 * @param mod		pointer to simulation model object
 * @param opt		pointer to simulate_options object
 * @return		error status
 */
int make_simulator(simulator **sim_in, data *dat, model *mod, simulate_options *opt)
{
	simulator *sim;
	*sim_in = malloc(sizeof **sim_in);

	if (!*sim_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "simulator");

	sim = *sim_in;

	sim->mod = mod;

	sim->true_cluster_id = malloc(dat->fdata->n_reads
		* sizeof *sim->true_cluster_id);
	sim->true_cluster_size = malloc(opt->true_K
		* sizeof *sim->true_cluster_size);

	if (!sim->true_cluster_id || !sim->true_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"simulator.true_cluster_*");
	
	if (opt->modo->parameterization == ART_PARAMETERIZATION
		&& opt->simulation_random & QUALITIES
		&& opt->simulation_from_data & QUALITIES)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "remove 'quality"
			"' from -f argument: art profile cannot be estimated "
			"from data.\n");

	/* data has more quality scores that DADA2 error profile allows */
	if (opt->modo->parameterization == DADA2_PARAMETERIZATION
		&& opt->simulation_from_data & QUALITIES
		&& dat->n_quality > mod->n_quality)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "DADA2 cannot "
			"generate reads conditional on quality scores exceeding"
			" %u\n", mod->n_quality);

	/* data has less quality scores that DADA2 error profile allows */
	if (opt->modo->parameterization == DADA2_PARAMETERIZATION
		&& opt->simulation_random & QUALITIES
		&& opt->simulation_from_data & QUALITIES) {
		if (dat->n_quality < mod->n_quality)
			mmessage(WARNING_MSG, INTERNAL_ERROR, "Quality score "
				"distribution to be estimated from data with %u"
				" quality scores, but DADA2 profile hanldes %u "
				"quality scores.  Highest quality scores are "
				"truncated.\n", dat->n_quality, mod->n_quality);
		mod->lambda = malloc(dat->max_read_position * mod->n_quality * sizeof *mod->lambda);
		if (!mod->lambda)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model::lambda");
	}
	
	/* randomizing qualities in a model that conditions on qualities */
	if ((opt->modo->parameterization == QUALITY_PARAMETERIZATION
		|| opt->modo->parameterization == MLOGIT_PARAMETERIZATION)
		&& opt->simulation_random & QUALITIES) {

		if (opt->simulation_from_data & QUALITIES)
			mod->n_quality = dat->n_quality;
		else
			mod->n_quality = MAX_ILLUMINA_QUALITY_SCORE
				- MIN_ILLUMINA_QUALITY_SCORE + 1;

		mod->lambda = malloc(dat->max_read_position * mod->n_quality * sizeof *mod->lambda);
		if (!mod->lambda)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model::lambda");
	}

	return NO_ERROR;
} /* make_simulator */


/**
 * Simulate data.  Set simulation parameters, simulate data, compute log
 * likelihood of data under true parameter values, write results file, and
 * write fastq data file.
 *
 * @param sim	simulator object
 * @param mod	model model
 * @param dat	data object
 * @param ini	initializer object
 * @param opt	options object
 * @return	error status
 */
int do_simulation(simulator *sim, model *mod, data *dat, initializer *ini,
						simulate_options *opt)
{
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
	int err = NO_ERROR;

	/* Simulate parameter values, storing them in default slots of
	 * model structure, such as model.pi.
	 */
	if ((err = set_simulation_parameters(dat, mod, ini, opt)))
		return err;

	debug_call(DEBUG_II <= fxn_debug, fxn_debug, fprint_fasta(stderr,
		mod->haplotypes, opt->true_K, dat->max_read_length, "H"));

	/* Simulate data using these parameters: sets
	 * \par simulation_options.true_cluster_id and
	 * \par simulation_options.true_cluster_size
	 */
	simulate_data(sim, dat, mod, opt);	/* overwrites data in memory */

	/* update no. quality scores in data */

	dat->n_quality = dat->fdata->max_quality - dat->fdata->min_quality + 1;
	dat->fdata->empty = 0;

	/* update no. quality scores in model if now exceeded */
/*
	if (mod->n_quality < dat->fdata->max_quality - dat->fdata->min_quality + 1
		&& (err = realloc_quality_information(mod, dat, opt->modo,
		dat->fdata->max_quality - dat->fdata->min_quality + 1)))
		return err;
*/

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "cluster sizes: ");
	debug_call(DEBUG_I <= fxn_debug, fxn_debug, fprint_uints(stderr,
		sim->true_cluster_size, opt->true_K, 4, 1));

	/* write simulated data as fastq file */
	if (opt->fastq_outfile) {
		fastq_options fop = {.outfile = opt->fastq_outfile,
			.read_encoding = XY_ENCODING};
		write_fastq(dat->fdata, &fop);
	}

	return err;
}/* do_simulation */


/**
 * Choose simulation parameters and place in model::* parameters.
 * Simulations may use real fastq files and some partition of the reads
 * to estimate simulation parameters that mimic real data.  To adjust
 * these simulations to control the difficulty of clustering, several
 * strategies are available, including providing one's own haplotypes in
 * a file, simulating haplotypes from JC69, shifting the quality scores
 * by a fixed integer amount, adjusting the errors by a constant factor.
 *
 * There is also code to generate de novo simulation parameters, for all
 * but the AMPLICLUST_MODEL, from some prior distributions.
 *
 * @param mod	model object pointer
 * @param dat	data object pointer
 * @param ini	initializer object pointer
 * @param opt	simulation options object pointer
 * @return	error status
 */
int set_simulation_parameters(data *dat, model *mod, initializer *ini,
					simulate_options *opt)
{
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
	int err = NO_ERROR;

	/* A silly, first attempt to control cluster spread, linearly adjust
	 * integer quality scores by fixed integer amount.
	 * [KSD, TODO] Does not work if there is no data (data::fdata::empty),
	 * but I don't think we ever use this simulation mechanism any more.
	 */
	if (opt->simulation_control & ADJUST_QUALITY)
		for (size_t i = 0; i < dat->sample_size; ++i)
			for (size_t j = 0; j < dat->lengths[i]; ++j) {
				if (opt->quality_change > 0)
					dat->qmat[i][j] = MIN(MAX_QUALITY_SCORE,
						dat->qmat[i][j]
						+ opt->quality_change);
				else
					dat->qmat[i][j] = MAX(MIN_QUALITY_SCORE,
						dat->qmat[i][j]
						+ opt->quality_change);
			}

	if (!dat->fdata->empty && opt->simulation_from_data == NOTHING)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "You have "
				"provided an input fastq file, but have not asked "
				"to estimate anything from it: -f -s --art -m\n");

	/* If user provides no fastq file, but just requests number and length
	 * of reads from command line, then there is no real data from which to
	 * estimate model parameters.  This code simulates parameters from
	 * "prior" distributions.  Alternatively, we could get parameters from
	 * user. [TODO,KSD]
	 */
	if (dat->fdata->empty && opt->simulation_from_data == NOTHING) {

		/* check for invalid input */
		if (dat->fdata->empty && !opt->sim_n_reads)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Must "
				"provide amount of data to simulate: -s -f.\n");

		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
				"simulate ampliclust model without fastq input "
				"file because lambda0|1 unknown: -s -m -f\n");

		/* everything is random and nothing is estimated from data */
		opt->simulation_random = MIXING | QUALITIES | NUCLEOTIDES
								| HAPLOTYPES;
		/* simulate model::pi (mixing proportions) */
		simulate_pi(mod, opt);

		/* simulate model::gamma (substitution generation) */
		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
			|| opt->modo->parameterization == QUALITY_PARAMETERIZATION
			|| opt->modo->parameterization == ART_PARAMETERIZATION)
			simulate_gamma(mod, opt);

		/* "simulate" model::beta multinomial logit parameters */
		if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION)
			simulate_beta(mod, opt);

		/* simulate model::lambda (quality score generation) */
		if (opt->modo->parameterization == QUALITY_PARAMETERIZATION
			|| opt->modo->parameterization == MLOGIT_PARAMETERIZATION
			|| opt->modo->parameterization == DADA2_PARAMETERIZATION)
			simulate_lambda(mod, dat, opt);

		if (!opt->ancestor)
			simulate_ancestor(&opt->ancestor, dat->max_read_length);

		/* [TODO, KSD] allow model::mut_positions */
		if ((err = simulate_haplotypes(dat->max_read_length,
			opt->true_K, NULL, opt->haplotype_spread, opt->ancestor,
							mod->haplotypes)))
			return err;

		return NO_ERROR;

	/* If user requests simulation parameters to be estimated from data,
	 * then we need real data to estimate those parameters.  However, if
	 * the user only requests quality scores from data and not under
	 * ampliclust model, then a partition is not necessary.
	 */
	} else if ((opt->simulation_from_data & QUALITIES)
					== opt->simulation_from_data &&
		opt->modo->parameterization != DEFAULT_PARAMETERIZATION) {


		/* check for invalid input */
		if (dat->fdata->empty)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Simulation settings require user to input "
				"fastq file: see -s or -f.\n");

		/* It is not possible to simulate when we are not to obtain
		 * quantities from data and they are not to be randomly
		 * generated.
		 */
		if (!(opt->simulation_random & QUALITIES)
			&& !(opt->simulation_from_data & QUALITIES))
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Simulation settings require fastq input "
				"file for qualities: see "
				"options.simulation_random and "
				"options.simulation_from_data\n");

		if (!(opt->simulation_random & NUCLEOTIDES)
			&& !(opt->simulation_from_data & NUCLEOTIDES))
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Simulation settings require fastq input "
				"file for nucleotides: see "
				"options.simulation_random and "
				"options.simulation_from_data\n");

		if (!(opt->simulation_random & MIXING)
			&& !(opt->simulation_from_data & MIXING))
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Simulation settings require fastq input "
				"file for mixing proportions: see "
				"options.simulation_random and "
				"options.simulation_from_data\n");

		if (!(opt->simulation_random & HAPLOTYPES)
			&& !(opt->simulation_from_data & HAPLOTYPES)
			&& !opt->simulation_infile)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Simulation settings require fastq input "
				"file for haplotypes: see "
				"options.simulation_random and "
				"options.simulation_from_data\n");

		/* Cannot estimate dada2 error profiles from data. [TODO] */
		if (opt->modo->parameterization == DADA2_PARAMETERIZATION
			&& opt->simulation_random & NUCLEOTIDES
			&& opt->simulation_from_data & NUCLEOTIDES)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Cannot estimate DADA2 substitution profiles "
				"from data: see -s or -f.\n");

		/* Cannot estimate art quality profiles from data. [TODO] */
		if (opt->modo->parameterization == ART_PARAMETERIZATION
			&& opt->simulation_random & QUALITIES
			&& opt->simulation_from_data & QUALITIES)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Cannot estimate ART substitution profiles "
				"from data: see -s or -f.\n");

		/* Cannot estimate lambda0|1 of ampliclust model from data. */
		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
			&& opt->simulation_random & QUALITIES
			&& !(opt->simulation_from_data & QUALITIES))
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Cannot simulate ampliclust quality score "
				"distributions: see -s or -f.\n");

		/* Cannot estimate delta of ampliclust model from data. */
		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
			&& opt->simulation_random & NUCLEOTIDES
			&& !(opt->simulation_from_data & NUCLEOTIDES))
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Cannot simulate ampliclust error "
				"distributions delta: see -s or -f.\n");

		/* simulate model::gamma (substitution generation) */
		if (opt->simulation_random & NUCLEOTIDES
			&& !(opt->simulation_from_data & NUCLEOTIDES)) {
			if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
				|| opt->modo->parameterization == QUALITY_PARAMETERIZATION
				|| opt->modo->parameterization == ART_PARAMETERIZATION)
				simulate_gamma(mod, opt);

			/* "simulate" model::beta multinomial logit parameters */
			else if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION)
				simulate_beta(mod, opt);


			/* [TODO, KSD] verify read_dada2_error_profile */
		}

		/* simulate model::pi (mixing proportions) */
		if (opt->simulation_random & MIXING
			&& !(opt->simulation_from_data & MIXING))
			simulate_pi(mod, opt);

		/* simulate model::lambda (quality score generation) */
		if (opt->simulation_random & QUALITIES
			&& (opt->modo->parameterization == QUALITY_PARAMETERIZATION
			|| opt->modo->parameterization == MLOGIT_PARAMETERIZATION
			|| opt->modo->parameterization == DADA2_PARAMETERIZATION)) {
			if (!(opt->simulation_from_data & QUALITIES))
				simulate_lambda(mod, dat, opt);
			else if (opt->modo->parameterization != DEFAULT_PARAMETERIZATION)
				estimate_lambda(mod, dat);
		}

		if (opt->simulation_random & HAPLOTYPES) {

			if (opt->simulation_from_data & HAPLOTYPES)
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"Not implemented.  See "
					"options:simulation_random and "
					"options:simulation_data_from\n");

			if (!opt->ancestor)
				simulate_ancestor(&opt->ancestor,
						dat->max_read_length);

			/* [TODO, KSD] allow model::mut_positions */
			simulate_haplotypes(dat->max_read_length, opt->true_K,
				NULL, opt->haplotype_spread, opt->ancestor,
				mod->haplotypes);
		} else {
			if ((err = read_initialization_file(
				opt->simulation_infile, dat, opt->inio, ini)))
				return err;
			/* technically may not need partition */
			if ((err = pull_bootstraps(dat, opt->inio, ini)))
				return err;
			for (unsigned char k = 0;k < opt->true_K; ++k)
				memcpy(&mod->haplotypes[k * dat->max_read_length],
					ini->best_modes[k], dat->max_read_length
						* sizeof *mod->haplotypes);
		}

		return NO_ERROR;
	}

	/* Otherwise, existing data produce simulation parameter estimates, with
	 * possible override of some parameters.
	 */

	/* To estimation simulation parameters, we must have a partition
	 * of the data.  This code reads a partition, set of seeds, or
	 * both from user-provided file.
	 */
//	if ((err = read_initialization_file(opt->simulation_infile,
//						dat, opt->inio, ini)))
//		return err;

	/* restore intended initialization method */
//	opt->simulation_initialization = opt->inio->initialization_method;
//	opt->inio->initialization_method = init_method;
	/* don't initialize estimation w/ wrong info */
//	opt->inio->partition_file = NULL;

//	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Set simulation "
//		"initialization to %d\n", opt->simulation_initialization);

	/* We need a partition of the data and haplotypes in order to
	 * estimate simulation parameters from data, but not all
	 * initialization methods set these, so this code makes
	 * sure \par initializer.best_cluster_id and \par initializer.best_modes
	 * are set.
	 */
//	if ((err = pull_bootstraps(dat, opt->inio, ini)))
//		return err;

	/* Estimate parameters from this partition of data: again
	 * we will use the same code that estimates/initializes
	 * parameters in the estimation stage, so we need to protect
	 * the requested settings.					// [TODO]  Well, not if we are using the model ONLY to do simulation, as would be the case for use in k-haplotypes.
	 */
	
/* [WORKING] make sure param_estimate is set correctly for the assumed parameterization
*/
	opt->modo->param_estimate = 0;

	if (opt->simulation_from_data & MIXING)
		opt->modo->param_estimate |= PARAM_PI;

	if (opt->simulation_from_data & HAPLOTYPES)
		opt->modo->param_estimate |= PARAM_HAPLOTYPE;

	if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION) {

		if (opt->simulation_from_data & NUCLEOTIDES)
			opt->modo->param_estimate |= PARAM_DELTA | PARAM_GAMMA;

		if (opt->simulation_from_data & QUALITIES)
			opt->modo->param_estimate |= PARAM_LAMBDA;
	}

	if ((opt->modo->parameterization == QUALITY_PARAMETERIZATION
		|| opt->modo->parameterization == ART_PARAMETERIZATION)
		&& opt->simulation_random & NUCLEOTIDES)
			opt->modo->param_estimate |= PARAM_GAMMA;

	if (fxn_debug) {
		if (opt->modo->param_estimate & PARAM_PI)
			mmessage(INFO_MSG, NO_ERROR, "Estimating pi from "
				"data.\n");
		if (opt->modo->param_estimate & PARAM_HAPLOTYPE)
			mmessage(INFO_MSG, NO_ERROR, "Estimating haplotypes "
				"from data.\n");
		if (opt->modo->param_estimate & PARAM_DELTA)
			mmessage(INFO_MSG, NO_ERROR, "Estimating delta "
				"from data.\n");
		if (opt->modo->param_estimate & PARAM_GAMMA)
			mmessage(INFO_MSG, NO_ERROR, "Estimating gamma "
				"from data.\n");
		if (opt->modo->param_estimate & PARAM_LAMBDA)
			mmessage(INFO_MSG, NO_ERROR, "Estimating lambda "
				"from data.\n");
	}

	/* trick initializer to set true parameter values
	 * [TODO?] Should true parameter values be in simulator object?
	 * [TODO?] Should we not use any parameters from model (see point_params_at_simulation())?  (Now it uses the model:n* variants when available.
	 */
//	point_params_at_simulation(mod, dat, opt, POINT_TRUTH);
//	initialize_model_parameters(mod, dat, opt->inio, ini, BEST_VALUES);
//	point_params_at_simulation(mod, dat, opt, UNPOINT_TRUTH);
//	opt->modo->param_estimate = param_est;

	if (opt->simulation_random & MIXING && opt->modo->K != opt->true_K
					&& opt->simulation_from_data & MIXING)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
			"obtain K=%u mixing proportions from data and "
			"simulate with K=%u clusters: see -ktrue, -sk, "
			"and -f.\n", opt->modo->K, opt->true_K);

	int store_val = opt->modo->em_info;
	opt->modo->em_info = ABSOLUTE_SILENCE;
	cluster_amplicons(dat, mod, ini, opt->modo, opt->inio, simple_update,
					sim_simple_initialize, (void *) opt);
	opt->modo->em_info = store_val;


	/* overwrite selected parameters by simulating them */
	/* simulate model::gamma (substitution generation) */
	if (opt->simulation_random & NUCLEOTIDES
		&& !(opt->simulation_from_data & NUCLEOTIDES)) {
		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
			|| opt->modo->parameterization == QUALITY_PARAMETERIZATION
			|| opt->modo->parameterization == ART_PARAMETERIZATION)
			simulate_gamma(mod, opt);

		/* "simulate" model::beta multinomial logit parameters */
		else if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION)
			simulate_beta(mod, opt);


		/* [TODO, KSD] verify read_dada2_error_profile */
	}

	/* mixing proportions complicated as they may differ between
	 * estimation model and simulation model
	 */
	if (!(opt->simulation_random & MIXING) && opt->modo->K != opt->true_K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
			"obtain membership in K=%u clusters from data and "
			"simulate from K=%u clusters: see -ktrue and -sk\n",
			opt->modo->K, opt->true_K);
	/* renormalize to exclude background model component */
	if (opt->simulation_random & MIXING && opt->modo->background_model
		&& opt->simulation_from_data & MIXING) {
		double sum = 0;
		for (unsigned int k = 0; k < opt->true_K; ++k)
			sum += mod->pi[k];
		for (unsigned int k = 0; k < opt->true_K; ++k)
			mod->pi[k] /= sum;
	}

	/* simulate model::pi (mixing proportions) */
	if (opt->simulation_random & MIXING
		&& !(opt->simulation_from_data & MIXING)) {
		if (opt->background_model)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Cannot simulate data with background model\n");
		/* we can estimate parameters with background model or any
		 * number K that does not match simulation model: reset
		 * simulation model here to the one desired.
		 */
		if (opt->modo->K + opt->modo->background_model
			!= opt->true_K + opt->background_model) {
			mod->n_mix = opt->true_K + opt->background_model;
			double *tmp = realloc(mod->pi, mod->n_mix);
			if (!tmp)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"model::pi");
			mod->pi = tmp;
			/* [NOTE] model::npi and model::best_pi not updated */
		}
		simulate_pi(mod, opt);
	}

	/* simulate model::lambda (quality score generation) */
	if (opt->simulation_random & QUALITIES
		&& (opt->modo->parameterization == QUALITY_PARAMETERIZATION
		|| opt->modo->parameterization == DADA2_PARAMETERIZATION)) {
		if (!(opt->simulation_from_data & QUALITIES))
			simulate_lambda(mod, dat, opt);
		else if (opt->modo->parameterization != DEFAULT_PARAMETERIZATION)
			estimate_lambda(mod, dat);
	}

	/* simulate haplotypes, overwriting those estimated from data */
	if (opt->simulation_random & HAPLOTYPES) {

		if (opt->simulation_from_data & HAPLOTYPES)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Not implemented.  See "
				"options:simulation_random and "
				"options:simulation_data_from\n");

		if (!opt->ancestor)
			simulate_ancestor(&opt->ancestor,
					dat->max_read_length);

		/* [TODO, KSD] allow model::mut_positions */
		simulate_haplotypes(dat->max_read_length, opt->true_K,
			NULL, opt->haplotype_spread, opt->ancestor,
			mod->haplotypes);
	}

	/* The rest of the code overwrites parameters as requested by the
	 * user, with a goal to control the difficulty of clustering.
	 */

	/* [TODO, KSD] This is only a quick solution to overwrite model::pi.  */
	if (opt->simulation_control & CONTROL_ABUNDANCE) {

		/* read the abundance file */
		/*
		FILE *file = fopen("abundance.csv","rb");
		if(!file)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,"abundance.csv");

		double rate;

		for (unsigned int k = 0; k < opt->true_K; k++){
			rate = 0.;
			fscanf(file,"%lf,",&rate);
			mod->pi[k] = rate;
		}
		*/
		for(unsigned int k = 0 ; k < opt->true_K; k++)
			mod->pi[k] = 1. / opt->true_K;
	}

	/* [TODO, KSD] This is a quick and dirty hard-coded solution to overwrite
	 * model::haplotypes from file.  It assumes 'haplotypes.fsa' exists
	 * and contains options::K sequences with the right length.  Many
	 * opportunities for error.
	 */
	if (opt->simulation_control & EXTRA_HAPLOTYPE) {

		char hap_file[] = "haplotypes.fsa";

		if ((err = read_initialization_file(hap_file, dat,
							opt->inio, ini)))
			return err;
		for (unsigned char k = 0;k < opt->true_K; ++k)
			memcpy(&mod->haplotypes[k * dat->max_read_length],
				ini->best_modes[k], dat->max_read_length
						* sizeof *mod->haplotypes);
	}


	/* Overwrite haplotypes with simulated haplotypes under JC69 with
	 * user-specified evolutionary separation.
	 */
	if (opt->simulation_control & SIMULATE_HAPLOTYPES) {

		/* simulate a center if one has not been provided: assumed
		 * iid Uniform(A, C, G, T).
		 */
		if (!opt->ancestor)
			simulate_ancestor(&opt->ancestor, dat->max_read_length);

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Simulating "
			"haplotypes from ancestor %s\n", display_sequence(
			opt->ancestor, dat->max_read_length, XY_ENCODING));

		/* optionally read locations where haplotypes may mutate */
		/* [TODO, KSD] Move to setup stage: make_model() */
		if (opt->hap_control_file) {
			mod->mut_position = malloc(dat->max_read_length
				* sizeof *mod->mut_position);
			if (!mod->mut_position)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"model::mut_position");
			if ((err = read_uints(opt->hap_control_file,
				mod->mut_position, dat->max_read_length)))
				return err;
		}

		/* simulate options::K haplotypes */
		simulate_haplotypes(dat->max_read_length, opt->true_K,
			mod->mut_position, opt->haplotype_spread,
			opt->ancestor, mod->haplotypes);

		fprint_fasta(stderr, mod->haplotypes, opt->true_K,
			dat->max_read_length, "H");
	}

	/* Overwrite delta (probability of error) according to requested
	 * amount of error (assumes AMPLICLUST_MODEL).
	 */
	if (opt->simulation_control & ADJUST_ERROR
		&& opt->simulation_random & NUCLEOTIDES) {

		if (!opt->modo->parameterization != DEFAULT_PARAMETERIZATION)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
				"use command line option -v without -m "
				"ampliclust");

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Adjusting error "
			"rate by %f.\n", opt->cluster_spread);

		if (opt->cluster_spread < 1)
			for (size_t i = 0; i < dat->max_read_position; ++i)
				mod->delta[i] = 1 - opt->cluster_spread
						* (1 - mod->delta[i]);
		else
			for (size_t i = 0; i < dat->max_read_position; ++i)
				mod->delta[i] = (opt->cluster_spread - 1.0)
							* mod->delta[i];
	}

	return err;
} /* set_simulation_parameters */

void simple_update(data *dat, model *mod, initializer *ini, void *obj)
{
	UNUSED(ini);
	simulate_options *sopt = (simulate_options *) obj;

	if (mod->ll >= mod->best_ll) {
		mod->ll = e_step(dat, mod, sopt->modo, FIRST);
		fprintf(stderr, "%15.3f (%15.3f): ", mod->ll, mod->best_ll);
		assign_clusters(mod->eik, mod->n_mix, dat->sample_size,
			ini->best_cluster_size, ini->best_cluster_id, 0);
		param_update(mod, dat, sopt->modo, OPTIMAL_VALUES);
		fprint_uints(stderr, ini->best_cluster_size, mod->n_mix, 3, 1);
	} else {
		assign_clusters(mod->eik, mod->n_mix, dat->sample_size,
			ini->cluster_size, ini->cluster_id, 0);
		fprintf(stderr, "%15.3f (%15.3f): ", mod->ll, mod->best_ll);
		fprint_uints(stderr, ini->cluster_size, mod->n_mix, 3, 1);
	}
} /* simple_update */

int sim_simple_initialize(data *dat, model *mod, initializer *ini, void *obj)
{
	simulate_options *sopt = (simulate_options *) obj;
//	return simple_initialize(dat, mod, ini, (void *) stop->inio);
	return initialize(mod, dat, ini, sopt->inio);
} /* sim_simple_initialize */


/**
 * Simulate coefficients of multinomial logistic regression model.
 * Actually, we just have one simulation value now.
 *
 * @param mod	pointer to model object
 * @param opt	pointer to simulate_options object
 */
void simulate_beta(model *mod, simulate_options *sopt)
{
	UNUSED(sopt);

	/* [TODO] hard-coded simulation values! */
	for (unsigned int i = 0; i < mod->n_beta_coef; ++i) {
		/* A->N for all N */
		mod->beta[i] = 0;

		/* N->A for all N */
		if (!(i % mod->n_predictors))
			mod->beta[i] = -4.5;//-8;
		/* N->N for all N!=A */
		else if ((i + mod->n_predictors) / mod->n_predictors
			== (i % mod->n_predictors))
			mod->beta[i] = 4.5;//8;
	}
} /* simulate_beta */



/**
 * Simulate mixing proportions.
 *
 * @param mod	pointer to model object
 * @param opt	pointer simulate_options object
 */
void simulate_pi(model *mod, simulate_options *opt)
{
	double sum = 0;
	for (unsigned int k = 0; k < opt->true_K; ++k) {
		mod->pi[k] = rgamma(opt->sim_pi_alpha, 1.0);
		sum += mod->pi[k];
	}
	for (unsigned int k = 0; k < opt->true_K; ++k)
		mod->pi[k] /= sum;
} /* simulate_pi */


/**
 * Simulate conditional substitution probabilities.
 *
 * @param mod	pointer to model object
 * @param opt	pointer to options object
 */
void simulate_gamma(model *mod, simulate_options *opt)
{
	double sum;

	for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
		sum = 0;
		for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b) {
			if (b == h) {
				mod->gamma[h * NUM_NUCLEOTIDES + b] = 0;
				continue;
			}
			mod->gamma[h * NUM_NUCLEOTIDES + b]
					= rgamma(opt->sim_gamma_alpha, 1.0);
			sum += mod->gamma[h * NUM_NUCLEOTIDES + b];
		}
		for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b)
			mod->gamma[h * NUM_NUCLEOTIDES + b] /= sum;
	}
} /* simulate_gamma */


/**
 * Simulate lambda, quality score pmf.
 *
 * This is silly model that generates roughly triangular pmf with
 * peak moving from right to left along read length.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @return	error status
 */
int simulate_lambda(model *mod, data *dat, simulate_options *opt)
{
	UNUSED(opt);
	int fxn_debug = ABSOLUTE_SILENCE;// DEBUG_II;//
	double sum;
	double *tmp = malloc(mod->n_quality * sizeof *tmp);

	if (!tmp)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"temporary memory");

	double max = 200, min = 1;
	unsigned int width = 10;
	unsigned int first_posn = dat->max_read_position - 20;
	unsigned int last_posn = dat->max_read_position;
	unsigned char min_qual = 20;
	unsigned char max_qual = mod->n_quality;
	double decrement = (max - min) / width;
	unsigned int shift = (last_posn - first_posn) / mod->n_quality;
	unsigned int remainder = (last_posn - first_posn) % mod->n_quality;
	unsigned int cnt = 0, times = 0;
	unsigned int peak_shift = (max_qual - min_qual) / (last_posn - first_posn);
	unsigned int peak_remainder = (max_qual - min_qual) % (last_posn - first_posn);
	unsigned int peak = max_qual - 1;
	unsigned int peak_times = 0;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "n_quality = %u; remainder "
		"= %u; shift = %u; peak_remainder =%u; peak_shift = %u\n",
		mod->n_quality, remainder, shift, peak_remainder, peak_shift);

	for (unsigned int j = 0; j < dat->max_read_position; ++j) {
		int k = peak;
		tmp[k] = max;
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "peak=%d, times=%u, "
			"peak_times=%u, posn=%u: ", peak, times, peak_times, j);
		while (++k < mod->n_quality)
			tmp[k] = MAX(min, max - (k - peak) * decrement);
		k = peak;
		while (--k > 0)
			tmp[k] = MAX(min, max - (peak - k) * decrement);
		tmp[0] = 0;
		sum = 0;
		for (k = 0; k < mod->n_quality; ++k) {
fprintf(stderr, "tmp[%u]: %f\n", k, tmp[k]);
			mod->lambda[j * mod->n_quality + k]
							= rgamma(tmp[k], 1.0);
			sum += mod->lambda[j * mod->n_quality + k];
		}
		for (k = 0; k < mod->n_quality; ++k) {
			mod->lambda[j * mod->n_quality + k] /= sum;
			debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, " %f",
				mod->lambda[j * mod->n_quality + k]);
		}
		debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "\n");
		if (j < first_posn || j >= last_posn)
			continue;
		if (shift > 1 && ++cnt >= shift + 1 && times < remainder) {
			++times;
			peak = MAX(min_qual,  (int) peak - 1);
			cnt = 0;
		} else if (shift > 1 && cnt >= shift && times >= remainder) {
			peak = MAX(min_qual,  (int) peak - 1);
			cnt = 0;
		} else if (shift == 0 && peak_times < peak_remainder) {
			peak = MAX(min_qual, (int) peak - (int) peak_shift - 1);
			++peak_times;
		} else if (shift == 0 && peak_times >= peak_remainder) {
			peak = MAX(min_qual, (int) peak - (int) peak_shift);
		}
	}

	free(tmp);

	return NO_ERROR;
} /* simulate_lambda */


/**
 * Estimate model::lambda from data.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 */
void estimate_lambda(model *mod, data *dat)
{
	/* estimate pmf from data */
	for (unsigned int j = 0; j < dat->max_read_position; ++j)
		for (unsigned char q = 0; q < mod->n_quality; ++q)
			mod->lambda[j*mod->n_quality + q] = 1e-12;

	for (size_t i = 0; i < dat->sample_size; ++i) {
		unsigned int j_off = dat->offset ? dat->offset[i] : 0;
		for (unsigned int j = 0; j < dat->lengths[i]; ++j)
			++mod->lambda[(j + j_off) * mod->n_quality
				+ dat->qmat[i][j]];
	}
	for (unsigned int j = 0; j < dat->max_read_position; ++j) {
		double tmp = 0.;
		for (unsigned char q = 0; q < mod->n_quality; ++q)
			tmp += mod->lambda[j*mod->n_quality + q];
		for (unsigned char q = 0; q < mod->n_quality; ++q)
			mod->lambda[j*mod->n_quality + q] /= tmp;
	}
} /* estimate_lambda */


/**
 * Simulate data. It uses the model::* parameters in, and stores the true
 * cluster ids and sizes in \ref simulation_options.true_cluster_id and
 * \ref simulation_options.true_cluster_size.  It stores the simulated data in
 * the \ref data object, possibly overwriting the previous data that was read
 * from * \ref opt.fastq_file.
 *
 * @param sim	simulator object pointer
 * @param dat	data object pointer
 * @param mod	model object pointer
 * @param opt	options object pointer
 * @return	error status
 */
int simulate_data(simulator *sim, data *dat, model *mod, simulate_options *opt)
{
	double r, cdf;
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//

	/* initialize empty clusters */
	if (opt->simulation_random & MIXING)
		for (size_t i = 0; i < opt->true_K; ++i)
			sim->true_cluster_size[i] = 0;

	/* compute transition probabilities once for mlogit model */
	if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION)
		compute_transition_probabilities(mod, dat, TRUE_VALUES);

	/* simulate each read */
	for (size_t i = 0; i < dat->sample_size; ++i) {

		/* choose cluster id */
		if (opt->simulation_random & MIXING) {
			r = rand() / (RAND_MAX + 1.);

			sim->true_cluster_id[i] = 0;
			cdf = mod->pi[sim->true_cluster_id[i]];
			while (r > cdf) {
				cdf += mod->pi[++sim->true_cluster_id[i]];
				if (sim->true_cluster_id[i] == opt->true_K - 1)
					break;
			}

			sim->true_cluster_size[sim->true_cluster_id[i]]++;

			debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Assigning "
						"read %d to cluster %d.\n", i,
						sim->true_cluster_id[i]);

		}

		/* simulate each position */
		for (size_t j = 0; j < dat->lengths[i]; ++j)
			simulate_position(sim, dat, mod, opt, i, j);

/*
		fprintf(stderr, "Simulate read %u:\n", i);
		fprintf(stderr, "Ancestor: %s\n", display_sequence(&mod->haplotypes[sim->true_cluster_id[i]*dat->max_read_length], dat->max_read_length, XY_ENCODING));
		fprintf(stderr, "    Read: %s\n", display_sequence(dat->dmat[i], dat->max_read_length, XY_ENCODING));
		fprintf(stderr, "   Quals: %s\n", display_quals(&dat->qmat[i], dat->max_read_length, dat->fdata->min_quality));
*/
	}

	return NO_ERROR;
} /* simulate_data */


/**
 * Simulate a single nucleotide and quality score of a read.
 *
 * @param sim	pointer to simulator object
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param opt	pointer to options object
 * @param i	index of read
 * @param j	index of position
 * @return	void
 */
void simulate_position(simulator *sim, data *dat, model *mod, simulate_options *opt,
						unsigned int i, unsigned int j)
{
	int fxn_debug = ABSOLUTE_SILENCE;// DEBUG_I;//
	data_t sim_nuc = 4, sim_q;	/* simulated nuc, quality pair */
	unsigned char hnuc		/* true haplotype nucleotide */
		= mod->haplotypes[sim->true_cluster_id[i]
						* dat->max_read_length + j];
	unsigned int j_off =		/* read position offset */
		dat->offset ? dat->offset[i] : 0;
	double cdf, ep, r;		/* auxiliary */
	double *pmf;
	unsigned int bjhq_idx, jhq_idx;

	if (opt->modo->parameterization == QUALITY_PARAMETERIZATION
		|| opt->modo->parameterization == DADA2_PARAMETERIZATION
		|| opt->modo->parameterization == ART_PARAMETERIZATION
		|| opt->modo->parameterization == MLOGIT_PARAMETERIZATION) {

		/* simulate quality score */
		if (opt->simulation_random & QUALITIES && !mod->art_profile) {
			r = rand() / (RAND_MAX + 1.);
			cdf = mod->lambda[(j + j_off) * mod->n_quality];
			sim_q = 0;
			while (r > cdf) {
				cdf += mod->lambda[(j + j_off)
					* mod->n_quality + ++sim_q];
				if (sim_q == mod->n_quality - 1)
					break;
			}

			dat->qmat[i][j] = sim_q;
			debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "\t%d: "
					"%d (%f:%f, %d)\n", j, sim_q, r, cdf,
					mod->n_quality);
		} else if (opt->simulation_random & QUALITIES
							&& mod->art_profile) {
			pmf = opt->modo->art_separate
				? &mod->art_profile[(j + j_off) * mod->n_hq_combos + hnuc * mod->n_quality]
				: &mod->art_profile[(j + j_off) * mod->n_hq_combos];

			r = rand() / (RAND_MAX + 1.);
			cdf = pmf[0];
			sim_q = 0;
			while (r > cdf) {
				cdf += pmf[++sim_q];
				if (sim_q == mod->n_quality)
					break;
			}

			dat->qmat[i][j] = sim_q;
			debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "\t%d: "
					"%d (%f:%f, %d)\n", j, sim_q, r, cdf,
					mod->n_quality);
		} else {
			sim_q = dat->qmat[i][j] + dat->fdata->min_quality;
		}

		if (sim_q + mod->min_quality < dat->fdata->min_quality)
			dat->fdata->min_quality = sim_q + mod->min_quality;
		if (sim_q + mod->min_quality > dat->fdata->max_quality)
			dat->fdata->max_quality = sim_q + mod->min_quality;

		/* given quality score, simulate nucleotide */
		r = rand() / (RAND_MAX + 1.);

		if (opt->modo->parameterization == QUALITY_PARAMETERIZATION
			|| opt->modo->parameterization == ART_PARAMETERIZATION) {

			/* quality score is literal */
			ep = error_prob(dat->fdata, sim_q);
			if (r <= ep) {	/* with error */
				r /= ep;
				sim_nuc = !hnuc;
				cdf = mod->gamma[hnuc * NUM_NUCLEOTIDES
								+ sim_nuc];
				while (r > cdf) {
					++sim_nuc;
					if (sim_nuc == hnuc)
						++sim_nuc;
					cdf += mod->gamma[hnuc
						* NUM_NUCLEOTIDES + sim_nuc];
					if (sim_nuc == NUM_NUCLEOTIDES - 1)
						break;
				}
			} else {	/* without error */
				sim_nuc = hnuc;
			}

			debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "\t%d: "
					"%d (%f:%f)\n", j, (int)sim_nuc, r, ep);
		} else if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION) {
			jhq_idx =  (j + j_off) * mod->n_hq_combos
						+ hnuc * mod->n_quality + sim_q;
			cdf = 0;
			for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b) {
				bjhq_idx = b * mod->n_jhq_combos + jhq_idx;
				cdf += mod->p_bjhq[bjhq_idx];
				if (r <= cdf) {
					sim_nuc = b;
					break;
				}
			}
//if (sim_nuc != hnuc) fprintf(stderr, "Simulate: %c from %c (%f)\n", nuc(dat->fdata, sim_nuc),  nuc(dat->fdata, hnuc), r);
		} else if (opt->modo->parameterization == DADA2_PARAMETERIZATION) {
			double vep[NUM_NUCLEOTIDES] = {0, 0, 0, 0};

			/* error profile was read in the order of A (0),C(1),G(3),T(2) */
			if (hnuc == XY_A) { //0
				vep[0] = mod->dada2_profile[sim_q]; //A->A
				vep[1] = mod->dada2_profile[mod->n_quality + sim_q];  //A->C
				vep[2] = mod->dada2_profile[mod->n_quality*3 + sim_q]; //A->T
				vep[3] = mod->dada2_profile[mod->n_quality*2 + sim_q]; //A->G
			} else if (hnuc == XY_C) {  //1
				vep[0] = mod->dada2_profile[mod->n_quality*4 + sim_q]; //C->A
				vep[1] = mod->dada2_profile[mod->n_quality*5 + sim_q];  //C->C
				vep[2] = mod->dada2_profile[mod->n_quality*7 + sim_q]; //C->T
				vep[3] = mod->dada2_profile[mod->n_quality*6 + sim_q]; //C->G
			} else if (hnuc == XY_G) {  //3
				vep[0] = mod->dada2_profile[mod->n_quality*8 + sim_q]; //G->A
				vep[1] = mod->dada2_profile[mod->n_quality*9 + sim_q];  //G->C
				vep[2] = mod->dada2_profile[mod->n_quality*11 + sim_q]; //G->T
				vep[3] = mod->dada2_profile[mod->n_quality*10 + sim_q]; //G->G
			} else {  //2, T
				vep[0] = mod->dada2_profile[mod->n_quality*12 + sim_q]; //T->A
				vep[1] = mod->dada2_profile[mod->n_quality*13 + sim_q];  //T->C
				vep[2] = mod->dada2_profile[mod->n_quality*15 + sim_q]; //T->T
				vep[3] = mod->dada2_profile[mod->n_quality*14 + sim_q]; //T->G
			}

			cdf = vep[0];
			sim_nuc = 0;
			while (r > cdf) {
				cdf += vep[++sim_nuc];
				if (sim_nuc == NUM_NUCLEOTIDES - 1)
					break;
			}
//if (sim_nuc != hnuc) fprintf(stderr, "Simulate: %c from %c (%f)\n", nuc(dat->fdata, sim_nuc),  nuc(dat->fdata, hnuc), r);
		}
		dat->dmat[i][j] = (data_t) sim_nuc;
	/* AMPLICLUST MODEL */
	} else if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION) {

		r = rand() / (RAND_MAX + 1.);

		/* no error */
		if (r <= mod->delta[j + j_off]) {

			dat->dmat[i][j] = hnuc;

			/* simulate quality score */
			r = rand() / (RAND_MAX + 1.);
			sim_q = 0;
			cdf = mod->lambda0[(j + j_off) * mod->n_quality];

			while (r > cdf) {
				cdf += mod->lambda0[(j + j_off) * mod->n_quality
								+ ++sim_q];
				if (sim_q == mod->n_quality - 1)
					break;
			}

			dat->qmat[i][j] = sim_q;

		/* error */
		} else {

			/* simulate nucleotide */
			r = rand() / (RAND_MAX + 1.);

			sim_nuc = !hnuc;
			cdf = mod->gamma[hnuc * NUM_NUCLEOTIDES + sim_nuc];
			while (r > cdf) {
				++sim_nuc;
				if (sim_nuc == hnuc)
					++sim_nuc;
				cdf += mod->gamma[hnuc * NUM_NUCLEOTIDES
								+ sim_nuc];
				if (sim_nuc == NUM_NUCLEOTIDES - 1)
					break;
			}

			dat->dmat[i][j] = sim_nuc;

			/* simulate quality score */
			r = rand() / (RAND_MAX + 1.);

			sim_q = 0;
			cdf = mod->lambda1[(j + j_off) * mod->n_quality];

			while (r > cdf)
				cdf += mod->lambda1[(j + j_off)
					* mod->n_quality + ++sim_q];

			dat->qmat[i][j] = sim_q;

		}

		if (sim_q + mod->min_quality < dat->fdata->min_quality)
			dat->fdata->min_quality = sim_q + mod->min_quality;
		if (sim_q + mod->min_quality > dat->fdata->max_quality)
			dat->fdata->max_quality = sim_q + mod->min_quality;

	} else {
		mmessage(ERROR_MSG, INTERNAL_ERROR, "Code not ready!\n");
	}
}/* simulate_position */


/**
 * Simulate one read with length len based on the current model.
 *
 *
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param len		simulate read of this length
 * @param simu_read	pointer to nucleotide sequences
 * @param simu_qual	pointer to quality score sequences
 * @param off		simulate read of this offset
 * @return		error status
 */
int simulate_read(data *dat, model *mod, unsigned int len,
	unsigned char *simu_read, unsigned char *simu_qual, unsigned int off)
{
	int fxn_debug = ABSOLUTE_SILENCE;

	if (!simu_qual || !simu_read)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"simulate read:space error");

	unsigned int true_clus_id = 0;
	double r, cdf;

	if (len > dat->max_read_length)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "simulate read:length error");

	/* choose cluster id read */
	r = rand() / (RAND_MAX + 1.);
	cdf = mod->best_pi[0];
	while (r > cdf){
		cdf += mod->best_pi[++true_clus_id];
		if (true_clus_id == (mod->K -1))
			break;
	}
	/* now simulate read and quality scores */
	/* the sequences use XY coding */

	for (size_t j = 0; j < len; ++j) {
		/* true nucleotide */
		char nuc = mod->best_haplotypes[
			true_clus_id * dat->max_read_length + j];

		/* read contains error */
		r = rand() / (RAND_MAX + 1.);
		if (r > mod->best_delta[j + off]) {
			/* choose new nucleotide */
			r = rand() / (RAND_MAX + 1.);
			cdf = mod->best_gamma[nuc * NUM_NUCLEOTIDES];
			char l = 0;
			while(r > cdf){
				cdf += mod->best_gamma[nuc * NUM_NUCLEOTIDES + ++l];
				if (l == (NUM_NUCLEOTIDES-1))
					break;
			}
			simu_read[j] = l;
		} else {
			simu_read[j] = nuc;
		}

		/* choose quality score */
		r = rand() / (RAND_MAX + 1.);

		/* no error */
		if (simu_read[j] == nuc) {
			char q = 0;
			cdf = mod->best_lambda0[(j + off) * mod->n_quality];
			while (r > cdf){
				cdf += mod->best_lambda0[mod->n_quality
					* (j + off) + ++q];
				if(q == mod->n_quality - 1)
					break;
			}
			simu_qual[j] = q;
		/* error */
		} else {
			char q = 0;
			//cdf = mod->best_lambda1_simple[0];
			cdf = mod->best_lambda1[(j + off) * mod->n_quality];
			while (r > cdf){
				cdf += mod->best_lambda1[mod->n_quality
					* (j + off) + ++q];
				if(q == mod->n_quality-1)
					break;
			}
				//cdf += mod->best_lambda1_simple[++q];
			simu_qual[j] = q;
		}
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "%s\n",
		display_sequence(simu_read, len, XY_ENCODING));

	return NO_ERROR;
} /* simulate_read */

/**
 * Simulate haplotypes for given ancester, rate, and (optionally) mutable
 * positions.
 *
 * [KSD] Thoroughly changed this code as it was not simulating under JC69
 * as advertised.
 *
 * [NOTE] If the model is changed to be nucleotide-specific (anything but
 * JC69), then we have to insure that it simulates xy_t type nucleotides.
 *
 * @param length	length of haplotypes
 * @param K		number of haplotypes
 * @param mut_position	mutable positions (1/0 indicator)
 * @param hap_spread	evolution rate, expected no. mutations / site
 * @param anc		pointer to ancestor sequence
 * @param hap		pointer to haplotypes, to be updated
 * @return		error status
 */
int simulate_haplotypes(unsigned int length, unsigned int K,
	unsigned int *mut_position, double hap_spread, unsigned char *anc,
	unsigned char *hap)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;// 
	int err = NO_ERROR;
	unsigned int k = 0;	/* current haplotype */
	unsigned int same;	/* test if proposed hap. matches prev. hap. */
	unsigned int itr = 0;	/* no. iterations to get distinct haplotype */
	unsigned int n_mutable_sites = 0;	/* number of mutatable sites */
	unsigned int max = 1000;/* max no attempts at unique hap. */
	double rate = 0.;	/* expected no. mutations per site */
	double r;
	size_t dim = NUM_NUCLEOTIDES * NUM_NUCLEOTIDES;
	double *a = malloc(dim * sizeof *a);
	double *Pt = NULL;	/* JC69 transition probability matrix */

	if (!a)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "rate matrix");

	if (mut_position)
		for (unsigned int j = 0; j < length; ++j) {
			if (mut_position[j])
				n_mutable_sites++;
		}
	else
		n_mutable_sites = length;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Ancestor: %s\n",
		display_sequence(anc, length, XY_ENCODING));

	/* set up JC69 rate matrix: rate is exp. no. changes / site */
	rate = hap_spread * length / n_mutable_sites;
	for (size_t i = 0; i < dim; ++i)
		a[i] = 0.25 * rate;	/* https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969) */
	for (size_t i = 0; i < NUM_NUCLEOTIDES; ++i)
		a[i * NUM_NUCLEOTIDES + i] = - 0.75 * rate;

	Pt = r8mat_expm1(NUM_NUCLEOTIDES, a);

	if (!Pt) {
		free(a);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "transition "
							"probability matrix");
	}

	while (k < K) {
		for (unsigned int j = 0; j < length; ++j) {
			if (mut_position && !mut_position[j]) {
				hap[k * length + j] = anc[j];
				continue;
			}
			r = rand() / (RAND_MAX + 1.);
			xy_t l = 0;
			for (double dsum = Pt[anc[j] * NUM_NUCLEOTIDES];
				dsum < r; dsum += Pt[anc[j]
					* NUM_NUCLEOTIDES + ++l]);
			hap[k * length + j] = l;
		}

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "%s\n",
			display_sequence(&hap[k * length], length, XY_ENCODING));//

		/* check if the haplotype is identical with previous haplotypes */
		same = 0;
		for (unsigned int h = 0; h < k; ++h)
			if (hamming_char_dis((char *) &hap[h*length],
					(char *) &hap[k*length], length) == 0) {
				same = 1;
				break;
			}

		if (!same)
			k++;

		if (++itr > max) {
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Sorry, "
				"cannot simulate %u distinct haplotypes with "
				"current settings.  See -c.\n", K);
			break;
		}
	}

	free(a);
	free(Pt);

	return err;
}/* simulate_haplotypes */

/**
 * Simulate a sequence under assumption iid Uniform(A, C, G, T) (XY_ENCODING)
 *
 * [KSD] Was pointer bug.
 *
 * @param seq	pointer to the sequence
 * @param len	length of the sequence
 * @return	error status
 */
int simulate_ancestor(unsigned char **seq, unsigned int len)
{
	*seq = malloc(len * sizeof **seq);
	if (!*seq)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"options.ancestor");
	for (size_t j = 0; j < len; ++j) {
		double r = rand() / (RAND_MAX + 1.);
		if (r <= 0.25)
			(*seq)[j] = XY_A;
		else if (r <= 0.5)
			(*seq)[j] = XY_C;
		else if (r <= 0.75)
			(*seq)[j] = XY_T;
		else
			(*seq)[j] = XY_G;
	}
	return NO_ERROR;
}/* simulate_ancestor */


/**
 * Initialize estimation model from true partition.
 *
 * @param mod	pointer to model object for estimation
 * @param dat	pointer to data object
 * @param ini	pointer to initializer object for estimation
 * @param sim	pointer to simulator object
 * @param sopt	pointer to simulation options object
 * @param mopt	pointer to model options object
 * @return	error status
 */
int initialize_from_true_partition(model *mod, data *dat, initializer *ini,
			simulator *sim, simulate_options *sopt, initialize_options *iopt)
{
	int err = NO_ERROR;

	/* copy true partition to initializer object */
	memcpy(ini->best_cluster_id, sim->true_cluster_id,
		dat->sample_size * sizeof *ini->cluster_id);
	for (unsigned int k = 0; k < sopt->true_K; ++k)
		memcpy(ini->best_modes[k], &sim->mod->haplotypes[dat->max_read_length * k],
			dat->max_read_length * sizeof *sim->mod->haplotypes);

	/* initialize parameters of estimation model */
	err = initialize_model_parameters(mod, dat, iopt, ini, BEST_VALUES);
	param_update(mod, dat, iopt->modo, DEFAULT_VALUES);

	return err;
} /* initialize_from_true_partition */


/**
 * Initialize from true parameters.
 *
 * @param mod		estimation model object pointer
 * @param sim_mod	simulation model object pointer
 * @param dat		data object pointer
 * @param mopt		estimation model options
 * @param sopt		simulate_options object
 * @return		error status
 */
int initialize_from_true_parameters(model *mod, model *sim_mod, data *dat,
		model_options *mopt, simulate_options *sopt)
{
	if (mopt->K != sopt->true_K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot initialize"
			" with true parameters if simulation model does not "
			"match estimation model: see -k, -sk\n");

	if (mopt->background_model != sopt->modo->background_model)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot initialize"
			" with true parameters when simulation does not "
			"match estimation model: see -b\n");

	if (mopt->parameterization != sopt->modo->parameterization)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot initialize"
			" with true parameters when simulation model does not "
			"match estimation model: see -p, -sp\n");

	memcpy(mod->haplotypes, sim_mod->haplotypes,
		dat->max_read_length * sopt->true_K * sizeof *mod->haplotypes);

	memcpy(mod->pi, sim_mod->pi, mod->n_mix * sizeof *mod->pi);

	if (mopt->parameterization == MLOGIT_PARAMETERIZATION) {

		memcpy(mod->beta, sim_mod->beta, mod->n_beta_coef * sizeof *mod->beta);

	} else if (mopt->parameterization == QUALITY_PARAMETERIZATION
		|| mopt->parameterization == DADA2_PARAMETERIZATION
		|| mopt->parameterization == ART_PARAMETERIZATION) {

		memcpy(mod->gamma, sim_mod->gamma, NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof *mod->gamma);

	} else {
		memcpy(mod->delta, sim_mod->delta, dat->max_read_position * sizeof *mod->pi);
		memcpy(mod->gamma, sim_mod->gamma, NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof *mod->gamma);
		memcpy(mod->lambda0, sim_mod->lambda0, dat->max_read_position * mod->n_quality * sizeof *mod->lambda0);
		memcpy(mod->lambda1, sim_mod->lambda1, dat->max_read_position * mod->n_quality * sizeof *mod->lambda1);
	}

	return NO_ERROR;
} /* initialize_from_true_parameters */
