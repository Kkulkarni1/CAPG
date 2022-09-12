/**
 * @file model.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 *
 * Manipulate model object.
 *
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#include "model.h"
#include "ampliclust.h"
#include "data.h"
#include "model_options.h"
#include "simulate.h"
#include "io.h"
#include "fastq.h"
#include "statistics.h"
#include "math.h"

const unsigned int NUM_NUCLEOTIDES_SQUARED = NUM_NUCLEOTIDES * NUM_NUCLEOTIDES;

int read_art_profile(model *mod, model_options *opt, data *dat);
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);

/**
 * Create model object.
 *
 * @param mod	model object to create
 * @param dat	pointer to data object
 * @param opt	pointer to model_options object
 * @param delay	delay sync
 * @return	error status
 */
int make_model(model **mod, data *dat, model_options *opt, int delay)
{
	int err = NO_ERROR;
	model *rm;
	*mod = malloc(sizeof **mod);

	if (*mod == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model");

	rm = *mod;

	rm->synced = 0;
	rm->mopt = opt;

	/* K */
	rm->n_mix = opt->K + opt->background_model;
	rm->K = opt->K;

	/* pi */
	rm->pi = malloc(rm->n_mix * sizeof *rm->pi);
	rm->npi = malloc(rm->n_mix * sizeof *rm->npi);
	rm->best_pi = malloc(rm->n_mix * sizeof *rm->best_pi);

	if (!rm->pi || !rm->npi || !rm->best_pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model.pi");

	/* gamma */
	if (opt->parameterization == DEFAULT_PARAMETERIZATION
		|| opt->parameterization == QUALITY_PARAMETERIZATION
		|| opt->parameterization == ART_PARAMETERIZATION) {
		rm->gamma = malloc(NUM_NUCLEOTIDES_SQUARED * sizeof *rm->gamma);
		rm->ngamma = malloc(NUM_NUCLEOTIDES_SQUARED
							* sizeof *rm->ngamma);
		rm->best_gamma = malloc(NUM_NUCLEOTIDES_SQUARED
						* sizeof *rm->best_gamma);

		if (!rm->gamma || !rm->ngamma || !rm->best_gamma)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model.gamma");

	} else {
		rm->gamma = rm->ngamma = rm->best_gamma = NULL;
	}


	/* delta */
	if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
		rm->delta = malloc(dat->max_read_position
							* sizeof *rm->delta);
		rm->ndelta = malloc(dat->max_read_position
							* sizeof *rm->ndelta);
		rm->best_delta = malloc(dat->max_read_position
							* sizeof *rm->delta);

		if (!rm->delta || !rm->ndelta || !rm->best_delta)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model.delta");
	} else {
		rm->delta = rm->ndelta = rm->best_delta = NULL;
	}


	/* pi in the background model */
	if (opt->background_model) {
		rm->bg_pi = malloc(NUM_NUCLEOTIDES * sizeof *rm->bg_pi);
		rm->nbg_pi = malloc(NUM_NUCLEOTIDES * sizeof *rm->bg_pi);
		rm->best_bg_pi = malloc(NUM_NUCLEOTIDES * sizeof *rm->bg_pi);

		if (!rm->bg_pi || !rm->nbg_pi || !rm->best_bg_pi)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"model.bg_pi");

	} else {
		rm->bg_pi = rm->nbg_pi = rm->best_bg_pi = NULL;
	}

	/* An attempt at reducing model complexity by using
	 * a coarser discretization on quality scores.
	 * [TODO] If this is revived, need to check that data is loaded; it will
	 * not be for some simulation cases.  Check data::fdata->empty
	 */
	/*
	if (opt->q_model == FIVE_BY_SIX) {
		for (size_t i = 0; i < dat->sample_size; ++i)
			for (size_t j = 0; j < dat->lengths[i]; ++j) {
				unsigned char *dptr = &dat->qmat[i][j];
				if (*dptr < 31)
					*dptr = *dptr / 5;
				else
					*dptr = *dptr - 24;
			}

		rm->n_quality = 6 + (dat->n_quality > 30
			? dat->n_quality - 30 : 0);
	} else
		rm->n_quality = dat->n_quality;
	*/

	/* by default: quality ranges from data */
	rm->n_quality = dat->n_quality;
	rm->min_quality = dat->fdata->min_quality;
	rm->max_quality = dat->fdata->max_quality;

//	rm->parameterization = opt->parameterization;
	rm->lambda0 = NULL;
	rm->nlambda0 = NULL;
	rm->best_lambda0 = NULL;

	rm->lambda1 = NULL;
	rm->nlambda1 = NULL;
	rm->best_lambda1 = NULL;

	rm->beta = NULL;
	rm->nbeta = NULL;
	rm->best_beta = NULL;
	rm->beta_diff = NULL;
	rm->p_bjhq = NULL;
//	rm->p_jhq = NULL;
	rm->e_bjhq = NULL;
	rm->e_jhq = NULL;
	rm->gradient = NULL;
	rm->hessian = NULL;
	rm->px = NULL;

	/* malloc space for both lambda0 and lambda1 */
	if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
		rm->lambda0 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->lambda0);
		rm->nlambda0 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->nlambda0);
		rm->best_lambda0 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->best_lambda0);

		rm->lambda1 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->lambda1);
		rm->nlambda1 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->nlambda1);
		rm->best_lambda1 = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->best_lambda1);

		if (!rm->lambda0 || !rm->nlambda0 || !rm->best_lambda0
			|| !rm->lambda1 || !rm->nlambda1 || !rm->best_lambda1)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::lambda");
	}

	/* background model */
	if (opt->background_model) {
		rm->bg_lambda = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->bg_lambda);
		rm->nbg_lambda = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->bg_lambda);
		rm->best_bg_lambda = malloc(dat->max_read_position
			* rm->n_quality * sizeof *rm->bg_lambda);

		if (!rm->bg_lambda || !rm->nbg_lambda || !rm->best_bg_lambda)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"model.bg_lambda");

	}

	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
		/* WORKING */
		rm->eq_inter_q_power = 1;	/* error * quality */
		rm->ep_inter_p_power = 0;	/* error * position */
//		rm->eh_inter = 0;		/* error * haplotype */
		rm->position_power = 0;		/* hard-coded */
		rm->quality_power = 0;		/* polynomial */
		rm->pq_inter_p_power = 0;	/* model */
		rm->pq_inter_q_power = 0;	/* [TODO] allow CL control */
		rm->n_predictors = 
			1			/* intercept */
			+ NUM_NUCLEOTIDES - 1	/* source nucleotide */
			+ rm->position_power 	/* terms for position */
			+ rm->quality_power	/* terms for quality */
			/* terms for interactions */
			+ rm->eq_inter_q_power + rm->ep_inter_p_power // + rm->eh_inter
			+ rm->pq_inter_p_power * rm->pq_inter_q_power;

		/* one set for each outcome nucleotide, except first:
		 * order of predictors: haplotype source nucleotide, 
		 * position powers, quality powers, interactions, ...
		 */
		rm->n_beta_coef = rm->n_predictors * (NUM_NUCLEOTIDES - 1);
		rm->beta = malloc(rm->n_beta_coef * sizeof *rm->beta);
		rm->nbeta = malloc(rm->n_beta_coef * sizeof *rm->beta);
		rm->best_beta = malloc(rm->n_beta_coef * sizeof *rm->beta);
		rm->beta_diff = malloc(rm->n_beta_coef * sizeof *rm->beta);

		if (!rm->beta || !rm->nbeta || !rm->best_beta || !rm->beta_diff)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model::beta");

		/* initialize beta with no effects */
		/* [KSD, TODO] setting true value here for simulation */
		for (unsigned int i = 0; i < rm->n_beta_coef; ++i) {
			rm->beta[i] = 0;

/*
			if (!(i % rm->n_predictors))
				rm->beta[i] = -8;
			else if ((i % rm->n_predictors) == (i + rm->n_predictors) / rm->n_predictors)
				rm->beta[i] = 8;
			if ((i % rm->n_predictors) >= rm->n_predictors)
				rm->beta[i] = 1;
*/
		}

		/* transition probabilities are available for every possible
		 * position, quality, haplotype nucleotide, and read nucleotide
		 * combination, for easy lookup in uncompressed read data
		 */

		/* read nucleotide, site, quality, haplotype combos */
		rm->n_satisfaction_indices = rm->n_quality
			* dat->max_read_position * NUM_NUCLEOTIDES_SQUARED;

		/* site, quality score, and haplotype nucleotide combos */
		rm->n_jhq_combos = dat->max_read_position * NUM_NUCLEOTIDES
							* rm->n_quality;

		rm->p_bjhq = malloc(rm->n_satisfaction_indices
						* sizeof *rm->p_bjhq);

//		rm->p_jhq = malloc(rm->n_jhq_combos * sizeof *rm->p_jhq);
		rm->max_sat_jhq = malloc(rm->n_jhq_combos * sizeof
							*rm->max_sat_jhq);

		if (!rm->p_bjhq || !rm->max_sat_jhq)// || !rm->p_jhq)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::p_bjhq");

		/* quality score, haplotype nucleotide combos */
		rm->n_hq_combos = NUM_NUCLEOTIDES * rm->n_quality;

		/* site, quality score, haplotype nucleotide, predictor */
		rm->n_jhqp_combos = dat->max_read_position * NUM_NUCLEOTIDES
			* rm->n_quality * rm->n_predictors;

		/* the expected number of transitions similarly available */
		rm->e_bjhq = malloc(rm->n_satisfaction_indices
							* sizeof *rm->e_bjhq);
		rm->e_jhq = malloc(rm->n_jhq_combos * sizeof *rm->e_jhq);

		if (!rm->e_bjhq || !rm->e_jhq)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::e_bjhq");

                /* initialize to 0 */
                for (unsigned int j = 0; j < dat->max_read_position; ++j) {
			unsigned int jhq_index = j * rm->n_hq_combos;
			for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
				jhq_index += h * rm->n_quality;
				for (unsigned int q = 0; q < rm->n_quality; ++q)
					rm->e_jhq[jhq_index + q] = 0;
				jhq_index -= h * rm->n_quality;
			}
		}

		rm->gradient = malloc(rm->n_beta_coef * sizeof *rm->gradient);
		rm->hessian = malloc(rm->n_beta_coef * rm->n_beta_coef
							* sizeof *rm->hessian);
		if (!rm->gradient || !rm->hessian)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"model::gradient or model::hessian");

		rm->px = malloc(rm->n_beta_coef * sizeof *rm->px);

		if (!rm->px)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::px");

	} else if (opt->parameterization == ART_PARAMETERIZATION) {
		if (opt->art_separate)
			rm->n_hq_combos = NUM_NUCLEOTIDES * rm->n_quality;
		else
			rm->n_hq_combos = rm->n_quality;
	}

	/* haplotypes */
	rm->haplotypes = malloc(dat->max_read_length * opt->K
		* sizeof *rm->haplotypes);
	rm->nhaplotypes = malloc(dat->max_read_length * opt->K
		* sizeof *rm->nhaplotypes);
	rm->best_haplotypes = malloc(dat->max_read_length * opt->K
		* sizeof *rm->best_haplotypes);

	if (!rm->haplotypes || !rm->nhaplotypes || !rm->best_haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"model::haplotypes");

	/* eik: allocate and initialize to all 0 */
	rm->eik = malloc(dat->sample_size * rm->n_mix * sizeof *rm->eik);

	if (rm->eik == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::eik");


	/* haplotypes simulation control */
	rm->mut_position = NULL;

	rm->ll = 0;
	rm->best_ll = opt->previous_ll;	/* best AECM log likelihood */
	rm->best_init_ll = -INFINITY;	/* best RND-EM criterion */
	rm->JC_ll = -INFINITY;

	rm->aic = INFINITY;
	rm->bic = INFINITY;

	/* JC69 approximate prior on haplotypes */
	rm->distance = NULL;
	rm->est_ancestor = NULL;

	if (opt->JC69_model) {
		rm->distance = malloc(rm->K * sizeof *rm->distance);
		rm->est_ancestor = malloc(dat->max_read_length
					* sizeof * rm->est_ancestor);

		if (!rm->distance || !rm->est_ancestor)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model::JC69");
	}

	rm->qual_list = NULL;
	rm->qual_per_list = NULL;

	if ((!dat->fdata->empty || delay)
		&& (err = sync_model(rm, dat, opt)))
		return err;

        rm->dada2_profile = NULL;

	/* allocate space for DADA2 profile; overwrite model::n_quality */
	if (opt->parameterization == DADA2_PARAMETERIZATION
		&& (err = read_dada2_error_profile(rm, dat, opt)))
			return err;

	rm->art_profile = NULL;

	/* allocate space for ART profile; overwrite model::n_quality */
	if (opt->parameterization == ART_PARAMETERIZATION
		&& (err = read_art_profile(rm, opt, dat)))
		return err;

	return err;
} /* make_model */


/**
 * Record information about quality scores observed at each site if mlogit
 * model selected.  Called once, not efficient.  Also count number of
 * parameters.
 *
 * Since most quality scores are not observed at every site, it is inefficient
 * to compute expected transitions of mlogit model for every possible
 * combination of site, quality score, and haplotype nucleotide.  This code
 * records a shortened list of quality scores observed in each site of
 * the entire dataset (assuming the data object has not subsampled the
 * fastq file).  Then, expected transitions will be computed only for
 * these observed combinations.
 *
 * [NOTE] Use data::n_quality here, not model::n_quality.  This is an accounting
 * of quality scores observed in the data, not in the model.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to model_options object
 * @return	error status
 */
int sync_model(model *mod, data *dat, model_options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;
	int err = NO_ERROR;

	if (mod->qual_list)
		free(mod->qual_list);
	if (mod->qual_per_list)
		free(mod->qual_per_list);

	/* allocate space enough for all quality scores observed at all sites */
	mod->qual_list = calloc(dat->max_read_length * dat->n_quality,
						sizeof *mod->qual_list);
	mod->qual_per_list = malloc(dat->max_read_length
						* sizeof *mod->qual_per_list);
	if (!mod->qual_list || !mod->qual_per_list) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::qual_list");
		goto SYNC_MODEL_EXIT;
	}

	/* for each read */
	unsigned char *qptr = dat->fdata->quals;
	for (unsigned int i = 0; i < dat->fdata->n_reads; ++i) {

		unsigned int off = dat->offset ? dat->offset[i] : 0;
		unsigned int len = read_length(dat->fdata, i);
		/* for each read position */
		for (unsigned int j = 0; j < len; ++j) {
			mod->qual_list[(j + off) * dat->n_quality + *qptr] = 1;
			++qptr;
		}
	}

	/* actual number of observed combinations of site and quality score */
	unsigned int n_jq_combos = 0;

	for (unsigned int j = 0; j < dat->max_read_length; ++j) {

		/* count and record quality scores used */
		unsigned n_used_qual = 0;
		for (unsigned char q = 0; q < dat->n_quality; ++q)
			if (mod->qual_list[j * dat->n_quality + q])
				mod->qual_list[n_jq_combos + n_used_qual++] = q;
		mod->qual_per_list[j] = n_used_qual;
		n_jq_combos += n_used_qual;
	}

	/* reallocate only necessary data */
	unsigned char *qtmp = realloc(mod->qual_list, n_jq_combos
						* sizeof *mod->qual_list);
	if (!qtmp) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::qual_list");
		goto SYNC_MODEL_EXIT;
	}
	mod->qual_list = qtmp;

	if (fxn_debug) {
		unsigned int l = 0;
		for (unsigned int j = 0; j < dat->max_read_length; ++j) {
/*
			debug_call(DEBUG_I <= fxn_debug, fxn_debug,
				fprint_uchars(stderr, &mod->qual_list[l], 
				mod->qual_per_list[j], 2, 1));
*/
			l += mod->qual_per_list[j];
		}
	}

	/* count parameters [TODO] update for new models */
	mod->n_param = opt->K * dat->max_read_length	/* haplotypes */
		+ opt->K - 1				/* <- pi, \/ gamma */
		+ (opt->parameterization != MLOGIT_PARAMETERIZATION && opt->parameterization != DADA2_PARAMETERIZATION ? 12 : 0)
		+ (opt->background_model ? NUM_NUCLEOTIDES - 1 : 0);/* bg_pi */

	/* nucleotide emission parameters */
	mod->n_param_lambda = 0;

	/* Multinomial logistic regression model has already counted the
	 * number of parameters above.
	 */
	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
		mod->n_param_lambda = mod->n_beta_coef;

		if (opt->background_model)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
							"not implemented yet");

	/* Quality model has one delta parameter for each read position and
	 * one emission probability for every observed quality at each
	 * position.
	 */
	} else if (opt->parameterization == DEFAULT_PARAMETERIZATION
		&& opt->model_quality) {

		mod->n_param_lambda += dat->max_read_position;	/* delta */

		mod->n_param_lambda += 2 * n_jq_combos;	/* lambda0, lambda1 */

		/* background model uses one lambda, as error status unknown */
		if (opt->background_model)
			mod->n_param_lambda += n_jq_combos;

	}

	mod->n_param += mod->n_param_lambda;	/* lambda */

	if (opt->parameterization != MLOGIT_PARAMETERIZATION) {
		free(mod->qual_per_list);
		mod->qual_per_list = NULL;
		free(mod->qual_list);
		mod->qual_list = NULL;
	}

	mmessage(INFO_MSG, NO_ERROR, "Number of parameters: %u\n",
							mod->n_param);

	mod->synced = 1;

	return err;

SYNC_MODEL_EXIT:
	if (mod->qual_per_list)
		free(mod->qual_per_list);
	if (mod->qual_list)
		free(mod->qual_list);

	return err;
} /* sync_model */

/**
 * Reallocate quality-related values.  When the number of distinct quality
 * values in the data changes, we need to update parameters that range over
 * all those values.
 *
 * @param mod	pointer to model object to update
 * @param dat	pointer to data object
 * @param mopt	pointer model_options object
 * @return	error status
 */
int realloc_quality_information(model *mod, data *dat, model_options *mopt, unsigned int nquality)
{
	double *nptr = NULL;
	mod->n_quality = nquality;

	/* malloc space for both lambda0 and lambda1 */
	if (mopt->parameterization == DEFAULT_PARAMETERIZATION) {
		nptr = realloc(mod->lambda0, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::lambda0");
		mod->lambda0 = nptr;
		nptr = realloc(mod->nlambda0, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::nlambda0");
		mod->nlambda0 = nptr;
		nptr = realloc(mod->best_lambda0, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::best_lambda0");
		mod->best_lambda0 = nptr;

		nptr = realloc(mod->lambda1, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::lambda1");
		mod->lambda1 = nptr;
		nptr = realloc(mod->nlambda1, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::nlambda1");
		mod->nlambda1 = nptr;
		nptr = realloc(mod->best_lambda1, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::best_lambda1");
		mod->best_lambda1 = nptr;
	}

	/* background model */
	if (mopt->background_model) {
		nptr = realloc(mod->bg_lambda, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::bg_lambda");
		mod->bg_lambda = nptr;
		nptr = realloc(mod->nbg_lambda, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::nbg_lambda");
		mod->nbg_lambda = nptr;
		nptr = realloc(mod->best_bg_lambda, dat->max_read_position
							* mod->n_quality);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::best_bg_lambda");
		mod->best_bg_lambda = nptr;
	}

	if (mopt->parameterization == MLOGIT_PARAMETERIZATION) {
		/* read nucleotide, site, quality, haplotype combos */
		mod->n_satisfaction_indices = mod->n_quality
			* dat->max_read_position * NUM_NUCLEOTIDES_SQUARED;

		/* site, quality score, and haplotype nucleotide combos */
		mod->n_jhq_combos = dat->max_read_position * NUM_NUCLEOTIDES
							* mod->n_quality;

		nptr = realloc(mod->p_bjhq, mod->n_satisfaction_indices);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::p_bjhq");
		mod->p_bjhq = nptr;

		nptr = realloc(mod->max_sat_jhq, mod->n_jhq_combos);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::max_sat_jhq");
		mod->max_sat_jhq = nptr;

		/* quality score, haplotype nucleotide combos */
		mod->n_hq_combos = NUM_NUCLEOTIDES * mod->n_quality;

		/* site, quality score, haplotype nucleotide, predictor */
		mod->n_jhqp_combos = dat->max_read_position * NUM_NUCLEOTIDES
			* mod->n_quality * mod->n_predictors;

		/* the expected number of transitions similarly available */
		nptr = realloc(mod->e_bjhq, mod->n_satisfaction_indices);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::e_bjhq");
		mod->e_bjhq = nptr;

		nptr = realloc(mod->e_jhq, mod->n_jhq_combos);
		if (!nptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model::e_jhq");
		mod->e_jhq = nptr;

                /* initialize to 0 */
                for (unsigned int j = 0; j < dat->max_read_position; ++j) {
			unsigned int jhq_index = j * mod->n_hq_combos;
			for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
				jhq_index += h * mod->n_quality;
				for (unsigned int q = 0; q < mod->n_quality; ++q)
					mod->e_jhq[jhq_index + q] = 0;
				jhq_index -= h * mod->n_quality;
			}
		}
	} else if (mopt->parameterization == ART_PARAMETERIZATION) {
		if (mopt->art_separate)
			mod->n_hq_combos = NUM_NUCLEOTIDES * mod->n_quality;
		else
			mod->n_hq_combos = mod->n_quality;
	}

	return NO_ERROR;
} /* realloc_quality_information */


/**
 * Read ART profile.
 *
 * @param mod	pointer to model object
 * @param opt	pointer to model_options object
 * @param dat	pointer to data object
 * @return	error status
 */
int read_art_profile(model *mod, model_options *opt, data *dat)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int err = NO_ERROR;
	unsigned char *observed_quals = NULL;
	double sum;
	int data_line = 0;
	unsigned int j, max_j = 0, cnt, l;
	unsigned char q, b = 0;
	unsigned char min_q = MAX_QUALITY_SCORE, max_q = MIN_QUALITY_SCORE;
	unsigned int jh_idx;
	char c;
	FILE *fp = fopen(opt->art_file, "r");

	if (!fp)
		return mmessage(ERROR_MSG, FILE_NOT_FOUND, opt->art_file);

	/* find min and max quality score */
	do {
		fgetc(fp);
		fscanf(fp, "%u", &j);

		/* fast forward through positions we don't need */
		if (j >= dat->max_read_length)
			do {
				c = fgetc(fp);
			} while (!feof(fp) && c != '\n');

		/* read quality scores */
		while (!feof(fp) && fscanf(fp, "%hhu", &q)) {
			if (q < min_q)
				min_q = q;
			if (q > max_q)
				max_q = q;
		}

		/* skip data line */
		do {
			c = fgetc(fp);
		} while (!feof(fp) && c != '\n');

	} while (!feof(fp));

	if (j + 1 < dat->max_read_length)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Error profile is "
			"for reads of length %u, but your data have reads up "
			"to %ubp long.\n", j + 1, dat->max_read_length);

	if (dat->n_quality > max_q - min_q + 1)
		mmessage(WARNING_MSG, INTERNAL_ERROR, "Error profile "
			"qualities (%u) do not match data file qualities "
			"(%u).\n", max_q - min_q + 1, dat->n_quality);

	/* update model::n_quality */
	if (mod->n_quality != max_q + 1 && (err =
		realloc_quality_information(mod, dat, opt, max_q + 1)))
		return err;

	/* assumption about art error profiles */
	mod->min_quality = MIN_ASCII_QUALITY_SCORE;
	mod->max_quality = mod->min_quality + max_q;
	
	observed_quals = malloc(mod->n_quality * sizeof *observed_quals);

	if (!observed_quals)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "observed_quals");
	
	if (opt->art_separate)
		mod->art_profile = malloc(NUM_NUCLEOTIDES * mod->n_quality *
			dat->max_read_length * sizeof *mod->art_profile);
	else
		mod->art_profile = malloc(dat->max_read_length *
			mod->n_quality * sizeof *mod->art_profile);

	if (!mod->art_profile) {
		free(observed_quals);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model:art_profile");
	}


	rewind(fp);

	do {
		c = fgetc(fp);		/* source haplotype */
		fscanf(fp, "%u", &j);	/* site */

		if (max_j < j)
			max_j = j;

		/* fast forward through positions we don't need */
		if (j >= dat->max_read_length)
			do {
				c = fgetc(fp);
			} while (!feof(fp) && c != '\n');

		/* fast forward through uninteresting lines */
		if ((opt->art_separate && c == '.')
				|| (!opt->art_separate && c != '.'))
			do {
				c = fgetc(fp);
			} while (!feof(fp) && c != '\n');

		if (c == '\n')
			continue;

		/* ignore cases where true base is N */
		if (c == 'N')
			break;
		else if (c == '.')
			b = 0;
		else if (c == 'A')
			b = XY_A;
		else if (c == 'C')
			b = XY_C;
		else if (c == 'G')
			b = XY_G;
		else if (c == 'T')
			b = XY_T;

		/* read in distribution */
		if (data_line) {
			l = 0;
			jh_idx = j * mod->n_hq_combos + b * mod->n_quality;
			for (unsigned char q = 0; q < mod->n_quality; ++q)
				mod->art_profile[jh_idx + q] = 0;
			sum = 0;
			while (fscanf(fp, "%u", &cnt)) {
				mod->art_profile[jh_idx + 
						observed_quals[l++]] = cnt;
				sum += cnt;
			}
			for (unsigned char q = 0; q < mod->n_quality; ++q) {
				mod->art_profile[jh_idx + q] /= sum;
				debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "art_profile[%u][%u][%u] = %f (%d)\n", mod->art_profile[jh_idx + q], j, b, q, opt->art_separate);
			}
			data_line = 0;

		/* read in observed quality scores */
		} else {
			l = 0;
			while (!feof(fp) && fscanf(fp, "%hhu", &q)) {
				if (q > MAX_QUALITY_SCORE)
					return mmessage(ERROR_MSG,
						INTERNAL_ERROR, "ART error "
						"profile has quality scores "
						"above %u: contact "
						"programmer\n",
						MAX_QUALITY_SCORE);
				observed_quals[l++] = q;
			}
			data_line = 1;
		}
	} while (!feof(fp));

	if (max_j + 1 < dat->max_read_length)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "ART profile is "
			"only for read lenths of %u or less, and you have "
			"requested %u.\n", max_j + 1, dat->max_read_length);

	free(observed_quals);

	return err;
} /* read_art_profile */


/**
 * Calculate aic and bic modified by approximate JC69 hierarchical model on
 * haplotypes.
 *
 * @param hap			haplotypes
 * @param est_anc		ancester sequence
 * @param distance		distance from haplotypes to ancestor sequence
 * @param best_ll		current log likelihood from data
 * @param K			number of haplotypes
 * @oaram JC_ll			log likelihood from JC69 model
 * @param n_aic			pointer to aic, value updated
 * @param n_bic			pointer to bic, value updated
 * @param n_param		number of parameters in current model
 * @param max_read_length	length of haplotypes
 * @param sample_size		sample size
 *
 * return			err status
 **/
int modified_ic(unsigned char *hap, unsigned char *est_anc, double *distance,
	double best_ll, unsigned int K, double *JC_ll, double *n_aic,
	double *n_bic, unsigned int n_param, unsigned int max_read_length,
	size_t sample_size)
{
	int param_change = 0;

	m_JC69(hap, est_anc, distance, K, max_read_length);
	*JC_ll = e_JC69(hap, est_anc, distance, K, max_read_length);

	/* K branch lengths, ancestral haplotype, but no haplotypes estimated */
	param_change = K - max_read_length * (K - 1);

	*n_aic = aic(best_ll + *JC_ll, n_param + param_change);
	*n_bic = bic(best_ll + *JC_ll, n_param + param_change, sample_size);

	return NO_ERROR;
}/* modified_ic */


/**
 * get MMEs of all parameters in the JC69 model. 
 * 
 * @param hap	pointer to haplotypes
 * @param anc	ancestor haplotype, to be calculated (tbc)
 * @param dist	expected no. changes/site b/w ancestor & haplotype, tbc
 * @param K	number of haplotypes
 * @param len	length of reads
 * @return err
 **/
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	unsigned int K, unsigned int len)
{
	
	int err = NO_ERROR;
	unsigned int count[NUM_NUCLEOTIDES];
	unsigned int max_count;

	/* most common nucleotide across haplotypes is the estimated ancestor */
	for (unsigned int j = 0; j < len; j ++){
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; n++)
			count[n] = 0;
		for (unsigned int k = 0; k < K; k++)
			count[hap[k * len + j]]++;
		max_count = 0;
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; ++n)
			if (count[n] > max_count) {
				max_count = count[n];
				anc[j] = n; 
			}
	}

	/* [KSD, BUG] Was the expected number of OBSERVED changes per site. */
	/* estimate the expected number of changes per site of all haplotypes */
	for (unsigned int k = 0; k < K; k++) {
		if (0) {	/* bug-free version */
			double tmp = (double) hamming_char_dis( (char *)
				&hap[k*len], (char *) anc, (size_t) len) / len;
			dist[k] = -0.75 * log(1 - tmp / 0.75);
		} else {	/* buggy version */
			dist[k] = (double) hamming_char_dis( (char *)
				&hap[k*len], (char *) anc, (size_t) len) / len;
		}
	}

	return err;
}/* m_JC69 */


/**
 * Calculate the log likelihood of all K haplotypes under JC69 model.
 * 
 * @param hap   haplotype sequences
 * @param anc   ancestor sequence
 * @param dis	expected no. changes/site for each haplotype
 * @param K	number of haplotypes
 * @param len	length of the haplotypes
 * @return	err status
 **/
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	unsigned int K, unsigned int len)
{
	double ll = 0;

	for (unsigned int k = 0; k < K; k++)
		for (unsigned int j = 0; j < len; j++)
			if (anc[j] == hap[k * len + j])
				ll += log(0.25 + 0.75 * exp(-dist[k] / 0.75));
			else
				ll += log(0.25 - 0.25 * exp(-dist[k] / 0.75));

	return ll;
}/* e_JC69 */


/**
 * Compute parametric transition probabilities given coefficients beta.
 *
 * The order of predictors is: intercept (0), haplotype (source) nucleotide
 * (1-3), position powers (4-), quality powers, interactions...
 * beta coefficients grouped by read nucleotide.
 *
 * The transition probabilities (model::trans_probs) are ordered by read
 * nucleotide, position, haplotype nucleotide, quality, and includes all
 * possibilities (unlike beta) so lookup is easy.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param type	which beta coefficients (DEFAULT_VALUES|BEST_VALUES)
 */
void compute_transition_probabilities(model *mod, data *dat, int type)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_II;//DEBUG_I;//
#endif
	double sat_bjh;	/* haplotype-, position-, and read-related satisfaction */
	double sat_q;	/* quality-related satisfaction */
	double *bptr = type == BEST_VALUES ?  mod->best_beta : mod->beta;
	double site_power, qual_power;
	double sum;
	unsigned int bjhq_idx = 0, jhq_idx = 0;
#ifdef DEBUGGING_ON
	debug_call(fxn_debug >= DEBUG_I, fxn_debug, mmessage(INFO_MSG, NO_ERROR,
								"beta: "));
	debug_call(fxn_debug >= DEBUG_I, fxn_debug, fprint_doubles(stderr, bptr,
						mod->n_beta_coef, 3, 1));
#endif
	/* the constraint implied for transition to A */
	for (unsigned int j = 0; j < dat->max_read_position; ++j) {
		for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
			for (unsigned char q = 0; q < mod->n_quality; ++q) {
				jhq_idx = j * mod->n_hq_combos + h * mod->n_quality + q;
				mod->p_bjhq[jhq_idx] = 1;
				//mod->p_jhq[jhq_idx] = 1;
				mod->max_sat_jhq[jhq_idx] = -INFINITY;
			}
		}
	}


	/* for each possible read nucleotide other than A */
	for (unsigned char b = 1; b < NUM_NUCLEOTIDES; ++b) {

		/* for each possible read position */
		for (unsigned int j = 0; j < dat->max_read_position; ++j) {

			sat_bjh = 0;
			/* satisfaction associated with read position */
			if (mod->position_power && j) {
				site_power = (double) j / dat->max_read_position;
				sat_bjh += bptr[NUM_NUCLEOTIDES] * site_power;
				for (unsigned int i = 1; i < mod->position_power; ++i) {
					site_power *= (double) j / dat->max_read_position;
					sat_bjh += bptr[NUM_NUCLEOTIDES + i] * site_power;
				}
			}

			/* for each observed quality score */
			for (unsigned char q = 0; q < mod->n_quality; ++q) {

				/* satisfaction associated with quality */
				sat_q = 0;
				if (mod->quality_power && q) {
					qual_power = (double) q / mod->n_quality;
					sat_q += bptr[NUM_NUCLEOTIDES + mod->position_power] * qual_power;
					for (unsigned int i = 1; i < mod->quality_power; ++i) {
						qual_power *= (double) q / mod->n_quality;
						sat_q += bptr[NUM_NUCLEOTIDES + mod->position_power + i] * qual_power;
					}
				}

				/* haplotype effect */
				/* for each possible haplotype source nucleotide */
				for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {

					sat_bjh += bptr[h];	/* source nucleotide */
					jhq_idx =  j * mod->n_hq_combos + h * mod->n_quality + q;
					bjhq_idx = b * mod->n_jhq_combos + jhq_idx;
//if (jhq_idx == 38 && fxn_debug >= DEBUG_II) mmessage(INFO_MSG, NO_ERROR, "sat_bhj = %f, sat_q = %f, p[j=%u,q=%u,h=%c,b=%c]: %f\n", sat_bjh, sat_q, j, q + MIN_ASCII_QUALITY_SCORE, std_to_char[h], std_to_char[b], exp(sat_bjh + sat_q));
					mod->p_bjhq[bjhq_idx] = sat_bjh + sat_q;
//					mod->p_jhq[jhq_idx] += mod->p_bjhq[bjhq_idx];
					if (mod->p_bjhq[bjhq_idx] > mod->max_sat_jhq[jhq_idx])
						mod->max_sat_jhq[jhq_idx] = mod->p_bjhq[bjhq_idx];
					sat_bjh -= bptr[h];
				}

			}
		}

		/* advance to next block: read nucleotide */
		bptr += mod->n_predictors;
	}
#ifdef DEBUGGING_ON
	jhq_idx = 338;
	debug_call(fxn_debug >= DEBUG_II, fxn_debug, mmessage(INFO_MSG, NO_ERROR, "p_bjhq[.][%u][%c][%u]: %.3f %.3f %.3f %.3f\n",
		jhq_idx / mod->n_hq_combos,
		std_to_char[(jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos) / mod->n_quality],
		(jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos - ((jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos) / mod->n_quality) * mod->n_quality) % (dat->max_read_position * mod->n_hq_combos + NUM_NUCLEOTIDES * mod->n_quality) + MIN_ASCII_QUALITY_SCORE,
		mod->p_bjhq[jhq_idx], mod->p_bjhq[jhq_idx + mod->n_jhq_combos], mod->p_bjhq[jhq_idx + 2*mod->n_jhq_combos], mod->p_bjhq[jhq_idx + 3*mod->n_jhq_combos]));
#endif
	for (unsigned int j = 0; j < dat->max_read_position; ++j) {
		jhq_idx = j * mod->n_hq_combos;
		for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
			jhq_idx += h * mod->n_quality;
			for (unsigned char q = 0; q < mod->n_quality; ++q) {
				jhq_idx += q;
				sum = 0;
				for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b) {
					bjhq_idx = b * mod->n_jhq_combos;
					mod->p_bjhq[bjhq_idx + jhq_idx] = exp(
						mod->p_bjhq[bjhq_idx + jhq_idx]
						- mod->max_sat_jhq[jhq_idx]);
					sum += mod->p_bjhq[bjhq_idx + jhq_idx];
				}
				for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b) {
					bjhq_idx = b * mod->n_jhq_combos;
					mod->p_bjhq[bjhq_idx + jhq_idx] /= sum;
				}
				jhq_idx -= q;
			}
			jhq_idx -= h * mod->n_quality;
		}
	}
#ifdef DEBUGGING_ON
	jhq_idx = 338;
	debug_call(fxn_debug >= DEBUG_II, fxn_debug, mmessage(INFO_MSG, NO_ERROR, "p_bjhq[.][%u][%c][%u]: %.3f %.3f %.3f %.3f\n",
		jhq_idx / mod->n_hq_combos,
		std_to_char[(jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos) / mod->n_quality],
		(jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos - ((jhq_idx - (jhq_idx / mod->n_hq_combos) * mod->n_hq_combos) / mod->n_quality) * mod->n_quality) % (dat->max_read_position * mod->n_hq_combos + NUM_NUCLEOTIDES * mod->n_quality) + MIN_ASCII_QUALITY_SCORE,
		mod->p_bjhq[jhq_idx], mod->p_bjhq[jhq_idx + mod->n_jhq_combos], mod->p_bjhq[jhq_idx + 2*mod->n_jhq_combos], mod->p_bjhq[jhq_idx + 3*mod->n_jhq_combos]));
#endif

} /* compute_transition_probabilities */


/**
 * Reallocate model struct for a different K or different sample with different size.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to model_options objet
 * @return	error status
 */
int realloc_model(model *mod, data *dat, model_options *opt)
{
	mod->n_mix = opt->K + opt->background_model;
	mod->K = opt->K;

	/* pi */
	if (!mod->n_mix)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.model.pi");

	double *pi = realloc(mod->pi, mod->n_mix * sizeof *mod->pi);
	double *npi = realloc(mod->npi, mod->n_mix * sizeof *mod->npi);
	double *best_pi = realloc(mod->best_pi,
		mod->n_mix * sizeof *mod->best_pi);

	if (!pi || !npi || !best_pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.model.pi");
	mod->pi = pi;
	mod->npi = npi;
	mod->best_pi = best_pi;

	/* haplotypes */
	if (!dat->max_read_length || !mod->K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.model.haplotypes");

	unsigned char *haplotypes = realloc(mod->haplotypes,
		dat->max_read_length * mod->K * sizeof *mod->haplotypes);
	unsigned char *nhaplotypes = realloc(mod->nhaplotypes,
		dat->max_read_length * mod->K * sizeof *mod->nhaplotypes);
	unsigned char *best_haplotypes = realloc(mod->best_haplotypes,
		dat->max_read_length * mod->K * sizeof *mod->best_haplotypes);

	if (!haplotypes || !nhaplotypes || !best_haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.model.haplotypes");
	mod->haplotypes = haplotypes;
	mod->nhaplotypes = nhaplotypes;
	mod->best_haplotypes = best_haplotypes;

	/* e_ik */
	if (!dat->sample_size || !mod->n_mix)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model.eik");
	double *eik = realloc(mod->eik, dat->sample_size
					* mod->n_mix * sizeof *mod->eik);
	if (!eik)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model.eik");
	mod->eik = eik;

	mod->best_ll = opt->previous_ll;	/* best AECM log likelihood */
	mod->best_init_ll = -INFINITY;	/* best RND-EM criterion */

	/* Set bins for quality score */
	/* keep this in case the sample data has been changed (multi_stage_method) */
	/*
	if (opt->q_model == FIVE_BY_SIX) {
		for (size_t i = 0; i < dat->sample_size; ++i)
			for (size_t j = 0; j < dat->lengths[i]; ++j) {
				unsigned char *dptr = &dat->qmat[i][j];
				if (*dptr < 31)
					*dptr = *dptr / 5;
				else
					*dptr = *dptr - 24;
			}

		mod->n_quality = 6 + (dat->n_quality > 30
			? dat->n_quality - 30 : 0);
	} else {
		mod->n_quality = dat->n_quality;
	}
	*/

	mod->aic = INFINITY;
	mod->bic = INFINITY;

	if(opt->JC69_model)
		mod->distance = realloc(mod->distance, mod->K * sizeof *mod->distance);

	return sync_model(mod, dat, opt);

}/* realloc_model */

/**
 * Overwrite one model with another.
 *
 * @param des_mod	model to overwrite
 * @param ori_mod	source model
 * @param dat		pointer to data object
 * @param partial	??
 * @return		error status
 */
int copy_model(model *des_mod, model *ori_mod, data *dat, int partial)
{
	/* check if K and n_mix are same between two models */
	if (des_mod->K != ori_mod->K || des_mod->n_mix != ori_mod->n_mix)
		return mmessage(ERROR_MSG, INTERNAL_MISMATCH, "copy_model()");
	
	if (des_mod->mopt->parameterization != ori_mod->mopt->parameterization)
		return mmessage(ERROR_MSG, INTERNAL_MISMATCH,
						"mismatched parameterization");
	if (des_mod->mopt->background_model != ori_mod->mopt->background_model)
		return mmessage(ERROR_MSG, INTERNAL_MISMATCH,
						"mismatched parameterization");
	if (des_mod->mopt->JC69_model != ori_mod->mopt->JC69_model)
		return mmessage(ERROR_MSG, INTERNAL_MISMATCH,
						"mismatched parameterization");

	des_mod->iter = ori_mod->iter;
	des_mod->n_param = ori_mod->n_param;

	/* copy best solution */
	des_mod->aic = ori_mod->aic;
	des_mod->bic = ori_mod->bic;
	des_mod->best_ll = ori_mod->best_ll;
	des_mod->best_init_ll = ori_mod->best_init_ll;
	des_mod->JC_ll = ori_mod -> JC_ll;

	memcpy(des_mod->best_pi, ori_mod->best_pi, des_mod->n_mix * sizeof *des_mod->best_pi);
	memcpy(des_mod->best_haplotypes, ori_mod->best_haplotypes, dat->max_read_length * des_mod->K * sizeof *des_mod->best_haplotypes);

	if (des_mod->mopt->parameterization == DEFAULT_PARAMETERIZATION) {
		memcpy(des_mod->best_delta, ori_mod->best_delta, dat->max_read_position * sizeof *des_mod->best_delta);
		memcpy(des_mod->best_gamma, ori_mod->best_gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->best_gamma);
		memcpy(des_mod->best_lambda0, ori_mod->best_lambda0, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->best_lambda0);
		memcpy(des_mod->best_lambda1, ori_mod->best_lambda1, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->best_lambda1);
	} else if (des_mod->mopt->parameterization == QUALITY_PARAMETERIZATION) {
		memcpy(des_mod->best_gamma, ori_mod->best_gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->best_gamma);
	} else if (des_mod->mopt->parameterization == ART_PARAMETERIZATION) {
		memcpy(des_mod->best_gamma, ori_mod->best_gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->best_gamma);
	} else if (des_mod->mopt->parameterization == MLOGIT_PARAMETERIZATION) {
		memcpy(des_mod->best_beta, ori_mod->best_beta, des_mod->n_beta_coef * sizeof *des_mod->best_beta);
	}

	memcpy(des_mod->eik, ori_mod->eik, dat->sample_size * des_mod->n_mix * sizeof *des_mod->eik);

	if (des_mod->mopt->background_model) {
		memcpy(des_mod->best_bg_pi, ori_mod->best_bg_pi, NUM_NUCLEOTIDES * sizeof *des_mod->best_bg_pi);
		if (des_mod->mopt->parameterization == DEFAULT_PARAMETERIZATION)
			memcpy(des_mod->best_bg_lambda, ori_mod->best_bg_lambda, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->best_bg_lambda);
	}
	if (des_mod->mopt->JC69_model) {
		memcpy(des_mod->distance, ori_mod->distance, des_mod->K * sizeof *des_mod->distance);
		memcpy(des_mod->est_ancestor, ori_mod->est_ancestor, dat->max_read_length * sizeof *des_mod->est_ancestor);
	}

	/* Only copy when needed. partial=0 */
	if (!partial) {
		/* copy log likelihood */
		des_mod->ll = ori_mod->ll;
		des_mod->nll = ori_mod->nll;
		des_mod->init_ll = ori_mod->init_ll;

		/* copy parameters in EM */

		memcpy(des_mod->pi, ori_mod->pi, des_mod->n_mix * sizeof *des_mod->pi);
		memcpy(des_mod->npi, ori_mod->npi, des_mod->n_mix * sizeof *des_mod->npi);

		memcpy(des_mod->haplotypes, ori_mod->haplotypes, dat->max_read_length * des_mod->K * sizeof *des_mod->haplotypes);
		memcpy(des_mod->nhaplotypes, ori_mod->nhaplotypes, dat->max_read_length * des_mod->K * sizeof *des_mod->nhaplotypes);

		if (des_mod->mopt->parameterization == DEFAULT_PARAMETERIZATION) {
			memcpy(des_mod->delta, ori_mod->delta, dat->max_read_position * sizeof *des_mod->delta);
			memcpy(des_mod->ndelta, ori_mod->ndelta, dat->max_read_position * sizeof *des_mod->ndelta);

			memcpy(des_mod->gamma, ori_mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->gamma);
			memcpy(des_mod->ngamma, ori_mod->ngamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->ngamma);

			memcpy(des_mod->lambda0, ori_mod->lambda0, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->lambda0);
			memcpy(des_mod->lambda1, ori_mod->lambda1, des_mod->n_quality * sizeof *des_mod->lambda1);
			memcpy(des_mod->nlambda0, ori_mod->nlambda0, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->nlambda0);
			memcpy(des_mod->nlambda1, ori_mod->nlambda1, des_mod->n_quality * sizeof *des_mod->nlambda1);
		} else if (des_mod->mopt->parameterization == QUALITY_PARAMETERIZATION) {
			memcpy(des_mod->gamma, ori_mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->gamma);
			memcpy(des_mod->ngamma, ori_mod->ngamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->ngamma);
		} else if (des_mod->mopt->parameterization == ART_PARAMETERIZATION) {
			memcpy(des_mod->gamma, ori_mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->gamma);
			memcpy(des_mod->ngamma, ori_mod->ngamma, NUM_NUCLEOTIDES_SQUARED * sizeof *des_mod->ngamma);
		} else if (des_mod->mopt->parameterization == MLOGIT_PARAMETERIZATION) {
			memcpy(des_mod->beta, ori_mod->beta, des_mod->n_beta_coef * sizeof *des_mod->beta);
			memcpy(des_mod->nbeta, ori_mod->nbeta, des_mod->n_beta_coef * sizeof *des_mod->nbeta);
		}

		if (des_mod->mopt->background_model) {
			memcpy(des_mod->bg_pi, ori_mod->bg_pi, NUM_NUCLEOTIDES * sizeof *des_mod->bg_pi);
			memcpy(des_mod->nbg_pi, ori_mod->nbg_pi, NUM_NUCLEOTIDES * sizeof *des_mod->nbg_pi);
			if (des_mod->mopt->parameterization == DEFAULT_PARAMETERIZATION) {
				memcpy(des_mod->bg_lambda, ori_mod->bg_lambda, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->bg_lambda);
				memcpy(des_mod->nbg_lambda, ori_mod->nbg_lambda, dat->max_read_position * des_mod->n_quality * sizeof *des_mod->nbg_lambda);
			}
		}
	}

	return NO_ERROR;
} /* copy_model */


/**
 * Read DADA2 error profile from hard-coded filename.  Note, that DADA2
 * reports error rates multiplied by 1000.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param mopt	pointer to model_options object
 * @return	error status
 */
int read_dada2_error_profile(model *mod, data *dat, model_options *mopt)
{
	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	unsigned int nqual;
	double rate;
	const char *error_profile_file = mopt->dada2_file;

	FILE *fp = fopen(error_profile_file, "r");
	if(!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, error_profile_file);
	
	unsigned int cnt = 0;
	while (fscanf(fp, "%lf", &rate) == 1) {
		++cnt;
		fgetc(fp);	/* comma */
	}

	if (cnt % NUM_NUCLEOTIDES_SQUARED)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
							error_profile_file);

	mod->dada2_profile = malloc(cnt * sizeof *mod->dada2_profile);

	if (!mod->dada2_profile)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"model:error_profile");

	nqual = cnt / NUM_NUCLEOTIDES_SQUARED;

	if (dat->n_quality > nqual)
		mmessage(WARNING_MSG, INVALID_USER_INPUT, "DADA2 error profile "
			"allows only %u quality values, but data uses %u\n",
			nqual, dat->n_quality);

	/* update model::n_quality */
	if (mod->n_quality != nqual
		&& (err = realloc_quality_information(mod, dat, mopt, nqual)))
		return err;

	/* min and max are determined by DADA2 error profile */
	mod->min_quality = MIN_ASCII_QUALITY_SCORE;
	mod->max_quality = mod->min_quality + mod->n_quality - 1;
	
	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Number of quality scores: "
		"%u\n", mod->n_quality);

	rewind(fp);

	for (unsigned int r = 0; r < NUM_NUCLEOTIDES_SQUARED; r++)
		for (unsigned int q = 0; q < mod->n_quality ; q++) {
			fscanf(fp, "%lf,", &rate);
			mod->dada2_profile[r * mod->n_quality + q] = rate / 1000;
			debug_msg(!q && DEBUG_I <= fxn_debug, fxn_debug, "r: %i,q:"
				"%i,rate: %f\n", r, q, rate/1000);
		}

	fclose(fp);

	return err;
}/* read_dada2_error_profile */


/**
 * Extract error probability from error profile.
 *
 * [TODO] Consider making this an inline function.
 * @param error_profile	the error profile
 * @param n_quality	the number of possible quality scores
 * @param hap_nuc	true nucleotide
 * @param obser_nuc	observed nucleotide (with error, possibility)
 * @param qual		observed quality score
 * @return		error probability
 */
double dada2_error(double *error_profile, unsigned char n_quality,
	unsigned char hap_nuc, unsigned char obser_nuc, unsigned char qual)
{

	double lp = 0.;

     //A (0),C(1),G(3),T(2)
	if (hap_nuc == XY_A) {
		if (obser_nuc == XY_A)
			lp = error_profile[qual];			//A->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality + qual];		//A->C
		else if (obser_nuc == XY_T)
			lp = error_profile[n_quality * 3 + qual];	//A->T
		else
			lp = error_profile[n_quality * 2 + qual];	//A->G
	} else if (hap_nuc == XY_C) {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality * 4 + qual];	//C->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality * 5 + qual];	//C->C
		else if (obser_nuc == XY_T)
			lp = error_profile[n_quality * 7 + qual];	//C->T
		else
			lp = error_profile[n_quality * 6 + qual];	//C->G
	} else if (hap_nuc == XY_G) {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality * 8 + qual];	//G->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality * 9 + qual];	//G->C
		else if (obser_nuc == XY_T)
			lp = error_profile[n_quality * 11 + qual];	//G->T
		else
			lp = error_profile[n_quality * 10 + qual];	//G->G
	} else {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality * 12 + qual];	//T->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality * 13 + qual];	//T->C
		else if (obser_nuc == XY_T)
			lp = error_profile[n_quality* 15 + qual];	//T->T
		else
			lp = error_profile[n_quality* 14 + qual];	//T->G
	}

	return lp;
}/* dada2_error */

/**
 * Zero expected counts in MLOGIT_PARAMETERIZATION.
 *
 * @param mod	pointer to model struct
 * @param dat	pointer to data struct
 */
void zero_e(model *mod, data *dat)
{
	/* reset expected counts */
	unsigned char *q_idx = mod->qual_list;
	unsigned int idx;
	for (unsigned int j = 0; j < dat->max_read_length; ++j) {
		idx = j * mod->n_hq_combos;
		for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
			idx += h * mod->n_quality;
			for (unsigned int q = 0; q < mod->qual_per_list[j]; ++q) {
				idx += *(q_idx + q);
				mod->e_jhq[idx] = 0;
				for (unsigned char b = 1; b < NUM_NUCLEOTIDES; ++b)
					mod->e_bjhq[b * mod->n_jhq_combos + idx] = 0;
				idx -= *(q_idx + q);
			}
			idx -= h * mod->n_quality;
		}
		q_idx += mod->qual_per_list[j];
	}
} /* zero_e */

/**
 * Delete model object.
 *
 * @param mod	pointer to model object to delete
 */
void free_model(model *mod)
{
	if (!mod)
		return;
	if (mod->pi) free(mod->pi);
	if (mod->npi) free(mod->npi);
	if (mod->best_pi) free(mod->best_pi);
	if (mod->gamma) free(mod->gamma);
	if (mod->ngamma)free(mod->ngamma);
	if (mod->best_gamma) free(mod->best_gamma);
	if (mod->lambda0) free(mod->lambda0);
	if (mod->nlambda0) free(mod->nlambda0);
	if (mod->best_lambda0) free(mod->best_lambda0);
	if (mod->lambda1) free(mod->lambda1);
	if (mod->nlambda1) free(mod->nlambda1);
	if (mod->best_lambda1) free(mod->best_lambda1);
	if (mod->delta) free(mod->delta);
	if (mod->ndelta) free(mod->ndelta);
	if (mod->best_delta) free(mod->best_delta);
	if (mod->qual_list) free(mod->qual_list);
	if (mod->qual_per_list) free(mod->qual_per_list);
	if (mod->haplotypes) free(mod->haplotypes);
	if (mod->nhaplotypes) free(mod->nhaplotypes);
	if (mod->best_haplotypes) free(mod->best_haplotypes);
	if (mod->eik) free(mod->eik);
	if (mod->distance) free(mod->distance);
	if (mod->est_ancestor) free(mod->est_ancestor);
	if (mod->mut_position) free(mod->mut_position);
	if (mod->beta) free(mod->beta);
	if (mod->nbeta) free(mod->nbeta);
	if (mod->best_beta) free(mod->best_beta);
	if (mod->beta_diff) free(mod->beta_diff);
	if (mod->p_bjhq) free(mod->p_bjhq);
//	if (mod->p_jhq) free(mod->p_jhq);
	if (mod->gradient) free(mod->gradient);
	if (mod->hessian) free(mod->hessian);
	if (mod->px) free(mod->px);
	if (mod->e_bjhq) free(mod->e_bjhq);
	if (mod->e_jhq) free(mod->e_jhq);
	free(mod);
	mod = NULL;
} /* free_model */
