/**
 * @file aecm.c
 * @author K. S. Dorman
 *
 * [TODO] currently changing options structure: this file is done when all references to options are variables are via opt->
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef YUDI
#include <cblas.h>
#else
#include <cblas/cblas.h>
#endif

#include "ampliclust.h"
#include "error.h"
#include "aecm.h"
#include "math.h"
#include "simulate.h"
#include "io.h"
#include "model_options.h"
#include "model.h"

void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
void dpotri_(char *UPLO, int *N, double *A, int *LDA, int *INFO);

void m_haplotype(data *dat, model *mod, model_options *opt);
void m_haplotype_mlogit(data *dat, model *mod, model_options *opt);
int m_beta(data *dat, model *mod, model_options *opt);

int aecm(data *dat, model *mod, model_options *opt)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	double pll;
#endif
	double rdelta;
	mod->iter = 0;
	mod->ll = -INFINITY;

	do {
		/* first E step */
		mod->nll = e_step(dat, mod, opt, FIRST);
#ifdef DEBUGGING_ON
		if (mod->nll < mod->ll)
			return mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (I. < %5.3f)****\n", mod->iter, mod->nll, mod->ll);
		pll = mod->nll;
#endif

		/* first CM step for haplotype */
		if (opt->param_estimate & PARAM_HAPLOTYPE) {
#ifdef DEBUGGING_ON
			cc_msg(global_wp, opt->em_info >= SILENT
				|| fxn_debug >= SILENT, fxn_debug, "%5d: "
				"%15.3f (e-step #1)\n", mod->iter, mod->nll);
#endif

			if (opt->parameterization == MLOGIT_PARAMETERIZATION)
				m_haplotype_mlogit(dat, mod, opt);
			else
				m_haplotype(dat, mod, opt);
#ifdef DEBUGGING_ON
			debug_call(fxn_debug >= DEBUG_I, fxn_debug,
				 haplotype_delta_hamming_matrix(mod, dat, opt,
						NEW_COPY, DEFAULT_COPY));
#endif

			/* second E step */
			mod->nll = e_step(dat, mod, opt, SECOND);
#ifdef DEBUGGING_ON
			if (mod->nll < pll)
				return mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (II. < %5.3f)****\n", mod->iter, mod->nll, pll);
#endif
		}

		/* second CM step for all other parameters */
		m_other(dat, mod, opt);

		rdelta = (mod->ll - mod->nll)/mod->ll;

#ifdef DEBUGGING_ON
		if (mod->nll < mod->ll)
			return mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (%5.3e)****\n", mod->iter, mod->nll, rdelta);
#endif

		cc_msg(global_wp, opt->em_info >= SILENT, QUIET,
			"%5d: %15.3f (%5.3e)\n", mod->iter, mod->nll, rdelta);

		/* store latest parameter updates in default slots */
		param_update(mod, dat, opt, DEFAULT_VALUES);

		/* terminate iterations */
		mod->iter++;
		if (mod->iter >= opt->n_iter)
			return MESSAGE(global_wp, WARNING_MSG,
				AECM_EXCEED_MAX_ITERATIONS, "%s: %d\n",
				aecm_error_message(AECM_EXCEED_MAX_ITERATIONS),
				opt->n_iter);

		if (isfinite(mod->ll) && rdelta < opt->epsilon)
			return NO_ERROR;

#ifdef USE_CURSES
		char key = getch();
		if (key == 27)
			return MESSAGE(global_wp, WARNING_MSG, NO_ERROR,
				 "aborting run\n");
#endif
	} while (1);
} /* aecm */


/**
 * Copy parameters from one location to another.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to optinos object
 * @param which	what do we want to copy where
 * 	FROM_BEST: copy model::best_* locations to default locations
 *	BEST_VALUES|OPTIMAL_VALUES: copy default locations to model::best_* locations
 *	DEFAULT_VALUES: copy update locations (model::n*) to default locations
 */
int param_update(model *mod, data *dat, model_options *opt, int which)
{
	if (which == FROM_BEST) {	/* set DEFAULT_VALUES from BEST_VALUES */
		memcpy(mod->haplotypes, mod->best_haplotypes, dat->max_read_length * opt->K * sizeof *mod->haplotypes);
		memcpy(mod->pi, mod->best_pi, mod->n_mix * sizeof *mod->pi);
		if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
			memcpy(mod->delta, mod->best_delta, dat->max_read_position * sizeof *mod->delta);
			memcpy(mod->gamma, mod->best_gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
			memcpy(mod->lambda0, mod->best_lambda0, dat->max_read_position * mod->n_quality * sizeof *mod->lambda0);
			memcpy(mod->lambda1, mod->best_lambda1, dat->max_read_position * mod->n_quality * sizeof *mod->lambda1);
		} else if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
			memcpy(mod->beta, mod->best_beta, mod->n_beta_coef * sizeof *mod->beta);
		} else if ((opt->parameterization == ART_PARAMETERIZATION && opt->art_separate)
			|| opt->parameterization == QUALITY_PARAMETERIZATION) {
			memcpy(mod->gamma, mod->best_gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
		}
		if (opt->background_model) {
			memcpy(mod->bg_pi, mod->best_bg_pi, NUM_NUCLEOTIDES * sizeof *mod->bg_pi);
			memcpy(mod->bg_lambda, mod->best_bg_lambda, dat->max_read_length * mod->n_quality * sizeof *mod->bg_lambda);
		}
	} else if (which) {	/* BEST_VALUES | OPTIMAL_VALUES: same thing for model parameters */
		mod->best_ll = mod->ll;
		memcpy(mod->best_haplotypes, mod->haplotypes, dat->max_read_length * opt->K * sizeof *mod->haplotypes);
		memcpy(mod->best_pi, mod->pi, mod->n_mix * sizeof *mod->pi);
		if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
			memcpy(mod->best_delta, mod->delta, dat->max_read_position * sizeof *mod->pi);
			memcpy(mod->best_gamma, mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
			memcpy(mod->best_lambda0, mod->lambda0, dat->max_read_position * mod->n_quality * sizeof *mod->lambda0);
			memcpy(mod->best_lambda1, mod->lambda1, dat->max_read_position * mod->n_quality * sizeof *mod->lambda1);
		} else if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
			memcpy(mod->best_beta, mod->beta, mod->n_beta_coef * sizeof *mod->beta);
		} else if ((opt->parameterization == ART_PARAMETERIZATION && opt->art_separate)
			|| opt->parameterization == QUALITY_PARAMETERIZATION) {
			memcpy(mod->best_gamma, mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
		}
		if (opt->background_model) {
			memcpy(mod->best_bg_pi, mod->bg_pi, NUM_NUCLEOTIDES * sizeof *mod->bg_pi);
			memcpy(mod->best_bg_lambda, mod->bg_lambda, dat->max_read_length * mod->n_quality * sizeof *mod->bg_lambda);
		}
	} else {		/* DEFAULT_VALUES */
		mod->ll = mod->nll;
		memcpy(mod->haplotypes, mod->nhaplotypes, dat->max_read_length * opt->K * sizeof *mod->haplotypes);
		memcpy(mod->pi, mod->npi, mod->n_mix * sizeof *mod->pi);
		if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
			memcpy(mod->delta, mod->ndelta, dat->max_read_position * sizeof *mod->delta);
			memcpy(mod->gamma, mod->ngamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
			memcpy(mod->lambda0, mod->nlambda0, dat->max_read_position * mod->n_quality * sizeof *mod->lambda0);
			memcpy(mod->lambda1, mod->nlambda1, dat->max_read_position * mod->n_quality * sizeof *mod->lambda1);
		} else if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
			memcpy(mod->beta, mod->nbeta, mod->n_beta_coef * sizeof *mod->beta);
		} else if ((opt->parameterization == ART_PARAMETERIZATION && opt->art_separate)
			|| opt->parameterization == QUALITY_PARAMETERIZATION) {
			memcpy(mod->gamma, mod->ngamma, NUM_NUCLEOTIDES_SQUARED * sizeof *mod->gamma);
		}
		if (opt->background_model) {
			memcpy(mod->bg_pi, mod->nbg_pi, NUM_NUCLEOTIDES * sizeof *mod->bg_pi);
			memcpy(mod->bg_lambda, mod->nbg_lambda, dat->max_read_length * mod->n_quality * sizeof *mod->bg_lambda);
		}
	}
	return NO_ERROR;
} /* param_update */


/**
 * Copy parameters from one model to another.
 *
 * @param to_mod	pointer to model to copy to
 * @param from_mod	pointer to model to copy from
 * @param dat		pointer to data object
 * @param model_options	pointer to model options
 * @return		error status
 */
int param_copy(model *to_mod, model *from_mod, data *dat, model_options *mopt)
{
	memcpy(to_mod->haplotypes, from_mod->haplotypes, dat->max_read_length * mopt->K * sizeof *to_mod->haplotypes);
	memcpy(to_mod->pi, from_mod->pi, to_mod->n_mix * sizeof *to_mod->pi);
	if (mopt->parameterization == DEFAULT_PARAMETERIZATION) {
		memcpy(to_mod->delta, from_mod->delta, dat->max_read_position * sizeof *to_mod->delta);
		memcpy(to_mod->gamma, from_mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *to_mod->gamma);
		memcpy(to_mod->lambda0, from_mod->lambda0, dat->max_read_position * to_mod->n_quality * sizeof *to_mod->lambda0);
		memcpy(to_mod->lambda1, from_mod->lambda1, dat->max_read_position * to_mod->n_quality * sizeof *to_mod->lambda1);
	} else if (mopt->parameterization == MLOGIT_PARAMETERIZATION) {
		memcpy(to_mod->beta, from_mod->beta, to_mod->n_beta_coef * sizeof *to_mod->beta);
	} else if ((mopt->parameterization == ART_PARAMETERIZATION && mopt->art_separate)
		|| mopt->parameterization == QUALITY_PARAMETERIZATION) {
		memcpy(to_mod->gamma, from_mod->gamma, NUM_NUCLEOTIDES_SQUARED * sizeof *to_mod->gamma);
	}
	if (mopt->background_model) {
		memcpy(to_mod->bg_pi, from_mod->bg_pi, NUM_NUCLEOTIDES * sizeof *to_mod->bg_pi);
		memcpy(to_mod->bg_lambda, from_mod->bg_lambda, dat->max_read_length * to_mod->n_quality * sizeof *to_mod->bg_lambda);
	}

	return NO_ERROR;
} /* param_copy */


/**
 * E step.
 *
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param opt		pointer to model_options object
 * @param second	second call after haplotypes updated
 * @return		log likelihood (from previous parameter estimates)
 */
double e_step(data *dat, model *mod, model_options *opt, int second)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//QUIET;//SILENT;//mod->iter > 0;//DEBUG_III;//DEBUG_I;//
#endif
	unsigned int i, j, k;
	unsigned int jhq_idx, bjhq_idx;
	double sum, ll = 0., max, ep = 0.;
	int model_quality = opt->parameterization == DEFAULT_PARAMETERIZATION
		|| opt->parameterization == ART_PARAMETERIZATION
					? opt->model_quality : 0;
	unsigned char *haps = second == SECOND ? mod->nhaplotypes : mod->haplotypes;	/* [KSD, WAS BUG] ??? */

	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
#ifdef DEBUGGING_ON
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "beta: ");
		debug_call(fxn_debug >= DEBUG_I, fxn_debug, fprint_doubles(
			stderr, mod->beta, mod->n_beta_coef, 3, 1));
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "pi: ");
		debug_call(fxn_debug >= DEBUG_I, fxn_debug, fprint_doubles(
			stderr, mod->pi, opt->K, 3, 1));
		compute_transition_probabilities(mod, dat, DEFAULT_VALUES);
#endif
		zero_e(mod, dat);
	}

	for (i = 0; i < dat->sample_size; ++i) {

		unsigned int off = dat->offset ? dat->offset[i] : 0;
		max = -INFINITY;	/* to hold max_k e_{ik} */

		for (k = 0; k < opt->K; ++k) {
			mod->eik[i*mod->n_mix + k] = log(mod->pi[k]);	/* probability of cluster k */
			for (j = 0; j < dat->lengths[i]; ++j) {
				if (opt->parameterization == MLOGIT_PARAMETERIZATION) {

					unsigned char b = dat->dmat[i][j];
					unsigned char h = haps[k*dat->max_read_length + j];
					jhq_idx = (j + off) * mod->n_hq_combos + h * mod->n_quality + dat->qmat[i][j];
					bjhq_idx = b * mod->n_jhq_combos + jhq_idx;
					mod->eik[i*mod->n_mix + k] += log(mod->p_bjhq[bjhq_idx]);						/* prob. h_{kj} -> x_{ij} */
				} else if (opt->parameterization == QUALITY_PARAMETERIZATION) {
					ep = error_prob(dat->fdata, dat->qmat[i][j]);
					if (dat->dmat[i][j] != haps[k*dat->max_read_length + j])
						mod->eik[i*mod->n_mix + k] += log(ep)								/* prob. error [TODO] precompute */
							+ log(mod->gamma[haps[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]]);	/* prob. h_{kj} -> x_{ij} given error */
					else
						mod->eik[i*mod->n_mix + k] += log(1 - ep);							/* prob. error [TODO] precompute  */
				} else if (opt->parameterization == DADA2_PARAMETERIZATION) {
					mod->eik[i*mod->n_mix + k] += log(dada2_error(								/* prob. h_{kj} -> x_{ij} */
						mod->dada2_profile, mod->n_quality,
						haps[k*dat->max_read_length + j],
						dat->dmat[i][j], dat->qmat[i][j]));
				} else if (opt->parameterization == ART_PARAMETERIZATION) {
					unsigned char h = opt->art_separate ? haps[k*dat->max_read_length + j] : 0;
					jhq_idx = (j + off) * mod->n_hq_combos + h * mod->n_quality + dat->qmat[i][j];
					if (model_quality)
						mod->eik[i*mod->n_mix + k] += log(mod->art_profile[jhq_idx]);					/* prob. q_{ij} [TODO] precompute */
					if (!opt->art_separate && dat->dmat[i][j] != haps[k*dat->max_read_length + j])
						mod->eik[i*mod->n_mix + k] += 
							+ log(mod->gamma[haps[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]]);	/* prob. h_{kj} -> x_{ij} given error */
				} else {
					if (dat->dmat[i][j] != haps[k*dat->max_read_length + j])
						mod->eik[i*mod->n_mix + k] += log(1 - mod->delta[j + off])					/* prob. error */
						+ (model_quality ? log(mod->lambda1[(j + off)*mod->n_quality + dat->qmat[i][j]]) : 0)		/* prob. q_{ij} given error */
						+ log(mod->gamma[haps[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]]);		/* prob. h_{kj} -> x_{ij} given error */
					else
						mod->eik[i*mod->n_mix + k] += log(mod->delta[j + off])						/* prob. no error */
							+ (model_quality ? log(mod->lambda0[(j + off)*mod->n_quality + dat->qmat[i][j]]) : 0);	/* prob q_{ij} given no error */
				}
#ifdef DEBUGGING_ON
				if (!isfinite(mod->eik[i*mod->n_mix + k])) {
					mmessage(ERROR_MSG, INTERNAL_ERROR, "eik[%u][%u] = %f at read position %u, iteration %u, round %d (h=%c r=%c q=%u)\n",
						i, k, mod->eik[i*mod->n_mix + k], j + off, mod->iter, second,
						nuc(dat->fdata, (int) haps[k*dat->max_read_length+j]), nuc(dat->fdata, (int) dat->dmat[i][j]), dat->qmat[i][j]);
					exit(0);
				}
#endif

			}
#ifdef DEBUGGING_ON
			debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "e[%zu][%u]: %f\n", i, k, mod->eik[i*mod->n_mix + k]);
#endif
			if (max < mod->eik[i*mod->n_mix + k])
				max = mod->eik[i*mod->n_mix + k];
		}

		if (opt->background_model) {
			/* bg_piA */
			/* bg_piA*bg_lambda */
			mod->eik[i*mod->n_mix + k] = log(mod->pi[k]);
			for (j = 0; j < dat->lengths[i]; ++j)
				mod->eik[i*mod->n_mix + k] += log(mod->bg_pi[(int) dat->dmat[i][j]])
					+ (model_quality 
					? log(mod->bg_lambda[(j + off)*mod->n_quality + dat->qmat[i][j]]) : 0);
			if (max < mod->eik[i*mod->n_mix + k])
				max = mod->eik[i*mod->n_mix + k];
		}

		/* normalize */
		sum = 0.;
		for (k = 0; k < mod->n_mix; ++k) {
			mod->eik[i*mod->n_mix + k] = exp(mod->eik[i*mod->n_mix + k] - max);	/* e_{ik} = e^{-max}\pi_k\prod_{j=1}^J ... */
			sum += mod->eik[i*mod->n_mix + k];					/* e^{-max}\sum_{k=1}^(K+1) \pi_k \prod_{j=1}^J ... */
		}

		/* called to report log likelihood of reads ([TODO] store?) */
		if (second > SECOND) {
			fprintf(active_fp, " %.3f", log(sum) + max);
			continue;
		}

		/* compute previous log likelihood */
		ll += log(sum) + max;							/* log(\sum_{k=1}^(K+1) \pi_k \prod_{j=1}^J ...) - max + max */

#ifdef DEBUGGING_ON
		if (!isfinite(ll)) {
			mmessage(ERROR_MSG, INTERNAL_ERROR, "ll(%u): log(%f) + %f\n", i, sum, max);
			exit(0);
		}
		double tmp = 0;
#endif
		for (k = 0; k < mod->n_mix; ++k) {
			mod->eik[i*mod->n_mix + k] /= sum;				/* e_{ik} = (e^{-max}\pi_k\prod_{j=1}^J ...) / (\sum_{k=1}^(K+1) e^{-max}\pi_k\prod_{j=1}^J...) */

			if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
				double eik = mod->eik[i * mod->n_mix + k];
				for (unsigned int j = 0; j < dat->lengths[i]; ++j) {
					unsigned char b = dat->dmat[i][j];
					unsigned char h = haps[k * dat->max_read_length + j];
					jhq_idx = (j + off) * mod->n_hq_combos + h * mod->n_quality + dat->qmat[i][j];
					bjhq_idx = b * mod->n_jhq_combos + jhq_idx;
					mod->e_bjhq[bjhq_idx] += eik;
					mod->e_jhq[jhq_idx] += eik;
				}
			}

#ifdef DEBUGGING_ON
			tmp += mod->eik[i*mod->n_mix + k];
#endif
		}
#ifdef DEBUGGING_ON
		if (abs(tmp - 1.) > 1e-6) {
			message(stderr, __FILE__, __func__, __LINE__, INFO_MSG, NO_ERROR, "eik[%u] sum != 1.0 (%e)\n", i, tmp = 1.);
			exit(0);
		}
#endif
	}

	return ll;
} /* e_step */


/**
 * M step for multinomial logistic regression.
 *
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param opt	pointer to model_options object
 */
void m_haplotype_mlogit(data *dat, model *mod, model_options *opt)
{
#ifdef DEBUGGING_ON
//	int fxn_debug = ABSOLUTE_SILENCE;//QUIET;//SILENT;//mod->iter > 0;//
#endif
	size_t i;
	unsigned int k, j;
	unsigned int off, bjhq_idx;
	unsigned char l, max = XY_A;
	double sum, dmax;

	for (k = 0; k < opt->K; ++k) {
		for (j = 0; j < dat->max_read_length; ++j) {
			dmax = -INFINITY;			/* maximum likelihood of all data dependent on h_{kj} */
			for (l = 0; l < NUM_NUCLEOTIDES; ++l) {	/* consider all possible haplotype nucleotides */
				sum = 0.;
				for (i = 0; i < dat->sample_size; ++i) {
					off = dat->offset ? dat->offset[i] : 0;

					/* current read i does not cover this position */
					if (j >= dat->lengths[i])
						continue;

					/* transition probabilities are updated as per last E step */

					bjhq_idx = dat->dmat[i][j] * mod->n_jhq_combos + (j + off) * mod->n_hq_combos + l * mod->n_quality + dat->qmat[i][j];
					sum += mod->eik[i*mod->n_mix + k] * log(mod->p_bjhq[bjhq_idx]);
				}
				if (sum > dmax) {
					dmax = sum;
					max = l;
				}
			}
			mod->nhaplotypes[k*dat->max_read_length + j] = max;
		}
	}

} /* m_haplotype_mlogit */


/**
 * M step for original parameterization.
 *
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param opt	pointer to model_options object
 */
void m_haplotype(data *dat, model *mod, model_options *opt)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//QUIET;//SILENT;//mod->iter > 0;//
#endif
	size_t i;
	unsigned int k, j, l, jhq_idx;
	double sum, dmax, ep;
	data_t max = 0;
	int model_quality = opt->parameterization == DEFAULT_PARAMETERIZATION
		|| opt->parameterization == ART_PARAMETERIZATION
						? opt->model_quality : 0;
	unsigned int off;

	for (k = 0; k < opt->K; ++k) {
		for (j = 0; j < dat->max_read_length; ++j) {
			dmax = -INFINITY;	/* maximum likelihood of all data dependent on h_{kj} */
			if (opt->parameterization == DEFAULT_PARAMETERIZATION && mod->delta[j] >= 1.0) {	/* no variation at this site */
				/* assign ancestor from  first descendent covering position since all identical */
				i = 0;
				while (dat->lengths[i] <= j) ++i;
				mod->nhaplotypes[k*dat->max_read_length + j] = dat->dmat[i][j];
				continue;
			}
			for (l = 0; l < NUM_NUCLEOTIDES; ++l) {		/* consider all possible haplotype nucleotides */
				sum = 0.;
				for (i = 0; i < dat->sample_size; ++i) {
					off = dat->offset ? dat->offset[i] : 0;
#ifdef DEBUGGING_ON
					debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "%s:%d: haplotype %u, site %u, read %u", __func__, __LINE__, k, j, i);
#endif

					/* current read i does not cover this position */
					if (j >= dat->lengths[i])
						continue;

					if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
						if (dat->dmat[i][j] == l) {	/* no change */
#ifdef DEBUGGING_ON
							debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " %f * %f * %f", mod->eik[i*mod->n_mix + k], model_quality ? log( mod->lambda0[(j + off)*mod->n_quality + dat->qmat[i][j]]):0, mod->delta[j + off]);
#endif
							sum += mod->eik[i*mod->n_mix + k] * (log(mod->delta[j + off])					/* prob. no error */
								+ (model_quality ? log( mod->lambda0[(j + off)*mod->n_quality + dat->qmat[i][j]]):0));	/* prob q_{ij} given no error */
#ifdef DEBUGGING_ON
							debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " = %f\n", sum);
#endif
						} else {
#ifdef DEBUGGING_ON
							debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " %f * %f * %f * %f", mod->eik[i*mod->n_mix + k], 1 - mod->delta[j + off], mod->gamma[l*NUM_NUCLEOTIDES + dat->dmat[i][j]], model_quality ? log(mod->lambda1[(j + off)*mod->n_quality + dat->qmat[i][j]]) : 0);
#endif
							sum += mod->eik[i*mod->n_mix + k] * (log(1 - mod->delta[j + off])				/* prob. error */
								+ (model_quality ? log(mod->lambda1[(j + off)*mod->n_quality + dat->qmat[i][j]]):0)	/* prob q_{ij} given error */
								+ log(mod->gamma[l*NUM_NUCLEOTIDES + dat->dmat[i][j]]));				/* prob l -> x[i][j] given error */
#ifdef DEBUGGING_ON
								debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " = %f\n", sum);
#endif
						}
					} else if (opt->parameterization == DADA2_PARAMETERIZATION) {
						sum += mod->eik[i*mod->n_mix + k] * dada2_error(					/* prob. h_{kj} -> x_{ij} */
							mod->dada2_profile, mod->n_quality, l,
							dat->dmat[i][j], dat->qmat[i][j]);
					} else if (opt->parameterization == ART_PARAMETERIZATION) {
						jhq_idx = (j + off) * mod->n_hq_combos + (opt->art_separate ? l : 0) * mod->n_quality + dat->qmat[i][j];
						ep = error_prob(dat->fdata, dat->qmat[i][j]);
						if (dat->dmat[i][j] != l)
							sum += mod->eik[i*mod->n_mix + k] * (log(ep)							/* prob. error */
								+ (model_quality ? log(mod->art_profile[jhq_idx]) : 0)					/* prob q_{ij} */
								+ (opt->art_separate ? 0 : log(mod->gamma[l*NUM_NUCLEOTIDES + dat->dmat[i][j]])));	/* prob. h_{kj} -> x_{ij} given error */
						else
							sum += mod->eik[i*mod->n_mix + k] * (log(1 - ep)						/* prob. error */
								+ (model_quality ? log(mod->art_profile[jhq_idx]) : 0));				/* prob q_{ij} */
					} else {	/* QUALITY_PARAMETERIZATION */
						ep = error_prob(dat->fdata, dat->qmat[i][j]);
						if (dat->dmat[i][j] != l)
							sum += mod->eik[i*mod->n_mix + k] * (log(ep)							/* prob. error */
								+ log(mod->gamma[l*NUM_NUCLEOTIDES + dat->dmat[i][j]]));				/* prob. h_{kj} -> x_{ij} given error */
						else
							sum += mod->eik[i*mod->n_mix + k] * log(1 - ep);						/* prob. error */
					}

#ifdef DEBUGGING_ON
					if (!isfinite(sum)) {
						message(stderr, __FILE__, __func__, __LINE__, INFO_MSG, NO_ERROR, "sum = %f (%u %u)\n", sum, dat->dmat[i][j], l);
						exit(0);
					}
#endif
				}
				if (sum > dmax) {
					dmax = sum;
					max = l;
				}
			}
			mod->nhaplotypes[k*dat->max_read_length + j] = max;
#ifdef DEBUGGING_ON
			if (mod->nhaplotypes[k*dat->max_read_length + j] != XY_A
				&& mod->nhaplotypes[k*dat->max_read_length + j] != XY_C
				&& mod->nhaplotypes[k*dat->max_read_length + j] != XY_G
				&& mod->nhaplotypes[k*dat->max_read_length + j] != XY_T) {
				debug_msg(ERROR_MSG, 0, "model::nhapltoypes[%u][%u] is %c (%d)\n", k, j, nuc(dat->fdata, mod->nhaplotypes[k*dat->max_read_length + j]), (int)mod->nhaplotypes[k*dat->max_read_length + j]);
				exit(0);
			}
#endif
		}
	}
} /* m_haplotype */


/**
 * M step for non-haplotype parameters.
 *
 */
void m_other(data *dat, model *mod, model_options *opt)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	double sum2;
#endif
	size_t i;
	unsigned int j, k, l;
	double sum, tmp1, tmp2, tmp3 = 0;
	unsigned int off;
	double epsilon = 1e-8;	/* [TODO, KSD] put in model_options control */
	double *dptr;

	if (opt->param_estimate & PARAM_PI)
		for (k = 0; k < mod->n_mix; ++k)
			mod->npi[k] = epsilon;

	/* use pseudocounts to avoid boundaries */
	if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
		if (opt->param_estimate & PARAM_DELTA)
			for (j = 0; j < dat->max_read_length; ++j)
				mod->ndelta[j] = epsilon;

		if (opt->param_estimate & PARAM_LAMBDA) {
			dptr = mod->nlambda0;
			for (j = 0; j < dat->max_read_length; ++j)
				for (l = 0; l < mod->n_quality; ++l) {
					*dptr = epsilon;
					dptr++;
				}

			dptr = mod->nlambda1;
			for (j = 0; j < dat->max_read_length; ++j)
				for (l = 0; l < mod->n_quality; ++l) {
					*dptr = epsilon;
					dptr++;
				}
		}
	}

	if ((opt->parameterization == DEFAULT_PARAMETERIZATION
		|| opt->parameterization == QUALITY_PARAMETERIZATION
		|| (opt->parameterization == ART_PARAMETERIZATION && opt->art_separate))
				&& opt->param_estimate & PARAM_GAMMA) {
		dptr = mod->ngamma;
		for (i = 0; i < NUM_NUCLEOTIDES; ++i)
			for (j = 0; j < NUM_NUCLEOTIDES; ++j) {
				*dptr = 0.;
				dptr++;
			}
	}
	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
		/* initialize with the last estimate */
		memcpy(mod->nbeta, mod->beta, mod->n_beta_coef * sizeof *mod->beta);
	}

	if (opt->background_model) {
		if (opt->param_estimate & PARAM_LAMBDA) {
			dptr = mod->nbg_lambda;
			for (j = 0; j < dat->max_read_length; ++j)
				for (l = 0; l < mod->n_quality; ++l) {
					*dptr = epsilon;
					dptr++;
				}
		}

		if (opt->param_estimate & PARAM_BG_PI)
			for (i = 0; i < NUM_NUCLEOTIDES; ++i)
				mod->nbg_pi[i] = epsilon;
	}

	tmp1 = 0;
	for (i = 0; i < dat->sample_size; ++i) {
		off = dat->offset ? dat->offset[i] : 0;
		if (opt->param_estimate & PARAM_PI) {
#ifdef DEBUGGING_ON
			sum = 0;
#endif
			for (k = 0; k < mod->n_mix; ++k) {
#ifdef DEBUGGING_ON
				sum += mod->eik[i*mod->n_mix + k];
#endif
				mod->npi[k] += mod->eik[i*mod->n_mix + k];
			}
#ifdef DEBUGGING_ON
		if (fabs(sum - 1.) > 1e-6) {
			message(stderr, __FILE__, __func__, __LINE__, INFO_MSG, NO_ERROR, "eik[%u] sum != 1.0 (%e)\n", i, sum - 1.);
			exit(0);
		}
#endif
		}

		for (j = 0; j < dat->lengths[i]; ++j)
			if (opt->background_model) {
				if (opt->param_estimate & PARAM_LAMBDA)
					mod->nbg_lambda[(j + off)*mod->n_quality + dat->qmat[i][j]] += mod->eik[i*mod->n_mix + opt->K];
				if (opt->param_estimate & PARAM_BG_PI) {
					mod->nbg_pi[(int) dat->dmat[i][j]] += mod->eik[i*mod->n_mix + opt->K];
					tmp1 += mod->eik[i*mod->n_mix + opt->K];	/* [CAUTION] assumes only one background model */
				}
			}
//fprintf(stderr, "paramterization == %d\n", opt->parameterization);
		if (opt->parameterization == MLOGIT_PARAMETERIZATION)
			continue;

		if (opt->parameterization == DEFAULT_PARAMETERIZATION) {
			for (j = 0; j < dat->lengths[i]; ++j)
				for (k = 0; k < opt->K; ++k)
					if (mod->nhaplotypes[k*dat->max_read_length + j] == dat->dmat[i][j]) {
						if (opt->param_estimate & PARAM_DELTA)
							mod->ndelta[j + off] += mod->eik[i*mod->n_mix + k];
						if (opt->param_estimate & PARAM_LAMBDA)
							mod->nlambda0[(j + off)*mod->n_quality + dat->qmat[i][j]] += mod->eik[i*mod->n_mix + k];
					} else {
						if (opt->param_estimate & PARAM_GAMMA)
							mod->ngamma[mod->nhaplotypes[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]] += mod->eik[i*mod->n_mix + k];
						if (opt->param_estimate & PARAM_LAMBDA)
							mod->nlambda1[(j + off)*mod->n_quality + dat->qmat[i][j]] += mod->eik[i*mod->n_mix + k];
					}
		} else if ((opt->parameterization == QUALITY_PARAMETERIZATION
			|| (opt->parameterization == ART_PARAMETERIZATION && opt->art_separate))
			&& opt->param_estimate & PARAM_GAMMA) {
			for (j = 0; j < dat->lengths[i]; ++j)
				for (k = 0; k < opt->K; ++k)
					mod->ngamma[mod->nhaplotypes[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]] += mod->eik[i*mod->n_mix + k];
		}
	}

	if (opt->param_estimate & PARAM_PI)
		for (k = 0; k < mod->n_mix; ++k) {
			mod->npi[k] /= dat->sample_size + epsilon * mod->n_mix;
#ifdef DEBUGGING_ON
			if (mod->npi[k] < 0 || mod->npi[k] > 1) {
				message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "pi[%u]: %f not in [0,1]\n", k, mod->npi[k]);
				exit(EXIT_FAILURE);
			}
#endif
		}

	if (opt->background_model && opt->param_estimate & PARAM_BG_PI)
		for (k = 0; k < NUM_NUCLEOTIDES; ++k)
			mod->nbg_pi[k] /= tmp1 + 4*epsilon;

	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
		if (m_beta(dat, mod, opt))
			exit(EXIT_FAILURE);
#ifdef DEBUGGING_ON
		debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "nbeta: ");
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_doubles(stderr, mod->nbeta, mod->n_beta_coef, 3, 1));
#endif
		return;
	}

/* [BUG] assumes model_options::background_model only ever combined with DEFAULT_PARAMETERIZATION */
	if (opt->parameterization == DEFAULT_PARAMETERIZATION
		&& opt->param_estimate & PARAM_LAMBDA) {
		for (j = 0; j < dat->max_read_length; ++j) {
			/* Using mod->ndelta[j] is numerically unstable because
			 * more error accumulates summing e_{ik} than \lambda_{jl}
			 * terms.  So we recompute the sums here and use them
			 * to normalize. */
			tmp1 = tmp2 = 0;
			for (l = 0; l < mod->n_quality; ++l) {
				tmp1 += mod->nlambda0[j*mod->n_quality + l];
				tmp2 += mod->nlambda1[j*mod->n_quality + l];
			}
#ifdef DEBUGGING_ON
			sum = sum2 = 0;
#endif

			for (l = 0; l < mod->n_quality; ++l) {
				mod->nlambda0[j*mod->n_quality + l] /= tmp1;	//mod->ndelta[j] + mod->n_quality * epsilon;
#ifdef DEBUGGING_ON
				sum2 += mod->nlambda0[j*mod->n_quality + l];
				if (mod->nlambda0[j*mod->n_quality + l] <= 0 || mod->nlambda0[j*mod->n_quality + l] > 1) {
					message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "lambda0[%u][%u]: %f not in [0,1]\n", j, l, mod->nlambda0[j*mod->n_quality + l]);
					exit(EXIT_FAILURE);
				}
#endif
				mod->nlambda1[j*mod->n_quality + l] /= tmp2;	//dat->fdata->n_reads - mod->ndelta[j] + mod->n_quality * epsilon;
#ifdef DEBUGGING_ON
				sum += mod->nlambda1[j*mod->n_quality + l];
				if (mod->nlambda1[j*mod->n_quality + l] <= 0 || mod->nlambda1[j*mod->n_quality + l] > 1) {
//					if (mod->nlambda1[j*mod->n_quality + l] < -epsilon || mod->nlambda1[j*mod->n_quality + l] > 1 + epsilon) {
						message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "lambda1[%u][%u]: %f not in [0,1]\n", j, l, mod->nlambda1[j*mod->n_quality + l]);
						exit(EXIT_FAILURE);
//					}
					if (mod->nlambda1[j*mod->n_quality + l] < 0)
						mod->nlambda1[j*mod->n_quality + l] = 0;
					else if (mod->nlambda1[j*mod->n_quality + l] > 1)
						mod->nlambda1[j*mod->n_quality + l] = 1;
				}
#endif

			}
#ifdef DEBUGGING_ON
			if (fabs(sum - 1.) > 1e-6) {
				message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "lambda1 does not sum to 1 for position %u (%f; delta = %f, n_reads = %u, should be %f)\n", j, sum, mod->ndelta[j], dat->fdata->n_reads, 1.0 / mod->n_quality);
				exit(EXIT_FAILURE);
			}
			if (fabs(sum2 - 1.) > 1e-6) {
				message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "lambda0 does not sum to 1 for position %u (%f; delta = %f, n_reads = %u, should be %f)\n", j, sum, mod->ndelta[j], dat->fdata->n_reads, 1.0 / mod->n_quality);
				exit(EXIT_FAILURE);
			}
#endif
		}
	}
	if (opt->background_model && opt->param_estimate & PARAM_LAMBDA) {
		for (j = 0; j < dat->max_read_length; ++j) {
			tmp3 = 0;
			for (l = 0; l < mod->n_quality; ++l)
				tmp3 += mod->nbg_lambda[j*mod->n_quality + l];
			for (l = 0; l < mod->n_quality; ++l) {
				mod->nbg_lambda[j*mod->n_quality + l] /= tmp3;	//dat->fdata->n_reads - mod->ndelta[j] + mod->n_quality * epsilon;
#ifdef DEBUGGING_ON
				debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "Position %u, quality %u: %f\n", j, l, mod->nbg_lambda[j*mod->n_quality + l]);
				if (mod->nbg_lambda[j*mod->n_quality + l] <= 0 || mod->nbg_lambda[j*mod->n_quality + l] > 1) {
					message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "bg_lambda[%u][%u]: %f not in [0,1]\n", j, l, mod->nbg_lambda[j*mod->n_quality + l]);
					exit(EXIT_FAILURE);
				}
#endif
			}
		}
	}
	if (opt->parameterization == DEFAULT_PARAMETERIZATION
		&& opt->param_estimate & PARAM_DELTA) {
		for (j = 0; j < dat->max_read_length; ++j) {
			mod->ndelta[j] /= dat->coverage[j] + 2*epsilon; //fdata->n_reads; // + 2*epsilon;
			if (mod->ndelta[j] < 0) {
				mmessage(WARNING_MSG, INTERNAL_ERROR, "delta[%u]: %.12e <= 0 (%u)\n", j, mod->ndelta[j], dat->coverage[j]);
				mod->ndelta[j] = 0;
#ifdef DEBUGGING_ON
				if (mod->ndelta[j] > -opt->tolerance) {
					message(stderr, __FILE__, __func__, __LINE__, WARNING_MSG, INTERNAL_ERROR, "delta[%u]: %.12e < 0\n", j, mod->delta[j]);
				} else {
					message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "delta[%u]: %.12f (%.12e) not in [0,1]\n", j, mod->delta[j], mod->delta[j] - 1.);
					exit(EXIT_FAILURE);
				}
#endif
			}
			if (mod->ndelta[j] > 1) {
				mmessage(WARNING_MSG, INTERNAL_ERROR, "delta[%u]: %.12e >= 1 (%u)\n", j, mod->ndelta[j], dat->coverage[j]);
				mod->ndelta[j] = 1.;
#ifdef DEBUGGING_ON
				if (mod->ndelta[j] < 1 + opt->tolerance) {
					message(stderr, __FILE__, __func__, __LINE__, WARNING_MSG, INTERNAL_ERROR, "delta[%u]: %.12f > 1 (delta - 1 = %.12e)\n", j, mod->delta[j], mod->delta[j] - 1.);
				} else {
					message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "delta[%u]: %.12f (%.12e) not in [0,1]\n", j, mod->delta[j], mod->delta[j] - 1.);
					exit(EXIT_FAILURE);
				}
#endif
			}
		}
	}


	if ((opt->parameterization == DEFAULT_PARAMETERIZATION
		|| opt->parameterization == QUALITY_PARAMETERIZATION
		|| (opt->parameterization == ART_PARAMETERIZATION && opt->art_separate))
		&& opt->param_estimate & PARAM_GAMMA) {
		for (i = 0; i < NUM_NUCLEOTIDES; ++i) {
			sum = 0;
			for (j = 0; j < NUM_NUCLEOTIDES; ++j)
				if (i != j)
					sum += mod->ngamma[i*NUM_NUCLEOTIDES + j];
			for (j = 0; j < NUM_NUCLEOTIDES; ++j)
				if (i != j) {
					mod->ngamma[i*NUM_NUCLEOTIDES + j] /= sum;
#ifdef DEBUGGING_ON
					if (mod->ngamma[i*NUM_NUCLEOTIDES + j] <= 0 || mod->ngamma[i*NUM_NUCLEOTIDES + j] >= 1) {
						message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "gamma[%u][%u]: %f not in (0,1)\n", i, j, mod->ngamma[i*NUM_NUCLEOTIDES + j]);
						exit(EXIT_FAILURE);
					}
#endif
				}
		}
	}

} /* m_other */


/**
 * M step for multinomial logistic regression coefficients.
 *
 * For models that are not overparameterized, the multinomial logistic
 * regression likelihood is concave everywhere and thus can be
 * numerically obtained by Newton-Raphson.
 *
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param opt	pointer to model_options object
 * @return	error status
 */
int m_beta(data *dat, model *mod, model_options *opt)
{
#ifdef DEBUGGING_ON
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int condn = DEBUG_I <= fxn_debug, lcondn;
#endif
	int err = NO_ERROR;
	unsigned int grad_idx = 0;	/* index of gradient */
	unsigned int bjhq_idx, jhq_idx;	/* indices of transitions */
	unsigned int q_idx = 0;		/* quality list index */
	unsigned int idx, px_idx;
	double site_power, site_power2;	/* x's for site */
	double qual_power, qual_power2;	/* x's for quality */
	int nbeta = mod->n_beta_coef;	/* first intercept is 0 */
	unsigned int nbeta_sq = nbeta * nbeta;
	double diff;			/* y_{ij} - P_{ij} [mlogit] */
	char uplo = 'U';		/* lapack: use upper triangle */
	int info;			/* lapack: error status */
	unsigned int max_iter = 10000;	/* hard-coded max. iterations */
	unsigned int niter = 0;		/* no. of iterations */
	double enorm_diff, enorm_vec;	/* Euclidean distances */

	/* The transition probabilities model.p_bjhq have been computed from
 	 * the previous estimates of beta, while model.e_bjhq and model.e_jhq have
	 * been computed during the last E step.
	 */
#ifdef DEBUGGING_ON
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "beta (initial): ");
	debug_call(fxn_debug >= DEBUG_I, fxn_debug, fprint_doubles(stderr,
					mod->beta, mod->n_beta_coef, 3, 1));
#endif

	/* loop until convergence */
	do {
		/* reset */
		q_idx = 0;
		if (niter)
			memcpy(mod->beta, mod->nbeta, mod->n_beta_coef * sizeof *mod->beta);

		compute_transition_probabilities(mod, dat, DEFAULT_VALUES);

		/* reset gradient and hessian */
		for (int l = 0; l < nbeta; ++l)
			mod->gradient[l] = 0;
		for (unsigned int l = 0; l < nbeta_sq; ++l)
			mod->hessian[l] = 0;
#ifdef DEBUGGING_ON
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_doubles(stderr, mod->p_bjhq, mod->n_satisfaction_indices, 3, 1));
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_vectorized_sq_matrix(stderr, mod->hessian, nbeta, 0));
#endif

	/* gradient and hessian are sums over j, h, q, and b */
	for (unsigned int j = 0; j < dat->max_read_length; ++j) {
		for (unsigned char h = 0; h < NUM_NUCLEOTIDES; ++h) {
			for (unsigned char q = 0; q < mod->qual_per_list[j]; ++q) {

				/* px is a sum over b */
				for (unsigned char l = 0; l < nbeta; ++l)
					mod->px[l] = 0;

				jhq_idx = j * mod->n_hq_combos + h * mod->n_quality + mod->qual_list[q_idx + q];

				/* compute contribution of combo {j, h, q, b} to gradient and and vector \sum_b gamma_{jhqp}x_{jhqp} */
				for (unsigned char b = 1; b < NUM_NUCLEOTIDES; ++b) {

					grad_idx = mod->n_predictors * (b - 1);
					bjhq_idx = b * mod->n_jhq_combos + jhq_idx;

					diff = mod->e_bjhq[bjhq_idx] - mod->e_jhq[jhq_idx] * mod->p_bjhq[bjhq_idx];
#ifdef DEBUGGING_ON					
					lcondn = condn && 0;
					debug_msg(lcondn, fxn_debug, "j=%u, h=%u, q=%u, b=%u: diff = %f (%f[%u] - %f[%u] * %f)\n", j, h, mod->qual_list[q_idx + q], b, diff, mod->e_bjhq[bjhq_idx], bjhq_idx, mod->e_jhq[jhq_idx], jhq_idx, mod->p_bjhq[bjhq_idx]);
#endif
					/* intercept and haplotype effect: x=1 */
					mod->gradient[grad_idx + h] += diff;
					mod->px[grad_idx + h] += mod->p_bjhq[bjhq_idx];

					/* site effect: j, j^2, j^3, ... */
					site_power = (double) j / dat->max_read_length;
					for (unsigned int l = 0; l < mod->position_power; ++l) {
						mod->gradient[grad_idx + NUM_NUCLEOTIDES + l] += diff * site_power;
						mod->px[grad_idx + NUM_NUCLEOTIDES + l] += mod->p_bjhq[bjhq_idx] * site_power;
						site_power *= (double) j / dat->max_read_length;
					}

					/* quality effect: q, q^2, q^3, ... */
					qual_power = (double) mod->qual_list[q_idx + q] / dat->n_quality;
					for (unsigned int l = 0; l < mod->quality_power; ++l) {
						mod->gradient[grad_idx + NUM_NUCLEOTIDES + mod->position_power + l] += diff * qual_power;
						mod->px[grad_idx + NUM_NUCLEOTIDES + mod->position_power + l] += mod->p_bjhq[bjhq_idx] * qual_power;
						qual_power *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
					}
				}

#ifdef DEBUGGING_ON
				debug_msg(j == 5 && DEBUG_II <= fxn_debug, fxn_debug, "j=%u, h=%u, q=%u: ", j, h, mod->qual_list[q_idx + q]);
				debug_call(j == 5 && DEBUG_II <= fxn_debug, fxn_debug, fprint_doubles(stderr, mod->px, nbeta, 3, 1));
#endif

				/* only non-zero contributions to b-block */

				for (unsigned char b = 1; b < NUM_NUCLEOTIDES; ++b) {
					//lcondn = condn && 0;

					grad_idx = mod->n_predictors * (b - 1);
					bjhq_idx = b * mod->n_jhq_combos + jhq_idx;
					double ep = mod->e_jhq[jhq_idx] * mod->p_bjhq[bjhq_idx];
#ifdef DEBUGGING_ON
					debug_call(lcondn, fxn_debug, fprintf(stderr, "e[j=%u][q=%u][h=%u][b=%u] * p[j=%u][q=%u][h=%u][b=%u] (%.3f * %.3f = %.3f): ", j, mod->qual_list[q_idx + q], h, b, j, mod->qual_list[q_idx + q], h, b, mod->e_jhq[jhq_idx], mod->p_bjhq[bjhq_idx], ep));
#endif
					if (ep == 0)
						continue;

#ifdef DEBUGGING_ON
					if (lcondn) {
						fprintf(stderr, "x:");
						for (unsigned char r = 1; r < NUM_NUCLEOTIDES; ++r) {
							if (r == b) {
								for (unsigned char r2 = 0; r2 < NUM_NUCLEOTIDES; ++r2)
									fprintf(stderr, " %.3f", (double)(r2 == h));
								site_power = (double) j / dat->max_read_length;
								for (unsigned int l = 0; l < mod->position_power; ++l) {
									fprintf(stderr, " %.3f", site_power);
									site_power *= (double) j / dat->max_read_length;
								}
								qual_power = (double) mod->qual_list[q_idx + q] / dat->n_quality;
								for (unsigned int l = 0; l < mod->quality_power; ++l) {
									fprintf(stderr, " %.3f", qual_power);
									qual_power *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
								}
							} else {
								for (unsigned int l = 0; l < mod->n_predictors; ++l)
								fprintf(stderr, " 0.000");
							}
						}
						fprintf(stderr, "\n");
					}
#endif

					/* First four rows, intercept and haplotype
					 * effects.
					 */
					for (unsigned char m = 0; m < NUM_NUCLEOTIDES; ++m) {
						/* [grad_idx + m][grad_idx + h] */
						idx = grad_idx + m + nbeta * (grad_idx + h);
						mod->hessian[idx] -= - ep * ((m == h) - mod->px[grad_idx + m]);
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(m=%u, h=%u) idx=%u %f (ep=%f)\n", m, h, idx, - ep * ((m == h) - mod->px[m]), mod->hessian[idx]);
#endif

						site_power = (double) j / dat->max_read_length;
						for (unsigned int l1 = 0; l1 < mod->position_power; ++l1) {
							/* [grad_idx + m][grad_idx + NUM_NUCLEOTIDES + l1] */
							idx = grad_idx + m + nbeta * (grad_idx + NUM_NUCLEOTIDES + l1);
							mod->hessian[idx] -= - ep * ((m == h) - mod->px[grad_idx + m]) * site_power;
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(m=%u, l1=%u) idx=%u %f (ep=%f)\n", m, NUM_NUCLEOTIDES + l1, idx, - ep * ((m == h) - mod->px[grad_idx + m]) * site_power, mod->hessian[idx]);
#endif
							site_power *= (double) j / dat->max_read_length;
						}
						qual_power = (double) mod->qual_list[q_idx + q] / dat->n_quality;
						for (unsigned int l2 = 0; l2 < mod->quality_power; ++l2) {
							/* [grad_idx + m][grad_idx + NUM_NUCLEOTIDES + mod->position_power + l2] */
							idx = grad_idx + m + nbeta * (grad_idx + NUM_NUCLEOTIDES + mod->position_power + l2);
							mod->hessian[idx] -= - ep * ((m == h) - mod->px[grad_idx + m]) * qual_power;
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(m=%u, l2=%u) idx=%u %f (ep=%f)\n", m, NUM_NUCLEOTIDES + l2, idx, - ep * ((m == h) - mod->px[grad_idx + m]) * qual_power, mod->hessian[idx]);
#endif
							qual_power *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
						}
					}

					/* rows for site effect */
					site_power = (double) j / dat->max_read_position;
					px_idx = NUM_NUCLEOTIDES;
					for (unsigned int l = 0; l < mod->position_power; ++l) {
						site_power2 = (double) j / dat->max_read_position;
						for (unsigned int l1 = 0; l1 < mod->position_power; ++l1) {
							/* [grad_idx + NUM_NUCLEOTIDES + l][grad_idx + NUM_NUCLEOTIDES + l1] */
							idx = grad_idx + NUM_NUCLEOTIDES + l + nbeta * (grad_idx + NUM_NUCLEOTIDES + l1);
							mod->hessian[idx] -= - ep * (site_power - mod->px[grad_idx + px_idx + l]) * site_power2;
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(l=%u, l1=%u) idx=%u %f (ep=%f)\n", NUM_NUCLEOTIDES + l, NUM_NUCLEOTIDES + l1, idx, - ep * (site_power - mod->px[grad_idx + px_idx + l]) * site_power2, mod->hessian[idx]);
#endif
							site_power2 *= (double) j / dat->max_read_position;
						}
						qual_power = (double) mod->qual_list[q_idx + q] / dat->n_quality;
						for (unsigned int l2 = 0; l2 < mod->quality_power; ++l2) {
							/* [grad_idx + NUM_NUCLEOTIDES + l][grad_idx + NUM_NUCLEOTIDES + mod->position_power + l2] */
							idx = grad_idx + NUM_NUCLEOTIDES + l + nbeta * (grad_idx + NUM_NUCLEOTIDES + mod->position_power + l2);
							mod->hessian[idx] -= - ep * (site_power - mod->px[grad_idx + px_idx + l]) * qual_power;
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(l=%u, l2=%u) idx=%u %f (ep=%f)\n", NUM_NUCLEOTIDES + l, NUM_NUCLEOTIDES + mod->position_power + l2, idx, - ep * (site_power - mod->px[grad_idx + px_idx + l]) * qual_power, mod->hessian[idx]);
#endif
							qual_power *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
						}
						site_power *= (double) j / dat->max_read_position;
					}

					/* rows for quality effect */
					qual_power = (double) mod->qual_list[q_idx + q] / dat->n_quality;
					px_idx = NUM_NUCLEOTIDES + mod->position_power;
					for (unsigned int l = 0; l < mod->quality_power; ++l) {
						qual_power2 = (double) mod->qual_list[q_idx + q] / dat->n_quality;
						for (unsigned int l2 = 0; l2 < mod->quality_power; ++l2) {
							/* [grad_idx + NUM_NUCLEOTIDES + mod->position_power + l][grad_idx + NUM_NUCLEOTIDES + mod->position_power + l2] */
							idx = grad_idx + px_idx + l + nbeta * (grad_idx + px_idx + l2);
							mod->hessian[idx] -= - ep * (qual_power - mod->px[grad_idx + px_idx + l]) * qual_power2;
#ifdef DEBUGGING_ON
if (lcondn) fprintf(stderr, "(l=%u, l2=%u) idx=%u %f (ep=%f)\n", l + px_idx, l2 + px_idx, idx, - ep * (qual_power - mod->px[grad_idx + px_idx + l]) * qual_power2, mod->hessian[idx]);
#endif
							qual_power2 *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
						}
						qual_power *= (double) mod->qual_list[q_idx + q] / dat->n_quality;
					}
				}

			}
		}
		q_idx += mod->qual_per_list[j];
	}
#ifdef DEBUGGING_ON
		debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "Gradient: ");
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_doubles(stderr, mod->gradient, nbeta, 3, 1));
		if (fxn_debug >= DEBUG_II && 1)
			fprintf(stderr, "Hessian:\n");
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_vectorized_sq_matrix(stderr, mod->hessian, nbeta, 0));
#endif

		/* invert hessian: it is symmetric and positive definite */
		/* Cholesky factorization */
		dpotrf_(&uplo, &nbeta, mod->hessian, &nbeta, &info);
		if (info)
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "dpotrf_ "
							"error: %d\n", info);
#ifdef DEBUGGING_ON
		if (condn && 0)
			fprintf(stderr, "LU:\n");
		debug_call(condn && 0, fxn_debug, fprint_vectorized_sq_matrix(stderr, mod->hessian, nbeta, 0));
#endif
		/* matrix inversion */
		dpotri_(&uplo, &nbeta, mod->hessian, &nbeta, &info);
		if (info)
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "dpotri_ "
							"error: %d\n", info);

		for (unsigned int l = 0; l < (unsigned int)nbeta; ++l)
			for (unsigned int j = l; j < (unsigned int)nbeta; ++j)
				mod->hessian[l*nbeta + j]
						= mod->hessian[j*nbeta + l];
#ifdef DEBUGGING_ON
		if (condn && 0)
			fprintf(stderr, "Inverse:\n");
		debug_call(condn && 0, fxn_debug, fprint_vectorized_sq_matrix(stderr, mod->hessian, nbeta, 0));
#endif
		/* symmetric matrix times vector */
		cblas_dsymv(CblasColMajor, CblasUpper, nbeta, 1, mod->hessian,
				nbeta, mod->gradient, 1, 0, mod->nbeta, 1);

		/* Newton-Raphson update */
		for (unsigned int i = 0; i < mod->n_beta_coef; ++i) {
			mod->nbeta[i] = mod->beta[i] + mod->nbeta[i];	/* reverse sign */
			mod->beta_diff[i] = mod->nbeta[i] - mod->beta[i];
		}

		enorm_diff = cblas_dnrm2(mod->n_beta_coef, mod->beta_diff, 1);
		enorm_vec = cblas_dnrm2(mod->n_beta_coef, mod->beta, 1);
		if (!enorm_vec)
			enorm_vec = 1.;
		++niter;
#ifdef DEBUGGING_ON
		debug_call(condn, fxn_debug, mmessage(DEBUG_MSG, NO_ERROR, " beta: "));
		debug_call(condn, fxn_debug, fprint_doubles(stderr, mod->beta, mod->n_beta_coef, 3, 1));
		debug_call(condn, fxn_debug, mmessage(DEBUG_MSG, NO_ERROR, "nbeta: "));
		debug_call(condn, fxn_debug, fprint_doubles(stderr, mod->nbeta, mod->n_beta_coef, 3, 1));
		debug_call(condn, fxn_debug, mmessage(DEBUG_MSG, NO_ERROR, "diff (%2u): %g / %g (%g)\n", niter, enorm_diff, enorm_vec, sqrt(enorm_diff / enorm_vec)));
#endif

	} while (sqrt(enorm_diff / enorm_vec) >= opt->beta_epsilon && niter < max_iter);

	if (niter >= max_iter) {
		mmessage(ERROR_MSG, EXCEED_ITERATIONS, "%u\n", max_iter);
		err = EXCEED_ITERATIONS;
	}

	return err;
} /* m_beta */


/**
 * Return description of AECM error.
 *
 * @param err_no	error number
 * @return		string describing error
 */
const char *aecm_error_message(int err_no) {
	if (err_no == AECM_EXCEED_MAX_ITERATIONS)
		return "Exceed maximum allowed iterations";
	else if (err_no == AECM_ASCENT_VIOLATION)
		return "Log likelihood decline";
	else
		return "";
} /* aecm_error_message */


/**
 * Display Hamming distance matrix between sets of haplotypes.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to model_options object	[TODO]
 * @param rows	which copy of haplotypes for rows
 * @param cols	which copy of haplotypes for columns
 */
void haplotype_delta_hamming_matrix(model *mod, data *dat, model_options *opt,
	int rows, int cols)
{
	unsigned char *rhaps = rows == NEW_COPY ? mod->nhaplotypes : mod->haplotypes;
	unsigned char *chaps = cols == NEW_COPY ? mod->nhaplotypes : mod->haplotypes;

	//fprint_fasta(stderr, mod->nhaplotypes, opt->K, dat->max_read_length, "E");
	for (unsigned int j = 0; j < opt->K; ++j) {
/* Nice idea, but a function in simulate.c should to do it.
		if (opt->do_simulation)
			for (unsigned int k = 0; k < opt->true_K; ++k)
				fprintf(stderr, " %2lu", hamming_char_dis(
					(char *) &mod->nhaplotypes[j*dat->max_read_length],
					(char *) &mod->true_haplotypes[k*dat->max_read_length],
					dat->max_read_length));
*/
		for (unsigned int k = 0; k < opt->K; ++k)
			fprintf(stderr, " %2lu", hamming_char_dis(
				(char *) &rhaps[j*dat->max_read_length],
				(char *) &chaps[k*dat->max_read_length],
				dat->max_read_length));
		fprintf(stderr, "\n");
	}
} /* haplotype_delta_hamming_matrix */
