/**
 * @file model.h
 * @author Karin S. Dorman
 *
 * Header file for model struct.
 */

#ifndef __H_MODEL__
#define __H_MODEL__

#include "data.h"
#include "model_options.h"

typedef struct _model model;

extern const unsigned int NUM_NUCLEOTIDES_SQUARED;

/**
 * Parameter sets.
 */
enum {
	PARAM_HAPLOTYPE = 1,
	PARAM_DELTA = 2,
	PARAM_LAMBDA = 4,
	PARAM_GAMMA = 8,
	PARAM_PI = 16,
	PARAM_BG_PI = 32,
	PARAM_BETA = 64,
};

/**
 * Model for bins
 * */
enum {
	NO_BINS = 0,
	BINS_LAMBDA0 = 1,
  	BINS_LAMBDA1 = 2,
	BINS_QUALITY = 4,
};

/**
 * Model parameterizations for nucleotides and quality scores.
 * Notice that MLOGIT, QUALITY, and DADA2 need a quality score
 * generator for simulation, which can be lambda or ART.
 */
enum {
	DEFAULT_PARAMETERIZATION,	/* H, pi, delta, lambda0|1, gamma */
	MLOGIT_PARAMETERIZATION,	/* H, pi, beta[, lambda] (q covariates) */
	QUALITY_PARAMETERIZATION,	/* H, pi, gamma[, lambda] (q literal) */
	DADA2_PARAMETERIZATION,		/* H, pi[, lambda] (q covariates) */
	ART_PARAMETERIZATION,		/* H, pi, gamma, lambdaA, lambdaC, ... */
	NUM_PARAMETERIZATIONS
};

enum{
	EQUAL_LENGTH = 2,
	EXPECTATION_MODEL = 6,  /* for lambda0 only */
};

 

struct _model {
	double ll;		/*<! current log likelihood */
	double nll;		/*<! next log likelihood */
	double best_ll;		/*<! best log likelihood so far */
	double init_ll;		/*<! best initialization log likelihood [KSD: move to initialization.h] */
	double best_init_ll;	/*<! best initialization log like. across initializations [KSD: move to initialization.h] */
	double JC_ll;      	/*<! log likelihood of haplotypes of JC69 model */
	model_options *mopt;	/*<! pointer to model_options object [KLUDGE] */

	/* status */
	size_t iter;		/*<! no. of EM iterations */
	unsigned synced;	/*<! has this model been synced? */

	unsigned char n_quality;/*<! no. of quality scores optionally
				 * post-compression; must be >=
				 * \ref data:n_quality.
				 */
	unsigned char min_quality;
	unsigned char max_quality;

	double *eik;		/*<! expectation: P(C_i=k|r_i,q_i) */
	double *e_bjhq;		/*<! expectation: E[h_j->b_j|q_j) */
	double *e_jhq;		/*<! expectation: E[h_j|q_j) */

//	int parameterization;		/*<! see Model parameterizations */
	unsigned int n_param;		/*<! number of parameters */
	unsigned int n_param_lambda;	/*<! no. of parameters in lambda0|1 */
	unsigned int n_mix;		/*<! number of mixing components */
	unsigned int K; 		/*<! number of clusters */

	/* mlogit model */
	unsigned int n_predictors;	/*<! number of predictors */
	unsigned int n_beta_coef;	/*<! number of beta coefficients */
	unsigned int n_jhq_combos;	/*<! site, haplotype, quality */
	unsigned int n_hq_combos;	/*<! haplotype, quality */
	unsigned int n_jhqp_combos;	/*<! site, haplotype, quality, beta */
	unsigned int position_power;	/*<! max. power on the read position */
	unsigned int quality_power;	/*<! max. power on quality score */
	/* error x quality score [TODO] */
	unsigned int eq_inter_q_power;	/*<! max. power on quality */
	unsigned int ep_inter_p_power;	/*<! max. power on quality */
	/* position x quality score */
	unsigned int pq_inter_p_power;	/*<! max. power on position */
	unsigned int pq_inter_q_power;	/*<! max. power on quality */
	/* haplotype x position and quality */
	unsigned int hp_inter_p_power;	/*<! max. power on position */
	unsigned int hq_inter_q_power;	/*<! max. power on quality */
	double *gradient;		/*<! gradient for mlogit NR */
	double *hessian;		/*<! hessian for mlogit NR */
	double *px;			/*<! used in mlogit NR */

	/* parameters */
	double *pi;		/*<! kx1 mixing proportions */
	double *npi;		/*<! updated kx1 mixing proportions */
	double *delta;		/*<! lx1 prob. error per position */
	double *ndelta;		/*<! updated lx1 prob. error per position */
	double *gamma;		/*<! 4x4 substitution probabilities */
	double *ngamma;		/*<! updated 4x4 substitution probabilities */
	double *lambda0;	/*<! qxl quality score pmfs */
	double *nlambda0;	/*<! updated qxl quality score pmfs */
	double *lambda1;	/*<! qx1 quality score pmfs */
	double *nlambda1;	/*<! updated qx1 quality score pmfs */
	double *lambda;		/*<! marginal pmf for q */
	double *beta;		/*<! coefficients in mlogit model */
	double *nbeta;		/*<! update coefficients in mlogit model */

	/* dada2 parameterization */
	double *dada2_profile;		/*<! DADA2-style error profile */

	/* ART parameterization */
	double *art_profile;	/*<! store art profile */

	unsigned char *haplotypes;	/*<! haplotypes */
	unsigned char *nhaplotypes;	/*<! updated haplotypes */


	double *bg_pi;		/*<! background nucleotide proportions */
	double *nbg_pi;		/*<! next background nucleotide proportions */
	double *bg_lambda;	/*<! background quality score distn */
	double *nbg_lambda;	/*<! next background quality score distn */

	/* best solution */
	unsigned char *best_haplotypes;
	double *best_pi;
	double *best_delta;
	double *best_gamma;
	double *best_lambda0;
	double *best_lambda1;
	double *best_beta;
	double *best_bg_pi;
	double *best_bg_lambda;

	/* simulation control */
	unsigned int *mut_position;	/*<! l x 1 positions of mut in hap */

	/* for model comparison */
	double aic;
	double bic;

	/* for the JC69 model */
	double * distance; /*<! Kx1 distances from haplotypes to the ancestor */
	unsigned char * est_ancestor; /*<! The estimated ancestor */

	/* mlogit model: efficient representation
	 * Reads will always be independent, sites are currently independent.
	 * Memory is arranged as reads (not by sites), but sites within a read
	 * are predictable and quality scores are not.  While loading a read
	 * into memory and working on its sites is efficient, loading the
	 * corresponding transition probabilities (all possible sites, any
	 * possible quality score, any possible haplotype nucleotide) is not.
	 * To avoid swapping cache lines, the solution is to rearrange the
	 * data by appropriate covariate, for example by site, haplotype
	 * nucleotide, and quality score.  If we wanted to model dependence
	 * between sites, the solution is to rearrange on site, context,
	 * haplotype nucleotide, then quality score.
	 *
	 * Also, L1 cache size is 32Kb, which fits 4000 ints or doubles.  There
	 * are slow downs at this limit and at 256Kb L2 cache also.
	 */
	unsigned int n_satisfaction_indices;	/*<! no. V_{ij}, mlogit.pdf */
	double *p_bjhq;		/*<! transition probability h->b given j,q */
	double *p_jhq;		/*<! sum transition counts h->b given j,q */
	double *max_sat_jhq;		/*<! internal to avoid overflow */
	unsigned int *qual_per_list;	/*<! no. quality scores observed / site */
	unsigned char *qual_list;	/*<! list of observed quality scores */
	double *beta_diff;	/*<! change in beta coefficients */

}; /* model */

int make_model(model **mod, data *dat, model_options *opt, int delay);
int sync_model(model *mod, data *dat, model_options *opt);
int realloc_quality_information(model *mod, data *dat, model_options *opt, unsigned int nquality);
int modified_ic(unsigned char *hap, unsigned char *est_anc, double *distance, double best_ll, unsigned int K, double *JC_ll, double *n_aic, double *n_bic, unsigned int n_param, unsigned int max_read_length, size_t sample_size);
int realloc_model(model *mod, data *dat, model_options *opt);
int copy_model(model *des_mod, model *ori_mod,data *dat, int partial);
void compute_transition_probabilities(model *mod, data *dat, int type);
void zero_e(model *mod, data *dat);
int read_dada2_error_profile(model *mod, data *dat, model_options *opt);
double dada2_error(double *error_profile, unsigned char n_quality,unsigned char hap_nuc,unsigned char obser_nuc, unsigned char qual);
void free_model(model *mod);

#endif
