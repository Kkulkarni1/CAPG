/**
 * @file model_options.h
 * @author K. S. Dorman
 */

#ifndef __MODEL_OPTIONS_H__
#define __MODEL_OPTIONS_H__

typedef struct _model_options model_options;

struct _model_options {
	/* model */
	int parameterization;	/*<! model parameterization */
	int q_model;		/*<! various q-models */
	int model_quality;	/*<! ignore quality during estimation */
	int background_model;	/*<! add a background category: iid nucs */
	int JC69_model;		/*<! use approx. JC69 hierarchical model */

	char const *dada2_file;	/*<! dada2 file for error profile */

	char const *art_file;	/*<! art file for error profiles */
	int art_separate;	/*<! separate error profiles */

	/* number of clusters */
	unsigned int K;		/*<! number of clusters */
	unsigned int max_K;	/*<! max number of clusters considered */
	unsigned int min_K;	/*<! min number of clusters considered */

	/* estimation */
	int use_aic;		/*<! use aic to do model selection */
	double previous_ll;	/*<! best log like.: save only if better */
	int em_info;		/*<! level of information to output */
	char param_estimate;	/*<! choose params. to estimate (model.h) */

	/* estimation run info */
	double epsilon;		/*<! change in relative log likelihood */
	double tolerance;	/*<! small number allowed numerical error */
	double beta_epsilon;	/*<! relative change in beta for mlogit */
	unsigned int n_iter;	/*<! max. number of AECM iterations */
}; /* model_options */

int make_model_options(model_options **);
int process_parameterization_option(const char *str, const int argc, const char **argv, int i, model_options *modo);
int process_estimation_option(const char **argv, int i, model_options *modo);


#endif
