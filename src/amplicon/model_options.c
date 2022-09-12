/**
 * @file model_options.c
 * @author Karin S. Dorman
 *
 * Command-line options handling.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>

#include "model_options.h"
#include "model.h"
#include "error.h"

int make_model_options(model_options **modo_in)
{
	*modo_in = malloc(sizeof **modo_in);

	if (!*modo_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model_options");

	model_options *modo = *modo_in;

	modo->K = 0;	/* invalid option */

	modo->parameterization = DEFAULT_PARAMETERIZATION;
	modo->q_model = DEFAULT_Q_MODEL;
	modo->model_quality = 1;
	modo->background_model = 0;
	modo->JC69_model = 1;
	modo->dada2_file = "Illumina_profiles/DADA2_mock_error_x1000_file.csv";
	modo->art_file = NULL;
	modo->art_separate = 0;

	modo->max_K = 1000;
	modo->min_K = 2;

	modo->param_estimate = PARAM_HAPLOTYPE | PARAM_DELTA | PARAM_PI
		| PARAM_LAMBDA | PARAM_GAMMA | PARAM_BG_PI | PARAM_BETA;

	modo->use_aic = 0;	/* [TODO?] add separate copy to initialization */

	modo->previous_ll = -INFINITY;

	modo->epsilon = 1e-6;
	modo->tolerance = 1e-12;
	modo->beta_epsilon = 1e-2;//0.1;//
	modo->em_info = SILENT;
	modo->n_iter = 100;

	return NO_ERROR;
} /* make_model_options */

/**
 * Process command line option -p or -sp.
 *
 * @param argv	command line arguments and options
 * @param i	index of -p or -sp on command line
 * @param modo	model_options object pointer
 * @return	index of last option to command line argument
 */
int process_parameterization_option(const char *str, const int argc,
				const char **argv, int i, model_options *modo)
{
	if (!strncmp(argv[i + 1], "mlogit", 6)) {
		modo->parameterization = MLOGIT_PARAMETERIZATION;
		mmessage(INFO_MSG, NO_ERROR, "%s parameterization: mlogit\n",
									str);
	} else if (!strncmp(argv[i + 1], "qualit", 6)) {
		modo->parameterization = QUALITY_PARAMETERIZATION;
		mmessage(INFO_MSG, NO_ERROR, "%s parameterization: quality\n",
									str);
	} else if (!strncmp(argv[i + 1], "dada2", 5)) {
		modo->parameterization = DADA2_PARAMETERIZATION;
		mmessage(INFO_MSG, NO_ERROR, "%s parameterization: dada2\n",
									 str);
		if (i + 2 < argc && argv[i + 2][0] != '-') {
			++i;
			modo->dada2_file = argv[i + 1];
			mmessage(INFO_MSG, NO_ERROR, "dada2 error profile "
					"file: %s\n", modo->dada2_file);
		}
	} else if (!strncmp(argv[i + 1], "art", 3)) {
		modo->parameterization = ART_PARAMETERIZATION;
		mmessage(INFO_MSG, NO_ERROR, "%s parameterization: art\n",
									str);
		if (i + 2 < argc && argv[i + 2][0] != '-') {
			++i;
			modo->art_file = argv[i + 1];
			mmessage(INFO_MSG, NO_ERROR, "ART error profile "
					"file: %s\n", modo->art_file);
			if (i + 2 < argc && argv[i + 2][0] != '-') {
				++i;
				modo->art_separate = 1;
			}
		}
	} else if (!strncmp(argv[i + 1], "defaul", 6)) {
		mmessage(INFO_MSG, NO_ERROR, "%s parameterization: default\n",
									str);
	} else {
		return i;
	}

	return i + 1;
} /* process_parameterization_option */


/**
 * Process command line option -e or -se.
 *
 * @param argv	command line arguments and options
 * @param i	index of -p or -sp on command line
 * @param modo	model_options object pointer
 * @return	index of last option to command line argument
 */
int process_estimation_option(const char **argv, int i, model_options *modo)
{
	++i;
	if (!strcmp(argv[i], "delta")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_DELTA;
		mmessage(INFO_MSG, NO_ERROR, "Estimate delta? %s\n",
			(modo->param_estimate & PARAM_DELTA) ? "yes" : "no");
	} else if (!strcmp(argv[i], "gamma")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_GAMMA;
		mmessage(INFO_MSG, NO_ERROR, "Estimate gamma? %s\n",
			(modo->param_estimate & PARAM_GAMMA) ? "yes" : "no");
	} else if (!strcmp(argv[i], "lambda")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_LAMBDA;
		mmessage(INFO_MSG, NO_ERROR, "Estimate lambda? %s\n",
			(modo->param_estimate & PARAM_LAMBDA) ? "yes" : "no");
	} else if (!strcmp(argv[i], "beta")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_BETA;
		mmessage(INFO_MSG, NO_ERROR, "Estimate beta? %s\n",
			(modo->param_estimate & PARAM_BETA) ? "yes" : "no");
	} else if (!strcmp(argv[i], "haplotype")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_HAPLOTYPE;
		mmessage(INFO_MSG, NO_ERROR, "Estimate haplotypes? %s\n",
			(modo->param_estimate & PARAM_HAPLOTYPE) ? "yes" : "no");
	} else if (!strcmp(argv[i], "pi")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_PI;
		mmessage(INFO_MSG, NO_ERROR, "Estimate pi? %s\n",
			(modo->param_estimate & PARAM_PI) ? "yes" : "no");
	} else if (!strcmp(argv[i], "bg_pi")) {
		modo->param_estimate = modo->param_estimate & ~PARAM_BG_PI;
		mmessage(INFO_MSG, NO_ERROR, "Estimate pi? %s\n",
			(modo->param_estimate & PARAM_PI) ? "yes" : "no");
	} else {
		mmessage(ERROR_MSG, NO_ERROR, "Fail to recognize -e argument "
							"'%s'\n", argv[i]);
		return i - 1;
	}
	return i;
} /* process_estimation_option */
