#include <stdlib.h>

#include "initialize_options.h"
#include "initialize.h"
#include "cmdline.h"
#include "error.h"

/**
 * Create initialize_options struct.
 *
 * @param inio_in	address of pointer to initialize_options object
 * @return		error status
 */
int make_initialize_options(initialize_options **inio_in, model_options *modo)
{

	initialize_options *io;
	*inio_in = malloc(sizeof **inio_in);

	if (!*inio_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initialization_options");

	io = *inio_in;

	io->modo = modo;
	io->estimate_K = 0;
	io->n_init = 1;
	io->initialization_method = INIT_RANDOM_SEEDS;	//INIT_KMODES;
	io->beta_epsilon = 0.1;	/* model_options:beta_epsilon is too stringent
				 * and very slow during initialization, which
				 * is why we have this second copy of the 
				 * parameter.simulation_model
				 */

	/* nested initializations */
	io->run_method = EXHAUSTIVE_EM;
	io->n_inner_init = 1;

	/* k-modes initialization */
	io->n_kmodes_iter = 100;	// 10000;
	io->kmodes_initialization_method = KMODES_INIT_RANDOM_SEEDS;
	io->kmodes_opt = NULL;
	io->kmodes_huang97 = 1;

	/* user initialization */
	io->initialization_file = NULL;
	io->partition_file = NULL;

	/* turning initialization into parameter estimates */
	io->assignment_method = HAMMING;

	/* for impro_deblur and ampliCI */
	io->stochastic_deplur = 0;
	io->convergence_deplur = 1;
	io->error_adjust = 1.;
	io->evaluate_initialization = 0;
	io->check_false_positive = 1;
	io->remove_gamma = 0;
	io->use_error_profile = 0;
	io->error_rate = 0.005;
	io->low_bound = 1.5;

	return NO_ERROR;
} /* make_initialize_options */


/**
 * Process command line option -i (or -si).
 *
 * @param argv	command line options and arguments
 * @param i	index of -i or -si
 * @param opt	pointer to master options object
 * @param inio	pointer to initialize_options object
 * @return	index of last argument of option
 */
int process_initialization_option(int argc, const char **argv, int i, void *opt,
						initialize_options *inio)
{
	/* option argument starts with a digit */
	if (argv[i+1][0] >= 48 && argv[i+1][0] <= 57) {
		inio->n_init = read_uint(argc, argv, ++i, opt);
		debug_msg(1, 0, "Number initializations: %u\n", inio->n_init);
	} else if (!strcmp(argv[i + 1], "kmodes-old")) {
		inio->initialization_method = INIT_HW_KMODES;
		++i;
	} else if (!strcmp(argv[i + 1], "kmodes")) {
		inio->initialization_method = INIT_HW_KMODES;
		++i;
	} else if (!strcmp(argv[i + 1], "true")) {
		inio->initialization_method = INIT_TRUE_VALUES;
		++i;
	} else if (!strncmp(argv[i + 1], "truep", 5)) {
		inio->initialization_method = INIT_TRUE_PARTITION;
		++i;
	} else if (!strcmp(argv[i + 1], "clb09")) {
		inio->initialization_method = INIT_CAO;
		++i;
	} else if (!strcmp(argv[i + 1], "deblur")) {
		inio->initialization_method = INIT_DEBLUR;
		++i;
	} else if (!strcmp(argv[i + 1], "imdeblur")) {
		inio->initialization_method = INIT_IMDEBLUR;
		++i;
	} else if (!strcmp(argv[i + 1], "quality")) {
		inio->assignment_method = QUALITY;
		++i;
	} else if (!strcmp(argv[i + 1], "hamming")) {
		inio->assignment_method = HAMMING;
		++i;
	} else if (!strncmp(argv[i + 1], "random", 6)) {
		inio->initialization_method = INIT_RANDOM_SEEDS;
		inio->kmodes_initialization_method = KMODES_INIT_RANDOM_SEEDS;
		inio->n_kmodes_iter = 0;
		inio->kmodes_huang97 = 0;
		++i;
	} else {
		++i;
		for (size_t j = 0; j < strlen(argv[i]); ++j)
			if (argv[i][j] == ':') {
				inio->initialization_file = malloc((j + 1)
					* sizeof *inio->initialization_file);
				strncpy(inio->initialization_file, argv[i], j);
				inio->partition_file = &argv[i][j+1];
				debug_msg(1, 0, "partition file: '%s'\n",
							inio->partition_file);
				break;
			}
		if (!inio->initialization_file) {
			inio->initialization_file = malloc((strlen(argv[i]) + 1)
					* sizeof *inio->initialization_file);
			strcpy(inio->initialization_file, argv[i]);
			mmessage(INFO_MSG, NO_ERROR, "Initialization file:  "
				"'%s'\n", inio->initialization_file);
		}
	}

	return i;
} /* process_initialization_option */
