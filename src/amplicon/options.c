/**
 * @file options.c
 * @author Karin S. Dorman
 *
 * Command-line options handling.
 *
 * DONE
 * X take a set of seeds to initialize (-i <seed_file_as_fasta>)
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>



#include "options.h"
#include "initialize.h"
#include "simulate.h"
#include "model.h"
#include "kmodes.h"
#include "cmdline.h"
#include "io.h"

#ifdef USE_CURSES
WINDOW *global_wp;
#else
FILE *global_wp = NULL;	/* never changed */
#endif
FILE *active_fp;

/**
 * Setup options object.
 */
int make_options(options **opt) {
	int err = NO_ERROR;
	options *op;
	*opt = malloc(sizeof **opt);

	if (*opt == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");

	op = *opt;

	op->seed = 0;
	op->fastq_file = NULL;
	op->outfile = NULL;
	op->outfile_k=NULL;
	op->seed_file = NULL;

	op->multi_stage_method = 0;
	op->sample_size = 1000;
	op->simu_size = 10000;
	op->proportion = 1.00;
	op->epsilon_s = 1;

	op->modo = NULL;	/* estimation model options */
	op->inio = NULL;	/* estimation initialization options */
	op->simo = NULL;	/* simulation options */
	op->sim_modo = NULL;	/* simulation model options */
	op->sim_inio = NULL;	/* simulation initialization options */

	if ((err = make_model_options(&op->modo)))
		return err;

	if ((err = make_model_options(&op->sim_modo)))
		return err;

	if ((err = make_initialize_options(&op->inio, op->modo)))
		return err;

	if ((err = make_initialize_options(&op->sim_inio, op->sim_modo)))
		return err;

	if ((err = make_simulate_options(&op->simo, op->sim_inio, op->sim_modo)))
		return err;

	op->do_estimation = 1;
	op->do_simulation = 0;
	op->estimate_K = 0;
	op->reference_file = NULL;
	op->offset_file = NULL;
	op->use_curses = 0;
	
	return NO_ERROR;
} /* make_options */

/**
 * Free options object.
 *
 * @param opt	pointer to options object
 */
void free_options(options *opt)
{
	if(opt) {
		if (opt->inio->partition_file && opt->inio->initialization_file)
			free(opt->inio->initialization_file);
		if (opt->simo->ancestor)
			free(opt->simo->ancestor);
		if (opt->simo->simulation_infile)
			free(opt->simo->simulation_infile);
		free(opt);
	}
} /* free_options */

/**
 * Parse command-line.
 */
int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'a':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				if (!strncmp(&argv[i][j], "anc", 3)) {
					opt->simo->ancestor = malloc((strlen(
						argv[++i]) + 1) * sizeof
						*opt->simo->ancestor);
					if (!opt->simo->ancestor) {
						err = mmessage(ERROR_MSG,
							MEMORY_ALLOCATION,
							"options.ancestor");
						goto CMDLINE_ERROR;
					}
					strcpy((char *) opt->simo->ancestor, argv[i]);
				} else if (!strncmp(&argv[i][j], "api", 3)) {
					opt->simo->sim_pi_alpha = read_cmdline_double(
						argc, argv, ++i, (void *)opt);
				} else if (!strncmp(&argv[i][j], "aga", 3)) {
					opt->simo->sim_gamma_alpha = read_cmdline_double(
						argc, argv, ++i, (void *)opt);
				} else if (!strncmp(&argv[i][j], "art", 3)) {
					opt->modo->parameterization = ART_PARAMETERIZATION;
					opt->modo->art_file = argv[++i];
					debug_msg(1, 0, "art file: %s\n",
						opt->modo->art_file);
					if (i + 1 < argc && argv[i + 1][0] != '-') {
						opt->modo->art_separate = 1;
						++i;
					}
				} else {
					opt->reference_file = argv[++i];
				}
				break;
			case 'b':
				opt->modo->background_model = 1;
				break;
			case 'n':
				opt->estimate_K = 1;
				break;
			case 'y':
				opt->multi_stage_method = 1;
				break;
			case 'x':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				else if (argv[i+1][0] >= 48
					&& argv[i+1][0] <= 57)
					opt->sample_size = read_uint(argc, argv, ++i,
						(void *)opt);
				break;
			case 'k':
				if (i == argc - 1) {
					goto CMDLINE_ERROR;
				} else if (!strcmp(&argv[i][j], "kmin")) {
					if (argv[i + 1][0] >= 48
							&& argv[i + 1][0] <= 57)
						opt->modo->min_K = read_uint(
							argc, argv, ++i,
								(void *)opt);
				} else if (!strcmp(&argv[i][j], "kmax")) {
					if (argv[i + 1][0] >= 48
							&& argv[i + 1][0] <= 57)
						opt->modo->max_K = read_uint(
							argc, argv, ++i,
							(void *)opt);
				} else if (!strcmp(&argv[i][j], "ktrue")) {
					if (argv[i + 1][0] >= 48
							&& argv[i + 1][0] <= 57)
						opt->simo->true_K = read_uint(argc,
							argv, ++i, (void *)opt);
					debug_msg(1, 0, "true K: %u\n",
								opt->simo->true_K);
				} else if (argv[i+1][0] >= 48
					&& argv[i+1][0] <= 57) {
					opt->modo->K = read_uint(argc, argv,
							++i, (void *)opt);
				} else if (!strcmp(argv[i + 1], "CLB09")) {
					opt->inio->kmodes_initialization_method
						= KMODES_INIT_CLB09_RANDOM;
					++i;
				} else if (!strcmp(argv[i + 1], "H97")) {
					opt->inio->kmodes_initialization_method
						= KMODES_INIT_H97_RANDOM;
					++i;
				} else if (!strcmp(argv[i + 1], "HD17")) {
					opt->inio->kmodes_initialization_method
						= KMODES_INIT_HD17;
					++i;
				}
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'l':
				opt->modo->previous_ll = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
				break;
			case 'c':
				opt->simo->simulation_control |= SIMULATE_HAPLOTYPES;
				if (!strcmp(&argv[i][j], "control"))
					opt->simo->hap_control_file = argv[++i];
				else
					opt->simo->haplotype_spread =
						read_cmdline_double(argc, argv,
							 ++i, (void *)opt) / 2;
				break;
			case 'v':
				opt->simo->simulation_control |= ADJUST_ERROR;
				opt->simo->cluster_spread = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
				if (opt->simo->cluster_spread <= 0 ||
					opt->simo->cluster_spread >= 2) {
					mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						"-v argument must be positive.\n");
					goto CMDLINE_ERROR;
				}
				break;
			case 'g':
				opt->simo->simulation_control |= ADJUST_QUALITY;
				opt->simo->quality_change = read_int(argc, argv,
					++i, (void *)opt);
				break;
			case 'd':
				if (i == argc - 1) {
					goto CMDLINE_ERROR;
				} else if (!strncmp(&argv[i][j], "dada2", 5)) {
					opt->modo->parameterization
						= DADA2_PARAMETERIZATION;
					opt->modo->dada2_file = argv[++i];
					debug_msg(1, 0, "dada2 file: %s\n",
						opt->modo->dada2_file);
					break;
				} else {
					opt->simo->fastq_outfile = argv[++i];
					debug_msg(1, 0, "fastq outfile: %s\n",
						opt->simo->fastq_outfile);
				}
				break;
			case 'e':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				if (i == process_estimation_option(argv, i,
								opt->modo))
					goto CMDLINE_ERROR;
				else
					++i;
				break;
			case 'i':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				i = process_initialization_option(argc, argv, i,
							(void *)opt, opt->inio);
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'f':
				if (i == argc - 1)
					err = INVALID_CMD_OPTION;
				opt->fastq_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Input fastq file:"
					" '%s'\n", opt->fastq_file);
				while (i + 1 < argc && argv[i + 1][0] != '-') {
					++i;
					if (!strncmp(argv[i], "mix", 3))
						opt->simo->simulation_from_data |= MIXING;
					else if (!strncmp(argv[i], "qual", 4))
						opt->simo->simulation_from_data |= QUALITIES;
					else if (!strncmp(argv[i], "nuc", 3))
						opt->simo->simulation_from_data |= NUCLEOTIDES;
					else if (!strncmp(argv[i], "hap", 3))
						opt->simo->simulation_from_data |= HAPLOTYPES;
					else
						--i;
				}
				mmessage(INFO_MSG, NO_ERROR, "Mixing "
					"proportions %s.\n",
					opt->simo->simulation_from_data & MIXING
						? "from data" : "simulated");
				mmessage(INFO_MSG, NO_ERROR, "Quality "
					"scores %s.\n",
					opt->simo->simulation_from_data &
					QUALITIES ? "from data" : "simulated");
				mmessage(INFO_MSG, NO_ERROR, "Nucleotide "
					"emission parameters %s.\n",
					opt->simo->simulation_from_data &
					NUCLEOTIDES ? "from data" : "simulated");
				mmessage(INFO_MSG, NO_ERROR, "Haplotypes %s.\n",
					opt->simo->simulation_from_data &
					HAPLOTYPES ? "from data" : "simulated");
				break;
			case 'o':
				if (i == argc - 1)
					err = INVALID_CMD_OPTION;
				opt->outfile = argv[++i];
				break;
			case 'p':
				if (i == argc - 1)
					err = INVALID_CMD_OPTION;

				i = process_parameterization_option("model",
						argc, argv, i, opt->modo);
				break;
			case 'q':
				opt->modo->model_quality = 0 ; // MODIFY IT LATER
				break;
			case 'r':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				else if (argv[i+1][0] >= 48
					&& argv[i+1][0] <= 57) {
					opt->seed = read_uint(argc, argv, ++i,
						(void *)opt);
					srand(opt->seed);
					debug_msg(1, 0, "seed: %lu\n", opt->seed);
				} else if (!strncmp(argv[i + 1], "ini-em", 6)) {
					opt->inio->run_method = INI_EM;
					++i;
					opt->inio->n_inner_init = read_uint(argc,
							argv, ++i, (void *)opt);
				} else if (!strncmp(argv[i + 1], "rnd-em", 6)) {
					opt->inio->run_method = RND_EM;
					++i;
					opt->inio->n_inner_init = read_uint(argc,
							argv, ++i, (void *)opt);
				} else if (!strncmp(argv[i + 1], "em-em", 5)) {
					err = INVALID_CMD_OPTION;
					mmessage(ERROR_MSG, err, "'em-em' "
						"initialization not yet implemented");
					++i;
					goto CMDLINE_ERROR;
				} else {
					while (i < argc && argv[i + 1][0] != '-') {
						if (!strncmp(argv[i + 1], "nuc", 3)) {
							opt->simo->simulation_random ^= NUCLEOTIDES;
							mmessage(INFO_MSG, NO_ERROR, "Simulate "
								"random nucleotides: %s.\n",
								opt->simo->simulation_random
								& NUCLEOTIDES ? "yes" : "no");
							++i;
						} else if (!strncmp(argv[i + 1], "mix", 3)) {
							opt->simo->simulation_random ^= MIXING;
							mmessage(INFO_MSG, NO_ERROR, "Simulate "
								"random mixing proportions: "
								"%s.\n",
								opt->simo->simulation_random
								& MIXING ? "yes" : "no");
							++i;
						} else if (!strncmp(argv[i + 1], "qual", 4)) {
							opt->simo->simulation_random ^= QUALITIES;
							mmessage(INFO_MSG, NO_ERROR, "Simulate "
								"random qualities: %s.\n",
								opt->simo->simulation_random
								& QUALITIES ? "yes" : "no");
							++i;
						} else if (!strncmp(argv[i + 1], "hap", 3)) {
							opt->simo->simulation_random ^= HAPLOTYPES;
							mmessage(INFO_MSG, NO_ERROR, "Simulate "
								"random haplotypes: %s.\n",
								opt->simo->simulation_random
								& HAPLOTYPES ? "yes" : "no");
							++i;
						}
					}
				}
				break;
			case 's':
				if (i == argc - 1)	/* requires >0 args */
					goto CMDLINE_ERROR;

				opt->do_simulation = 1;

				/* initialization for simulation */
				if (!strncmp(&argv[i][j], "si", 2)) {
					i = process_initialization_option(argc,
							argv, i, (void *) opt,
								opt->sim_inio);
					if (errno)
						goto CMDLINE_ERROR;
					break;
				} else if (!strncmp(&argv[i][j], "sp", 2)) {

					i = process_parameterization_option(
						"simulation", argc, argv, i,
								opt->sim_modo);
					if (opt->sim_modo->parameterization
						== ART_PARAMETERIZATION) {
						opt->simo->simulation_from_data
							&= ~QUALITIES;
						opt->simo->simulation_random
							|= QUALITIES;
					}
					break;
				} else if (!strncmp(&argv[i][j], "sk", 2)) {
					opt->sim_modo->K = read_uint(argc, argv,
							++i, (void *)opt);
					mmessage(INFO_MSG, NO_ERROR,
						"Simulation K: %u\n",
							opt->sim_modo->K);
					break;
				}

				++i;
				for (size_t j = 0; j < strlen(argv[i]); ++j)
					if (argv[i][j] == ':') {
						opt->simo->simulation_infile = malloc((j + 1)
							* sizeof *opt->simo->simulation_infile);
						strncpy(opt->simo->simulation_infile, argv[i], j);
						opt->simo->simulation_infile[j] = '\0';
						opt->inio->partition_file = &argv[i][j+1];
						debug_msg(1, 0, "partition "
							"file: '%s'\n",
							opt->inio->partition_file);
						break;
					}
				if (is_numeric(argv[i])) {
					if (i == argc - 1)	/* requires 2 args */
						goto CMDLINE_ERROR;
					opt->simo->sim_n_reads = read_uint(argc, argv,
							i, (void *)opt);
					opt->simo->sim_read_length = read_uint(argc,
							argv, ++i, (void *)opt);
					if (i + 1 < argc && argv[i + 1][0] != '-')
						opt->simo->simulation_outfile
								= argv[++i];
					debug_msg(1, 0, "Simulate %u reads, %u "
						"positions, outputting results "
						"in '%s'\n", opt->simo->sim_n_reads,
						opt->simo->sim_read_length,
						opt->simo->simulation_outfile);
				} else if (!strcmp(argv[i], "quality")) {
					opt->sim_modo->parameterization
						= QUALITY_PARAMETERIZATION;
					++i;
				} else if (!strcmp(argv[i], "ampliclust")) {
					opt->sim_modo->parameterization
						= DEFAULT_PARAMETERIZATION;
					++i;
				} else if (!strcmp(argv[i], "dada2")){
					opt->sim_modo->parameterization
						= DADA2_PARAMETERIZATION;
					++i;
				} else if (!strcmp(argv[i], "art")){
					opt->sim_modo->parameterization
						= ART_PARAMETERIZATION;
					++i;
				} else if (!strcmp(argv[i], "mlogit")){
					opt->sim_modo->parameterization
						= MLOGIT_PARAMETERIZATION;
					++i;
				} else if (!opt->simo->simulation_infile) {
					if (i == argc - 1
						|| argv[i + 1][0] == '-') {
						opt->simo->simulation_outfile = argv[i];
					} else {
						opt->simo->simulation_infile = malloc(
							(strlen(argv[i]) + 1) *
									sizeof
							*opt->simo->simulation_infile);
						strcpy(opt->simo->simulation_infile,
									argv[i]);
						debug_msg(1, 0, "simulation infile: "
							"'%s'\n",
							opt->simo->simulation_infile);
						opt->simo->simulation_outfile = argv[++i];
					}
					debug_msg(1, 0, "simulation outfile: "
						"'%s'\n",
						opt->simo->simulation_outfile);
				} else {
					opt->simo->simulation_outfile = argv[++i];
				}
				break;
			case 't':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				opt->offset_file = argv[++i];
				break;
			case 'w':
				opt->use_curses = 1;
				break;
			case 'h':
				fprint_usage(stderr, argv[0], opt);
				free_options(opt);
				exit(EXIT_SUCCESS);
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}

	/* user is requesting a specific initialization */
	/* I think this is no longer relevant
	if (opt->inio->initialization_file)
		opt->inio->n_inner_init = 1;
	*/
	if (!opt->inio->n_init)
		opt->do_estimation = 0;

	/* cannot initialize with truth if not simulating */
	if (!opt->do_simulation) {
		if (opt->inio->initialization_method == INIT_TRUE_VALUES
			|| opt->inio->initialization_method == INIT_TRUE_PARTITION)
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
				"initialize with truth when not simulating.\n");
	} else {
		if (!opt->simo->true_K)
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"You must set true_K for the simulation!\n");
	}

	if (opt->inio->initialization_method == INIT_TRUE_VALUES
		|| opt->inio->initialization_method == INIT_TRUE_PARTITION) {
		if (opt->inio->n_init > 1)
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Repeated"
				" initializations senseless when initializating"
				" with truth.\n");
		if (!opt->modo->K) {
			opt->modo->K = opt->simo->true_K;
			debug_msg(1, 0, "Setting K = %u\n", opt->modo->K);
		}
		if (opt->simo->true_K != opt->modo->K)
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
				"initialize with truth when simulation K (%u) "
				"not same as estimation K (%u).\n", opt->simo->true_K,
				opt->modo->K);
	}

	if (opt->do_simulation && opt->modo->parameterization == ART_PARAMETERIZATION
		&& !opt->modo->art_file && opt->sim_modo->art_file)
		opt->modo->art_file = opt->sim_modo->art_file;

	/* If the user does not give a K and does not use option -n */
	if (!opt->modo->K && !opt->estimate_K && opt->inio->n_init > 0)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"You must choose K to use during estimation or a "
			"method to estimate K (see option -n)!\n");

	/* if user does not provide estimation K for simulation */
	if (opt->do_simulation && !opt->sim_modo->K)
		opt->sim_modo->K = opt->simo->true_K;

	/* If the user gives a K but uses option -n by mistake */
	if (opt->modo->K && opt->estimate_K)
		opt->estimate_K = 0;

	if(opt->modo->min_K > opt->modo->max_K && opt->estimate_K)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
					 "Invalid --kmin and --kmax.\n");
	
	if(opt->modo->min_K == opt->modo->max_K) {
		opt->estimate_K = 0;
		opt->modo->K = opt->modo->min_K;
	}

	/* K identified in different stages should not be constant */
	if (opt->multi_stage_method)
		opt->estimate_K = 1;

	/* user choose not to model quality */
	if (!opt->modo->model_quality)
		opt->modo->param_estimate = PARAM_HAPLOTYPE | PARAM_DELTA
				| PARAM_PI | PARAM_GAMMA | PARAM_BG_PI;

	if (!opt->modo->background_model)	/* [KSD] I don't think this is necessary. */
		opt->modo->param_estimate = opt->modo->param_estimate & ~PARAM_BG_PI;
		
	/* If the user does not choose a K and want to find a true K */
	if (!opt->modo->K && opt->estimate_K)
		opt->modo->K = opt->modo->min_K;

	if (opt->do_simulation && opt->simo->simulation_random & MIXING
		&& !(opt->simo->simulation_random & NUCLEOTIDES))
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You cannot "
			"randomize the mixing proportions without randomizing "
			"the reads.\n");

	if (opt->do_simulation && opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION &&
		!(opt->simo->simulation_random & QUALITIES))
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You cannot "
			"simulated under the ampliclust model without "
			"simulating quality scores.\n");

	if (opt->do_simulation && opt->sim_modo->parameterization == QUALITY_PARAMETERIZATION
		&& !(opt->simo->simulation_random & NUCLEOTIDES))
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You cannot "
			"simulate under the quality model and not simulate "
			"reads.\n");

	if (opt->do_simulation
		&& opt->sim_modo->parameterization == ART_PARAMETERIZATION
		&& !opt->sim_modo->art_file)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You must provide"
			" an ART profile file to simulate with ART model.\n");

	if (opt->do_simulation && opt->simo->simulation_control & ADJUST_ERROR
		&& (!(opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION)
		|| !(opt->simo->simulation_random & NUCLEOTIDES)))
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You cannot "
			"simulate while controlling error rates without "
			"using ampliclust model and simulating reads.\n");

	if (opt->inio->initialization_method == INIT_KMODES) {
		if ((err = make_kmodes_options(&opt->inio->kmodes_opt)))
			goto CMDLINE_ERROR;
		opt->inio->kmodes_opt->weighted = 0;
		opt->inio->kmodes_opt->init_update = opt->inio->kmodes_huang97;
		opt->inio->kmodes_opt->use_qtran = 0;
		opt->inio->kmodes_opt->use_hartigan = 0;
	}

	if (opt->inio->initialization_method == INIT_TRUE_VALUES
		&& opt->modo->parameterization != opt->sim_modo->parameterization)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "You cannot "
			"initialize with true values (-i true) if simulation "
			"and estimation models do not match.\n");

	if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION
		&& opt->modo->background_model)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot combine "
			"mlogit model and background model (yet).\n");
	
	if (opt->fastq_file && opt->do_simulation && opt->simo->sim_n_reads)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot input "
			"fastq file and specify number and length of reads to "
			"simulate; it will use the number and lengths in the "
			"fastq file.\n");

	if (opt->do_simulation && !err)
		set_seed(rand(), rand());

	return err;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

/**
 * Usage statement.
 */
void fprint_usage(FILE *fp, const char *cmdname, void *obj)
{
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - cluster amplicon sequences\n",
		&cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [-a <astr> -r <rulong> -i <istr> "
		"--run <runstr> <runuint> --kinit <kinitstr> "
		"-s <sstr1> [<sstr2>] -l <ldbl> "
		"-c <cdbl> -v <vdbl> -o <ostr> -d <dstr> -t <tstr>] "
		"-k <kuint> -f <fstr>...\n",
		&cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s clusters amplicon sequences presented "
		"as reads in fastq file <fstr> into <kuint> clusters.  It "
		"records the best solution over <iuint> initializations "
		"in file <ostr> if that solution achieves better log "
		"likelihood than than <ldbl>."
		"\n\n\tBy default, it uses k-modes initialization, itself "
		"initialized with method <kinitstr>, but you can request random"
		" seeds initialization by setting <initstr>=random.  It can "
		"also use RND-EM with <inituint> random initializations if "
		"<initstr>=rnd-em.  The randomness of initialization is "
		"controlled by setting random number seed <rulong>."
		"\n\n\tIt can simulate from partitioned data in file "
		"<sstr>, writing fastq to <dstr>, before running regularly.\n",
		&cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-a <astr>\n\t\tAlign reads to this \"ancestor\" sequence [implementation not complete].\n");
	fprintf(fp, "\t--art <art_str> [separate]\n\t\tUse ART model (no indels), where <art_str> is ART profile file. Implies -p art. [DEFAULT: none]\n");
	fprintf(fp, "\t\tOptional second argument if selected ART profile uses different profile for each haplotype nucleotide.\n");
	fprintf(fp, "\t--api <adbl>\n\t\tSimulate mixing proportions as Dirichlet(<adbl>, <adbl>, ...).\n");
	fprintf(fp, "\t--agamma <adbl>\n\t\tSimulate substitution proportions as Dirichlet(<adbl>, <adbl>, ...).\n");
	fprintf(fp, "\t--anc <astr>\n\t\tUse this \"ancestor\" sequence for simulation.\n");
	fprintf(fp, "\t-b\n\t\tInclude a background scatter model.  Specifically, an unknown proportion of reads are assumed to be generated by an iid Unif(piN) model, where piN is estimated.\n");
	fprintf(fp, "\t-c <cdbl>\n\t\tSimulate haplotype centers with average expected number of substitutions per site <cdbl> from each other (currently assumes JC69 model on a star phylogeny).\n");
	fprintf(fp, "\t--control <cstr>\n\t\tSimulate haplotype centers with changes only at positions marked 0 (fixed) / 1 (mutable) in file <cstr>.\n\t\tMust use with -c <cdbl>.\n");
	fprintf(fp, "\t-v <vdbl>\n\t\tScale within-cluster \"variance\" during simulation with respect to a fastq file (see -f).\n");
	fprintf(fp, "\t\tIf 0 < <vdbl> < 1, scale probability of error at each site by this proportional amount.  If 1 < <vdbl> < 2, scale within-cluster probability of no error by <vdbl> -  1.\n");
	fprintf(fp, "\t-g <gint>\n\t\tRough scaling of within-cluster \"variance\" during simulation with respect to a fastq file (see -f).\n");
	fprintf(fp, "\t\tShift quality scores by this integer amount, while respecting boundaries [0,40].\n");
	fprintf(fp, "\t-d <dstr>\n\t\tFastq output file for simulated data.\n");
	fprintf(fp, "\t-e <estr>\n\t\tUse multiple times to list parameters to NOT estimate (delta|gamma|lambda|beta|haplotype|pi|bg_pi).\n");
	fprintf(fp, "\t-f <fstr> [<fstr1> <fstr2> ...]\n\t\tThe fastq input file. [OPTIONAL]\n");
	fprintf(fp, "\t\tOptional arguments toggle simulation parameters to be estimated from data (mixing|nucleotides|qualities|haplotypes) [DEFAULT: qualities)\n");
	fprintf(fp, "\t\tmixing\tmixing proportions are not taken from data by default\n");
	fprintf(fp, "\t\tnucleotides\tnucleotide emission parameters are not taken from data by default\n");
	fprintf(fp, "\t\tqualities\tqualities are taken from data by default\n");
	fprintf(fp, "\t\thaplotypes\thaplotypes are not taken from data by default\n");
	fprintf(fp, "\t-i, --init <istr>\n\t\tSet desired number of initializations [OPTIONAL; DEFAULT: %u] or set up type of initialization...\n"
		"\t\t\tidentify partition file or seed file or haplotype file or haplotype_file:partition_file or ...\n"
		"\t\t\tchoose \"kmodes\" for k-modes, \"true\" for true parameter values, \"truepartition\" for true partition, \"clb09\" for Cao, \"deblur\" for deblur, \"imdeblur\" for imdeblur, or \"random\" for random [OPTIONAL; DEFAULT: k-modes] or ...\n"
		"\t\t\tchoose \"hamming\" for Hamming distance or \"quality\" for quality model to assign reads to modes [OPTIONAL; DEFAULT: hamming].\n",
		 opt->inio->n_init);
	fprintf(fp, "\t\tA partition file has the 0-based index of the cluster to which each read belongs on one or separate lines (see the output of data::optimal_cluster_id on \"assignments\" line of output file [-o argument]).\n");
	fprintf(fp, "\t\tA seed file lists K 0-based indices of reads that will serve as the haplotypes (cluster centers).\n");
	fprintf(fp, "\t\tA haplotype file is a fasta-formatted file with the K haplotypes.\n");
	fprintf(fp, "\t\tYou can combine a haplotype file and partition file by placing a ':' (colon) with no spaces between their names.\n");
	fprintf(fp, "\t--run <runstr> <runuint>\n\t\tRun method and number of inner initializations for ini-em|rnd-em|em-em (exhaustive|ini-em|rnd-em|em-em). [DEFAULT: exhaustive; em-em not yet implemented]\n");
	fprintf(fp, "\t-k <kuint>\n\t\tSet number of clusters K to estimate.\n");
	fprintf(fp, "\t--ktrue <ktrueuint>\n\t\tSet simulated number of clusters K.\n");
	fprintf(fp, "\t--kinit <kinitstr>\n\t\tK-modes initialization method (rnd|H97|CBL09|HD17). [DEFAULT: rnd]\n");
	fprintf(fp, "\t-l <ldbl>\n\t\tPrevious best log likelihood solution that a new solution is required to beat in order to save the results in outfile (-o argument).  [DEFAULT: none]\n");
	fprintf(fp, "\t-o <ostr>\n\t\tOutput file to record results of best solution.  [DEFAULT: none]\n");
	fprintf(fp, "\t-p <pstr>\n\t\tSet estimation model parameterization (default|mlogit|dada2|art|quality) [DEFAULT: default].\n");
	fprintf(fp, "\t\t\tampliclust: model with delta, gamma, lambda0 and lambda1.\n");
	fprintf(fp, "\t\t\tmlogit: model with beta.\n");
	fprintf(fp, "\t\t\tquality: model with gamma and quality scores taken literally.\n");
	fprintf(fp, "\t\t\tdada2: model with (hard-coded) error profile taken from DADA2.\n");
	fprintf(fp, "\t\t\tart: model with (hard-coded) quality profile taken from ART, gamma and quality scores taken literally.\n");
	fprintf(fp, "\t-q\n\t\tModel quality scores in ampliclust or art model.  [DEFAULT: not selected]\n");
	fprintf(fp, "\t-r <rulong>\n\t\tSet random number seed. [OPTIONAL]\n");
	fprintf(fp, "\t-r <rstr1> <rstr2> ...\n\t\tRandomize these data during simulation. [DEFAULT: mixing|nucleotides|haplotypes]\n");
	fprintf(fp, "\t\t\tmixing\t\tmixing proportions\n");
	fprintf(fp, "\t\t\tnucleotides\tread nucleotides\n");
	fprintf(fp, "\t\t\tqualities\tread qualities\n");
	fprintf(fp, "\t\t\thaplotypes\thaplotypes\n");
	fprintf(fp, "\t-s <sstr1> <sstr2>\n\t\tSimulate data based on <fstr>, using initialization in <sstr1> (see -i for allowed format), output simulation information in <sstr2>, simulated data in <dstr>.  [DEFAULT: none]\n");
	fprintf(fp, "\t-s <sstr1>\n\t\tSimulate data based on <fstr>, output simulation information in <sstr1>, simulated data in <dstr>.  [DEFAULT: none]\n");
	fprintf(fp, "\t-s <sint1> <sint2>[ <sstr3>]\n\t\tSimulate <sint1> reads of length <sint2> under model <sstr> with K=<ktrueuint>, output simulation information in <sstr3>, simulated data in <dstr>.\n");
	fprintf(fp, "\t-sk <skuint>\n\t\tSet number of clusters K to estimate during simulation initialization.\n");
	fprintf(fp, "\t-sp <spstr>[ <soptional>]\n\t\tSet simulation model parameterization (default|mlogit|dada2|art|quality) [DEFAULT: default].\n");
	fprintf(fp, "\t\tCertain argument values take additional optional arguments:\n");
	fprintf(fp, "\t\t\t-sp sdada2[ <sdada2_file>]\n\t\t\t\tSelect dada2 model with error profile given in <sdada2_file> [DEFAULT: %s].\n", opt->sim_modo->dada2_file);
	fprintf(fp, "\t\t\t-sp sart[ <sart_file>]\n\t\t\t\tSelect ART model with error profile given in <sart_file> [DEFAULT: %s].\n", opt->sim_modo->art_file);
	fprintf(fp, "\t\t\t\tUse of this argument implies '-r qualities' to toggle simulation of quality scores.\n");
	fprintf(fp, "\t-si ...\n\t\tSet initialization method when estimating parameters from data for simulation. See -i.\n");
	fprintf(fp, "\t-t <tstr>\n\t\tFile with offsets for differentially trimmed reads.  [DEFAULT: none]\n");
	fprintf(fp, "\t-w\n\t\tUse curses, which allows user to hit <esc> key to abort the current AECM iterations and proceed to next initialization. [DEFAULT: not set]\n");
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
