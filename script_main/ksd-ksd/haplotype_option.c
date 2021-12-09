/**
 * @file haplotype_option.c
 * @author Yudi Zhang
 */

#include "haplotype_option.h"
#include "error.h"
#include "util.h"
#include "math.h"
#include "fastq.h"
#include "io.h"
#include "io_kmodes.h"
#include "cmdline.h"
#include "myfun.h"

int make_options(options **opt) {
	*opt = malloc(sizeof **opt);
	
	if (*opt == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");
	
	(*opt)->K = 4;	/* invalid value */
	(*opt)->true_column = UINT_MAX;
	(*opt)->run_with_quals = 1;
	(*opt)->true_cluster = NULL;
	(*opt)->true_cluster_size = NULL;
	(*opt)->true_K = 0;
	(*opt)->n_init = 3; /* Used to test the running time */
	(*opt)->n_inner_init = 1;
	(*opt)->n_max_iter = 10000;
	(*opt)->info = QUIET;
	(*opt)->subtract_one = 0;
	(*opt)->shuffle = 0;
	(*opt)->kmodes_algorithm = FASTQ_HW;
	(*opt)->init_method = KMODES_INIT_USER_SEEDS;
	(*opt)->estimate_k = 0;
	(*opt)->kopt = NULL;
	(*opt)->update_modes = 0;
#ifdef __KMODES_NO_QTRANS__
	(*opt)->use_qtran = 0;
#else
	(*opt)->use_qtran = 1;
#endif
	(*opt)->use_hartigan = 0;
	(*opt)->target = 0;
	(*opt)->seed = 3;
	(*opt)->seed_idx = NULL;
	(*opt)->seed_set = NULL;
	(*opt)->n_sd_idx = 0;
	(*opt)->n_seed_set = 0;
	(*opt)->datafile = "sim3.1.fastq";
	(*opt)->data_outfile = NULL;
	(*opt)->ini_file = NULL;
	(*opt)->soln_file = "zoo_Ud_k5_llO.txt";
	(*opt)->pfile = NULL;
	(*opt)->sfile = NULL;
	(*opt)->mfile = NULL;
	(*opt)->mfile_out = NULL;
	(*opt)->weight = KMODES_NO_WEIGHTING;
	(*opt)->quiet = MINIMAL;
	(*opt)->seconds = 0;
	(*opt)->simulate = 0;
	(*opt)->sim_K = 0;
	(*opt)->require_sim_K = 1;
	(*opt)->sim_between_t = 0;
	(*opt)->sim_within_t = 0;
	(*opt)->sim_n_observations = 0;
	(*opt)->sim_n_coordinates = 0;
	(*opt)->sim_n_categories = 0;
	(*opt)->sim_pi = NULL;
	(*opt)->sim_alpha = NULL;
	(*opt)->sim_cluster = NULL;
	(*opt)->sim_modes = NULL;
	(*opt)->result_files = NULL;
	(*opt)->n_result_files = 0;
	(*opt)->true_modes = NULL;
	
	return NO_ERROR;
} /* make_options */

void free_options(options *opt) {
	if (opt) {
		if (opt->result_files) {
			for (unsigned int i = 0; i < opt->n_k; ++i)
				free(opt->result_files[i]);
			free(opt->result_files);
			opt->result_files = NULL;
		}
		if (opt->sim_cluster) {
			free(opt->sim_cluster);
			opt->sim_cluster = NULL;
		} else if (opt->true_cluster) {
			free(opt->true_cluster);
			opt->true_cluster = NULL;
		}
		if (opt->true_cluster_size) {
			free(opt->true_cluster_size);
			opt->true_cluster_size = NULL;
		}
		if (opt->seed_set) {
			if (opt->seed_set[0])
				free(opt->seed_set[0]);
			free(opt->seed_set);
			opt->seed_set = NULL;
		}
		if (opt->sim_pi) {
			free(opt->sim_pi);
			opt->sim_pi= NULL;
		}
		if (opt->sim_alpha) {
			free(opt->sim_alpha);
			opt->sim_alpha= NULL;
		}
		if (opt->sim_modes) {
			if (opt->sim_modes[0])
				free(opt->sim_modes[0]);
			free(opt->sim_modes);
			opt->sim_modes = NULL;
		}
		/* opt->seed_idx free'd via dat->seed_idx */
		free(opt);
	}
} /* free_options */

int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;
	
	opt->n_k = 0;
	opt->min_k = UINT_MAX;
	opt->max_k = 0;
	for (i = 1; i < argc; ++i) {
		j = 0;
		if (argv[i][j] == '-') {
			a = argv[i][++j];
			while (a == '-' && ++j < (int) strlen(argv[i]))
				a = argv[i][j];
			if (argv[i][j] == 'k'
				&& isdigit((unsigned char)argv[i][j+1])) {
				unsigned int k = strtoul(&argv[i][j+1], NULL, 0);
				++opt->n_k;
				if (k < opt->min_k)
					opt->min_k = k;
				else if (k > opt->max_k)
					opt->max_k = k;
			}
		}
	}
	if (opt->n_k && opt->n_k != opt->max_k - opt->min_k + 1)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "-k[0-9] "
				"arguments must identify contiguous values of K.\n");
	if (opt->n_k) {
		debug_msg(MINIMAL, opt->quiet, "%u -k[0-9] arguments\n",
								opt->n_k);
		opt->result_files = malloc(opt->n_k
					* sizeof *opt->result_files);
		opt->n_result_files = malloc(opt->n_k
					* sizeof *opt->n_result_files);
		if (opt->result_files == NULL
					|| opt->n_result_files == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"options:result_files");
	}
	
	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'c':
				if (!strncmp(&argv[i][j], "cont", 4)) {
					opt->continue_run = 1;
					break;
				} else if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->true_column = read_uint(argc, argv, ++i, (void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet, "true column = "
						"%u\n", opt->true_column);
				break;
			case 'j':
				opt->estimate_k = 1;
				break;
			case 'k':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				if (strlen(&argv[i][j]) > 1) {
					unsigned int k = strtoul(&argv[i][j+1],
								 NULL, 0);
					opt->n_result_files[k - opt->min_k]
						= read_cmdline_strings(argc,
						argv, i + 1, &opt->result_files[
						k - opt->min_k], (void *)opt);
					i += opt->n_result_files[k - opt->min_k];//
				} else {
					opt->K = read_uint(argc, argv, ++i,
								(void *)opt);
					if (errno)
						goto CMDLINE_ERROR;
					debug_msg(MINIMAL, opt->quiet,
							"K = %u.\n", opt->K);
				}
				break;
			case 'l':
				opt->kmodes_algorithm = FASTQ_LLOYDS;
				debug_msg(QUIET, opt->quiet,
					"Using Lloyd's algorithm for fastq.\n");
				break;
			case 'e':
				opt->kmodes_algorithm = FASTQ_LLOYDS_EFFICIENT;
				debug_msg(QUIET, opt->quiet, "Using efficient "
					"Lloyd's algorithm for fastq.\n");
				break;
			case 'm':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->mfile = argv[++i];
				if (access(opt->mfile, F_OK) == -1) {
					opt->mfile = NULL;
					opt->target = read_cmdline_double(argc,
							argv, i, (void *)opt);
				} else if (i + 1 < argc && argv[i + 1][0] != '-')
					opt->mfile_out = argv[++i];
				break;
			case 'w':
				opt->kmodes_algorithm = KMODES_HARTIGAN_WONG;
				debug_msg(QUIET, opt->quiet,
					"Using Hartigan and Wong algorithm.\n");
				break;
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->datafile = argv[++i];
				if (i + 1 < argc && argv[i+1][0] != '-') {
					opt->data_outfile = argv[++i];
					debug_msg(QUIET, opt->quiet, "Will "
						"write data, after possible "
						"adjustments, to file = %s\n",
						opt->data_outfile);
				}
				debug_msg(QUIET, opt->quiet, "Data file = %s\n",
								opt->datafile);
				break;
			case 'o':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->soln_file = argv[++i];
				if (i + 1 < argc && argv[i+1][0] != '-')
					opt->ini_file = argv[++i];
				break;
			case 'i':
				if (i + 1 == argc || argv[i + 1][0] == '-')
					goto CMDLINE_ERROR;
				if (!strcmp(argv[i+1], "rnd"))
					opt->init_method = KMODES_INIT_RANDOM_SEEDS;
				else if (!strcmp(argv[i+1], "h97"))
					opt->init_method = KMODES_INIT_H97;
				else if (!strcmp(argv[i+1], "h97rnd"))
					opt->init_method = KMODES_INIT_H97_RANDOM;
				else if (!strcmp(argv[i+1], "hd17"))
					opt->init_method = KMODES_INIT_HD17;
				else if (!strcmp(argv[i+1], "clb09"))
					opt->init_method = KMODES_INIT_CLB09;
				else if (!strcmp(argv[i+1], "clb09rnd"))
					opt->init_method = KMODES_INIT_CLB09_RANDOM;
				else if (!strcmp(argv[i+1], "av07"))
					opt->init_method = KMODES_INIT_AV07;
				else if (!strcmp(argv[i+1], "av07grd"))
					opt->init_method = KMODES_INIT_AV07_GREEDY;
				else if (!strcmp(argv[i+1], "rndp"))
					opt->init_method = KMODES_INIT_RANDOM_FROM_PARTITION;
				else if (!strcmp(argv[i+1], "rnds"))
					opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
				
				/* assume seed indices are being provided */
				else if (access(argv[i+1], F_OK) == -1) {
					
					opt->init_method = KMODES_INIT_USER_SEEDS;
					
					opt->n_sd_idx = read_cmdline_uints(argc,
						argv, ++i, &opt->seed_idx,
								(void *)opt);
					if (errno || !opt->n_sd_idx)
						goto CMDLINE_ERROR;
					debug_msg(MINIMAL, opt->quiet, "Seed "
								"indices:");
					debug_call(MINIMAL, opt->quiet,
						fprint_uints(stderr,
							opt->seed_idx,
							opt->n_sd_idx, 1, 1));
					i += opt->n_sd_idx - 1;
					
					/* assume seeds are provided in a file */
				} else {
					opt->sfile = argv[i+1];
					//opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
				}
				debug_msg(QUIET, opt->quiet, "Using %s "
					"initialization.\n",
					kmodes_init_method(opt->init_method));
				++i;
				break;
			case 'n':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->n_init = read_uint(argc, argv, ++i,
							(void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				debug_msg(MINIMAL, opt->quiet,
					"Initializations = %u\n", opt->n_init);
				break;
			case '1':
				opt->subtract_one = 1;
				break;
			case 'p':
				if (i + 1 == argc || (err =
					process_arg_p(argc, argv, &i, j, opt)))
					goto CMDLINE_ERROR;
				break;
			case 'r':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				if (!strncmp(&argv[i][j], "run", 3)) {
					opt->n_inner_init = read_uint(argc,
						argv, ++i, (void *)opt);
					debug_msg(MINIMAL, opt->quiet,
						"Inner initializations; %u\n",
						opt->n_inner_init);
				} else {
					opt->seed = read_uint(argc, argv, ++i,
								(void *)opt);
					srand(opt->seed);
					debug_msg(MINIMAL, opt->quiet, "Seed: "
							"%lu\n", opt->seed);
				}
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 's':	/*-s <sn> <sp> <sc> <st1> <st2>*/
				/* no appropriate arguments */
				if (!strcmp(&argv[i][j], "shuffle")) {
					opt->shuffle = 1;
					debug_msg(MINIMAL, opt->quiet, "Data "
						"will be shuffled.\n");
					break;
				}
				if (i + 5 >= argc || argv[i + 1][0] == '-')
					goto CMDLINE_ERROR;
				opt->simulate = 1;
				opt->sim_n_observations = read_ulong(argc, argv,
							++i, (void *)opt);
				opt->sim_n_coordinates = read_ulong(argc, argv,
							++i, (void *)opt);
				unsigned int tst = read_uint(argc, argv, ++i,
								(void *)opt);
				opt->sim_between_t = read_cmdline_double(argc,
						 argv, ++i, (void *)opt);
				opt->sim_within_t = read_cmdline_double(argc,
						argv, ++i, (void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				if (tst > pow(2, 8*sizeof(data_t))) {
					mmessage(ERROR_MSG, INVALID_USER_INPUT,
						 "Cannot simulate data in more "
						 "than %u categories.\n",
						 (unsigned int) pow(2,
							8*sizeof(data_t)));
					goto CMDLINE_ERROR;
				} else
					opt->sim_n_categories = tst;
				debug_msg(MINIMAL, opt->quiet, "Simulation:\n"
					"\t%u observations\n\t%u coordinates\n"
					"\t%u categories\n"
					"\t%f between variance\n"
					"\t%f within variance\n",
					opt->sim_n_observations,
					opt->sim_n_coordinates,
					opt->sim_n_categories,
					opt->sim_between_t, opt->sim_within_t);
				
				break;
			case 't':
				opt->seconds = read_cmdline_double(argc, argv,
							++i, (void *)opt);
				debug_msg(QUIET, opt->quiet, "Running for "
					"%.0fs.\n", opt->seconds);
				break;
			case 'u':
				opt->update_modes = 1;
				debug_msg(MINIMAL, opt->quiet,
					"Update modes: on\n");
				break;
			case 'q':
#ifndef __KMODES_NO_QTRANS__
				if (!strncmp(&argv[i][j], "qt", 2))
					opt->use_qtran = 0;
				else
#endif
					opt->quiet = QUIET;
				break;
			case 'h':
				if (!strcmp(&argv[i][j], "h97")) {
					opt->kmodes_algorithm = KMODES_HUANG;
					debug_msg(QUIET, opt->quiet,
						"Using Huang's algorithm.\n");
				} else if (!strncmp(&argv[i][j], "hart",
						strlen("hart"))) {
					opt->use_hartigan = 1;
					debug_msg(QUIET, opt->quiet,
						"Use Hartigan's update.\n");
				} else {
					fprint_usage(stderr, argv[0], opt);
					free_options(opt);
					exit(EXIT_SUCCESS);
				}
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}
	
	if (opt->seconds)
		opt->n_init = 1;
	if (opt->K == 1) {
		opt->n_init = 1;
		opt->n_inner_init = 1;
	}
	if (opt->n_sd_idx && opt->n_sd_idx != opt->K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Incorrect usage"
				"of -k and -s.\n");
	else if (opt->n_sd_idx || opt->pfile) {
		if (opt->n_init > 1)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Resetting to"
				 " one initialization.\n");
		opt->n_init = 1;
		opt->n_inner_init = 1;
		if (opt->seconds > 0)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Resetting to"
				 " run 0 seconds.\n");
		opt->seconds = 0.;
		if (opt->shuffle)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Resetting to"
				 " not shuffle data.\n");
		opt->shuffle = 0;
	} else if (opt->result_files) {
		if (opt->K > 0)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "No sense "
				 "setting K if estimating K. Ignoring -k %u.\n",
				 opt->K);
		opt->K = 0;
		if (opt->n_init > 1)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Resetting to"
				 " zero initializations.\n");
		opt->n_init = 0;
		if (opt->seconds > 0)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Resetting to"
				 " run 0 seconds.\n");
		opt->seconds = 0.;
		
		if (opt->pfile)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Ignoring "
				 "partition file '%s'.\n", opt->pfile);
		opt->pfile = NULL;
		
		if (opt->sfile)
			mmessage(WARNING_MSG, INVALID_USER_INPUT, "Ignoring "
				 "seed file '%s'.\n", opt->sfile);
		opt->sfile = NULL;
	}
	
	if (opt->continue_run && access(opt->soln_file, F_OK) == -1)
		return mmessage(ERROR_MSG, FILE_NOT_FOUND, opt->soln_file);
	else if (opt->continue_run && opt->ini_file && access(opt->ini_file,
							F_OK) == -1)
		return mmessage(ERROR_MSG, FILE_NOT_FOUND, opt->ini_file);
	if (opt->init_method == KMODES_INIT_RANDOM_FROM_PARTITION
				&& opt->true_column == UINT_MAX)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "To use rndp "
				"initialization, you must know the true cluster "
				"assignments: see option -c.\n");
	
	if (!opt->shuffle && opt->n_init > 1)
		mmessage(WARNING_MSG, NO_ERROR, "Failing to shuffle the data "
			 "can produce input order-determined behavior.\n");
	
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
 * Process command-line arguments starting with p.
 *
 * @param argc	number of arguments
 * @param argv	command-line arguments
 * @param i	current argument
 * @param j	first character of argument after -
 * @param opt	options object
 */
int process_arg_p(int argc, char const **argv, int *i, int j, options *opt)
{
	int use_dirichlet = 0;
	
	if (!strcmp(&argv[*i][j], "pi")) {
		if (!strcmp(argv[*i+1], "dir")) {
			use_dirichlet = 1;
			++(*i);
		}
		if (opt->sim_K && (opt->sim_K != read_cmdline_doubles(argc,
			argv, *i + 1, &opt->sim_pi, (void *)opt) || errno))
			return INVALID_CMD_OPTION;
		else {
			opt->sim_K = read_cmdline_doubles(argc, argv,
					*i + 1, &opt->sim_pi, (void *)opt);
			if (errno)
				return INVALID_CMD_OPTION;
			debug_msg(MINIMAL, opt->quiet, "Simulation K = %u\n",
								opt->sim_K);
		}
		(*i) += opt->sim_K;
		if (use_dirichlet) {
			opt->sim_alpha = malloc(opt->sim_K
						* sizeof *opt->sim_alpha);
			if (!opt->sim_alpha)
				return(mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"options::sim_alpha"));
			memcpy(opt->sim_alpha, opt->sim_pi, opt->sim_K
						* sizeof *opt->sim_pi);
			debug_msg(MINIMAL, opt->quiet, "Simulation alpha:");
			debug_call(MINIMAL, opt->quiet, fprint_doubles(stderr,
					opt->sim_alpha, opt->sim_K, 2, 1));
		} else {
			double sum = 0;
			for (unsigned int j = 1; j < opt->sim_K; ++j)
				sum += opt->sim_pi[j];
			if (sum >= 1.)
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"-pi <pdbl1> ... <pdblK> arguments must"
								"sum to 1.0.\n");
			opt->sim_pi[0] = 1 - sum;
			debug_msg(MINIMAL, opt->quiet, "Simulation pi:");
			debug_call(MINIMAL, opt->quiet, fprint_doubles(stderr,
						opt->sim_pi, opt->sim_K, 2, 1));
		}
	} else if (access(argv[*i + 1], F_OK) == -1)
		opt->n_effective_coordinates = read_cmdline_double(argc, argv,
							++(*i), (void *)opt);
	else
		opt->pfile = argv[++(*i)];
	
	return NO_ERROR;
} /* process_arg_p */
