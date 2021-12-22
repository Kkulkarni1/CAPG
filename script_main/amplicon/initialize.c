/**
 * @file initialize.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 * 
 * Initialize ampliclust.
 *
 * TODO
 * - initialization cannot handle background_model
 * - initialization clusters haplotypes of length \par data.min_read_length;
 *   initializes remaining positions in a foolish way
 *
 * DONE
 * X initialize from set of seeds
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <string.h>
#include <stdlib.h>

#include "ampliclust.h"
#include "initialize.h"
#include "simulate.h"
#include "model.h"
#include "data.h"
#include "initialize_options.h"
#include "aecm.h"
#include "math.h"
#include "kmodes.h"
#include "util.h"
#include "io.h"
#include "fastq.h"

int initial_partition(model *mod, data *dat, initialize_options *opt);
double stash_initialization(model *mod, data *dat, initialize_options *opt, initializer *ini, double total);
int randomize_seeds(data *dat, initializer *ini, size_t K, unsigned int len);
void seeds_to_haplotypes(data *dat, initialize_options *opt, initializer *ini);
int partition_model_free(data *dat, initialize_options *opt, initializer *ini, int best);
int evaluate_ini(data *dat, initialize_options *opt, initializer *ini);
double iterate_expected_true_abundance(data *dat, double *abun_true, double *e_trans, hash *unit, unsigned int *H, unsigned int current_k, unsigned int K, unsigned int obs_abun, double true_abun);
int expected_true_abundance(initialize_options *opt, data *dat, double *e_trans, double *abun_true, unsigned int *H, size_t *idx, unsigned int count_i, unsigned int select, unsigned int K, unsigned int i, int iterate, double max);

/**
 * Create initializer object.
 *
 * @param ini	initializer object
 * @param dat	data object
 * @param opt	initialize_options object
 * @paraa mod	pointer to model object
 * @param tbd	read nucleotides to be determined
 * @return	error status
 */
int make_initializer(initializer **ini, data *dat, model *mod,
					initialize_options *opt, int tbd)
{
	//int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	initializer *in;
	*ini = malloc(sizeof **ini);

	if (*ini == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer object");
	in = *ini;

	in->synced = 0;
	in->cluster_id = NULL;
	in->best_cluster_id = NULL;

	in->seed_idx = calloc(opt->modo->K, sizeof *in->seed_idx);
	in->best_seed_idx = calloc(opt->modo->K, sizeof *in->best_seed_idx);

	if (!in->seed_idx || !in->best_seed_idx)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "initializer::seed_idx");

	in->criterion = malloc(opt->modo->K * sizeof *in->criterion);
	in->best_criterion = malloc(mod->n_mix * sizeof *in->best_criterion);

	if (!in->criterion || !in->best_criterion)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer::criterion");
	in->cluster_size = malloc(mod->n_mix * sizeof *in->cluster_size);

	in->best_cluster_size = malloc(mod->n_mix
					* sizeof *in->best_cluster_size);

	if (!in->cluster_size || !in->best_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer::cluster_size");

	in->best_total = INFINITY;
	in->optimal_total = INFINITY;
	in->n_inits = 0;

	/* allocate room for k-modes initialization */
	in->seeds = malloc(opt->modo->K * sizeof *in->seeds);
	in->best_modes = malloc(opt->modo->K * sizeof *in->best_modes);
	in->seed_lengths = calloc(opt->modo->K, sizeof *in->seed_lengths);

	if (!in->seeds || !in->best_modes || !in->seed_lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "initializer.seeds");

	/* [TODO] allocate seeds in one block: easier to free and realloc:
	 * data_t *dptr = malloc(dat->max_read_length * opt->modo->K * sizeof **in_seeds);
	 * for (size_t i = 1; i < in->K; ++i) {
	 *	in->seeds[i] = dptr;
	 *	dptr += dat->max_read_length;
	 * }
	 * //... use in->seeds ... //
	 * if (in->seeds) {
	 *	if (in->seeds[0])
	 *		 free(in->seeds[0])
	 *	free(in->seeds);
	 * }
	 * */
	for (size_t i = 0; i < opt->modo->K; i++) {
		in->seeds[i] = calloc(dat->max_read_length, sizeof **in->seeds);
		in->best_modes[i] = calloc(dat->max_read_length,		// [TODO] initializer needs to know length of observations in data
						sizeof **in->best_modes);
		if (in->seeds[i] == NULL || in->best_modes[i] == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"initializer.seeds");
	}

	/* allocate room for cluster assignments */
	in->cluster_id = calloc(dat->sample_size, sizeof *in->cluster_id);	// [TODO] initializer needs to know number of observations from data
	in->best_cluster_id = calloc(dat->sample_size,
		sizeof *in->best_cluster_id);

	if (!in->cluster_id || !in->best_cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer.cluster_id");

	in->abun_true = NULL;
	in->p = NULL;
	in->H = NULL;
	in->e_trans = NULL;
	in->uniq_seq_idx = NULL;
	in->uniq_seq_count = NULL;

        if (opt->use_error_profile && !mod->dada2_profile
		&& opt->initialization_method == INIT_IMDEBLUR
                && (err = read_dada2_error_profile(mod, dat, mod->mopt)))
                        return err;


	if (!dat->fdata->empty && !tbd)
		err = sync_initializer(in, dat, opt);

	return err;
} /* make_initializer */


/**
 * Finish initializer setup once data available.
 *
 * @param ini	initializer pointer
 * @param dat	data pointer
 * @param opt	pointer to initialize_options object
 * @return	error status
 */
int sync_initializer(initializer *ini, data *dat, initialize_options *opt)
{
	if (opt->initialization_method != INIT_DEBLUR
		&& opt->initialization_method != INIT_IMDEBLUR)
		return NO_ERROR;

	/* allocate room for index and count array of unique sequences */
	ini->uniq_seq_idx = calloc(dat->hash_length, sizeof *ini->uniq_seq_idx);
	ini->uniq_seq_count = calloc(dat->hash_length,
						sizeof *ini->uniq_seq_count);

	if (!ini->uniq_seq_idx || !ini->uniq_seq_count)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer.uniq_seq\n");

	/* hash table should be sorted in an order of decreasing */
	if (store_index(dat->seq_count, ini->uniq_seq_idx, dat->hash_length)
		|| store_count(dat->seq_count, ini->uniq_seq_count,
							dat->hash_length))
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"store_index() or store_count()\n");

	ini->synced = 1;

	return NO_ERROR;
} /* sync_initializer */


/**
 * Reallocate initializer struct for a different K or sample_size.  It is
 * assumed that initialize_options::K and data::sample_size contain the new numbers.
 *
 * @param ini	pointer to initialization object
 * @param dat	pointer to data object
 * @param opt	pointer to initialize_options object
 * @return	error status
 */
int realloc_initializer(initializer *ini, data *dat, initialize_options *opt)
{
	/* reallocate the room for a new K*/
	/* seed_idx */
	/* claculate the previous K */
	unsigned int pre_K = opt->modo->K;

	if(pre_K < 2)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"realloc.initializer\n");
	size_t *seed_idx = realloc(ini->seed_idx,
					opt->modo->K * sizeof *ini->seed_idx);
	size_t *best_seed_idx = realloc(ini->best_seed_idx,
					opt->modo->K * sizeof *ini->best_seed_idx);

	if (!seed_idx || !best_seed_idx )
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.seed_idx\n");

	ini->seed_idx = seed_idx;
	ini->best_seed_idx = best_seed_idx;

	/* criterion */
	double *criterion = realloc(ini->criterion,
					opt->modo->K * sizeof *ini->criterion);
	double *best_criterion = realloc(ini->best_criterion,
					opt->modo->K * sizeof *ini->best_criterion);

	if (!criterion || !best_criterion)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.criterion\n");

	ini->criterion = criterion;
	ini->best_criterion = best_criterion;

	/* seeds, best_modes,seed_lengths*/
	/* free previous points */
	/* [TODO] easier if you allocate ini->{seeds, best_modes, etc} in one block (see other TODOs) */

	if (ini->seeds || ini->best_modes) {
		for (unsigned int i = 0; i < pre_K; ++i){
			if (ini->seeds[i])
				free(ini->seeds[i]);
			if (ini->best_modes[i])
				free(ini->best_modes[i]);
		}
		free(ini->seeds);
		free(ini->best_modes);
	}
	ini->seeds = NULL;
	ini->best_modes = NULL;

	data_t **seeds = malloc(opt->modo->K * sizeof *ini->seeds);
	data_t **best_modes = malloc(opt->modo->K * sizeof *ini->best_modes);

	/*
	for (unsigned int i = 0; i < opt->modo->K - K_change; i++)
	{
		if (ini->seeds[i])
			free(ini->seeds[i]);
		if (ini->best_modes[i])
			free(ini->best_modes[i]);
	}

	data_t **seeds = realloc(ini->seeds, opt->modo->K * sizeof *ini->seeds);
	data_t **best_modes = realloc(ini->best_modes, opt->modo->K * sizeof *ini->best_modes);
	*/

	unsigned int *seed_lengths = realloc(ini->seed_lengths,
				opt->modo->K * sizeof *ini->seed_lengths);

	if (!seeds || !best_modes || !seed_lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.initializer.seeds\n");

	ini->seeds = seeds;
	ini->best_modes = best_modes;
	ini->seed_lengths = seed_lengths;

	for (unsigned int i = 0; i < opt->modo->K; ++i) {
		if (!dat->max_read_length)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds\n");
		ini->seeds[i] = calloc(dat->max_read_length, sizeof *ini->seeds);
		ini->best_modes[i] = calloc(dat->max_read_length, sizeof *ini->best_modes);
		if (ini->seeds[i] == NULL || ini->best_modes[i] == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds\n");
	}

	/* cluster_size */
	unsigned int *cluster_size = realloc(ini->cluster_size, opt->modo->K
						* sizeof *ini->cluster_size);
	unsigned int *best_cluster_size = realloc(ini->best_cluster_size,
				opt->modo->K * sizeof *ini->best_cluster_size);

	if (!cluster_size || !best_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.cluster_size\n");

	ini->cluster_size = cluster_size;
	ini->best_cluster_size = best_cluster_size;

	/* cluster_id */
	unsigned int *cluster_id = realloc(ini->cluster_id,
		dat->sample_size * sizeof *ini->cluster_id);
	unsigned int *best_cluster_id = realloc(ini->best_cluster_id,
		dat->sample_size * sizeof *ini->best_cluster_id);

	if (!cluster_id || !best_cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.cluster_id\n");

	ini->cluster_id = cluster_id;
	ini->best_cluster_id = best_cluster_id;

	/* index and count array of unique sequences */
	if (ini->uniq_seq_idx) {
		size_t *uniq_seq_idx = realloc(ini->uniq_seq_idx,
			dat->hash_length * sizeof *ini->uniq_seq_idx);

		unsigned int *uniq_seq_count = realloc(ini->uniq_seq_count,
			dat->hash_length * sizeof *ini->uniq_seq_count);

		if (!uniq_seq_idx || !uniq_seq_count)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"realloc.initializer.uniq_seq\n");

		ini->uniq_seq_idx = uniq_seq_idx;
		ini->uniq_seq_count = uniq_seq_count;
	}

	if (ini->abun_true)
		free(ini->abun_true);
	if (ini->p)
		free(ini->p);
	if (ini->H)
		free(ini->H);
	if (ini->e_trans)
		free(ini->e_trans);

	ini->abun_true = NULL;
	ini->p = NULL;
	ini->H = NULL;
	ini->e_trans = NULL;

	/* hash table should have been sorted in an order of decreasing */
	if (store_index(dat->seq_count, ini->uniq_seq_idx, dat->hash_length)
		|| store_count(dat->seq_count, ini->uniq_seq_count,
							dat->hash_length))
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"store_index() or store_count()\n");

	ini->best_total = INFINITY;
	ini->optimal_total = INFINITY;
	ini->n_inits = 0;

	return NO_ERROR;
}/* realloc_initializer */


/**
 * Initialize parameters.
 *
 * User determines initialization by -i command line parameter.  When this
 * function is complete, all model parameters, including \ref model.haplotypes,
 * \ref model.pi, \ref model.gamma, \ref model.delta, \ref model.lambda0,
 * and \ref model.lambda1 are set.
 *
 * Note, this same code is used to estimate simulation parameters from data to
 * simulate data.  Please see set_simulation_paramters() for more information.
 *
 * @param mod	model object pointer
 * @param dat	data object pointer
 * @param ini	initializer object pointer
 * @param opt	initialize_options object pointer
 * @return	error status
 */
int initialize(model *mod, data *dat, initializer *ini, initialize_options *opt)
{
	int fxn_debug = DEBUG_III;//DEBUG_I;//ABSOLUTE_SILENCE;//
	int err = NO_ERROR;

	debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "entering with "
		"initialization method %d\n", opt->initialization_method);

	if (opt->initialization_method == INIT_TRUE_VALUES
		|| opt->initialization_method == INIT_TRUE_PARTITION)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Cannot use "
			"initialize() to initialize with truth.\n");

	/* user has provided initialization file of some kind */
	/* [KSD: must reconsider if makes sense inside iterative version.] */
	if (opt->initialization_file) {

		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "processing "
			"initialization file '%s'\n", opt->initialization_file);

		/* read cluster, seed file, or haplotypes file */
		if ((err = read_initialization_file(opt->initialization_file,
			dat, opt, ini)))
			return err;

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Initialization "
			"method: %d\n", opt->initialization_method);

		if (fxn_debug == DEBUG_III && opt->initialization_method
			== INIT_SEEDS)
			for (unsigned int k = 0; k < opt->modo->K; ++k)
				debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"Seed %u: %u\n", k, ini->seed_idx[k]);

		/* to initialize model parameters, we need a partition of the
		 * data and haplotypes, so this computes \par
		 * initializer.best_cluster_id or \par initializer.best_modes
		 * based on initialization */
		if ((err = pull_bootstraps(dat, opt, ini)))
			return err;
		//pull_bootstraps(dat, opt, ini);

/*
	} else if (opt->initialization_method == INIT_IMDEBLUR) {

		err = impro_deblur(opt, dat, ini, mod);

		for (size_t k = 0; k < opt->modo->K; ++k)
			memcpy(ini->best_modes[k], ini->seeds[k],
				ini->seed_lengths[k] * sizeof **ini->best_modes);
		memcpy(ini->best_seed_idx, ini->seed_idx,
				opt->modo->K * sizeof *ini->best_seed_idx);
		memcpy(ini->best_cluster_size, ini->cluster_size, opt->modo->K
					* sizeof *ini->cluster_size);
		memcpy(ini->best_cluster_id, ini->cluster_id, dat->sample_size
					* sizeof *ini->cluster_id);
		//mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->ll);
*/

	/* otherwise do some kind of random initialization */
	} else if (opt->initialization_method == INIT_RANDOM_SEEDS
		|| opt->initialization_method == INIT_KMODES
		|| opt->initialization_method == INIT_HW_KMODES
		|| opt->initialization_method == INIT_RANDOM_SEEDS
		|| opt->initialization_method == INIT_DEBLUR
		|| opt->initialization_method == INIT_IMDEBLUR
		|| opt->initialization_method == INIT_CAO ) {

		ini->best_total = INFINITY;	/* if initialize_options.n_inner_init>1 */
		mod->init_ll = -INFINITY;	/* RND-EM */					// [TODO] model
		int null_cluster = 1;		/* repeat if null cluster */
		int exceed_iter = 0;		/* unconverged initializer */
		double avg_total = 0;
		double avg_ll = 0;

		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "Random "
			"initialization to set initializer.best_cluster_id and "
			"initializer.best_modes.\n");

		/* repeat \par initialize_options.n_inner_init times or as long as the
		 * initializer provides a null cluster (no members) */
		for (long j = 0; null_cluster || j < (int) opt->n_inner_init;
			++j) {
			unsigned int iter;

			err = NO_ERROR;

			/* set \par initializer.seeds and \par
			 * initializer.cluster_id */

			/* initialize haplotypes by randomly selecting K
			 * observations and assigning reads to closest
			 * observation by model-free distance based on
			 * quality scores */
			if (opt->initialization_method == INIT_RANDOM_SEEDS) {

				/* choose seeds: data.seeds & data.seed_idx */
				randomize_seeds(dat, ini, opt->modo->K, 0);

				/* partition based on chosen seeds */
				err = partition_model_free(dat, opt, ini,
								DEFAULT_VALUES);


			/* initialize haplotypes by running Cao's method */
			} else if (opt->initialization_method == INIT_CAO) {

				err = kmodes_init_clb09(dat->dmat,
					dat->sample_size, dat->min_read_length,
					opt->modo->K, 0, 0, 1, ini->seeds,
						ini->seed_idx, _SIZE_T_P);  

				if (!err) {

					for (size_t k = 0; k < opt->modo->K; ++k)
						ini->seed_lengths[k]
							= dat->min_read_length;

					err = partition_model_free(dat, opt,
							ini, DEFAULT_VALUES);
				}

			/* initialize haplotypes by running deblur method */
			} else if (opt->initialization_method == INIT_DEBLUR) {

				err = deblur(dat, ini->uniq_seq_idx,
					ini->uniq_seq_count, opt->error_rate,
					opt->modo->K, ini->seeds, ini->seed_idx, 
					ini->seed_lengths, ini->p,
					ini->abun_true, ini->H);

				if (!err)
					err = partition_model_free(dat, opt,
							ini, DEFAULT_VALUES);

			/* initialize haplotypes by running modified deblur method */
			} else if (opt->initialization_method == INIT_IMDEBLUR) {
				err = impro_deblur(opt, dat, ini, mod);				// [TODO] model

				if (!err)
					err = partition_model_free(dat, opt,
							ini, DEFAULT_VALUES);

			/* initialize haplotypes by running k-modes */
			} else if (opt->initialization_method == INIT_KMODES) {
				/* choose seeds: \par initializer.seeds */
				kmodes_init(dat->dmat, dat->sample_size,
					dat->min_read_length, opt->modo->K, 0,
					ini->seeds, (unsigned int *)ini->seed_idx,
					opt->kmodes_initialization_method,
					opt->kmodes_opt->weighted);       /* bugs here*/

				for (size_t k = 0; k < opt->modo->K; ++k)
					ini->seed_lengths[k] = dat->min_read_length;

				kmodes_huang(dat->dmat, dat->sample_size,
					dat->min_read_length, ini->seeds,
					opt->modo->K, ini->cluster_id, 
					ini->cluster_size, opt->n_kmodes_iter,
					ini->criterion, &err, &iter,
					opt->kmodes_opt); //opt->distance, 1);

			/* initialize haplotypes by running buggy k-modes */
			} else {	/* INIT_HW_KMODES */

				/* choose seeds: data.seeds */
				randomize_seeds(dat, ini, opt->modo->K,
							dat->min_read_length);

				kmodes_hw(dat->dmat, dat->sample_size,
					dat->min_read_length, ini->seeds,
					opt->modo->K, ini->cluster_id,
					ini->cluster_size, opt->n_kmodes_iter,
					ini->criterion, &err, &iter,
					opt->kmodes_opt);// 0,1,0); 
			}


			/* compute sum of Hamming distances criterion */
			null_cluster = 0;
			double total = 0.;
			for (size_t k = 0; k < opt->modo->K; ++k) {
				total += ini->criterion[k];
				if (ini->cluster_size[k] == 0) {
					MESSAGE(global_wp, WARNING_MSG, NO_ERROR,
						"kmodes() null cluster:  "
						"retrying...\n");
					null_cluster = 1;
					break;
				}
			}
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
							"total: %d \n", total);

			/* check for error in k-modes */
			if (opt->initialization_method == INIT_KMODES) {
				if (err && err != KMODES_NULL_CLUSTER_ERROR
					&& err != KMODES_EXCEED_ITER_WARNING)
					MESSAGE(global_wp, ERROR_MSG, err, "%s (%d)\n",
						kmodes_error(err), err);
				if (err == KMODES_EXCEED_ITER_WARNING)
					exceed_iter = err;
			} else if (err) {
				MESSAGE(global_wp, ERROR_MSG, err, "");
			}

			/* null cluster error: try again */
			if (null_cluster) {
				--j;
				continue;
			}

			avg_total += total;
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "total: "
				"%d; avg_total: %d\n", total, avg_total / j);

			if (opt->run_method == RND_EM)
				MESSAGE(global_wp, DEBUG_MSG, NO_ERROR,
					"Inner initialization %d (%.0f):", j, total);
			else if (opt->run_method == INI_EM)
				MESSAGE(global_wp, DEBUG_MSG, NO_ERROR,
					"Inner initialization %d: %.0f\n", j, total);

			/* store initialization if better in \par
			 * initializer.best_* (also \par model.init_ll for
			 * RND_EM) */
			avg_ll += stash_initialization(mod, dat, opt, ini, total);		// [TODO] model tracing ...
		}

		if (opt->run_method == INI_EM)
			MESSAGE(global_wp, DEBUG_MSG, NO_ERROR,
				"Average criterion in %zu attempts: %.0f\n",
				opt->n_inner_init,
				avg_total / opt->n_inner_init);
		else if (opt->run_method == RND_EM)
			MESSAGE(global_wp, DEBUG_MSG, NO_ERROR,
				"Average criterion, ll in %zu attempts: "
				"%.0f, %.3f\n", opt->n_inner_init,
				avg_total / opt->n_inner_init,
				avg_ll / opt->n_inner_init);

		/* record better initial criterion */
		if (mod->best_init_ll < mod->init_ll)						// [TODO] model
			mod->best_init_ll = mod->init_ll;					// [TODO] model
		if (ini->best_total < ini->optimal_total)
			ini->optimal_total = ini->best_total;
		if (exceed_iter)
			MESSAGE(global_wp, ERROR_MSG, err, "%s\n",
				kmodes_error(exceed_iter));

		free_kmodes();
		reset_k();

	}

	/* use best initial partition and best modes from \par
	 * initializer.best_* to initialize parameters;
	 * lower stringency needed for fast initialization of mlogit model
	 */
	double eps = opt->modo->beta_epsilon;
	opt->modo->beta_epsilon = opt->beta_epsilon;
	initialize_model_parameters(mod, dat, opt, ini, BEST_VALUES);			// [TODO] model
	param_update(mod, dat, opt->modo, DEFAULT_VALUES);				// [TODO] model
	opt->modo->beta_epsilon = eps;

	//mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->ll);


	return err;
} /* initialize */


/**
 * Read initialization information.  Initialization file can be a fasta file
 * containing haplotypes, a file containing \ref initialize_options.K seed indices, a file
 * containing a partition (\ref data.sample_size unsigned ints in {0, ...,
 * initialize_options.K}), or a file with initial values for all parameters [TODO].  After
 * the call, either \ref initializer.best_cluster_id or \ref
 * initializer.best_seed_idx or \ref initializer.best_modes, each at least of
 * length \ref data.min_read_length, contains the necessary initialization
 * information and \ref initialize_options.initialization_method* indicates what kind of
 * initialization has been set up.
 *
 * [TODO] parameter initialization file not yet implemented
 *
 * @param filename	name of file containing initialization information
 * @param dat		data object pointer
 * @param opt		initialize_options object pointer
 * @param ini		initializer object pointer
 * @return		error status
 */
int read_initialization_file(char const * const filename, data *dat,
	initialize_options *opt, initializer *ini)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_III;//
	int err = NO_ERROR;
	FILE *fp = fopen(filename, "r");

	debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "Opening file '%s'\n",
								filename);

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	char c = fgetc(fp);

	rewind(fp);

	debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "First char '%c'\n", c);

	/* not fasta file */
	if (c != '>') {
		/* try reading partition file */
		err = fread_uints(fp, ini->best_cluster_id,
			dat->sample_size);

		/* failed: maybe this is a seed index file */
		if (err == FILE_FORMAT_ERROR) {
			rewind(fp);
			err = fread_size_ts(fp, ini->seed_idx, opt->modo->K);	// [TODO] opt, by the way this will break if initializing for simulation and initialize_options::true_K != initialize_options::K!
			if (!err) {
				mmessage(INFO_MSG, NO_ERROR, "recovered from "
					"previous error\n");
				opt->initialization_method = INIT_SEEDS;	// [TODO] opt
			}
		} else {
			opt->initialization_method = INIT_PARTITION;		// [TODO] opt
		}
	
	/* fasta file containing possible haplotypes */
	} else {
		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "entering fasta "
								"read\n");

		/* read in fasta-formatted haplotypes */
		fastq_data *fqd = NULL;
		fastq_options fop = {.read_encoding = XY_ENCODING};
		if ((err = fread_fastq(fp, &fqd, &fop))) {
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "err=%d\n",
									err);
			fclose(fp);
			return err;
		}

		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "finished reading "
			"%u sequences in %s format\n", fqd->n_reads,
			fqd->file_type == FASTA_FILE ? "fasta" : "fastq");

		if (fqd->n_lengths || fqd->n_max_length
						!= dat->max_read_length) {
			fclose(fp);
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"haplotypes in '%s' must be same length.\n",
				filename);
		}

		/* exactly K haplotypes provided: copy them to best_modes */
		if (opt->modo->K == fqd->n_reads) {					// [TODO] opt
			for (unsigned int k = 0; k < opt->modo->K; ++k)		// [TODO] opt
				memcpy(ini->best_modes[k],
					&fqd->reads[k * fqd->n_max_length],
					fqd->n_max_length * sizeof *fqd->reads);

		/* too many haplotypes provided, choose K randomly */
		} else if (opt->modo->K < fqd->n_reads) {				// [TODO] opt
			/* [TODO] should not be random here */
			/* [KSD] Why not? */
			unsigned int k = 0, t = 0;
			while (k < opt->modo->K) {					// [TODO] opt
				double r = rand() / (RAND_MAX + 1.);
				if ((fqd->n_reads - t) * r < opt->modo->K - k)	// [TODO] opt
					memcpy(ini->best_modes[k++],
						&fqd->reads[fqd->n_max_length * t],	//
						fqd->n_max_length * sizeof *fqd->reads);//
				++t;
			}
		/* too few haplotypes, choose the rest randomly from reads */
		} else {
			unsigned int k = 0;

			/* copy all available haplotypes over */
			for (k = 0; k < fqd->n_reads; ++k)
				memcpy(ini->best_modes[k],
					&fqd->reads[fqd->n_max_length * k],
					fqd->n_max_length * sizeof *fqd->reads);

			/* choose remaining at random from reads */
			unsigned int t = k;
			while (k < opt->modo->K) {					// [TODO] opt
				double r = rand() / (RAND_MAX + 1.);
				if ((dat->sample_size - t) * r < opt->modo->K - k	// [TODO] opt
					&& dat->lengths[t] == fqd->n_max_length)
					memcpy(ini->best_modes[k++],
						dat->dmat[t], dat->lengths[t]
						* sizeof **dat->dmat);
				++t;
			}
		}
		free_fastq(fqd);

		if (fxn_debug >= DEBUG_III)
			fprint_alignment2(stderr, ini->best_modes, opt->modo->K,	// [TODO] opt
				dat->max_read_length);

		opt->initialization_method = INIT_HAPLOTYPES;			// [TODO] opt

		/* evaluate hash table */
		if (opt->evaluate_initialization				// [TODO] opt
			&& (err = evaluate_ini(dat, opt, ini)))			// [TODO] opt tracing...
				return err;

		if (opt->partition_file) {					// [TODO] opt
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"partition file: '%s'\n",
					opt->partition_file);			// [TODO] opt
			err = read_uints(opt->partition_file,			// [TODO] opt
				ini->best_cluster_id,
				dat->sample_size);
			if (!err)
				opt->initialization_method			// [TODO] opt
					= INIT_HAPLOTYPES_AND_PARTITION;
		}
	}

	fclose(fp);
	return err;
} /* read_initialization_file */


/**
 * Pull ourselves up by the bootstraps.  Before calling this function, either
 * \ref data.best_seed_idx or \ref initializer.best_cluster_id is set.  At the
 * end of this function, \ref initializer.best_cluster_id and \ref
 * initializer.best_modes should be set indicating a partition and an estimated
 * haplotype.  We could have directly set \ref model.haplotypes, but the other
 * initializer, kmodes, does not use the vectorized matrix form, so
 * initialize_model_parameters() sets \ref model.haplotypes from \ref
 * initializer.best_modes.  Bad planning?
 *
 * @param dat	data object pointer
 * @param opt	initialize_options object pointer
 * @param ini	initializer object pointer
 * @return	error status
 */
int pull_bootstraps(data *dat, initialize_options *opt, initializer *ini)
{
	int fxn_debug = DEBUG_III;//ABSOLUTE_SILENCE;//
	int err = NO_ERROR;

	/* need to produce haplotype estimates in \ref initializer.best_modes,
	 * each of length \ref data.max_read_length */
	if (opt->initialization_method == INIT_PARTITION) {			// [TODO] opt

		/* [TODO] code is rather sloppy in here */

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Initializing with "
								"partition.\n");

		/* find maximum character in each column */
		size_t max_nchar = 0;
		size_t *nchar = malloc(dat->max_read_length * sizeof *nchar);
		if (!nchar)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nchar");
		for (size_t j = 0; j < dat->max_read_length; ++j) {
			size_t max = 0;
			for (size_t i = 0; i < dat->sample_size; ++i) {
				if (j >= dat->lengths[i])
					break;
				if (dat->dmat[i][j] > max)
					max = dat->dmat[i][j];
			}
			nchar[j] = max + 1;
			if (max + 1 > max_nchar)
				max_nchar = max + 1;
		}

		/* allocate space to keep character counts */
		size_t *cnt_chars = malloc(opt->modo->K * max_nchar * sizeof *cnt_chars);	// [TODO] opt

		if (!cnt_chars)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "cnt_chars");

		/* find the mode of each cluster */
		for (size_t j = 0; j < dat->max_read_length; ++j) {
			/* initialize counts */
			for (size_t k = 0; k < opt->modo->K; ++k)				// [TODO] opt
				for (size_t l = 0; l < nchar[j]; ++l)
					cnt_chars[k * max_nchar + l] = 0;

			/* count */
			for (size_t i = 0; i < dat->sample_size; ++i) {
				if (j >= dat->lengths[i])
					break;
				cnt_chars[ini->best_cluster_id[i] * max_nchar
					+ dat->dmat[i][j]]++;
			}

			/* find max: assumes all positions observed [TODO] */
			for (size_t k = 0; k < opt->modo->K; ++k) {				// [TODO] opt
				size_t max = 0;
				for (size_t l = 0; l < nchar[j]; ++l)
					if (cnt_chars[l] > max) {
						max = cnt_chars[k * max_nchar + l];
						ini->best_modes[k][j] = l;
					}
			}
		}
		for (size_t k = 0; k < opt->modo->K; ++k)					// [TODO] opt
			ini->seed_lengths[k] = dat->max_read_length;

		free(nchar);
		free(cnt_chars);

	/* need to sets \ref initializer.best_cluster_id from \par
	 * initializer.seed_idx */
	} else if (opt->initialization_method == INIT_SEEDS				// [TODO] opt
		|| opt->initialization_method == INIT_RANDOM_SEEDS) {			// [TODO] opt

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Initializing with "
								"seeds.\n");

		seeds_to_haplotypes(dat, opt, ini);					// [TODO] opt tracing ... done
		err = partition_model_free(dat, opt, ini, BEST_VALUES);			// [TODO] opt tracing ... done

	/* need to produce \ref initializer.best_cluster_id from \par
	 * initializer.best_modes */
	} else if (opt->initialization_method == INIT_HAPLOTYPES) {			// [TODO] opt
fprintf(stderr, "model: %u\n", opt->modo->parameterization);
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Using haplotypes "
			"to compute initializer.best_cluster_id using quality "
			"scores and uniform error model.\n");

		err = partition_model_free(dat, opt, ini, BEST_VALUES);			// [TODO] opt

	} else
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "No bootstraps to "
			"pull: initializer.best_modes and "
			"initializer.best_cluster_id specified by user.\n");

	return err;
} /* pull_bootstraps */


/**
 * initialize_func-compliant function to do initialization.  Use this one if 
 * there can be no option to initialize from truth, such as when you are
 * simulating.
 *
 * @param dat	data pointer
 * @param mod	model pointer
 * @param ini	initializer pointer
 * @param obj	obscured initialize_options pointer
 * @return	error status
 */
int simple_initialize(data *dat, model *mod, initializer *ini, void *obj)
{
	initialize_options *inio = (initialize_options *) obj;
	return initialize(mod, dat, ini, inio);
} /* simple_initialize */


/**
 * Convert seeds into haplotypes
 *
 * @param dat	data object pointer
 * @param opt	initialize_options object pointer
 * @param ini	initializer object pointer
 */
void seeds_to_haplotypes(data *dat, initialize_options *opt, initializer *ini)
{
	int fxn_debug = ABSOLUTE_SILENCE;

	/* set up the haplotypes in \par initializer.best_modes */
	for (unsigned int k = 0; k < opt->modo->K; ++k) {					// [TODO] opt
		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "Seed %u: %u", k,
			ini->seed_idx[k]);

		ini->seed_lengths[k] = dat->lengths[ini->seed_idx[k]];

		debug_msg_cont(DEBUG_III <= fxn_debug, fxn_debug, " %s\n",
			display_sequence(dat->dmat[ini->seed_idx[k]],
			ini->seed_lengths[k], XY_ENCODING));

		for (size_t j = 0; j < ini->seed_lengths[k]; ++j)
			ini->best_modes[k][j] = dat->dmat[ini->seed_idx[k]][j];

		/* DOES NOT WORK: SERIOUS BUG I STILL DON'T UNDERSTAND
		memcpy(ini->best_modes[k], dat->dmat[ini->seed_idx[k]],
			dat->min_read_length * **ini->seeds);*/

		debug_msg_cont(DEBUG_III <= fxn_debug, fxn_debug, "\t%s\n",
			display_sequence(ini->best_modes[k],
			ini->seed_lengths[k], XY_ENCODING));
	}

	if (fxn_debug >= DEBUG_I)
		fprint_alignment2(stderr, ini->best_modes, opt->modo->K,			// [TODO] opt
			dat->min_read_length);

} /* seeds_to_haplotypes */


/**
 * Partition data given true haplotypes.  Assumes quality scores are the
 * true probability of error and all substitutions are equally likely.
 *
 * @param dat	data object pointer
 * @param opt	initialize_options object pointer
 * @param ini	initializer object pointer
 * @param best	use \par initializer.best_modes and \par
 *		initializer.best_cluster_id
 * @return	error status
 */
int partition_model_free(data *dat, initialize_options *opt, initializer *ini,
								int best)
{
	int fxn_debug = DEBUG_III;//DEBUG_I;//ABSOLUTE_SILENCE;//DEBUG_II;//
	int err = NO_ERROR;
	unsigned int iter = 0;


	unsigned int *id = best ? ini->best_cluster_id : ini->cluster_id;
	data_t **haplotypes = best ? ini->best_modes : ini->seeds;
	unsigned int * size = best ? ini->best_cluster_size : ini->cluster_size;
	double *criterion = best ? ini->best_criterion : ini->criterion;


	if (opt->assignment_method == HAMMING) {					// [TODO] opt
		/* 0 iterations of k-modes without Huang97 mode updates
		 * assigns reads to clusters using Hamming distance */
		kmodes_options kopt = {0, 0, 0, 0};
		kmodes_huang(dat->dmat, dat->sample_size,
			dat->min_read_length, haplotypes,
			opt->modo->K, id, (unsigned int *)size, 0,	/* no iterations */	// [TODO] opt
			criterion, &err, &iter, &kopt);
		return err;
	}


	for (size_t k = 0; k < opt->modo->K; ++k) {						// [TODO] opt
		size[k] = 0;
		criterion[k] = 0.;
	}


	/* [WORKING] the Hamming distance version does not work the
	 * same as the call to kmodes_huang above: one difference
	 * is the lengths used here, but there must be more */
	for (size_t i = 0; i < dat->sample_size; ++i) {
		double mll = -INFINITY;
		double best_criterion = INFINITY;
		for (size_t k = 0; k < opt->modo->K; ++k) {					// [TODO] opt
			double ll = 0;
			double k_criterion = 0.;
			unsigned int len = MIN(ini->seed_lengths[k],
				dat->lengths[i]);

			/* no contribution from the read positions
			 * not represented by the haplotype */
			for (unsigned int j = 0; j < len; ++j) {
				double eprob = 0.;
				if (opt->assignment_method)				// [TODO] opt
					eprob = error_prob(dat->fdata,
						dat->qmat[i][j]);
				if (opt->assignment_method && haplotypes[k][j]		// [TODO] opt
					== dat->dmat[i][j])
					ll += log(1 - eprob);
				else if (opt->assignment_method) {			// [TODO] opt
					ll += log(eprob / 3);
					k_criterion += 1.;
				} else
					k_criterion += 1.;
			}
			if (opt->assignment_method && ll > mll) {			// [TODO] opt
				mll = ll;
				id[i] = k;
				best_criterion = k_criterion;
			} else if (k_criterion < best_criterion) {
				id[i] = k;
				best_criterion = k_criterion;
			}
			if (opt->assignment_method)					// [TODO] opt
				debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"Read %u, cluster %u: ll=%.3f hd=%.0f "
					"(ll=%.3f hd=%.0f k=%u)\n", i, k, ll,
					k_criterion, mll, best_criterion,
					id[i]);
			else
				debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"Read %u, cluster %u: hd=%.0f (hd=%.0f "
					"k=%u)\n", i, k, k_criterion,
					best_criterion, id[i]);
		}
		size[id[i]]++;
		if (opt->assignment_method)						// [TODO] opt
			criterion[id[i]] += mll;
		else
			criterion[id[i]] += best_criterion;

		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Assign read %u", i);//
		debug_msg_cont(DEBUG_III <= fxn_debug, fxn_debug, " (%s)",
			display_sequence(dat->dmat[i], dat->lengths[i],
								XY_ENCODING));
		debug_msg_cont(DEBUG_II <= fxn_debug, fxn_debug, " to cluster "
				"%u with score %f\n", id[i],
				opt->assignment_method ? mll : best_criterion);		// [TODO] opt
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Cluster sizes: ");
	debug_call(DEBUG_I <= fxn_debug, fxn_debug, fprint_uints(stderr, size,
								opt->modo->K, 4, 1));		// [TODO] opt	
	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Criteria: ");
	debug_call(DEBUG_I <= fxn_debug, fxn_debug, fprint_doubles(stderr,
			criterion, opt->modo->K, opt->assignment_method ? 3 : 0, 1));		// [TODO] opt


	return err;
} /* partition_model_free */


/**
 * Stash better initialization.  Judgement of better is determined by
 * initialize_options.run_method.  If using RND-EM, then the log likelihood is computed
 * on the intial parameters, computed from \par initializer.best_modes and
 * \par initializer.best_cluster_ids.  Otherwise, the k-modes criterion is used.
 * If a better solution is found, it is stored in \par initializer.best_seed_idx,
 * \par initializer.best_cluster_id, \par initializer.best_cluster_size, and
 * \par initializer.best_criterion.  \par initializer.best_total stores the
 * sum criterion and model.init_ll stores the best likelihood observed during
 * initialization.
 *
 * @param mod	model object pointer
 * @param dat	data object pointer
 * @param opt	initialize_options object pointer
 * @param ini	initializer object pointer
 * @param total	current initialization criterion (Hamming distances sum now)
 */
double stash_initialization(model *mod, data *dat, initialize_options *opt,
	initializer *ini, double total)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int better = 0;
	double ll = -INFINITY;

	/* RND_EM computes log likelihood to judge quality of initialization */
	if (opt->run_method == RND_EM) {
		/* initialize model parameters with \par data.cluster_id and
		 * \par data.seeds as set by initializer */
		initialize_model_parameters(mod, dat, opt, ini, DEFAULT_VALUES);	// [TODO] model
		param_update(mod, dat, opt->modo, DEFAULT_VALUES);			// [TODO] model

		/* calculate ll under quality model ? */
		//opt->literal_qualities = 1;
		ll = e_step(dat, mod, opt->modo, FIRST);					// [TODO] model tracing...
		//opt->literal_qualities = 0;

		MESSAGE_CONT(global_wp, " %f\n", ll);
		better = ll > mod->init_ll;						// [TODO] model
	/* other methods use the sum of Hamming distances criterion */
	} else {
		better = total < ini->best_total;
	}

	/* better solution to be stored in \par initializer.best_* variables */
	if (better) {
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Storing improved "
			"initialization with total %.0f\n", total);
		ini->best_total = total;
		mod->init_ll = ll;							// [TODO] model
		for (size_t k = 0; k < opt->modo->K; ++k)
			memcpy(ini->best_modes[k], ini->seeds[k],
				ini->seed_lengths[k] * sizeof **ini->best_modes);
		if (fxn_debug >= DEBUG_I)
			fprint_alignment2(stderr, ini->best_modes, opt->modo->K,
				dat->min_read_length);
		memcpy(ini->best_seed_idx, ini->seed_idx, opt->modo->K
			* sizeof *ini->best_seed_idx);
		memcpy(ini->best_cluster_id, ini->cluster_id,
			dat->sample_size * sizeof *ini->best_cluster_id);
		/* [TODO] This is recomputed in initialize_model_parameters. */
		memcpy(ini->best_cluster_size, ini->cluster_size, mod->n_mix
			* sizeof *ini->best_cluster_size);
		memcpy(ini->best_criterion, ini->criterion, mod->n_mix
			* sizeof *ini->best_criterion);
	}

	return opt->run_method == RND_EM ? ll : total;
} /* stash_initialization */


/**
 * Estimate model update parameters (model::n*) from given haplotypes and
 * partition of dataset.
 * 
 * After preparing a partition in \ref data.best_cluster_id or \ref
 * data.cluster_id and a modal haplotype in \ref data.best_modes or \ref
 * data.seeds [see initialize()], estimate model parameters by calling the
 * second CM-step to update all parameters but the haplotypes.
 *
 * It is very important to notice that you need to call \ref param_update(mod,
 * dat, opt, DEFAULT_VALUES) to actually initialize the default parameter
 * blocks.  Also, if you only want to estimate some parameters, then set
 * \ref model_initialize_options::param_estimate to indicate the ones to estimate.
 *
 * This function assumes that model::nhaplotypes is allocated, and as of now
 * it always is.
 *
 * @param mod	model object pointer
 * @param dat	data object pointer
 * @param opt 	initialize_options object pointer
 * @param ini	initialize object pointer
 * @param type	which partition to use (BEST_VALUES|DEFAULT_VALUES)
 * @return	error status
 */
int initialize_model_parameters(model *mod, data *dat, initialize_options *opt,
	initializer *ini, int type)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_II;//DEBUG_I;//
	unsigned int i, j, k;
	unsigned int *cid =		/* chosen clustering */
		type == BEST_VALUES ? ini->best_cluster_id : ini->cluster_id;
	unsigned char **hap = 		/* chosen haplotype */
		type == BEST_VALUES ? ini->best_modes : ini->seeds;

	debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Initializing from "
		"initializer.%scluster_id and initializer.%s.\n",
		type == BEST_VALUES ? "best_" : "",
		type == BEST_VALUES ? "best_modes" : "seeds");

	/* Initialize model.nhaplotypes since this is used in m_other(). */
	for (k = 0; k < opt->modo->K; ++k) {
		ini->best_cluster_size[k] = 0;
		for (j = 0; j < ini->seed_lengths[k]; ++j) {
//fprintf(stderr, "k = %u, j = %u\n", k, j);
			mod->nhaplotypes[k * dat->max_read_length + j] 			// [TODO] model
								= hap[k][j];
}
	}
	/* [TODO] initializer is unaware of background model */
	for (k = opt->modo->K; k < mod->n_mix; ++k)
		ini->best_cluster_size[k] = 0;

	debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_fasta(stderr,
		mod->nhaplotypes, opt->modo->K, dat->max_read_length, "S"));		// [TODO] model
	debug_call(fxn_debug >= DEBUG_I, fxn_debug,
		haplotype_delta_hamming_matrix(mod, dat, opt->modo, NEW_COPY,
								 NEW_COPY));

	/* [TODO] fix this silly kludge to fill last positions of haplotypes */
	for (k = 0; k < opt->modo->K; ++k) {
		for (j = ini->seed_lengths[k]; j < dat->max_read_length; ++j) {
			mod->nhaplotypes[k * dat->max_read_length + j] = XY_A;		// [TODO] model
			for (i = 0; i < dat->sample_size; ++i)
				if (j < dat->lengths[i])
					mod->nhaplotypes[cid[i]				// [TODO] model
						* dat->max_read_length + j]
						= dat->dmat[i][j];
		}
	}

	/* mini E-step from \ref initializer.best_cluster.id or \ref
	 * initializer.cluster_id; model.eik previously initialize to 0
	 */
	for (i = 0; i < dat->sample_size; ++i) {
		for (k = 0; k < opt->modo->K; ++k)
			mod->eik[i * mod->n_mix + k] = 0.;				// [TODO] model: move to estimator within model?  No, we need the whole model object anyway.
		mod->eik[i * mod->n_mix + cid[i]] = 1.;					// [TODO] model
		++ini->best_cluster_size[cid[i]];
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Assigning read %u "
					"to cluster %u (%u)\n", i, cid[i],
					 ini->best_cluster_size[cid[i]]);
	}

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "cluster sizes: ");
	debug_call(fxn_debug >= DEBUG_I, fxn_debug, fprint_uints(stderr,
		ini->best_cluster_size, opt->modo->K, 3, 1));

	/* this is the pseudo E-step for mlogit parameterization */
	if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION
		&& opt->modo->param_estimate & PARAM_BETA) {
		unsigned int jhq_index;

		zero_e(mod, dat);							// [TODO] model

		/* stupid initialization: nnet uses Unif(-r,r), where r is such
		 * that r*max(|x| \approx 1.
		 */
		for (i = 0; i < mod->n_beta_coef; ++i) {				// [TODO] model
			mod->beta[i] = 0;						// [TODO] model
                        if (!(i % mod->n_predictors))					// [TODO] model
                                mod->beta[i] = -1;					// [TODO] model
                        else if ((i % mod->n_predictors) == (i + mod->n_predictors) / mod->n_predictors)	// [TODO] model
                                mod->beta[i] = 1;					// [TODO] model
		}

		/* count expected number of h->b with q at position j;
		 *  model.e_bjhq and model.ejhq previous set to 0
		 */
		unsigned char *rptr = dat->fdata->reads;
		unsigned char *qptr = dat->fdata->quals;
		for (i = 0; i < dat->fdata->n_reads; ++i) {
			unsigned int len = read_length(dat->fdata, i);
			unsigned int off = dat->offset ? dat->offset[i] : 0;

			for (j = 0; j < len; ++j) {
				unsigned char h = mod->nhaplotypes[cid[i]		// [TODO] model
					* dat->max_read_length + j + off];
				unsigned char b = *(rptr + j);
				unsigned char q = *(qptr + j);
				jhq_index = (j + off) * mod->n_hq_combos		// [TODO] model
					+ h * mod->n_quality + q;			// [TODO] model
				mod->e_bjhq[b * mod->n_jhq_combos + jhq_index] += 1;	// [TODO] model
				mod->e_jhq[jhq_index] += 1;				// [TODO] model
			}
			rptr += len;
			qptr += len;
		}
		compute_transition_probabilities(mod, dat, DEFAULT_VALUES);		// [TODO] model tracing... No, I'm starting to conclude we need the whole model object, the whole data object, and our own initialize_options object in k-haplotypes.
	}

	m_other(dat, mod, opt->modo);							// [TODO] opt tracing ...

	return NO_ERROR;
} /* initialize_model_parameters */


/**
 * Select K random observations as seeds, setting \par initializer.seed_lengths
 * \par initializer.seeds, \par initializer.seeds_idx.
 *
 * @param dat	data object pointer
 * @param ini	initializer object pointer
 * @param K	number of seeds to select
 * @param len	fixed length (0 to use chosen read lengths)
 * @return	error status
 */
int randomize_seeds(data *dat, initializer *ini, size_t K, unsigned int len)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//DEBUG_II;//
	int same;

	if (dat->sample_size < K)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Requesting %u clusters with only %u reads.\n", K,
			dat->sample_size);

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Choosing seeds:");
	debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "\n");

	for (size_t i = 0; i < K; ++i) {
		do {
			same = 0;
			ini->seed_idx[i] = (size_t) ((double) rand()
				/ RAND_MAX * dat->sample_size);
			ini->seed_lengths[i] = !len
				? dat->lengths[ini->seed_idx[i]] : len;

			for (size_t j = 0; j < i; ++j)
				/* sample without replacement */
				/* check for same or equal seeds */
				if (ini->seed_idx[i] == ini->seed_idx[j]
					|| !amplicon_read_compare(dat,
					ini->seed_idx[i], ini->seed_idx[j])) {
					same = 1;
					break;
				}
		} while (same);
		/*debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "\n");
		debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "i:%i", (int)i);*/
		debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, " %4lu",
							ini->seed_idx[i]);
		memcpy(ini->seeds[i], dat->dmat[ini->seed_idx[i]],
			ini->seed_lengths[i] * sizeof **ini->seeds);
		debug_msg_cont(DEBUG_II <= fxn_debug, fxn_debug, "\t%s\n",
			display_sequence(ini->seeds[i],
			ini->seed_lengths[i], XY_ENCODING));
	}
	debug_msg_cont(DEBUG_I <= fxn_debug, fxn_debug, "\n");
	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "finish randomize_seeds\n");

	return NO_ERROR;
} /* randomize_seeds */


/**
 * Implement deblur method as a good initialization method (given K).
 * 
 * @param dat		data struct pointer
 * @param idx		index of uniq sequences 
 * @param count		abundances of uniq sequences
 * @param error_rate	error rate per nuc
 * @param K		num of sequences selected
 * @param seeds		seeds, to be set
 * @param seed_idx	index of seeds, to be set
 * @param seed_lengths	lengths of seeds, to be set
 * @param p		??
 * @param abun_true	??
 * @param H		haplotypes
 * 
 * @return		error status
 * 
 * */
int deblur(data *dat, size_t *idx, unsigned int *count, double error_rate,
	unsigned int K, data_t **seeds, size_t *seed_idx,
	unsigned int *seed_lengths, double *p, double *abun_true,
	unsigned int *H)
{

	int err = NO_ERROR;
	int fxn_debug = DEBUG_II;
	unsigned int cprime, d, ord, select;
	double alpha, cdf, r, sum=0.;
	double adj = 1.;

	if (!abun_true)
		abun_true = malloc(dat->hash_length * sizeof *abun_true);
	if (!abun_true)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"deblur.copy_count");

	if (!H)
		H = malloc(K * sizeof *H);
	if (!H)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "deblur.H");

	/* update the count table*/
	//if (store_count(dat->seq_count,count,dat->hash_length))
	//	return mmessage(ERROR_MSG,INTERNAL_ERROR,"deblur.uniqe_seq.store");

	/* initialize true abundance to observed counts */
	for (unsigned int i = 0; i < dat->hash_length; i++)
		abun_true[i] = count[i];

	/* give beta from the reference */
	/* error upper bound for 1 to 11 hamming distance and indel */
	double beta[13] = {0,0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001,
						0.001, 0.001, 0.0005, 0.01};

	/* caculate the alpha */
	double itr = 1;

	for (unsigned int j = 0; j < dat->min_read_length; j++)
		itr *= 1 - error_rate;
	alpha = 1 - itr;

	/* main algorithm from deblur */
	for (unsigned int i = 0; i < dat->hash_length; i++) {

		if (abun_true[i] <= 0)
			continue;

		/* adjusted num of the sequence */
		cprime = (unsigned int) abun_true[i] / (1 - alpha);

		unsigned int next = i + 1;

		for (unsigned int j = next; j < dat->hash_length; j++) {
			d = hamming_char_dis((char *)dat->dmat[idx[i]],
					(char *) dat->dmat[idx[j]],
					dat->max_read_length);

			if (d > 11) /* TODO: think about how to deal with indels */
				continue;
			/*
			if (count[j] - cprime * beta[d] < 0)
				count[j] = 0;
			else
				count[j] = count[j] - cprime * beta[d];
			*/
			/* could be lower than 0 */
			abun_true[j] = abun_true[j] - cprime * beta[d] * adj;
		}
	}

	/* check out bugs */
	for (unsigned int i = 0; i < dat->hash_length; i++)
		if (abun_true[i] > -10)
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"count%i:%i\n", i, abun_true[i]);

	/* copy K sequences of results to the seeds for a good initialization */
	/* select sequences from high to low observed abundance */
	select = 0;
	for (unsigned int i = 0; i < dat->hash_length; i++) {

		if (abun_true[i] > 0) {
			H[select] = i;
			seed_idx[select] = idx[i];
			seed_lengths[select] = dat->lengths[seed_idx[select]];

			/* TODO: bugs here. Have not checked the length of seeds */
			memcpy(seeds[select], dat->dmat[seed_idx[select]],
				seed_lengths[select] * sizeof ** seeds);
			select++;
		}
		if (select == K)
			break;
	}

	/* randomly choose sequences based on the count */
	if (select < K) {
		if (!p)
			p = malloc(dat->hash_length * sizeof *p);
		if (!p)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "deblur.p");

		for (unsigned int i = 0; i < dat->hash_length; i++) {
			if (abun_true[i] > 0)
				p[i] = 0.;
			else
				p[i] = exp(-((int)count[i] - abun_true[i])*3.
					/ ((int)count[i]*(int)count[i]));	/* [KSD] Why casting to int? */
				//p[i] = exp(count[i]/2.);
			sum += p[i];
		}

		/* Let the sum = 1 */
		for (unsigned int i = 0; i < dat->hash_length; i++)
			p[i] /= sum;

		for (unsigned int i = 0; i < dat->hash_length; i++)
			if (p[i] > 0.001)
				debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
							"p%i:%f\n", i, p[i]);

		sum = 1;
		while (select < K) {
			r = (double) rand() /RAND_MAX;
			r = r * sum; // scale r

			ord = 0;
			cdf = p[ord];
			while(r > cdf){
				cdf += p[++ord];
				if (ord == (dat->hash_length-1))
					break;
			}

			/* TODO: bugs here. Have not checked the length of seeds */

			/* check the identity */
			/*
			unsigned int same = 0;
			for (unsigned int h = 0; h < select; h++) {
				if (idx[ord ]== seed_idx[h]){
					same = 1;
					break;
				}
			}
			if (same == 1)
				continue;
			*/
			/* find a better way */
			sum = sum-p[ord];
			p[ord] = 0.;
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug,
					"seed_idx:%i:%i\n", ord, idx[ord]);
			H[select] = ord;
			seed_idx[select]=idx[ord];
			seed_lengths[select]=dat->lengths[idx[ord]];
			memcpy(seeds[select], dat->dmat[idx[ord]],
				dat->max_read_length * sizeof **seeds);
			select++;
		}

	}

	for(unsigned int k = 0; k < K; k++) {
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "True abundance of "
			"the %ith haplotype (%ith unique seq id): %f\n", k,
			H[k], abun_true[H[k]]);
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "observed "
			"abundance of the %ith haplotype (%i th unique seq id):"
			" %i\n", k, H[k], count[H[k]]);
	}


	/* Randomly choose rest of seeds */
	/*
	while (select < K) {
		size_t s = (size_t) ((double) rand()
				/ RAND_MAX * dat->sample_size);

		// check the length
		if (dat->lengths[s] != dat->max_read_length)
			continue;

		// check the identity of haplotypes 
		unsigned int same = 0;
		for (unsigned int h = 0; h < select; h++)
			if (hamming_char_dis((char *)seeds[h],(char *)dat->dmat[s],
				dat->max_read_length) == 0)
						same = 1;
		if (same == 1)
			continue;

		seed_idx[select] = s;
		seed_lengths[select] = dat->lengths[s];
		memcpy(seeds[select], dat->dmat[s],
			dat->max_read_length * sizeof **seeds);
		select++;
	}*/

	return err;
}/* deblur */

/**
 * An improved version of deblur when given K. Uses similar idea. 
 * 
 * @param dat	data struct pointer
 * @param opt	initialize_options struct pointer
 * @param ini	initializer struct pointer
 * 
 * @return	error status
 **/
int impro_deblur(initialize_options *opt, data *dat, initializer *ini, model *mod)
{

	int err = NO_ERROR;
	int fxn_debug = DEBUG_I;	//DEBUG_II;	//

	unsigned int select = 0;
	unsigned int n_candidate = dat->hash_length;
	double r, cdf;
	int false_positive = 0;
	double pre_bic = INFINITY;
	double pre_aic = INFINITY;
	double up_bound = INFINITY;	/* allow skip of false positives */
	double low_bound = opt->low_bound;		/* haplotypes have abundance > this */
	//unsigned int *array_fp = NULL;
	double sum;

	//double epsilon = 1e-8;

	/* eliminate singletons: hash sorted by decreasing abundance */
	for (unsigned int i = 0; i < dat->hash_length; i++)
		if (ini->uniq_seq_count[i] < low_bound) {
			n_candidate = i;
			break;
		}

	if (!n_candidate)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "All reads "
						"unique: too much error!\n");

	/* malloc arrays to store true abundance */
	if (!ini->abun_true)
		ini->abun_true = malloc(n_candidate * sizeof *ini->abun_true);
	if (!ini->abun_true)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::abun_true");

	/* store probabilities for stochastic approach */
	if (opt->stochastic_deplur) {
		if (!ini->p)
			ini->p = malloc(n_candidate * sizeof *ini->p);
		if (!ini->p)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"initializer::p");
	}

	/* malloc arrays to store the idx of haplotypes */
	if (!ini->H)
		ini->H = malloc(opt->modo->K * sizeof *ini->H);
	if (!ini->H)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "initializer::H");

	/* malloc e matrix (log value) */
	if (!ini->e_trans)
		ini->e_trans = malloc(dat->sample_size * (opt->modo->K + 1)
						* sizeof *ini->e_trans);
	if (!ini->e_trans)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::e_trans");

	/* initialize estimated abundance with observed abundance */
	for (unsigned int i = 0; i < n_candidate; i++) {
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Unique sequence "
			"%u abundance: %u\n", i, ini->uniq_seq_count[i]);
		ini->abun_true[i] = ini->uniq_seq_count[i];
		if (opt->stochastic_deplur)
			ini->p[i] = 0.;
	}

	/* choose the most abundant unique sequence as the first haplotype */
	unsigned int ord = 0;
	select = 0;

	ini->seed_idx[select] = ini->uniq_seq_idx[ord];  // idx in dmat and qmat
	ini->seed_lengths[select] = dat->lengths[ini->uniq_seq_idx[ord]];
	memcpy(ini->seeds[select], dat->dmat[ini->uniq_seq_idx[ord]],
		dat->max_read_length * sizeof **ini->seeds);
	ini->H[select] = ord ; // idx in unique sequence table

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Selecting %d with estimated "
		"true abundance %.3f\n", ord, ini->abun_true[ord]);

	/* calculate the expected correct reads (last row of the matrix ) */
	cal_e_trans(dat, opt, ini->e_trans, opt->modo->K,
				mod->dada2_profile, NULL, 1);

	/* calculate expected misreads from H to all reads */
	cal_e_trans(dat, opt, ini->e_trans, select, mod->dada2_profile,
						ini->seeds[select], 0);

	select++;

	for (unsigned int k = 1; k < opt->modo->K; k++) {

		double max = 0.;    // maximum estimated true abundance
		sum = 0.;

		/* deterministic version chooses next haplotype to be one with
		 * maximum estimated true abundance as long as it exceeds lower
		 * bound.
		 */
		for (unsigned int i = 0; i < n_candidate; i++) {

			/* [TODO] should find a more efficient way */
			if (find_uint(ini->H, i, select))
				continue;

			/* update the true abundance i */
			if (ini->abun_true[i] > low_bound
						|| opt->stochastic_deplur) {

				/* check for bugs later */
				/*
				if ((err = expected_true_abundance(opt, dat,
					ini->e_trans, ini->abun_true, ini->H,
					ini->uniq_seq_idx,
					ini->uniq_seq_count[i], select, opt->modo->K,
					i, opt->convergence_deplur, max)))
					mmessage(WARNING_MSG, INTERNAL_ERROR,
						"warning about no convergence generated when updating true abundance of %i in %ith step \n", i, k + 1);
				*/
				if((e_TrueAbun(opt, dat, ini->e_trans,
					ini->abun_true, ini->H,
					ini->uniq_seq_idx,
					ini->uniq_seq_count[i], select,
					opt->modo->K, i,
					opt->convergence_deplur)))
					mmessage(WARNING_MSG, INTERNAL_ERROR,
						"warning about no convergence generated when updating true abundance of %i in %ith step \n", i, k + 1);
				
				debug_msg(DEBUG_II <= fxn_debug, fxn_debug,
					"Estimated abundance of %d: %.3f\n", i,
					ini->abun_true[i]);
				if (ini->abun_true[i] > max
					&& ini->abun_true[i] <= up_bound) {
					ord = i;
					max = ini->abun_true[i];
				}
				if (opt->stochastic_deplur && ini->abun_true[i] > 1)
					sum += ini->abun_true[i];
			}
		}

		/* stochastic version chooses next haplotype in proportion to
		 * expected true abundance
		 */
		if (opt->stochastic_deplur) {

			/* calculate prob for each uniq seq */
			for (unsigned int i = 0 ;i < n_candidate; i++)
				if (!find_uint(ini->H, i, select)) {
					if (ini->abun_true[i] > 1)
						ini->p[i] = ini->abun_true[i] / sum;
					else
						ini->p[i] = 0.;
				}


			/* randomly choose a seed based on the count */
			if (sum > 0.) {

				r = (double) rand() / RAND_MAX;

				ord = 0;
				cdf = ini->p[ord];
				while (r > cdf) {
					cdf += ini->p[++ord];
					if (ord == n_candidate - 1)
						break;
				}
			}
		}

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Selecting %d with estimated true"
				" abundance %.3f\n", ord, ini->abun_true[ord]);
		/* update seed */
		ini->seed_idx[select] = ini->uniq_seq_idx[ord];  // idx in dmat and qmat
		ini->seed_lengths[select] = dat->lengths[ini->uniq_seq_idx[ord]];
		memcpy(ini->seeds[select], dat->dmat[ini->uniq_seq_idx[ord]],
			dat->max_read_length * sizeof **ini->seeds);
		ini->H[select] = ord;  //idx in unique sequence table

		/* calculate the transition probability from new haplotype */
		cal_e_trans(dat, opt, ini->e_trans, select, mod->dada2_profile,
							ini->seeds[select], 0);

		/* check for false positive 
		 * only consider the deterministic version here
		 */
		if (!opt->stochastic_deplur && opt->check_false_positive) {

			/* No further searching if we have reached the lower bound */
			if (max > low_bound)
				false_positive = 1; 

			do {
				
				/* calculate bic or other statistic to check 
				 * if there is support for new haplotype */
				false_positive = detect_false_positive(dat, mod,
					ini, &pre_bic, &pre_aic, select,
					opt->modo->use_aic,
						opt->use_error_profile, 0);

				/* Not a false positive */
				if (!false_positive) {
fprintf(stderr, "Accepting %u (%u)\n", ord, select);
					break;
				}
fprintf(stderr, "Rejecting %u (%u)\n", ord, select);

				up_bound = max;  // set the upper bound for further searching 

				max = 0.;

				/* search for next most abundant sequence */
				for (unsigned int i = 0; i < n_candidate; i++) {
					if (ini->abun_true[i] > max &&
						ini->abun_true[i] < up_bound) {
						ord = i;
						max = ini->abun_true[i];
					}
				}

				/* reach low bound, no more haplotypes */
				if (max < low_bound)
					break;

				debug_msg(DEBUG_I, fxn_debug, "Selecting %d "
					"with estimated true abundance %f\n",
					ord, max);
				/* update seed */
				ini->seed_idx[select] = ini->uniq_seq_idx[ord];  // idx in dmat and qmat
				ini->seed_lengths[select]
					= dat->lengths[ini->uniq_seq_idx[ord]];
				memcpy(ini->seeds[select],
					dat->dmat[ini->uniq_seq_idx[ord]],
					dat->max_read_length * sizeof **ini->seeds);
				ini->H[select] = ord;  //idx in unique sequence table

				/* calculate the transition probability */
				cal_e_trans(dat, opt, ini->e_trans, select,
					mod->dada2_profile, ini->seeds[select], 0);


			} while (false_positive);
		}

		/* reach low bound, no more haplotypes */
		if (max < low_bound)
			break;

		/* do not choose the same haplotype again */
		if (opt->stochastic_deplur)
			ini->p[ord] = 0.;

		/* to exclude detected false positive */
		up_bound = max; // set to true abundance of last haplotype
		select ++;
	}

	if (!opt->estimate_K) {
		opt->modo->K = select;
	} else if (select < opt->modo->K) {
		do {
			sum = 0;
			for (unsigned int i = 0 ;i < n_candidate; i++)
				if (!find_uint(ini->H, i, select))
					sum += ini->abun_true[i];
			for (unsigned int i = 0 ;i < n_candidate; i++)
				if (!find_uint(ini->H, i, select))
					ini->p[i] = ini->abun_true[i] / sum;
				else
					ini->p[i] = 0;

			/* randomly choose a seed based on the count */
			r = (double) rand() / RAND_MAX;

			ord = 0;
			cdf = ini->p[ord];
			while (r > cdf) {
				cdf += ini->p[++ord];
				if (ord == n_candidate - 1)
					break;
			}
			debug_msg(DEBUG_I, fxn_debug, "Randomly selecting %d "
				"with estimated true abundance %.3f\n", ord,
							 ini->abun_true[ord]);
			ini->seed_idx[select] = ini->uniq_seq_idx[ord];  // idx in dmat and qmat
			ini->seed_lengths[select] = dat->lengths[ini->uniq_seq_idx[ord]];
			memcpy(ini->seeds[select], dat->dmat[ini->uniq_seq_idx[ord]],
				dat->max_read_length * sizeof **ini->seeds);
			ini->H[select] = ord;  //idx in unique sequence table
			select++;
		} while (select < opt->modo->K);
	}

	/* the following call sets ini->cluster_id, ini->cluster_size,
	 * and mod->ll using up-to-date haplotypes
	 */
	detect_false_positive(dat, mod, ini, &pre_bic, &pre_aic, select, 0, opt->use_error_profile, 1);

	return err;
}/* impro_deblur */

/**
 * Expected true abundace of unique sequence i
 * 
 * @param dat         data object
 * @param e_trans     log expected transition prob of reads from H (imcomplete)
 * @param abun_true    updated true abun of unique sequences
 * @param H            idx of hap in unique sequence table
 * @param count_i      observed abundance of the unique sequence i 
 * @param select       Num of hap has been selected 
 * @param i            idx of the unique sequence in unique sequence table 
 * @param conve        convergence 1 or not 0
 * 
 * 
 **/
int e_TrueAbun(initialize_options *opt, data *dat, double *e_trans, double *abun_true, unsigned int *H, size_t *idx,
	unsigned int count_i,unsigned int select, unsigned int K,unsigned int i, int conve){

	//double true_abun_i_init = abun_true[i];
	double delta;
	unsigned int iter = 0;
	//double adj = opt->error_adjust;
	double sum_pro_H;  // expected prob the observed sequence from H
	//double epsilon = 1e-6;  
	//unsigned int n_iter = 100;
	
	/* find the read idx in the sample set of this unique sequence */
	hash *unit = NULL;
	unsigned char *seq = dat->dmat[idx[i]];
	unsigned int length = dat->lengths[idx[i]];
	HASH_FIND(hh, dat->seq_count, seq, length * sizeof *seq, unit);

	if (!unit)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Do not find the sequence in hash table");

	/* update true abundance of unique sequence i */
	double true_abun_i = count_i;
	for (unsigned int r = 0; r < count_i; r++){
		/* prob rmi <- H */
		sum_pro_H = 0.;
		for (unsigned int k =0 ; k < select; k++)
			sum_pro_H += exp(e_trans[k*dat->sample_size+unit->idx_array[r]])*abun_true[H[k]];
		true_abun_i -= sum_pro_H / (sum_pro_H + exp(e_trans[K*dat->sample_size+unit->idx_array[r]])*abun_true[i]);
	}
	iter++;

	if(conve){
		do{
			delta = (abun_true[i] - true_abun_i)/abun_true[i];
			//delta = abun_true[i]-true_abun_i;
			/* updated abundance should be smaller than current abundance */
			
			if((delta < 0) || (true_abun_i < 0) )
				return mmessage(WARNING_MSG, INTERNAL_ERROR,"True abundnace increase or under 0\n");

			abun_true[i] = true_abun_i;

			if (delta < opt->modo->epsilon || true_abun_i < 1){
				return NO_ERROR;
			}

			true_abun_i = count_i;
			for (unsigned int r = 0; r < count_i; r++){
				/* prob rmi <- H */
				sum_pro_H = 0.;
				for (unsigned int k =0 ; k < select; k++)
					sum_pro_H += exp(e_trans[k*dat->sample_size+unit->idx_array[r]])*abun_true[H[k]];
				true_abun_i -= sum_pro_H / (sum_pro_H + exp(e_trans[K*dat->sample_size+unit->idx_array[r]])*abun_true[i]);
			}
			iter ++;

			if(iter > opt->modo->n_iter)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
					"exceed max interations,%i\n", iter);

		}while(1);
	}else
		abun_true[i] = true_abun_i;

	
	return NO_ERROR;
}/* e_TrueAbun */

/**
 * Compute expected true abundance of ith unique sequence.  Terminates
 * early if deterministic version without check for false positives is
 * activated and the iterations drop below the current maximum.
 * 
 * @param dat		data object
 * @param e_trans	log expected misreads of each unique sequence from each
 			haplotype
 * @param abun_true	estimated true abundances (updated at index i)
 * @param H		indices of haplotypes in unique sequence table
 * @param idx		indices of reads in data matrix
 * @param count_i	observed abundance of the unique sequence i 
 * @param select	number of current haplotype
 * @param K		total number of haplotypes
 * @param i		idx of the unique sequence
 * @param iterate	iterate to convergence
 * @param max		maximum expected true abundance so far
 * @return		error status
 **/
int expected_true_abundance(initialize_options *opt, data *dat, double *e_trans,
	double *abun_true, unsigned int *H, size_t *idx, unsigned int count_i,
	unsigned int select, unsigned int K, unsigned int i, int iterate,
	double max)
{
	UNUSED(max);
//	int can_abort = !opt->check_false_positive && !opt->stochastic_deplur;
	double delta;
	unsigned int iter = 0;

	/* find the read idx in the sample set of this unique sequence */
	hash *unit = NULL;
	unsigned char *seq = dat->dmat[idx[i]];
	unsigned int length = dat->lengths[idx[i]];

	HASH_FIND(hh, dat->seq_count, seq, length * sizeof *seq, unit);

	if (!unit)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"Do not find the sequence in hash table");

/* we may not use this strategy here */
/*
	if (can_abort && abun_true[i] < max)
		return NO_ERROR;
*/
	/* abun_true[i] is initialized to count_i, but gets decremented for
	 * each haplotype added.  Since additional haplotypes will only
	 * decrease the expected abundance, the relative change should be
	 * negative even on the first iteration.
	 */
	if (iterate) {
		do {
			/* one iteration: abun_true[i] -> true_abun_i */
			double true_abun_i = iterate_expected_true_abundance(
					dat, abun_true, e_trans, unit, H,
					select, K, count_i, abun_true[i]);
			/*
			if (can_abort && true_abun_i < max)
				return NO_ERROR;
			*/
			/* relative change */
			delta = (abun_true[i] - true_abun_i) / abun_true[i];

			/* iterates should decrease */
			if (delta < 0 || true_abun_i < 0)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
					"True abundance increase or under 0\n");

			/* end if convergence or estimate is below 1 */
			if (delta < opt->modo->epsilon || true_abun_i < 1)
				return NO_ERROR;

			abun_true[i] = true_abun_i;

			iter++;

			if (iter > opt->modo->n_iter)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
					"exceed max interations,%i\n", iter);

		} while (1);
	} else {
		abun_true[i] = iterate_expected_true_abundance(dat, abun_true,
			e_trans, unit, H, select, K, count_i, abun_true[i]);
	}


	return NO_ERROR;
}/* expected_true_abundance */

/**
 * One iterate of the fixed point iteration to compute expected abundance of
 * a unique sequence.
 *
 * @param dat		pointer to data object
 * @param abun_true	estimated true abundances, updating entry
 * @param e_trans	log probability of misread to each read from existing
 *			haplotypes
 * @param unit		unique sequence, containing indices of individual reads
 * @param H		indices of current haplotypes, length current_k - 1
 * @param current_k	this is next haplotype to find
 * @param K		total number of haplotypes
 * @param obs_abun	observed abundance of unique sequence
 * @param true_abun	estimated expected true abundance of unique sequence, x
 * @return		new estimated expected true abundance, f(x)
 */
double iterate_expected_true_abundance(data *dat, double *abun_true,
	double *e_trans, hash *unit, unsigned int *H, unsigned int current_k,
	unsigned int K, unsigned int obs_abun, double true_abun)
{
	double sum_pro_H, tmp;

	double true_abun_new = obs_abun;
	for (unsigned int r = 0; r < obs_abun; r++) {
		/* prob rmi <- H */
		sum_pro_H = 0.;
		for (unsigned int k = 0; k < current_k; k++) {
			tmp = exp(e_trans[k * dat->sample_size
					+ unit->idx_array[r]]);
			sum_pro_H += tmp * abun_true[H[k]];
		}
		tmp = exp(e_trans[K * dat->sample_size + unit->idx_array[r]]);
		true_abun_new -= sum_pro_H / (sum_pro_H + tmp * true_abun);
	}
	return true_abun_new;
} /* iterate_expected_true_abundance */


/**
 * Detect whether the chosen seed is a false positive.
 *
 * Use a simple model to assess whether there is support in the data for
 * the most recently chosen haplotype.  Here we compute the likelihood
 * assuming the quality scores are true or given some kind of pre-computed
 * error profile, the same model used to compute the expected misreads.
 *
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param bic		previously computed bic, updated to new bic
 * @param aic		previously computed aic, updated to new aic
 * @param select	size of current haplotype set
 * @param use_aic	use aic or -1 for do not compute ll or ic
 * 
 * @return		1/0 false positive or not
 * */
int detect_false_positive(data *dat, model *mod, initializer *ini,
	double *bic, double *aic, unsigned int select, int use_aic,
	int error_profile, int final)
{
		
	//int err = NO_ERROR;
	int false_positive = 1;
	double new_bic, pre_bic = *bic; 
	double new_aic, pre_aic = *aic;
	double JC_ll;
	double max, sum;
	double epsilon = 1e-8;
	unsigned int K = final ? select : (select + 1);	/* [TODO, BUG] it was select + 1, but this causes segfaults */
	unsigned int class;
	unsigned int n_param =  K * dat->max_read_length + K - 1 // haplotypes and pi 
		+ (error_profile ? dat->n_quality * NUM_NUCLEOTIDES_SQUARED: 0); 
   
	/* initialization */
	for (unsigned int k = 0; k < K; k++)
		ini->cluster_size[k] = 0;

	/* assign each read to haplotype */
	for (unsigned int i = 0; i < dat->sample_size; i++) {
		max = -INFINITY;
		class = 0;
		for (unsigned int k = 0; k < K; ++k) {
			if (max < ini->e_trans[k * dat->sample_size + i]) {
				max = ini->e_trans[k * dat->sample_size + i];
				class = k;
			}
		}
		ini->cluster_size[class]++;
		ini->cluster_id[i] = class;
	}

	//for (unsigned int k = 0; k < K; ++k)
	//	mmessage(INFO_MSG, NO_ERROR, "%i.\n",ini->cluster_size[k]);

	/*
	kmodes_huang(dat->dmat, dat->sample_size,
		dat->min_read_length, ini->seeds,
		K, ini->cluster_id, (unsigned int *) ini->cluster_size, 0,
		ini->criterion, &err, &iter, opt->kmodes_opt);	// BUGGY
	reset_k();
	*/

	if (use_aic < 0)
		return 0;

	/* estimate haplotypes and mixing proportions */
	for (unsigned int k = 0; k < K; ++k) {
		//mmessage(INFO_MSG, NO_ERROR, "%i.\n",ini->cluster_size[k]);
		for (unsigned int j = 0; j < ini->seed_lengths[k]; ++j)
			mod->haplotypes[k * dat->max_read_length + j]
				= ini->seeds[k][j];
		mod->pi[k] = (double)ini->cluster_size[k]/dat->sample_size;
		if(mod->pi[k] < epsilon)
			mod->pi[k] = epsilon;    // avoid bugs
	}

	/* calculate ll */
	mod->ll = 0.;
	for (unsigned int i = 0; i < dat->sample_size; i++) {
		max = -INFINITY;
		for (unsigned int k = 0; k < K; k++) {
			mod->eik[k * dat->sample_size + i] = log(mod->pi[k])
				+ ini->e_trans[k * dat->sample_size + i];
			if (max < mod->eik[k*dat->sample_size + i])
				max = mod->eik[k*dat->sample_size + i];
		}
		sum = 0.;
		for (unsigned int k = 0; k < K; ++k) {
			mod->eik[k * dat->sample_size + i]
				= exp(mod->eik[k * dat->sample_size + i] - max);
			sum += mod->eik[k * dat->sample_size + i];
		}

		mod->ll += log(sum) + max;

		for (unsigned k = 0; k < K; ++k)
					mod->eik[k * dat->sample_size + i] /= sum;
	}
	mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->ll);

	/* calculate modified aic and bic */
	modified_ic(mod->haplotypes, mod->est_ancestor, mod->distance, mod->ll,
				K, &JC_ll, &new_aic, &new_bic, n_param,
				dat->max_read_length, dat->sample_size);

	/* better bic means the chosen seed is not a false positive */
	if ((use_aic ? new_aic : new_bic) < (use_aic ? pre_aic : pre_bic)) {
		false_positive = 0;
		*bic = new_bic;
		*aic = new_aic;
	}
	mmessage(INFO_MSG, NO_ERROR, "pre_bic : %f. new_bic: %f\n", pre_bic, new_bic);
	mmessage(INFO_MSG, NO_ERROR, "pre_aic : %f. new_aic: %f\n", pre_aic, new_aic);

	return false_positive;
}/* detect_false_positive */


/**
 * Calculate expected number of misreads from selected haplotype.
 *
 * [TODO, BUG] This code uses data::n_quality for dimensions of dada2 error
 * profile, but the error profile may consider more qualities than the data.
 * Use model::n_quality.
 *
 * @param dat		pointer to data object
 * @param opt		pointer to initialize_options object
 * @param e_trans	expected counts
 * @param select	index of source haplotype
 * @param error_profile	optional error probabilities
 * @param seq		sequence we are generating counts of (NULL if self true)
 * @param self		is the sequence the same as haplotype
 * @return		error status
 */
int cal_e_trans(data *dat, initialize_options *opt, double *e_trans, unsigned int select,
			double *error_profile, unsigned char *seq, int self)
{

	unsigned int n = dat->sample_size;
	int remove_gamma = opt->remove_gamma;
	double adj = opt->error_adjust;
	double l1third = remove_gamma ? 0 : log(1./3);

	for (unsigned int r = 0; r < n; r++) {
		e_trans[n*select + r] = 0;
		for (unsigned int j = 0; j < dat->lengths[r]; j++) {
			if (opt->use_error_profile) {  //A C G T
				if (self)
					e_trans[n*select + r] +=
						dada2_error(error_profile,
						dat->n_quality, dat->dmat[r][j],
						dat->dmat[r][j], dat->qmat[r][j]);
				else
					e_trans[n*select + r] +=
						dada2_error(error_profile,
						dat->n_quality, seq[j], 
						dat->dmat[r][j], dat->qmat[r][j]);

			/* assume quality scores are valid */
			} else {
				double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);
				if (self)
					e_trans[n*select + r] += log(1 - ep);
				else if (dat->dmat[r][j] != seq[j])
					e_trans[n*select + r] += log(ep) + l1third;
				else
					e_trans[n*select + r] += log(1 - ep);
			}
		}
	}

	return NO_ERROR;
}/* cal_e_trans */


/**
 *  Compare true haplotypes with duplicated unique sequences table. 
 * 	Output the result for training the parameters.
 * 
 * [XP] May not be useful since there is a better way using updated hash table.
 */
int evaluate_ini(data *dat, initialize_options *opt, initializer *ini)
{
	UNUSED(opt);
	int err = NO_ERROR;
	int find = 0;
	unsigned int *d_matrix = NULL;
	int fxn_debug = DEBUG_I;

	FILE *fp = fopen("training_file.txt", "w");
	if (!fp)
		return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR,			// [TODO] opt: problem calling global_wp!  Solution: return error status and output message there
							"training_file");

	/* output the idx and abundance of unique sequences */
	fprintf(fp, "uniq seq idx: ");
	fprint_size_ts(fp, ini->uniq_seq_idx, dat->hash_length, 4, 1);
	fprintf(fp, "uniq seq abundance: ");
	fprint_uints(fp, ini->uniq_seq_count, dat->hash_length, 2, 1);

	/* compare the true seeds unique sequences */
	for (unsigned int i = 0; i < opt->modo->K; i++) {
		for (unsigned int j = 0; j < dat->hash_length; j++)
			if (hamming_char_dis((char *)ini->best_modes[i],
				(char *)dat->dmat[ini->uniq_seq_idx[j]],
				dat->max_read_length) == 0) {
				ini->best_seed_idx[i] = ini->uniq_seq_idx[j];
				find = 1;
				debug_msg(DEBUG_I <= fxn_debug, fxn_debug,
					"The %ith unique sequence is the %ith "
					"haplotype with abundance %i\n", j, i,
						 ini->uniq_seq_count[j]);
				break;
			}

		if (!find)
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "The true "
				"sequence %i does not exist in dataset\n", i);
	}

	/* output the idx of true seeds */
	fprintf(fp, "true seeds: ");
	fprint_size_ts(fp, ini->best_seed_idx, opt->modo->K, 4, 1);

	/* compare the sample dataset with unique sequences table */
	/* TODO: need a better way */
	for (unsigned int i = 0; i < dat->sample_size; i++)
		for (unsigned int j = 0; j < dat->hash_length; j++)
			if (hamming_char_dis((char *)dat->dmat[i],
				(char *)dat->dmat[ini->uniq_seq_idx[j]],
				dat->max_read_length) == 0) {
				ini->best_cluster_id[i]
					= (unsigned int) ini->uniq_seq_idx[j];
				break;
			}

	/* output the duplicated dataset */
	fprintf(fp, "duplicated: ");
	fprint_uints(fp, ini->best_cluster_id, dat->sample_size,4,1);

	/* distance matrix */
	d_matrix = malloc(dat->hash_length * dat->hash_length
							* sizeof *d_matrix);
	if (!d_matrix)
		return mmessage(ERROR_MSG,MEMORY_ALLOCATION,"evalu_ini");

	for (unsigned int i = 0; i < dat->hash_length; i++) {
		for (unsigned int j = i; j < dat->hash_length; j++) {
			unsigned int d = hamming_char_dis(
				(char *)dat->dmat[ini->uniq_seq_count[i]],
				(char *)dat->dmat[ini->uniq_seq_count[j]],
							dat->max_read_length);
			d_matrix[j + i * dat->hash_length] = d;
			d_matrix[i + j * dat->hash_length] = d;

		}
	}

	fprintf(fp, "distance matrix \n");
	fprint_vectorized_uintmatrix(fp, d_matrix, dat->hash_length,
						dat->hash_length, 1);

	fclose(fp);
	if (d_matrix)
		free(d_matrix);

	return err;
}/* evaluate_ini */


void free_initializer(initializer *ini, initialize_options *opt)
{
	if (ini) {
		if (ini->cluster_id)
			free(ini->cluster_id);
		if (ini->best_cluster_id)
			free(ini->best_cluster_id);

		if (ini->seed_idx)
			free(ini->seed_idx);
		if (ini->best_seed_idx)
			free(ini->best_seed_idx);

		if (ini->criterion)
			free(ini->criterion);
		if (ini->best_criterion)
			free(ini->best_criterion);

		if (ini->cluster_size)
			free(ini->cluster_size);
		if (ini->best_cluster_size)
			free(ini->best_cluster_size);

		if (ini->seeds) {
			for (size_t i = 0; i < opt->modo->K; ++i)
				if(ini->seeds[i])
					free(ini->seeds[i]);
			free(ini->seeds);
		}

		if (ini->seed_lengths)
			free(ini->seed_lengths);

		if (ini->best_modes) {
			for (size_t i = 0; i < opt->modo->K; ++i)
				if (ini->best_modes[i])
					free(ini->best_modes[i]);
			free(ini->best_modes);
		}
		if (ini->uniq_seq_count)
			free(ini->uniq_seq_count);
		if (ini->uniq_seq_idx)
			free(ini->uniq_seq_idx);
		if (ini->abun_true)	
			free(ini->abun_true);
		if (ini->p)
			free(ini->p);
		if (ini->H)
			free(ini->H);
		if(ini->e_trans)
			free(ini->e_trans);

		free(ini);
	}
} /* free_initializer */
