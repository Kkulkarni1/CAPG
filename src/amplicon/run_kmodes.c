/**
 * @file run_kmodes.c
 * @author Karin S. Dorman
 *
 */

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

//#define USE_CURSES
//#define __KMODES_NO_QTRANS__

#include "run_kmodes.h"
#include "cmdline.h"
#include "error.h"
#include "io.h"
#include "io_kmodes.h"
#include "cluster.h"
#include "timing_mach.h"
#include "matrix_exponential.h"

int make_options(options **opt);
int parse_options(options *opt, int argc, const char **argv);
int process_arg_p(int argc, char const **argv, int *i, int j, options *opt);
void free_options(options *opt);
int make_data(data **data, options *opt);
int allocate_data_for_k(data *dat, unsigned int k);
int finish_make_data(data *dat, options *opt);
int simulate_data(data *dat, options *opt);
int read_data(data *dat, options *opt);
int initialize_outfiles(data *dat, options * opt, FILE **in_fps, FILE **in_fpi);
int estimate_k(data *dat, options *opt);
void compute_costs(double *crit, double *var, double *size, unsigned int n, unsigned int K, double *cost, double *rrcost, double *krcost, double *sd, double *rsd, double *ksd);
void compute_jump_stats(double cost, double rrcost, double krcost, double pcost, double prrcost, double pkrcost, double *jump, double *rjump, double *kjump, double Y, int reset);
int restore_state(data *dat, options *opt);
int shuffle_data(data *dat, options *opt);
static inline int initialize(data *dat, options *opt);
void write_status(data *dat, options *opt, FILE *fps, FILE *fpi, unsigned int i, TIME_STRUCT *start);
void write_best_solution(data *dat, options *opt, FILE *fps);
static inline void stash_state(data *dat, options *opt);
static inline double hd(data_t *x, data_t *y, unsigned int p);
void free_data_for_k(data *dat, options *opt);
void free_data(data *dat);

int main(int argc, const char **argv)
{
	int err = NO_ERROR;		/* error code */
	options *opt = NULL;		/* run options */
	data *dat = NULL;		/* data object */
	TIME_STRUCT start;		/* timing */

	if ((err = make_options(&opt)))
		goto CLEAR_AND_EXIT;

	if ((err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;

	if ((err = make_data(&dat, opt)))
		goto CLEAR_AND_EXIT;

	if (opt->simulate) {
		if ((err = simulate_data(dat, opt)))
			goto CLEAR_AND_EXIT;
	} else if ((err = read_data(dat, opt)))
		goto CLEAR_AND_EXIT;

	if ((err = finish_make_data(dat, opt)))
		goto CLEAR_AND_EXIT;

	if (opt->estimate_k) {
		err = estimate_k(dat, opt);
		goto CLEAR_AND_EXIT;
	}

	if (opt->continue_run && (err = restore_state(dat, opt)))
		goto CLEAR_AND_EXIT;

	FILE *fpi = NULL, *fps = NULL;	/* output, solution */
	if (opt->soln_file && (err = initialize_outfiles(dat, opt, &fps, &fpi)))
		goto CLEAR_AND_EXIT;

	if (opt->n_init && (opt->kmodes_algorithm == KMODES_HUANG
			|| opt->kmodes_algorithm == KMODES_HARTIGAN_WONG)) {
		err = make_kmodes_options(&opt->kopt);
		if (err)
			goto CLEAR_AND_EXIT;
		opt->kopt->weighted = opt->weight;
		opt->kopt->init_update = opt->update_modes;
		opt->kopt->use_qtran = opt->use_qtran;
		opt->kopt->use_hartigan = opt->use_hartigan;
	}

	MARK_TIME(&start);
	for (unsigned int i = 0; i < opt->n_init; ++i) {

		/* shuffle data to avoid input order artifacts */
		if (opt->shuffle && opt->K > 1 && (err = shuffle_data(dat, opt)))
			goto CLEAR_AND_EXIT;

		/* initialization */
		if ((err = initialize(dat, opt))) {
			mmessage(ERROR_MSG, CUSTOM_ERROR, "%s\n",
				kmodes_error(err));
			goto CLEAR_AND_EXIT;
		}

		/* run one of k-modes algorithms */
		if (opt->kmodes_algorithm == KMODES_HUANG)
			dat->total = kmodes_huang(dat->dmat,
				dat->n_observations, dat->n_coordinates,
				dat->use_ini ? dat->ini_seeds : dat->seeds,
				opt->K, dat->cluster_id, dat->cluster_size,
				opt->n_max_iter, dat->criterion, &err,
				&dat->iter, opt->kopt);//opt->weight, opt->update_modes);
		else if (opt->kmodes_algorithm == KMODES_LLOYD)
			dat->total = kmodes_lloyd(dat->dmat,
				dat->use_ini ? dat->ini_seeds : dat->seeds,
				dat->cluster_size, dat->cluster_id,
				dat->n_observations, dat->n_coordinates, opt->K,
				opt->n_max_iter, dat->criterion, &err,
				&dat->iter, opt->weight);
		else if (opt->kmodes_algorithm == KMODES_HARTIGAN_WONG)
			dat->total = kmodes_hw(dat->dmat, dat->n_observations,
				dat->n_coordinates,
				dat->use_ini ? dat->ini_seeds : dat->seeds,
				opt->K, dat->cluster_id, dat->cluster_size,
				opt->n_max_iter, dat->criterion, &err,
				&dat->iter, opt->kopt); //opt->weight, opt->update_modes, opt->use_qtran);
		else {
			mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Invalid algorithm.\n");
			goto CLEAR_AND_EXIT;
		}

		/* repeat if obtain null cluster */
		if (err == KMODES_NULL_CLUSTER_ERROR) {
			if (i) --i;
			continue;
		} else if (err) {
			mmessage(ERROR_MSG, CUSTOM_ERROR, "%s\n",
				kmodes_error(err));
			if (err != KMODES_EXCEED_ITER_WARNING)
				goto CLEAR_AND_EXIT;
		}

		/* counts function call time ... oh well */
		write_status(dat, opt, fps, fpi, i, &start);

		/* record best solution */
		if ((!err || err == KMODES_EXCEED_ITER_WARNING)
			&& dat->total < dat->best_total) {
			stash_state(dat, opt);
		} else if (err)
			fprintf(stderr, "[ERROR] %s (%d)\n",
				kmodes_error(err), err);

		if (opt->seconds > 0 && ELAP_TIME(&start) > opt->seconds)
			break;
		else if (opt->seconds > 0) {
			i = 0;
			opt->n_init++;
		}
	}
	dat->seconds += ELAP_TIME(&start) - dat->uncounted_seconds;
	dat->n_init += opt->n_init;

	write_best_solution(dat, opt, fps);

CLEAR_AND_EXIT:

	free_kmodes();
	free_cluster_statics();

	if (dat)
		free_data(dat);
	if (opt)
		free_options(opt);

	return(err);
} /* main */


static inline
int initialize(data *dat, options *opt)
{
	int err = NO_ERROR;
	double inner_total = INFINITY, dtmp;
	data_t **seeds = dat->seeds;
	unsigned int *seed_idx = dat->seed_idx;

	dat->use_ini = 0;
	for (unsigned int j = 0; j < opt->n_inner_init; ++j) {
		if (!opt->n_sd_idx && !opt->pfile) {

			if (opt->init_method
				== KMODES_INIT_RANDOM_FROM_PARTITION)
				kmodes_init_random_from_partition(dat->dmat,
					dat->n_observations, dat->n_coordinates,
					opt->K, seeds, seed_idx,
					opt->true_cluster);
			else if (opt->init_method
				== KMODES_INIT_RANDOM_FROM_SET)
				kmodes_init_random_from_set(opt->K,
					dat->n_coordinates, opt->n_seed_set,
					seeds, opt->seed_set);
			else if (opt->shuffle && opt->init_method
				== KMODES_INIT_RANDOM_SEEDS) {
				unsigned int k = 0, l = 0, m;
				do {
						/*
					memcpy(seeds[k], dat->dmat[k],
						dat->n_coordinates * **seeds);
						*/
					seed_idx[k] = l;
					for (m = 0; m < dat->n_coordinates; ++m)
						seeds[k][m] = dat->dmat[l][m];

					for (m = 0; m < k; ++m)
						if (!hd(seeds[k], seeds[m], dat->n_coordinates))
							break;
					if (!k || m == k) ++k;
					++l;
				} while (k < opt->K);
			} else
				kmodes_init(dat->dmat, dat->n_observations,
					dat->n_coordinates, opt->K,
					opt->n_seed_set,
					seeds, seed_idx,
					opt->init_method, opt->weight);
		} else if (opt->n_sd_idx) {	/* will only happen once */
			for (unsigned int k = 0; k < opt->K; ++k) {
				for (unsigned int j = 0; j < dat->n_coordinates;
					++j)
					seeds[k][j] =
						dat->dmat[seed_idx[k]][j];
				/* retaining: WHY DOES THIS NOT WORK???
				memcpy(dat->seeds[k],
					dat->dmat[seed_idx[k]],
					dat->n_coordinates * *dat->seeds[k]);
				*/
			}
		} else if (opt->pfile) {	/* will only happen once */
			kmodes_init_from_partition(dat->dmat,
				dat->n_observations, dat->n_coordinates, opt->K,
				opt->weight, seeds, dat->cluster_id);
		}
		if (opt->n_inner_init > 1) {
			if (opt->kmodes_algorithm == KMODES_HUANG)
				dtmp = kmodes_huang(dat->dmat,
					dat->n_observations, dat->n_coordinates,
					seeds, opt->K, dat->cluster_id,
					dat->cluster_size, 0, dat->criterion,
					&err, &dat->iter, opt->kopt);//opt->weight, opt->update_modes);
			else if (opt->kmodes_algorithm == KMODES_HARTIGAN_WONG)
				dtmp = kmodes_hw(dat->dmat, dat->n_observations,
					dat->n_coordinates, seeds, opt->K,
					dat->cluster_id, dat->cluster_size, 0,
					dat->criterion, &err, &dat->iter,
					opt->kopt); //opt->weight, opt->update_modes, opt->use_qtran);
			else
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"Invalid algorithm in ini-kmodes.\n");

			if (opt->quiet > MINIMAL)
				fprintf(stdout, "Inner initialization %*u of "
					"%u: %f (%f)\n", (int)
					(log10(opt->n_inner_init) + 1), j,
					opt->n_inner_init, dtmp, inner_total);

			/* repeat if obtain null cluster */
			if (err == KMODES_NULL_CLUSTER_ERROR) {
				if (j) --j;
				continue;
			} else if (err)
				return err;

			if (dtmp < inner_total) {
				inner_total = dtmp;
				seeds = seeds == dat->seeds
					? dat->ini_seeds : dat->seeds;
				seed_idx = seed_idx == dat->seed_idx
					? dat->ini_seed_idx : dat->seed_idx;
				dat->use_ini = !dat->use_ini;
			}
		}
	}

	if (opt->n_inner_init > 1)
		dat->use_ini = !dat->use_ini;
	return err;
} /* initialize */

/**
 * Stash current state in best_* slot.
 *
 * @param dat	data object pointer
 * @param opt	options object pointer
 */
static inline void stash_state(data *dat, options *opt)
{
	dat->best_total = dat->total;
	for (unsigned int k = 0; k < opt->K; ++k)
		if (dat->use_ini)
			memcpy(dat->best_modes[k], dat->ini_seeds[k],
				dat->n_coordinates * sizeof **dat->best_modes);
		else
			memcpy(dat->best_modes[k], dat->seeds[k],
				dat->n_coordinates * sizeof **dat->best_modes);

	if (dat->use_ini)
		memcpy(dat->best_seed_idx, dat->ini_seed_idx, opt->K * sizeof
			*dat->best_seed_idx);
	else
		memcpy(dat->best_seed_idx, dat->seed_idx, opt->K * sizeof
			*dat->best_seed_idx);

	memcpy(dat->best_cluster_id, dat->cluster_id, dat->n_observations
		* sizeof *dat->best_cluster_id);
	memcpy(dat->best_cluster_size, dat->cluster_size, opt->K
		* sizeof *dat->best_cluster_size);
	memcpy(dat->best_criterion, dat->criterion, opt->K
		* sizeof *dat->best_criterion);
	if (opt->shuffle)
		memcpy(dat->best_obsn_idx, dat->obsn_idx, dat->n_observations
			* sizeof *dat->obsn_idx);
} /* stash_state */


/**
 * Restore state from previous run.
 *
 * @param dat		data object
 * @param opt		options object
 * @return		error status
 *
 */
int restore_state(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	FILE *fp = fopen(opt->soln_file, "r");
	if (fp == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->soln_file);

	int found_best_modes = 0;
	unsigned int i, k;
	char temp[80];
	long int lline = 0;
	double req_sec = 0;
	while (fgets(temp, 80, fp) != NULL) {
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "restore_state: "
								"%s\n", temp);
		if (!strncmp(temp, "Best optimized criterion:",	/* 25 */
			strlen("Best optimized criterion:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Best optimized criterion: %lf",
				&dat->best_total) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
			for (k = 0; k < opt->K; ++k)
				if (fscanf(fp, "%lf", &dat->best_criterion[k])
					!= 1)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
		} else if (!strncmp(temp, "Best cluster sizes:", /* 19 */
			strlen("Best cluster sizes:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Best cluster sizes: %u",
				&dat->best_cluster_size[0]) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
			for (k = 1; k < opt->K; ++k)
				if (fscanf(fp, "%u", &dat->best_cluster_size[k])
					!= 1)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
		} else if (!strncmp(temp, "Best solution originating",/* 25 */
			strlen("Best solution originating"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Best solution originating seeds: %u",
				&dat->best_seed_idx[0]) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
			for (k = 1; k < opt->K; ++k)
				if (fscanf(fp, "%u", &dat->best_seed_idx[k])
					!= 1)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
		} else if (!strncmp(temp, "Best solution cluster assignments:",/* 34 */
			strlen("Best solution cluster assignments:"))) {
			/* data on newline */
			for (i = 0; i < dat->n_observations; ++i)
				if (fscanf(fp, "%u", &dat->best_cluster_id[i])
					!= 1)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
		} else if (!strncmp(temp, "Best solution indexing:",/* 23 */
			strlen("Best solution indexing:"))) {
			/* data on newline */
			if (!dat->best_obsn_idx) {
				dat->best_obsn_idx = malloc(dat->n_observations
					* sizeof *dat->best_obsn_idx);
				if (dat->best_obsn_idx == NULL)
					return mmessage(ERROR_MSG,
						MEMORY_ALLOCATION,
						"data:best_obsn_idx.\n");
			}
			for (i = 0; i < dat->n_observations; ++i)
				if (fscanf(fp, "%u", &dat->best_obsn_idx[i])
					!= 1)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
		} else if (!strncmp(temp, "Best modes:",/* 11 */
			strlen("Best modes:"))) {
			found_best_modes = 1;
			debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Found "
							"\"Best modes:\".\n");
			/* data on newline */
			for (k = 0; k < opt->K; ++k) {
				if (fscan_data_ts(fp, dat->best_modes[k],
					dat->n_coordinates))
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR,
						opt->soln_file);
				debug_msg(DEBUG_I <= fxn_debug, fxn_debug,
					"Best modes, read mode %u:\n", k);
				debug_call(DEBUG_I <= fxn_debug, fxn_debug,
					fprint_data_ts(stderr,
						dat->best_modes[k],
						dat->n_coordinates, 0, 1));
			}
		} else if (!strncmp(temp, "Maximum AR:",/* 11 */
			strlen("Maximum AR:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Maximum AR: %lf", &dat->best_rand) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
		} else if (!strncmp(temp, "Time to better:",/* 15 */
			strlen("Time to better:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Time to better: %lf +/- %lf; times=%u",
				&dat->avg_time, &dat->sd_time, &dat->ntimes)
				!= 3)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
		} else if (!strncmp(temp, "Run for",/* 7 */
			strlen("Run for"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			fscanf(fp, "Run for %lf", &req_sec);
			if (fgetc(fp) != 's')
				req_sec = 0;
		} else if (!strncmp(temp, "Cost scaling factor:",/* 20 */
			strlen("Cost scaling factor:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Cost scaling factor: %lf",
				&dat->first_cost) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
		} else if (!strncmp(temp, "Maximum cost:",/* 13 */
			strlen("Maximum cost:"))) {
			fseek(fp, -strlen(temp), SEEK_CUR);
			if (fscanf(fp, "Maximum cost: %lf", &dat->worst_cost)
				!= 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->soln_file);
		}

		/* fast forward through current line */
		lline = strlen(temp);
		char c = temp[strlen(temp) - 1];
		while (c != '\n' && !feof(fp)) {
			c = fgetc(fp);
			++lline;
		}
	}
	if (!found_best_modes)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "Missing "
			"\"Best modes:\" line in '%s'.\n", opt->soln_file);
	fseek(fp, -lline, SEEK_CUR);
	if (fscanf(fp, "%*f %lf %lf %lf %lf", &dat->avg_iter, &dat->sd_iter,
		&dat->avg_cost, &dat->sd_cost) != 4)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, opt->soln_file);
	if ((opt->simulate || opt->true_cluster)
		&& fscanf(fp, "%lf %lf %lf %lf %lf %lf", &dat->avg_ar,
			&dat->sd_ar, &dat->avg_mi, &dat->sd_mi, &dat->avg_vi,
			&dat->sd_vi) != 6)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				opt->soln_file);

	if (fscanf(fp, "%lf", &dat->seconds) != 1)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, opt->soln_file);
	if (fscanf(fp, "%u", &dat->n_init) != 1)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, opt->soln_file);
	fclose(fp);

	/* cumulative sums for cost are scaled to avoid overflow */
	if (req_sec > 0) {
		dat->avg_cost /= dat->first_cost;
		dat->sd_cost /= dat->first_cost;
	} else {
		dat->avg_cost /= (dat->n_init + opt->n_init);
		dat->sd_cost /= (dat->n_init + opt->n_init);
	}

	/* invert calculation of mean and sample standard deviation */
	dat->sd_cost = (pow(dat->sd_cost, 2) * dat->n_init
		+ dat->avg_cost * dat->avg_cost) * dat->n_init;
	dat->avg_cost *= dat->n_init;

	dat->sd_iter = (pow(dat->sd_iter, 2) * dat->n_init
		+ dat->avg_iter * dat->avg_iter) * dat->n_init;
	dat->avg_iter *= dat->n_init;

	/* optionally invert mean and sample standard deviation of ar, mi, vi */
	if (opt->simulate || opt->true_cluster) {
		dat->sd_ar = (pow(dat->sd_ar, 2) * dat->n_init
			+ dat->avg_ar * dat->avg_ar) * dat->n_init;
		dat->avg_ar *= dat->n_init;
		dat->sd_mi = (pow(dat->sd_mi, 2) * dat->n_init
			+ dat->avg_mi * dat->avg_mi) * dat->n_init;
		dat->avg_mi *= dat->n_init;
		dat->sd_vi = (pow(dat->sd_vi, 2) * dat->n_init
			+ dat->avg_vi * dat->avg_vi) * dat->n_init;
		dat->avg_vi *= dat->n_init;
	}

	return NO_ERROR;
} /* restore_state */


/**
 * Open and write initial state to output files.
 *
 * @param dat	data object
 * @param opt	options object
 * @param fps	solution file
 * @param fpi	initialization file
 * @return	error status
 */
int initialize_outfiles(data *dat, options * opt, FILE **in_fps, FILE **in_fpi)
{
	/* open solution file */
	if (opt->continue_run)	/* continuing a previous run */
		*in_fps = fopen(opt->soln_file, "a");
	else			/* starting a fresh run */
		*in_fps = fopen(opt->soln_file, "w");
	if (*in_fps == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->soln_file);
	FILE *fps = *in_fps;

	/* open initialization file */
	if (opt->ini_file) {
		if (opt->continue_run)
			*in_fpi = fopen(opt->ini_file, "a");
		else
			*in_fpi = fopen(opt->ini_file, "w");
		if (*in_fpi == NULL)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					opt->ini_file);
	}
	FILE *fpi = *in_fpi;

	/* write initial state to solution file */
	if (opt->n_init || opt->seconds > 0) {
		if (opt->continue_run)
			fprintf(fps, "***Continuing run\n");
		if (opt->K > 0) {
			fprintf(fps, "Seed: %lu\n", opt->seed);
			fprintf(fps, "Algorithm: %s%s\n",
				kmodes_algorithm(opt->kmodes_algorithm),
					opt->update_modes ? " with mode updates"
					: "");
			fprintf(fps, "Initialization: %s%s\n",
				kmodes_init_method(opt->init_method),
				opt->shuffle ? " with shuffling." : ".");
			if (opt->seconds > 0)
				fprintf(fps, "Run for %fs\n", opt->seconds);
			else
				fprintf(fps, "Run for %u initializations\n",
					opt->n_init);
		}

		fprintf(fps, "Number of clusters to estimate: %u\n", opt->K);
		if (fpi)
			fprintf(fps, "Initialization output file: %s\n",
				opt->ini_file);
	}
	/* record run information in outfile */
	if (opt->simulate) {
		fprintf(fps, "Number of observations: %u\n",
			dat->n_observations);
		fprintf(fps, "Number of coordinates: %u\n", dat->n_coordinates);
		fprintf(fps, "Number of categories: %u\n",
			opt->sim_n_categories);
		fprintf(fps, "Number of true clusters: %u (%u realized)\n",
			opt->sim_K, opt->true_K);
		fprintf(fps, "CTMC times: %f %f\n", opt->sim_between_t,
			opt->sim_within_t);
		fprintf(fps, "CTMC probabilities: %f %f\n",
			opt->sim_between_prob, opt->sim_within_prob);
		if (opt->sim_alpha) {
			fprintf(fps, "Dirichlet prior (alpha):\n");
			fprint_doubles(fps, opt->sim_alpha, opt->sim_K, 6, 1);
		}
		fprintf(fps, "Mixing proportions:\n");
		fprint_doubles(fps, opt->sim_pi, opt->sim_K, 6, 1);
		fprintf(fps, "Simulated modes:\n");
		for (unsigned int k = 0; k < opt->sim_K; ++k)
			fprint_data_ts(fps, opt->sim_modes[k],
				dat->n_coordinates, (int) (log10(
				opt->sim_n_categories
				- (opt->sim_n_categories>1))+1), 1);
		fprintf(fps, "Mode pairwise distances:\n");
		for (unsigned int k = 0; k < opt->sim_K - 1; ++k) {
			for (unsigned int j = 1; j <= k; ++j)
				fprintf(fps, " %*s",
					(int)(log10(dat->n_coordinates) + 1),
					" ");
			for (unsigned int j = k + 1; j < opt->sim_K; ++j)
				fprintf(fps, " %*.0f",
					(int)(log10(dat->n_coordinates) + 1),
					hd(opt->sim_modes[k], opt->sim_modes[j],
					dat->n_coordinates));
			fprintf(fps, "\n");
		}
		fprintf(fps, "Simulated cluster assignments:\n");
		fprint_uints(fps, opt->sim_cluster, dat->n_observations,
			(int)(log10(opt->sim_K - (opt->sim_K > 1)) + 1), 1);
		fprintf(fps, "Simulated cluster sizes:\n");
		fprint_uints(fps, opt->true_cluster_size, opt->sim_K,
			(int)(log10(dat->n_observations - opt->sim_K + 1) + 1), 1);
		if (opt->datafile)
			fprintf(fps, "Data written to file: %s\n",
				opt->datafile);

	} else if (opt->true_cluster) {
		fprintf(fps, "Number of true clusters: %u\n", opt->true_K);
		if (opt->true_modes) {
			fprintf(fps, "True modes:\n");
			for (unsigned int k = 0; k < opt->true_K; ++k)
				fprint_data_ts(fps, opt->true_modes[k],
					dat->n_coordinates, (int) (log10(
					opt->sim_n_categories
					- (opt->sim_n_categories>1))+1), 1);
			fprintf(fps, "Mode pairwise distances:\n");
			for (unsigned int k = 0; k < opt->true_K; ++k) {
				for (unsigned int j = k + 1; j
					< opt->true_K; ++j)
					fprintf(fps, " %.0f", hd(
						opt->true_modes[k],
						opt->true_modes[j],
						dat->n_coordinates));
				fprintf(fps, "\n");
			}
		}
		fprintf(fps, "True cluster assignments:\n");
		fprint_uints(fps, opt->true_cluster, dat->n_observations,
			(int)(log10(opt->true_K - (opt->true_K > 1)) + 1), 1);
		fprintf(fps, "True cluster sizes:\n");
		fprint_uints(fps, opt->true_cluster_size, opt->true_K,
			(int)(log10(dat->n_observations - opt->true_K + 1) + 1), 1);
	} else {
		fprintf(fps, "Input file: %s\n", opt->datafile);
		fprintf(fps, "Number of observations: %u\n",
			dat->n_observations);
		fprintf(fps, "Number of coordinates: %u\n", dat->n_coordinates);
		fprintf(fps, "Number of categories: ");
		fprint_data_ts(fps, dat->n_categories, dat->n_coordinates,
			(int)(log10(dat->max_n_categories -
			(dat->max_n_categories>1))+1), 1);
		fprintf(fps, "Total categories: %zu\n", dat->tot_n_categories);
	}

	return NO_ERROR;
} /* initialize_outfiles */


/**
 * Write status to stdout and output files.
 *
 * @param dat	data object
 * @param opt	options object
 * @param fps	solution file pointer
 * @param fpi	initialization file pointer
 * @param i	current initialization
 * @param start	epoch program started
 * @return	number of seconds spend in this function
 */
void write_status(data *dat, options *opt, FILE *fps, FILE *fpi,            /**/
	unsigned int i, TIME_STRUCT *start)
{
	TIME_STRUCT stop;	/* timing */

	MARK_TIME(&stop);

	/* output information about found solution to stdout */
	if (opt->quiet >= MINIMAL)
		fprintf(stdout, "Init %*u (%*u iterations):",
			(int)(log10(opt->n_init) + 1), (opt->seconds > 0
			? opt->n_init : i),
			(int)(log10(dat->max_iter) + 1), dat->iter);

	/* internal settings */
	if (dat->iter > dat->max_iter)
		dat->max_iter = dat->iter;

	for (unsigned int k = 0; k < opt->K; ++k) {
		if (opt->quiet >= MINIMAL)
			fprintf(stdout, " %*.0f", (int)(log10(dat->worst_cost)
				+ 1), dat->criterion[k]);
	}

	if (dat->worst_cost < dat->total)
		dat->worst_cost = dat->total;
	if (opt->seconds > 0 && dat->first_cost == 1)
		dat->first_cost = dat->total;

	if (opt->quiet >= MINIMAL)
		fprintf(stdout, ": %*.0f (%*.0f)", (int)(log10(dat->worst_cost)
			+ 1), dat->total, (int)(log10(dat->worst_cost) + 1),
			dat->best_total);
	dat->avg_cost += dat->total / (opt->seconds > 0 ? dat->first_cost
			: (dat->n_init + opt->n_init));
	dat->sd_cost += dat->total * dat->total / (opt->seconds > 0
		? dat->first_cost * dat->first_cost
		: (dat->n_init + opt->n_init) * (dat->n_init + opt->n_init));
	dat->avg_iter += dat->iter;
	dat->sd_iter += (double) dat->iter * dat->iter;

	unsigned int *tc = opt->true_cluster;
	unsigned int tk = opt->true_K;
	double ar = 0, mi = 0, vi = 0;
	if (tc) {
		mi = mutual_information(dat->cluster_id, tc,
			dat->n_observations, opt->K, tk,
			NORMALIZED_MUTUAL_INFORMATION, MAX_SCALING);
		vi = mutual_information(dat->cluster_id, tc,
			dat->n_observations, opt->K, tk,
			NORMALIZED_VARIATION_OF_INFORMATION, NO_MI_SCALING);
		if (opt->quiet >= MINIMAL)
			fprintf(stdout, " %.3f %.3f", mi, vi);
		ar = cluster_index(dat->cluster_id, tc, dat->n_observations,
			opt->K, tk, ADJUSTED_RAND_INDEX);
		dat->avg_mi += mi;
		dat->avg_vi += vi;
		dat->avg_ar += ar;
		dat->sd_mi += mi * mi;
		dat->sd_vi += vi * vi;
		dat->sd_ar += ar * ar;
	}

	if ((opt->simulate || opt->true_cluster) && opt->quiet >= MINIMAL)
		fprintf(stdout, " %6.3f (%.3f)", ar, dat->best_rand);
	if (opt->quiet >= MINIMAL)
		fprintf(stdout, "\n");

	/* output initialization information */
	if (fpi || fps) {
		fprintf(fpi ? fpi : fps, "%u %u", opt->seconds > 0 ? opt->n_init
			: i, dat->iter);
		fprint_doubles(fpi ? fpi : fps, dat->criterion, opt->K, 0, 0);
		fprintf(fpi ? fpi : fps, " %.0f %.0f", dat->total,
			dat->best_total);
		if (opt->simulate || opt->true_cluster)
			fprintf(fpi ? fpi : fps, " %f %f %f %f", ar,
				dat->best_rand, mi, vi);
	}

	/* record information about times beating target */
	++dat->ctime;
	if (dat->total <= opt->target) {
		dat->avg_time += dat->ctime;
		dat->sd_time += dat->ctime * dat->ctime;
		dat->ctime = 0;
		dat->ntimes++;
	}

	if ((opt->simulate || opt->true_cluster) && ar > dat->best_rand)
		dat->best_rand = ar;

	/* the rest is necessary bookkeeping */
	dat->uncounted_seconds += ELAP_TIME(&stop);

	if (fpi || fps)
		fprintf(fpi ? fpi : fps, " %f\n", ELAP_TIME(start) - dat->uncounted_seconds);

} /* write_status */


/**
 * Create the options object and set defaults.
 *
 * @param opt	options object (to be allocated)
 * @return	error status
 */
int make_options(options **opt) {
	*opt = malloc(sizeof **opt);

	if (*opt == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");

	(*opt)->K = 0;	/* invalid value */
	(*opt)->true_column = UINT_MAX;
	(*opt)->true_cluster = NULL;
	(*opt)->true_cluster_size = NULL;
	(*opt)->true_K = 0;
	(*opt)->n_init = 1;
	(*opt)->n_inner_init = 1;
	(*opt)->n_max_iter = 10000;
	(*opt)->info = QUIET;
	(*opt)->subtract_one = 0;
	(*opt)->shuffle = 0;
	(*opt)->kmodes_algorithm = KMODES_HUANG;
	(*opt)->init_method = KMODES_INIT_RANDOM_SEEDS;
	(*opt)->continue_run = 0;
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
	(*opt)->seed = 0;
	(*opt)->seed_idx = NULL;
	(*opt)->seed_set = NULL;
	(*opt)->n_sd_idx = 0;
	(*opt)->n_seed_set = 0;
	(*opt)->datafile = NULL;
	(*opt)->data_outfile = NULL;
	(*opt)->ini_file = NULL;
	(*opt)->soln_file = NULL;
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

int parse_options(options *opt, int argc, const char **argv)
{
	int fxn_debug = ABSOLUTE_SILENCE;
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
		debug_msg(MINIMAL <= fxn_debug, opt->quiet, "%u -k[0-9] "
						"arguments\n", opt->n_k);
		opt->result_files = malloc(opt->n_k *
			sizeof *opt->result_files);
		opt->n_result_files = malloc(opt->n_k *
			sizeof *opt->n_result_files);
		if (opt->result_files == NULL || opt->n_result_files == NULL)
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
				debug_msg(MINIMAL <= fxn_debug, opt->quiet,
					"true column = %u\n", opt->true_column);
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
					i += opt->n_result_files[k - opt->min_k];
				} else {
					opt->K = read_uint(argc, argv, ++i,
						(void *)opt);
					if (errno)
						goto CMDLINE_ERROR;
					debug_msg(MINIMAL <= fxn_debug,
						opt->quiet, "K = %u.\n", opt->K);//
				}
				break;
			case 'l':
				opt->kmodes_algorithm = KMODES_LLOYD;
				debug_msg(QUIET <= fxn_debug, opt->quiet,
					"Using Lloyd's algorithm.\n");
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
				debug_msg(QUIET <= fxn_debug, opt->quiet,
					"Using Hartigan and Wong algorithm.\n");
				break;
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->datafile = argv[++i];
				if (i + 1 < argc && argv[i+1][0] != '-') {
					opt->data_outfile = argv[++i];
					debug_msg(QUIET <= fxn_debug, opt->quiet,
						"Will write data, after possible "
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
				else if (access(argv[i+1],  F_OK) == -1) {

					opt->init_method = KMODES_INIT_USER_SEEDS;

					opt->n_sd_idx = read_cmdline_uints(argc,
						argv, ++i, &opt->seed_idx,
						(void *)opt);
					if (errno || !opt->n_sd_idx)
						goto CMDLINE_ERROR;
					debug_msg(MINIMAL <= fxn_debug,
						opt->quiet, "Seed indices:");
					debug_call(MINIMAL <= fxn_debug,
						opt->quiet, fprint_uints(stderr,
						opt->seed_idx, opt->n_sd_idx,
						1, 1));
					i += opt->n_sd_idx - 1;

				/* assume seeds are provided in a file */
				} else {
					opt->sfile = argv[i+1];
					//opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
				}
				debug_msg(QUIET <= fxn_debug, opt->quiet,
					"Using %s initialization.\n",
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
				debug_msg(MINIMAL <= fxn_debug, opt->quiet,
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
					debug_msg(MINIMAL <= fxn_debug,
						opt->quiet, "Inner "
						"initializations; %u\n",
						opt->n_inner_init);
				} else {
					opt->seed = read_uint(argc, argv, ++i,
						(void *)opt);
					srand(opt->seed);
					debug_msg(MINIMAL <= fxn_debug,
						opt->quiet, "Seed: "
						"%lu\n", opt->seed);
				}
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 's':	/*-s <sn> <sp> <sc> <st1> <st2>*/
				/* no appropriate arguments */
				if (!strcmp(&argv[i][j], "shuffle")) {
					opt->shuffle = 1;
					debug_msg(MINIMAL <= fxn_debug,
						opt->quiet, "Data will be "
						"shuffled.\n");
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
				debug_msg(MINIMAL <= fxn_debug, opt->quiet,
					"Simulation:\n\t%u observations\n\t%u "
					"coordinates\n\t%u categories\n\t%f "
					"between variance\n\t%f within "
					"variance\n", opt->sim_n_observations,
					opt->sim_n_coordinates,
					opt->sim_n_categories,
					opt->sim_between_t, opt->sim_within_t);

				break;
			case 't':
				opt->seconds = read_cmdline_double(argc, argv,
					++i, (void *)opt);
				debug_msg(QUIET <= fxn_debug, opt->quiet,
					"Running for %.0fs.\n", opt->seconds);
				break;
			case 'u':
				opt->update_modes = 1;
				debug_msg(MINIMAL <= fxn_debug, opt->quiet,
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
					debug_msg(QUIET <= fxn_debug, opt->quiet,//
						"Using Huang's algorithm.\n");
				} else if (!strncmp(&argv[i][j], "hart",
							strlen("hart"))) {
					opt->use_hartigan = 1;
					debug_msg(QUIET <= fxn_debug, opt->quiet,
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
				"setting K if estimating K.  Ignoring -k %u.\n",
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
	int fxn_debug = ABSOLUTE_SILENCE;
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
			debug_msg(MINIMAL <= fxn_debug, opt->quiet,
					"Simulation K = %u\n", opt->sim_K);
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
			debug_msg(MINIMAL <= fxn_debug, opt->quiet,
							"Simulation alpha:");
			debug_call(MINIMAL <= fxn_debug, opt->quiet,
				fprint_doubles(stderr, opt->sim_alpha,
							opt->sim_K, 2, 1));
		} else {
			double sum = 0;
			for (unsigned int j = 1; j < opt->sim_K; ++j)
				sum += opt->sim_pi[j];
			if (sum >= 1.)
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"-pi <pdbl1> ... <pdblK> arguments must"
					"sum to 1.0.\n");
			opt->sim_pi[0] = 1 - sum;
			debug_msg(MINIMAL <= fxn_debug, opt->quiet,
							"Simulation pi:");
			debug_call(MINIMAL <= fxn_debug, opt->quiet,
				fprint_doubles(stderr, opt->sim_pi,
							opt->sim_K, 2, 1));
		}
	} else if (access(argv[*i + 1], F_OK) == -1)
		opt->n_effective_coordinates = read_cmdline_double(argc, argv,
			++(*i), (void *)opt);
	else
		opt->pfile = argv[++(*i)];

	return NO_ERROR;
} /* process_arg_p */


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

int make_data(data **dat, options *opt) {
	int err = NO_ERROR;
	*dat = malloc(sizeof **dat);

	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");

	(*dat)->data = NULL;
	(*dat)->dmat = NULL;
	(*dat)->n_observations = 0;
	(*dat)->n_coordinates = 0;
	(*dat)->n_categories = NULL;
	(*dat)->max_n_categories = 0;
	(*dat)->seeds = NULL;
	(*dat)->ini_seeds = NULL;
	(*dat)->seed_idx = NULL;
	(*dat)->best_seed_idx = NULL;
	(*dat)->ini_seed_idx = NULL;
	(*dat)->cluster_id = NULL;
	(*dat)->obsn_idx = NULL;
	(*dat)->best_obsn_idx = NULL;
	(*dat)->best_cluster_id = NULL;
	(*dat)->best_modes = NULL;
	(*dat)->criterion = NULL;
	(*dat)->best_criterion = NULL;
	(*dat)->cluster_size = NULL;
	(*dat)->best_cluster_size = NULL;
	if ((err = allocate_data_for_k(*dat, opt->K)))
		return err;

	(*dat)->best_total = INFINITY;
	(*dat)->best_rand = -INFINITY;
	(*dat)->n_init = 0;
	(*dat)->iter = 0;
	(*dat)->max_iter = 1;
	(*dat)->seconds = 0.;
	(*dat)->uncounted_seconds = 0.;
	(*dat)->first_cost = 1;
	(*dat)->worst_cost = 1;
	(*dat)->avg_cost = (*dat)->sd_cost = 0;
	(*dat)->avg_iter = (*dat)->sd_iter = 0;
	(*dat)->avg_time = (*dat)->sd_time = 0;
	(*dat)->ctime = 0;
	(*dat)->ntimes = 0;
	(*dat)->avg_ar = (*dat)->sd_ar = 0;
	(*dat)->avg_mi = (*dat)->sd_mi = 0;
	(*dat)->avg_vi = (*dat)->sd_vi = 0;
	(*dat)->n_categories = NULL;
	(*dat)->tot_n_categories = 0;

	return NO_ERROR;
} /* make_data */

int allocate_data_for_k(data *dat, unsigned int k)
{
	if (!k)
		return NO_ERROR;

	dat->criterion = malloc(k * sizeof *dat->criterion);
	dat->best_criterion = malloc(k * sizeof *dat->best_criterion);
	dat->cluster_size = malloc(k * sizeof *dat->cluster_size);
	dat->best_cluster_size = malloc(k * sizeof *dat->best_cluster_size);

	if (dat->criterion == NULL || dat->best_criterion == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.criterion");
	if (dat->cluster_size == NULL || dat->best_cluster_size == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"data.cluster_size");

	/* WARNING: not allocating seed_idx assuming no new initializations */
	dat->best_seed_idx = malloc(k * sizeof *dat->best_seed_idx);

	if (dat->best_seed_idx == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"data.best_seed_idx");

	/* WARNING: not allocating seeds assuming no new initializations */

	if (dat->n_coordinates) {
		dat->best_modes = malloc(k * sizeof *dat->best_modes);
		data_t *tmp = malloc(k*dat->n_coordinates * sizeof *tmp);
		if (dat->best_modes == NULL || tmp == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data.best_modes");
		for (unsigned int i = 0; i < k; ++i) {
			dat->best_modes[i] = tmp;
			tmp += dat->n_coordinates;
		}
	}

	return NO_ERROR;
} /* allocate_data_for_k */

/**
 * Free data that depends on options::K.
 *
 * @param dat	data object pointer
 */
void free_data_for_k(data *dat, options *opt)
{
	if (dat->criterion) {
		free(dat->criterion);
		dat->criterion = NULL;
	}
	if (dat->best_criterion) {
		free(dat->best_criterion);
		dat->best_criterion = NULL;
	}
	if (dat->cluster_size) {
		free(dat->cluster_size);
		dat->cluster_size = NULL;
	}
	if (dat->best_cluster_size) {
		free(dat->best_cluster_size);
		dat->best_cluster_size = NULL;
	}

	if (dat->seed_idx && !opt->seed_idx) {
		free(dat->seed_idx);
		dat->seed_idx = NULL;
	}
	if (dat->best_seed_idx) {
		free(dat->best_seed_idx);
		dat->best_seed_idx = NULL;
	}
	if (dat->seeds) {
		if (dat->seeds[0])
			free(dat->seeds[0]);
		free(dat->seeds);
		dat->seeds = NULL;
	}
	if (dat->best_modes) {
		if (dat->best_modes[0])
			free(dat->best_modes[0]);
		free(dat->best_modes);
		dat->best_modes = NULL;
	}
} /* free_data_for_k */

int read_data(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	unsigned int j = 0, i;
	char c, pc = 0;
	FILE *fp = fopen(opt->datafile, "r");
	data_t *dptr;

	if (fp == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->datafile);

	dat->n_coordinates = 1;
	dat->n_observations = 0;

	/* count number of columns, assuming space-separated fields,
	 * allowing possibility of no terminal newline */
	while (!feof(fp)) {
		c = fgetc(fp);
		/* new line starting */
		if (!feof(fp) && c != '\n' && (!pc || pc == '\n'))
		       dat->n_observations++;
		/* new column starting */
		if (dat->n_observations == 1 && pc && c == ' ' && pc != ' ')
		       dat->n_coordinates++;
		pc = c;
	}

	debug_msg(MINIMAL <= fxn_debug, opt->quiet, "Data %u x %u\n",
				dat->n_observations, dat->n_coordinates);

	if (opt->true_column < dat->n_coordinates) {
		--dat->n_coordinates;
		opt->true_cluster = malloc(dat->n_observations
			* sizeof *opt->true_cluster);
		if (!opt->true_cluster) {
			err = mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"options::true_cluster");
			goto ABORT_READ_DATA;
		}
	}

	dat->n_categories = calloc(dat->n_coordinates,
		sizeof *dat->n_categories);
	dat->data = malloc(dat->n_coordinates * dat->n_observations
		* sizeof *dat->data);
	if (dat->data == NULL || dat->n_categories == NULL) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::data");
		goto ABORT_READ_DATA;
	}

	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations
		* sizeof *dat->best_cluster_id);

	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"data::cluster_id");

	rewind(fp);

	/* read in data */
	/* assume categories are 0, 1, 2, ... without skips */
	dptr = dat->data;
	if (opt->true_column > dat->n_coordinates) {
		unsigned int nc = dat->n_coordinates;
		while (fscanf(fp, "%" AMPLICLUST_SCNu_data_t, dptr) == 1) {
			if (opt->subtract_one) -- (*dptr);
			if (dat->n_categories[j] < *dptr + 1u)
				dat->n_categories[j] = *dptr + 1u;
			++ dptr;
			j = (j + 1u) % nc;
		}
	} else {
		unsigned int nc = dat->n_coordinates + 1u;
		i = 0;
		while (fscanf(fp, "%" AMPLICLUST_SCNu_data_t, dptr) == 1) {
			if (opt->subtract_one) -- (*dptr);
			if (j == opt->true_column) {
				opt->true_cluster[i ++] = *dptr;
				if (*dptr + 1u > opt->true_K)
					opt->true_K = *dptr + 1u;
			} else {
				if (*dptr + 1u > dat->n_categories[j
					- (j > opt->true_column)])
					dat->n_categories[j - (j
						> opt->true_column)]
						= *dptr + 1u;
				++ dptr;
			}
			j = (j + 1u) % nc;
		}
	}
	fclose(fp);
	fp = NULL;

	/* true modes provided: use to estimate options:true_K */
	if (opt->mfile) {

		/* open file */
		fp = fopen(opt->mfile, "r");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
				"temporary file");

		/* count number of lines if file: options:true_K */
		unsigned int mode_K = 0;
		do {
			char c = fgetc(fp);

			/* new line with content */
			if (c != '\n' && !feof(fp))
				++mode_K;

			/* fast-forward through line */
			while (c != '\n' && !feof(fp)) c = fgetc(fp);
		} while (!feof(fp));

		if (!opt->true_K)
			opt->true_K = mode_K;

		if (mode_K != opt->true_K) {
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Number of modes does not match number of "
				"clusters in data file '%s'.\n", opt->datafile);
			goto ABORT_READ_DATA;
		}

		rewind(fp);

		/* allocate memory */
		data_t *tmp = malloc(opt->true_K * dat->n_coordinates * sizeof *tmp);
		if (!tmp) {
			err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "tmp");
			goto ABORT_READ_DATA;
		}
		opt->true_modes = malloc(opt->true_K * sizeof *opt->true_modes);
		if (!opt->true_modes) {
			err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "tmp");
			goto ABORT_READ_DATA;
		}

		/* read modes */
		for (unsigned int k = 0; k < opt->true_K; ++k) {
			opt->true_modes[k] = tmp;
			tmp += dat->n_coordinates;
			if (fscan_data_ts(fp, opt->true_modes[k],
				dat->n_coordinates)) {
				err = mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->mfile);
				goto ABORT_READ_DATA;
			}
		}
		fclose(fp);
		fp = NULL;

		if (opt->true_K > dat->n_observations) {
			err = mmessage(ERROR_MSG, INTERNAL_ERROR, "Your "
				"friendly programmer has assumed the number "
				"of clusters is smaller than the number of "
				"observations.\n");
			goto ABORT_READ_DATA;
		}

		/* count true number of clusters (modes may be identical) */
		opt->sim_n_categories = 0;	/* use to store categories */
		unsigned int distinct_k = 1;
		dat->cluster_id[0] = 0;/* tmp use: map old cluster to new */
		for (unsigned int k = 1; k < opt->true_K; ++k) {
			unsigned int same = k;
			for (unsigned int j = 0; j < k; ++j) {
				same = j;
				for (unsigned int i = 0; i < dat->n_coordinates;
					++i) {
					if (!j && opt->true_modes[k][i] + 1u
						> opt->sim_n_categories)
						opt->sim_n_categories =
							opt->true_modes[k][i] + 1u;
					if (opt->true_modes[k][i]
						!= opt->true_modes[j][i]) {
						same = k;
						break;
					}
				}
				if (same < k)
					break;
			}
			if (same == k) {
				/* IMPORTANT: opt->true_modes[0] NOT reassigned,
				 * otherwise memory true_modes memory block would
				 * become unreachable */
				if (distinct_k < k)
					opt->true_modes[distinct_k] = opt->true_modes[k];
				dat->cluster_id[k] = distinct_k++;
			} else
				dat->cluster_id[k] = dat->cluster_id[same];
		}

		if (distinct_k < opt->true_K)
			mmessage(WARNING_MSG, NO_ERROR, "Mode file contains %u "
				"modes, but only %u are unique.\n", opt->true_K,
				distinct_k);
		opt->true_K = distinct_k;
	} else if (opt->true_cluster)	/* can only guess all modes are distinct */
		for (unsigned int k = 0; k < opt->true_K; ++k)
			dat->cluster_id[k] = k;	/* temporary use */

	/* count cluster sizes : one more opportunity to lose a cluster */
	if (opt->true_cluster) {
		opt->true_cluster_size = calloc(opt->true_K, sizeof
			*opt->true_cluster_size);
		if (!opt->true_cluster_size)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"options:true_cluster_size");

		/* count cluster sizes using (new) indices */
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			++opt->true_cluster_size[dat->cluster_id[
				opt->true_cluster[i]]];

		/* count clusters with members */
		unsigned int nonempty_k = 0;
		for (unsigned int k = 0; k < opt->true_K; ++k) {
			opt->true_cluster_size[nonempty_k]
				= opt->true_cluster_size[k];
			if (opt->true_cluster_size[k]) {
				if (opt->true_modes && !nonempty_k && k)
					memcpy(opt->true_modes[0],
						opt->true_modes[k],
					       	dat->n_coordinates
					       	* sizeof **opt->true_modes);
				else if (opt->true_modes)
					opt->true_modes[nonempty_k]
						= opt->true_modes[k];
				dat->cluster_id[k] = nonempty_k++;
			} else
				dat->cluster_id[k] = nonempty_k;
		}
		if (nonempty_k < opt->true_K) {
			mmessage(WARNING_MSG, NO_ERROR, "Mode file contained "
				"%u distinct modes, but only %u have data.\n",
				opt->true_K, nonempty_k);
			opt->true_K = nonempty_k;
			for (unsigned int i = 0; i < dat->n_coordinates; ++i)
				opt->true_cluster[i] =
					dat->cluster_id[opt->true_cluster[i]];
		}
	}

	/* fix assumption of 0, 1, 2, ... categories */
	/* SLOW: should write an option to bypass it! */
	for (j = 0; j < dat->n_coordinates; ++j) {
		unsigned int *present = calloc(dat->n_categories[j], sizeof *present);
		for (i = 0; i < dat->n_observations; ++i)
			++present[dat->data[dat->n_coordinates*i + j]];
		unsigned int cnt = 0;
		for (i = 0; i < dat->n_categories[j]; ++i) {
			if (present[i]) ++cnt;
			present[i] = (i ? present[i-1] : 0) + (present[i] == 0);
		}
		if (cnt == dat->n_categories[j]) {
			free(present);
			continue;
		}
		mmessage(WARNING_MSG, NO_ERROR, "Coordinate %u uses only %u "
			"categories, but has %u categories.  To get correct "
			"category counts, use double arguments to -o and maybe "
			"-m, and then run with corrected files.\n", j,
			cnt, dat->n_categories[j]);
		dat->n_categories[j] = cnt;
		for (i = 0; i < dat->n_observations; ++i)
			dat->data[dat->n_coordinates*i + j] -= present[dat->data[dat->n_coordinates*i + j]];
		if (opt->true_modes)
			for (unsigned int k = 0; k < opt->true_K; ++k)
				opt->true_modes[k][j] -= present[opt->true_modes[k][j]];
		free(present);
	}
	for (j = 0; j < dat->n_coordinates; ++j) {
		dat->tot_n_categories += dat->n_categories[j];
		if (dat->n_categories[j] > dat->max_n_categories)
			dat->max_n_categories = dat->n_categories[j];
	}

	/* write modes (presumably after conversion above) */
	if (opt->mfile_out) {
		fp = fopen(opt->mfile_out, "w");
		for (unsigned int k = 0; k < opt->true_K; ++k)
			fprint_data_ts(fp, opt->true_modes[k],
				dat->n_coordinates, 0, 1);
		fclose(fp);
		fp = NULL;
	}

ABORT_READ_DATA:
	if (fp) fclose(fp);

	return err;
} /* read_data */

int simulate_data(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	unsigned int k;
	double dsum;
	double factor = (double) opt->sim_n_categories * (opt->sim_n_categories - 1);
	size_t dim = opt->sim_n_categories * opt->sim_n_categories;	/* WARNING/BUG: may be very large! */
	double *Pt = NULL;
	double *a = malloc(dim * sizeof *a);

	dat->n_observations = opt->sim_n_observations;
	dat->n_coordinates = opt->sim_n_coordinates;

	/* allocate simulated modes and data */
	data_t *tmp = malloc(opt->sim_K * dat->n_coordinates
		* sizeof **opt->sim_modes);
	opt->sim_modes = malloc(opt->sim_K * sizeof *opt->sim_modes);
	/* allocate space for data */
	dat->data = malloc(dat->n_coordinates * dat->n_observations
		* sizeof *dat->data);
	dat->n_categories = calloc(dat->n_coordinates,
		sizeof *dat->n_categories);
	opt->sim_cluster = malloc(dat->n_observations
		* sizeof *opt->sim_cluster);

	if (!a || !tmp || !opt->sim_modes || !dat->data
		|| !opt->sim_cluster || !dat->n_categories) {
		err = MEMORY_ALLOCATION;
		goto ABORT_SIMULATE_DATA;
	}
	opt->true_modes = opt->sim_modes;
	opt->true_cluster = opt->sim_cluster;

	for (k = 0; k < opt->sim_K; ++k) {
		opt->sim_modes[k] = tmp;
		tmp += dat->n_coordinates;
		dat->n_categories[k] = opt->sim_n_categories;
		dat->cluster_size[k] = 0;
	}

	/* allocate instantaneous rate matrix */
	for (size_t i = 0; i < dim; ++i)
		a[i] = opt->sim_between_t / factor;
	for (unsigned int i = 0; i < opt->sim_n_categories; ++i)
		a[i*opt->sim_n_categories + i] = - (double)
			(opt->sim_n_categories - 1) * opt->sim_between_t
			/ factor;

	/* compute transition probability matrix */
	Pt = r8mat_expm1(opt->sim_n_categories, a);
	if (!Pt) {
		err = MEMORY_ALLOCATION;
		goto ABORT_SIMULATE_DATA;
	}

	debug_msg(QUIET <= fxn_debug, opt->quiet, "Probability of no change: "
								"%f\n", Pt[0]);
	opt->sim_between_prob = Pt[0];

	do {
		/* simulate modes */
		for (unsigned int j = 0; j < dat->n_coordinates; ++j) {
			data_t c_ancestor = (data_t) ((double) rand() / RAND_MAX
				* opt->sim_n_categories);
			/* simulate ancestor character */
			for (k = 0; k < opt->sim_K; ++k) {
				double r = (double) rand() / RAND_MAX;
				data_t l = 0;
				for (dsum = Pt[c_ancestor*opt->sim_n_categories];
					dsum < r; dsum += Pt[c_ancestor
						* opt->sim_n_categories + ++l]);
				opt->sim_modes[k][j] = l;
			}
		}

		/* count number of distinct modes */
		opt->true_K = 1;
		for (unsigned int k = 1; k < opt->sim_K; ++k) {
			int same = 1;
			for (unsigned int j = 0; j < k; ++j) {
				same = 1;
				for (unsigned int i = 0; i < dat->n_coordinates;
					++i)
					if (opt->sim_modes[k][i]
						!= opt->sim_modes[j][i]) {
						same = 0;
						break;
					}
				if (same)
					break;
			}
			if (!same)
				++opt->true_K;
		}
	} while (opt->require_sim_K && opt->true_K < opt->sim_K);

	if (opt->true_K < opt->sim_K) {
		mmessage(WARNING_MSG, NO_ERROR, "Asked to simulate %u modes, "
			"but only %u are unique (possibility not ready).\n",
			opt->sim_K, opt->true_K);
		goto ABORT_SIMULATE_DATA;
	}

	/* allocate space for true cluster sizes */
	opt->true_cluster_size = malloc(opt->true_K *
		sizeof *opt->true_cluster_size);
	if (!opt->true_cluster_size) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"options:true_cluster_size");
		goto ABORT_SIMULATE_DATA;
	}

	/* reset instantaneous rate matrix and transition probability matrix */
	for (size_t i = 0; i < dim; ++i)
		a[i] = opt->sim_within_t / factor;
	for (unsigned int i = 0; i < opt->sim_n_categories; ++i)
		a[i*opt->sim_n_categories + i] = - (double)
			(opt->sim_n_categories - 1) * opt->sim_within_t
			/ factor;

	free(Pt);
	Pt = r8mat_expm1(opt->sim_n_categories, a);
	if (!Pt) {
		err = MEMORY_ALLOCATION;
		goto ABORT_SIMULATE_DATA;
	}
	debug_msg(QUIET <= fxn_debug, opt->quiet, "Probability of no change: "
								"%f\n", Pt[0]);
	opt->sim_within_prob = Pt[0];

	/* verify that the simulated mode is the theoretical mode */
	for (data_t l = 0; l < opt->sim_n_categories; ++l) {
		dsum = 0;
		for (data_t j = 0; j < opt->sim_n_categories; ++j) {
			dsum += Pt[l*opt->sim_n_categories + j];
			if (l != j && Pt[l*opt->sim_n_categories + l]
				< Pt[l*opt->sim_n_categories + j]) {
				err = mmessage(ERROR_MSG,
					INVALID_USER_INPUT, "Within "
					"cluster time (%f) is too "
					"large\n", opt->sim_within_t);
				goto ABORT_SIMULATE_DATA;
			}
		}
		if (fabs(dsum - 1.) > 1e-6) {
			err = mmessage(ERROR_MSG, INTERNAL_ERROR,
				"matrix exponentiation failed\n");
			goto ABORT_SIMULATE_DATA;
		}
	}

	unsigned int true_k = 0;
	do {
		/* simulate pi */
		if (opt->sim_alpha) {
			/* this is the only place the R rng is used */
			set_seed(rand(), rand());
			double sum = 0;
			for (k = 0; k < opt->sim_K; ++k) {
				opt->sim_pi[k] = rgamma(opt->sim_alpha[k], 1.0);
				sum += opt->sim_pi[k];
			}
			for (k = 0; k < opt->sim_K; ++k)
				opt->sim_pi[k] /= sum;
			debug_msg(QUIET <= fxn_debug, opt->quiet, "Simulated "
									"pi: ");
			if (QUIET <= opt->quiet)
				fprint_doubles(stderr, opt->sim_pi, opt->sim_K,
					3, 1);
		}

		for (unsigned int k = 0; k < opt->true_K; ++k)
			opt->true_cluster_size[k] = 0;

		/* simulate data */
		for (unsigned int i = 0; i < dat->n_observations; ++i) {

			/* choose cluster */
			double r = (double) rand() / RAND_MAX;
			for (k = 0, dsum = opt->sim_pi[0];
				dsum < r; dsum += opt->sim_pi[++k]);
			opt->sim_cluster[i] = k;
			++opt->true_cluster_size[k];
			for (unsigned int j = 0; j < dat->n_coordinates; ++j) {
				data_t c_ancestor = opt->sim_modes[k][j];

				/* choose coordinate */
				r = (double) rand() / RAND_MAX;
				data_t l = 0;
				for (dsum = Pt[c_ancestor*opt->sim_n_categories];
					dsum < r; dsum += Pt[c_ancestor
						* opt->sim_n_categories + ++l]);
				dat->data[i*dat->n_coordinates + j] = l;

			}
		}

		true_k = 0;
		for (unsigned int k = 0; k < opt->true_K; ++k)
			if (opt->true_cluster_size[k]) ++true_k;
	} while (opt->require_sim_K && true_k < opt->true_K);

	if (true_k < opt->true_K) {
		mmessage(WARNING_MSG, NO_ERROR, "Simulated %u distinct modes, "
			"but only %u produced data (not implemented).\n",
			opt->true_K, true_k);
		goto ABORT_SIMULATE_DATA;
	}

	FILE *fp = NULL;
	if (opt->datafile) {
		fp = fopen(opt->datafile, "w");
		if (!fp) {
			err = mmessage(ERROR_MSG, FILE_OPEN_ERROR,
				opt->datafile);
			goto ABORT_SIMULATE_DATA;
		}

		for (unsigned int i = 0; i < dat->n_observations; ++i) {
			/* first column is true assignments */
			fprintf(fp, "%u", opt->sim_cluster[i]);
			fprint_data_ts(fp, &dat->data[i*dat->n_coordinates],
				dat->n_coordinates, 0, 1);
		}
		fclose(fp);
	}


ABORT_SIMULATE_DATA:
	if (Pt) free(Pt);
	if (a) free(a);

	return err;
} /* simulate_data */

int finish_make_data(data *dat, options *opt) {
	int fxn_debug = ABSOLUTE_SILENCE;

	dat->dmat = malloc(dat->n_observations * sizeof *dat->dmat);
	data_t *rptr = dat->data;
	for (unsigned int i = 0; i < dat->n_observations; i++) {
		dat->dmat[i] = rptr;
		rptr += dat->n_coordinates;
	}
	if (fxn_debug)
		mmessage(DEBUG_MSG, NO_ERROR, "Allocated %dx%d data matrix\n",
			dat->n_observations, dat->n_coordinates);

	if (opt->data_outfile) {
		FILE *fp = fopen(opt->data_outfile, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
				opt->data_outfile);
		for (unsigned int i = 0; i < dat->n_observations; ++i) {
			if (opt->true_column < UINT_MAX)
				fprintf(fp, "%u", opt->true_cluster[i]);
			for (unsigned int j = 0; j < dat->n_coordinates; ++j)
				if (opt->true_column < UINT_MAX || j)
					fprintf(fp, " %u", dat->dmat[i][j]);
				else
					fprintf(fp, "%u", dat->dmat[i][j]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		if (fxn_debug)
			mmessage(DEBUG_MSG, NO_ERROR, "Wrote data to file"
				"'%s'\n", opt->data_outfile);
	}

	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations
		* sizeof *dat->best_cluster_id);

	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"data::cluster_id");

	if (opt->shuffle) {
		dat->obsn_idx = malloc(dat->n_observations
			* sizeof *dat->obsn_idx);
		dat->best_obsn_idx = malloc(dat->n_observations
			* sizeof *dat->obsn_idx);
		if (!dat->obsn_idx || !dat->best_obsn_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data:obsn_idx");
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			dat->obsn_idx[i] = i;
	}

	if (!opt->K)
		return NO_ERROR;

	if (!opt->seed_idx) {
		dat->seed_idx = malloc(opt->K * sizeof *dat->seed_idx);
		if (!dat->seed_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::seeds_idx");
	} else
		dat->seed_idx = opt->seed_idx;

	dat->best_seed_idx = malloc(opt->K * sizeof *dat->best_seed_idx);
	dat->seeds = malloc(opt->K * sizeof *dat->seeds);
	dat->best_modes = malloc(opt->K * sizeof *dat->best_modes);

	if (dat->best_seed_idx == NULL || dat->seeds == NULL
		|| dat->best_modes == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");

	if (opt->n_inner_init > 1) {
		dat->ini_seeds = malloc(opt->K * sizeof *dat->ini_seeds);
		dat->ini_seed_idx = malloc(opt->K * sizeof *dat->ini_seed_idx);
		if (!dat->ini_seeds || !dat->ini_seed_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::ini_*");
	}

	data_t *tmp1 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp2 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp3 = NULL;
	if (!tmp1 || !tmp2)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");
	if (opt->n_inner_init > 1) {
		tmp3 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
		if (!tmp3)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::ini_seeds");
	}
	for (unsigned int i = 0; i < opt->K; i++) {
		dat->seeds[i] = tmp1;
		dat->best_modes[i] = tmp2;
		tmp1 += dat->n_coordinates;
		tmp2 += dat->n_coordinates;
		if (opt->n_inner_init > 1) {
			dat->ini_seeds[i] = tmp3;
			tmp3 += dat->n_coordinates;
		}
	}

	if (opt->pfile) {
		FILE *fp = fopen(opt->pfile, "r");
		for (unsigned int i = 0; i < dat->n_observations; ++i) {
			if (fscanf(fp, "%u", &dat->cluster_id[i]) != 1) {
				fclose(fp);
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->pfile);
			}
			if (dat->cluster_id[i] >= opt->K) {
				fclose(fp);
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"Partition file assigns to more than %u"
					" clusters.", opt->K);
			}
		}
		fclose(fp);

	}

	if (opt->sfile) {
		FILE *fp = fopen(opt->sfile, "r");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->sfile);

		opt->n_seed_set = 0;
		do {
			char c = fgetc(fp);

			/* new line with content */
			if (c != '\n' && !feof(fp))
				++opt->n_seed_set;

			/* fast-forward through line */
			while (c != '\n' && !feof(fp)) c = fgetc(fp);
		} while (!feof(fp));

		data_t **seeds = dat->seeds;
		/* allocate space for seed set */
		if (opt->n_seed_set > opt->K) {
			data_t *tmp = malloc(opt->n_seed_set
				* dat->n_coordinates * sizeof *tmp);
			opt->seed_set = malloc(opt->n_seed_set * sizeof
				*opt->seed_set);
			if (!tmp || !opt->seed_set)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"options:seed_set");
			for (unsigned int k = 0; k < opt->n_seed_set; ++k) {
				opt->seed_set[k] = tmp;
				tmp += dat->n_coordinates;
			}
			seeds = opt->seed_set;
			opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
		} else if (opt->n_seed_set < opt->K)
			mmessage(WARNING_MSG, NO_ERROR,
				"Requesting %u clusters, but only %u seeds in "
				"seed file '%s'.  Will generate remaining seeds"
				" with chosen initialization method.\n", opt->K,
				opt->n_seed_set, opt->sfile);
		else {	/* exactly options::K seeds provided */
			opt->init_method = KMODES_INIT_USER_SEEDS;
			if (opt->n_init > 1)
				mmessage(WARNING_MSG, INVALID_USER_INPUT,
					"Resetting to one initialization.\n");
			opt->n_init = 1;
			opt->n_inner_init = 1;
		}

fprintf(stderr, "Found %u seeds\n", opt->n_seed_set);
		rewind(fp);
		for (unsigned int k = 0; k < opt->n_seed_set; ++k) {
			if (fscan_data_ts(fp, seeds[k], dat->n_coordinates)) {
				fclose(fp);
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					opt->sfile);
			}
			if (opt->subtract_one)
				for (unsigned int j = 0; j < dat->n_coordinates;
					++j) {
					if (seeds[k][j] == 0)
						return mmessage(ERROR_MSG,
							INVALID_USER_INPUT,
							"Command option -1 "
							"requested, but -i "
							"<istr> has 0-based "
							"data.\n");
					else
						seeds[k][j] -= 1;
				}
			if (dat->ini_seeds && opt->n_seed_set  < opt->K)
				memcpy(dat->ini_seeds[k], seeds[k],
					dat->n_coordinates * sizeof **seeds);
		}
		fclose(fp);
	}

	/*
	for (unsigned int i = 0; i < dat->n_observations; ++i) {
		unsigned int n_identical = 0;
		for (unsigned int j = 0; j < dat->n_observations; ++j) {
			unsigned int l = 0;
			for (; l < dat->n_coordinates; ++l)
				if (dat->data[i*dat->n_observations + l] != dat->data[j*dat->n_observations + l]) {
					break;
				}
			if (l == dat->n_coordinates)
				n_identical++;
		}
		if (n_identical > dat->n_observations - 7)
			fprintf(stderr, "Read %u identical to %u (of %u [%u]) others\n", i, n_identical, dat->n_observations, dat->n_observations - 7);
	}
	exit(0);
	*/


	return NO_ERROR;
} /* finish_make_data */

/**
 * Shuffle data input order.  Since the k-modes algorithms take the observations
 * in input order for initialization, and many of the real datasets are ordered
 * by class, it is important to randomize this aspect of the data.
 *
 * See https://stackoverflow.com/questions/3343797/is-this-c-implementation-of-fisher-yates-shuffle-correct
 * and https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 *
 * @param dat	data pointer
 * @return	error status
 */
int shuffle_data(data *dat, options *opt)
{
	if (dat->n_observations > RAND_MAX)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Data set too big!");
	long n = (long) dat->n_observations;
	long i, j, lim;
	data_t *dptr;

	for (i = n - 1; i > 0; --i) {
		lim = RAND_MAX - RAND_MAX % (i + 1);
		do {
			j = rand();
		} while (j >= lim);
		j = j % (i + 1);

		dptr = dat->dmat[j];
		dat->dmat[j] = dat->dmat[i];
		dat->dmat[i] = dptr;

		if (dat->obsn_idx) {
			unsigned int ui = dat->obsn_idx[j];
			dat->obsn_idx[j] = dat->obsn_idx[i];
			dat->obsn_idx[i] = ui;
		}

		/* options:true_cluster and options:sim_cluster same */
		if (opt->sim_cluster) {
			unsigned int ui = opt->sim_cluster[j];
			opt->sim_cluster[j] = opt->sim_cluster[i];
			opt->sim_cluster[i] = ui;
		} else if (opt->true_cluster) {
			unsigned int ui = opt->true_cluster[j];
			opt->true_cluster[j] = opt->true_cluster[i];
			opt->true_cluster[i] = ui;
		}
	}

	return NO_ERROR;
} /* shuffle_data */


/**
 * Write best solution
 *
 * @param dat	data object pointer
 * @param opt	options object pointer
 * @param fps	solution file pointer
 * @param fpi	initialization file pointer
 */
void write_best_solution(data *dat, options *opt, FILE *fps)
{
	/* report best solution to stdout and outfile */
	if (dat->n_init) {
		if (opt->seconds > 0) {
			if (opt->quiet >= QUIET)
				fprintf(stdout, "Cost scaling factor: %f\n",
					dat->first_cost);
			else
				fprintf(fps, "Cost scaling factor: %f\n",
					dat->first_cost);
		}
		if (opt->quiet >= QUIET)
			fprintf(stdout, "Maximum cost: %.0f\n",
				dat->worst_cost);
		else
			fprintf(fps, "Maximum cost: %.0f\n", dat->worst_cost);
		if (opt->quiet >= QUIET) {
			fprintf(stdout, "Best optimized criterion: %.0f",
				dat->best_total);
			fprint_doubles(stdout, dat->best_criterion, opt->K, 0,
				1);
		}
		if (fps) {
			fprintf(fps, "Best optimized criterion: %.0f\n",
				dat->best_total);
			fprint_doubles(fps, dat->best_criterion, opt->K, 0, 1);
		}
		if (opt->simulate || opt->true_cluster) {
			if (opt->quiet >= QUIET)
				fprintf(stdout, "Maximum AR: %f\n", dat->best_rand);
			if (fps)
				fprintf(fps, "Maximum AR: %f\n", dat->best_rand);
		}
		if (opt->quiet >= QUIET && opt->K > 0) {
			fprintf(stdout, "Best cluster sizes:");
			fprint_uints(stdout, dat->best_cluster_size, opt->K,
				(int)(log10(dat->n_observations) + 1), 1);
		}
		if (fps && opt->K > 0) {
			fprintf(fps, "Best cluster sizes:");
			fprint_uints(fps, dat->best_cluster_size, opt->K, 0, 1);
		}
		if (opt->quiet >= QUIET && opt->K > 0 &&
			opt->init_method < KMODES_INIT_NUMBER_RANDOM_METHODS) {

			fprintf(stdout, "Best solution originating seeds:");
			fprint_uints(stdout, dat->best_seed_idx, opt->K,
				(int)(log10(dat->n_observations) + 1), 1);
		}
		if (fps && opt->K > 0 && opt->init_method
			< KMODES_INIT_NUMBER_RANDOM_METHODS) {

			fprintf(fps, "Best solution originating seeds:");
			fprint_uints(fps, dat->best_seed_idx, opt->K, 0, 1);
		}
		if (fps && opt->K > 0) {
			fprintf(fps, "Best solution cluster assignments:\n");
			fprint_uints(fps, dat->best_cluster_id,
				dat->n_observations, 0, 1);
			fprintf(fps, "\n");
		}
		if (fps && opt->K > 0 && opt->shuffle) {
			fprintf(fps, "Best solution indexing:\n");
			fprint_uints(fps, dat->best_obsn_idx,
				dat->n_observations, 0, 1);
			fprintf(fps, "\n");
		}
		if (opt->quiet >= QUIET)
			fprintf(stdout, "Best modes:\n");
		for (unsigned int k = 0; k < opt->K; ++k) {
			for (unsigned int j = 0; j < dat->n_coordinates; ++j) {
				if (opt->subtract_one)
					++dat->best_modes[k][j];
				if (opt->quiet >= QUIET)
					fprintf(stdout, " %*"
						AMPLICLUST_PRIu_data_t, (int)
						(log10(dat->n_categories[j] -
						(dat->n_categories[j]>1)) + 1),
						(unsigned int)
						dat->best_modes[k][j]);
			}
			if (opt->quiet >= QUIET)
				fprintf(stdout, "\n");
		}
		if (fps) {
			fprintf(fps, "Best modes:\n");
			for (unsigned int k = 0; k < opt->K; ++k)
				fprint_data_ts(fps, dat->best_modes[k],
					dat->n_coordinates, 0, 1);
		}
		if (opt->true_modes) {
			if (opt->quiet >= QUIET)
				fprintf(stdout, "Distance matrix between true "
					"and best modes:\n");
			if (fps)
				fprintf(fps, "Distance matrix between true and "
					"best modes:\n");
			for (unsigned int k = 0; k < opt->true_K; ++k) {
				for (unsigned int i = 0; i < opt->K; ++i) {
					unsigned int cnt = 0;
					for (unsigned int j = 0; j
						< dat->n_coordinates; ++j)
						cnt += (dat->best_modes[i][j]
							!= opt->true_modes[k][j]);
					if (opt->quiet >= QUIET)
						fprintf(stdout, " %*u", (int)
							ceil(log10(dat->n_coordinates + 1)),
							cnt);
					if (fps)
						fprintf(fps, " %u", cnt);
				}
				if (opt->quiet >= QUIET)
					fprintf(stdout, "\n");
				if (fps)
					fprintf(fps, "\n");
			}
		}
	}
	if (dat->ntimes > 1) {
		dat->avg_time /= dat->ntimes;
		dat->sd_time = sqrt((dat->sd_time / dat->ntimes
			- dat->avg_time * dat->avg_time) / dat->ntimes);
		if (opt->quiet >= QUIET)
			fprintf(stdout, "Time to better: %f +/- %f; times=%u\n",
				dat->avg_time, dat->sd_time, dat->ntimes);
		if (fps)
			fprintf(fps, "Time to better: %f +/- %f; times=%u\n",
				dat->avg_time, dat->sd_time, dat->ntimes);
	}

	if (dat->n_init) {
		dat->avg_iter /= dat->n_init;
		dat->sd_iter = sqrt((dat->sd_iter / dat->n_init
			- dat->avg_iter*dat->avg_iter)
			/ dat->n_init);
		dat->avg_cost /= dat->n_init;
		dat->sd_cost = sqrt((dat->sd_cost / dat->n_init
			- dat->avg_cost*dat->avg_cost)
			/ dat->n_init);
		/* scaled by data::n_init or first_cost to avoid overflow */
		if (opt->seconds > 0) {
			dat->avg_cost *= dat->first_cost;
			dat->sd_cost *= dat->first_cost;
		} else {
			dat->avg_cost *= dat->n_init;
			dat->sd_cost *= dat->n_init;
		}
		if (opt->quiet > ABSOLUTE_SILENCE)
			fprintf(stdout, "%.0f %f %f %f %f", dat->best_total,
				dat->avg_iter, dat->sd_iter, dat->avg_cost,
				dat->sd_cost);
		if (fps)
			fprintf(fps, "%.0f %f %f %f %f", dat->best_total,
				dat->avg_iter, dat->sd_iter, dat->avg_cost,
				dat->sd_cost);
		if (opt->simulate || opt->true_cluster) {
			dat->avg_ar /= dat->n_init;
			dat->sd_ar = sqrt((dat->sd_ar / dat->n_init
				- dat->avg_ar * dat->avg_ar) / dat->n_init);
			dat->avg_mi /= dat->n_init;
			dat->sd_mi = sqrt((dat->sd_mi / dat->n_init
				- dat->avg_mi * dat->avg_mi) / dat->n_init);
			dat->avg_vi /= dat->n_init;
			dat->sd_vi = sqrt((dat->sd_vi / dat->n_init
				- dat->avg_vi * dat->avg_vi) / dat->n_init);
			if (opt->quiet > ABSOLUTE_SILENCE)
				fprintf(stdout, " %f %f %f %f %f %f",
					dat->avg_ar, dat->sd_ar, dat->avg_mi,
					dat->sd_mi, dat->avg_vi, dat->sd_vi);
			if (fps)
				fprintf(fps, " %f %f %f %f %f %f", dat->avg_ar,
					dat->sd_ar, dat->avg_mi, dat->sd_mi,
					dat->avg_vi, dat->sd_vi);
		}
		if (opt->quiet > ABSOLUTE_SILENCE)
			fprintf(stdout, " %f", dat->seconds);
		if (fps)
			fprintf(fps, " %f", dat->seconds);
		if (opt->quiet > ABSOLUTE_SILENCE)
			fprintf(stdout, " %u\n", dat->n_init);
		if (fps)
			fprintf(fps, " %u\n", dat->n_init);
	}
} /* write_best_solution */

/**
 * Estimate K.
 *
 * @param dat	data object pointer
 * @param opt	options object pointer
 * @return	error status
 */
int estimate_k(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int err = NO_ERROR;
	double Y = opt->n_effective_coordinates / 2.;
	double cost, rrcost, krcost;	/* rescaled costs (or distortions) */
	double sd, rsd, ksd;
	double pcost = 0., prrcost = 0., pkrcost = 0.;	/* previous r-costs */
	double jump, rjump, kjump;
	double pjump = 0, prjump = 0, pkjump = 0, pKL = 0, prKL = 0, pkKL = 0;
	double max_jump = 0, max_rjump = 0, max_kjump = 0, max_KL = 0, max_rKL = 0, max_kKL = 0;
	double max_jump2 = 0, max_rjump2 = 0, max_kjump2 = 0;
	unsigned int max_K_jump = 0, max_K_rjump = 0, max_K_kjump = 0, max_K_KL = 0, max_K_rKL = 0, max_K_kKL = 0;
	unsigned int max_K_jump2 = 0, max_K_rjump2 = 0, max_K_kjump2 = 0;
	unsigned int peak_K_jump = 0, peak_K_rjump = 0, peak_K_kjump = 0, peak_K_KL = 0, peak_K_rKL = 0, peak_K_kKL = 0;
	double Ck, pCk, rCk, kCk, prCk, pkCk, KL, rKL, kKL;
	double max_cost = 0;
	int too_small = 0;

	if (!opt->result_files)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Use -k<kuint1> "
			"... arguments to provides output files.\n");

	double *asize = malloc(opt->max_k * sizeof *asize);
	double *pasize = malloc(opt->max_k * sizeof *pasize);
	double *var = malloc(opt->max_k * sizeof *var);
	double *pcrit = malloc(opt->max_k * sizeof *pcrit);
	unsigned int *ridx = malloc(dat->n_observations * sizeof *ridx);
	double *obsn_hd = malloc(dat->n_observations * sizeof *obsn_hd);
	unsigned int **obsn_cnt = malloc(dat->n_observations * sizeof *obsn_cnt);
	unsigned int *obsn_cnt_sum = malloc(dat->n_observations * sizeof *obsn_cnt_sum);
	double *pobsn_hd = malloc(dat->n_observations * sizeof *pobsn_hd);
	unsigned int **pobsn_cnt = malloc(dat->n_observations * sizeof *pobsn_cnt);
	unsigned int *pobsn_cnt_sum = malloc(dat->n_observations * sizeof *pobsn_cnt_sum);
	if (asize == NULL || var == NULL || pasize == NULL
		|| pcrit == NULL || ridx == NULL
		|| obsn_hd == NULL || obsn_cnt == NULL || obsn_cnt_sum == NULL
		|| pobsn_hd == NULL || pobsn_cnt == NULL || pobsn_cnt_sum == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "rewrite this");
	unsigned int *tmp1 = malloc(dat->n_observations * opt->max_k * sizeof *tmp1);
	unsigned int *tmp2 = malloc(dat->n_observations * opt->max_k * sizeof *tmp2);
	if (tmp1 == NULL || tmp2 == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "rewrite this");
	for (unsigned int i = 0; i < dat->n_observations; ++i) {
		obsn_cnt[i] = tmp1;
		pobsn_cnt[i] = tmp2;
		tmp1 += opt->max_k;
		tmp2 += opt->max_k;
	}

	opt->n_bootstrap = 1000;
	for (unsigned int l = 0; l < opt->n_k; ++l) {
		opt->K = opt->min_k + l;
		allocate_data_for_k(dat, opt->K);

		/* find best solution among all files with this k */
		unsigned best_j = 0;
		double best_total = INFINITY;
		for (unsigned int j = 0; j < opt->n_result_files[l]; ++j) {
			opt->soln_file = opt->result_files[l][j];
			restore_state(dat, opt);
			if (best_total > dat->best_total) {
				best_j = j;
				best_total = dat->best_total;
			}
		}

		/* load best solution */
		opt->soln_file = opt->result_files[l][best_j];
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Best solution is "
			"%.0f from file '%s'.\n", best_total, opt->soln_file);
		restore_state(dat, opt);
		max_cost = dat->best_total;

		/* set up reverse index lookup */
		if (opt->K > 1 && dat->best_obsn_idx)
			for (unsigned int i = 0; i < dat->n_observations; ++i)
				ridx[dat->best_obsn_idx[i]] = i;
		else
			for (unsigned int i = 0; i < dat->n_observations; ++i)
				ridx[i] = i;

		/* initialize */
		for (unsigned int j = 0; j < opt->K; ++j) {
			dat->criterion[j] = 0.;
			var[j] = 0.;
			asize[j] = 0.;
			if (dat->best_cluster_size[j] == 1) {
				mmessage(ERROR_MSG, INTERNAL_ERROR, "Cluster %u"
					" has just 1 member (K=%u).\n", j,
					opt->K);
				too_small = 1;
				break;
			}
		}

		if (too_small)
			break;

		for (unsigned int i = 0; i < dat->n_observations; ++i) {
			unsigned int j = dat->best_cluster_id[ridx[i]];
			debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Observation "
						"%u assigned to %u\n", i, j);
			obsn_hd[i] = hd(dat->dmat[i], dat->best_modes[j],
				dat->n_coordinates);
			//fprintf(stderr, "Observation %u in cluster %u/%u, distance %.0f\n", i, j, opt->K, obsn_hd[i]);

			/* allocate observations equally to equidistant modes */
			obsn_cnt[i][j] = 1;
			obsn_cnt_sum[i] = 1;
			for (unsigned m = 0; m < opt->K; ++m) {
				if (m == j) continue;
				obsn_cnt[i][m] = 0;
				if (hd(dat->dmat[i], dat->best_modes[m],
					dat->n_coordinates) == obsn_hd[i]) {
					obsn_cnt[i][m] = 1;
					++obsn_cnt_sum[i];
				}
			}
			double tmp1 = obsn_hd[i] / obsn_cnt_sum[i];
			double tmp2 = 1. / obsn_cnt_sum[i];
			for (unsigned k = 0; k < opt->K; ++k)
				if (obsn_cnt[i][k]) {
					asize[k] += tmp2;
					dat->criterion[k] += tmp1;
					var[k] += tmp1 * tmp1;
				}
		}

		compute_costs(dat->criterion, var, asize, dat->n_observations, opt->K, &cost, &rrcost, &krcost, &sd, &rsd, &ksd);
		if (!l && opt->K > 1) {
			pcost = cost;
			prrcost = rrcost;
			pkrcost = krcost;
			continue;
		}

		/* we can compute jump statistics: these are the observed */
		compute_jump_stats(cost, rrcost, krcost, pcost, prrcost, pkrcost, &jump, &rjump, &kjump, Y, !l);

		Ck = pow(opt->K - 1, 1/Y) * pcost - pow(opt->K, 1/Y) * cost;
		rCk = pow(opt->K - 1, 1/Y) * prrcost - pow(opt->K, 1/Y) * rrcost;
		kCk = pow(opt->K - 1, 1/Y) * pkrcost - pow(opt->K, 1/Y) * krcost;
		if (l > 1) {
			KL = fabs(pCk/Ck);
			rKL = fabs(prCk/rCk);
			kKL = fabs(pkCk/kCk);
		} else
			KL = rKL = kKL = 0.;

		/* bootstrap hd */
		double jmu = 0, rjmu = 0, kjmu = 0;
		double jvar = 0, rjvar = 0, kjvar = 0;
		for (unsigned i = 0; i < opt->n_bootstrap; ++i) {
			debug_msg(DEBUG_I <= fxn_debug, fxn_debug,
							"BOOTSTRAP %u\n", i);
			double cost_bs, rrcost_bs, krcost_bs;
			double pcost_bs = 0., prrcost_bs = 0., pkrcost_bs = 0.;
			for (unsigned int k = 0; k < opt->K; ++k) {
				dat->criterion[k] = 0.;
				pcrit[k] = 0.;
				var[k] = 0.;
				asize[k] = 0.;
				pasize[k] = 0.;
			}
			for (unsigned int j = 0; j < dat->n_observations; ++j) {
				long lim = RAND_MAX - RAND_MAX % dat->n_observations;
				long rnd;
				do {
					rnd = rand();
				} while (rnd >= lim);
				rnd = rnd % dat->n_observations;
				double tmp1 = obsn_hd[rnd] / obsn_cnt_sum[rnd];
				double tmp2 = 1. / obsn_cnt_sum[rnd];
				for (unsigned k = 0; k < opt->K; ++k)
					if (obsn_cnt[rnd][k]) {
						asize[k] += tmp2;
						dat->criterion[k] += tmp1;
					}
				if (l) {
					tmp1 = pobsn_hd[rnd] / pobsn_cnt_sum[rnd];
					tmp2 = 1. / pobsn_cnt_sum[rnd];
					for (unsigned int k = 0; k < opt->K - 1; ++k)
						if (pobsn_cnt[rnd][k]) {
							pasize[k] += tmp2;
							pcrit[k] += tmp1;
						}
				}
			}

			/* compute previous costs */
			if (l)
				compute_costs(pcrit, NULL, pasize, dat->n_observations, opt->K - 1, &pcost_bs, &prrcost_bs, &pkrcost_bs, NULL, NULL, NULL);

			/* compute current costs and jump statistics */
			compute_costs(dat->criterion, NULL, asize, dat->n_observations, opt->K, &cost_bs, &rrcost_bs, &krcost_bs, NULL, NULL, NULL);

			double jump_bs, rrjump_bs, krjump_bs;
			compute_jump_stats(cost_bs, rrcost_bs, krcost_bs, pcost_bs, prrcost_bs, pkrcost_bs, &jump_bs, &rrjump_bs, &krjump_bs, Y, 0);
			jmu += jump_bs;
			rjmu += rrjump_bs;
			kjmu += krjump_bs;
			jvar += jump_bs * jump_bs;
			rjvar += rrjump_bs * rrjump_bs;
			kjvar += krjump_bs * krjump_bs;
		}
		jmu /= opt->n_bootstrap;
		rjmu /= opt->n_bootstrap;
		kjmu /= opt->n_bootstrap;
		jvar = (jvar/opt->n_bootstrap - jmu * jmu);// / opt->n_bootstrap;
		rjvar = (rjvar/opt->n_bootstrap - rjmu * rjmu);// / opt->n_bootstrap;
		kjvar = (kjvar/opt->n_bootstrap - kjmu * kjmu);// / opt->n_bootstrap;
		if (jump > max_jump) {
			max_K_jump = opt->K;
			max_jump = jump;
		}
		if (jump > max_jump2 && jmu - 1.96*sqrt(jvar) > 0) {
			max_K_jump2 = opt->K;
			max_jump2 = jump;
		}
		if (rjump > max_rjump) {
			max_K_rjump = opt->K;
			max_rjump = rjump;
		}
		if (rjump > max_rjump2 && rjmu - 1.96*sqrt(rjvar) > 0) {
			max_K_rjump2 = opt->K;
			max_rjump2 = rjump;
		}
		if (kjump > max_kjump) {
			max_K_kjump = opt->K;
			max_kjump = kjump;
		}
		if (kjump > max_kjump2 && kjmu - 1.96*sqrt(kjvar) > 0) {
			max_K_kjump2 = opt->K;
			max_kjump2 = kjump;
		}
		if (KL > max_KL) {
			max_K_KL = opt->K - 1;
			max_KL = KL;
		}
		if (rKL > max_rKL) {
			max_K_rKL = opt->K - 1;
			max_rKL = rKL;
		}
		if (kKL > max_kKL) {
			max_K_kKL = opt->K - 1;
			max_kKL = kKL;
		}
		if (l && opt->K > 1) {
			if (!peak_K_jump && jump < pjump)
				peak_K_jump = opt->K - 1;
			if (!peak_K_rjump && rjump < prjump)
				peak_K_rjump = opt->K - 1;
			if (!peak_K_kjump && kjump < pkjump)
				peak_K_kjump = opt->K - 1;
			if (!peak_K_KL && KL < pKL)
				peak_K_KL = opt->K - 2;
			if (!peak_K_rKL && rKL < prKL)
				peak_K_rKL = opt->K - 2;
			if (!peak_K_kKL && kKL < pkKL)
				peak_K_kKL = opt->K - 2;
		}

		fprintf(stdout, "%s d=%*.0f (%*.0f) rd=%*.0f (%*.0f) kd=%*.0f (%*.0f)\n", opt->datafile,
			(int)(log10(max_cost-1) + 1), cost, (int)log10(max_cost-1), sd, (int)(log10(max_cost-1) + 1), rrcost, (int)log10(max_cost-1), rsd, (int)(log10(max_cost-1) + 1), krcost, (int)log10(max_cost-1), ksd);
		fprintf(stdout, "%s J[%2u] = %9.3g %9.3g (%9.3g: %9.3g - %9.3g); rJ[%2u] = %9.3g %9.3g (%9.3g: %9.3g - %9.3g); kJ[%2u] = %9.3g %9.3g (%9.3g: %9.3g - %9.3g)\n",
			opt->datafile,
			opt->K, jump, jmu, sqrt(jvar), jmu - 1.96*sqrt(jvar), jmu + 1.96*sqrt(jvar),
			opt->K, rjump, rjmu, sqrt(rjvar), rjmu - 1.96*sqrt(rjvar), rjmu + 1.96*sqrt(rjvar),
			opt->K, kjump, kjmu, sqrt(kjvar), kjmu - 1.96*sqrt(kjvar), kjmu + 1.96*sqrt(kjvar));
		fprintf(stdout, "%s KL[%2u] = %9.3g; rKL[%2u] = %9.3g; kKL[%2u] = %9.3g\n",
			opt->datafile, opt->K - 1, KL, opt->K - 1, rKL, opt->K - 1, kKL);
		pcost = cost;
		prrcost = rrcost;
		pkrcost = krcost;
		pjump = jump;
		prjump = rjump;
		pkjump = kjump;
		pKL = KL;
		prKL = rKL;
		pkKL = kKL;
		pCk = Ck;
		prCk = rCk;
		pkCk = kCk;
		memcpy(pasize, asize, opt->K * sizeof(*asize));
		memcpy(pcrit, dat->criterion, opt->K * sizeof(*asize));
		memcpy(pobsn_cnt_sum, obsn_cnt_sum, dat->n_observations * sizeof *obsn_cnt_sum);
		memcpy(pobsn_hd, obsn_hd, dat->n_observations * sizeof *obsn_hd);
		for (unsigned int i = 0; i < dat->n_observations; ++i) {
			memcpy(pobsn_cnt[i], obsn_cnt[i], opt->K * sizeof *obsn_cnt[i]);
			for (unsigned int j = 0; j < opt->K; ++j)
				if (pobsn_cnt[i][j] != obsn_cnt[i][j]) exit(0);
		}
	}

	fprintf(stdout, "%s     Maxima: J = %2u, rJ = %2u, kJ = %2u, J2 = %2u, rJ2 = %2u, kJ2 = %2u, KL = %2u, rKL = %2u, kKL = %2u\n", opt->datafile,
		max_K_jump, max_K_rjump, max_K_kjump, max_K_jump2, max_K_rjump2, max_K_kjump2, max_K_KL, max_K_rKL, max_K_kKL);
	fprintf(stdout, "%s First peak: J = %2u, rJ = %2u, kJ = %2u, KL = %2u, rKL = %2u, kKL = %2u\n", opt->datafile,
		peak_K_jump, peak_K_rjump, peak_K_kjump, peak_K_KL, peak_K_rKL, peak_K_kKL);

	free(asize);
	free(pasize);
	free(pcrit);
	free(var);
	free(obsn_cnt_sum);
	free(pobsn_cnt_sum);
	free(obsn_hd);
	free(pobsn_hd);
	free(obsn_cnt[0]);
	free(obsn_cnt);
	free(pobsn_cnt[0]);
	free(pobsn_cnt);

	return err;
} /* estimate_k */

void compute_costs(double *crit, double *var, double *size, unsigned int n, unsigned int K,
	double *cost, double *rrcost, double *krcost, double *sd, double *rrsd, double *krsd)
{
	double sum = 0.;

	*rrcost = *cost = 0;
	if (var)
		*sd = *rrsd = 0;
	for (unsigned int j = 0; j < K; ++j) {
		if (!size[j])
			continue;
		*cost += crit[j];
		crit[j] /= size[j];
if (isnan(crit[j])) {
	printf("Exiting here... size=%g\n", size[j]);
	exit(0);
}
		*rrcost += crit[j];
		sum += 1. / size[j];
		if (var) {
			var[j] = (var[j]/size[j] - crit[j]*crit[j]);
			*rrsd += var[j] / size[j];
			*sd += var[j] * size[j];
		}
		//fprintf(stderr, "%f %f (%f)\n", *rrcost, *rrsd, size[j]);
	}
	*krcost = *rrcost * n / K;
	*rrcost *= K / sum;

	if (var) {
		*sd = sqrt(*sd);
		*krsd = sqrt(*rrsd) * n / K;
		*rrsd = K * sqrt(*rrsd) / sum;
	}
/*
	fprintf(stderr, "w[%u] = %.0f (%.1f); rrw[%u] = %.3f (%.3f); krw[%u] = %.3f (%.3f)\n",
		K, *cost, sqrt(*sd), K, *rrcost, *rrsd, K, *krcost, *krsd);
*/
} /* compute_costs */

void compute_jump_stats(double cost, double rrcost, double krcost,
	double pcost, double prrcost, double pkrcost,
	double *jump, double *rjump, double *kjump, double Y, int reset)
{
	static double factor, rfactor, kfactor;
	*jump = *rjump = *kjump = 0;
	if (reset) {
		factor = cost/2.;
		rfactor = rrcost/2.;
		kfactor = krcost/2.;
		fprintf(stdout, "factor=%f, rfactor=%f, kfactor=%f\n", factor, rfactor, kfactor);
	}
	*jump = pow(cost/factor, -Y) - (pcost > 0 ? pow(pcost/factor, -Y) : 0);
	*rjump = pow(rrcost/rfactor, -Y) - (prrcost > 0 ? pow(prrcost/rfactor, -Y) : 0);
	*kjump = pow(krcost/kfactor, -Y) - (pkrcost > 0 ? pow(pkrcost/kfactor, -Y) : 0);
	/*fprintf(stderr, "j = %0.3f; rj = %6.3f; kj = %6.3f\n", *jump, *rjump, *kjump);*/

} /* compute_jump_stats */

void free_data(data *dat) {
	if (dat) {
		if (dat->dmat) free(dat->dmat);
		if (dat->seeds) {
			free(dat->seeds[0]);
			free(dat->seeds);
		}
		if (dat->best_modes) {
			free(dat->best_modes[0]);
			free(dat->best_modes);
		}
		if (dat->seed_idx) free(dat->seed_idx);
		if (dat->best_seed_idx) free(dat->best_seed_idx);
		if (dat->cluster_id) free(dat->cluster_id);
		if (dat->best_cluster_id) free(dat->best_cluster_id);
		if (dat->criterion) free(dat->criterion);
		if (dat->best_criterion) free(dat->best_criterion);
		if (dat->cluster_size) free(dat->cluster_size);
		if (dat->best_cluster_size) free(dat->best_cluster_size);
		if (dat->data) free(dat->data);
		if (dat->n_categories) free(dat->n_categories);
		free(dat);
	}
} /* free_data */

void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i) fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - cluster observations with categorical predictors\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [-r <rulong> -n <nuint> -h97|-l|-w -i <istr>|<iuint1...k> -p <pfile> -o <ofile>] -k <kuint> -f <ffile> ...\n", &cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s clusters observations found in file <ffile> into <kuint> clusters.  It randomly initialize <iuint> times after setting random number seed <sulong>.\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-c <cint>\n\t\tSet the column containing the truth.\n");
	fprintf(fp, "\t--cont\n\t\tContinue previous run.\n");
	fprintf(fp, "\t-k <kuint>\n\t\tSet the desired number of clusters K.\n");
	fprintf(fp, "\t-k<kuint> <kstr1> <kstr2> ...\n\t\tThe names of output files for K=<kuint>, limited by POSIX ARG_MAX.\n");
	fprintf(fp, "\t-j\n\t\tEstimate K from multiple -k<kuint1> ... -k<kuint2> ... arguments.\n");
	fprintf(fp, "\t-l\n\t\tRun Lloyd's algorithm (cannot combine with -u).\n");
	fprintf(fp, "\t--hartigan\n\t\tUse hartigan updates with k-modes.\n");
	fprintf(fp, "\t-h97\n\t\tRun Huang's algorithm (combine with -u to replicate klaR).\n");
	fprintf(fp, "\t-w\n\t\tRun Hartigan and Wong algorithm (can combine with -u).\n");
	fprintf(fp, "\t-1\n\t\tSubtract 1 from the observation categories.\n");
	fprintf(fp, "\t-f <ffile> [<ffile2>]\n\t\tSet the input filename (with -s, simulate data are written to this file).\n");
	fprintf(fp, "\t\tIf <ffile2> specified, then write the data to this file: same as <ffile> to overwrite.\n");
	fprintf(fp, "\t-o <ofile> [<ofile2>]\n\t\tSet the output filename.  If second argument given, split information into first.\n");
	fprintf(fp, "\t-i <istr>\n\t\tSet initialization method (rnd|h97|h97rnd|hd17|clb09|clb09rnd|av07|av07grd|rndp|rnds).\n");
	fprintf(fp, "\t   <suint1> ... <suintk>\n\t\tSet the indices of the seeds.\n");
	fprintf(fp, "\t   <ifile>\n\t\tProvide file with possible seeds.\n");
	fprintf(fp, "\t\tIf there are more than <kuint> seeds in <ifile>, then method is 'rnds'.\n");
	fprintf(fp, "\t\tIf there are <kuint> seeds in <ifile>, then deterministic initialization with these seeds.\n");
	fprintf(fp, "\t\t'rndp' randomly selects seeds from given partitions (use with -c option).\n");
	fprintf(fp, "\t\t'rnds' randomly selects seeds from given seed set (repeat -i to provide <ifile> or just use -i <ifile> or -m <mfile>).\n");
	fprintf(fp, "\t-n <nuint>\n\t\tSet the desired number of initializations [OPTIONAL; DEFAULT: %u].\n", opt->n_init);
	fprintf(fp, "\t-pi <pflt1> ... <pfltk>\n\t\tThe cluster proportions in simulation or:\n");
	fprintf(fp, "\t-pi dir <pflt1> ... <pfltk>\n\t\tThe alpha for Dirichlet prior on pi.\n");
	fprintf(fp, "\t-p <pdbl>\n\t\tEffective number of independent coordinates, use with -k<kuint> ... arguments.\n");
	fprintf(fp, "\t-p <pfile>\n\t\tPartition file for initialization (overrides -i).\n");
	fprintf(fp, "\t-r <rulong>\n\t\tSet random number seed. [OPTIONAL]\n");
	fprintf(fp, "\t--run <ruint>\n\t\tNumber of inner initializations.\n");
	fprintf(fp, "\t-m <mdbl>|<mfile> [<mfile2>]\n\t\tTarget minimum or mode file.\n");
	fprintf(fp, "\t\tIf <mfile2> specified with <mfile>, then write modified mode information to file <mfile2>.\n");
	fprintf(fp, "\t--shuffle\n\t\tShuffle the data input order on each initialization.\n");
	fprintf(fp, "\t-s <sn> <sp> <sc> <st1> <st2>\n\t\tSimulation size (observations by coordinates by categories) and times.\n");
	fprintf(fp, "\t\t<st1> is the time separating the centers.\n");
	fprintf(fp, "\t\t<st2> is the time separating the observations from their centers.\n");
	fprintf(fp, "\t-t <tsec>\n\t\tNumber of seconds to run initializations.\n");
	fprintf(fp, "\t-u\n\t\tUse mode updates (combine with -m to replicate klaR).\n");
#ifndef __KMODES_NO_QTRANS__
	fprintf(fp, "\t-qt\n\t\tTurn off Hartigan & Wong quick-transfer stage.\n");
#endif
	fprintf(fp, "\t-q\n\t\tQuiet\n");
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\nUse <ofile2> (-o) and <mfile2> (-m) to convert from format output by fqmorph to run_kmodes format, where categories use contiguous categories 0, 1, 2, ..., without skipping.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i) fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */


static inline double hd(data_t *x, data_t *y, unsigned int p)
{
	double d = 0;
	for (unsigned int j = 0; j < p; ++j)
		d += (x[j] != y[j]);
	return d;
} /* hd */
