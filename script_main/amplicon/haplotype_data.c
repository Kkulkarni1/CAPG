/**
 * @file haplotype_data.c
 * @author Yudi Zhang
 */

#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>

#include "haplotype_data.h"
#include "array.h"
#include "kmodes.h"
#include "cluster_lloyds.h"
#include "error.h"
#include "util.h"
#include "math.h"


/**
 * Setup data object.
 */
int make_data(data **dat)
{
	*dat = malloc(sizeof **dat);
	
	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");
	
	(*dat)->fdata = NULL;
	(*dat)->dmat = NULL;
	(*dat)->qmat = NULL;
	(*dat)->n_quality = 0;
	(*dat)->prob_f = NULL;
	(*dat)->prob_t = NULL;
	(*dat)->data = NULL;
	(*dat)->v_ik = NULL;
	(*dat)->e_kjN = NULL;
	(*dat)->last_cent = NULL;
	(*dat)->last_assign = NULL;
	(*dat)->last_cost = NULL;
	(*dat)->last_vik = NULL;
	(*dat)->last_ekjn = NULL;
	(*dat)->ic2 = NULL;
	(*dat)->live = NULL;
	(*dat)->index = NULL;
	
	(*dat)->read_idx = NULL;
	(*dat)->max_read_length = 0;
	(*dat)->min_read_length = 0;
	(*dat)->lengths = NULL;
	
	(*dat)->true_cluster_id = NULL;
	(*dat)->true_cluster_size = NULL;
	
	(*dat)->categories = NULL;
	(*dat)->n_observations = 0;
	(*dat)->n_coordinates = 0;
	(*dat)->n_categories = NULL;
	(*dat)->max_n_categories = 0;
	(*dat)->n_categories = NULL;
	(*dat)->tot_n_categories = 0;
	
	(*dat)->seeds = NULL;
	(*dat)->ini_seeds = NULL;
	(*dat)->seed_idx = NULL;
	(*dat)->best_seed_idx = NULL;
	(*dat)->ini_seed_idx = NULL;
	
	(*dat)->total = 0;
	(*dat)->cluster_id = NULL;
	(*dat)->obsn_idx = NULL;
	(*dat)->best_obsn_idx = NULL;
	(*dat)->best_cluster_id = NULL;
	(*dat)->best_modes = NULL;
	(*dat)->criterion = NULL;
	(*dat)->best_criterion = NULL;
	(*dat)->cluster_size = NULL;
	(*dat)->best_cluster_size = NULL;
	(*dat)->best_seed_idx = NULL;
	
	(*dat)->best_total = -INFINITY;
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
	
	
	return NO_ERROR;
}/* make_data */

/**
 * Sample n from N objects without replacemnt.
 *
 * @param N	number of objects to sample from
 * @param n	number of objects to sample
 * @param idx	0-based indices of n selected objects from set {0,1,...,N-1}
 */
void sample_better(unsigned int N, unsigned int n, unsigned int *idx)
{
	unsigned int t = 0, m = 0, lim, i;
	
	while (m < n) {
		/* note integer division */
		lim = (N - t) * (RAND_MAX / (N - t));
		do {
			i = rand();
		} while (i >= lim);
		i %= (N - t);
		if (i < n - m)
			idx[m++] = t++;
		else
			++t;
	}
} /* sample_better */


/**
 * Sync state of data object.
 *
 * @param dat	pointer to data object
 * @param opt	pointer to option object
 * @return	error status
 */
int sync_state_data(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	
	// [TODO, BUG] Now temporarily shrink the data: 100 * 11
	dat->n_observations = dat->fdata->n_reads - 2900;
	dat->max_read_length = dat->fdata->n_max_length - 241;
	dat->min_read_length = dat->fdata->n_min_length - 241;

	// [TODO, KSD, BUG] Since you assume data::n_coordinates is nonzero
	// below, an inequality here must be a bug.  Raise and handle error.
	if (dat->max_read_length == dat->min_read_length)
		dat->n_coordinates = dat->min_read_length;
	
	dat->n_quality = dat->fdata->max_quality - dat->fdata->min_quality + 1;
	dat->lengths = malloc(dat->n_observations * sizeof *dat->lengths);
	
	if (!dat->lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.lengths");
	
	if (dat->fdata->n_lengths)
		memcpy(dat->lengths, dat->fdata->n_lengths, dat->n_observations
							* sizeof *dat->lengths);
	else
		for (size_t i = 0; i < dat->n_observations; ++i)
			dat->lengths[i] = dat->max_read_length;
	
	/* allocate the index array of reads */
	dat->read_idx = malloc(dat->fdata->n_reads * sizeof *dat->read_idx);
	
	if (!dat->read_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.read_idx");
	
	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		dat->read_idx[i] = i;
	
	/* make reads matrix */
	dat->dmat = malloc(dat->n_observations * sizeof *dat->dmat);
	
	if (!dat->dmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.dmat");
	
	unsigned char *rptr = dat->fdata->reads;
	for (size_t i = 0; i < dat->n_observations; ++i) {
		dat->dmat[i] = rptr;
		rptr += dat->lengths[i];
	}
	
	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) sequence matrix\n",
				dat->n_observations, dat->max_read_length);
	
	/* make quals matrix */
	/* allocate short-cut pointers to quality sequences */
	dat->qmat = malloc(dat->n_observations * sizeof *dat->qmat);
	
	if (!dat->qmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.qmat");
	
	unsigned char *qptr = dat->fdata->quals;
	for (size_t i = 0; i < dat->n_observations; i++) {
		dat->qmat[i] = qptr;
		qptr += dat->lengths[i];
	}
	
	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) quality matrix\n",
				dat->n_observations, dat->max_read_length);
	
	/* Find the number of category and count it */
	//dat->categories = malloc(dat->tot_n_categories * sizeof *(dat->categories));
	//if (!dat->categories)
	//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.cat");

	// [KSD, TODO] Move this earlier where you seamlessly handle variable
	// length reads.  Also get rid of data::max_read_length and
	// data::min_read_length if you do not want variable length reads.
	if (dat->fdata->n_lengths)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "k-haplotypes does "
					"not work with variable length reads!");
	
	dat->categories = find_unique(dat->fdata->reads, dat->n_observations
				* dat->max_read_length, &dat->tot_n_categories);
	
	/* Allocate for prob_t and prob_f and compute the values */
	MAKE_2ARRAY(dat->prob_f, dat->n_observations, dat->max_read_length);
	MAKE_2ARRAY(dat->prob_t, dat->n_observations, dat->max_read_length);
	
	compute_pij(dat, dat->prob_t, dat->prob_f);
	
	/* initialize the auxiliary variables */
	MAKE_2ARRAY(dat->v_ik, dat->n_observations, opt->K);
	SETZERO_MATRIX(dat->v_ik, dat->n_observations, opt->K);
	MAKE_3ARRAY(dat->e_kjN, opt->K, dat->n_coordinates, dat->tot_n_categories);
	SETZERO_3DARRAY(dat->e_kjN, opt->K, dat->n_coordinates, dat->tot_n_categories);
	
	MAKE_1ARRAY(dat->last_assign, dat->n_observations);
	MAKE_2ARRAY(dat->last_cent, opt->K, dat->n_coordinates);
	MAKE_1ARRAY(dat->last_cost, opt->K);
	MAKE_2ARRAY(dat->last_vik, dat->n_observations, opt->K);
	MAKE_3ARRAY(dat->last_ekjn, opt->K, dat->n_coordinates, dat->tot_n_categories);
	MAKE_1ARRAY(dat->ic2, dat->n_observations);
	MAKE_1ARRAY(dat->live, opt->K);
	MAKE_1ARRAY(dat->index, dat->n_observations);
	
	///* Sample seeds from reads, use initialization in kmodes */
	//unsigned int idx[opt->K];
	//if (RAND_SEED) {
	//	srand((unsigned int) time(NULL));
	//	sample_better(dat->n_observations, opt->K, idx);
	//} else {
	//	idx[0] = 111; idx[1] = 2222; idx[2] = 1000; // input seed
	//}
	//MAKE_2ARRAY(dat->seeds, opt->K, dat->max_read_length);
	//for (int i = 0; i < opt->K; ++i)
	//	COPY_1ARRAY(dat->seeds[i], dat->dmat[idx[i]], dat->max_read_length);
	//

	/* Allocate for writing results */
	dat->cluster_id = malloc(dat->n_observations * sizeof *dat->cluster_id);
	dat->best_cluster_id = malloc(dat->n_observations
						* sizeof *dat->best_cluster_id);
	
	if (dat->cluster_id == NULL || dat->best_cluster_id == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_id");
	
	//MAKE_2ARRAY(dat->best_modes, opt->K, dat->max_read_length);
	
	dat->cluster_size = malloc(opt->K * sizeof *dat->cluster_size);
	dat->best_cluster_size = malloc(opt->K * sizeof *dat->best_cluster_size);
	
	if (dat->cluster_size == NULL || dat->best_cluster_size == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::cluster_size");
	
	dat->criterion = malloc(opt->K * sizeof *dat->criterion);
	dat->best_criterion = malloc(opt->K * sizeof *dat->best_criterion);
	
	if (dat->criterion == NULL || dat->best_criterion == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data::criterion");
	
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
	
	return err;
}/* sync_state_data */

/**
 * Set up seeds
 
 * @param dat data
 * @param opt option
 * @return
 */
int make_seeds(data *dat, options *opt)
{
	
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
	
	data_t *tmp1 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp2 = malloc(opt->K * dat->n_coordinates * sizeof **dat->seeds);
	data_t *tmp3 = NULL;
	if (!tmp1 || !tmp2)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::seeds");

	// [KSD, BUG] data::ini_seeds is not allocated, so if
	// options::n_inner_init > 1, there will be a segmentation fault.
	if (opt->n_inner_init > 1) {
		tmp3 = malloc(opt->K * dat->n_coordinates
							* sizeof **dat->seeds);
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
	
	return NO_ERROR;
} /* make_seeds */

void free_data(data *dat)
{
	if (!dat) return;
	
	if (dat->fdata) free_fastq(dat->fdata);
	if (dat->dmat) free(dat->dmat);
	if (dat->qmat) free(dat->qmat);
	if (dat->read_idx) free(dat->read_idx);
	if (dat->lengths) free(dat->lengths);
	if (dat->true_cluster_id) free(dat->true_cluster_id);
	if (dat->true_cluster_size) free(dat->true_cluster_size);
	FREE_2ARRAY(dat->prob_f);
	FREE_2ARRAY(dat->prob_t);
	FREE_2ARRAY(dat->v_ik);
	FREE_3ARRAY(dat->e_kjN);
	FREE_1ARRAY(dat->last_assign);
	FREE_2ARRAY(dat->last_cent);
	FREE_1ARRAY(dat->last_cost);
	FREE_2ARRAY(dat->last_vik);
	FREE_3ARRAY(dat->last_ekjn);
	FREE_1ARRAY(dat->live);
	FREE_1ARRAY(dat->ic2);
	FREE_VECTOR(dat->index);
	//if (dat->coverage) free (dat->coverage);
	
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
	if (dat->n_categories) free(dat->n_categories);
	if (dat->categories) free(dat->categories);
	free(dat);
	
} /* free_data */

/**
 * Write solution
 *
 * @param dat		data object pointer
 * @param opt		options object pointer
 * @param in_fps	solution file pointer
 */
int write_solution(data *dat, options *opt, FILE **in_fps)
{
	/* report best solution to stdout and outfile */
	
	if (opt->continue_run)	/* continuing a previous run */
		*in_fps = fopen(opt->soln_file, "a");
	else			/* starting a fresh run */
		*in_fps = fopen(opt->soln_file, "w");
	if (*in_fps == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->soln_file);
	FILE *fps = *in_fps;
	
	//for (unsigned int k = 0; k < opt->K; ++k)
	//	memcpy(dat->best_modes[k], dat->seeds[k],
	//		dat->max_read_length * sizeof **dat->best_modes);
	//COPY_1ARRAY(dat->best_cluster_id, dat->cluster_id, dat->n_observations);
	//COPY_1ARRAY(dat->best_cluster_size, dat->cluster_size, opt->K);
	//COPY_1ARRAY(dat->best_criterion, dat->criterion, opt->K);
	
	if (fps) {
		fprintf(fps, "Best optimized criterion: %.0f\n",
			dat->best_total);
		fprint_doubles(fps, dat->best_criterion, opt->K, 0, 1);
	}
	
	if (fps && opt->K > 0) {
		fprintf(fps, "Best cluster sizes:");
		fprint_uints(fps, dat->best_cluster_size, opt->K, 0, 1);
	}
	
	if (fps && opt->K > 0) {
		fprintf(fps, "Best solution cluster assignments:\n");
		fprint_uints(fps, dat->best_cluster_id,
						dat->n_observations, 0, 1);
		fprintf(fps, "\n");
	}
	
	if (fps) {
		fprintf(fps, "Best modes:\n");
		for (unsigned int k = 0; k < opt->K; ++k)
			fprint_data_ts(fps, dat->best_modes[k],
						dat->n_coordinates, 0, 1);
	}
	
	if (fps) {
		fprintf(fps, "Time cost:\n");
		
	}
	fclose(fps);
	return NO_ERROR;
} /* write_solution */

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
	data_t *dptr, *qptr;
	
	for (i = n - 1; i > 0; --i) {
		lim = RAND_MAX - RAND_MAX % (i + 1);
		do {
			j = rand();
		} while (j >= lim);
		j = j % (i + 1);
		
		dptr = dat->dmat[j];
		dat->dmat[j] = dat->dmat[i];
		dat->dmat[i] = dptr;
		
		qptr = dat->qmat[j];
		dat->qmat[j] = dat->qmat[i];
		dat->qmat[i] = qptr;
		
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

/* Read true seeds */

void read_fsa(const char *filename, int n, int p, data_t **true_seed) {
	
	FILE *fp = fopen(filename, "r");
	
	if (fp == NULL)
		fprintf(stderr, "Unable to open file %s.\n", filename);
	
	char temp;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < p; ++j) {
			fscanf(fp, "%c", &temp);
			printf("%c", temp);
			
			switch (temp) {
				case 'A':
					true_seed[i][j] = 0;
					break;
					
				case 'C':
					true_seed[i][j] = 1;
					break;
					
				case 'G':
					true_seed[i][j] = 3;
					break;
					
				case 'T':
					true_seed[i][j] = 2;
					break;
			}
			printf("%d", true_seed[i][j]);
		}
	fclose(fp);
}
