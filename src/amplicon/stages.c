/**
 * @file stages.c
 * @author Xiyu Peng
 *
 * Functions for multi-stage ampliclust.
 *
 * [TODO, BUG, KSD] I believe there is a bug in the use of data::n_quality and
 * model::n_quality (it used to be called model::n_quality also).  It will only
 * manifest when quality scores are compressed.  They are used  interchangeably
 * throughout this code, but only one is correct.
 */

#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#include "error.h"
#include "run_ampliclust.h"
#include "ampliclust.h"
#include "model.h"
#include "options.h"
#include "data.h"
#include "simulate.h"
#include "stages.h"
#include "fastq.h"
#include "math.h"
#include "io.h"
#include "aecm.h"
#include "hash.h"
#include "sample.h"


/**
 * If the data size is too large, for example, larger than the choosen size of
 * sample, we will perform a multi-stage algorithm to classify those reads.
 *
 * @param mod	model object pointer
 * @param dat	data object pointer
 * @param opt	options object pointer
 * @param ini	initializer object pointer
 * @param sim	simulator object pointer
 * @param ri	run_info object pointer
 * @param stag	stages object pointer
 * @return	error status
 */
int multi_stage(model *mod, data *dat, options *opt, initializer *ini,
	simulator *sim, run_info *ri, stages *stag)
{
	UNUSED(sim);

	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;  /* DEBUG_I*/
	size_t next_D_size = 0;

	next_D_size = stag->D_size[0];

	if (!next_D_size)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"Size of the initial dataset is 0");

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "%i \n", next_D_size);

	FILE *fps = fopen("sample.txt", "w");
	if (!fps)
		return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR,"sample.txt");

	while (next_D_size) {
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug,"next_D_size:%i \n",
								next_D_size);

		unsigned int current = stag->current_stage; 	/* [KSD] Why do you need a local variable distinct from stage::current_stage? */
		/* it is used in line 203. keep it unless I find a better strategy */
		mmessage(INFO_MSG, NO_ERROR, "Current Stage: %i \n", stag->current_stage);

		mmessage(INFO_MSG, NO_ERROR, "Number of the Reads in  Stage "
			"%i: %i \n", stag->current_stage, stag->D_size[stag->current_stage]);

		
		if (stag->D_size[stag->current_stage] > 100) {	/* [TODO] do not hard-code 100 */

			/* opt->sample_size is larger than the remaining available data (happens in the final stage) */
			if (opt->sample_size >= stag->D_size[stag->current_stage]) {
				stag->sample_size[stag->current_stage]= stag->D_size[stag->current_stage];
				memcpy(stag->sample_idx, stag->D_idx,stag->D_size[stag->current_stage]
				* sizeof *stag->sample_idx);
			}else{
				stag->sample_size[stag->current_stage]=opt->sample_size;
				/* find sample_idx at this stage */
				if ((err = random_sample(stag->D_size[stag->current_stage],
					stag->sample_size[stag->current_stage],stag->D_idx,stag->sample_idx)))
					return err;
			}
			fprintf(fps, "S[%i]: ",stag->current_stage);
			fprint_size_ts(fps, stag->sample_idx, stag->sample_size[stag->current_stage], 2, 1);
		
		/* if size of D_s is too small, merge remaining reads with 
		 * previous sample */
		} else {

			/* If size of D_s is too small */
			if (stag->D_size[stag->current_stage]< 10 )
				break;

			/* reverse effect of last iteration */
			--stag->current_stage;
			stag->sum_K = stag->sum_K - stag->all_K[stag->current_stage];

			/* create a merged sample (change sample_size and sample_idx */
			if ((err = reconstruct_sample(stag, opt)))
				return err;
		}

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "sample_size: %i\n",
			stag->sample_size[stag->current_stage]);

		/* update data object preparing ampliclustK */
		if ((err = update_data(dat, opt,stag->sample_size[stag->current_stage],
			stag->sample_idx)))
			return err;
	
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "size of D_i: %i\n",
			stag->D_size[stag->current_stage]);

		/* Do ampliclustK() on the sample set */
		if (opt->estimate_K) {
			if ((err = ampliclustK(dat, mod, opt, ini, ri)))
				return err;
		} else { /* codes here just for test */
			if ((err = realloc_model(mod, dat, opt->modo)))
				return err;
			if ((err = realloc_run_info(ri, dat->sample_size,
				mod->n_mix)))
				return err;
			if ((err = realloc_initializer(ini, dat, opt->inio)))
				return err;
			if ((err = cluster_amplicons(dat, mod, ini, ri->opt->modo,
	                        ri->opt->inio, do_per_iterate, do_initialize,
       				(void *) ri)))
				return err;
			mmessage(INFO_MSG, NO_ERROR, "aic:%f; bic:%f.\n",
			mod->aic, mod->bic);
		}

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "size of D_i: "
			"%i,%i,%i\n", stag->D_size[0], stag->D_size[1],
			stag->D_size[2]);

		/* Code below would not be used, since we do hypothesis test on the D_s */
		/* move index of reads in S_s to end of index of reads in D_s */
		/*
		if ((err = reorder_D_idx(stag,SAMPLE)))
			return err;
		*/

		/* If D_s\S_s = NULL, it is the end of the algorithm */
		size_t remain_size = stag->D_size[stag->current_stage]
			- stag->sample_size[stag->current_stage];

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "current stage:%i; "
			"check remaining size: %i - %i \n", stag->current_stage,
			stag->D_size[stag->current_stage],
			stag->sample_size[stag->current_stage]);
		
		if (remain_size > 0) {
			
			/* generate p_values of reads in D_s */
			if ((err = pvalue_generator(dat, stag, mod, opt)))
				return err;

			/* Note: alpha values decrease with order */
			/* 0.05,0.025,0.010,0.005,0.001 */
			double alpha[ALPHA_SIZE] = {(double) ALPHA_1/ALPHA_NU,
				(double) ALPHA_2/ALPHA_NU, (double) ALPHA_3/ALPHA_NU,
				(double)ALPHA_4/ALPHA_NU, (double)ALPHA_5/ALPHA_NU};
			unsigned int count[ALPHA_SIZE] = {0, 0, 0, 0, 0};

			/* count the size of each alpha */
			for (size_t i = 0; i < stag->D_size[stag->current_stage]; ++i) {
				for (int a = 0; a < ALPHA_SIZE; ++a) {
					if (stag->pvalue[i] < alpha[a])
						count[a]++;
				}
			}

			debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "alpha "
				"%f,%f,%f,%f,%f\n", alpha[0], alpha[1],
				alpha[2], alpha[3], alpha[4]);

			next_D_size = 0;

			for (int a = 0; a < ALPHA_SIZE; ++a) {
				if (count[a] > stag->D_size[stag->current_stage] * alpha[a]) {
					stag->confi_level[stag->current_stage] = (double) alpha[a];
					next_D_size = count[a];

					/* [TODO] maximum stages=100 */

					stag->D_size[stag->current_stage + 1]
						= next_D_size;
					if ((err = reorder_D_idx(stag, D_SET)))
						return err;
					++stag->current_stage;
					if (stag->current_stage >= 99)	/* [TODO] do not hard-code */
						return mmessage(ERROR_MSG, INTERNAL_ERROR,
							"exceed maximum number of stages");
					mmessage(INFO_MSG, NO_ERROR, "next "
						"stage %i has %i reads\n",
						stag->current_stage,
						next_D_size);

					break;
				}
			}
		} else {
			next_D_size = 0; /* the end */
		}

		/* store the result from model object in stages object after ampliclust K */
		if ((err = update_stages(stag, dat, mod, opt, current))) 
			return err;

	}

	/* rescale stag->all_pi (optional, if we have dropped several reads in the final) */
	double sum = 0.;
	for (unsigned int k = 0; k < stag->sum_K; ++k)
		sum += stag->all_pi[k];
	for (unsigned int k = 0; k < stag->sum_K; ++k)
		stag->all_pi[k] /= sum;

	/* reestimate cluster here */
	if ((err = reestimate_parameters(dat, stag, opt)))
		return err;

	/* final assignment of reads */
	if ((err = final_assignment(dat, stag, opt, ri)))
		return err;

	/* print the current result */
	mmessage(INFO_MSG, NO_ERROR, "print the result.....\n");

	unsigned int max_read_position = dat->fdata->n_max_length
		+ dat->max_offset_ori;
	
	FILE *fp = fopen(opt->outfile, "w");
	if (!fp)
		return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR,
			opt->outfile);
	fprintf(fp, "All K: %i\n ", stag->sum_K);
	fprint_fasta(fp, stag->all_haplotypes, stag->sum_K,
			dat->fdata->n_max_length, "H");
	fprintf(fp, "All pi: ");
	fprint_doubles(fp, stag->all_pi, stag->sum_K, 6, 1);

	unsigned int pre_K = 0;
	for (unsigned int i = 0; i <= stag->current_stage; ++i) {
		fprintf(fp, "stage: %u\n", i);

		fprintf(fp,"Number of reads in stage %u ( D_s ): %zu\n", i,
						stag->D_size[i]);

		fprintf(fp, "best_K: %u\n", stag->all_K[i]);

		fprintf(fp, "pi: ");
		fprint_doubles(fp, stag->all_pi + pre_K, stag->all_K[i], 6, 1);

		fprint_fasta(fp, stag->all_haplotypes + pre_K
			* dat->fdata->n_max_length, stag->all_K[i],
			dat->fdata->n_max_length, "H");

		fprintf(fp, "lambda0:\n");
		fprint_vectorized_matrix(fp, stag->lambda0_s
			+ i * max_read_position * dat->n_quality,
			max_read_position, dat->n_quality, 1);
		fprintf(fp, "lambda1:\n");
		fprint_vectorized_matrix(fp, stag->lambda1_s
			+ i * max_read_position * dat->n_quality,
			max_read_position, dat->n_quality, 1);
		fprintf(fp, "delta: \n");
		fprint_doubles(fp, stag->delta_s + i * max_read_position,
					  max_read_position, 6, 1);

		fprintf(fp, "gamma: \n");
		fprint_vectorized_sq_matrix(fp, stag->gamma_s
			+ i * NUM_NUCLEOTIDES * NUM_NUCLEOTIDES,
					NUM_NUCLEOTIDES, 1);
		pre_K += stag->all_K[i];
	}

	/* print the best result */
	mmessage(INFO_MSG, NO_ERROR, "now print the best result.....\n");

	fprintf(fp, "best K: %i\n ", stag->best_K);
	fprint_fasta(fp, stag->best_haplotypes, stag->best_K,
		dat->fdata->n_max_length, "H");
	fprintf(fp, "best pi: ");
	fprint_doubles(fp, stag->best_pi, stag->best_K, 6, 1);

	fprintf(fp, "assignments: ");
	fprint_uints(fp, ri->optimal_cluster_id, dat->fdata->n_reads, 2, 1);

	fprintf(fp, "log likelihood: %f\n", stag->ll);

	fprintf(fp, "D[0]: ");
	fprint_size_ts(fp, stag->D_idx, dat->fdata->n_reads, 2, 1);	

	mmessage(INFO_MSG, NO_ERROR, "finished\n");

	fclose(fp);
	fclose(fps);	

	return err;
} /* multi_stage */

/**
 * Initialize the stages object
 *
 * @param stag	pointer to stages object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @return	error status
 */
int make_stages(stages **stag, data *dat, options *opt)
{

	int fxn_debug = ABSOLUTE_SILENCE;
	int maximum = 100; /* maximum number of stages */

	stages *rsta;
	*stag = malloc(sizeof **stag);

	if (*stag == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages");

	rsta = *stag;

	rsta->D_idx = calloc(dat->fdata->n_reads,sizeof * rsta->D_idx);

	if (!rsta->D_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.Didx");

	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		rsta->D_idx[i] = dat->read_idx[i];

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "D_idx");

	/* suppose stages are less than maximum=100 */
	rsta->D_size = calloc(maximum, sizeof *rsta->D_size);
	rsta->sample_size = calloc(maximum, sizeof *rsta->sample_size);

	if(!rsta->D_size || !rsta->sample_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.size");

	rsta->current_stage = 0; /* initial value */

	rsta->D_size[0] = dat->fdata->n_reads;

	rsta->sample_idx = calloc(opt->sample_size, sizeof *rsta->sample_idx);
	if (!rsta->sample_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.sample");

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "sample_idx");


	rsta->max_offset = NULL;
	rsta->max_read_length = NULL;
	rsta->max_read_position = NULL;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "offset");


	/* rsta->pvalue=calloc((rsta->D_size[0]-rsta->sample_size[0]),sizeof *rsta->pvalue);*/
	rsta->pvalue = calloc(rsta->D_size[0], sizeof *rsta->pvalue);

	if(!rsta->pvalue)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.pvalue");

	rsta->confi_level = calloc(maximum, sizeof *rsta->confi_level);

	if (!rsta->confi_level)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.confi_level");


	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "pvalue");

	rsta->simu_read = malloc(dat->fdata->n_max_length
		* sizeof *rsta->simu_read);
	rsta->simu_qual = malloc(dat->fdata->n_max_length
		* sizeof *rsta->simu_qual);

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "simulate_read_H0");

	/* we generate 10000 reads to get the distribution of the test_stat*/
	rsta->simu_test_stat = calloc(opt->simu_size
		* (dat->max_offset_ori + 1), sizeof * rsta->simu_test_stat);

	if (!rsta->simu_test_stat || !rsta->simu_read || !rsta->simu_qual)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.simu");

	rsta->all_K = calloc(maximum,sizeof * rsta->all_K);
	rsta->sum_K = 0;
	rsta->best_K = 0;

	if(!rsta->all_K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.K");

	rsta->all_haplotypes = NULL;
	rsta->all_pi = NULL;
	rsta->best_haplotypes = NULL;
	rsta->best_pi = NULL;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "K,haplotypes,pi");

	/* to save the memory space, just malloc room of delta_s for one stage */

	unsigned int max_read_position = dat->max_offset_ori
		+ dat->fdata->n_max_length;

	rsta->delta_s = malloc(max_read_position * sizeof *rsta->delta_s);
	rsta->best_delta = malloc(max_read_position * sizeof *rsta->best_delta);

	if(!rsta->delta_s || !rsta->best_delta)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.delta");

	rsta->gamma_s = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
		* sizeof *rsta->gamma_s);
	rsta->best_gamma = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
		* sizeof *rsta->best_gamma);

	if (!rsta->gamma_s || !rsta->best_gamma)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.gamma");

	rsta->lambda0_s = malloc(max_read_position * dat->n_quality
		* sizeof *rsta->lambda0_s);
	rsta->best_lambda0 = malloc(max_read_position * dat->n_quality
		* sizeof *rsta->best_lambda0);

	if(!rsta->lambda0_s || !rsta->best_lambda0)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.lambda0");

	rsta->lambda1_s = malloc(max_read_position * dat->n_quality
		* sizeof *rsta->lambda1_s);
	rsta->best_lambda1 = malloc(max_read_position * dat->n_quality
		* sizeof *rsta->best_lambda1);

	if (!rsta->lambda1_s || !rsta->best_lambda1)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "stages.lambda1");

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Dispersion parameter");

	rsta->eik = NULL;
	rsta->ll = -INFINITY;

	return NO_ERROR;
}/* make_stages */


/**
 * Store the result from model object in stages object.
 * need to increase size of some components in stages object (realloc)
 *
 * [KSD] to verify
 *
 * @param stag		pointer to stages object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param opt		pointer to options object
 * @param current	current iteration
 * @return err		error status
 */
int update_stages(stages *stag, data *dat, model *mod, options *opt,
	unsigned int current)
{
	int err = NO_ERROR;
	unsigned int n_stages = current + 1;

	/* K */
	stag->all_K[current] = mod->K;
	stag->sum_K = stag->sum_K + mod->K;

	/* haplotypes */
	if (!dat->fdata->n_max_length || !stag->sum_K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.haplotypes");
	unsigned char *all_haplotypes = realloc(stag->all_haplotypes,
		dat->fdata->n_max_length * stag->sum_K
		* sizeof *stag->all_haplotypes);
	if (!all_haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.haplotypes");

	stag->all_haplotypes = all_haplotypes;
	unsigned char *start_h = stag->all_haplotypes + dat->fdata->n_max_length
		* (stag->sum_K - mod->K);
	memcpy(start_h, mod->best_haplotypes, dat->fdata->n_max_length
		* mod->K * sizeof *stag->all_haplotypes);

	/* pi */
	if (!mod->n_mix)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.pi");
	double *all_pi = realloc(stag->all_pi, sizeof *stag->all_pi *
		(stag->sum_K + n_stages * opt->modo->background_model));
	if (!all_pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.pi");
	stag->all_pi = all_pi;
	double *start_pi = stag->all_pi + stag->sum_K
		+ n_stages * opt->modo->background_model - mod->n_mix;

	/* Need to scale pi */
	double scale = (double) (stag->D_size[current]
		- stag->D_size[current + 1]) / stag->D_size[0];

	for (unsigned int i = 0; i < mod->n_mix; ++i)
		mod->best_pi[i] = mod->best_pi[i] * scale;

	memcpy(start_pi, mod->best_pi, mod->n_mix * sizeof *stag->all_pi);

	/* delta */

	unsigned int max_read_position = dat->fdata->n_max_length
		+ dat->max_offset_ori;

	double *delta_s = realloc(stag->delta_s,
		max_read_position * n_stages * sizeof *stag->delta_s);
	if (!delta_s)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.delta");
	stag->delta_s = delta_s;
	double *start_d = stag->delta_s + max_read_position * (n_stages - 1);
	memcpy(start_d, mod->best_delta,
		max_read_position * sizeof *stag->delta_s);

	/* gamma */
	double *gamma_s = realloc(stag->gamma_s, n_stages * NUM_NUCLEOTIDES
		* NUM_NUCLEOTIDES * sizeof *stag->gamma_s);
	if (!gamma_s)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.gamma");
	stag->gamma_s = gamma_s;

	double *start_g = stag->gamma_s + (n_stages - 1)
		* NUM_NUCLEOTIDES * NUM_NUCLEOTIDES;
	memcpy(start_g, mod->best_gamma,
		NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof *stag->gamma_s);

	/* lambda */
	if (!mod->n_quality)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.lambda");
	double *lambda0_s = realloc(stag->lambda0_s, n_stages
		* max_read_position * mod->n_quality * sizeof *stag->lambda0_s);
	double *lambda1_s = realloc(stag->lambda1_s, n_stages
		* max_read_position * mod->n_quality * sizeof *stag->lambda1_s);
	if (!lambda0_s || !lambda1_s)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.realloc.lambda");

	stag->lambda0_s = lambda0_s;
	stag->lambda1_s = lambda1_s;

	double *start_l0 = stag->lambda0_s + (n_stages - 1)
			* max_read_position * mod->n_quality;
	double *start_l1 = stag->lambda1_s + (n_stages - 1)
			* max_read_position * mod->n_quality;

	memcpy(start_l0, mod->best_lambda0, max_read_position
			* mod->n_quality * sizeof *stag->lambda0_s);
	memcpy(start_l1, mod->best_lambda1, max_read_position
			* mod->n_quality* sizeof *stag->lambda1_s);

	return err;
} /* update_stages */

/**
 * Free stages object
 *
 * @param stag	pointer to stages object
 */
void free_stages(stages *stag)
{
	if(stag) {
		if (stag->D_idx) free(stag->D_idx);
		if (stag->D_size) free(stag->D_size);
		if (stag->sample_size) free(stag->sample_size);

		if (stag->sample_idx) free(stag->sample_idx);

		if (stag->max_offset) free(stag->max_offset);
		if (stag->max_read_length) free(stag->max_read_length);
		if (stag->max_read_position) free(stag->max_read_position);

		if (stag->pvalue) free(stag->pvalue);
		if (stag->confi_level) free(stag->confi_level);
		if (stag->simu_read) free(stag->simu_read);
		if (stag->simu_qual) free(stag->simu_qual);
		if (stag->simu_test_stat) free(stag->simu_test_stat);

		if (stag->all_K) free(stag->all_K);
		if (stag->all_haplotypes) free(stag->all_haplotypes);
		if (stag->all_pi) free(stag->all_pi);
		if (stag->delta_s)free(stag->delta_s);
		if (stag->gamma_s)free(stag->gamma_s);
		if (stag->lambda0_s)free(stag->lambda0_s);
		if (stag->lambda1_s)free(stag->lambda1_s);

		if (stag->best_pi) free(stag->best_pi);
		if (stag->best_haplotypes) free(stag->best_haplotypes);
		if (stag->best_delta)free(stag->best_delta);
		if (stag->best_gamma)free(stag->best_gamma);
		if (stag->best_lambda0) free(stag->best_lambda0);
		if (stag->best_lambda1) free(stag->best_lambda1);

		if (stag->eik) free(stag->eik);
		free(stag);
	}
}/* free_stages */

/**
 * Merge penultimate sample and remainder when final stage D_s is too small.
 *
 * Note, we will have to repeat the previous stage....
 *
 * @param stag	pointer to stages object
 * @param opt	pointer to options object
 * @return	error status
 */
int reconstruct_sample(stages *stag, options *opt)
{
	/* current is previous stage */
	unsigned int current = stag->current_stage;
	unsigned int ind;
	size_t num=0;

	/* size of the new sample: current + 1 is aborted stage */
	stag->sample_size[current] = opt->sample_size
		+ stag->D_size[current + 1];

	/* realloc the sample idx */
	size_t *sample_idx = realloc(stag->sample_idx,
		stag->sample_size[current] * sizeof *stag->sample_idx);
	if (!sample_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.sample_idx");
	stag->sample_idx = sample_idx;

	/* integrate reads in small D_s into previous sample set (should avoid duplicates) */
	for(size_t i=0; i <stag->D_size[current+1]; i++){
		ind = 0;
		for(size_t j=0; j <stag->sample_size[current];j++){
			if (stag->D_idx[i]==stag->sample_idx[j]){
				ind = 1;
				break;
			}
		}
		if(ind==0){
			size_t idx = opt->sample_size+num;
			stag->sample_idx[idx]= stag->D_idx[i];
			num++;
		}
	}
	stag->sample_size[current] = opt->sample_size +num;

	/*
	memcpy(stag->sample_idx + opt->sample_size, stag->D_idx,
		stag->D_size[current + 1] * sizeof *stag->sample_idx); */
			/* [BUG] This looks like a bug to me, where you copy past end of array: should be "sizeof *sample_idx". */

	return NO_ERROR;
}/* reconstruct_sample */

/**
 * reorder D_idx for storing D_{s+1} at the beginning of D_s (type D_SET)
 * or storing S_s at the end of D_s (type SAMPLE) .
 *
 * [KSD] not yet verified
 *
 * @param stag	pointer stages object
 * @param type	SAMPLE or D_SET
 * @return	error status
 */
int reorder_D_idx(stages *stag, int type)
{
	int fxn_debug = ABSOLUTE_SILENCE;  /* DEBUG_I */ /* DEBUG_II */
	unsigned int current = stag->current_stage;
	int count_s = 0;
	size_t next_D_size = 0;

	/* Move the index in the S_s at the end of D_s array*/
	if (type == SAMPLE) {/* Not be used any more   */
		for (size_t i = 0; i < stag->sample_size[current]; ++i) {
			for (size_t j = 0; j < stag->D_size[current]; ++j){
				if (stag->sample_idx[i] == stag->D_idx[j]) {
					debug_msg(DEBUG_I <= fxn_debug,
						fxn_debug, "before:sample_idx: "
						"%zu; D_idx(change): %zu\n",
						stag->D_idx[j], stag->D_idx[
						stag->D_size[current] - i - 1]);
					size_t trf = stag->D_idx[
						stag->D_size[current] - i - 1];
					stag->D_idx[
						stag->D_size[current] - i - 1]
						= stag->D_idx[j];
					stag->D_idx[j] = trf;
					debug_msg(DEBUG_I <= fxn_debug,
						fxn_debug, "after:sample_idx: "
						"%zu; D_idx(change): %zu\n",
						stag->D_idx[j], stag->D_idx[
						stag->D_size[current] - i - 1]);
					++count_s;
					break;
				}

			}
		}
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "sample_size:%i, "
			"count: %i\n", stag->sample_size[current], count_s);

	/* Based on the final confidence level chosen and p-values for all reads,
	 * store the index of D_{s+1} at the beginning of the D_idx.
	 */
	} else if(type == D_SET) {
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug,
						"current stages: %i", current);
		size_t size_p = stag->D_size[current];
		double alpha = stag->confi_level[current];
		for (size_t i = 0; i < size_p; ++i) {
			if (stag->pvalue[i] < alpha) {
				debug_msg(DEBUG_II <= fxn_debug, fxn_debug,
					"before:d_i+1dx: %zu (location %zu); "
					"D_idx(change): %zu (location %zu)\n",
					stag->D_idx[i], i,
					stag->D_idx[next_D_size],
					next_D_size);
				size_t trf = stag->D_idx[next_D_size];
				stag->D_idx[next_D_size] = stag->D_idx[i];
				stag->D_idx[i] = trf;
				debug_msg(DEBUG_II <= fxn_debug, fxn_debug,
					"after:d_i+1dx: %zu (location %i); "
					"D_idx(change): %zu (location %zu)\n",
					stag->D_idx[i], i,
					stag->D_idx[next_D_size], next_D_size);
				++next_D_size;
			}
		}
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "alpha:%f, "
			"next_D_size: %zu\n", alpha, next_D_size);
	}

	return NO_ERROR;
} /* reorder_D_idx */

/**
 * Generate an array of p-values for reads in D_s/S_s
 * Note that the size of p-values array is different at every stage,
 * and it is determined by current D_size and current sample_size.
 *
 * [TODO] There is probably some opportunity for speed-up in this code.
 *
 * @param dat	pointer to data object
 * @param stag	pointer to stages object
 * @param mod	pointer to model object
 * @param opt	pointer to options object
 * @return	    error status
 */
int pvalue_generator(data *dat, stages *stag, model *mod, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	unsigned int count_offset=dat->max_offset_ori+1; /* max_offset starts with 0 */

	/* Simulate the dsitribution of the test statistics */
	for (unsigned int i = 0; i < count_offset; ++i) {
		for (size_t j = 0; j < opt->simu_size; ++j) {
			/* [TODO] handle variable lengths and trimmed lengths */
			int len = dat->max_read_length;

			if ((err = simulate_read(dat, mod, len,
					stag->simu_read, stag->simu_qual, i)))
				return err;

			stag->simu_test_stat[i*opt->simu_size + j] =
				likelihood_ratio(mod, stag->simu_read,
						stag->simu_qual, len, i,dat->max_read_length);

			/* log(ratio) should smaller than 0 */
			if (stag->simu_test_stat[i*opt->simu_size + j] >= 0)
				return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"simulate.likelihood_ratio");

		}
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "simulation size:%i\n",
								opt->simu_size);

	/* Generate the p-value for each read in the D_s */
	size_t pvalue_size = stag->D_size[stag->current_stage];

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "pvalue size:%i\n",
								pvalue_size);
	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "D size:%i\n",
					stag->D_size[stag->current_stage]);

	for (size_t i = 0; i < pvalue_size; ++i) {
		size_t r_idx = stag->D_idx[i];
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "r_idx:%i\n", r_idx);

		unsigned char *rptr = dat->fdata->reads;
		unsigned char *qptr = dat->fdata->quals;

		/* [TODO] This is slow. */
		for (size_t i = 0; i < r_idx; ++i) {
			rptr += read_length(dat->fdata, i);
			qptr += read_length(dat->fdata, i);
		}

		unsigned int off = dat->offset_ori ? dat->offset_ori[r_idx] : 0;
		double ratio = likelihood_ratio(mod, rptr, qptr,
			read_length(dat->fdata, r_idx), off,dat->max_read_length);

		if (ratio >= 0)  /* log(ratio) should smaller than 0 */
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
				"reads.likelihood_ratio");

		unsigned int count = 0;
		for (size_t j = 0; j < opt->simu_size; ++j)
			if (stag->simu_test_stat[off*opt->simu_size + j] <= ratio)
				++count;

		stag->pvalue[i] = (double)count / opt->simu_size;

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "pvalue:%f\n",
							stag->pvalue[i]);
	}
	return NO_ERROR;
}/* pvalue_generator */


/**
 * This function calculate the likelihood ratio as the test statistic.
 *
 * [KSD] I need to check this code since I did not do it on 1/15/18.
 *
 * @param mod	           model object
 * @param read	           tested nucleotide sequence
 * @param qual	           tested quality score sequence
 * @param len	           length of the sequence
 * @param off	           the trimmed length of the sequence (usually 0)
 * @param max_read_length  maximum read length in data object
 * @return	log likelihood ratio
 */
double likelihood_ratio(model *mod, unsigned char *read,
	unsigned char *qual, unsigned int len, unsigned int off,unsigned int max_read_length)
{
	double log_llr = 0., sum = 0., p = 0.;
	double max = -INFINITY;
	unsigned char *haps = mod->best_haplotypes;

	double *sum_K = malloc(mod->K * sizeof *sum_K);	/* [KSD,BUG] was sizeof(double *) */

	if (!sum_K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.likelihood_ratio.sum_K");

	/* initialization */
	for (unsigned int k = 0; k < mod->K; ++k)
		sum_K[k] = 0.;

	/* [TODO] consider the background model */

	/* calculate log value of the numerator of the test statistic */

	for (unsigned int k = 0; k < mod->K; ++k) {
		sum_K[k] = (double)log(mod->best_pi[k]); /* probability of cluster k */
		for (unsigned int j = 0; j < len; ++j) {
			if (read[j] != haps[k * max_read_length + j]) {
				sum_K[k] += log(1 - mod->best_delta[j + off])	 /* prob. error */
					+ log(mod->best_lambda1[(j + off)
					* mod->n_quality + qual[j]])		/* prob. of quality score */
					+ log(mod->best_gamma[haps[
					k * max_read_length + j]
					* NUM_NUCLEOTIDES + read[j]]);		/* prob. h_{kj} -> x_{ij} given error */
			} else {
				sum_K[k] += log(mod->best_delta[j + off])	/* prob. no error */
					+ log(mod->best_lambda0[(j + off)
					* mod->n_quality + qual[j]]);		/* prob. of quality score */
			}
		}
		if (max < sum_K[k])
			max = sum_K[k];
	}
	for (unsigned int k = 0; k < mod->K; ++k) {
		sum_K[k] = exp(sum_K[k] - max);
		sum += sum_K[k];
	}
	log_llr = log(sum) + max;

	/* calculate log value of the denominator of the test statistic */

	/* [KSD,TODO] avoid allocating memory inside compute intensive code */
	unsigned char *haps_hat = malloc(max_read_length * sizeof * read);

	if(!haps_hat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"stages.likelihood_ratio.haps_hat");

	for (unsigned int j = 0; j < len; ++j) {
		/* check we need to reestimate hat{H_{s,K_s+1,j}} */
		double ratio_lambda = 
			mod->best_lambda0[(j + off) * mod->n_quality
				+ qual[j]] / mod->best_lambda1[(j + off)
					* mod->n_quality + qual[j]];
		double ratio_delta = (1 - mod->best_delta[j + off])
			/ mod->best_delta[j + off];

		/* use r_j */
		if (ratio_lambda > ratio_delta) {
			haps_hat[j] = read[j];

			log_llr -= log(mod->best_delta[j + off])
				+ log(mod->best_lambda0[(j + off)
				* mod->n_quality + qual[j]]);


			/* I guess the case below is very rare */
		} else {

			double current_max = -INFINITY;

			/* estimated a best nucleotide at each position */
			for (int i = 0; i < NUM_NUCLEOTIDES; ++i) {
				if (i == read[j])
					p = mod->best_delta[j + off]
						* mod->best_lambda0[(j + off)
						* mod->n_quality + qual[j]];
				else
					p = (1 - mod->best_delta[j + off])
						* mod->best_lambda1[(j + off)
						* mod->n_quality + qual[j]]
						* mod->best_gamma[
						i * NUM_NUCLEOTIDES + read[j]];
				if (p > current_max) {
					haps_hat[j] = i;
					current_max = p;
				}
			}

			/* Since all gamma are smaller than 1, read[j] == haps_hat[j] still could happen here */
			if (read[j] != haps_hat[j])
				log_llr -= log(1 - mod->best_delta[j + off])
					+ log(mod->best_lambda1[(j + off)
					* mod->n_quality + qual[j]])
					+ log(mod->best_gamma[haps_hat[j]
					* NUM_NUCLEOTIDES + read[j]]);
			else
				log_llr -= log(mod->best_delta[j + off]) /* prob. no error */
					+ log(mod->best_lambda0[(j + off)
					* mod->n_quality + qual[j]]);
		}

	}
	free(sum_K);	/* [TODO] repeated allocation/free is inefficient */
	free(haps_hat);

	return log_llr;
} /* likelihood_ratio */

/**
 * Final estimates of K and pi by reclassification
 *
 * [KSD] to verify
 *
 * @param dat	pointer to data object
 * @param stag	pointer to stages object
 * @param opt	pointer to options object
 * @return	error status
 */
int reestimate_parameters(data *dat, stages *stag, options *opt)
{
	/*int fxn_debug = ABSOLUTE_SILENCE;*/
	int err = NO_ERROR;

	size_t sample_size = (size_t) (dat->fdata->n_reads * opt->proportion);

	/* realloc stag->sample_idx for the given sample_size */
	if (!sample_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"reestimate_parameters.realloc.sample_idx");
	size_t *sample_idx = realloc(stag->sample_idx,
			sample_size * sizeof *stag->sample_idx);
	if (!sample_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"reestimate_parameters.realloc.sample_idx");
	stag->sample_idx = sample_idx;

	/* update the sample_idx and data object */
	if(opt->proportion==1){
		for (size_t i = 0; i< sample_size ; ++i)
			stag->sample_idx[i] = i;
	} else {
		if ((err = random_sample(dat->fdata->n_reads, sample_size, stag->D_idx,stag->sample_idx)))
			return err;
	}

	if ((err = update_data(dat, opt, sample_size,stag->sample_idx)))
		return err;

	/* select or estimate dispersion parameters from the pool */
	/* Below I just choose to use dispersion parameters from the first stage */
	unsigned int max_read_position = dat->max_offset_ori
		+ dat->fdata->n_max_length;

	memcpy(stag->best_delta, stag->delta_s,
		max_read_position * sizeof * stag->best_delta);
	memcpy(stag->best_gamma, stag->gamma_s,
		NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof * stag->best_gamma);
	memcpy(stag->best_lambda0, stag->lambda0_s,
		max_read_position * dat->n_quality * sizeof *stag->best_lambda0);
	memcpy(stag->best_lambda1, stag->lambda1_s,
		max_read_position * dat->n_quality * sizeof *stag->best_lambda1);

	/* generate matrix E */
	if ((err = generate_E(stag, dat, stag->sum_K,NULL)))
		return err;

	/* reestimate K and haplotypes */
	if ((err = eliminate_nullcluster(stag, dat, opt)))
		return err;

	return err;
}/* reestimate_parameters */

/**
 * Reassign all reads in the complete data set to final clusters.
 *
 * [KSD] to verify
 *
 * @param dat	pointer to data object
 * @param stag	pointer to stages object
 * @param opt	pointer to options object
 * @param ri	pointer to run_info object
 * @return	    error status
 */
int final_assignment(data *dat, stages *stag, options *opt, run_info *ri)
{
	int err = NO_ERROR;
	FILE *fp = fopen("ll_assignment.txt", "w");
	if (!fp)
		return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR,"ll_assignment.txt");

	/* if we just use part of the complete data set to reestimate
	 * dispersion parameters and eliminate null clusters, here we should realloc space for
	 * a larger data size */

	if (opt->proportion < 1) {

		size_t sample_size = dat->fdata->n_reads; /* complete data set */

		/* realloc stag->sample_idx for the given sample_size */
		size_t *sample_idx = realloc(stag->sample_idx,
			sample_size * sizeof *stag->sample_idx);
		if (!sample_idx)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"final_assignment.realloc.sample_idx");
		stag->sample_idx = sample_idx;

		/* update the sample_idx and data object */
		for (size_t i = 0; i<sample_size ; ++i)
			stag->sample_idx[i] = i;

		if ((err = update_data(dat, opt, sample_size,stag->sample_idx)))
			return err;
	}

	/* regenerate E matrix for an adjusted K and the complete data set */
	if ((err = generate_E(stag, dat, stag->best_K,fp)))
		return err;

	if ((err = realloc_run_info(ri, dat->fdata->n_reads, stag->best_K)))
		return err;

	if ((err = assign_clusters(stag->eik, stag->best_K, dat->fdata->n_reads,
		ri->optimal_cluster_size, ri->optimal_cluster_id,0)))
		return err;
	fclose(fp);

	return err;
}/* final_assignment */

/**
 * reidentify the number of clusters K and hyplotypes, and reestimate pi
 *
 * @param stag		pointer to stages object
 * @param dat		pointer to data object
 * @param opt		pointer to options object
 * @return		    error status
 */
int eliminate_nullcluster(stages *stag, data *dat, options *opt)
{
	/* int fxn_debug = ABSOLUTE_SILENCE;*/
	int err = NO_ERROR;
	unsigned int k, a, skm1 = stag->sum_K - 1;
	size_t i, hd = 0;

	double *pro_k = malloc(stag->sum_K * sizeof *pro_k);				/* [BUG] used to say sizeof (double *) */
	unsigned int* ind_eliminate = calloc(stag->sum_K, sizeof *ind_eliminate);	/* [BUG] Used to say sizeof (int *) */
	stag->best_pi = malloc(stag->sum_K * sizeof *stag->best_pi);
	stag->best_haplotypes = malloc(dat->fdata->n_max_length * stag->sum_K
		* sizeof *stag->best_haplotypes);

	/* Find the null cluster in the existing sum_K clusters  */
	for (k = 0; k < stag->sum_K; ++k) {
		pro_k[k] = 0.;
		for (i = 0; i < dat->sample_size; ++i)
			pro_k[k] += stag->eik[i*stag->sum_K + k];
		if (pro_k[k] < opt->epsilon_s)
			ind_eliminate[k] = 1;
	}

	/* check if any two haplotypes are same and keep only the former */
	for (k = 0; k < skm1; ++k) {	/* [TODO] compiler may optimize, but avoid stag->sum_K - 1 calculation per iteration */
		for (a = k + 1; a < stag->sum_K ; ++a) {
			hd = hamming_char_dis((char *) &stag->all_haplotypes[
				k * dat->fdata->n_max_length],
				(char * ) &stag->all_haplotypes[
				a * dat->fdata->n_max_length],
				dat->fdata->n_max_length);
			if (hd == 0) {
				ind_eliminate[a] = 1;
				stag->all_pi[k] += stag->all_pi[a];
			}
		}
	}

	/* remove null clusters and duplicated hyplotypes */
	stag->best_K = 0;
	for (k = 0; k < stag ->sum_K; ++k) {
		if (ind_eliminate[k] == 0) {
			memcpy(&stag->best_haplotypes[
				stag->best_K * dat->fdata->n_max_length],
				&stag->all_haplotypes[k * dat->fdata->n_max_length],
				dat->fdata->n_max_length);
			stag->best_pi[stag->best_K] = stag->all_pi[k];
			stag->best_K += 1 ;
		}
	}

	free(pro_k);
	free(ind_eliminate);

	return err;
}/* eliminate_nullcluster */

/**
 * E step to estimate e_ik for every reads in the chosen sample
 *
 * @param stag		pointer to stages object
 * @param dat		pointer to data object
 * @param K         K clusters
 * @param fp        pointer to file
 * @return		error status
 */
int generate_E(stages *stag, data *dat,unsigned int K,FILE *fp)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;

	double sum, max, ll = 0.;
	unsigned int length;
	size_t i, j, k;

	unsigned char *haplotypes;
	double *pi;

	/* realloc space for different K and different sample_size */
	if (!stag->eik) {
		stag->eik = malloc(dat->sample_size * K * sizeof *stag->eik);
		if (!stag->eik)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"generate_E.stag::eik");
	} else {
		double *eik = realloc(stag->eik, dat->sample_size * K
			* sizeof *stag->eik);
		if (!eik)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"generate_E.stag::eik");
		stag->eik = eik;
	}

	if (K < stag->sum_K) {
		/* if best_K < sum_K */
		haplotypes = stag->best_haplotypes;
		pi = stag->best_pi;
	} else {
		haplotypes = stag->all_haplotypes;
		pi = stag->all_pi;
	}

	/* Below are modified from function e_step() in aecm.c */
	for (i = 0; i < dat->sample_size; ++i) {
		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Read %u\n");
		max = -INFINITY;	/* to hold max_k e_{ik} */
		for (k = 0; k < K; ++k) {
			stag->eik[i* K + k] = log(pi[k]);	/* probability of cluster k */
			length = dat->lengths[i];
			for (j = 0; j < length; ++j) {

				if (dat->dmat[i][j] != haplotypes[k*dat->max_read_length + j]) {

					stag->eik[i*K + k] += log(1 - stag->best_delta[j + (dat->offset ? dat->offset[i] : 0)])				/* prob. error */
						+ log(stag->best_lambda1[(j + (dat->offset ? dat->offset[i] : 0))*dat->n_quality + dat->qmat[i][j]])	/* prob. q_{ij} given error */
						+ log(stag->best_gamma[haplotypes[k*dat->max_read_length + j]*NUM_NUCLEOTIDES + dat->dmat[i][j]]);	/* prob. h_{kj} -> x_{ij} given error */
				} else {
					stag->eik[i* K + k] += log(stag->best_delta[j + (dat->offset ? dat->offset[i] : 0)])				/* prob. no error */
						+ log(stag->best_lambda0[(j + (dat->offset ? dat->offset[i] : 0))*dat->n_quality + dat->qmat[i][j]]);	/* prob q_{ij} given no error */
				}
			}
			if (max < stag->eik[i*K + k])
				max = stag->eik[i*K + k];
		}

		sum = 0.;
		for (k = 0; k < K; ++k) {
			stag->eik[i*K + k] = exp(stag->eik[i*K + k] - max);	/* e_{ik} = e^{-max}\pi_k\prod_{j=1}^J ... */
			sum += stag->eik[i*K + k];				/* e^{-max}\sum_{k=1}^(K+1) \pi_k \prod_{j=1}^J ... */
		}
		if(fp)
			fprintf(fp, " %.3f", log(sum) + max);
	
		ll += log(sum) + max;

		for (k = 0; k < K; ++k)
			stag->eik[i*K + k] /= sum;

	}
	stag->ll = ll;
	return err;
}/* generate_E */


