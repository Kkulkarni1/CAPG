/**
 * @file run_ampliclust.c
 * @author Karin S. Dorman
 *
 * Cluster amplicon sequences.
 *
 * TODO
 * - debug: why does partition_model_free()'s method of initializing clusters not produce the same result as equivalent call to kmodes_macqueen (maybe b/c log(error_prob)/3 intead of log(error_prob/3), but not tested)
 * - initialization clusters haplotypes of length \par data.fdata.n_min_length[initialization.c]
 *   initializes remaining positions in a foolish way [initialization.c]
 * - handle error indel in reads: very rare in Illumina
 * - unsigned char data_t (previous) VS uint8_t data_t (current);
 *
 * DONE
 * X initialization is unaware of \par options.background_model	[initialization.c]
 * X rewrite cluster_amplicons() so that writing a file of results is optional
 * X created initializer object to store initialization information separate from data and model
 * X change default to NOT use curses: use -w to turn on curses
 * X split into separate files for options, model, data, initalization, simulation, etc.
 * X take a set of seeds to initialize (-i <seed_file_as_fasta>)
 * X estimate error rates (& compare to estimated from true error data)
 * X allow different length reads
 * X rewrite message handling to use global debug level
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <curses.h>

#include "ampliclust.h"
#include "run_ampliclust.h"
#include "aecm.h"
#include "data.h"
#include "model.h"
#include "options.h"
#include "initialize.h"
#include "fastq.h"
#include "simulate.h"
#include "io.h"
#include "stages.h"
#include "cluster.h"
#include "statistics.h"
#include "math.h"

int make_run_info(run_info **ri, data *dat, simulator *sim, options *opt);
int do_cluster_amplicons(data *dat, model *mod, initializer *ini, run_info *ri);
void compute_differences(data *dat, model *mod, model_options *opt, double *dpi, double *ddelta, double *dgamma, double *dlambda0, double *dlambda1, double *dbeta, double *davg_mhd);
int modified_ic(unsigned char* hap, unsigned char *est_anc, double *distance, double best_ll, unsigned int K, double *JC_ll, double *n_aic, double *n_bic, unsigned int n_param, unsigned int max_read_length,size_t sample_size);
int write_results(char const * const outfile, data *dat, model *mod, options *opt, initializer *ini, simulator *sim, void *ri, int which);


int main(int argc, const char **argv)
{
	int err = NO_ERROR;		/* error state */
	options *opt = NULL;		/* run options */
	data *dat = NULL;		/* data object */
	model *mod = NULL;		/* estimation model, fitted by AECM */
	model *sim_mod = NULL;		/* simulation model */
	stages *stag = NULL;       	/* multi_stages clustering */
	fastq_options *fqo = NULL;	/* fastq file options */
	initializer *ini = NULL;	/* initializer */
	initializer *sim_ini = NULL;	/* simulation initializer */
	simulator *sim = NULL;		/* simulator */
	run_info *ri = NULL;		/* run_info object */

	/* parse command line */
	if ((err = make_options(&opt))
		|| (err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;

	/* make data object */
	if ((err = make_data(&dat, opt)))
		goto CLEAR_AND_EXIT;

	/* setup fastq options */
	if ((err = make_fastq_options(&fqo)))
		goto CLEAR_AND_EXIT;

	/* encode nucleotides in 2-bits: error raised if ambiguous bases
	 * NOTE:  changing this has implications because we assume data.dmat
	 * and model.haplotypes contain xy_t type data.
	 */
	fqo->read_encoding = XY_ENCODING;

	/* read sequence data */
	if (opt->fastq_file && (err = read_fastq(opt->fastq_file,
						&dat->fdata, fqo)))
		goto CLEAR_AND_EXIT;

	/* or allocate space for data */
	else if (!opt->fastq_file && opt->simo->sim_n_reads && (err =
			allocate_empty_fastq(&dat->fdata, fqo,
			opt->simo->sim_n_reads, opt->simo->sim_read_length)))
		goto CLEAR_AND_EXIT;

	/* with data loaded, can polish off data object,
	 * although indexing of the reads is postponed
	 * if they are still to be simulated
	 */
	if ((err = load_data(dat, opt, opt->do_simulation
		&& opt->simo->simulation_random & NUCLEOTIDES)))
		goto CLEAR_AND_EXIT;

/* [WORKING] data::n_quality set to range of quality scores observed in loaded data */
debug_msg(1, 1, "data::n_quality: %u\n", dat->n_quality);

	/* simulate data */
	if (opt->do_simulation) {

		/* make simulation model: synchronize now */
		if ((err = make_model(&sim_mod, dat, opt->sim_modo, 0)))
			goto CLEAR_AND_EXIT;
/* [WORKING] model::n_quality set to range of quality scores observed (or predicted) in loaded data or to the number of qualities use in error profile profile; if data loaded, model::qual_list as per data and model::n_param counted */

		if ((err = make_simulator(&sim, dat, sim_mod, opt->simo)))
			goto CLEAR_AND_EXIT;

		/* make initializer for simulator: synchronize now */
		if ((err = make_initializer(&sim_ini, dat, sim_mod,
			opt->sim_inio, 0)))
			goto CLEAR_AND_EXIT;

		if ((err = do_simulation(sim, sim_mod, dat, sim_ini, opt->simo)))
			goto CLEAR_AND_EXIT;

		if (!sim_mod->synced && (err = sync_model(sim_mod, dat, opt->sim_modo)))
			goto CLEAR_AND_EXIT;
debug_msg(1, 1, "data::n_quality: %u\n", dat->n_quality);
/* [WORKING] data::n_quality may be updated and if it exceeds model::n_quality, update that and memory allocated for quality-related parameters in model */

		/* compute loglikelihood with true parameteris if possible */
		sim_mod->ll = e_step(dat, sim_mod, opt->sim_modo, FIRST);
		debug_msg(1, 0, "Loglikelihood under true model: %f\n",
								sim_mod->ll);

		/* write parameters and true membership */
		if (opt->simo->simulation_outfile)
			write_results(opt->simo->simulation_outfile, dat,
				sim_mod, opt, sim_ini, sim, NULL, TRUE_VALUES);

		if (sim_ini)
			free_initializer(sim_ini, opt->inio);

		/* data available: finish some setup */
		if ((err = sync_data(dat)))
			goto CLEAR_AND_EXIT;
/*
		if (!mod->synced && (err = sync_model(mod, dat, opt->modo)))
			goto CLEAR_AND_EXIT;
		if (!ini->synced
			&& (err = sync_initializer(ini, dat, opt->inio)))
			goto CLEAR_AND_EXIT;
*/
	}

	if (opt->do_estimation) {

		/* create estimation model */
		if ((err = make_model(&mod, dat, opt->modo, 0)))
			goto CLEAR_AND_EXIT;

		/* [TODO] allow haplotype indels: not implemented */
		if (opt->reference_file && (err = pw_align_reads(dat->fdata,
							 opt->reference_file)))
			goto CLEAR_AND_EXIT;

		/* make initializer */
		if ((err = make_initializer(&ini, dat, mod, opt->inio, 0)))
			goto CLEAR_AND_EXIT;

		/* create run_info object */
		if ((err = make_run_info(&ri, dat, sim, opt)))
			goto CLEAR_AND_EXIT;

		/* main algorithm(s) */

		/* multi-stage method */
		if (opt->multi_stage_method && (dat->fdata->n_reads > opt->sample_size )) {
			/* create stages object */
			if ((err = make_stages(&stag, dat, opt)))
				goto CLEAR_AND_EXIT;

			if ((err = multi_stage(mod, dat, opt, ini, sim, ri, stag)))
				goto CLEAR_AND_EXIT;
		/* estimate K */
		} else if (opt->estimate_K) {
			if ((err = ampliclustK(dat, mod, opt, ini, ri)))
				goto CLEAR_AND_EXIT;

		/* estimate model given K */
		} else if ((err = do_cluster_amplicons(dat, mod, ini, ri))) {
			goto CLEAR_AND_EXIT;
		}
	} else if (opt->outfile) {
		write_results(opt->outfile, dat, NULL, opt, NULL, NULL, NULL, 0);
	}


CLEAR_AND_EXIT:

	if (dat)
		free_data(dat);
	if (ini)
		free_initializer(ini, opt->inio);
	if (mod)
		free_model(mod);
	if (ri)
		free_run_info(ri);
	if (opt)
		free_options(opt);
	if (stag)
		free_stages(stag);
	if (sim_mod)
		free_model(sim_mod);

	if (fqo)
		free_fastq_options(fqo);

	return (EXIT_FAILURE);
} /* main */


/**
 * Cluster the amplicons while interacting with user.
 *
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param ini	pointer to initializion object
 * @param ri	pointer to run information object
 * @return	error status
 */
int do_cluster_amplicons(data *dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;

	/* optionally start curses interface to allow interruption of AECM */
#ifdef USE_CURSES
	if (ri->opt->use_curses) {
		global_wp = initscr();
		noecho();
		timeout(0);
		scrollok(global_wp, TRUE);
		cbreak();
	}
#endif

	if ((err = cluster_amplicons(dat, mod, ini, ri->opt->modo,
			ri->opt->inio, do_per_iterate, do_initialize,
			(void *) ri)))
		return err;

#ifdef USE_CURSES
	if (ri->opt->use_curses) {
		echo();
		endwin();
	}
#endif

	/* though AECM not run, can still output results from initialization */
	if (ri->opt->inio->n_init == 0 && ri->opt->do_simulation) {

		ri->opt->inio->initialization_method = INIT_TRUE_PARTITION;
		if ((err = do_initialize(dat, mod, ini, (void *) ri)))
			return err;
	//	mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->ll);

		/* \par dat.best_modes and \par dat.best_cluster_id set and
		 * can be used to compute log likelihood under agnostic model
		 * or ampliclust model with initial parameter estimates */
		if (ri->opt->inio->initialization_method == INIT_IMDEBLUR) {
			/* initializer computes mod->ll and mod->K */
			mod->n_param =  ri->opt->modo->K * dat->max_read_length
							+ ri->opt->modo->K - 1;
			assign_clusters(mod->eik, ri->opt->modo->K,
				dat->sample_size, ini->best_cluster_size,
				ini->best_cluster_id, 1);
		} else {
			mod->ll = e_step(dat, mod, ri->opt->modo, FIRST);
			assign_clusters(mod->eik, mod->n_mix, dat->sample_size,
				ini->best_cluster_size, ini->best_cluster_id, 0);//
		}
	//	mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->ll);
		param_update(mod, dat, ri->opt->modo, BEST_VALUES);
		mod->best_ll = mod->ll;

		memcpy(ri->optimal_cluster_id, ini->best_cluster_id,
				dat->sample_size * sizeof *ini->cluster_id);
		memcpy(ri->optimal_cluster_size, ini->best_cluster_size,
				mod->n_mix * sizeof *ini->cluster_size);

		/* and the seeds that beget this great solution */
		if (ri->opt->inio->initialization_method != INIT_TRUE_PARTITION
			&& ri->opt->inio->initialization_method != INIT_TRUE_VALUES)
			memcpy(ri->optimal_seed_idx, ini->seed_idx,
				ri->opt->modo->K * sizeof *ini->seed_idx);

		/* calculate aic and bic */
		if (ri->opt->modo->JC69_model) {
			modified_ic(mod->haplotypes, mod->est_ancestor,
				mod->distance, mod->best_ll, ri->opt->modo->K,
				&mod->JC_ll, &mod->aic, &mod->bic, mod->n_param,
				dat->max_read_length, dat->sample_size);
		} else {
				mod->aic = aic(mod->best_ll, mod->n_param);
				mod->bic = bic(mod->best_ll, mod->n_param,
						dat->sample_size);
		}
	//	mmessage(INFO_MSG, NO_ERROR, "ll: %f.\n", mod->best_ll);
	//	mmessage(INFO_MSG, NO_ERROR, "bic : %f. aic: %f\n", bic, aic);

		if (ri->opt->outfile && !ri->opt->estimate_K)
			write_results(ri->opt->outfile, dat, mod, ri->opt, ini,
							ri->sim, ri, BEST_VALUES);
	}

	free_cluster_statics();

	return err;
} /* do_cluster_amplicons */


/**
 * A call_per_iterate() function.  Produce the summary displays that go to 
 * screen when running ampliclust.
 *
 * @param dat	pointer to data object
 * @param mod	pointer to model object
 * @param ini	pointer to initializer object
 * @param obj	hidden run_info object
 */
void do_per_iterate(data *dat, model *mod, initializer *ini, void *obj)
{
	run_info *ri = (run_info *) obj;
	options *opt = ri->opt;
	simulator *sim = ri->sim;
	double dpi = -1, dgamma = -1;           /* variables to report status */
	double dlambda0 = -1, dlambda1 = -1, ddelta = -1, dbeta = -1;
	double ari = 0, davg_mhd = -1;
	int better = 0;

	/* given previous stored solution, calculate differences */
	if (isfinite(mod->best_ll) || isfinite(opt->modo->previous_ll))//
		compute_differences(dat, mod, opt->modo, &dpi, &ddelta,
			&dgamma, &dlambda0, &dlambda1, &dbeta, &davg_mhd);

	if (mod->ll <= mod->best_ll) {

		assign_clusters(mod->eik, mod->n_mix, dat->sample_size,
				ini->cluster_size, ini->cluster_id, 0);
		ari = cluster_index(ini->cluster_id, ri->optimal_cluster_id,
			 dat->sample_size, mod->n_mix, mod->n_mix,
						ADJUSTED_RAND_INDEX);
	/* better solution */
	} else {

		better = 1;

		mod->ll = e_step(dat, mod, opt->modo, FIRST);
		assign_clusters(mod->eik, mod->n_mix, dat->sample_size,
			ini->best_cluster_size, ini->best_cluster_id, 0);
		ari = cluster_index(ini->best_cluster_id, ri->optimal_cluster_id,
			 dat->sample_size, mod->n_mix, mod->n_mix,
						ADJUSTED_RAND_INDEX);


/* verbose output I now find to wordy (holding for debugging)
	if (opt->do_simulation) {
		model *smod = obj->sim_mod;
		PRINT(global_wp, "Comparing to truth:\n");
		PRINT(global_wp, "ARI: %f\n",
			cluster_index(ini->best_cluster_id,
			sim->true_cluster_id, dat->sample_size,
			mod->n_mix, mod->n_mix, ADJUSTED_RAND_INDEX));
		PRINT(global_wp, "\\Delta pi: %f\n", euc_dis(
			mod->pi, smod->pi, mod->n_mix)
			/ mod->n_mix);
		PRINT(global_wp, "\\Delta delta: %f\n", euc_dis(
			mod->delta, smod->delta, dat->max_read_position)
			/ dat->max_read_position);
		PRINT(global_wp, "\\Delta gamma: %f\n", euc_dis(
			mod->gamma, smod->gamma,
			NUM_NUCLEOTIDES*NUM_NUCLEOTIDES)/(NUM_NUCLEOTIDES*(NUM_NUCLEOTIDES-1)));
		PRINT(global_wp, "haplotype distances (%u x %u):\n", opt->modo->K, opt->true_K);
		for (size_t j = 0; j < opt->modo->K; ++j) {
			for (size_t k = 0; k < opt->true_K; ++k)
				PRINT(global_wp, " %2u",
				hamming_char_dis(
				&mod->haplotypes[j*dat->max_read_position],
				&smod->haplotypes[k*dat->max_read_position],
				dat->max_read_position));
			PRINT(global_wp, "\n");
		}
	} else {
		PRINT(global_wp, "Comparing to previous best:\n");
		PRINT(global_wp, "ARI: %f\n", ari);
		PRINT(global_wp, "\\Delta pi: %f", dpi);
		PRINT(global_wp, "; delta: %f", ddelta);
		PRINT(global_wp, "; gamma: %f", dgamma);
		PRINT(global_wp, "; lambda0: %f", dlambda0);
		PRINT(global_wp, "; lambda1: %f\n", dlambda1);
		if (mod->best_ll > opt->modo->previous_ll) {
			PRINT(global_wp, "haplotype distances:\n");
			for (size_t j = 0; j < opt->modo->K; ++j) {
				for (size_t k = 0; k < opt->modo->K; ++k)
					PRINT(global_wp, " %2u",
					hamming_char_dis(
					&mod->haplotypes[j*dat->max_read_position],
					&mod->best_haplotypes[k*dat->max_read_position],
					dat->max_read_position));
				PRINT(global_wp, "\n");
			}
		}
	}
*/

		/* store newly identified MLEs in model.best_* slots */
		param_update(mod, dat, opt->modo, OPTIMAL_VALUES);

		/* store newly identified hard clustering solution */
		memcpy(ri->optimal_cluster_id, ini->best_cluster_id,
			dat->sample_size * sizeof *ini->cluster_id);
		memcpy(ri->optimal_cluster_size, ini->best_cluster_size,
			mod->n_mix * sizeof *ini->cluster_size);

		/* and the seeds that beget this great solution */
		if (opt->inio->initialization_method != INIT_TRUE_PARTITION
			&& opt->inio->initialization_method != INIT_TRUE_VALUES)
			memcpy(ri->optimal_seed_idx, ini->seed_idx,
				opt->modo->K * sizeof *ini->seed_idx);

	/* output information about better solution */
/* again: I find this output too wordy now
	PRINT(global_wp, "assignments: ");
	wprint_uints(global_wp, ri.optimal_cluster_id,
		dat->sample_size, 2, 0);
	PRINT(global_wp, "pi: ");
	if (global_wp)
		wprint_doubles(global_wp, mod->best_pi, mod->n_mix, 3, 1);
	else
		fprint_doubles(stderr, mod->best_pi, mod->n_mix, 3, 1);
	wprintw(global_wp, "delta: ");
	wprint_doubles(global_wp, mod->best_delta, dat->max_read_position, 3, 1);
	wprintw(global_wp, "gamma:\n");
	wprint_vectorized_sq_matrix(global_wp, mod->best_gamma, NUM_NUCLEOTIDES,
								ROW_ORDER);
	wprint_fasta(global_wp, mod->best_haplotypes, opt->modo->K,
				dat->max_read_position, "H");
	PRINT(global_wp, "sizes: ");
	if (global_wp)
		wprint_size_t(global_wp, ri.optimal_cluster_size,
			mod->n_mix, 3, 1);
	else
		fprint_size_ts(stderr, ri.optimal_cluster_size, mod->n_mix,
									3, 1);
*/

		/* for model comparison, calculate aic and bic */
		if (opt->modo->JC69_model) {
			modified_ic(mod->haplotypes, mod->est_ancestor, mod->distance,
				mod->best_ll, mod->K, &mod->JC_ll, &mod->aic, &mod->bic,
				mod->n_param, dat->max_read_length, dat->sample_size);
		} else {
			mod->aic = aic(mod->best_ll, mod->n_param);
			mod->bic = bic(mod->best_ll, mod->n_param, dat->sample_size);
		}

		/* overwrite output file with this solution */
		if (opt->outfile && !opt->estimate_K)
			write_results(opt->outfile, dat, mod, opt,
				ini, sim, ri, OPTIMAL_VALUES);
	}

	/* friendly output about the quality of current solution */
	cc_msg(global_wp, SILENT, SILENT, " %3u", ini->n_inits);
	PRINT(global_wp, " | pi=%4.2f", dpi);
	if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION) {
		PRINT(global_wp, " beta=%4.2f", dbeta);
	} else if ((opt->modo->parameterization == ART_PARAMETERIZATION
						&& opt->modo->art_separate)
		|| opt->modo->parameterization == QUALITY_PARAMETERIZATION) {
		PRINT(global_wp, " gamma=%4.2f", dgamma);
	} else {
		PRINT(global_wp, " delta=%4.2f", ddelta);
		PRINT(global_wp, " gamma=%4.2f", dgamma);
		PRINT(global_wp, " lambda0=%4.2f", dlambda0);
		PRINT(global_wp, " lambda1=%4.2f", dlambda1);
	}
	PRINT(global_wp, " ari=%5.3f", ari);
	PRINT(global_wp, " ahd=%4.2f", davg_mhd);
	PRINT(global_wp, " |");
	/* output: cluster sizes */
#ifdef USE_CURSES
	if (global_wp)
		wprint_uints(global_wp, better ? ri->optimal_cluster_size
			: ini->cluster_size, mod->n_mix, 3, 0);
	else
#endif
		fprint_uints(stderr, better ? ri->optimal_cluster_size
			: ini->cluster_size, mod->n_mix, 3, 0);
	/* initialization criterion (best) & log likelihood (best) */
	PRINT(global_wp, " | %.0f (%.0f) | ll=%.3f (ll=%.3f)\n",
		ini->best_total, ini->optimal_total, mod->ll,
		mod->best_ll);

} /* do_per_iterate */


/**
 * Do initialization, either from truth or by pulling on bootstraps.
 *
 * @param dat	data object
 * @param mod	model object
 * @param ini	initialization object
 * @param obj	options object
 * @return	error status
 */
int do_initialize(data *dat, model *mod, initializer *ini, void *obj)
{
	run_info *ri = (run_info *) obj;

	if (ri->opt->inio->initialization_method == INIT_TRUE_VALUES)
		return initialize_from_true_parameters(mod, ri->sim->mod, dat,
						ri->opt->modo, ri->opt->simo);
	else if (ri->opt->inio->initialization_method == INIT_TRUE_PARTITION)
		return initialize_from_true_partition(mod, dat, ini, ri->sim,
						ri->opt->simo, ri->opt->inio);
	else
		return initialize(mod, dat, ini, ri->opt->inio);
} /* do_initialize */


/**
 * Compute distances between previous best solution and current solution.
 *
 * [TODO] The average Hamming distance computed between current haplotypes
 * and previous best haplotypes is not quite right since there is no
 * correct way to match the haplotypes.  It currently just finds the HD
 * to the closest haplotype, which means that multiple haplotypes could
 * be mapped to the same haplotype.  Thus the average HD is underestimated.
 *
 * @param mod		model pointer
 * @param dat		data pointer
 * @param opt		model_options pointer
 * @param dpi		to store Euclidean difference in pi parameter
 * @param ddelta	to store Euclidean distance in delta parameter
 * @param dgamma	to store Euclidean distance in gamma parameter
 * @param dlambda0	to store Euclidean distance in lambda0 parameter
 * @param dlambda1	to store Euclidean distance in lambda1 parameter
 * @param dbeta		to store Euclidean distance in beta parameter
 * @param davg_mhd	to store average Hamming distance among haplotypes
 * @return		return values are in the last 6 arguments
 */
void compute_differences(data *dat, model *mod, model_options *opt, double *dpi,
	double *ddelta, double *dgamma, double *dlambda0, double *dlambda1,
	double *dbeta, double *davg_mhd)
{
	*dpi = euc_dis(mod->pi, mod->best_pi, mod->n_mix) / mod->n_mix;
	if (opt->parameterization == MLOGIT_PARAMETERIZATION) {
		*dbeta = euc_dis(mod->beta, mod->best_beta,
			mod->n_beta_coef)/mod->n_beta_coef;
	} else if ((opt->parameterization == ART_PARAMETERIZATION
						&& opt->art_separate)
		|| opt->parameterization == QUALITY_PARAMETERIZATION) {
		*dgamma = frobenius_norm(mod->gamma, mod->best_gamma,
			NUM_NUCLEOTIDES, 0)
			/ (NUM_NUCLEOTIDES*(NUM_NUCLEOTIDES-1));
	} else {
		*ddelta = euc_dis(mod->delta, mod->best_delta,
			dat->max_read_position)/dat->max_read_position;
		*dgamma = frobenius_norm(mod->gamma, mod->best_gamma,
			NUM_NUCLEOTIDES, 0)
			/ (NUM_NUCLEOTIDES*(NUM_NUCLEOTIDES-1));
		*dlambda0 = euc_dis(mod->lambda0, mod->best_lambda0,
			dat->max_read_position * mod->n_quality)
			/(dat->max_read_position * mod->n_quality);
		*dlambda1 = euc_dis(mod->lambda1, mod->best_lambda1,
			dat->max_read_position * mod->n_quality)
			/(dat->max_read_position * mod->n_quality);
	}
	*davg_mhd = 0;
	for (size_t j = 0; j < opt->K; ++j) {
		size_t nl = dat->max_read_length;
		double min_hd = nl;
		for (size_t k = 0; k < opt->K; ++k) {
			double hd = hamming_char_dis((char *)
				&mod->haplotypes[j*nl], (char *)
				&mod->best_haplotypes[k*nl], nl);
			if (hd < min_hd)
				min_hd = hd;
		}
		*davg_mhd += min_hd;
	}
	*davg_mhd /= opt->K;
} /* compute_differences */


/**
 * Cluster amplicon reads without a given K
 *
 * @param dat	data object
 * @param best_mod	model object (best model identified below )
 * @param opt	options object
 * @param ini   initializer object
 * @param ri    run_info object
 * @return	error status
 */
int ampliclustK(data *dat, model *best_mod, options *opt, initializer *ini,
	run_info *ri)
{
	int err = NO_ERROR;
	unsigned int K_ini = opt->modo->min_K;
	opt->modo->K = opt->modo->min_K;

	/* model object mod for storing best result under different K */
	model *mod = NULL;
	if ((err = make_model(&mod, dat, opt->modo, 0)))
		return err;

	/* Set a lower bound when finding K */
	/* TODO set up the lowerbound when initialize the data object */
	unsigned int abun_count = (unsigned int) dat->sample_size * 0.01;
	unsigned int lowerbound = count_sequences(dat->seq_count,
		abun_count > 10? abun_count:10);
	// lowerbound = 0;

	/* use lowerbound only if the min_K is the default value */
	if (lowerbound > opt->modo->min_K && opt->modo->min_K == 2) {
		K_ini = lowerbound;
		opt->modo->K = lowerbound;
		if ((err = realloc_model(mod, dat, opt->modo)))
			return err;
		if ((err = realloc_initializer(ini, dat, opt->inio)))
			return err;
		if ((err = realloc_run_info(ri, dat->sample_size,
			 opt->modo->K + opt->modo->background_model)))
			return err;
	} else if(opt->multi_stage_method) {

		if ((err = realloc_model(best_mod, dat, opt->modo)))
			return err;
		if ((err = realloc_initializer(ini, dat, opt->inio)))
			return err;
		if ((err = realloc_run_info(ri, dat->sample_size,
				opt->modo->K + opt->modo->background_model)))
			return err;
	}

	for (unsigned int k = K_ini; k <= opt->modo->max_K; ++k) {
		opt->modo->K = k;

		/*
		When increasing k ,maybe we should also increase number of initialization.
		Suppose the number of initialization required is alpha * k * n * p
		alpha to be determined here
		*/
		/* If number of initialization is not given in the command */
		/* TODO may think a better way to do this */
		/*
		if(opt->inio->n_init == 1)
			opt->inio->n_init = (unsigned int) 0.002 * k * dat->max_read_position * dat->sample_size;
		*/

		/* If it is not the first time of the loop or we change the opt->min_K, we need reallocate some objects based on different K */
		if (k != K_ini) {
			if ((err = realloc_model(mod, dat, opt->modo)))
				return err;
			if ((err = realloc_initializer(ini, dat, opt->inio)))
				return err;
			if ((err = realloc_run_info(ri, dat->sample_size,
				opt->modo->K + opt->modo->background_model)))
				return err;
		}

		/* mod is the pointer to the best model under k, returned from cluster_amplicons() */
		if ((err = cluster_amplicons(dat, mod, ini, opt->modo,
			opt->inio, do_per_iterate, do_initialize, (void *) opt)))
			return err;

		mmessage(INFO_MSG, NO_ERROR, "mod.aic:%f; mod.bic:%f.\n",
			mod->aic, mod->bic);

		/* currently, we use bic or aic to do model selection   */

		if ((opt->modo->use_aic ? mod->aic : mod->bic)
			< (opt->modo->use_aic ? best_mod->aic : best_mod->bic)) {

			/* check if we need reallocate struct best_mod for a different K */
			if (best_mod->K != mod->K)
				if ((err = realloc_model(best_mod, dat, opt->modo)))
					return err;

			/* copy current optimal model */
			if ((err = copy_model(best_mod, mod, dat, 1)))
				return err;
		}

		/* Maybe the previous K is a better solution */
		if ((opt->modo->use_aic ? mod->aic : mod->bic)
			> (opt->modo->use_aic ? best_mod->aic : best_mod->bic)) {
			/* If we do not choose an appropriate lowbound for K, we may have problem here */
			if (k == K_ini+1) {
				mmessage(INFO_MSG, NO_ERROR, "The best K is "
					"unlikely in the range (%i,%i] "
					"and the algorithm will be stopped and output the best model K = %i.\n",
					K_ini, opt->modo->max_K, K_ini);
				break;
			} else {
				mmessage(INFO_MSG, NO_ERROR, "The best K is %i.\n",
					best_mod->K);
				break;
			}
		} else if (k == opt->modo->max_K) {
			mmessage(INFO_MSG, NO_ERROR, "The best K has "
				"not been detected in the range [%i,%i]"
				" and you may choose a larger range.\n",
				K_ini, opt->modo->max_K);
			break;
		}

		mmessage(INFO_MSG, NO_ERROR, "Best_mod.aic:%f; Best_mod.bic:%f.\n",
			best_mod->aic, best_mod->bic);

	}

	/* reassign reads under the best model */
	if ((err = realloc_run_info(ri, dat->sample_size, best_mod->n_mix)))
		return err;

	assign_clusters(best_mod->eik, best_mod->n_mix, dat->sample_size,
		ri->optimal_cluster_size, ri->optimal_cluster_id, 0);

	/* print the result of the best model */

	/* just copy and modify part of write_result().
	rewrite the below later */
	/* TODO move them to one appropriate function later */

	if (!opt->multi_stage_method ) {
		FILE *fp = fopen(opt->outfile, "w");
		if (!fp)
			return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR,
				opt->outfile);

		fprintf(fp, "best_K: %i\n", best_mod->K);
		if (opt->modo->model_quality) {
			fprintf(fp, "lambda0:");
			fprint_vectorized_matrix(fp, best_mod->best_lambda0,
			 	dat->max_read_position, best_mod->n_quality, 1);
			fprintf(fp, "lambda1:");
			fprint_vectorized_matrix(fp, best_mod->best_lambda1,
				dat->max_read_position, best_mod->n_quality, 1);
		}
		fprintf(fp, "pi: ");
		fprint_doubles(fp, best_mod->best_pi, best_mod->K, 6, 1);
		fprintf(fp, "delta: ");
		fprint_doubles(fp, best_mod->best_delta, dat->max_read_position,
			6, 1);
		fprintf(fp, "gamma: ");
		fprint_vectorized_sq_matrix(fp, best_mod->best_gamma,
							NUM_NUCLEOTIDES, 1);
		if (opt->modo->background_model) {
			fprintf(fp, "bg_pi: ");
			fprint_doubles(fp, best_mod->best_bg_pi,
				NUM_NUCLEOTIDES, 6, 1);
			fprintf(fp, "bg_lambda: ");
			fprint_vectorized_matrix(fp, best_mod->best_bg_lambda,
				 dat->max_read_position, best_mod->n_quality,
									 1);
		}
		fprintf(fp, "assignments: ");
		fprint_uints(fp, ri->optimal_cluster_id, dat->sample_size,
								2, 1);
		fprintf(fp, "sizes: ");
		fprint_uints(fp, ri->optimal_cluster_size, best_mod->n_mix, 3, 1);

		fprint_fasta(fp, best_mod->best_haplotypes, best_mod->K,
					 dat->max_read_length, "H");
		if (opt->modo->JC69_model){
			fprintf(fp, "Estimated common ancestor: \n");
			fprint_fasta(fp, best_mod->est_ancestor, 1,
					 dat->max_read_length, "A");
			fprintf(fp, "Evolution_rate: ");
			fprint_doubles(fp, best_mod->distance, best_mod->K,
				6, 1);
			fprintf(fp, "log likelihood from JC69 model: %f\n",
							best_mod->JC_ll);
		}
		fprintf(fp, "log likelihood: %f\n", best_mod->best_ll);
		fprintf(fp, "aic: %f\n", best_mod->aic);
		fprintf(fp, "bic: %f\n", best_mod->bic);

		fclose(fp);
	}

	if(mod)
		free_model(mod);

	return err;
}/* ampliclust_K */

/**
 * Write results from both model and data to file.
 *
 * @param outfile	name of file
 * @param dat		data object pointer
 * @param mod		model object pointer
 * @param opt		options object pointer
 * @param sim		simulator object pointer
 * @param ini		initializer object pointer
 * @param vri		obscured run_info object pointer
 * @param which		which results to output (OPTIMAL_VALUES|BEST_VALUES|TRUE_VALUES)
 * 	OPTIMAL_VALUES	after running AECM from an initialization or just initializing
 * @return		error status
 */
int write_results(char const * const outfile, data *dat, model *mod,
	options *opt, initializer *ini, simulator *sim, void *vri, int which)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "entering\n");

	if (which == OPTIMAL_VALUES && !vri)
		return MESSAGE(global_wp, ERROR_MSG, INTERNAL_ERROR,
			"missing run_info object");

	FILE *fp = fopen(outfile, "w");
	if (!fp)
		return MESSAGE(global_wp, ERROR_MSG, FILE_OPEN_ERROR, outfile);

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "file '%s' opened\n", outfile);

	fprintf(fp, "rng seed: %lu\n", opt->seed);
	fprintf(fp, "quality range: %u %u\n", dat->fdata->min_quality,
						dat->fdata->max_quality);

	fprintf(fp, "expected number of errors:");
	for (size_t i = 0; i < dat->sample_size; ++i) {
		double eecnt = 0;
		for (size_t j = 0; j < dat->lengths[i]; ++j)
			eecnt += error_prob(dat->fdata, dat->qmat[i][j]);
		fprintf(fp, " %.3f", eecnt);
	}
	fprintf(fp, "\n");

	/* simulation data: \ref simulator.true_cluster_(id|size) are truth */
	if (which == TRUE_VALUES) {
		fprintf(fp, "pi (truth): ");
		/* [TODO] simulation cannot include background model */
		fprint_doubles(fp, sim->mod->pi, opt->simo->true_K, 6, 1);
		if (opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION) {
			fprintf(fp, "delta (truth): ");
			fprint_doubles(fp, sim->mod->delta,
						dat->max_read_position, 6, 1);
			fprintf(fp, "lambda0 (truth):");
			fprint_vectorized_matrix(fp, sim->mod->lambda0,
				dat->max_read_position, mod->n_quality, 1);
			fprintf(fp, "lambda1 (truth):");
			fprint_vectorized_matrix(fp, sim->mod->lambda1,
				dat->max_read_position, mod->n_quality, 1);
		} else if (opt->sim_modo->parameterization != ART_PARAMETERIZATION
				&& opt->simo->simulation_random & QUALITIES) {
			fprintf(fp, "lambda (truth):");
			fprint_vectorized_matrix(fp, sim->mod->lambda,
				dat->max_read_position, mod->n_quality, 1);
		}
		if (mod->mut_position) {
			fprintf(fp, "mutable positions:");
			fprint_uints(fp, mod->mut_position,
						dat->max_read_length, 2, 1);
		}
		if (opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION
			|| opt->sim_modo->parameterization == ART_PARAMETERIZATION
			|| opt->sim_modo->parameterization == QUALITY_PARAMETERIZATION) {
			fprintf(fp, "gamma (truth): ");
			fprint_vectorized_sq_matrix(fp, sim->mod->gamma,
							NUM_NUCLEOTIDES, 1);
		} else if (opt->sim_modo->parameterization == MLOGIT_PARAMETERIZATION) {
			fprintf(fp, "beta (truth): ");
			fprint_doubles(fp, sim->mod->beta,
						mod->n_beta_coef, 6, 1);
		}
		fprintf(fp, "True number of errors:");
		for (size_t i = 0; i < dat->sample_size; ++i) {
			unsigned int cnt = 0;
			for (size_t j = 0; j < dat->lengths[i]; ++j)
				if (dat->dmat[i][j] != sim->mod->haplotypes[
					sim->true_cluster_id[i]*dat->max_read_position
					+ j + (dat->offset ? dat->offset[i] : 0)])
					++cnt;
			fprintf(fp, " %u", cnt);
		}
		fprintf(fp, "\n");
		fprintf(fp, "assignments (truth): ");
		fprint_uints(fp, sim->true_cluster_id, dat->sample_size, 2, 1);
		/* [TODO] simulation cannot include background model */
		fprintf(fp, "sizes (truth): ");
		fprint_uints(fp, sim->true_cluster_size, opt->simo->true_K, 3, 1);
		fprint_fasta(fp, sim->mod->haplotypes, opt->simo->true_K,
			dat->max_read_position, "True H");
		fprintf(fp, "log likelihood (true values): %f\n", sim->mod->ll);
		fprintf(fp, "aic (true values): %f\n", aic(sim->mod->ll, sim->mod->n_param));
		fprintf(fp, "bic (true values): %f\n", bic(sim->mod->ll, sim->mod->n_param,
						dat->sample_size));

		fclose(fp);
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "file '%s' closed\n", outfile);

		goto WRITE_RESULTS_RETURN;
	}

	if (!opt->inio->n_init)
		goto WRITE_RESULTS_RETURN;

	/* now ini is defined */

	fprintf(fp, "initialization: %u\n", ini->n_inits);

	run_info *ri = (run_info *) vri;

	size_t *osi = which == OPTIMAL_VALUES ? ri->optimal_seed_idx
						: ini->best_seed_idx;
	fprintf(fp, "seeds: ");
	fprint_size_ts(fp, osi, opt->modo->K, 4, 1);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "seeds output\n");

	double *pi = which == OPTIMAL_VALUES ? mod->best_pi : mod->pi;
	if (opt->modo->param_estimate & PARAM_PI) {
		fprintf(fp, "pi: ");
		fprint_doubles(fp, pi, opt->modo->K, 6, 1);
	}

	double *beta = which == OPTIMAL_VALUES ? mod->best_beta : mod->beta;
	double *lambda1 = which == OPTIMAL_VALUES ? mod->best_lambda1 : mod->lambda1;
	double *lambda0 = which == OPTIMAL_VALUES ? mod->best_lambda0 : mod->lambda0;
	double *delta = which == OPTIMAL_VALUES ? mod->best_delta : mod->delta;
	double *gamma = which == OPTIMAL_VALUES ? mod->best_gamma : mod->gamma;
	double *bg_pi = which == OPTIMAL_VALUES ? mod->best_bg_pi : mod->bg_pi;
	double *bg_lambda = which == OPTIMAL_VALUES ? mod->best_bg_lambda
							: mod->bg_lambda;
	if (opt->modo->parameterization == MLOGIT_PARAMETERIZATION) {
		if (opt->modo->param_estimate & PARAM_BETA) {
			fprintf(fp, "beta:");
			fprint_doubles(fp, beta, mod->n_beta_coef, 6, 1);
		}
	} else {
		if (opt->modo->model_quality &&
			opt->modo->parameterization == DEFAULT_PARAMETERIZATION
				&& opt->modo->param_estimate & PARAM_LAMBDA) {
			fprintf(fp, "lambda0:");
			fprint_vectorized_matrix(fp, lambda0,
				dat->max_read_position, mod->n_quality, 1);
			fprintf(fp, "lambda1:");
			fprint_vectorized_matrix(fp, lambda1,
				dat->max_read_position, mod->n_quality, 1);
		}
		if (opt->modo->parameterization == DEFAULT_PARAMETERIZATION
				&& opt->modo->param_estimate & PARAM_DELTA) {
			fprintf(fp, "delta: ");
			fprint_doubles(fp, delta, dat->max_read_position, 6, 1);
		}
		if (opt->modo->param_estimate & PARAM_GAMMA) {
			fprintf(fp, "gamma: ");
			fprint_vectorized_sq_matrix(fp, gamma,
							NUM_NUCLEOTIDES, 1);
		}
		if (opt->modo->background_model) {
			if (opt->modo->param_estimate & PARAM_BG_PI) {
				fprintf(fp, "bg_pi: ");
				fprint_doubles(fp, bg_pi, NUM_NUCLEOTIDES, 6, 1);//
			}
			if (opt->modo->model_quality) {
				fprintf(fp, "bg_lambda: ");
				fprint_vectorized_matrix(fp, bg_lambda,
					dat->max_read_position, mod->n_quality, 1);//
			}
		}
	}

	unsigned int *cid = which == OPTIMAL_VALUES ? ri->optimal_cluster_id
							: ini->best_cluster_id;
	fprintf(fp, "assignments: ");
	fprint_uints(fp, cid, dat->sample_size, 2, 1);

	unsigned char *haps = which == OPTIMAL_VALUES ? mod->best_haplotypes
						: mod->haplotypes;
	size_t *qprob = malloc(mod->n_quality * sizeof *qprob);
	size_t *nprob = malloc(mod->n_quality * sizeof *nprob);
	for (size_t q = 0; q < mod->n_quality; ++q) {
		qprob[q] = 0;
		nprob[q] = 0;
	}
	for (size_t i = 0; i < dat->sample_size; ++i)
		for (size_t j = 0; j < dat->lengths[i]; ++j) {
			if (dat->dmat[i][j] != haps[
				cid[i] * dat->max_read_position + j
				+ (dat->offset ? dat->offset[i] : 0)])
				++qprob[(int) dat->qmat[i][j]];
			nprob[(int) dat->qmat[i][j]]++;
		}
	fprintf(fp, "     PHRED error probabilities:");	/* TODO: assumes certain probs */
	fprint_error_probs(fp, dat->fdata);
	fprintf(fp, "\nempirical error probabilities:");
	for (size_t q = MIN_ASCII_QUALITY_SCORE;
		q <= MAX_ASCII_QUALITY_SCORE; ++q) {
		if (q < dat->fdata->min_quality)
			continue;
		else if (q > dat->fdata->max_quality)
			break;
		else
			fprintf(fp, " %f", (double) qprob[q
				- dat->fdata->min_quality]
				/ nprob[q - dat->fdata->min_quality]);
	}
	fprintf(fp, "\n");

	free(qprob);
	free(nprob);

	/* report log likelihood of each read in most likely assignment */
	fprintf(fp, "read ll under quality model (estimated haplotypes):");
	for (size_t i = 0; i < dat->sample_size; ++i) {
		size_t k = cid[i];
		double ll = 0;
		for (size_t j = 0; j < dat->lengths[i]; ++j) {
			double eprob = error_prob(dat->fdata, dat->qmat[i][j]);
			if (haps[k*dat->max_read_position + j
				+ (dat->offset ? dat->offset[i] : 0)]
							== dat->dmat[i][j])
					ll += log(1 - eprob);
				else
					ll += log(eprob)/3;
		}
		fprintf(fp, " %.3f",  ll);
	}
	fprintf(fp, "\n");

	fprintf(fp, "read ll under model:");
	if (which == OPTIMAL_VALUES)
		param_update(mod, dat, opt->modo, FROM_BEST);	/* [KSD, TODO] This is dangerous because it overwrites DEFAULT slots */
	active_fp = fp;
	e_step(dat, mod, opt->modo, PRINT_ESTD_READ_LLS);
	fprintf(fp, "\n");

	unsigned int *cs = which == OPTIMAL_VALUES ? ri->optimal_cluster_size
						: ini->best_cluster_size;
	fprintf(fp, "sizes: ");
	fprint_uints(fp, cs, mod->n_mix, 3, 1);
	fprint_fasta(fp, haps, opt->modo->K, dat->max_read_position, "H");

	if (opt->modo->JC69_model) {	/* these parameters are not really estimated */
		fprintf(fp, "Estimated common ancestor: \n");
		fprint_fasta(fp, mod->est_ancestor, 1,
				 dat->max_read_length, "A");
		fprintf(fp, "Evolution_rate: ");
		fprint_doubles(fp, mod->distance, mod->K, 6, 1);
		fprintf(fp, "log likelihood from JC69 model: %f\n", mod->JC_ll);
	}
	fprintf(fp, "log likelihood: %f\n", which == OPTIMAL_VALUES ? mod->best_ll
								: mod->ll);
	fprintf(fp, "aic: %f\n", mod->aic);
	fprintf(fp, "bic: %f\n", mod->bic);

	if (opt->do_simulation) {
		fprintf(fp, "Comparing to truth:\n");
		fprintf(fp, "ARI: %f\n",
			cluster_index(cid, sim->true_cluster_id,
			dat->sample_size, mod->n_mix, mod->n_mix,
			ADJUSTED_RAND_INDEX));
		if (opt->modo->param_estimate & PARAM_PI)
			fprintf(fp, "\\Delta pi: %f\n", euc_dis(pi,
				sim->mod->pi, mod->n_mix) / mod->n_mix);
		if (opt->sim_modo->parameterization == MLOGIT_PARAMETERIZATION
			&& opt->modo->parameterization == MLOGIT_PARAMETERIZATION) {
			if (opt->modo->param_estimate & PARAM_BETA)
				fprintf(fp, "\\Delta beta: %f\n", euc_dis(beta,
					sim->mod->beta, mod->n_beta_coef)
					/ mod->n_beta_coef);
		} else if (opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION
			&& opt->modo->parameterization == DEFAULT_PARAMETERIZATION) {
			if (opt->modo->param_estimate & PARAM_DELTA)
				fprintf(fp, "\\Delta delta: %f\n", euc_dis(
					delta, sim->mod->delta,
					dat->max_read_position)
					/ dat->max_read_position);
			if (opt->modo->param_estimate & PARAM_LAMBDA) {
				fprintf(fp, "\\Delta lambda0: %f\n", euc_dis(
					lambda0, sim->mod->lambda0,
					dat->max_read_position * mod->n_quality)
					/ (dat->max_read_position
							* mod->n_quality));
				fprintf(fp, "\\Delta lambda1: %f\n",  euc_dis(
					lambda1, sim->mod->lambda1,
					dat->max_read_position * mod->n_quality)
					/ (dat->max_read_position
							* mod->n_quality));
			}
		}
		if ((opt->sim_modo->parameterization == DEFAULT_PARAMETERIZATION
			|| opt->sim_modo->parameterization == ART_PARAMETERIZATION
			|| opt->sim_modo->parameterization == QUALITY_PARAMETERIZATION)
			&& opt->modo->parameterization == DEFAULT_PARAMETERIZATION
			&& opt->modo->param_estimate & PARAM_GAMMA)
			fprintf(fp, "\\Delta gamma: %f\n", euc_dis(
				gamma, sim->mod->gamma, NUM_NUCLEOTIDES
				* NUM_NUCLEOTIDES) / (NUM_NUCLEOTIDES
					* (NUM_NUCLEOTIDES - 1)));
		if (opt->modo->param_estimate & PARAM_HAPLOTYPE) {
			fprintf(fp, "haplotype distances (%u x %u):\n", opt->modo->K,
							opt->simo->true_K);
			for (size_t j = 0; j < opt->modo->K; ++j) {
				for (size_t k = 0; k < opt->simo->true_K; ++k)
					fprintf(fp, " %2lu", hamming_char_dis(
					(char *) &haps[j*dat->max_read_length],
					(char *) &sim->mod->haplotypes[k*dat->max_read_length],//
					dat->max_read_length));
				fprintf(fp, "\n");
			}
		}
	}

	fclose(fp);

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "file '%s' closed\n", outfile);

WRITE_RESULTS_RETURN:
	return NO_ERROR;
} /* write_results */


/**
 * Create run_info object.
 *
 * @param ri	pointer to run_info object
 * @param dat	pointer to data object
 * @param sim	pointer to simulator object
 * @param opt	pointer to options object
 */
int make_run_info(run_info **ri, data *dat, simulator *sim, options *opt)
{
	run_info *rio;
	*ri = malloc(sizeof **ri);

	if (*ri == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "run_info");

	rio = *ri;
	rio->opt = opt;
	rio->sim = sim;

	/* initialize status variables*/
	rio->optimal_cluster_id = calloc(dat->sample_size,
		sizeof *rio->optimal_cluster_id);

	if (!rio->optimal_cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_cluster_id");

	rio->optimal_seed_idx = calloc(opt->modo->K, sizeof *rio->optimal_seed_idx);		// opt->modo

	if (!rio->optimal_seed_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_seed_idx");

	rio->optimal_cluster_size = calloc(opt->modo->K + opt->modo->background_model,		// opt->modo
		sizeof *rio->optimal_cluster_size);

	if (!rio->optimal_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_cluster_size");

	return NO_ERROR;
}/* make_run_info*/

/**
 * Reallocate run_info objfect for a different K or sample_size.
 *
 * [KSD, TODO] Should you not reallocate if there is no change?  Is there cost
 * to calling realloc for the same size?
 *
 * @param ri		pointer to run_info object
 * @param sample_size	size of new sample
 * @param K		new K
 * @return		error status
 */
int realloc_run_info(run_info *ri, size_t sample_size, unsigned int K)
{
	size_t *optimal_seed_idx = realloc(ri->optimal_seed_idx,
		K * sizeof *ri->optimal_seed_idx);

	if (!ri->optimal_seed_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.run_info.optimal_seed_idx");

	ri->optimal_seed_idx = optimal_seed_idx;

	unsigned int *optimal_cluster_size = realloc(ri->optimal_cluster_size,
			K * sizeof *ri->optimal_cluster_size);

	if (!ri->optimal_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.run_info.optimal_cluster_size");

	ri->optimal_cluster_size = optimal_cluster_size;

	unsigned int *optimal_cluster_id = realloc(ri->optimal_cluster_id,
		sample_size * sizeof *ri->optimal_cluster_id);

	if (!ri->optimal_cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.run_info.optimal_cluster_id");

	ri->optimal_cluster_id = optimal_cluster_id;


	return NO_ERROR;
}/* realloc_run_info */


/**
 * Delete run_info object
 */
void free_run_info(run_info *ri)
{
	if (ri) {
		if (ri->optimal_cluster_id)
			free(ri->optimal_cluster_id);
		if (ri->optimal_seed_idx)
			free(ri->optimal_seed_idx);
		if (ri->optimal_cluster_size)
			free(ri->optimal_cluster_size);
		free(ri);
	}
}/*free_run_info*/
