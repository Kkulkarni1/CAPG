/**
 * @file run_khaplotype.c
 * @author Yudi Zhang
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>

#include "array.h"
#include "cluster_lloyds.h"
#include "myfun.h"
#include "fastq.h"
#include "haplotype_data.h"
#include "haplotype_option.h"
#include "effici.h"
#include "cmdline.h"
//#include "test_func.h"
//#include "test.h"
#include "initialize_kmodes.h"

static inline void stash_state(data *dat, options *opt);

int main(int argc, const char **argv)
{
	int err = NO_ERROR;
	data *dat = NULL;
	options *opt = NULL;
	fastq_options *fqo = NULL;	/* fastq file options */
	
	/* parse command line */
	if ((err = make_options(&opt)))
		goto CLEAR_AND_EXIT;
	
	if ((err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;
	
	/* make data object */
	if ((err = make_data(&dat)))
		goto CLEAR_AND_EXIT;
	
//	if (opt->run_with_quals) {
		
		if ((err = make_fastq_options(&fqo)))
			goto CLEAR_AND_EXIT;
		
		/* encode bases as xy_t: error raised if ambiguous bases. */
		fqo->read_encoding = XY_ENCODING;
		
		/* read sequence data */
		if (opt->datafile && (err = read_fastq(opt->datafile,
						&dat->fdata, fqo)))
			goto CLEAR_AND_EXIT;
		
		/* with data now loaded, can polish off data object */
		if ((err = sync_state_data(dat, opt)))
			goto CLEAR_AND_EXIT;

//	} else {
		//	goto CLEAR_AND_EXIT;
		//	[KSD] I don't think this works yet, as it is buried in run_kmodes.c.
		//	To test if I have the same results with kmodes_lloyds
//		if ((err = read_data(dat, opt)))
//			goto CLEAR_AND_EXIT;
//		
//		if ((err = finish_make_data(dat, opt)))
//			goto CLEAR_AND_EXIT;
//	}

	if ((err = make_seeds(dat, opt)))
		goto CLEAR_AND_EXIT;

	if ((err = kmodes_ini(dat, opt)))	// [KSD] Where is this function?
		goto CLEAR_AND_EXIT;
	
	FILE *fps = NULL;
	
	//FILE *fp1 = fopen("dmat.txt", "w"), *fp2 = fopen("qmat.txt", "w");
	//
	//for (int i = 0; i < dat->n_observations; ++i) {
	//	for (int j = 0; j < dat->n_coordinates; ++j) {
	//		fprintf(fp1, "%d\t", dat->dmat[i][j]);
	//		fprintf(fp2, "%d\t", dat->qmat[i][j]);
	//	}
	//	fprintf(fp1, "\n");
	//	fprintf(fp2, "\n");
	//}
	//
	//fclose(fp1); fclose(fp2);
	//exit(123);
	
	/* Record time */
	clock_t begin, end;
	double cost;
	begin = clock();
	
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
		
		/* Run of the algorithm */
		if (opt->kmodes_algorithm == FASTQ_LLOYDS)
			dat->total = cluster_lloyds(dat, dat->use_ini ?
				dat->ini_seeds : dat->seeds, dat->cluster_size,
				dat->cluster_id, dat->n_observations,
				dat->n_coordinates, opt->K, opt->n_max_iter,
				dat->criterion, &err, &dat->iter,
				fastq_lloyds_step1, fastq_lloyds_step2,
						compute_criterion_hap);
//		else if (opt->kmodes_algorithm == LLOYDS_OLD) {
//			dat->total = lloyd_fastq(dat,
//						 dat->use_ini ? dat->ini_seeds : dat->seeds, dat->cluster_size,
//						 dat->cluster_id, dat->n_observations, dat->n_coordinates,
//						 opt->K, opt->n_max_iter, dat->criterion, &err, &dat->iter, compute_criterion_hap);
//		}
		
		else if (opt->kmodes_algorithm == FASTQ_LLOYDS_EFFICIENT)
		/* [KSD] This function needs to be written to use the same interface
		 * as cluster_lloyds(), while taking different step1 and step2
		 * functions.  All the auxiliary variables needed to make it more
		 * efficient should be available in dat and allocated as part of
		 * the data setup, say in sync_state_hap().
		 */
			dat->total = cluster_lloyds2 (dat, dat->use_ini ?
				dat->ini_seeds : dat->seeds, dat->cluster_size,
				dat->cluster_id, dat->n_observations,
				dat->n_coordinates, opt->K, opt->n_max_iter,
				dat->criterion, &err, &dat->iter,
				fastq_lloyds_efficient_step1,
				fastq_lloyds_efficient_step2, compute_cost);
		
		else if (opt->kmodes_algorithm == FASTQ_MACQUEEN)
			dat->total = cluster_macqueen(dat, dat->n_observations,
				dat->n_coordinates, dat->use_ini ? dat->ini_seeds
				: dat->seeds, opt->K, dat->cluster_id,
				dat->cluster_size, opt->n_max_iter,
				dat->criterion, &err, &dat->iter,
				fastq_macqueen_ini, fastq_macqueen_iter,
								compute_cost);
		else if (opt->kmodes_algorithm == FASTQ_HW)
			dat->total = cluster_hw(dat, dat->use_ini ? dat->ini_seeds : dat->seeds, dat->cluster_size, dat->cluster_id, dat->n_observations, dat->n_coordinates, opt->K, opt->n_max_iter, dat->criterion, &err, &dat->iter, hw_init);
		
		else {
			mmessage(ERROR_MSG, INVALID_USER_INPUT,
				 "Invalid algorithm.\n");
			goto CLEAR_AND_EXIT;
		}
		
		/* repeat if obtain null cluster */
		if (err == KMODES_NULL_CLUSTER_ERROR) {
			if (i)
				--i;
			continue;
		} else if (err) {
			mmessage(ERROR_MSG, CUSTOM_ERROR, "%s\n",
				 kmodes_error(err));
			if (err != KMODES_EXCEED_ITER_WARNING)
				goto CLEAR_AND_EXIT;
		}
		
		/* record best solution */
		if ((!err || err == KMODES_EXCEED_ITER_WARNING)
				&& dat->total > dat->best_total) {
			stash_state(dat, opt);
		} else if (err)
			fprintf(stderr, "[ERROR] %s (%d)\n",
				kmodes_error(err), err);
	}
	
	end = clock();
	cost = (double)(end - begin)/CLOCKS_PER_SEC;
	printf("Time cost is: %lf secs\n", cost);
	
	/* Write solution */
	write_solution(dat, opt, &fps);
	
CLEAR_AND_EXIT:
	
	if (dat)
		free_data(dat);
	if (opt)
		free_options(opt);
	
	return(err);
} /* main */

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
	
	COPY_1ARRAY(dat->best_cluster_id, dat->cluster_id, dat->n_observations);
	COPY_1ARRAY(dat->best_cluster_size, dat->cluster_size, opt->K);
	COPY_1ARRAY(dat->best_criterion, dat->criterion, opt->K);
	
	if (opt->shuffle)
		memcpy(dat->best_obsn_idx, dat->obsn_idx, dat->n_observations
						* sizeof *dat->obsn_idx);
} /* stash_state */



/**
 * print command line usage
 *
 * @param fp file
 * @param cmdname command line
 * @param obj option struct
 */
void fprint_usage(FILE *fp, const char *cmdname, void *obj)
{
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
	fprintf(fp, "\t-e\n\t\tRun Lloyd's efficient algotithm for haplotype.\n");
	fprintf(fp, "\t-h97\n\t\tRun Macqueen(Huang)'s algorithm (combine with -u to replicate klaR).\n");
	fprintf(fp, "\t-w\n\t\tRun Hartigan and Wong algorithm (can combine with -u).\n");
	fprintf(fp, "\t-1\n\t\tSubtract 1 from the observation categories.\n");
	fprintf(fp, "\t-f <ffile> [<ffile2>]\n\t\tSet the input filename (with -s, simulate data are written to this file).\n");
	fprintf(fp, "\t\tIf <ffile2> specified, then write the data to this file: same as <ffile> to overwrite.\n");
	fprintf(fp, "\t-o <ofile> [<ofile2>]\n\t\tSet the output filename.  If second argument given, split information into first.\n");
	fprintf(fp, "\t-i <istr>\n\t\tSet initialization method (rnd|h97|h97rnd|hd17|clb09|clb09rnd|av07|av07grd|rndp|rnds).\n");
	fprintf(fp, "\t\t<suint1> ... <suintk>\n\t\tSet the indices of the seeds.\n");
	fprintf(fp, "\t\t<ifile>\n\t\tProvide file with possible seeds.\n");
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
