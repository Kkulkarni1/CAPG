/**
 * @file haplotype_data.h
 * @author Yudi Zhang
 */

#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H

#include "fastq.h"
#include "constants.h"
#include "haplotype_option.h"
#include "io_kmodes.h"
#include "io.h"
#include "myfun.h"

#define RAND_SEED 0

typedef struct _data data;


/**
 * Data structure for clustering fastq data.
 */
struct _data {
	
	/* fastq data */
	fastq_data *fdata;
	data_t **qmat;		/*<! quality score matrix*/
	unsigned char n_quality;/*<! number quality scores [min, max] */
	double **prob_t;	/*<! probability of no error: log(1-p_ij) */
	double **prob_f;	/*<! probability of error: log(p_ij/3) */
	
	unsigned int max_read_length;	/*<! maximum observed read length */
	unsigned int min_read_length;	/*<! minimum observed read length */
	unsigned int *lengths;		/*<! length of reads */
	
	/* index */
	size_t *read_idx;		/*<! index array of all reads* */
	
	/* truth when data simulated */
	unsigned int *true_cluster_id;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */
	
	data_t *data;			/*<! data */
	data_t **dmat;			/*<! data as matrix */
	unsigned int n_observations;	/*<! number of observations */
	unsigned int n_coordinates;	/*<! no. of coordinates (all categorical) */
	data_t *n_categories;		/*<! no. of categories per coordinate */
	data_t max_n_categories;	/*<! max. number categories per coordinate */
	data_t tot_n_categories;	/*<! total number of categories */
	data_t *categories;		/*<! total no. of the category */
	
	/* auxiliary variables */
	double ***e_kjN;	/* likelihood for each category of each site per cluster */
	double **v_ik;		/* likelihood of each observation to every cluster */
	int **last_cent;	/* last centers */
	int *last_assign;	/* last assignment of each observation */
	double *last_cost;
	double **last_vik;
	double ***last_ekjn;
	int *index;		/* index of the observations in one perticular cluster */
	unsigned int *ic2;
	unsigned int *live;
	
	
	/* initialization */
	data_t **seeds;		/*<! chosen seeds */
	data_t **ini_seeds;	/*<! trial seeds */
	unsigned int *seed_idx;	/*<! seed indices */
	unsigned int *ini_seed_idx;	/*<! trial seed indices */
	unsigned int n_init;	/*<! number of initializations done */
	uint8_t use_ini;	/*<! using ini_* versions or not */
	
	/* current solution */
	double total;			/*<! current criterion */
	unsigned int *cluster_id;	/*<! cluster assignments */
	unsigned int *obsn_idx;		/*<! used if shuffling */
	double *criterion;		/*<! criterion */
	unsigned int *cluster_size;	/*<! cluster sizes */
	unsigned int iter;		/*<! iterations */
	
	/* best solution */
	double best_total;		/*<! total criterion */
	double best_rand;		/*<! if simulated */
	unsigned int *best_seed_idx;	/*<! seeding */
	data_t **best_modes;		/*<! estimated modes */
	double *best_criterion;		/*<! criterion */
	unsigned int *best_cluster_id;	/*<! cluster assignments */
	unsigned int *best_cluster_size;/*<! cluster sizes */
	unsigned int *best_obsn_idx;	/*<! used if shuffling */
	
	/* summary statistics */
	double seconds;			/*<! seconds used */
	double avg_cost, sd_cost;	/*<! minimum criterion */
	double avg_iter, sd_iter;	/*<! iterations to convergence */
	double avg_ar, sd_ar;		/*<! adjusted rand (truth known) */
	double avg_mi, sd_mi;		/*<! mutual information (truth known) */
	double avg_vi, sd_vi;		/*<! variance of info (truth known) */
	double avg_time, sd_time;	/*<! inits to target */
	unsigned int ntimes;		/*<! times hit target */
	
	/* internal use */
	double uncounted_seconds;
	double first_cost;
	double worst_cost;
	unsigned int max_iter;
	unsigned int ctime;
}; /* data */

int make_data(data **data);
int sync_state_data(data *dat, options *opt);
void free_data(data *dat);
int write_solution(data *dat, options *opt, FILE **in_fps);
void sample_better(unsigned int N, unsigned int n, unsigned int *idx);
void read_fsa(const char *filename, int n, int p, data_t **true_seed);
int make_seeds(data *dat, options *opt);
int shuffle_data(data *dat, options *opt);

#endif /* HAPLOTYPE_H */
