#ifndef __ROSHAN_H__
#define __ROSHAN_H__

#include <stdio.h>
#include <ulimit.h>

#include "vcf.h"

#define ROSHAN_VERSION 1.0
#define ALLOW_DEBUGGING
#define N_SUBGENOMES 2
#define NA_ASSIGNMENT	UINT_MAX

typedef struct _options options;
typedef struct _input input;

/**
 * structure for information about amplici-related coverage tests
 * Let p1, p2, p3, p4, be the ordered true proportions of the four most abundant
 * haplotypes.
 * H04: p1 = p2, p3 = p4 (4 haplotypes, equal chromosomal coverage)
 * H04a: p1 = p2 = p3 = p4 (4 haplotypes, equal chromosomal/subgenomic coverage)
 * H03a: p1, p2 = p3 (3 haplotypes, equal chromosomal coverage)
 * H03b: p1 = p2, p3 (3 haplotypes, equal chromosomal coverage)
 * H02: p1, p2 (2 haplotypes)
 */
struct _input {

	unsigned int *assignment;	/*<! (0.. -> cluster_id) assigned cluster of each read */
	unsigned int *cluster_sizes;	/*<! (cluster_id -> int) cluster sizes */
	unsigned int K;			/*<! no. haplotypes (clusters) found by AmpliCI */
	unsigned int proptest;		/*<! levels of screening based on prop test: 0 means not screening, 1 means rejecting the null that has 4 haplotypes, 2 rejects the nulls that has 3 haplotyoes */
	unsigned int n_prop_exclude;

	/**
	 * coordinate systems:
	 * integers: 0, 1, 2, ...
	 * integers[n]: 0, 1, 2, ..., n-1
	 * cluster_id: 0, 1, 2, 3, 4, ..., K-1
	 * select_id (retained cluster_ids): 0, 2, 3
	 * include_id == integers[keep]: 0, 1, 2 for keep 3
	 * csome_id (A1, A2, B1, B2) == integers[4]: 0, 1, 2, 3 for allotetraploid
	 */
	size_t *sort_id;				/*<! (integers[K] -> cluster_id) abundance sorted clusters */
	unsigned int *excluded_clusters;		/*<! (cluster_id -> 0|1) indicate excluded clusters (of K) */
	unsigned int included_clusters[2*N_SUBGENOMES];	/*<! (include_id -> cluster_id) included clusters */
	double coverage_A[2*N_SUBGENOMES];		/*<! (include_id -> float) expected A coverage */
	double proportion_A[2*N_SUBGENOMES];		/*<! (include_id -> float) proportion A coverage */
//	size_t haplotype_id[2*N_SUBGENOMES];		/*<! (csome_id -> cluster_id) H0, H1, ... of each csome */
	unsigned int *cluster_to_csome;			/*<! (cluster_id -> csome_id) haplotype id H0, H1, ... to csome */
	unsigned int csome_to_inc[2*N_SUBGENOMES];	/*<! (csome_id -> include_id) include id of each csome, may repeat */
	unsigned int n_indels[N_SUBGENOMES];		/*<! number indels within subgenome alignment */
	unsigned char **ref_profile_alignments[N_SUBGENOMES];
	unsigned int reject_state;		/*<! see enum abundance_test */
	double mle_w;				/*<! estimated w */
	double pval_equal_subgenomic_coverage;	/*<! p-value for H04a */
	size_t n_excluded;			/*<! no. of excluded haplotypes */
	size_t n_hash_excluded;			/*<! no. already excluded reads */
	size_t n_observation;			/*<! no. of reads in fastq file */
};

enum abundance_test {
	ACCEPT_FOUR,	/* p1 = p2, p3 = p4 */
	ACCEPT_THREE_A,	/* p1, p2 = p3 */
	ACCEPT_THREE_B,	/* p1 = p2, p3 */
	ACCEPT_TWO	/* p1, p2 */
};

enum msa_tool {
	CLUSTAL_OMEGA,
	MAFFT
};


struct _options {
	unsigned char amplicon;		/*<! amplicon data? */
	unsigned char drop_unmapped;	/*<! drop reads that are unmapped
					 * in either alignment
					 */
	unsigned char drop_secondary;	/*<! drop secondary alignments */
	unsigned char display_alignment;/*<! display alignments to stderr */
	int drop_soft_clipped;		/*<! drop reads soft-clipped by this
					 * length of more in either alignment
					 */
	int drop_indel;			/*<! drop reads whose alignments
					 * contain indel this long or more
					 * in either alignment
					 */
	int proptest_screen;		/*<! perform subgenomic coverage tests */
	int amplici_complete;		/*<! always run AmpliCI to completion */
	double proptest_alpha;		/*<! critical value for rejection */
	unsigned int proptest_min_freq;	/*<! minimum read frequency */
	double coverage_screen;		/*<! drop sites w/ low coverage */
	double biallelic_screen;	/*<! drop sites w/ >2 alleles */
	double min_expected_coverage;	/*<! abort if subgenome coverage < this */
	double				/*<! min. posterior alignment prob. */
		min_posterior_alignment_prob;
	double max_misalignment_proportion;
					/*<! maximum allowed proportion of mis-
					 *   alignments, for use in checking
					 *   H03b.
					 */
	int posthoc_coverage_test;	/*<! perform post hoc coverage test */
	int equal_homolog_coverage_test;	/*<! perform equal homologous coverage test */
	double phc_min_genotype_pp;	/*<! min. posterior probability to call heterozygote */
	double phc_min_alignment_pp;	/*<! min. posterior probability of alignment */
	int min_length;			/*<! drop reads shorter than this */
	int max_length;			/*<! drop reads longer than this */
	double max_eerr;		/*<! drop reads with more expected
					 * errors than this
					 */
	double min_log_likelihood;	/*<! drop reads smaller log likelihood */
	int max_quality_score;		/*<! censory quality scores here */
	unsigned int n_sample;		/*<! number of Monte Carlo samples */
	const char *param_file;		/*<! error probabilities file */
	FILE *error_file;		/*<! error data file */
	unsigned int use_bam;		/*<! files are in bam format */
	const char *sbam_files[N_SUBGENOMES];/*<! sam/bam files */
	const char *fsa_files[N_SUBGENOMES];	/*<! fsa files (DEPRECATED) */
	const char *ref_names[N_SUBGENOMES];	/*<! name of references */
	const char *ref_fsas[N_SUBGENOMES];	/*<! reference fsa files */
	char *subref_fsas[N_SUBGENOMES];	/*<! selected reference fsa files */
	const char *vcf_files[N_SUBGENOMES];	/*<! vcf files */
	const char *ref_alignment;	/*<! alignment of references (SAM) */
	const char *subref_fsa_base;	/*<! base name for subref_fsas */
	vcf_options *vcf_opt;		/*<! output genotype likelihoods */
	const char *sample_name;	/*<! name of accession */
	int write_fastq_and_quit;	/*<! unexcluded reads to fastq */
	const char *fastq_after_filter;	/*<! after paralog filtering */
	const char *sam_after_filter
			[N_SUBGENOMES];	/*<! after paralog filtering */
	int genotype_by_clustering;	/*<! genotype by clustering */
	const char *amplici_command;	/*<! screen reads w/ ampliCI command */
	const char *ac_fastq_file;	/*<! ampliCI fastq filename */
	const char *ac_outfile;		/*<! ampliCI output filename */
	double ac_low_bound;		/*<! ampliCI lower bound */
	double ac_ll_threshold;		/*<! ampliCI log likelihood lower bound */
	int msa;
	const char *clustalo_command;	/*<! clustal omega command */
	const char *clustalo_infile;	/*<! infile for clustalo alignment */
	const char *clustalo_outfile;	/*<! outfile for clustalo alignment */
	const char *mafft_command;
	const char *mafft_infile;
	const char *mafft_outfile;
	const char *coverage_file;
	int debugging_site;		/*<! debug this site (0 no) */
	int karin_old;
};


#endif
