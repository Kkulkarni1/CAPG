#ifndef __CAPG_H__
#define __CAPG_H__

#include <stdio.h>
#include "uthash.h"
#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "vcf.h"

#define N_FILES 2
typedef struct _options options;
typedef struct _input input;

struct _input {
	unsigned int *exclude;		/*<! index of excluded cluster after cut-off */
	unsigned int *exclude_id;	/*<! index of excluded reads after cut-off */
	unsigned int *not_input;	/*<! index of already excluded reads in
					 *   merge hash, but not in fastq file */
	unsigned int *assignment;	/*<! assigned cluster of each read */
	unsigned int K;			/*<! no. haplotypes (clusters) */
	unsigned char proptest;		/*<! levels of screeing based on prop test: 0 means not screeing, 1 means rejecting the null that has 4 haplotypes, 2 rejects the null that has 3 haplotyoes */
	unsigned int n_prop_exclude;
	size_t n_excluded;		/*<! no. of excluded haplotypes */
	size_t n_excluded_id;		/*<! no. of excluded reads by cut-off */
	size_t n_hash_excluded;		/*<! no. already excluded reads */
	size_t n_observation;		/*<! no. of reads in fastq file */
};


struct _options {
	unsigned int drop_soft_clipped;	/*<! drop reads soft-clipped by this
					 * length of more in either alignment
					 */
	unsigned int drop_indel;	/*<! drop reads whose alignments
					 * contain indel this long or more
					 * in either alignment
					 */
	unsigned int min_length;	/*<! drop reads shorter than this */
	unsigned int max_length;	/*<! drop reads longer than this */
	double weight_penalty;		/*<! drop reads whose coverage */
	double biallelic_screen;	/*<! drop sites w/ >2 alleles */
	double min_expected_coverage;	/*<! abort if subgenome coverage < this */
	double min_genotype_post_prob;	/*<! min. posterior probability to call heterozygote */
	double min_alignment_post_prob;	/*<! min. posterior probability of alignment */
	double max_eerr;		/*<! drop reads with more expected
					 * errors than this
					 */
	double min_log_likelihood;	/*<! drop reads smaller log likelihood */
	int max_quality_score;		/*<! censory quality scores here */
	unsigned int n_sample;		/*<! number of Monte Carlo samples */
	const char *param_file;		/*<! error probabilities file */
	FILE *error_file;		/*<! error data file */
	unsigned int use_bam;		/*<! files are in bam format */
	const char *sbam_files[N_FILES];/*<! sam/bam files */
	const char *fsa_files[N_FILES];	/*<! fsa files */
	const char *ref_names[N_FILES];	/*<! name of references */
	const char *vcf_files[N_FILES];	/*<! vcf files */
	const char *output_file;	/*<! final output file */
	vcf_options *vcf_opt;		/*<! output genotype likelihoods */
	const char *extracted_rf; /*<! targeted reference fsa file name */
	int write_fastq_and_quit;	/*<! unexcluded reads to fastq */
	const char *ampliclust_command;	/*<! screen reads w/ ampliCI command */
	const char *ac_fastq_file;	/*<! ampliCI fastq filename */
	const char *ac_outfile;		/*<! ampliCI output filename */
	double ac_low_bound;		/*<! ampliCI lower bound */
	const char *coverage_file;
	char const *sam_file;		/*<! reference sam file */
	const char *sample_name;	/*<! name of accession */

	unsigned char display_alignment;/*<! display alignments to stderr */
	unsigned char drop_unmapped;	/*<! drop reads that are unmapped
					 * in either alignment
					 */
	unsigned char proptest_screen;	/*<! {0,1,2,3} drop reads according to abundance test */
	unsigned char drop_secondary;	/*<! drop secondary alignments */
	unsigned char equal_homolog_coverage_test; /*<! perform equal coverage test */
	unsigned char posthoc_coverage_test;	/*<! perform post hoc coverage test */
	/* this is how we ran CAPG: references are chrom:start-end,
	 * where start, end are 0-based [inclusive, exclusive), but 
	 * the default is now samtools format, where start, end are
	 * 1-based, [inclusive, inclusive].
	 */
	unsigned char legacy_region_specification;
};

typedef struct nuc_state {
	UT_hash_handle hh;		/* 56 bytes w/ pointer alignment */
	double prob;			/* 8 bytes */
	unsigned int pos;		/* 4 */
	data_t true_nuc, read_nuc;	/* 1, 1 bytes */
	data_t quality;			/* 1 */
	/* total: 71 bytes round up to
	 * 72 with 1 byte slop
	 */
} nuc_state;

typedef struct mlogit_stuff {
	nuc_state *ns;
	int pos;
} mlogit_stuff;

int read_ampliclust_results(FILE *fp, input *in);
int make_input(input **in);
void free_input(input *in);
int default_options(options *opt);
int parse_options_capg(options *opt, int argc, const char **argv);
nuc_state *read_param_file(char const *param_file);
double mlogit_sub_prob(data_t s, data_t r, data_t q, void *vptr);
double mlogit_sub_lprob(data_t s, data_t r, data_t q, void *vptr);
#endif
