//
//  simulation.h
//  Test
//
//  Created by Yudi Zhang on 5/14/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#ifndef simulation_h
#define simulation_h

#include <stdio.h>
#include "fastq.h"
#include "nuc.h"
#include "qual.h"

typedef struct _simu_options simu_options;
typedef struct _simu_data simu_dat;

struct _simu_options {
	unsigned int len_N;		/*<! length of sequence sampled from fsa */
	const char *out_file;		/*<! out_file (individual fasta) */
	const char *fsa_file;		/*<! fsa files (A & B seperately && together) */
	const char *ref_name;		/*<! random sample ref name */
	const char *out_sam;		/*<! alignment of A and B in a sam file */
	const char *extracted_rf;	/*<! subgenome A reference */
	const char *sub_ref_b;		/*<! already existing subgenome B reference */
	char *delim_ref;
	char *delim_len;
	unsigned long seed;
	double homo_rate;		/*<! rate of homologous SNPs rate */
	double heter_rate;		/*<! rate of homologous SNPs rate */
	double alpha;			/*<! beta distribution(proportion of the reference allele at the jth allelic SNP) */
	double beta;			/*<! beta distribution(proportion of the reference allele at the jth allelic SNP) */
	double substitution_rate;	/*<! snps substition rate */
	double prop_allele;		/*<! HWE */
	unsigned int num_ind;
	
	const char *ART_command;	/*<! art command */
	const char *error_file1;	/*<! error files */
	const char *error_file2;
	const char *fq_file;		/*<! simulated fastq file */
	int length;
	int coverage;
	unsigned int imperfect_ref;	/*<! whether to simulate imperfect ref A */
	double mismatch_prob;		/*<! probability of mismatch in subgenome A */
	
	const char *bwa_command;	/*<! bwa to make alignment */
	const char *samAB;		/*<! prefix of read-ref sam files */
};

struct _simu_data {
	char_t *seq_A;			/*<! perfect reference A */
	char_t *seq_B;			/*<! perfect reference B */
	char_t *seq_A2;			/*<! 2nd parental csome A */
	char_t *seq_B2;			/*<! 2nd parental csome B */
	char_t *ref_A;			/*<! actual reference A */
	char_t **ind;
	unsigned int *homo_loci;
	unsigned int *heter_loci;
	unsigned int *mm_loci;		/*<! mismatches between ref and true subgenome */
};

void make_simu_opt(simu_options *opt);
int parse_opt(simu_options *opt, int argc, const char **argv);
int make_data(simu_dat **dat);
int load_data(simu_dat *dat, simu_options *opt, fastq_data *fds);
void fprint_fsa(FILE *fp, char_t **data, size_t n, size_t p, char const * const prefix);
void fprint_seq(FILE *fp, char_t *data, size_t p, char const * const prefix);
void fprint_usage(FILE *fp, const char *cmdname, void *obj);
void write_sam(FILE *fp, char_t *B, size_t p, char *name_A, char *name_B);
void call_art(simu_options *opt, char *fq_out, char *fq_in);
void call_bwa(simu_options *opt, char *ref_in, char *reads1, char *reads2, char *sam);
#endif /* simulation_h */
