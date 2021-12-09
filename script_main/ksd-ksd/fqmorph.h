/**
 * @file fqmorph.h
 * @author Karin S. Dorman
 *
 * Header file for fqmorph structs, extern functions and defines.
 */

#ifndef __H_FQMORPH__
#define __H_FQMORPH__

#include <stdio.h>
#include <stddef.h>

#include "fastq.h"


typedef struct _options options;
typedef struct _data data;

enum {
	FQMORPH_NO_ERROR,
	FQMORPH_AMBIGUOUS_NUCLEOTIDE	/* example only */
};

/**
 * Output format types.
 */
enum {
	FASTQ_FORMAT,
	FASTA_FORMAT,
	TABLE_FORMAT,	/* format used by k-modes */
	NUMBER_FORMATS
};

/**
 * Run options.
 */
struct _options {
	/* data */
	char const *fastq_file;	/*<! name of fastq input file */
	char const *partition_file;	/*<! partition file */
	char const *partition_outfile;	/*<! partition output file */
	fastq_options *fqo;

	/* running */
	int reverse_complement;	/*<! request reverse complement */
	int split_file;		/*<! use partition file to split */
	int output_clusters;	/*<! output cluster assignments */
	int append;		/*<! append output to existing file */
	int output_format;	/*<! format of output */
	int expected_errors;	/*<! expected number of errors */
	int read_names;		/*<! retain fastq names */
	int qs_to_fasta;	/*<! output quality scores to fasta */
	unsigned long seed;	/*<! random number seed [srand()] */
	unsigned int cut_start;	/*<! start of read region to cut */
	unsigned int cut_end;	/*<! end of read region to cut */
	unsigned int nsample;	/*<! randomly sample this many reads */
	unsigned int nni_k;	/*<! k of k-nearest neighbor */
	double max_ee;		/*<! maximum allowed expected errors */
				/*<! max. avg. per site err. prob. */
	double leading_ee;	/*<! in 5' window */
	double trailing_ee;	/*<! in 3' window */
	char_t leading_qs;	/*<! as in trimmomatic LEADING */
	char_t trailing_qs;	/*<! as in trimmomatic TRAILING */
	unsigned int leading_nuc;	/*<! number of 5' nucs to trim */
	unsigned int trailing_nuc;	/*<! number of 3' nucs to trim */
	unsigned int max_length;	/*<! trim 3' to maximum read length */
	int s_class;		/*<! selected class */
	unsigned int n_name_match;	/*<! number name match patterns */
	char const **name_match;	/*<! pattern to match in name */
	unsigned char report_lengths;	/*<! report lengths of reads */
	unsigned char report_names;	/*<! report names of reads */

	unsigned char strip_slash;	/*<! strip /# from names */
	unsigned char strip_dash;	/*<! strip -.* from names */
	unsigned char casavize;	/*<! convert names to casava format */

	char const *outfile;	/*<! output file */
	char const *distn_file;	/*<! distribution file */
}; /* options */

/**
 * Data.
 */
struct _data {
	unsigned int *cluster_id;	/*<! store contents of partition file */
	/* could add pointers to reads here */
}; /* data */

const char * fqmorph_error_message(int err_no);


#endif
