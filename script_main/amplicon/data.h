/**
 * @file data.h
 * @author Karin S. Dorman
 *
 * Header file for data struct.
 */

#ifndef __H_DATA__
#define __H_DATA__

#include "fastq.h"
#include "options.h"
//#include "constants.h"
#include "hash.h"

typedef struct _data data;

/**
 * Data.  Store the fastq file.
 */
struct _data {

	/* data */
	fastq_data *fdata;		/*<! contents of fastq file */
	size_t *coverage;		/*<! coverage at each read position */
	unsigned char n_quality;	/*<! number quality scores [min, max] */

	/* [KSD] We assume that all reads are trimmed of technical sequence and
	 * begin at the same location in the target genome.  Nevertheless, each
	 * read may have been trimmed of variable amounts of technical sequence
	 * at the 5' end.
	 *
	 * The offset indicates the differences in the amount of 5' trimming,
	 * so a positive offset indicates the target sequence occurs at higher
	 * read positions: \par data.offset.
	 *
	 * Reads may be of variable length, because of the variable trimming
	 * or some other cause.  Under these assumptions, the amplicon length is
	 * the maximal length of observed amplicon sequence is the maximum
	 * trimmed read length: \par data.max_read_length.
	 *
	 * Finally, the maximum read position is the maximum of the read
	 * read lengths plus offset, since offset read positions are consumed
	 * in the 5' end technical sequence: \par data.max_read_position.
	 */
	unsigned int *offset;		/*<! differential 5' read trims */
	unsigned int max_offset;	/*<! maximum offset */
	unsigned int max_read_position;	/*<! maximum read length + offset */

	/* [KSD] The data duplication below is a tiny flaw in design.  The
	 * fastq_data struct has no use for max and min read lengths.  Some data
	 * objects, like this one, do.  Technically, the information should be
	 * stored here.  However, the most efficient way to obtain the max/min
	 * read lengths is during fastq file reading, which suggests passing
	 * in these variables to the reader function.  On the other hand, not
	 * all users of fastq_data need the max/min information.  To keep a
	 * simple interface, we would need wrapper functions, which will
	 * make the fastq code less transparent.  Getting this right is much
	 * easier in an object-oriented language.  Here, we decide to keep
	 * simplicity at the cost of two extra variables. */

	/* the following are duplicated from fastq_data struct */
	unsigned int max_read_length;	/*<! maximum observed read length */
	unsigned int min_read_length;	/*<! minimum observed read length */

	/* [KSD] There is more aggregious replication below.  The jagged array
	 * representation below is less efficient than in fastq_data, but far
	 * easier to use.  We make the fundamental assumption that clustering
	 * will never be done directly on large datasets that cannot be
	 * retained in memory. Instead, direct clustering will be done on
	 * subsets of data. Thus, this data layer allows users of fastq data
	 * to avoid knowing the gory innards of fastq_data, at the cost of
	 * extra, but not prohibitive, memory usage. */

	/* easy access to reads and quality scores */
	char_t **dmat;		/*<! nucleotide sequences as matrix */
	char_t **qmat;	 	/*<! quality sequences as matrix */
	unsigned int *lengths;	/*<! length of reads */
	size_t sample_size;	/*<! number of reads */

	unsigned int *offset_ori;    /*<! A copy of offset, unchanged among stages */
	unsigned int max_offset_ori; /*<! A copy of max_offset, unchanged among stages */

	/* index */
	size_t *read_idx;	/*<! index array of all reads* */

	/* hash table */
	hash *seq_count; /*<! frequency of unique sequences (hash table) */
	unsigned int hash_length;  /*<! num of unique sequences in hash table */
}; /* data */

int make_data(data **data, options *opt);
int load_data(data *dat, options *opt, int tbd);
int sync_data(data *dat);
void fprint_fasta(FILE *fp, char_t *data, size_t n, size_t p, char const * const prefix);
void fprint_alignment2(FILE *fp, char_t **data, size_t n, size_t p);
#ifdef USE_CURSES
void wprint_fasta(WINDOW *wp, char_t *data, size_t n, size_t p, char const * const prefix);
void wprint_alignment2(WINDOW *wp, char_t **data, size_t n, size_t p);
#endif
int amplicon_read_compare(data *dat, size_t i, size_t j);
int update_data(data *dat, options *opt, size_t sample_size,size_t *sample_idx);
void free_data(data *dat);


#endif
