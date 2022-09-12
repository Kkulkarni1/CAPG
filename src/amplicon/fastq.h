/**
 * @file fastq.c
 * @author Karin S. Dorman
 *
 * FASTA/FASTQ library.
 *
 * DATA ENCODING:
 * Because NGS data are big data, it can be important to consider how they
 * are internally represented.  The four nucleotides, A, C, G, and T, can be
 * stored in 2 bits, the IUPAC nucleotide codes in 4 (if we treat U and T as
 * equivalent).  FASTA/FASTQ files represent the nucleotides (and quality
 * scores) as 8 bit characters even though all 8 bits are not necessary, so
 * some of the characters do not correspond to valid codes.  Because the
 * printable ASCII characters are in the mid-range of the 8-bit integers, it is
 * also possible to subtract 'A' (the integer interpretation of character 'A')
 * from the FASTA/FASTQ characters to encode the nucleotides as integers in the
 * range 'A' - 'A' (0) to 'Y' - 'A' (24).
 *
 * When reading FASTA/FASTQ files, we convert to the 4-bit encoding (iupac_t),
 * or the 2-bit encoding (xy_t).  Only the former can encode ambiguous
 * nucleotides such as R, Y, or N.
 *
 * MORE EFFICIENT ENCODING FOR FUTURE:
 * Despite all this fancy stuff about encoding, the current library STILL STORES
 * the nucleotides and quality scores in 8 bits, thus wasting at least 4 bits of
 * unused memory.  However, see sequence.h and sequence.c for a memory-efficient
 * encoding of sequence data that we could use at some future date.  ANOTHER
 * IDEA: store nucleotides (2 bits) in the low bits, quality scores (6 bits) in
 * the high bits.
 *
 * VARIABLE READ LENGTHS:
 * This implementation handles variable length reads (fastq_data:n_lengths), but
 * does not provide pointers to the start of each read.  There is a slight
 * saving of memory because read lengths are assumed to be rather short
 * (unsigned int), while pointers are large.  Unless you store pointers to the
 * start of each read in fastq_data:reads or qualities in fastq:quals, you will
 * have to use pointer arithmetic to find individual read data.  If reads are
 * all the same length, then fastq_data:n_lengths is NULL and
 * fastq_data:n_max_length is the shared length.  This inconsistent interface
 * to the read length is a nuisance, so see read_length() for access.
 *
 * DATA TYPES AND LIMITS:
 * This implementation assumes the number and lengths of reads can be stored in
 * unsigned int, and there is no verification of this truth.  Be cautious.
 */

#ifndef __H_FASTQ_DATA__
#define __H_FASTQ_DATA__

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "nuc.h"
#include "qual.h"
#include "error.h"
//#include "constants.h"

/**
 * Type for character data (reads and qualities).
 */
typedef uint8_t char_t;

/**
 * Types of files.
 */
enum {
	FASTQ_FILE,
	FASTA_FILE,
	NUM_FILE_TYPES
};

extern char const *file_type_name[NUM_FILE_TYPES];

/**
 * Generated errors.
 */
enum {
	FASTQ_INVALID_READ_CHAR		/*!< invalid character in read */
		= NUM_ERRORS + 1,
	FASTQ_AMBIGUOUS_READ_CHAR,	/*!< ambiguous character in read */
	FASTQ_INVALID_QUALITY_CHAR,	/*!< invalid quality score (NOT USED) */
	FASTQ_INCOMPLETE_READ,		/*!< file ends before end of read */
	FASTQ_FILE_FORMAT_ERROR,	/*!< other file format error */
	FASTQ_PREMATURE_EOF,		/*!< premature end of fast[qa] file */
	FASTQ_EOF,			/*!< natural end of fast[qa] file */
	FASTQ_READ_TOO_SHORT,		/*!< read too short as per filter */
	FASTQ_READ_TOO_LONG,		/*!< read too long as per filter */
	FASTQ_NUM_ERRORS		/*!< number of possible errors */
};

typedef struct _fastq_data fastq_data;
typedef struct _fasta_data fasta_data;
typedef struct _fastq_options fastq_options;

/**
 * FASTQ data structure.
 */
struct _fastq_data {
	int empty;			/*!< whether there is data loaded */
	int file_type;			/*!< type of file */
	int read_encoding;		/*!< read encoding used */
	unsigned int n_reads;		/*!< number of reads */
	unsigned int n_max_length;	/*!< length of reads if identical */
	unsigned int n_min_length;	/*!< length of shortest read */
	size_t *index;			/*!< byte pointer of reads in file */
	unsigned int *n_lengths;	/*!< length of reads if not identical */
	unsigned int *name_lengths;	/*!< lengths of names */
	unsigned int *read_flag;	/*!< possible filtering read */
	unsigned char *site_flag;	/*!< possible filtering per site */
	char_t *reads;			/*!< sequence reads */
	char_t *quals;			/*!< quality score strings */
	char *names;			/*!< headers */
	char_t *reference_seq;		/*!< optional reference sequence */
	char_t max_quality;		/*!< maximum observed quality score */
	char_t min_quality;		/*!< minimum observed quality score */
};

/**
 * FASTQ options.
 */
struct _fastq_options {
	int drop_ambiguous_reads;	/*!< drop invalid reads */
	int drop_ambiguous_nucs;	/*!< drop nucleotides */
	char const *outfile;		/*!< outfile */
	int read_names;			/*!< read in read names */
	int paired;			/*!< paired with another file: keep pairing */
	int append;			/*!< append to outfile */
	int read_encoding;		/*!< which read encoding to request */
	int reverse_complement;		/*!< reverse complement the data */
	int fasta;			/*!< convert reads to fasta format */
	int qs_fasta;			/*!< convert qualities to fasta */
	unsigned int min_length;	/*!< min. len. to read/write sequence */
	unsigned int max_length;	/*!< max. len. to read/write sequence */
	unsigned char casavize;		/*!< append casava-style naming to names */
};

/* read fasta/fastq files */
int allocate_empty_fastq(fastq_data **in_fqd, fastq_options *fqo, unsigned int nreads, unsigned int read_length);
int fread_fastq(FILE *fp, fastq_data **fqd, fastq_options *fqo);
int read_fastq(const char *filename, fastq_data **fqd, fastq_options *fqo);
int read_read(FILE *fp, fastq_data *fqd, fastq_options *fqo, unsigned int *len, char_t *rptr, char_t *qptr, char *nptr, unsigned int *nlen);
unsigned int cnt_reads(char const * const filename);
unsigned int fcnt_reads(FILE *fp);
int allocate_read_flag(fastq_data *fqd);
int allocate_site_flag(fastq_data *fqd);
int findex_reads(FILE *fp, fastq_data *fqd);
int make_fastq_options(fastq_options **opt);


/* do stuff with reads */
int read_compare(fastq_data *fqd, unsigned int i, unsigned int j);
int pw_align_reads(fastq_data *fqd, char const * const rfile);
double read_distance(fastq_data *fqd, unsigned int i, unsigned int j);
double read_distance_ptr(fastq_data *fqd, unsigned int len, char_t *rptr1, char_t *rptr2, char_t *qptr1, char_t *qptr2);
unsigned char *reverse_complement(char_t const *const in_str, unsigned int len, int encoding);	/* caller must free return pointer */
void in_situ_reverse_complement(char_t *const str, unsigned int len, int encoding);
void in_situ_reverse(char_t * const str, unsigned int len);

/* output */
unsigned char *display_sequence(char_t const * const in_str, unsigned int len, int encoding);	/* caller must free return pointer */
unsigned char *display_sequence_trimmed(char_t const * const in_str, unsigned int len, int encoding, unsigned char *trim);	/* caller must free return pointer */
unsigned char *display_reverse_complement(char_t const * const in_str, unsigned int len, int encoding);	/* caller must free return pointer */
unsigned char *display_reverse_complement_trimmed(char_t const * const in_str, unsigned int len, int encoding, unsigned char *trim);	/* caller must free return pointer */
unsigned char *display_quals(char_t const * const in_str, unsigned int len, char_t min);	/* caller must free return pointer */
unsigned char *display_reverse_quals(char_t const * const in_str, unsigned int len, char_t min);	/* caller must free return pointer */
unsigned char *display_reverse_quals_trimmed(char_t const * const in_str, unsigned int len, char_t min, unsigned char *trim);	/* caller must free return pointer */
unsigned char *display_reverse_quals_trimmed(char_t const * const in_str, unsigned int len, char_t min, unsigned char *trim);	/* caller must free return pointer */
void write_sequence(FILE *fp, char_t const * const in_str, unsigned int len, int encoding);
void write_quals(FILE *fp, char_t const * const in_str, unsigned int len, unsigned char min);
int write_fastq(fastq_data *fqd, fastq_options *fqo);
int write_fastq_marked_trimmed(fastq_data *fqd, fastq_options *fqo, unsigned int *id, unsigned int selected_id, int (*trim)(fastq_data *, char_t *, char_t *, unsigned int, unsigned int, void *), void *obj);
int write_table(fastq_data *fqd, char const *filename);
int write_table_marked(fastq_data *fqd, char const *filename, unsigned int *id, unsigned int val);
char const *read_name(fastq_data *fqd, unsigned int j);
char const *next_read_name(fastq_data *fqd, char const **names, unsigned int j);
char *next_read_name_rw(fastq_data *fqd, char **names, unsigned int j);

const char *fastq_error_message(int err_no);

void free_fastq(fastq_data *fqd);
void free_fastq_options(fastq_options *opt);


/**
 * Convert quality code to error probability.
 *
 * @param q	quality code encoded as ASCII
 * @param fqd	fastq_data object pointer
 * @return	probability
 */
inline double error_prob(fastq_data *fqd, char_t q) {
	return pow(10, - (q + fqd->min_quality - 33) / 10.);
} /* error_prob */

/**
 * Print out true error probabilities.
 *
 * @param fp	FILE pointer
 * @param fqd	fastq_data object pointer
 */
inline void fprint_error_probs(FILE *fp, fastq_data *fqd)
{
	unsigned int n_quality = fqd->max_quality - fqd->min_quality + 1;
	for (unsigned int q = 0; q < n_quality; ++q)
		fprintf(fp, " %f", error_prob(fqd, q));
} /* fprint_error_probs */

/**
 * Return length of given read.
 *
 * @param fqd	fastq_data object pointer
 * @param i	index of read
 * @return	length of read
 */
inline unsigned int read_length(fastq_data *fqd, unsigned int i)
{
	return fqd->n_lengths ? fqd->n_lengths[i] : fqd->n_max_length;
} /* read_length */

/**
 * Write one read of a given length as a row in R table format.
 *
 * @param fp	file handle pointer
 * @param read	read to write
 * @param len	length of read
 */
inline void write_read_in_table(FILE *fp, char_t *read, unsigned int len)
{
	for (unsigned int j = 0; j < len; ++j) {
		if (j)
			fprintf(fp, " ");
		fprintf(fp, "%u", (unsigned int) read[j]);
	}
	fprintf(fp, "\n");
} /* write_read_in_table */

inline char nuc(fastq_data *fqd, char_t c)
{
	if (fqd->read_encoding == IUPAC_ENCODING)
		return iupac_to_char[c];
	else if (fqd->read_encoding == XY_ENCODING)
		return xy_to_char[c];
	else
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Call programmer!\n");
}/* nuc */

#endif
