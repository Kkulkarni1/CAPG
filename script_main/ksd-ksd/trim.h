/**
 * @file ampliclust.h
 * @author Karin S. Dorman
 *
 * Header file for trim structs, extern functions and defines.
 */

#ifndef __H_TRIM__
#define __H_TRIM__

#include <stdio.h>

#include "fastq.h"
#include "error.h"

/**
 * Get rid of unused parameter warnings.
 */
#define UNUSED(x) (void)(x)

typedef struct _options options;
typedef struct _data data;
typedef struct _model model;

/**
 * Where to cache results.  Is this just the current solution,
 * the best local solution, or the best global solution?
 */
enum {
	DEFAULT,	/*<! cluster_id, criterion, etc. */
	BEST,		/*<! best_cluster_id, best_criterion, etc. */
	OPTIMAL		/*<! optimal_cluster_id, etc. */
};

/**
 * EM cycle.
 */
enum {
	FIRST,
	SECOND
};

/**
 * Error list.
 */
enum {
	EM_NO_ERROR,			/*<! no EM errors */
	EM_EXCEED_MAX_ITERATIONS,	/*<! over maximum allowed iterations */
	EM_ASCENT_VIOLATION,		/*<! EM log likelihood decline */
};

enum {
	TRIM_NO_ERROR,
	TRIM_AMBIGUOUS_NUCLEOTIDE
};

enum {
	ESTIMATE_PI = 1,	/*<! estimate mixing proportions */
	ESTIMATE_ALL = 1	/*<! estimate everything */
};

/**
 * Run options.
 */
struct _options {
	/* data */
	char const *fastq_file;	/*<! name of fastq input file */
	char const *in_primer;	/*<! primer from command line */
	char_t *primer;		/*<! primer sequence, with ambiguities */
	char_t *cprimer;	/*<! current primer */
	size_t n_primers;	/*<! no. sequences allowed by options.primer */
	unsigned int primer_len;/*<! primer length */

	/* running */
	unsigned int n_init;	/*<! number of random initializations */
	unsigned int n_iter;	/*<! max. number of AECM iterations */
	unsigned long seed;	/*<! random number seed [srand()] */
	double epsilon;		/*<! change in relative log likelihood */
	double tolerance;	/*<! small number allowed numerical error */
	int info;		/*<! level of information to output */
	int estimate;		/*<! which parameters to estimate */

	/* model choice */
	unsigned int n_indices;	/*<! number of possible start indices */
	size_t n_states;	/*<! no. possible hidden states */

	const char *trimmed_outfile;	/*<! fastq outfile for trimmed reads */
	const char *trim_outfile;	/*<! fastq outfile for trimmed part */
	const char *outfile;	/*<! ... */
}; /* options */

/**
 * Data object
 */
struct _data {
	fastq_data *fqd;		/*<! fastq data */
	char_t **reads;			/*<! pointers to individual reads */
	unsigned int *index;		/*<! primer start index */
	unsigned int *best_index;
	unsigned int *optimal_index;
}; /* data */

/**
 * Model object
 */
struct _model {
	double ll;		/*<! current log likelihood */
	double nll;		/*<! next log likelihood */
	double best_ll;		/*<! best log likelihood so far */
	unsigned int iter;	/*<! no. of realized iterations */

	double *eik;		/*<! expectations */

	/* parameters */
	double *pi;		/*<! kx1 mixing proportions */
	double *delta;		/*<! lx1 prob. error per position */
	double *gamma;		/*<! 4x4 substitution probabilities */
	double *eta;		/*<! sample sequence nucleotide probabilities */
	double *beta;		/*<! barcode nucleotide probabilities */
	double *lambda0;	/*<! qx1 quality score pmfs */
	double *lambda1;	/*<! qx1 quality score pmfs */

	/* update parameters */
	double *npi;
	double *ndelta;
	double *ngamma;
	double *neta;
	double *nbeta;
	double *nlambda0;
	double *nlambda1;

	/* best solution */
	char *best_haplotypes;
	double *best_pi;
	double *best_delta;
	double *best_eta;
	double *best_beta;
	double *best_gamma;
	double *best_lambda0;
	double *best_lambda1;
}; /* model */

int em(data *dat, model *mod, options *opt);
double e_step(data *dat, model *mod, options *opt);
void m_step(data *dat, model *mod, options *opt);
int param_update(model *mod, options *opt, int best);
void assign_index(data *dat, model *mod, options *opt, int type);
char *em_error_message(int err_no);
char *trim_error_message(int err_no);


void print_double(FILE *fp, double *v, size_t k);
void print_gamma(FILE *fp, double *gamma, size_t n);
void print_lambda(FILE *fp, double *gamma, size_t n, size_t l);
void print_uints(FILE *fp, size_t *v, size_t n, int width);

/**
 * Compare primer and read nucleotides.  Primer nucleotide is encoded as
 * iupac_t, while read nucleotide may not be.  Primer nucleotide should
 * not be ambiguous nucleotide.
 *
 * @param nuc	primer nucleotide
 * @param dat	data object pointer
 * @param i	read index
 * @param j	read position index
 * @return	-1, 0, 1
 */
inline int nuccmp(iupac_t nuc, data *dat, size_t i, size_t j) {
	if (popcnt[(int)nuc] > 1)
		exit(mmessage(WARNING_MSG, TRIM_AMBIGUOUS_NUCLEOTIDE,
			"comparing read nucleotide to %s in promoter: %c (%d "
			"%u)\n", trim_error_message(TRIM_AMBIGUOUS_NUCLEOTIDE),
			iupac_to_char[(int)nuc], nuc, popcnt[(int)nuc]));
	if (dat->fqd->read_encoding == IUPAC_ENCODING) {
		if (nuc < dat->reads[i][j]) return -1;
		if (nuc > dat->reads[i][j]) return 1;
	} else {
		if (nuc < xy_to_iupac[(int)dat->reads[i][j]]) return -1;
		if (nuc > xy_to_iupac[(int)dat->reads[i][j]]) return 1;
	}
        return 0;
} /* nuccmp */


/**
 * Return index of transition from primer to read nucleotide.  Primer nucleotide
 * is encoded as iupac_t, while read nucleotide may not be.  Neither nucleotide
 * should be ambiguous.  Use this function to avoid remembering nucleotide order.
 * Also see nuc_idx().
 *
 * @param nuc	primer nucleotide
 * @param dat	data object pointer
 * @param i	read index
 * @param j	read position index
 * @return	index
 */
inline size_t err_idx(iupac_t nuc, data *dat, size_t i, size_t j) {
	if (dat->fqd->read_encoding == IUPAC_ENCODING)
		return iupac_to_xy[(int)nuc] * NUM_NUCLEOTIDES + iupac_to_xy[(int)dat->reads[i][j]];
	else
		return iupac_to_xy[(int)nuc] * NUM_NUCLEOTIDES + dat->reads[i][j];
} /* err_idx */

/**
 * Return index of read nucleotide.  Nucleotide should not be ambiguous.  Use
 * this function to avoid remembering nucleotide order.  Also see err_idx().
 *
 * @param dat	data object pointer
 * @param i	read index
 * @param j	read position index
 * @return	index
 */
inline size_t nuc_idx(data *dat, size_t i, size_t j) {
	if (dat->fqd->read_encoding == IUPAC_ENCODING)
		return iupac_to_xy[(int)dat->reads[i][j]];
	else
		return dat->reads[i][j];
} /* nuc_idx */

/*
 * Return nucleotide char of nucleotide in read i, position j.
 *
 * @param dat	data object pointer
 * @param i	read index
 * @param j	read position index
 * @return	index
 */
inline char nuc_char(data *dat, size_t i, size_t j) {
	if (dat->fqd->read_encoding == IUPAC_ENCODING)
		return xy_to_char[(int)iupac_to_xy[(int)dat->reads[i][j]]];
	else
		return xy_to_char[(int)dat->reads[i][j]];
} /* nuc_char */

/**
 * Return first unambiguous sequence matching ambiguous sequence.  Given \p seq,
 * populate \p cseq with the first nucleotide sequence that lacks any of the
 * ambiguities in \p seq.  Also see next_sequence().
 *
 * @param cseq	stores resulting unambiguous sequence
 * @param seq	stores original sequence with ambiguities
 * @param len	length of sequence
 */
inline void first_sequence(char_t * const cseq, char_t * const seq, size_t len) {
	for (size_t i = 0; i < len; ++i) {
		/* ambiguous nucleotides: assign first possible nucleotide */
		if (popcnt[(int)seq[i]] > 1 && popcnt[(int)seq[i]] < NUM_NUCLEOTIDES)
			for (size_t n = 0; n < NUM_NUCLEOTIDES; ++n) {
				if (1L << n & seq[i]) {
					cseq[i] = 1L << n;
					break;
				}
			}
		else
			cseq[i] = seq[i];
	}
} /* first_sequence */

/**
 * Return next unambiguous sequence matching ambiguous sequence.  Given \p seq,
 * and \p cseq, update \p cseq to contain the next unambiguous sequence 
 * consistent with \p seq.  Return 0 if cycle back around to first such
 * sequence.  See first_sequence().
 *
 * @param cseq	stores resulting unambiguous sequence
 * @param seq	stores original sequence with ambiguities
 * @param len	length of sequence
 * @return	1 if cycles back to first unambiguous sequence, ow 0
 */
inline int next_sequence(char_t * const cseq, char_t * const seq, size_t len) {
	size_t n;
	for (size_t i = len - 1; 1; --i) {
		if (popcnt[(int)seq[i]] == 1 || popcnt[(int)seq[i]] == NUM_NUCLEOTIDES) {
			if (!i)
				return 0;
			else
				continue;
		}
		/* assign to next allowable nucleotide */
		for (n = iupac_to_std[(int)cseq[i]] + 1; n < NUM_NUCLEOTIDES; ++n)
			if (1L << n & seq[i]) {
				cseq[i] = 1L << n;
				break;
			}
		/* if no more allowable nucleotides, assign to first */
		if (n == NUM_NUCLEOTIDES) {
			for (size_t n = 0; n < NUM_NUCLEOTIDES; ++n)
				if (1L << n & seq[i]) {
					cseq[i] = 1L << n;
					break;
				}
		/* a change has been wrought: return */
		} else
			break;

		/* we have reset to the original */
		if (!i)
			return 0;
	}
	return 1;
} /* next_sequence */

#endif
