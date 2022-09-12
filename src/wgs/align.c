#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "align.h"
#include "error.h"

enum {
	TRACEBACK_DIAGONAL = 1,
	TRACEBACK_LEFT = 2,
	TRACEBACK_UP = 3
};


unsigned const char bases[NUM_NUCLEOTIDES] = {'A', 'C', 'G', 'T'};
static double *qacgt = NULL;


/**
 * Read in quality-dependent null distribution from file.
 * Variable qacgt is allocated as row-order matrix containing probability of
 * each nucleotide (column) given each quality score (row), with values obtained
 * from the named file.
 *
 * @param filename	name of file containing data
 * @return		error status
 */
int setup_null(char const * const filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	unsigned int cnt_lines = 0;	/* aka number quality scores */
	int char_after = 0;
	char c = fgetc(fp);
	while (c != EOF) {
		char_after = 1;
		if (c == '\n') {
			++cnt_lines;
			c = fgetc(fp);
			if (!feof(fp))
				char_after = 0;
		}
		c = fgetc(fp);
	}
	if (!char_after)
		++cnt_lines;

	/* allocate matrix */
	qacgt = malloc(NUM_NUCLEOTIDES * cnt_lines * sizeof *qacgt);
	if (!qacgt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "qacgt");

	fseek(fp, 0, SEEK_SET);
	double *dptr = qacgt;
	for (unsigned int i = 0; i < cnt_lines; ++i) {
		unsigned int q;
		if (fscanf(fp, "%u %lf %lf %lf %lf", &q, dptr, dptr + 1,
						dptr + 2, dptr + 3) != 5)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"%s: line %u\n", filename, i);
		dptr += 4;
	}
	fclose(fp);

	return NO_ERROR;
} /* setup_null */

/**
 * Score alignment of a pair of read nucleotides.  Consider aligning nucleotides
 * s1 and s2 with quality scores q1 and q2 (equivalently error probabilities e1
 * and e2), then assuming the quality scores communicate the true probability of
 * error, an appropriate score is
 * S(s1,s2,q1,q2) 	= \log[ Pr(s1,s2|q1,q2)/Pr(s1,s2) ]
 *			= \log[ \sum_N Pr(s1,s2|q1,q1,N)P(N) ]
 *				- \log[ Pr(s1) ] - \log[ Pr(s2) ]
 * We assume P(N) = 0.25 and
 * P(N) = P(s1) = P(s2) = 0.25 when qacgt [see setup_null()] is not loaded.
 * For Pr(s1,s2|q1,q2,N), we assume reads are independent and uniform
 * substitutions upon error.
 *
 * WARNING: There is a global variable in play.  If qacgt is set, then we use
 * Pr(s1|q1) in place of Pr(s1), where Pr(s1|q1) is provided by the user, read
 * from file via setup_null().
 *
 * @param s1	read 1 nucleotide
 * @param s2	read 2 nucleotide
 * @param q1	read 1 quality score
 * @param q2	read 2 quality score
 * @param e1	read 1 error probability, most likely 10^(-q1/10)
 * @param e2	read 2 error probability, most likely 10^(-q1/10)
 * @param add	shift log score up by this factor
 * @param enc	read encoding
 * @return	the log score
 */
double dynamic_scoring(unsigned char s1, unsigned char s2, unsigned char q1,
		unsigned char q2, double e1, double e2, double add, int enc)
{
	double score = 0, lscore;

	/* sum over possible source nucleotides */
	for (unsigned char b = 0; b < NUM_NUCLEOTIDES; ++b)
		//score += (b == s1 ? 1 - e1 : e1 / 3) * (b == s2 ? 1 - e2 : e2 / 3);
		//score += sub_prob_given_q_with_encoding(b, s1, XY_ENCODING, enc, q1, 0) * sub_prob_given_q_with_encoding(b, s2, XY_ENCODING, enc, q2, 0);
		score += sub_prob_with_encoding(b, s1, XY_ENCODING, enc, e1, 0) * sub_prob_with_encoding(b, s2, XY_ENCODING, enc, e2, 0);

	if (!qacgt)
		lscore = log(score * 4) + add;
	else
		lscore = log(score / 4)
			- log(qacgt[q1 * NUM_NUCLEOTIDES + s1])
			- log(qacgt[q2 * NUM_NUCLEOTIDES + s2]);

//if (!isfinite(lscore)) fprintf(stderr, "S(%c, %c, %u, %u, %f, %f): %e = %e - %e - %e ( = %e)\n", xy_to_char[s1], xy_to_char[s2], q1, q2, e1, e2, lscore, log(score/4), log(qacgt[q1 * NUM_NUCLEOTIDES + s1]), log(qacgt[q2 * NUM_NUCLEOTIDES + s2]), log(score/4) - log(qacgt[q1 * NUM_NUCLEOTIDES + s1]) -  log(qacgt[q2 * NUM_NUCLEOTIDES + s2]));

	return lscore;
} /* dynamic_scoring */


/**
 * Compute Z(T) via forward-sum variant of Viterbi.
 *
 * @param s1	first sequence
 * @param s2	second sequence
 * @param len1	length of first sequence
 * @param len2	length of first sequence
 * @param score	4x4 matrix of substitution scores
 * @param gap_p	gap penalty
 * @param band	band, -1 for no band
 * @param semi	semi-global alignment
 * @param perr1	error probability in read 1
 * @param perr2	error probability in read 2
 * @param q1	quality scores in read 1
 * @param q2	quality scores in read 2
 * @param err	pointer to integer error status
 * @param alen	pointer to alignment length
 * @param ascr	pointer to alignment score
 * @param add	how much to shift up score per aligned position
 * @param enc	read encoding
 * @return	alignment
 */
unsigned char **nw_alignment(unsigned char const * const s1,
	unsigned char const * const s2, int len1, int len2,
	int score[4][4], double gap_p, int band, int semi, double const *perr1,
	double const *perr2, unsigned char const * const q1,
	unsigned char const * const q2, int *err, unsigned int *alen,
	double *ascr, double add, int enc)
{
	static size_t nnw = 0;
	int i, j;
	int l, r;
	int iband = band >= 0 ? band : 0;
	double diag, left, up;

	*err = NO_ERROR;

	unsigned int nrow = len1 + 1;
	unsigned int ncol = len2 + 1;

	/* allocate matrices, including traceback */
	double *d = (double *) malloc(nrow * ncol * sizeof(double)); 
	int *p = (int *) malloc(nrow * ncol * sizeof(int));
	if (d == NULL || p == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "d & p");
		return NULL;
	}

	// Fill out left columns of d, p.
	for (i = 0; i <= len1; i++) {
		d[i*ncol] = semi ? 0 : i * gap_p;
		p[i*ncol] = 3;
	}

	// Fill out top rows of d, p.
	for (j = 0; j <= len2; j++) {
		d[j] = semi ? 0 : j * gap_p;
		p[j] = 2;
	}

	// Calculate left/right-bands in case of different lengths
	int lband = iband, rband = iband;
	if (len2 > len1)
		rband = iband + len2 - len1;
	else if (len1 > len2)
		lband = iband + len1 - len2;

	// Fill out band boundaries of d.
	if (band >= 0 && (iband < len1 || iband < len2)) {
		for (i = 0; i <= len1; i++) {
			if (i - lband - 1 >= 0)
				d[i*ncol + i - lband - 1] = -9999;
			if (i + rband + 1 <= len2)
				d[i*ncol + i + rband + 1] = -9999;
		}
	}

	if (ascr)
		*ascr = 0;

	// Fill out the body of the DP matrix.
	for (i = 1; i <= len1; i++) {
		if (band >= 0) {
			l = i - lband;
			if (l < 1)
				l = 1;
			r = i + rband;
			if (r > len2)
				r = len2;
		} else {
			l = 1;
			r = len2;
		}

		for (j = l; j <= r; j++) {
			// Score for the left move.
			if (i == len1)
				left = (double) d[i*ncol + j - 1]
					+ (semi ? 0 : gap_p);
			else
				left = (double) d[i*ncol + j - 1] + gap_p;

			// Score for the up move.
			if (j == len2)
				up = d[(i-1)*ncol + j]
					+ (semi ? 0 : gap_p);
			else
				up = d[(i-1)*ncol + j] + gap_p;

			// Score for the diagonal move.
			diag = d[(i-1)*ncol + j-1] + (perr1
				? dynamic_scoring(s1[i-1], s2[j-1], q1[i-1],
				q2[j-1], perr1[i-1], perr2[j-1], add, enc)
				: score[(int) s1[i-1]][(int) s2[j-1]]);

			// Break ties and fill in d,p.
			if (up >= diag && up >= left) {
				d[i*ncol + j] = up;
				p[i*ncol + j] = 3;
				if (ascr)
					*ascr = up;
			} else if (left >= diag) {
				d[i*ncol + j] = left;
				p[i*ncol + j] = 2;
				if (ascr)
					*ascr = left;
			} else {
				d[i*ncol + j] = diag;
				p[i*ncol + j] = 1;
				if (ascr)
					*ascr = diag;
			}
		}
	}

	unsigned char *al0 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	unsigned char *al1 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	if (al0 == NULL || al1 == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al0 & al1");
		return NULL;
	}

	// Trace back over p to form the alignment.
	int len_al = 0;
	i = len1;
	j = len2;

	while ( i > 0 || j > 0 ) {
		switch ( p[i*ncol + j] ) {
			case 1:
				al0[len_al] = s1[--i];
				al1[len_al] = s2[--j];
				break;
			case 2:
				al0[len_al] = '-';
				al1[len_al] = s2[--j];
				break;
			case 3:
				al0[len_al] = s1[--i];
				al1[len_al] = '-';
				break;
			default:
				*err = OUT_OF_BAND;
				mmessage(ERROR_MSG, *err,
					"NW alignment out of range");
				return NULL;
		}
		len_al++;
	}
//	al0[len_al] = '\0';
//	al1[len_al] = '\0';
/*
	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j)
			fprintf(stderr, " %2d", d[i*ncol + j]);
		fprintf(stderr, "\n");
	}
*/


	// Allocate memory to alignment strings.
	unsigned char **al = (unsigned char **) malloc( 2 * sizeof(unsigned char *) ); //E
	if (al == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al");
		return NULL;
	}

	al[0] = (unsigned char *) malloc(len_al); //E
	al[1] = (unsigned char *) malloc(len_al); //E
	if (al[0] == NULL || al[1] == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al[]");
		return NULL;
	}

	// Reverse the alignment strings (since traced backwards).
	for (i = 0 ; i < len_al ; i++) {
		al[0][i] = al0[len_al-i-1];
		al[1][i] = al1[len_al-i-1];
	}
//	al[0][len_al] = '\0';
//	al[1][len_al] = '\0';

	// Free allocated memory
	free(d);
	free(p);
	free(al0);
	free(al1);

	*alen = len_al;

	nnw++;
	return al;
} /* nw_alignment */
