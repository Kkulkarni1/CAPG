/*
 * @file merger.c
 * @author K. S. Dorman
 *
 * Purpose: read merger
 * 
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MATHLIB_STANDALONE
#include <R/Rmath.h>

#include "align.h"
#include "error.h"
#include "fastq.h"

char *char_reverse_complement(char const *read)
{
	int len = strlen(read);
	char *rc = malloc((len + 1) * sizeof *read);

	if (!rc)
		return rc;

	rc[len] = '\0';
	for (int i = 0; i < len; ++i)
		if (read[len - 1 - i] == 'A')
			rc[i] = 'T';
		else if (read[len - 1 - i] == 'C')
			rc[i] = 'G';
		else if (read[len - 1 - i] == 'G')
			rc[i] = 'C';
		else if (read[len - 1 - i] == 'T')
			rc[i] = 'A';

	return rc;
} /* char_reverse_complement */

int main(int argc, const char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Invalid command line: ./merger fastq1 fastq2 "
			"out_root\n");
		exit(1);
	}
	int score[NUM_NUCLEOTIDES][NUM_NUCLEOTIDES] = {{2, -3, -3, -2},
			{-3, 2, -2, -3}, {-3, -2, 2, -3}, {-2, -3, -3, 2}};
	int err = NO_ERROR;
	double gap_penalty = -7;//-6.5;-INFINITY;//	/* gap open/extension penalty */
	double score_add = 0;//0.5;	/* shift scores by this much */
	int band = -1;			/* no banded alignment w/ semi-global */
	int semi_global = 1;		/* semi-global align.: don't change */
	int use_null = 0;		/* read and use null distn */
	int alignment_to_stderr = 1;	/* alignments to stderr */
	double avg_mm_ep = 0;		/* average mismatch penalty */
	double avg_match_ep = 0;	/* average match penalty */
	unsigned int n_mm_ep = 0, n_match_ep = 0;	/* counts for average */
	unsigned int n_indel = 0, n_alen = 0;
	unsigned int n_match = 0, n_mismatch = 0;
	fastq_data *fqd1 = NULL;	/* forward reads */
	fastq_data *fqd2 = NULL;	/* reverse reads */
	fastq_options *fqo = NULL;	/* options for fastq files */
	char const *fastq_file1 = argv[1];	/* forward read file */
	char const *fastq_file2 = argv[2];	/* reverse read file */
	unsigned int max_nreads = 100;	/* debugging: quit after these reads */
	char const *out_base = argv[3];		/* base name for output files */
	char *outfile1 = NULL, *outfile2 = NULL;
	FILE *out1 = NULL, *out2 = NULL;
	double *err1 = NULL, *err2 = NULL;	/* error probabilities */
	double *aquals = NULL;
	int read_encoding = XY_ENCODING;/* [TODO] XY_ENCODING assumed below */
	int degap_readthrough = 1;	/* automatically trim gaps  -----RRRRRRR
					 * from this configuration: RRRRRRRRR---
					 * as probable read-through.
					 */

	if (!strlen(out_base)) {
		fprintf(stderr, "Invalid command line: ./merger fastq1 fastq2 "
								"out_root\n");
		goto EXIT_NOW;
	}

	/* set up fastq options */
	if ((err = make_fastq_options(&fqo)))
		goto EXIT_NOW;

	fqo->read_names = 1;			/* store read names */
	fqo->read_encoding = read_encoding;	/* disallow ambiguous nucs */
	fqo->drop_ambiguous_nucs = 1;		/* by dropping them */
	fqo->paired = 1;			/* reads are paired */

	outfile1 = malloc((strlen(out_base) + 19) * sizeof *outfile1);
	outfile2 = malloc((strlen(out_base) + 19) * sizeof *outfile2);
	sprintf(outfile1, "%s_unmerged_R1.fastq", out_base);
	sprintf(outfile2, "%s_unmerged_R2.fastq", out_base);
	out1 = fopen(outfile1, "w");
	out2 = fopen(outfile2, "w");

	if (!out1 || !out2)
		goto EXIT_NOW;

	/* read in fastq files */
	if ((err = read_fastq(fastq_file1, &fqd1, fqo)))
		goto EXIT_NOW;
	if ((err = read_fastq(fastq_file2, &fqd2, fqo)))
		goto EXIT_NOW;

	/* read null distribution for unaligned positions */
	if (use_null && (err = setup_null("acgt_by_q.Rtxt")))
		goto EXIT_NOW;

	unsigned char *rptr1 = fqd1->reads;
	unsigned char *rptr2 = fqd2->reads;
	unsigned char *qptr1 = fqd1->quals;
	unsigned char *qptr2 = fqd2->quals;
	char *nptr = fqd1->names;

	err1 = malloc(fqd1->n_max_length * sizeof *err1);
	err2 = malloc(fqd2->n_max_length * sizeof *err2);

	if (!err1 || !err2)
		goto EXIT_NOW;

	for (unsigned int i = 0; i < fqd1->n_reads; ++i) {
		unsigned int len1 = read_length(fqd1, i);
		unsigned int len2 = read_length(fqd2, i);
		unsigned int alen;
		double ascore;

		/* reverse reverse reads */
		in_situ_reverse_complement(rptr2, len2, fqd2->read_encoding);
		in_situ_reverse(qptr2, len2);

		/* compute error probabilities */
		for (unsigned j = 0; j < len1; ++j)
			err1[j] = error_prob(fqd1, qptr1[j]);
		for (unsigned j = 0; j < len2; ++j)
			err2[j] = error_prob(fqd2, qptr2[j]);

#ifdef DEBUGGING
		unsigned char *str1 = display_sequence(rptr1, len1,
								XY_ENCODING);
		unsigned char *str2 = display_sequence(rptr2, len2,
								XY_ENCODING);
		unsigned char *qstr1 = display_quals(qptr1, len1,
							fqd1->min_quality);
		unsigned char *qstr2 = display_quals(qptr2, len2,
							fqd2->min_quality);
		fprintf(stderr, "=============================================="
			"======================================================"
			"\n@R1\n%s\n+R1\n%s\n@R2\n%s\n+R2\n%s\n================"
			"======================================================"
			"==============================\n", str1, qstr1, str2,
									qstr2);
		free(str1);
		free(str2);
		free(qstr1);
		free(qstr2);
#endif

		/* get alignment */
		unsigned char **aln = nw_alignment(rptr1, rptr2, len1, len2,
			score, gap_penalty, band, semi_global, err1, err2,
			qptr1, qptr2, &err, &alen, &ascore, score_add,
								read_encoding);

		/* not aligned */
		if (!ascore) {
			fprintf(out1, "@%.*s\n", fqd1->name_lengths[i], nptr);
			fprintf(out2, "@%.*s\n", fqd2->name_lengths[i], nptr);
			write_sequence(out1, rptr1, len1, read_encoding);
			write_sequence(out2, rptr2, len2, read_encoding);
			fprintf(out1, "\n+\n");
			fprintf(out2, "\n+\n");
			write_quals(out1, qptr1, len1, fqd1->min_quality);
			write_quals(out2, qptr2, len2, fqd2->min_quality);
			fprintf(out1, "\n");
			fprintf(out2, "\n");
			goto CONTINUE_LOOP;
		}

		aquals = malloc(alen * sizeof *aquals);
		if (!aquals)
			goto EXIT_NOW;

		if (alignment_to_stderr)
			fprintf(stderr, "@%.*s score=%.3e %f vs. %f vs. %f\n",
				fqd1->name_lengths[i], nptr, ascore, avg_mm_ep,
						gap_penalty, avg_match_ep);
		unsigned start = 0, end = 0;
		unsigned in1 = 0, in2 = 0;
		unsigned int i1 = 0, i2 = 0;
		int readthrough = 0;
		for (unsigned int j = 0; j < alen; ++j) {
			if (alignment_to_stderr)
				fprintf(stderr, "%c", aln[0][j] == '-' ? '-'
							: nuc(fqd1, aln[0][j]));
			if (!in1 && aln[0][j] != '-') {
				in1 = 1;
				if (in2) {
					readthrough = 1;
					start = j;
				}
			}
			if (!in2 && aln[1][j] != '-') {
				in2 = 1;
				if (in1)
					start = j;
			}
			if (!end && i1 == len1)
				end = j;
			else if (!end && i2 == len2)
				end = j;
			if (aln[0][j] != '-')
				++i1;
			if (aln[1][j] != '-')
				++i2;
		}

		if (alignment_to_stderr) {
			fprintf(stderr, "\n");
			for (unsigned int j = 0; j < alen; ++j)
				fprintf(stderr, "%c", aln[1][j] == '-' ? '-' :
					aln[1][j] == aln[0][j] ? '.' :
							nuc(fqd2, aln[1][j]));
			fprintf(stderr, "\n");
		}

		/* output fastq format */

		fprintf(stdout, "@%.*s score=%.3e (%u, %u) %s%f vs. %f vs. %f (%f = %u / %u; %f = %u / %u)\n",
			fqd1->name_lengths[i], nptr, ascore, start, end,
					readthrough ? "read-through " : "",
					avg_mm_ep, gap_penalty, avg_match_ep,
					(double) n_mismatch / (n_mismatch + n_match),
					n_mismatch, n_mismatch + n_match,
					(double) n_indel / n_alen, n_indel, n_alen);

		/* also compute overlap nucleotides and scores */
		i1 = i2 = 0;
		n_alen += 2 * (end - start);
		for (unsigned int j = 0; j < alen; ++j) {
			if (start <= j && j < end) {

				/* gap in alignment: Schirmer2016 indicates 
				 * deletions slightly more likelyx
				 */
				if (aln[0][j] == '-' || aln[1][j] == '-') {
					fprintf(stdout, "%c", nuc(fqd1,
						aln[0][j] == '-' ? aln[1][j]
								: aln[0][j]));
					aquals[j - start] = aln[0][j] == '-'
						? error_prob(fqd2, qptr2[i2])
						: error_prob(fqd1, qptr1[i1]);
					++n_indel;

				/* otherwise we must estimate true nucleotide */
				} else {
					aquals[j - start] = (!rptr1[i1] ?  1 -
						err1[i1] : err1[i1] / 3) *
						(!rptr2[i2] ? 1 - err2[i2]
								: err2[i2] / 3);
					double sum = aquals[j - start];
					unsigned char mpp_nuc = 0;
					for (unsigned int b = 1; b
						< NUM_NUCLEOTIDES; ++b) {
						double pp = (rptr1[i1] == b ?
							1 - err1[i1] : err1[i1]
							/ 3) * (rptr2[i2] == b ?
							1 - err2[i2] : err2[i2]
									/ 3);
						if (pp > aquals[j - start]) {
							aquals[j - start] = pp;
							mpp_nuc = b;
						}
						sum += pp;
					}
					fprintf(stdout, "%c",
							nuc(fqd1, mpp_nuc));
					aquals[j - start] = 1 - (rptr1[i1] ==
						mpp_nuc ? 1 - err1[i1] :
						err1[i1] / 3) * (rptr2[i2] ==
						mpp_nuc ? 1 - err2[i2] :
							err2[i2] / 3) / sum;

					/* read error */
					if (rptr1[i1] != rptr2[i2] &&
						(mpp_nuc == aln[0][j] || mpp_nuc
								== aln[1][j])) {
						/* equiv. to dynamic_scoring */
						avg_mm_ep = (n_mm_ep * avg_mm_ep
								+ log(sum * 4)
								+ score_add)
								/ (n_mm_ep + 1);
						++n_mm_ep;
					} else if (rptr1[i1] == rptr2[i2]
						&& mpp_nuc == rptr1[i]) {
						/* equiv. to dynamic_scoring */
						avg_match_ep = (n_match_ep *
								avg_match_ep
								+ log(sum * 4)
								+ score_add)
							/ (n_match_ep + 1);
						++n_match_ep;
					}
					if (rptr1[i1] == rptr2[i2])
						++n_match;
					else
						++n_mismatch;
/* appears to be the same
					if (rptr1[i1] != rptr2[i2] && aquals[j - start] > 0) {
						avg_mm_ep = (n_mm_ep*avg_mm_ep + log(aquals[j - start])) / (n_mm_ep + 1);
						++n_mm_ep;
					}
*/
				}
			} else if (j < start && !(degap_readthrough && readthrough)) {
				if (aln[0][j] != '-')
					fprintf(stdout, "%c", nuc(fqd1, aln[0][j]));
				else
					fprintf(stdout, "%c", nuc(fqd2, aln[1][j]));
			} else if (!(degap_readthrough && readthrough)) {
				if (aln[0][j] != '-')
					fprintf(stdout, "%c", nuc(fqd1, aln[0][j]));
				else
					fprintf(stdout, "%c", nuc(fqd2, aln[1][j]));
			}
			if (aln[0][j] != '-')
				++i1;
			if (aln[1][j] != '-')
				++i2;
		}
		fprintf(stdout, "\n+\n");

		i1 = i2 = 0;
		for (unsigned int j = 0; j < alen; ++j) {
			if (start <= j && j < end) {
				double qs = -10*log10(aquals[j - start]) + MIN_ASCII_QUALITY_SCORE;
				fprintf(stdout, "%c", (char)
					(qs > MAX_ASCII_QUALITY_SCORE ? MAX_ASCII_QUALITY_SCORE : qs));
			} else if (j < start && !(degap_readthrough && readthrough)) {
				if (aln[0][j] != '-')
					fprintf(stdout, "%c", (char)
						(qptr1[i1] + fqd1->min_quality));
				else
					fprintf(stdout, "%c", (char)
						(qptr2[i2] + fqd2->min_quality));
			} else if (!(degap_readthrough && readthrough)) {
				if (aln[0][j] != '-')
					fprintf(stdout, "%c", (char)
						(qptr1[i1] + fqd1->min_quality));
				else
					fprintf(stdout, "%c", (char)
						(qptr2[i2] + fqd2->min_quality));
			}
			if (aln[0][j] != '-')
				++i1;
			if (aln[1][j] != '-')
				++i2;
		}
		fprintf(stdout, "\n");

		free(aln[0]);
		free(aln[1]);
		free(aln);
		free(aquals);
CONTINUE_LOOP:
		rptr1 += len1;
		rptr2 += len2;
		qptr1 += len1;
		qptr2 += len2;
		nptr += fqd1->name_lengths[i];
		if (i + 1 == max_nreads)
			break;
	}

	err = EXIT_SUCCESS;

EXIT_NOW:

	if (out1)
		fclose(out1);
	if (out2)
		fclose(out2);
	if (outfile1)
		free(outfile1);
	if (outfile2)
		free(outfile2);

	if (fqd1)
		free_fastq(fqd1);
	if (fqd2)
		free_fastq(fqd2);
	if (fqo)
		free_fastq_options(fqo);
	if (err1)
		free(err1);
	if (err2)
		free(err2);
	
	exit(err);
} /* main */


/**
 * Likelihood ratio to evaluate the hypothesis of overlap.
 *
 * Assume iid sites conditional on quality scores with probability
 * Pr(X1=x1 | Q1=q1)Pr(X2=x2 | Q2=q2) when unaligned
 * where Pr(X=x | Q=q) = \sum_{n\in{A,C,G,T}} Pr(X=x | N=n,Q=q)Pr(N=n) != 1/4 if Pr(N=n) != 1/4.
 * Pr(X1=x1,X2=x2 | Q1=q1,Q2=q2) = \sum_{n\in{A,C,G,T}} Pr(X1=x1 | N=n,Q1=q1)Pr(X1=x1 | N=n,Q1=q1)Pr(N=n), but then we assumed Pr(N=n) = 1/4, which contradicts the above model.
 *
 * @param aln   the maximum likelihood alignment
 * @param s1    read 1
 * @param s2    read 2
 * @param q1    quality score 1
 * @param q2    quality score 21
 * @param alen  length of alignment
 * @param len1  length of read 1
 * @param len2  length of read 2
 */

