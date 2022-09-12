#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include "error.h"
#include "fastq.h"

/**
 * Nelder-Mead errors.
 */
enum {
	OUT_OF_BAND = NUM_ERRORS + 1	/* outside band in banded alignment */
};

int setup_null(char const * const filename);
unsigned char **nw_alignment(unsigned char const * const s1, unsigned char const * const s2, int len1, int len2, int score[4][4], double gap_p, int band, int ends_free, double const *perr1, double const *perr2, unsigned char const * const q1, unsigned char const * const q2, int *err, unsigned int *alen, double *ascore, double score_add, int enc);
double dynamic_scoring(unsigned char s1, unsigned char s2, unsigned char q1, unsigned char q2, double e1, double e2, double score_add, int enc);

#endif
