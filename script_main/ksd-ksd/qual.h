#ifndef __QUAL_H__
#define __QUAL_H__

#include "constants.h"

typedef unsigned char qual_t;

/*
extern const uint8_t _q_bits;
extern const uint8_t _q_mask;
*/
#define _q_bits	8U
#define _q_mask	255U

extern sequence_opt _qual_sequence_opt;

typedef double (*sub_prob_fxn)(data_t, data_t, data_t, void *);
typedef double (*sub_lprob_fxn)(data_t, data_t, data_t, void *);

extern sub_prob_fxn sub_prob;
extern sub_lprob_fxn sub_lprob;

#define MIN_QUALITY_SCORE	0	/*<! '!' - 33 */
#define MAX_QUALITY_SCORE	93	/*<! '~' - 33 */
#define MIN_ASCII_QUALITY_SCORE	33      /*<! '!' */
#define MAX_ASCII_QUALITY_SCORE	126	/*<! '~' */
#define MIN_ILLUMINA_QUALITY_SCORE	0	/*<! '!' - 33 */
#define MAX_ILLUMINA_QUALITY_SCORE	41	/*<! 'I' - 33 */

/**
 * Return probability of error at position \par i of read sequence \par s.
 */
inline double qual_to_prob(sequence *s, size_t i)
{
	return exp(- read_char(s, &_qual_sequence_opt, i) / 10. * log(10.));
} /* qual_to_prob */

/**
 * Return expected number of errors in quality sequence.
 */
inline double qual_to_exp_err(sequence *s)
{
	double ee = 0;
	for (size_t i = 0; i < s->len; ++i)
		ee += qual_to_prob(s, i);
	return ee;
} /* qual_to_exp_err */

inline double prob_error_free(sequence *s)
{
	double lef = 0;
	for (size_t i = 0; i < s->len; ++i)
		lef += log(1 - qual_to_prob(s, i));
	return exp(lef);
} /* prob_error_free */

/**
 * Return log probability of error at position \par i of read sequence \par s.
 */
inline double qual_to_lprob(sequence *s, size_t i)
{
	return - read_char(s, &_qual_sequence_opt, i) * log(10.) / 10.;
} /* qual_to_lprob */

/**
 * Probability of substitution from \par s to \par r with quality score \par q,
 * assuming XY_ENCODING.
 */
inline
double literal_sub_prob(data_t s, data_t r, data_t q, void *vptr)
{
	UNUSED(vptr);
	if (s == r)
		return 1 - exp( -q / 10. * log(10.) );
	else
		return exp( -q / 10. * log(10.) ) / 3;
} /* literal_sub_prob */

/**
 * Log probability of substitution from \par s to \par r with quality score
 * \par q assuming XY_ENCODING.
 */
inline
double literal_sub_lprob(data_t s, data_t r, data_t q, void *vptr)
{
	UNUSED(vptr);
	if (s == r)
		return log(1 - exp( -q * log(10) / 10 ));
	else
		return -q * log(10) / 10 - log(3);
} /* literal_sub_lprob */

double sub_prob_given_q_with_encoding(data_t s, data_t r, int es, int er, data_t q, int log, void *vptr);
double sub_prob_with_encoding(data_t s, data_t r, int es, int er, double e, int log);

/**
 * Convert character to quality.
 */
inline qual_t char_to_qual(unsigned char c)
{
	if (c >= MIN_ASCII_QUALITY_SCORE && c <= MAX_ASCII_QUALITY_SCORE)
{
//fprintf(stderr, "char_to_qual: %u\n", (c - MIN_ASCII_QUALITY_SCORE) & _q_mask);
		return (c - MIN_ASCII_QUALITY_SCORE) & _q_mask;
}
	return _q_mask;	/* impossible value */
} /* char_to_qual */

/**
 * Convert .
 */
inline qual_t char_to_censored_qual(unsigned char c, qual_t minq, qual_t maxq)
{
	if (((c - MIN_ASCII_QUALITY_SCORE) & _q_mask) < minq)
		return minq;
	if (((c - MIN_ASCII_QUALITY_SCORE) & _q_mask) > maxq)
		return maxq;

	return (c - MIN_ASCII_QUALITY_SCORE) & _q_mask;
} /* char_to_censored_qual */

/**
 * Get quality value at position \par i of sequence \par s.
 */
inline data_t get_qual(sequence *s, size_t i)
{
	return read_char(s, &_qual_sequence_opt, i);
} /* get_qual */

inline char qual_to_char(sequence *s, size_t i)
{
	return (char) read_char(s, &_qual_sequence_opt, i) + MIN_ASCII_QUALITY_SCORE;
}

/**
 * Write quality sequence to file pointer \par fp.
 */
inline void fwrite_qual_sequence(FILE *fp, sequence *s)
{
	for (unsigned int i = 0; i < s->len; ++i)
		fprintf(fp, "%c", (char) (read_char(s, &_qual_sequence_opt, i)
					+ MIN_ASCII_QUALITY_SCORE));
} /* fwrite_qual_sequence */

#endif
