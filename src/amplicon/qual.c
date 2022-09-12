#include <math.h>

#include "sequence.h"
#include "qual.h"
#include "nuc.h"

/*
const uint8_t _q_bits = 8;
const uint8_t _q_mask = 255;
#define _q_bits 8U
#define _q_mask 255U
*/

extern double literal_sub_prob(data_t s, data_t r, data_t q, void *);
extern double literal_sub_lprob(data_t s, data_t r, data_t q, void *);

/* default: substitution probabilities are */
sub_prob_fxn sub_prob = literal_sub_prob;
sub_lprob_fxn sub_lprob = literal_sub_lprob;

sequence_opt _qual_sequence_opt = {_q_bits, _q_mask};

/* inline functions defined in qual.h */
extern qual_t char_to_qual(unsigned char c);
extern qual_t char_to_censored_qual(unsigned char c, qual_t, qual_t);
extern void fwrite_qual_sequence(FILE *fp, sequence *s);
extern double qual_to_exp_err(sequence *s);
extern double prob_error_free(sequence *s);
extern double qual_to_prob(sequence *s, size_t i);
extern double qual_to_lprob(sequence *s, size_t i);
extern char qual_to_char(sequence *s, size_t i);
extern data_t get_qual(sequence *s, size_t i);

/**
 * Probability of substitution, allowing for different encodings and ambiguous
 * nucleotides.
 *
 * @param s		source nucleotide
 * @param r		read nucleotide
 * @param es		source nucleotide encoding
 * @param er		read nucleotide encoding
 * @param e		probability of error in read nucleotide
 * @param logged	return logged value
 * @return		substitution probability, possibly logged
 */
double sub_prob_with_encoding(data_t s, data_t r, int es, int er,
						double e, int logged)
{
	double ll = 0;
	if (es == XY_ENCODING && er == XY_ENCODING) {
		return s == r
			? logged ? log(1 - e) : 1 - e
			: logged ? log(e / 3) : e / 3;
	} else if (es == IUPAC_ENCODING && er == XY_ENCODING) {
		if (popcnt[s] == 1) {
			return iupac_to_xy[s] == r
				? logged ? log(1 - e) : 1 - e
				: logged ? log(e / 3) : e / 3;
		} else {
			for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
				if (s >> i & 1)
					ll += (iupac_to_xy[1 << i] == r
						? 1 - e : e / 3) / popcnt[s];
			if (logged)
				ll = log(ll);
		}
	} else if (es == XY_ENCODING && er == IUPAC_ENCODING) {
		if (popcnt[r] == 1)
			return s == iupac_to_xy[r]
				 ? logged ? log(1 - e) : 1 - 3
				 : logged ? log(e / 3) : e / 3;
		else
			for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
				if (r >> i & 1)
					ll += s == iupac_to_xy[1 << i]
						? 1 - e : e / 3;
		if (logged)
			ll = log(ll);
	} else {
		if (popcnt[s] == 1) {
			if (popcnt[r] == 1)
				return s == r
					? logged ? log(1 - e) : 1 - e
					: logged ? log(e / 3) : e / 3;
			else	/* read nucleotide ambiguous */
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
					if (r >> i & 1)
						ll += s == iupac_to_xy[1 << i]
							? 1 - e : e / 3;

			if (logged)
				ll = log(ll);
		} else {	/* source nucleotide ambiguous */
			if (popcnt[r] == 1) {
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
					if (s >> i & 1)
						ll += (iupac_to_xy[1 << i] == r
							? 1 - e : e / 3)
							/ popcnt[s];
				if (logged)
					ll = log(ll);
			} else {	/* read nucleotide ambiguous */
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i) {
					if (s >> i & 1) {
						double prob = 0;
						for (int j = 0; j
							< NUM_NUCLEOTIDES; ++j)
							if (r >> j & 1)
								prob += iupac_to_xy[1 << i] == iupac_to_xy[1 << j]
									? 1 - e : e / 3;
						ll += prob / popcnt[s];
					}
				}
				if (logged)
					ll = log(ll);
			}
		}
	}
	return ll;
} /* sub_prob_with_encoding */

/**
 * Probability of substitution, allowing for different encodings and ambiguous
 * nucleotides.
 *
 * @param s		source nucleotide
 * @param r		read nucleotide
 * @param es		source nucleotide encoding
 * @param er		read nucleotide encoding
 * @param q		quality score of read nucleotide
 * @param logged	return logged value
 * @param vptr		void pointer
 * @return		substitution probability, possibly logged
 */
double sub_prob_given_q_with_encoding(data_t s, data_t r, int es, int er,
					data_t q, int logged, void *vptr)
{
	double ll = 0;
	if (es == XY_ENCODING && er == XY_ENCODING) {
		return logged ? sub_lprob(s, r, q, vptr) : sub_prob(s, r, q, vptr);
	} else if (es == IUPAC_ENCODING && er == XY_ENCODING) {
		if (popcnt[s] == 1) {
			return logged ? sub_lprob(iupac_to_xy[s], r, q, vptr)
				: sub_prob(iupac_to_xy[s], r, q, vptr);
		} else {	/* assume uniform */
			for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
				if (s >> i & 1)
					ll += sub_prob(iupac_to_xy[1 << i], r, q, vptr)
								/ popcnt[s];
			if (logged)
				ll = log(ll);
		}
	} else if (er == IUPAC_ENCODING && es == XY_ENCODING) {
		if (popcnt[r] == 1)
			return logged ? sub_lprob(s, iupac_to_xy[r], q, vptr)
					 : sub_prob(s, iupac_to_xy[r], q, vptr);
		else 
			for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
				if (r >> i & 1)
					ll += sub_prob(s, iupac_to_xy[1 << i], q, vptr);
		if (logged)
			ll = log(ll);
	} else {
		if (popcnt[s] == 1) {
			if (popcnt[r] == 1)
				return logged ? sub_lprob(s, r, q, vptr) : sub_prob(s, r, q, vptr);
			else
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
					if (r >> i & 1)
						ll += sub_prob(s, iupac_to_xy[1 << i], q, vptr);
			if (logged)
				ll = log(ll);
		} else {
			if (popcnt[r] == 1) {
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i)
					if (s >> i & 1)
						ll += sub_prob(iupac_to_xy[1 << i], r, q, vptr) / popcnt[s];
				if (logged)
					ll = log(ll);
			} else {
				for (int i = 0; i < NUM_NUCLEOTIDES; ++i) {
					if (s >> i & 1) {
						double prob = 0;
						for (int j = 0; j < NUM_NUCLEOTIDES; ++j)
							if (r >> j & 1)
								prob += sub_prob(iupac_to_xy[1 << i], iupac_to_xy[1 << j], q, vptr);
						ll += prob / popcnt[s];
					}
				}
				if (logged)
					ll = log(ll);
			}
		}
	}
	return ll;
} /* sub_prob_given_q_with_encoding */

