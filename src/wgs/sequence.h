#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"

/**
 * The fundamental unit we use, which is bigger than all characters in
 * our sequences.
 */
//typedef uint8_t data_t;

/**
 * The new data types created here.
 */
typedef struct _sequence sequence;
typedef struct _sequence_opt sequence_opt;

/**
 * The size of our fundamental unit.
 */
extern size_t _data_size;
//size_t _data_size = 8 * sizeof(data_t);


/**
 * Structure to store sequence.
 */
struct _sequence {
	data_t *seq;	/* sequence */
	size_t len;	/* sequence length */
};

/**
 * Structure for sequence options.
 */
struct _sequence_opt {
	uint8_t bits;	/* bits for character in sequence */
	data_t mask;	/* mask for single character */
};

/**
 * Reach single character in sequence.
 *
 * @param s	sequence data
 * @param so	sequence options
 * @param i	index of character to retrieve
 * @return	the requested character as data_t type
 */
inline data_t read_char(sequence *s, sequence_opt *so, size_t i)
{
	size_t idx = i * so->bits;
	return (s->seq[ idx / _data_size ] >> ( idx % _data_size )) & so->mask;
} /* read_char */

/**
 * Write single character in sequence.
 *
 * @param s	sequence data
 * @param so	sequence options
 * @param i	index of character to overwrite
 * @param in	nucleotide to write
 */
/*@unused@*/ static inline void write_char(sequence *s, sequence_opt *so, size_t i, data_t in)
{
	size_t idx = i * so->bits;
	size_t shift = idx % _data_size;
//fprintf(stderr, "idx = %zu (%zu), shift = %u, data_size = %zu, in = %u\n", idx, idx / _data_size, shift, _data_size, in);
	idx /= _data_size;
	s->seq[idx] = (s->seq[idx] & ~(so->mask << shift)) | (in << shift);
} /* write_char */

/*@unused@*/ /*@null@*/ /*@out@*/ static inline data_t *sequence_alloc(size_t nchar, sequence_opt *so)
{
	return malloc((size_t) ceil((double) nchar / so->bits) * _data_size);
} /* sequence_alloc */

/*@unused@*/ static inline size_t seqlen(sequence *s, sequence_opt *so)
{
	return (size_t) ceil((double) s->len / so->bits);
} /* seqlen */

/*@unused@*/ static inline size_t seqbytes(sequence *s, sequence_opt *so)
{
	return (size_t) ceil((double) s->len / so->bits) * _data_size;
} /* seqbytes */

#endif
