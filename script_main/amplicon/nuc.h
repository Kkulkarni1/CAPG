/**
 * @file nuc.c
 * @author Karin S. Dorman
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
 * Beware: There must be internal consistency in the orderings encoding.  For
 * example, which bit of the 4-bit encoding represents A?  I use the order (high
 * bit) TGCA (low bit).
 * The xy_t encoding orders the nucleotides A=0, C=1, T=2, G=3, because it makes
 * for convenient conversion between the FASTA/FASTQ encoding and the 2-bit
 * encoding, but many people assume the alphabetic ordering A, C, G, T, or the
 * biochemical ordering G, A, C, T we used to use on sequencing gels (Why did we
 * do that exactly?).  The encodings are shown in gory detail in fastq.c.
 *
 */

#ifndef __H_NUC_DATA__
#define __H_NUC_DATA__

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <ctype.h>

#include "sequence.h"
#include "error.h"

/**
 * Types of nucleotide encodings.
 */
enum {
	DEFAULT_ENCODING,	/*!< default: xy_t unless need iupac_t */
	IUPAC_ENCODING,		/*!< iupac_t */
	XY_ENCODING,		/*!< xy_t */
	NUC_ENCODING		/*!< nuc_t: should not use */
};

/**
 * IUPAC symbols encoded in 4 lower bits.
 */
typedef unsigned char iupac_t;

/**
 * A, C, G, T as encoded by Xin Yin in 2 lower bits.
 */
typedef unsigned char xy_t;

/**
 * IUPAC symbols and 'X' encoded (with some gaps) as chars 0 to 24.
 */
typedef unsigned char nuc_t;

/**
 * Number of nucleotides: A, C, G, T.
 */
#define NUM_NUCLEOTIDES	4
#define NUM_IUPAC_SYMBOLS	16

/**
 * A, C, G, T as xy_t, iupac_t, and nuc_t.
 */
enum {
	XY_A = 0,
	XY_C = 1,
	XY_G = 3,
	XY_T = 2
};

enum {
	IUPAC_X = 0,
	IUPAC_A = 1,
	IUPAC_C = 2,
	IUPAC_G = 4,
	IUPAC_T = 8,
	IUPAC_U = 8,
	IUPAC_R = 5,
	IUPAC_Y = 10,
	IUPAC_S = 6,
	IUPAC_W = 9,
	IUPAC_K = 12,
	IUPAC_M = 3,
	IUPAC_B = 14,
	IUPAC_D = 13,
	IUPAC_H = 11,
	IUPAC_V = 7,
	IUPAC_N = 15
};
enum {
	STD_A,	/* 0 */
	STD_C,	/* 1 */
	STD_G,	/* 2 */
	STD_T	/* 3 */
};

/**
 * Minimum allowed nucleotide sequence character in ASCII is char 'A'.
 */
#define MIN_NUCLEOTIDE_ASCII	'A'

/**
 * Maximum number of letters in ASCII nucleotide alphabet: 'A' to 'Y'
 */
#define NUCLEOTIDE_ALPHABET_SIZE	('Z' - 'A')

/**
 * Convert xy_t, std_t, iupac_t to display char.
 */
extern unsigned char const xy_to_char[NUM_NUCLEOTIDES];
extern unsigned char const std_to_char[NUM_NUCLEOTIDES];
extern unsigned char const iupac_to_char[NUM_IUPAC_SYMBOLS];

/**
 * For sequences of various types.
 */
/*
extern const uint8_t _xy_bits;
extern const uint8_t _xy_mask;
*/
extern sequence_opt _xy_sequence_opt;
/*
extern const uint8_t _iupac_bits;
extern const uint8_t _iupac_mask;
*/
extern sequence_opt _iupac_sequence_opt;

/**
 * Convert among iupac_t, xy_t, std, and nuc_t.
 */
extern unsigned char const iupac_to_xy[NUM_IUPAC_SYMBOLS];
extern unsigned char const iupac_to_std[NUM_IUPAC_SYMBOLS];

extern unsigned char const xy_to_std[NUM_NUCLEOTIDES];
extern unsigned char const xy_to_iupac[NUM_NUCLEOTIDES];
extern unsigned char const std_to_xy[NUM_NUCLEOTIDES];
extern iupac_t const nuc_to_iupac[NUCLEOTIDE_ALPHABET_SIZE];
extern xy_t const nuc_to_xt[NUCLEOTIDE_ALPHABET_SIZE];

/**
 * Reverse complement the codes.
 */
extern xy_t const xy_to_rc[NUM_NUCLEOTIDES];
extern iupac_t const iupac_to_rc[NUM_IUPAC_SYMBOLS];

#define NUM_IUPAC_SYMBOLS	16
extern unsigned char iupac_symbols[NUM_IUPAC_SYMBOLS];
extern const unsigned char popcnt[NUM_IUPAC_SYMBOLS];

/**
 * Check whether character is valid IUPAC symbol.  Incoming character is
 * human-readable letter and output is 0/1 to indicate if it is one of the
 * allowed IUPAC symbols.
 *
 * @param c	human-readable character
 * @return	0|1
 */
inline int valid_iupac(char *c) {
	*c = toupper(*c);
	for (size_t i = 0; i < NUM_IUPAC_SYMBOLS; ++i)
		if (iupac_to_char[i] == *c)
			return 1;
	if (*c == 'X')
		return 1;
	return 0;
} /* valid_iupac */

/**
 * Check whether character is valid nucleotide.  Incoming character is
 * human-readable letter and output is 0/1 to indicate if it is one of
 * A, C, G, or T.
 *
 * @param c	human-readable character
 * @return	1|0
 */
inline int valid_nucleotide(char *c) {
	*c = toupper(*c);
	return *c == 'A' || *c == 'C' || *c == 'G' || *c == 'T' || *c == 'U';
} /* valid_nucleotide */

/**
 * Convert char to xy.
 *
 * @param c	ASCII char
 * @return	xy_t
 */
inline xy_t char_to_xy(char c)
{
	if (c == 'A' || c == 'a')
		return 'A' >> 1 & 3L;
	else if (c == 'C' || c == 'c')
		return 'C' >> 1 & 3L;
	else if (c == 'G' || c == 'g')
		return 'G' >> 1 & 3L;
	else if (c == 'T' || c == 't' || c == 'U' || c == 'u')
		return 'T' >> 1 & 3L;
	else
		return 1 << 7;	/* non-nuc */
} /* char_to_xy */

/**
 * Convert char to iupac_t.
 *
 * @apram c	ASCII char
 * @return	iupac_t
 */
inline iupac_t char_to_iupac(char c)
{
	if (c == 'A' || c == 'a')
		return 1L;
	else if (c == 'C' || c == 'c')
		return 2L;
	else if (c == 'G' || c == 'g')
		return 4L;
	else if (c == 'T' || c == 't' || c == 'U' || c == 'u')
		return 8L;
	else if (c == 'R' || c == 'r')
		return 5L;
	else if (c == 'Y' || c == 'y')
		return 10L;
	else if (c == 'S' || c == 's')
		return 6L;
	else if (c == 'W' || c == 'w')
		return 9L;
	else if (c == 'K' || c == 'k')
		return 12L;
	else if (c == 'M' || c == 'm')
		return 3L;
	else if (c == 'B' || c == 'b')
		return 14L;
	else if (c == 'D' || c == 'd')
		return 13L;
	else if (c == 'H' || c == 'h')
		return 11L;
	else if (c == 'V' || c == 'v')
		return 7L;
	else if (c == 'N' || c == 'N')
		return 15L;
	else
		return 0L;
} /* char_to_iupac */


/**
 * Convert char to encoding, returned as data_t.
 *
 * @param c		ASCII character
 * @param encoding	encoding to use
 * @return		encoded character as data_t
 */
inline data_t char_to_data(char c, int encoding)
{
	if (encoding == XY_ENCODING)
		return char_to_xy(c);
	else if (encoding == IUPAC_ENCODING)
		return char_to_iupac(c);
	else
		return (1 << (sizeof(data_t) - 1)) - 1;
} /* char_to_data */

/**
 * Write a xy-encoded nucleotide sequence to specified file stream.
 *
 * @param fp	file stream
 * @param s	sequence
 * @param sidx	start index (inclusive)
 * @param eidx	end index (not included)
 */
inline void fwrite_xy_sequence(FILE *fp, sequence *s, size_t sidx, size_t eidx)
{
	for (size_t i = sidx; i < eidx; ++i)
		fprintf(fp, "%c", xy_to_char[read_char(s,
					&_xy_sequence_opt, i)]);
} /* fwrite_xy_sequence */


/**
 * Write iupac-encoded nucleotide sequence to specified file stream.
 *
 * @param fp	file stream
 * @param s	sequence
 * @param sidx	start index (inclusive)
 * @param eidx	end index (not included)
 */
inline void fwrite_iupac_sequence(FILE *fp, sequence *s, size_t sidx,
							size_t eidx)
{
	for (size_t i = sidx; i < eidx; ++i)
		fprintf(fp, "%c", iupac_to_char[read_char(s,
						 &_iupac_sequence_opt, i)]);
} /* fwrite_iupac_sequence */

/**
 * Write a nucleotide sequence to specified file pointer.
 *
 * @param fp		file stream to write to (may be stdout or stderr)
 * @param s		sequence to write
 * @param encoding	how is this sequence encoded?
 * @return		error status
 */
inline int fwrite_nuc_sequence(FILE *fp, sequence *s, int encoding)
{
	if (encoding == XY_ENCODING)
		fwrite_xy_sequence(fp, s, 0, s->len);
	else if (encoding == IUPAC_ENCODING)
		fwrite_iupac_sequence(fp, s, 0, s->len);
	else
		return 1;
	return 0;
} /* fwrite_nuc_sequence */

/**
 * Write segment of nucleotide sequence to file.
 *
 * @param fp		file stream (may be stdout or stderr)
 * @param s		sequence to write
 * @param encoding	how is sequence encoded
 * @param sidx		start index
 * @param eidx		end index
 * @return		error status
 */
inline int fwrite_nuc_segment(FILE *fp, sequence *s, int encoding,
					size_t sidx, size_t eidx)
{
	if (encoding == XY_ENCODING)
		fwrite_xy_sequence(fp, s, sidx, eidx);
	else if (encoding == IUPAC_ENCODING)
		fwrite_iupac_sequence(fp, s, sidx, eidx);
	else
		return 1;
	return 0;
} /* fwrite_nuc_segment */

/**
 * Return appropiate sequence options for requested encoding.  Sequences
 * are passed around with options that specify their internals.  Rather
 * than ask the user to construct these options, we have options for the
 * standard encodings accessible via this function.
 *
 * @param encoding	desired encoding
 * @return		pointer sequence_opt object
 */
inline sequence_opt *nuc_sequence_opt(int encoding)
{
	if (encoding == XY_ENCODING)
		return &_xy_sequence_opt;
	else if (encoding == IUPAC_ENCODING)
		return &_iupac_sequence_opt;
	else
		return NULL;
} /* nuc_sequence_opt */

inline data_t get_nuc(sequence *s, int encoding, size_t i)
{
	return read_char(s, nuc_sequence_opt(encoding), i);
} /* get_nuc */

#endif
